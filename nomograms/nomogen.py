#!/usr/bin/env python3
"""
Nomogen - Numerical Nomogram Generator
===============================================================================
A high-performance automated generator for nomograms using numerical relaxation 
techniques.

Summary
-------
This module generates non-standard nomograms by solving a constrained optimization 
problem. It treats the nomogram scales as elastic threads and relaxes them into a 
configuration that minimizes alignment error (geometric accuracy) while maintaining 
smoothness (minimizing internal strain/curvature).

The generator is capable of creating:
1. Standard 3-scale alignment charts.
2. Curved scale nomograms (Type 9).
3. Ladder charts (Type 8 - implicitly supported via dual alignment).

Key Features
------------
* **Numerical Relaxation**: Uses a custom cost function combining alignment error 
    (the determinant of the alignment points) and smoothness constraints (derivatives).
* **Numba Optimization**: The core cost function and interpolation routines are 
    JIT-compiled using Numba to achieve C-like performance during the optimization 
    loop.
* **Chebyshev Nodes**: Uses Chebyshev spacing for discretization to minimize 
    Runge's phenomenon (oscillation at the edges of polynomial interpolation).
* **Barycentric Interpolation**: Implements stable barycentric Lagrange 
    interpolation for high-precision function approximation across the scales.

Dependencies
------------
* numpy: Matrix operations and array handling.
* scipy: Optimization routines (L-BFGS-B) and interpolation helpers.
* numba: JIT compilation for performance-critical kernels.
===============================================================================
"""

import sys
import math
import functools
import warnings
from typing import Dict, Any, Callable, Tuple, List, Optional, Union

import numpy as np
import scipy.interpolate
import scipy.optimize

# Check for Numba and import JIT decorators
try:
    from numba import njit, prange
    HAS_NUMBA = True
except ImportError:
    print("CRITICAL ERROR: Numba is required for this module.")
    print("Please install via: pip install numba")
    sys.exit(1)


# =============================================================================
#  SECTION 1: JIT-Compiled Numerical Kernels
# =============================================================================
#  These functions form the computational engine of the generator. They are 
#  compiled to machine code by Numba for maximum performance during the 
#  optimization loop.
# =============================================================================

@njit(cache=True)
def _get_barycentric_weights(nodes: np.ndarray) -> np.ndarray:
    """
    Calculate the weights for the second form of barycentric Lagrange interpolation.
    
    Time Complexity: O(N^2)
    
    Parameters:
        nodes (np.ndarray): The x-coordinates of the interpolation nodes.
        
    Returns:
        np.ndarray: Array of weights corresponding to each node.
    """
    n = len(nodes)
    weights = np.ones(n, dtype=np.float64)
    for j in range(n):
        for k in range(n):
            if j != k:
                weights[j] /= (nodes[j] - nodes[k])
    return weights


@njit(cache=True)
def _barycentric_interp_scalar(x: float, nodes: np.ndarray, weights: np.ndarray, y_vals: np.ndarray) -> float:
    """
    Perform barycentric interpolation for a single scalar value x using the 
    precomputed weights. Handles the singularity where x matches a node exactly.
    """
    epsilon = 1e-14
    # Check for exact node match to avoid division by zero
    for i in range(len(nodes)):
        if abs(x - nodes[i]) < epsilon:
            return y_vals[i]

    num = 0.0
    den = 0.0
    for i in range(len(nodes)):
        w_curr = weights[i] / (x - nodes[i])
        num += w_curr * y_vals[i]
        den += w_curr
    return num / den


@njit(parallel=True, cache=True)
def _barycentric_interp_grid(grid_flat: np.ndarray, nodes: np.ndarray, weights: np.ndarray, y_vals: np.ndarray) -> np.ndarray:
    """
    Perform barycentric interpolation for a flattened array of query points in parallel.
    
    Parameters:
        grid_flat (np.ndarray): Flattened array of x-coordinates to interpolate at.
        nodes (np.ndarray): The known x-nodes.
        weights (np.ndarray): Precomputed barycentric weights.
        y_vals (np.ndarray): The known y-values at the nodes.
        
    Returns:
        np.ndarray: Interpolated values corresponding to grid_flat.
    """
    n = len(grid_flat)
    res = np.empty(n, dtype=np.float64)
    
    for k in prange(n):
        x = grid_flat[k]
        val = 0.0
        match_idx = -1
        
        # Check for exact node match
        for i in range(len(nodes)):
            if abs(x - nodes[i]) < 1e-14:
                match_idx = i
                break
        
        if match_idx != -1:
            val = y_vals[match_idx]
        else:
            num = 0.0
            den = 0.0
            for i in range(len(nodes)):
                w = weights[i] / (x - nodes[i])
                num += w * y_vals[i]
                den += w
            val = num / den
        res[k] = val
    return res


@njit(parallel=True, cache=True)
def _calc_cost_core(
    lxu: np.ndarray, lyu: np.ndarray, 
    lxv: np.ndarray, lyv: np.ndarray, 
    unodes: np.ndarray, vnodes: np.ndarray, 
    txwa_flat: np.ndarray, tywa_flat: np.ndarray, 
    dxdwa_flat: np.ndarray, dydwa_flat: np.ndarray,
    fdxdu: np.ndarray, fdydu: np.ndarray, 
    fdxdv: np.ndarray, fdydv: np.ndarray,
    dwdu_flat: np.ndarray, dwdv_flat: np.ndarray,
    cost_pos_u: np.ndarray, cost_pos_v: np.ndarray,
    muShape: float, resolution: float, 
    umax: float, umin: float, vmax: float, vmin: float
) -> Tuple[float, float, float, float]:
    """
    The computationally intensive inner loop of the cost function.
    
    Calculates the error accumulators for:
    1. Alignment Error (Geometric validity): do the 3 points form a line?
    2. Derivative Error (Smoothness): are the scales changing abruptly?
    3. Shape Error (Constraint): enforces placement within the paper bounds.
    
    This function iterates over the grid of (u, v) pairs.
    """
    local_eAcc = 0.0
    local_eDeru = 0.0
    local_eDerv = 0.0
    local_eShape = 0.0

    n_u = len(unodes)
    n_v = len(vnodes)
    
    # Flattened loop for easier parallelization
    for i in prange(n_u * n_v):
        # Decode index back to 2D
        iu = i // n_v
        iv = i % n_v
        
        u = unodes[iu]
        v = vnodes[iv]

        # Fetch U scale properties at this node
        txu = lxu[iu]
        tyu = lyu[iu]
        dxdu = fdxdu[iu]
        dydu = fdydu[iu]
        
        # Fetch V scale properties at this node
        txv = lxv[iv]
        tyv = lyv[iv]
        dxdv = fdxdv[iv]
        dydv = fdydv[iv]
        
        # Fetch W scale interpolated properties (pre-calculated outside loop)
        txw = txwa_flat[i]
        tyw = tywa_flat[i]
        dxdw = dxdwa_flat[i]
        dydw = dydwa_flat[i]
        
        # --- Alignment Error (eAcc) ---
        # Calculate distance between U and V points
        tx = txu - txv
        ty = tyu - tyv
        td2 = tx * tx + ty * ty
        td = math.sqrt(td2)
        
        # Collinearity check (determinant form / distance)
        # If points are collinear, this value should be zero.
        e0 = (tx * (tyu - tyw) - (txu - txw) * ty) / td
        
        if td2 * resolution > 1.0:
            local_eAcc += e0 * e0

        # --- Derivative/Gradient Error Logic ---
        # Checks if the scales are smooth relative to the function change
        tmp = ty * dxdw - tx * dydw
        tuc = e0 * (tx * dxdu + ty * dydu) / td2
        tvc = e0 * (tx * dxdv + ty * dydv) / td2
        
        # Calculate contribution from dw/du
        dwdu = dwdu_flat[i]
        dedu = ((dwdu * tmp + (tyv - tyw) * dxdu - (txv - txw) * dydu)) / td - tuc
        local_eDeru += dedu * dedu

        # Calculate contribution from dw/dv
        dwdv = dwdv_flat[i]
        dedv = ((dwdv * tmp + (tyw - tyu) * dxdv - (txw - txu) * dydv)) / td + tvc
        local_eDerv += dedv * dedv
        
        # --- Shape Error (Scale Positioning Constraints) ---
        # Enforces the scales to stay roughly within the defined box
        if muShape != 0:
            t_u = 0.01 * (umax - umin)**2 + u * u
            uShape = ((dxdu * ty - dydu * tx))**2 * t_u
            
            t_v = 0.01 * (vmax - vmin)**2 + v * v
            vShape = ((dxdv * ty - dydv * tx))**2 * t_v
            
            inv_sum = 0.0
            if uShape > 1e-12: inv_sum += 1.0/uShape
            if vShape > 1e-12: inv_sum += 1.0/vShape
            
            tShape = td2 * inv_sum
            tShape *= (cost_pos_u[iu] + cost_pos_v[iv])
            local_eShape += tShape

    return local_eAcc, local_eDeru, local_eDerv, local_eShape


# =============================================================================
#  SECTION 2: Main Nomogen Class
# =============================================================================

class Nomogen:
    """
    Manages the numerical generation of a nomogram.
    
    This class handles the entire lifecycle:
    1. Parsing configuration and establishing scale transformations.
    2. Discretizing the domain using Chebyshev nodes.
    3. Precomputing derivatives.
    4. Running the Scipy optimization loop.
    5. Verifying results and injecting the resulting geometry back into parameters.
    """

    def __init__(self, func: Callable, main_params: Dict[str, Any]):
        """
        Initialize the Nomogram Generator.

        Args:
            func: The Python function f(u, v) solving for w.
            main_params: The dictionary configuration (standard PyNomo format).
        """
        self.func = func
        self.main_params = main_params
        
        # Optimization state tracking
        self.iteration_count = 0
        self.current_cost = 0.0
        self.old_cost = 0.0
        self.metrics = {'eAcc': 0.0, 'eDeru': 0.0, 'eDerv': 0.0, 'eShape': 0.0}

        # Main execution pipeline
        self._validate_and_setup_params()
        self._initialize_grids()
        self._precompute_derivatives()
        self._optimize()
        self._apply_and_verify_results()

    # -------------------------------------------------------------------------
    # Initialization & Validation
    # -------------------------------------------------------------------------

    def _validate_and_setup_params(self):
        """
        Validates inputs, sets up U/V/W limits, and creates log/linear wrappers
        for transformation functions.
        """
        # =====================================================================
        # FIXED: Define 'params' FIRST to avoid NameError
        # =====================================================================
        params = self.main_params

        # =====================================================================
        # FORCE A3 PAGE SIZE (29.7cm x 42.0cm)
        # =====================================================================
        # Overwrite whatever was passed in parameters to ensure consistent A3 output.
        # This is a hard override requested by the user.
        params['paper_width'] = 29.7
        params['paper_height'] = 42.0
        print(" >> FORCING A3 PAGE SIZE (29.7cm x 42.0cm)")
        # =====================================================================

        # Ensure correct block type
        if params['block_params'][0]['block_type'] != 'type_9':
            sys.exit(f"'type_9' block expected, found '{params['block_params'][0]['block_type']}'")

        # Setup N (Number of Chebyshev Points)
        if 'pdegree' in params:
            self.NN = params['pdegree']
            print("\nWARNING: 'pdegree' renamed to 'npoints'. Please update your code.")
        else:
            self.NN = params.get('npoints', 9)

        if not isinstance(self.NN, int):
            print(f"\n'npoints' must be integer, found {type(self.NN)} - ", end=' ')
            self.NN = 9
        elif self.NN < 3:
            sys.exit("npoints must be >= 3")
        print("Using", self.NN, "Chebyshev points")

        # Shape Parameter (regularization)
        self.muShape = 1.0e-14 if 'muShape' in params else 0
        if self.muShape != 0:
            print('Setting nomogram for best measureability, muShape is', self.muShape)

        # --- Configure U Scale ---
        self.params_u = params['block_params'][0]['f1_params']
        umin_user, umax_user = self.params_u['u_min'], self.params_u['u_max']
        if umax_user < umin_user:
            sys.exit(f"error: umax ({umax_user}) < umin ({umin_user})")

        if 'log' in self.params_u.get('scale_type', ''):
            if umin_user <= 0: sys.exit(f"error: umin {umin_user} <= 0 for log scale")
            self.u_user_trans = math.exp
            self.u_plot_trans = math.log
        else:
            self.u_user_trans = lambda t: t
            self.u_plot_trans = lambda t: t

        # --- Configure V Scale ---
        self.params_v = params['block_params'][0]['f3_params']
        vmin_user, vmax_user = self.params_v['u_min'], self.params_v['u_max']
        if vmax_user < vmin_user:
            sys.exit(f"error: vmax ({vmax_user}) < vmin ({vmin_user})")

        if 'log' in self.params_v.get('scale_type', ''):
            if vmin_user <= 0: sys.exit(f"error: vmin {vmin_user} <= 0 for log scale")
            self.v_user_trans = math.exp
            self.v_plot_trans = math.log
        else:
            self.v_user_trans = lambda t: t
            self.v_plot_trans = lambda t: t

        # --- Configure W Scale (Result) ---
        self.params_w = params['block_params'][0]['f2_params']
        wmin_in, wmax_in = self.params_w['u_min'], self.params_w['u_max']

        # Determine actual W range based on U and V corners
        corners = [self.func(u, v) for u in [umin_user, umax_user] for v in [vmin_user, vmax_user]]
        wval = [wmin_in, wmax_in] + corners
        wmin_user, wmax_user = min(wval), max(wval)

        if 'log' in self.params_w.get('scale_type', ''):
            if wmin_user <= 0: sys.exit(f"error: wmin {wmin_user} <= 0 for log scale")
            self.w_user_trans = math.exp
            self.w_plot_trans = math.log
        else:
            self.w_user_trans = lambda t: t
            self.w_plot_trans = lambda t: t

        # Store Transformed Limits (internal working coordinates)
        self.umin, self.umax = self.u_plot_trans(umin_user), self.u_plot_trans(umax_user)
        self.vmin, self.vmax = self.v_plot_trans(vmin_user), self.v_plot_trans(vmax_user)
        self.wmin, self.wmax = self.w_plot_trans(wmin_user), self.w_plot_trans(wmax_user)

        # Transformed Function Wrapper
        def w_func_impl(u, v):
            """Internal wrapper handling log/linear conversions."""
            return self.w_plot_trans(self.func(self.u_user_trans(u), self.v_user_trans(v)))
        self.w_func = w_func_impl

        # Physical Properties (Resolution used for error thresholds)
        self.width = 10 * params['paper_width']
        self.height = 10 * params['paper_height']
        self.resolution = 10 * 10 * self.width * self.height

    def _initialize_grids(self):
        """
        Sets up the Chebyshev grid nodes, barycentric weights, and generates
        the initial geometric estimates for the scales.
        """
        # Chebyshev Grid Generation
        # Grid clustering near ends reduces oscillation
        self.chebyGrid = (1 - np.cos(np.linspace(0, math.pi, self.NN))) / 2
        
        self.unodes = self.umin + (self.umax - self.umin) * self.chebyGrid
        self.vnodes = self.vmin + (self.vmax - self.vmin) * self.chebyGrid
        self.wnodes = self.wmin + (self.wmax - self.wmin) * self.chebyGrid

        # Barycentric Weights (Required for Numba interpolation)
        self.u_weights = _get_barycentric_weights(self.unodes)
        self.v_weights = _get_barycentric_weights(self.vnodes)
        self.w_weights = _get_barycentric_weights(self.wnodes)

        # --- Initial Geometry Estimation ---
        # We perform a rough layout to give the optimizer a valid starting point.
        w0, w1 = self.w_func(self.umin, self.vmin), self.w_func(self.umax, self.vmin)
        w2, w3 = self.w_func(self.umin, self.vmax), self.w_func(self.umax, self.vmax)

        # Determine vertical orientation based on w-values
        if (w1 > w0) == (w2 > w0):
            # vmax is at top
            wB, wE = w0, w3
            wG, wH = w1, w2
            yv0, yv2 = 0, 1
        else:
            # vmax is at bottom
            wB, wE = w2, w1
            wG, wH = w3, w0
            yv0, yv2 = 1, 0

        yw0, yw2 = (0, 1) if wE > wB else (1, 0)

        # Initialize Coordinate Arrays (Normalized 0-1)
        self.xu = np.zeros(self.NN)
        self.yu = self.chebyGrid.copy()
        self.xv = np.ones(self.NN)
        self.yv = yv0 + (yv2 - yv0) * self.chebyGrid

        # Initial W scale estimation
        if math.isclose(wG, wH):
            self.xw = np.full(self.NN, 0.5)
            yw_interp = scipy.interpolate.interp1d([wB, wG, wE], [0, 0.5, 1])
            self.yw = np.array([yw_interp(w) for w in self.wnodes])
        else:
            # Alpha approximation for crossed scales
            alphaxw = (wE - wG - wH + wB) / (wE - wB) / (wG - wH)
            xwB = (wH - wB) / (wE - wB) - alphaxw * (wH - wB)
            xwE = xwB + alphaxw * (wE - wB)
            xwB, xwE = np.clip([xwB, xwE], 0, 1)
            self.xw = xwB + self.chebyGrid * (xwE - xwB)
            self.yw = yw0 + (yw2 - yw0) * self.chebyGrid

    def _precompute_derivatives(self):
        """
        Calculates differentiation matrices and partial derivatives of W.
        These are static throughout the optimization.
        """
        
        def diffmat(x):
            """Generate differentiation matrix for a set of nodes x."""
            n = len(x)
            if n == 0: return 0
            XX = np.tile(x, (n, 1))
            c = np.empty((n, 1))
            c[::2] = 1
            c[1::2] = -1
            c[0] = 2; c[-1] *= 2
            D = c * (1 / c).T
            dXX = XX.T - XX + np.eye(n)
            D = D / dXX
            for ii in range(n):
                D[ii][ii] -= D[ii].sum()
            return D

        def derivative(f, x, x1, x0):
            """Numerical derivative using 4-point stencil."""
            try:
                h = 1.0e-4 * (x1 - x0)
                r = (-f(x + 2 * h) + 8 * f(x + h) - 8 * f(x - h) + f(x - 2 * h)) / (12 * h)
                if math.isnan(r):
                    print(f'Derivative NaN at x={x}, result={r}')
            except ValueError:
                r = math.nan
                print(f'Derivative exception at x={x}')
            return r

        # Generate Differentiation Matrices
        self.diffmatu = diffmat(self.unodes)
        self.diffmatv = diffmat(self.vnodes)
        self.diffmatw = diffmat(self.wnodes)

        # Pre-calculate W values grid
        self.w_values = np.array([[self.w_func(u, v) for v in self.vnodes] for u in self.unodes])
        self.w_values_flat = self.w_values.ravel()

        # Partial Derivatives (dw/du, dw/dv) calculation
        dwdu_values = np.array([[derivative(lambda uu: self.w_func(uu, v), u, self.umax, self.umin) 
                                 for v in self.vnodes] for u in self.unodes])
        
        dwdv_values = np.array([[derivative(lambda vv: self.w_func(u, vv), v, self.vmax, self.vmin) 
                                 for v in self.vnodes] for u in self.unodes])

        # Validate derivatives
        if np.isnan(dwdu_values).any() or np.isnan(dwdv_values).any():
            sys.exit("\nError: The scales cannot be represented by a finite polynomial in this range.\n"
                     "Derivatives contain NaNs. Please reduce limits on the scales and try again.")

        self.dwdu_flat = dwdu_values.ravel()
        self.dwdv_flat = dwdv_values.ravel()
        # Keep 2D versions for verification later
        self.dwdu_values = dwdu_values
        self.dwdv_values = dwdv_values

    # -------------------------------------------------------------------------
    # Optimization Loop
    # -------------------------------------------------------------------------

    def _cost_callback(self, xk):
        """Callback for Scipy Optimizer to print progress in-place."""
        self.iteration_count += 1
        
        # Shorten filename to first 10 chars (ASCII safe)
        fname = self.main_params.get('filename', 'nomogram')
        if len(fname) > 10:
            fname = fname[:9] + "~"  # Use tilde instead of ellipsis to prevent encoding errors
            
        # Compact formatting to prevent line-wrapping
        msg = (f"\r{fname} "
               f"It:{self.iteration_count:<4d} "
               f"Ac:{self.metrics['eAcc']:.1e} "
               f"dU:{self.metrics['eDeru']:.1e} "
               f"dV:{self.metrics['eDerv']:.1e}")
        
        if self.muShape != 0:
            msg += f" Sh:{self.metrics['eShape']:.1e}"
            
        # Write to stdout and flush immediately
        sys.stdout.write(msg)
        sys.stdout.flush()
        
        self.old_cost = self.current_cost

    def _calculate_cost(self, variables_vector):
        """The Objective Function."""
        if self.muShape != 0:
            lxu, lyu, lxv, lyv, lxw, lyw = np.array_split(variables_vector, 6)
            aa = 40
            cost_pos_u = np.exp(aa*(-0.3-lxu)) + np.exp(aa*(lxu-1.3)) + \
                         np.exp(aa*(-0.3-lyu)) + np.exp(aa*(lyu-1.3)) + 1
            cost_pos_v = np.exp(aa*(-0.3-lxv)) + np.exp(aa*(lxv-1.3)) + \
                         np.exp(aa*(-0.3-lyv)) + np.exp(aa*(lyv-1.3)) + 1
        else:
            lxu, lyu, lxv, lyv, lyw = self.xu.copy(), self.yu.copy(), self.xv.copy(), self.yv.copy(), self.yw.copy()
            lxw = self.xw 
            idx_splits = [self.NN-2, 2*self.NN-4, 3*self.NN-6, 4*self.NN-8, 5*self.NN-8]
            parts = np.array_split(variables_vector, idx_splits)
            lxu[1:-1], lyu[1:-1], lxv[1:-1], lyv[1:-1], lxw, lyw[1:-1] = parts
            cost_pos_u = np.zeros_like(lxu)
            cost_pos_v = np.zeros_like(lxv)

        # Spatial Derivatives
        fdxdu = self.diffmatu @ lxu
        fdydu = self.diffmatu @ lyu
        fdxdv = self.diffmatv @ lxv
        fdydv = self.diffmatv @ lyv
        fdxdw = self.diffmatw @ lxw
        fdydw = self.diffmatw @ lyw

        # Interpolate W scale properties
        txwa_flat = _barycentric_interp_grid(self.w_values_flat, self.wnodes, self.w_weights, lxw)
        tywa_flat = _barycentric_interp_grid(self.w_values_flat, self.wnodes, self.w_weights, lyw)
        dxdwa_flat = _barycentric_interp_grid(self.w_values_flat, self.wnodes, self.w_weights, fdxdw)
        dydwa_flat = _barycentric_interp_grid(self.w_values_flat, self.wnodes, self.w_weights, fdydw)

        # Numba Optimized Cost Loop
        eAcc, eDeru, eDerv, eShape = _calc_cost_core(
            lxu, lyu, lxv, lyv,
            self.unodes, self.vnodes,
            txwa_flat, tywa_flat, dxdwa_flat, dydwa_flat,
            fdxdu, fdydu, fdxdv, fdydv,
            self.dwdu_flat, self.dwdv_flat,
            cost_pos_u, cost_pos_v,
            self.muShape, self.resolution,
            self.umax, self.umin, self.vmax, self.vmin
        )

        eDeru *= (self.umax - self.umin) ** 2
        eDerv *= (self.vmax - self.vmin) ** 2
        eShape *= self.muShape

        self.metrics = {'eAcc': eAcc, 'eDeru': eDeru, 'eDerv': eDerv, 'eShape': eShape}
        
        mu = 1.0e-3
        total_cost = (eAcc + mu * (eDeru + eDerv) + eShape) / (len(self.unodes) * len(self.vnodes))
        
        self.current_cost = total_cost
        return total_cost

    def _optimize(self):
        """
        Runs the Scipy L-BFGS-B optimizer to minimize the cost function.
        Updates self.xu, self.yu, etc. with the optimized geometry.
        """
        if self.muShape != 0:
            x0 = np.concatenate([self.xu, self.yu, self.xv, self.yv, self.xw, self.yw])
        else:
            # If shape is not optimized, we lock the endpoints of U and V scales
            # to prevent rigid body rotation drift.
            x0 = np.concatenate([
                self.xu[1:-1], self.yu[1:-1], 
                self.xv[1:-1], self.yv[1:-1], 
                self.xw, self.yw[1:-1]
            ])

        # L-BFGS-B Settings
        # High maxiter/maxfun to ensure deep convergence.
        # Strict gtol/ftol for precision.
        res = scipy.optimize.minimize(
            self._calculate_cost, 
            x0, 
            method='L-BFGS-B',
            callback=self._cost_callback,
            options={
                'disp': False, 
                'gtol': 1e-15, 
                'ftol': 1e-15, 
                'maxiter': 500000,
                'maxfun': 500000
            }
        )

        print()
        print(f"cost function is {res.fun:.2e}, {res.message}")

        # Unpack results back into coordinate arrays
        if self.muShape != 0:
            self.xu, self.yu, self.xv, self.yv, self.xw, self.yw = np.array_split(res.x, 6)
        else:
            self.xu[1:-1], self.yu[1:-1], self.xv[1:-1], self.yv[1:-1], self.xw, self.yw[1:-1] = \
                np.array_split(res.x, [self.NN-2, 2*self.NN-4, 3*self.NN-6, 4*self.NN-8, 5*self.NN-8])

    # -------------------------------------------------------------------------
    # Results, Verification, and Injection
    # -------------------------------------------------------------------------

    def _apply_and_verify_results(self):
        """
        Performs post-optimization verification and injects the generated scale 
        functions back into the main parameters dictionary for plotting.
        """
        
        # --- Part 1: Derivative Error Verification ---
        self._verify_derivatives()

        # --- Part 2: Dense Geometric Verification ---
        self._verify_alignment_dense()

        # --- Part 3: Inject Functions into Params ---
        self._inject_results()
        
        # --- Part 4: Auto-configure Ticks ---
        self._configure_tick_placement()

        # --- Part 5: Handle Dual Scales ---
        self._configure_dual_scales()

    def _verify_derivatives(self):
        """Checks consistency between numerical derivatives and geometry."""
        fdxdu = self.diffmatu @ self.xu
        fdydu = self.diffmatu @ self.yu
        fdxdv = self.diffmatv @ self.xv
        fdydv = self.diffmatv @ self.yv
        fdxdw = self.diffmatw @ self.xw
        fdydw = self.diffmatw @ self.yw

        maxerr = 0.0

        # Create temporary grid interpolators
        xwcoorda = _barycentric_interp_grid(self.w_values_flat, self.wnodes, self.w_weights, self.xw).reshape(self.NN, self.NN)
        ywcoorda = _barycentric_interp_grid(self.w_values_flat, self.wnodes, self.w_weights, self.yw).reshape(self.NN, self.NN)
        dxdwa = _barycentric_interp_grid(self.w_values_flat, self.wnodes, self.w_weights, fdxdw).reshape(self.NN, self.NN)
        dydwa = _barycentric_interp_grid(self.w_values_flat, self.wnodes, self.w_weights, fdydw).reshape(self.NN, self.NN)

        for i in range(len(self.unodes)):
            xucoord = self.xu[i]
            yucoord = self.yu[i]
            dxdu = fdxdu[i]
            dydu = fdydu[i]

            for j in range(len(self.vnodes)):
                xvcoord = self.xv[j]
                yvcoord = self.yv[j]
                dxdv = fdxdv[j]
                dydv = fdydv[j]
                
                xwcoord = xwcoorda[i, j]
                ywcoord = ywcoorda[i, j]
                dxdw = dxdwa[i, j]
                dydw = dydwa[i, j]

                dwdu = self.dwdu_values[i][j]
                dwdv = self.dwdv_values[i][j]

                # Cross-product checks
                # Check 1
                lhs = dwdu * ((xucoord - xvcoord) * dydw - (yucoord - yvcoord) * dxdw)
                rhs = (xwcoord - xvcoord) * dydu - (ywcoord - yvcoord) * dxdu
                t = abs(lhs - rhs) * (self.umax - self.umin)
                if t > maxerr: maxerr = t

                # Check 2
                lhs = dwdv * ((xucoord - xvcoord) * dydw - (yucoord - yvcoord) * dxdw)
                rhs = (xucoord - xwcoord) * dydv - (yucoord - ywcoord) * dxdv
                t = abs(lhs - rhs) * (self.vmax - self.vmin)
                if t > maxerr: maxerr = t

                # Check 3
                lhs = dwdv * ((xwcoord - xvcoord) * dydu - (ywcoord - yvcoord) * dxdu)
                rhs = dwdu * ((xucoord - xwcoord) * dydv - (yucoord - ywcoord) * dxdv)
                t = abs(lhs - rhs) * (self.umax - self.umin) * (self.vmax - self.vmin) / (self.wmax - self.wmin)
                if t > maxerr: maxerr = t

        print("max derivative error is {:.2g}".format(maxerr))

    def _verify_alignment_dense(self):
        """Checks alignment error on a dense grid (every 1mm)."""
        d = 1  # check every 1 mm
        print("checking solution every ", d, "mm, ")
        ds = d / math.sqrt(self.height * self.width)
        
        # Prepare interpolators for derivatives
        fdxdu = self.diffmatu @ self.xu
        fdydu = self.diffmatu @ self.yu
        fdxdv = self.diffmatv @ self.xv
        fdydv = self.diffmatv @ self.yv

        baryfdxdu = scipy.interpolate.BarycentricInterpolator(self.unodes, fdxdu)
        baryfdydu = scipy.interpolate.BarycentricInterpolator(self.unodes, fdydu)
        baryfdxdv = scipy.interpolate.BarycentricInterpolator(self.vnodes, fdxdv)
        baryfdydv = scipy.interpolate.BarycentricInterpolator(self.vnodes, fdydv)

        # Generate dense U points based on arc length
        upoints = []
        u = self.umin
        while u <= self.umax:
            upoints.append(u)
            dx = baryfdxdu(u)
            dy = baryfdydu(u)
            dd = ds / math.hypot(dx, dy)
            if u + dd > self.umax and u < self.umax:
                u = self.umax
            else:
                u = u + dd

        # Generate dense V points
        vpoints = []
        v = self.vmin
        while v <= self.vmax:
            vpoints.append(v)
            dx = baryfdxdv(v)
            dy = baryfdydv(v)
            dd = ds / math.hypot(dx, dy)
            if v + dd > self.vmax and v < self.vmax:
                v = self.vmax
            else:
                v = v + dd

        # Batch Interpolate coordinates for these dense points
        upoints = np.array(upoints)
        vpoints = np.array(vpoints)
        
        xucoorda = _barycentric_interp_grid(upoints, self.unodes, self.u_weights, self.xu)
        yucoorda = _barycentric_interp_grid(upoints, self.unodes, self.u_weights, self.yu)
        xvcoorda = _barycentric_interp_grid(vpoints, self.vnodes, self.v_weights, self.xv)
        yvcoorda = _barycentric_interp_grid(vpoints, self.vnodes, self.v_weights, self.yv)

        maxdiff = 0
        
        def evaluate(val, nodes, weights, y_vals):
            return _barycentric_interp_scalar(val, nodes, weights, y_vals)

        # Iterate over all pairs
        for u, xucoord, yucoord in zip(upoints, xucoorda, yucoorda):
            for v, xvcoord, yvcoord in zip(vpoints, xvcoorda, yvcoorda):
                
                wvalue = self.w_func(u, v) 
                
                xwcoord = evaluate(wvalue, self.wnodes, self.w_weights, self.xw)
                ywcoord = evaluate(wvalue, self.wnodes, self.w_weights, self.yw)
                
                # Check linearity (normalized)
                difference = abs((xucoord - xvcoord) * (yucoord - ywcoord) - (xucoord - xwcoord) * (yucoord - yvcoord)) / \
                    math.hypot(xucoord - xvcoord, xvcoord - yvcoord)

                # Report range errors
                w_curr_user = self.w_user_trans(wvalue)
                wmin_user = self.w_user_trans(self.wmin)
                wmax_user = self.w_user_trans(self.wmax)
                
                if wvalue < self.wmin and not math.isclose(wvalue, self.wmin, rel_tol=0.01):
                    print( "scale range error, please check w scale min limits", \
                           "\nw({:g}, {:g}) = {:g} < wmin, {:g}".format(
                               self.u_user_trans(u), self.v_user_trans(v), w_curr_user, wmin_user))
                
                elif wvalue > self.wmax and not math.isclose(wvalue, self.wmax, rel_tol=0.01):
                    print("scale range error, please check w scale max limits", \
                          "\nw({:g}, {:g}) = {:g} > wmax, {:g}".format(
                              self.u_user_trans(u), self.v_user_trans(v), w_curr_user, wmax_user))

                if difference > maxdiff:
                    maxdiff = difference

        aler = maxdiff * math.sqrt(self.width * self.height)
        print("alignment error is estimated at less than {:5.2g} mm".format(aler))
        if aler > 0.2:
            print("alignment errors are possible - please check.")
            print(f"This nomogram used a polynomial defined with {self.NN} points ")
            print("Try increasing this, or reduce the range of one or more scales")

    def _inject_results(self):
        """Creates callables for the final geometry and attaches them to params."""
        def make_eval_func(nodes, weights, y_coeffs, plot_trans):
            return lambda val: _barycentric_interp_scalar(plot_trans(val), nodes, weights, y_coeffs)

        self.params_u.update({
            'f': make_eval_func(self.unodes, self.u_weights, self.xu, self.u_plot_trans),
            'g': make_eval_func(self.unodes, self.u_weights, self.yu, self.u_plot_trans),
            'h': lambda u: 1.0
        })
        self.params_v.update({
            'f': make_eval_func(self.vnodes, self.v_weights, self.xv, self.v_plot_trans),
            'g': make_eval_func(self.vnodes, self.v_weights, self.yv, self.v_plot_trans),
            'h': lambda v: 1.0
        })
        self.params_w.update({
            'f': make_eval_func(self.wnodes, self.w_weights, self.xw, self.w_plot_trans),
            'g': make_eval_func(self.wnodes, self.w_weights, self.yw, self.w_plot_trans),
            'h': lambda w: 1.0
        })

    def _configure_tick_placement(self):
        """Heuristic logic to place ticks on the correct side of the curve."""
        if self.yw[0] > self.yw[-1]:
            wEast, wWest = 'left', 'right'
        else:
            wEast, wWest = 'right', 'left'
        
        vEast = 'left' if self.yv[0] > self.yv[-1] else 'right'

        distuw = self.xw - self.xu
        distvw = self.xv - self.xw

        if 'tick_side' not in self.params_u and np.min(distuw) < 0.2:
            self.params_u['tick_side'] = 'left'

        if 'tick_side' not in self.params_v and np.min(distvw) < 0.2:
            self.params_v['tick_side'] = vEast
            self.params_v['turn_relative'] = True

        if 'tick_side' not in self.params_w:
            self.params_w['tick_side'] = wWest if np.linalg.norm(distuw) > np.linalg.norm(distvw) else wEast
            self.params_w['turn_relative'] = True
            print(f'Putting w scale ticks on {self.params_w["tick_side"]} side')

    def _configure_dual_scales(self):
        """Links secondary scales (Type 8 or Type 9) to the main generated scales."""
        for params in [self.params_u, self.params_v, self.params_w]:
            if 'tag' not in params: continue
            ltag = params['tag']
            
            # Search other blocks for matching tags
            for b in self.main_params['block_params'][1:]:
                if b['block_type'] == 'type_8':
                    laxis = b.get('f_params', {})
                    if laxis.get('tag') == ltag and 'function_x' not in laxis:
                        if 'align_func' not in laxis: sys.exit(f"dual axis '{ltag}' needs 'align_func'")
                        fal = laxis['align_func']
                        laxis['function_x'] = functools.partial(lambda u, p, f: p['f'](f(u)), p=params, f=fal)
                        laxis['function_y'] = functools.partial(lambda u, p, f: p['g'](f(u)), p=params, f=fal)
                        
                elif b['block_type'] == 'type_9':
                    for fp in ['f1_params', 'f2_params', 'f3_params']:
                        if fp in b:
                            laxis = b[fp]
                            if laxis.get('tag') == ltag and 'f' not in laxis:
                                if 'align_func' not in laxis: sys.exit(f"dual axis '{ltag}' needs 'align_func'")
                                fal = laxis['align_func']
                                laxis['f'] = functools.partial(lambda u, p, f: p['f'](f(u)), p=params, f=fal)
                                laxis['g'] = functools.partial(lambda u, p, f: p['g'](f(u)), p=params, f=fal)
                                laxis['h'] = params['h']