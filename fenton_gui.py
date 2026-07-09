# ==============================================================================
#  ENGINEERING TECHNICAL REFERENCE & THEORETICAL FORMULATION
# ==============================================================================
#  PROGRAM:      Fenton Nonlinear Wave Suite - Python Tkinter GUI solver
#  FILE:         fenton_gui.py
#  METHOD:       Stream-function / Fourier collocation method for steady,
#                two-dimensional, finite-amplitude gravity waves.
#  REFERENCE:    Rienecker-Fenton / Fenton numerical-wave formulation,
#                as implemented in this production solver.
# ==============================================================================
#
#  1. SCOPE OF THIS IMPLEMENTATION
#  -----------------------------------------------------------------------------
#  This program computes steady, periodic, finite-amplitude gravity waves in
#  finite water depth, with optional collinear current.  The calculation follows
#  the stream-function / Fourier reference formulation: the Fourier basis satisfies Laplace's
#  equation and the impermeable-bed condition analytically, while Newton iteration
#  solves the nonlinear free-surface streamline condition, Bernoulli condition,
#  wave-height constraint, period relation, current closure and global flux
#  identities.
#
#  The central engineering mapping is
#
#      L = F(H, T, d, U_c)
#
#  where H is crest-to-trough wave height, T is the apparent fixed-frame period,
#  d is still-water depth, U_c is the imposed collinear current under the selected
#  current convention, and L is the dynamically consistent wavelength.  In the
#  Fenton formulation a steady wave train is fundamentally defined by H, d and L;
#  when T is specified instead of L, the current or wave-speed convention is part
#  of the closure because the observed period is Doppler-shifted by current.
#
#  2. COORDINATES, NONDIMENSIONAL VARIABLES AND VELOCITIES
#  -----------------------------------------------------------------------------
#  The fixed-frame coordinate uses x in the direction of propagation, y = 0 at
#  still-water level and y = -d at the bed.  The moving wave-frame phase and the
#  nondimensional vertical coordinate are
#
#      X = k(x - c t),        Y = k y,        k = 2*pi/L.
#
#  One wavelength corresponds to X in [0, 2*pi], the bed is Y = -kd, and the
#  free surface is represented by
#
#      Y = zeta(X) = k eta(X).
#
#  The nondimensional stream function is
#
#      Psi = psi * sqrt(k^3/g),
#
#  and the wave-frame velocities are
#
#      U_hat = U * sqrt(k/g) = dPsi/dY,
#      V_hat = V * sqrt(k/g) = -dPsi/dX.
#
#  3. STREAM-FUNCTION / FOURIER REPRESENTATION
#  -----------------------------------------------------------------------------
#  The finite-depth Fourier representation used by the residual equations is
#
#      Psi(X,Y) = -Ubar*(Y + kd)
#                 + sum_{j=1..N} B_j sinh[j(Y+kd)]/cosh(jkd) cos(jX),
#
#  where Ubar* = Ubar sqrt(k/g).  At the bed, Y = -kd, both the linear term and
#  all hyperbolic sine terms vanish, so the bed condition is satisfied before the
#  nonlinear equations are assembled.
#
#  For numerical stability the hyperbolic ratios are evaluated in the stable hyperbolic-ratio
#  form
#
#      S_j(Y) = sinh(jY) + cosh(jY) tanh(jkd),
#      C_j(Y) = cosh(jY) + sinh(jY) tanh(jkd).
#
#  For large jkd, tanh(jkd) tends to 1 and both functions tend to exp(jY).
#
#  4. STATE VECTOR z[i]
#  -----------------------------------------------------------------------------
#  This file implements the Fenton-style 1-based state vector.  Index 0
#  is intentionally unused so that the residual equations map directly to the
#  z[i] notation:
#
#      z[1]  = kd
#      z[2]  = kH
#      z[3]  = T sqrt(gk)
#      z[4]  = c sqrt(k/g)
#      z[5]  = Eulerian current variable ubar_1 sqrt(k/g)
#      z[6]  = Stokes / mass-transport current variable ubar_2 sqrt(k/g)
#      z[7]  = mean wave-frame velocity Ubar sqrt(k/g)
#      z[8]  = q sqrt(k^3/g), with q = Ubar d - Q
#      z[9]  = r k/g, with r = R - g d
#      z[10] ... z[N+10]       = free-surface ordinates zeta_m = k eta_m
#      z[N+11] ... z[2N+10]    = Fourier coefficients B_j
#
#  Therefore, for Fourier order N, the active vector length is
#
#      num = 2N + 10.
#
#  With the production value N = 50, this gives 110 active 1-based state entries
#  and 110 nonlinear residual equations, matching the full nonlinear residual
#  system for the stream-function / Fourier collocation solver.
#
#  5. RESIDUAL SYSTEM IMPLEMENTED IN _eqns()
#  -----------------------------------------------------------------------------
#  At convergence every component of F(z) is zero.  The first eight residuals are
#  the global scalar constraints; the next N+1 enforce the free-surface streamline
#  condition; the final N+1 enforce Bernoulli's equation at the same collocation
#  nodes X_m = m*pi/N, m = 0..N.
#
#      r1 = z2 - z1(H/d)
#      r2 = z2 - Hs z3^2
#      r3 = z4 z3 - 2*pi
#      r4 = z5 + z7 - z4
#      r5 = z1(z6 + z7 - z4) - z8
#      r6 = z[c+4] - U_c sqrt(z1)
#      r7 = z10 + z[N+10] + 2 sum_{i=1..N-1} z[10+i]
#      r8 = z10 - z[N+10] - z2
#
#  For each free-surface node m = 0..N,
#
#      r[9+m] = psi_m - z8 - z7 z[10+m]
#
#  and
#
#      r[N+10+m] = 0.5*((-z7 + u_m)^2 + v_m^2) + z[10+m] - z9.
#
#  These equations solve the wave form, wave speed, current variables, mean-flow
#  flux, Bernoulli offset and Fourier coefficients as one coupled nonlinear
#  free-boundary problem.  They are not empirical fitting equations and are not a
#  post-processing correction to linear Airy theory.
#
#  6. CURRENT CONVENTIONS
#  -----------------------------------------------------------------------------
#  The solver keeps the current variables separate, as required by the Fenton current convention:
#
#      ubar_1 = c - Ubar,
#      ubar_2 = c - Q/d,
#      q      = Ubar d - Q.
#
#  The GUI uses the Eulerian-current criterion by default, but the state vector and
#  residual equation r6 preserve the Fenton current-selector convention.
#
#  7. NUMERICAL SOLUTION STRATEGY
#  -----------------------------------------------------------------------------
#  Wave height is introduced by continuation from a near-linear wave to the target
#  value.  At each continuation step, Newton iteration linearizes F(z) about the
#  current iterate and solves J dz = -F(z).  Numerical finite-difference Jacobians
#  and singular-value stabilized dense linear algebra are used because high-order
#  stream-function waves, near-limiting waves and current cases can make the
#  Jacobian ill-conditioned.
#
#  Linear Airy/Fenton-McKee estimates are used only as initialization and checking
#  aids.  The final result is the converged nonlinear stream-function / Fourier
#  collocation solution.
#
#  8. OUTPUT INTERPRETATION
#  -----------------------------------------------------------------------------
#  The report follows the distinction between k-based and d-based
#  nondimensional quantities: kd, kH, T sqrt(gk), c sqrt(k/g), current variables,
#  fluxes, Bernoulli constants, impulse, energy, momentum flux, radiation stress
#  and wave power are reported alongside their engineering depth-scaled forms.
#
#  The same calculation is run for no-current and, when requested, with-current
#  cases so the output can be compared under the selected current convention.
#
#  9. RUNNING FROM SOURCE
#  -----------------------------------------------------------------------------
#  Recommended development environment:
#
#      python -m venv build_env
#      build_env\Scripts\activate        (Windows)
#      python -m pip install --upgrade pip
#      python -m pip install -U "numpy>=1.20"
#      python -m pip install -U "numba>=0.57"   (optional acceleration)
#      python fenton_gui.py
#
#  If Numba is unavailable, the solver automatically uses the pure-NumPy path with
#  the same residual equations and state-vector convention.
#
#  10. COMPILING TO A STANDALONE EXECUTABLE
#  -----------------------------------------------------------------------------
#  To create a single-file Windows GUI executable:
#
#      python -m pip install pyinstaller
#      pyinstaller --onefile --noconsole fenton_gui.py
#
#  The generated executable is written to the dist folder.
# ==============================================================================

import os
import sys
import io
import tempfile
import warnings
import threading
import textwrap
import tkinter as tk
from tkinter import ttk, scrolledtext, messagebox, Menu

# ------------------------------------------------------------------------------
#  RUNTIME STABILITY (Windows / BLAS / JIT)
# ------------------------------------------------------------------------------
#  On some Windows 10 configurations, NumPy's underlying BLAS/LAPACK libraries
#  may oversubscribe CPU threads during repeated SVD solves, which can look like
#  a "hang" in GUI applications (high CPU usage with no visible progress).
#  The environment variables below cap common BLAS thread pools to 1 thread.
#  Users can override these by defining the variables before launching Python.
# ------------------------------------------------------------------------------
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("NUMEXPR_NUM_THREADS", "1")

# Ensure Numba cache is writable on Windows (avoids slowdowns/permission issues).
if os.name == "nt":
    try:
        _numba_cache_dir = os.path.join(tempfile.gettempdir(), "fenton_numba_cache")
        os.makedirs(_numba_cache_dir, exist_ok=True)
        os.environ.setdefault("NUMBA_CACHE_DIR", _numba_cache_dir)
    except Exception:
        pass

# GUI mode (pythonw / PyInstaller --noconsole) provides no stdout/stderr streams.
# We capture diagnostics in-memory and flush them into the GUI upon failure.
_ORIG_STDOUT = sys.stdout
_ORIG_STDERR = sys.stderr
_STDIO_BUFFER = io.StringIO()

if (sys.stdout is None) or (sys.stderr is None) or getattr(sys, "frozen", False):
    sys.stdout = _STDIO_BUFFER
    sys.stderr = _STDIO_BUFFER

import numpy as np

# NumPy compatibility:
# - NumPy 2.x removed the legacy alias `numpy.trapz`.
# - Prefer `numpy.trapezoid` when available, else fall back to `numpy.trapz`.
def _np_trapz(y, x):
    trapezoid = getattr(np, "trapezoid", None)
    if trapezoid is not None:
        return trapezoid(y, x)
    trapz = getattr(np, "trapz", None)
    if trapz is not None:
        return trapz(y, x)
    # Minimal fallback (should never trigger on normal NumPy builds)
    y = np.asarray(y)
    x = np.asarray(x)
    if y.size < 2:
        return 0.0
    dx = np.diff(x)
    return np.sum((y[:-1] + y[1:]) * 0.5 * dx)

# ==============================================================================
#  NUMBA ACCELERATION LAYER (OPTIONAL, SAME RESIDUAL SYSTEM)
# ==============================================================================
#
# The solver is dominated by repeated evaluations of the documented residual system
# _eqns() inside a finite-difference Jacobian (≈ O(num^2) residual calls per Newton
# step). Pure NumPy vectorisation creates many temporary arrays per collocation
# node; Numba JIT removes those temporaries by compiling tight scalar loops.
#
# IMPORTANT:
# - fastmath is NOT used (to avoid changing floating-point semantics).
# - These kernels implement the same residual equations as the NumPy path.
# - If Numba is unavailable, the code automatically falls back to the original
#   NumPy implementation with identical behaviour.
#
try:
    from numba import njit
    NUMBA_AVAILABLE = True
except Exception:  # pragma: no cover
    NUMBA_AVAILABLE = False

    def njit(*args, **kwargs):
        def _wrap(fn):
            return fn
        return _wrap


# ==============================================================================
#  GLOBAL CONSTANTS AND SOLVER CONFIGURATION
# ==============================================================================

# Physical constants used for dimensional conversion and reporting
G_STD = 9.80665         # Standard Gravity [m/s^2]
RHO   = 1025.0          # Density of Seawater [kg/m^3]

# Numerical configuration for the Fenton Fourier collocation system
DTYPE = np.float64      # Precision for floating point arithmetic
N_FOURIER = 50          # Fourier order N for the N+1 half-wave collocation nodes. 
                        # The active documented state-vector length is 2N+10 for this solver.
N_NUMBERS = 8           # formatting precision for output text

# Suppress benign numerical warnings from early Newton/continuation iterations
# where the nonlinear residual system can be temporarily ill-conditioned.
warnings.filterwarnings('ignore')

@njit(cache=True)
def _eqns_numba(z, rhs, coeff, Tanh, cos_nm, sin_nm, n, num, Hoverd, height, Current, Current_criterion):
    """Numba implementation of the finite-depth period residuals F(z)."""
    pi = np.pi

    # r1: depth-height relation, z2 = z1(H/d)
    rhs[1] = z[2] - z[1] * Hoverd

    # r2: period-input height relation, z2 = Hs*z3^2
    rhs[2] = z[2] - height * z[3] * z[3]

    # r3: apparent-period/celerity relation, z4*z3 = 2*pi
    rhs[3] = z[4] * z[3] - 2.0 * pi

    # r4: Eulerian current relation, ubar_1 = c - Ubar
    rhs[4] = z[5] + z[7] - z[4]

    # r5: Stokes-current/flux relation, q = Ubar*d - Q
    rhs[5] = z[1] * (z[6] + z[7] - z[4]) - z[8]

    # Fourier coefficients B_j and tanh(jkd) table for S_j and C_j
    kd = z[1]
    for i in range(1, n + 1):
        coeff[i] = z[n + i + 10]
        Tanh[i] = np.tanh(i * kd)

    # r6: selected current-criterion closure, z[c+4] = Uc*sqrt(kd)
    rhs[6] = z[Current_criterion + 4] - Current * np.sqrt(kd)

    # r7: mean-level convention over the half-wave collocation ordinates
    rhs7 = z[10] + z[n + 10]
    for i in range(1, n):
        rhs7 += 2.0 * z[10 + i]
    rhs[7] = rhs7

    # r8: crest-to-trough height definition, z10 - z[N+10] = z2
    rhs[8] = z[10] - z[n + 10] - z[2]

    # r[9+m] and r[N+10+m]: streamline and Bernoulli residuals at m=0..N
    for m in range(0, n + 1):
        zsurf = z[10 + m]  # zeta_m = k*eta_m at this free-surface node

        psi = 0.0
        u = 0.0
        v = 0.0

        for jj in range(1, n + 1):
            cj = coeff[jj]
            tj = Tanh[jj]

            x = jj * zsurf
            # Numerical safeguard for diverging trial free-surface ordinates
            if x > 60.0 or x < -60.0:
                rhs[1] = np.inf
                return np.inf
            e = np.exp(x)
            inv_e = 1.0 / e
            sinhkd = 0.5 * (e - inv_e)
            coshkd = 0.5 * (e + inv_e)

            # Stable hyperbolic rewrite: S_j and C_j using tanh(jkd)
            S = sinhkd + coshkd * tj
            C = coshkd + sinhkd * tj

            ccos = cos_nm[m, jj - 1]
            ssin = sin_nm[m, jj - 1]

            psi += cj * S * ccos
            jcj = jj * cj
            u += jcj * C * ccos
            v += jcj * S * ssin

        rhs[m + 9] = psi - z[8] - z[7] * z[m + 10]
        rhs[n + m + 10] = 0.5 * ((-z[7] + u) ** 2 + v * v) + z[m + 10] - z[9]

    # Objective value ||F(z)||^2 used only to monitor Newton progress
    ss = 0.0
    for i in range(1, num + 1):
        ss += rhs[i] * rhs[i]
    return ss


@njit(cache=True)
def _compute_Y_and_B_numba(z, B, Y, cosa, n):
    """Numba implementation of the post-convergence half-wave Fourier transform."""
    for i in range(0, Y.shape[0]):
        Y[i] = 0.0

    two_n = 2 * n

    for j in range(1, n + 1):
        B[j] = z[j + n + 10]

        # Half-wave cosine transform contribution from crest and trough ordinates
        sign = 1.0
        if (j % 2) == 1:
            sign = -1.0
        s = 0.5 * (z[10] + z[n + 10] * sign)

        for m in range(1, n):
            s += z[10 + m] * cosa[(m * j) % two_n]

        Y[j] = 2.0 * s / n


@njit(cache=True)
def _surface_keta_numba(Y, n, X):
    """Numba implementation of the surface reconstruction zeta(X)=k*eta(X)."""
    kEta = 0.0
    for j in range(1, n):
        kEta += Y[j] * np.cos(j * X)
    kEta += 0.5 * Y[n] * np.cos(n * X)
    return kEta


@njit(cache=True)
def _point_numba(X, Y, z, Tanh, B, n):
    """Numba implementation of the finite-depth velocity kernel."""
    kd = z[1]

    # Depth-scaled celerity/current variables derived from the k-based state vector
    c  = z[4] / np.sqrt(kd)
    ce = z[5] / np.sqrt(kd)
    R  = 1.0 + z[9] / kd

    # Local Fourier velocity sums in wave-scaled coordinates (X,Y)
    u = 0.0
    v = 0.0
    ux = 0.0
    vx = 0.0

    for j in range(1, n + 1):
        Cos = np.cos(j * X)
        Sin = np.sin(j * X)

        coshdelta = np.cosh(j * Y)
        sinhdelta = np.sinh(j * Y)
        C = coshdelta + sinhdelta * Tanh[j]
        S = sinhdelta + coshdelta * Tanh[j]

        Bj = B[j]
        u  += j * Bj * C * Cos
        v  += j * Bj * S * Sin
        ux += - (j * j) * Bj * C * Sin
        vx += (j * j) * Bj * S * Cos

    # Convert k-based derivatives to the d-based output scaling used in output tables
    inv_kd_sqrt = 1.0 / np.sqrt(kd)
    inv_kd_32   = 1.0 / (kd ** 1.5)

    u  *= inv_kd_sqrt
    v  *= inv_kd_sqrt
    ux *= np.sqrt(kd)
    vx *= np.sqrt(kd)

    # Add the Eulerian current component to the fixed-frame horizontal velocity
    u = ce + u

    # Time derivatives obtained from steady wave-frame fields and phase speed
    ut = -c * ux
    vt = -c * vx
    uy = vx
    vy = -ux

    dudt = ut + u * ux + v * uy

    return float(u), float(v), float(dudt)


# ==============================================================================
#  CORE SOLVER CLASS (FENTON z-VECTOR / FOURIER COLLOCATION)
# ==============================================================================
#
# This class implements the finite-depth, period-input stream-function /
# Fourier collocation solver.  The algebra follows the documented Fenton-style
# state vector and residual system, including:
#   - Unknown vector z[1..2N+10] with index 0 intentionally unused
#   - Residual system F(z)=0 with eight global equations and two collocation blocks
#   - Newton iteration with finite-difference Jacobian and SVD-stabilized solve
#   - Wave-height continuation from a near-linear seed to the target height
#
# Notes
# -----
# - The GUI (Tkinter) is a presentation layer; the solver core is the stream-function / Fourier collocation calculation.
# - This implementation keeps 1-based indexing for z[] so the equations map
#   directly to the documented residual notation.
#

class FentonStreamFunction:
    """
    Fenton steady-wave solver using the stream-function / Fourier
    collocation formulation with the documented 1-based z[i] state vector.

    Public API is preserved to keep the GUI unchanged.
    """

    # --------------------------- construction and state-vector allocation ---------

    def __init__(self, H, T, d, Uc=0.0):
        # Dimensional inputs: H, T, d and imposed Eulerian current Uc
        self.H_target = float(H)    # [m]
        self.T_target = float(T)    # [s]
        self.d        = float(d)    # [m]
        self.Uc       = float(Uc)   # [m/s] (Eulerian / lab-frame)

        # Physical and numerical constants
        self.g = G_STD
        self.N = N_FOURIER

        # Solver control: continuation and Newton iteration tolerances
        self.nstep  = 4        # continuation steps in wave height
        self.number = 40       # max Newton iterations per step
        self.crit   = 1.0e-8   # intermediate-step convergence factor (C++: crit)
        self.criter_final = 1.0e-10  # final-step convergence factor

        # Problem mode represented by this GUI: finite depth with period as input
        self.Depth = "Finite"
        self.Case  = "Period"

        # Current selector c: 1=Eulerian current, 2=Stokes/mass-transport current
        self.Current_criterion = 1

        # Derived nondimensional inputs used by the documented residual equations
        # MaxH=H/d, T_nd=T*sqrt(g/d), Height=H/(gT^2)
        self.MaxH = self.H_target / self.d if self.d > 0 else 0.0
        self.T_nd = self.T_target * np.sqrt(self.g / self.d) if self.d > 0 else 0.0
        self.Height = (self.MaxH / (self.T_nd * self.T_nd)) if self.T_nd > 0 else 0.0

        # Current input in depth scaling; r6 converts it to k-based z[c+4]
        self.Current = self.Uc / np.sqrt(self.g * self.d) if self.d > 0 else 0.0

        # ------------------------- dimensional outputs and diagnostics ----------------
        self.k = 0.0
        self.L = 0.0
        self.c = 0.0
        self.converged = False

        # Human-readable failure reason for GUI reporting when F(z) does not converge
        self.last_error = ""

        # Robustness for large ambient currents: allow more continuation/Newton attempts
        # without changing the documented residual equations.
        #
        if abs(self.Current) >= 1.0:
            self.nstep = max(self.nstep, 8)
            self.number = max(self.number, 80)
        self.Bj = np.zeros(self.N, dtype=DTYPE)             # B_1..B_N (0-based in Python)
        self.eta_nodes = np.zeros(self.N + 1, dtype=DTYPE)  # absolute z from bed [m]

        self.eta_crest = 0.0
        self.eta_trough = 0.0
        self.steepness = 0.0
        self.rel_depth = 0.0
        self.ursell = 0.0
        self.regime = ""

        self.breaking_index = 0.0
        self.is_breaking = False
        self.breaking_limit_miche = 0.0

        self.u_bed = 0.0
        self.tau_bed = 0.0
        self.acc_max = 0.0
        self.u_surf = 0.0
        self.w_max = 0.0
        self.asymmetry = 0.0
        self.ExcursionBed = 0.0

        self.Cg = 0.0
        self.Power = 0.0
        self.EnergyDensity = 0.0
        self.Sxx = 0.0
        self.Impulse = 0.0
        self.MassTransport = 0.0
        self.BernoulliR = 0.0

        # ------------------------- internal documented state arrays -----------------------
        self.n = int(self.N)
        self.num = 2 * self.n + 10

        # 1-based vectors; index 0 unused to preserve z[i] notation
        self.z    = np.zeros(self.num + 1, dtype=DTYPE)
        self.rhs1 = np.zeros(self.num + 1, dtype=DTYPE)
        self.rhs2 = np.zeros(self.num + 1, dtype=DTYPE)
        self.coeff = np.zeros(self.n + 1, dtype=DTYPE)    # coeff[1..n]
        self.Tanh  = np.zeros(self.n + 1, dtype=DTYPE)    # Tanh[1..n]
        self.B     = np.zeros(self.n + 1, dtype=DTYPE)    # B[1..n]
        self.Y     = np.zeros(self.num + 1, dtype=DTYPE)  # Y[0..n] used; keep size

        # Precomputed half-wave trigonometric tables for X_m = m*pi/N
        self.cosa = np.zeros(2 * self.n + 1, dtype=DTYPE)  # [0..2n]
        self.sina = np.zeros(2 * self.n + 1, dtype=DTYPE)

        # Precompute collocation lookup tables cos(jX_m) and sin(jX_m)
        k_idx = np.arange(0, 2 * self.n + 1, dtype=DTYPE)
        self.cosa[:] = np.cos(k_idx * np.pi / self.n)
        self.sina[:] = np.sin(k_idx * np.pi / self.n)

        self._j = np.arange(1, self.n + 1, dtype=DTYPE)
        self._j_int = np.arange(1, self.n + 1, dtype=np.int64)
        self._nm_map = (np.arange(0, self.n + 1, dtype=np.int64)[:, None] * self._j_int[None, :]) % (2 * self.n)
        self._cos_nm = self.cosa[self._nm_map]  # shape (n+1, n)
        self._sin_nm = self.sina[self._nm_map]

        # Continuation extrapolation storage for consecutive z[i] states
        self.sol = np.zeros((self.num + 1, 3), dtype=DTYPE)

        # Continuation-step variables Hs and H/d used in r1 and r2
        self.height = 0.0    # stepped 'height' (dimensionless)
        self.Hoverd = 0.0    # stepped H/d

    # --------------------------- solver helper methods ---------------------------

    def _init_linear(self):
        """
        Initialization for the finite-depth, period-input Fenton formulation.
        Produces the first continuation seed for the documented z[1..num] state.
        """
        n = self.n
        pi = np.pi

        # Finite-depth period-input initialization
        sigma = 2.0 * pi * np.sqrt(self.height / self.Hoverd) if self.Hoverd > 0 else 0.0

        # Linear/Fenton-McKee estimate used only to seed the nonlinear nonlinear system
        if sigma > 0:
            self.z[1] = (sigma * sigma) / (np.tanh(sigma ** 1.5) ** (2.0 / 3.0))
        else:
            # Very small-wave fallback seed; final solution still comes from F(z)=0
            self.z[1] = 2.0 * pi * max(self.height, 1e-12) / max(self.Hoverd, 1e-12)

        self.z[2] = self.z[1] * self.Hoverd
        self.z[4] = np.sqrt(np.tanh(self.z[1]))
        self.z[3] = 2.0 * pi / self.z[4]

        # Initial current variables for the selected Fenton current criterion
        if self.Current_criterion == 1:
            self.z[5] = self.Current * np.sqrt(self.z[2])
            self.z[6] = 0.0
        else:
            self.z[6] = self.Current * np.sqrt(self.z[2])
            self.z[5] = 0.0

        self.z[7] = self.z[4]
        self.z[8] = 0.0
        self.z[9] = 0.5 * self.z[7] * self.z[7]

        # Initial free-surface ordinates z[10..N+10] and Fourier coefficients B_j
        self.z[10] = 0.5 * self.z[2]
        for i in range(1, n + 1):
            self.z[n + i + 10] = 0.0
            self.z[i + 10] = 0.5 * self.z[2] * self.cosa[i]

        self.z[n + 11] = 0.5 * self.z[2] / self.z[7]

        # Store first continuation seed for later z[i] extrapolation
        for i in range(1, 10):
            self.sol[i, 1] = self.z[i]
        for i in range(10, self.num + 1):
            self.sol[i, 1] = 0.0

    def _eqns(self, rhs_out):
        """
        Evaluate the nonlinear residual vector F(z). Fills rhs_out[1..num] and
        returns ||F(z)||^2 for the finite-depth, period-input GUI mode.
        """
        # Numba-accelerated evaluation of the same residual equations.
        #
        # Robustness note:
        # If the JIT path produces NaN/Inf during a rejected Newton trial,
        # fall back to pure NumPy for clearer diagnostics; the algebra is unchanged.
        if NUMBA_AVAILABLE:
            ss = _eqns_numba(self.z, rhs_out, self.coeff, self.Tanh, self._cos_nm, self._sin_nm,
                             self.n, self.num, self.Hoverd, self.height, self.Current, self.Current_criterion)
            if np.isfinite(ss) and np.isfinite(rhs_out[1:self.num + 1]).all():
                return ss
            # else: continue into the NumPy implementation

        n = self.n
        num = self.num
        pi = np.pi
        z = self.z
        rhs = rhs_out

        # r1: depth-height relation, z2 = z1(H/d)
        rhs[1] = z[2] - z[1] * self.Hoverd

        # r2: period-input height relation, z2 = Hs*z3^2
        rhs[2] = z[2] - self.height * z[3] * z[3]

        # r3: apparent-period/celerity relation, z4*z3 = 2*pi
        rhs[3] = z[4] * z[3] - 2.0 * pi

        # r4: Eulerian current relation, ubar_1 = c - Ubar
        rhs[4] = z[5] + z[7] - z[4]

        # r5: Stokes-current/flux relation, q = Ubar*d - Q
        rhs[5] = z[1] * (z[6] + z[7] - z[4]) - z[8]

        # Fourier coefficients B_j and tanh(jkd) table for S_j and C_j
        for i in range(1, n + 1):
            self.coeff[i] = z[n + i + 10]
            self.Tanh[i] = np.tanh(i * z[1])

        # r6: selected current-criterion closure, z[c+4] = Uc*sqrt(kd)
        rhs[6] = z[self.Current_criterion + 4] - self.Current * np.sqrt(z[1])

        # r7: mean-level convention over the half-wave collocation ordinates
        rhs[7] = z[10] + z[n + 10]
        for i in range(1, n):
            rhs[7] += 2.0 * z[10 + i]

        # r8: crest-to-trough height definition, z10 - z[N+10] = z2
        rhs[8] = z[10] - z[n + 10] - z[2]

        # r[9+m] and r[N+10+m]: streamline and Bernoulli residuals at m=0..N
        j = self._j                      # shape (n,)
        coeff = self.coeff[1:n + 1]      # shape (n,)
        tanh = self.Tanh[1:n + 1]        # shape (n,)
        jcoeff = j * coeff               # shape (n,)

        for m in range(0, n + 1):
            zsurf = z[10 + m]  # zeta_m = k*eta_m at this free-surface node

            x = j * zsurf
            if np.any(x > 60.0) or np.any(x < -60.0):
                raise FloatingPointError("Divergence: exp(j*zsurf) out of safe range.")
            e = np.exp(x)
            inv_e = 1.0 / e
            sinhkd = 0.5 * (e - inv_e)
            coshkd = 0.5 * (e + inv_e)

            # Stable hyperbolic rewrite: S_j and C_j using tanh(jkd)
            S = sinhkd + coshkd * tanh
            C = coshkd + sinhkd * tanh

            cosnm = self._cos_nm[m]
            sinnm = self._sin_nm[m]

            psi = float(np.sum(coeff * S * cosnm))
            u   = float(np.sum(jcoeff * C * cosnm))
            v   = float(np.sum(jcoeff * S * sinnm))

            rhs[m + 9] = psi - z[8] - z[7] * z[m + 10]
            rhs[n + m + 10] = 0.5 * ((-z[7] + u) ** 2 + v * v) + z[m + 10] - z[9]

        return float(np.dot(rhs[1:num + 1], rhs[1:num + 1]))
    @staticmethod

    def _svd_solve(A, b):
        """
        Solve A x = b via SVD with Press et al. truncation:
          wmin = wmax * 1e-12

        Numerical safety:
        - Reject NaNs/Infs before entering LAPACK.
        - Fallback to least-squares if SVD fails to converge.
        """
        if (not np.isfinite(A).all()) or (not np.isfinite(b).all()):
            raise FloatingPointError("Non-finite values in Jacobian system (A or b).")

        try:
            U, s, Vt = np.linalg.svd(A, full_matrices=False)
            smax = np.max(s) if s.size else 0.0
            wmin = smax * 1.0e-12

            s_inv = np.zeros_like(s)
            mask = s > wmin
            s_inv[mask] = 1.0 / s[mask]

            return (Vt.T @ (s_inv * (U.T @ b)))

        except np.linalg.LinAlgError:
            # Conservative fallback: least-squares solution of the same linearised system.
            x, *_ = np.linalg.lstsq(A, b, rcond=1.0e-12)
            return x



    def _newton(self, iter_count):
        """
        Port of the C++ Newton(...) update with additional damping safeguards.

        The governing residual equations are unchanged. The only additions are:
        - Finite-difference step clamping (prevents extreme perturbations if z[i] diverges).
        - Backtracking line-search (reduces step when a full Newton step increases residuals).
        """
        n = self.n
        num = self.num

        # baseline residual
        ss0 = float(self._eqns(self.rhs1))
        if not np.isfinite(ss0):
            raise FloatingPointError("Non-finite residual norm at start of Newton step.")

        z0 = self.z.copy()

        A = np.zeros((num, num), dtype=DTYPE)
        b = np.zeros((num,), dtype=DTYPE)

        # finite-difference Jacobian (column-wise)
        for i in range(1, num + 1):
            h = 0.01 * z0[i]
            if abs(z0[i]) < 1.0e-4:
                h = 1.0e-5
            # clamp perturbation magnitude (purely numerical safeguard)
            if abs(h) > 1.0:
                h = np.copysign(1.0, h)

            self.z[i] = z0[i] + h
            self._eqns(self.rhs2)
            self.z[i] = z0[i]

            b[i - 1] = -self.rhs1[i]
            A[:, i - 1] = (self.rhs2[1:num + 1] - self.rhs1[1:num + 1]) / h

        dx = self._svd_solve(A, b)
        if not np.isfinite(dx).all():
            raise FloatingPointError("Non-finite Newton correction vector (dx).")

        # Backtracking: prefer alpha=1, reduce if it worsens residuals or violates kd>0
        alpha = 1.0
        ss_best = ss0
        z_best = z0

        while alpha >= 1.0e-4:
            z_try = z0.copy()
            z_try[1:num + 1] = z0[1:num + 1] + alpha * dx

            # Must keep kd positive and all values finite
            if (z_try[1] <= 0.0) or (not np.isfinite(z_try[1:num + 1]).all()):
                alpha *= 0.5
                continue

            self.z[:] = z_try
            ss1 = float(self._eqns(self.rhs2))

            if np.isfinite(ss1) and (ss1 <= ss_best):
                ss_best = ss1
                z_best = z_try
                # accept immediately if improvement is adequate
                if ss1 <= ss0:
                    break

            alpha *= 0.5

        # Commit best found (or revert if none acceptable)
        self.z[:] = z_best

        corr = float(np.mean(np.abs((z_best[10:n + 11] - z0[10:n + 11]))))
        return corr


    def _compute_Y_and_B(self):
        """
        Post-convergence half-wave cosine transform.
        Produces stream-function coefficients B[1..n] and surface coefficients Y[0..n].
        """
        # Numba-accelerated evaluation of the same residual equations.
        if NUMBA_AVAILABLE:
            _compute_Y_and_B_numba(self.z, self.B, self.Y, self.cosa, self.n)
            return

        n = self.n
        z = self.z
        self.Y[:] = 0.0

        for j in range(1, n + 1):
            self.B[j] = z[j + n + 10]

            s = 0.5 * (z[10] + z[n + 10] * ((-1.0) ** j))
            for m in range(1, n):
                s += z[10 + m] * self.cosa[(m * j) % (2 * n)]
            self.Y[j] = 2.0 * s / n

    def _surface_keta(self, X):
        """
        Reconstruct the free-surface ordinate zeta(X)=k*eta(X) at phase X.
        """
        # Numba-accelerated evaluation of the same residual equations.
        if NUMBA_AVAILABLE:
            return float(_surface_keta_numba(self.Y, self.n, float(X)))

        n = self.n
        kEta = 0.0
        for j in range(1, n):
            kEta += self.Y[j] * np.cos(j * X)
        kEta += 0.5 * self.Y[n] * np.cos(n * X)
        return float(kEta)

    def _point(self, X, Y):
        """
        Evaluate the finite-depth Fourier velocity field. Returns:
          u_dimless (w.r.t sqrt(g d))
          v_dimless (w.r.t sqrt(g d))
          dudt_dimless (w.r.t g)
        Input:
          X : phase in radians (k x)
          Y : vertical coordinate in wave scaling (k(z-d))
        """
        # Numba-accelerated evaluation of the same residual equations.
        if NUMBA_AVAILABLE:
            return _point_numba(float(X), float(Y), self.z, self.Tanh, self.B, self.n)

        n = self.n
        kd = float(self.z[1])

        # Depth-scaled celerity/current variables derived from the k-based state vector
        c  = float(self.z[4] / np.sqrt(kd))
        ce = float(self.z[5] / np.sqrt(kd))
        R  = float(1.0 + self.z[9] / kd)

        # Local Fourier velocity sums in wave-scaled coordinates (X,Y)
        u = 0.0
        v = 0.0
        ux = 0.0
        vx = 0.0

        for j in range(1, n + 1):
            Cos = np.cos(j * X)
            Sin = np.sin(j * X)

            coshdelta = np.cosh(j * Y)
            sinhdelta = np.sinh(j * Y)
            C = coshdelta + sinhdelta * self.Tanh[j]
            S = sinhdelta + coshdelta * self.Tanh[j]

            Bj = self.B[j]
            u  += j * Bj * C * Cos
            v  += j * Bj * S * Sin
            ux += - (j * j) * Bj * C * Sin
            vx += (j * j) * Bj * S * Cos

        # Convert k-based derivatives to the d-based output scaling used in output tables
        inv_kd_sqrt = 1.0 / np.sqrt(kd)
        inv_kd_32   = 1.0 / (kd ** 1.5)

        u  *= inv_kd_sqrt
        v  *= inv_kd_sqrt
        ux *= np.sqrt(kd)
        vx *= np.sqrt(kd)

        # Add the Eulerian current component to the fixed-frame horizontal velocity
        u = ce + u

        # Time derivatives obtained from steady wave-frame fields and phase speed
        ut = -c * ux
        vt = -c * vx
        uy = vx
        vy = -ux

        dudt = ut + u * ux + v * uy
        # Vertical acceleration is not printed but can be formed from the same derivatives:
        # dvdt = vt + u * vx + v * vy

        # Bernoulli pressure reconstruction follows the same dynamic condition:
        # y = 1.0 + Y / kd
        # Pressure = R - y - 0.5 * (((u - c) ** 2) + v * v)

        return float(u), float(v), float(dudt)

    # ----------------------------- public methods -----------------------------

    def solve(self):
        """
        Solve the steady nonlinear wave problem using Fenton's Fourier / stream-function method.

        The implementation follows the solver architecture:
        - Continuation in wave height (nstep)
        - Newton iterations on the full residual vector F(z)
        - Linear solve via SVD with Press-style truncation

        Robustness features:
        - Fail-fast on NaNs/Infs before calling LAPACK (prevents silent stalls).
        - Clear convergence state + message for GUI reporting.
        - Increased iteration budget automatically enabled for large |Uc|.
        """
        # Default outcome: failure unless we reach the end successfully
        self.converged = False
        self.last_error = ""

        # Basic input screening (physical requirements)
        if self.H_target <= 0.0 or self.T_target <= 0.0 or self.d <= 0.0:
            self.last_error = "Invalid inputs: H, T, and d must be > 0."
            return

        old_err = np.geterr()
        try:
            # Make numerical faults explicit. Underflow is benign in this context.
            np.seterr(over="raise", invalid="raise", divide="raise", under="ignore")

            # continuation step sizes
            dhe = self.Height / self.nstep
            dho = self.MaxH / self.nstep

            # height stepping
            for ns in range(1, self.nstep + 1):
                self.height = ns * dhe
                self.Hoverd = ns * dho

                # initial/extrapolated guess
                if ns == 1:
                    self._init_linear()
                else:
                    # z[i] = 2*sol[i][2] - sol[i][1]
                    self.z[1:self.num + 1] = 2.0 * self.sol[1:self.num + 1, 2] - self.sol[1:self.num + 1, 1]
                    # Fallback: if extrapolation produces an invalid start state, use the last converged state.
                    # This does not change any equations; it only prevents the Newton step from starting from NaN/Inf.
                    if (not np.isfinite(self.z[1:self.num + 1]).all()) or (self.z[1] <= 0.0):
                        self.z[1:self.num + 1] = self.sol[1:self.num + 1, 2]
                    if (not np.isfinite(self.z[1:self.num + 1]).all()) or (self.z[1] <= 0.0):
                        raise FloatingPointError("Invalid extrapolated start state for continuation step.")

                # Newton iterations
                step_converged = False
                for it in range(1, self.number + 1):
                    # Newton iteration. If the start state was poisoned by extrapolation,
                    # retry once from the last converged state (no extrapolation).
                    try:
                        err = self._newton(it)
                    except FloatingPointError:
                        if (ns > 1) and (it == 1):
                            self.z[1:self.num + 1] = self.sol[1:self.num + 1, 2]
                            err = self._newton(it)
                        else:
                            raise

                    if not np.isfinite(err):
                        raise FloatingPointError("Non-finite Newton correction.")

                    # -----------------------------------------------------------------
                    # Continuation storage must be refreshed before leaving the converged step.
                    # Otherwise sol[:,2] may hold an older z[i] state and corrupt the next
                    # continuation extrapolation.
                    # -----------------------------------------------------------------
                    if ns == 1:
                        self.sol[1:self.num + 1, 2] = self.z[1:self.num + 1]
                    else:
                        self.sol[1:self.num + 1, 1] = self.sol[1:self.num + 1, 2]
                        self.sol[1:self.num + 1, 2] = self.z[1:self.num + 1]

                    # Protect linear algebra calls on diverging states
                    if (not np.isfinite(self.z[1:self.num + 1]).all()) or (self.z[1] <= 0.0):
                        raise FloatingPointError("Divergence: non-finite/invalid state vector encountered.")

                    criter = self.criter_final if (ns == self.nstep) else self.crit
                    if (it > 1) and (err < criter * abs(self.z[1])):
                        step_converged = True
                        break


                if not step_converged:
                    self.last_error = (
                        f"Newton did not converge within {self.number} iterations "
                        f"at continuation step {ns}/{self.nstep}."
                    )
                    return

                # Update Fourier surface coefficients Y_j and stream-function coefficients B_j
                self._compute_Y_and_B()

            # ------------------------- dimensional post-processing -----------------------

            kd = float(self.z[1])
            if (not np.isfinite(kd)) or (kd <= 0.0):
                raise FloatingPointError("Invalid wavenumber (kd).")

            k_phys = kd / self.d
            L_phys = 2.0 * np.pi / k_phys
            c_phys = L_phys / self.T_target

            if (not np.isfinite(L_phys)) or (L_phys <= 0.0):
                raise FloatingPointError("Invalid wavelength.")
            if (not np.isfinite(c_phys)) or (c_phys <= 0.0):
                raise FloatingPointError("Invalid celerity.")

            # Surface nodes correspond to X_m=m*pi/N over the symmetric half-wave
            eta_nodes = np.zeros(self.n + 1, dtype=DTYPE)
            for m in range(0, self.n + 1):
                kEta = float(self.z[10 + m])             # zeta_m = k*eta_m at the free-surface node
                eta_nodes[m] = self.d * (1.0 + kEta / kd)

            self.eta_nodes = eta_nodes
            self.k = float(k_phys)
            self.L = float(L_phys)
            self.c = float(c_phys)

            z3_period = self.T_target * np.sqrt(self.g * k_phys)
            z4_period = c_phys * np.sqrt(k_phys / self.g)
            if (
                np.isfinite(z3_period)
                and np.isfinite(z4_period)
                and z3_period > 0.0
                and z4_period > 0.0
            ):
                self.z[3] = z3_period
                self.z[4] = z4_period
                if self.Current_criterion == 1:
                    self.z[5] = self.Uc * np.sqrt(k_phys / self.g)
                    self.z[7] = self.z[4] - self.z[5]
                    self.z[6] = self.z[8] / self.z[1] - self.z[7] + self.z[4]
                else:
                    self.z[6] = self.Uc * np.sqrt(k_phys / self.g)
                    self.z[7] = self.z[8] / self.z[1] - self.z[6] + self.z[4]
                    self.z[5] = self.z[4] - self.z[7]

            # Store B_j as a 0-based array for GUI/report access
            self.Bj = self.B[1:self.n + 1].copy()

            # Crest and trough elevations relative to still-water level
            self.eta_crest = float(self.eta_nodes[0] - self.d)
            self.eta_trough = float(self.eta_nodes[-1] - self.d)

            # Engineering nondimensional descriptors based on H, L and d
            self.steepness = self.H_target / self.L
            self.rel_depth = self.d / self.L
            self.ursell = (self.H_target * self.L * self.L) / (self.d ** 3)

            # Depth-regime classification using d/L
            if self.rel_depth < 0.05:
                self.regime = "Shallow"
            elif self.rel_depth < 0.5:
                self.regime = "Intermediate"
            else:
                self.regime = "Deep"

            # Miche limiting-height check used as an engineering warning
            self.breaking_limit_miche = float(0.142 * self.L * np.tanh(self.k * self.d))
            self.breaking_index = float(self.H_target / self.breaking_limit_miche) if self.breaking_limit_miche > 0 else 0.0
            self.is_breaking = bool(self.breaking_limit_miche > 0 and self.H_target > self.breaking_limit_miche)

            # Integral quantities from the k-based and d-based invariant scalings
            self._calc_integral_props_cpp()

            # Kinematic summary at bed, crest and trough phases
            self.u_bed, _, _ = self.get_kinematics(0.0, 0.0)

            # Quadratic bed-shear estimate is an engineering diagnostic, not a Fenton residual
            cf_est = 0.005
            self.tau_bed = 0.5 * RHO * cf_est * (self.u_bed ** 2)

            self.ExcursionBed = abs(self.u_bed) * self.T_target / (2.0 * np.pi)

            # Crest and trough surface velocities in the fixed frame
            self.u_surf, _, _ = self.get_kinematics(self.d + self.eta_crest, 0.0)
            u_trough, _, _ = self.get_kinematics(self.d + self.eta_trough, np.pi)
            self.asymmetry = abs(self.u_surf / u_trough) if abs(u_trough) > 0 else 0.0

            # Phase scan for maximum surface vertical velocity and horizontal acceleration
            scan_phases = np.linspace(0.0, np.pi, 40)
            max_ax = 0.0
            max_w  = 0.0
            for X in scan_phases:
                kEta = self._surface_keta(X)
                z_surf = self.d * (1.0 + kEta / kd)  # absolute from bed [m]
                _, w, ax = self.get_kinematics(z_surf, X)
                max_ax = max(max_ax, abs(ax))
                max_w  = max(max_w, abs(w))

            self.acc_max = float(max_ax)
            self.w_max   = float(max_w)

            self.converged = True

        except FloatingPointError as e:
            self.last_error = f"Floating point failure: {e}"
            self.converged = False
        except Exception as e:
            self.last_error = f"Solver error: {e}"
            self.converged = False
        finally:
            np.seterr(**old_err)

    def get_kinematics(self, z_bed, phase=0.0):
        """
        GUI-facing kinematics: (u_abs, w_abs, a_x) at a given vertical position.

        Parameters
        ----------
        z_bed : float
            Vertical coordinate from the bed [m]. Bed=0, mean level=d.
        phase : float
            Phase angle X = kx in radians (0 at crest, pi at trough for half-wave).

        Returns
        -------
        u_abs : float [m/s]
        w_abs : float [m/s]
        ax    : float [m/s^2]
        """
        kd = float(self.z[1])
        if kd <= 0.0 or self.d <= 0.0:
            return 0.0, 0.0, 0.0

        k_phys = kd / self.d
        X = float(phase)
        Y = float(k_phys * (float(z_bed) - self.d))  # wave scaling: k(z-d)

        u_nd, v_nd, dudt_nd = self._point(X, Y)

        u_abs = u_nd * np.sqrt(self.g * self.d)
        w_abs = v_nd * np.sqrt(self.g * self.d)
        ax = dudt_nd * self.g

        return float(u_abs), float(w_abs), float(ax)

    def _mean_square_bed_orbital_velocity(self, nph=720):
        """
        Mean square near-bed *orbital* horizontal velocity [m^2/s^2].

        Definition adopted (non-negative by construction):
            u_b^2 = < (u_b(t) - ū₁)^2 >
        where ū₁ is the Eulerian current (Uc), and <·> denotes averaging over
        one wave period (equivalently one wavelength for a steady progressive wave).

        Notes
        -----
        - This matches the standard RMS-orbital-velocity concept used in coastal
          engineering and sediment/force calculations: U_rms = sqrt(<u^2>).
        - Computed numerically by sampling phase uniformly.
        """
        if self.d <= 0.0 or self.T_target <= 0.0:
            return 0.0

        # Sample one full cycle in phase X. For a steady periodic wave,
        # spatial averaging over one wavelength is equivalent to temporal averaging.
        phases = np.linspace(0.0, 2.0 * np.pi, int(max(36, nph)), endpoint=False)

        ub2 = 0.0
        for ph in phases:
            u_abs, _, _ = self.get_kinematics(z_bed=0.0, phase=float(ph))  # bed: z_bed=0
            u_orb = u_abs - float(self.Uc)  # Remove imposed Eulerian current to isolate the orbital component
            ub2 += u_orb * u_orb

        return float(ub2 / len(phases))


    # ------------------------ integral quantities and depth-scaled output ---------

    def _momentum_flux_S_depth(self, phase=0.0, npts=1200):
        """
        Compute the depth-scaled momentum flux S/(ρ g d²) in the *moving frame*:
            S = ∫₀^{η} [ p + ρ (u-c)² ] dz   (per unit crest width)

        This is the quantity printed by Fenton in Solution.res as:
            Momentum flux      S/(ρ g d²)

        Notes
        -----
        - Uses the same non-dimensional pressure from Bernoulli as in Point().
        - Evaluated at a single phase; the result is invariant with phase for the steady solution.
        """
        kd = float(self.z[1])
        if kd <= 0.0:
            return 0.0

        c = float(self.z[4] / np.sqrt(kd))          # c/√(g d)
        R = float(1.0 + self.z[9] / kd)            # R/(g d)

        X = float(phase)
        kEta = float(self._surface_keta(X))
        eta_over_d = 1.0 + kEta / kd               # y = z/d at free surface

        ys = np.linspace(0.0, eta_over_d, int(max(50, npts)), dtype=DTYPE)
        integ = np.zeros_like(ys)
        for idx, y in enumerate(ys):
            Y = kd * (y - 1.0)                     # Y = k(z-d) = kd*(y-1)
            u_nd, v_nd, _ = self._point(X, Y)      # u,v scaled by √(g d)
            urel = u_nd - c
            # Pressure scaled by rho*g*d from the Bernoulli dynamic condition
            P = R - y - 0.5 * (urel * urel + v_nd * v_nd)
            integ[idx] = P + urel * urel

        return float(_np_trapz(integ, ys))

    def _calc_integral_props_cpp(self):
        """
        Compute integral quantities using the Fenton invariant scalings.
        Values are dimensionalised to match the GUI units.
        """
        kd = float(self.z[1])
        if kd <= 0.0:
            self.Power = self.EnergyDensity = self.Sxx = self.Impulse = self.Cg = 0.0
            self.MassTransport = 0.0
            self.BernoulliR = 0.0
            return

        # depth-scaled dimensionless bulk quantities
        c_dimless  = float(self.z[4] / np.sqrt(kd))
        ce_dimless = float(self.z[5] / np.sqrt(kd))
        cs_dimless = float(self.z[6] / np.sqrt(kd))
        ubar_dimless = float(self.z[7] / np.sqrt(kd))

        Q_dimless = float(ubar_dimless - self.z[8] / (kd ** 1.5))
        R_dimless = float(1.0 + self.z[9] / kd)

        # k-based invariants used by the Fenton output convention
        pulse = float(self.z[8] + kd * self.z[5])
        ke = 0.5 * (self.z[4] * pulse - self.z[5] * Q_dimless * (kd ** 1.5))

        pe = 0.0
        for i in range(1, self.n + 1):
            pe += 0.25 * (self.Y[i] ** 2)

        ub2 = float(2.0 * self.z[9] - self.z[4] * self.z[4])
        q_term = float(self.z[7] * kd - self.z[8])

        sxx = float(4.0 * ke - 3.0 * pe + ub2 * kd + 2.0 * self.z[5] * q_term)
        f = float(self.z[4] * (3.0 * ke - 2.0 * pe) + 0.5 * ub2 * (pulse + self.z[4] * kd) + self.z[4] * self.z[5] * q_term)

        # Convert to the depth-scaled nondimensional output system
        E_depth = float((ke + pe) / (kd ** 2))
        KE_depth = float(ke / (kd ** 2))
        PE_depth = float(pe / (kd ** 2))

        # Store d-based invariants for Solution-Flat-style reporting
        self.E_depth = E_depth
        self.KE_depth = KE_depth
        self.PE_depth = PE_depth
        Sxx_depth = float(sxx / (kd ** 2))
        F_depth = float(f / (kd ** 2.5))
        I_depth = float(pulse / (kd ** 1.5))

        self.Sxx_depth = Sxx_depth
        self.F_depth = F_depth
        self.I_depth = I_depth

        # Dimensionalise d-based quantities to engineering GUI units
        self.EnergyDensity = float(RHO * self.g * (self.d ** 2) * E_depth)         # [J/m^2]
        self.Sxx = float(RHO * self.g * (self.d ** 2) * Sxx_depth)                 # [N/m]
        self.Power = float(RHO * (self.g ** 1.5) * (self.d ** 2.5) * F_depth)      # [W/m]
        # Momentum flux S/(rho*g*d^2) in the moving-frame convention
        self.MomentumFluxDepth = self._momentum_flux_S_depth(phase=0.0, npts=1200)  # S/(ρ g d²)
        self.MomentumFlux = float(RHO * self.g * (self.d ** 2) * self.MomentumFluxDepth)             # [N/m]
        self.Impulse = float(RHO * np.sqrt(self.g * (self.d ** 3)) * I_depth)  # [kg/(m·s)] per unit crest width

        self.BernoulliR = float(R_dimless * self.g * self.d)                       # [m^2/s^2] head* g? (consistent scalar)
        self.MassTransport = float(cs_dimless * np.sqrt(self.g * self.d))          # [m/s] (Stokes current)

        # Values reported using Fenton current, flux and Bernoulli notation
        self.EulerianCurrent = float(self.Uc)                              # u1 [m/s]
        self.StokesCurrent = float(self.MassTransport)                      # u2 [m/s]
        self.MeanFluidSpeed = float(ubar_dimless * np.sqrt(self.g * self.d))  # Ū [m/s]
        self.VolumeFluxQ = float(Q_dimless * np.sqrt(self.g * (self.d ** 3))) # Q [m^2/s]
        self.WaveVolumeFlux_q = float(self.MeanFluidSpeed * self.d - self.VolumeFluxQ)  # q [m^2/s]
        self.BernoulliR_dimless = float(R_dimless)                          # R/(g d)
        self.Bernoulli_r = float((R_dimless - 1.0) * self.g * self.d)        # r = R - g d [m^2/s^2]
        self.KineticEnergy = float(RHO * self.g * (self.d ** 2) * self.KE_depth)     # [J/m^2]
        self.PotentialEnergy = float(RHO * self.g * (self.d ** 2) * self.PE_depth)  # [J/m^2]
        # Mean square bed orbital velocity for the Solution-Flat output row:
        # u_b^2 = <(u_b(t)-ubar_1)^2>, removing the Eulerian current component
        self.MeanSquareBedVelocity = float(self._mean_square_bed_orbital_velocity(nph=720))   # [m^2/s^2]

        self.Cg = float(self.Power / self.EnergyDensity) if abs(self.EnergyDensity) > 1e-12 else 0.0


# ==============================================================================
#  GUI APPLICATION CLASS
# ==============================================================================

class FentonApp(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("Fenton Wave Solver")
        sw = self.winfo_screenwidth()
        sh = self.winfo_screenheight()
        w = min(1600, sw - 80)
        h = min(950,  sh - 120)
        self.geometry(f"{w}x{h}")
        
        # --- Style Configuration ---
        style = ttk.Style(self)
        style.theme_use('clam')
        
        # --- Input Frame ---
        input_frame = ttk.LabelFrame(self, text="Wave Parameters", padding="10")
        input_frame.pack(fill="x", padx=10, pady=5)
        
        # Grid layout for inputs
        ttk.Label(input_frame, text="Wave Height (H) [m]:").grid(row=0, column=0, padx=5, pady=5, sticky="e")
        self.ent_H = ttk.Entry(input_frame)
        self.ent_H.grid(row=0, column=1, padx=5, pady=5)
        self.ent_H.insert(0, "3.0") # Default H=3
        
        ttk.Label(input_frame, text="Wave Period (T) [s]:").grid(row=0, column=2, padx=5, pady=5, sticky="e")
        self.ent_T = ttk.Entry(input_frame)
        self.ent_T.grid(row=0, column=3, padx=5, pady=5)
        self.ent_T.insert(0, "9.0") # Default T=9
        
        ttk.Label(input_frame, text="Water Depth (d) [m]:").grid(row=1, column=0, padx=5, pady=5, sticky="e")
        self.ent_d = ttk.Entry(input_frame)
        self.ent_d.grid(row=1, column=1, padx=5, pady=5)
        self.ent_d.insert(0, "5.0") # Default d=5
        
        ttk.Label(input_frame, text="Current Vel (Uc) [m/s]:").grid(row=1, column=2, padx=5, pady=5, sticky="e")
        self.ent_Uc = ttk.Entry(input_frame)
        self.ent_Uc.grid(row=1, column=3, padx=5, pady=5)
        self.ent_Uc.insert(0, "1.0") # Default U=1
        
        # Calculate Button
        self.btn_calc = ttk.Button(input_frame, text="CALCULATE", command=self.run_calculation)
        self.btn_calc.grid(row=0, column=4, rowspan=2, padx=20, sticky="ns")
        
        # --- Output Frame ---
        output_frame = ttk.LabelFrame(self, text="Output Console", padding="10")
        output_frame.pack(fill="both", expand=True, padx=10, pady=5)
        
        # Scrolled Text Widget for results with INCREASED FONT (14)
        self.txt_output = scrolledtext.ScrolledText(output_frame, font=("Courier New", 14), state='disabled')
        self.txt_output.pack(fill="both", expand=True)
        
        # --- Right-Click Context Menu ---
        self.context_menu = Menu(self, tearoff=0)
        self.context_menu.add_command(label="Copy", command=self.copy_to_clipboard)
        self.context_menu.add_command(label="Select All", command=self.select_all)
        
        # Bind Right Click (Button-3 on Windows/Linux, Button-2 on Mac often)
        self.txt_output.bind("<Button-3>", self.show_context_menu)
        
    def show_context_menu(self, event):
        self.context_menu.post(event.x_root, event.y_root)

    def copy_to_clipboard(self):
        try:
            self.clipboard_clear()
            text = self.txt_output.get("sel.first", "sel.last")
            self.clipboard_append(text)
        except tk.TclError:
            pass # Handle case where no text is selected

    def select_all(self):
        self.txt_output.tag_add("sel", "1.0", "end")

    def log(self, text, file_handle=None):
        """Writes text to the GUI text widget and optionally a file."""
        self.txt_output.config(state='normal')
        self.txt_output.insert(tk.END, text)
        self.txt_output.see(tk.END)
        self.txt_output.config(state='disabled')
        
        if file_handle:
            try:
                file_handle.write(text)
            except UnicodeEncodeError:
                # Fallback: never abort GUI just because the host encoding cannot represent symbols.
                enc = getattr(file_handle, "encoding", "utf-8") or "utf-8"
                file_handle.write(text.encode(enc, errors="replace").decode(enc, errors="replace"))
            except Exception:
                pass

    def _finalize_error(self, title: str, details: str):
        """Render a failure state safely in the GUI (main thread only)."""
        # Restore output area
        self.txt_output.config(state="normal")
        self.txt_output.delete(1.0, tk.END)
        self.txt_output.config(state="disabled")

        # Re-enable button so the user can try again
        self.btn_calc.config(state="normal")

        # Print the error
        self.log(f"{title}\n\n")
        if details:
            self.log(details + "\n")

        # If stdout/stderr were redirected (frozen/no-console), include buffered output.
        try:
            buf = _STDIO_BUFFER.getvalue()
            if buf.strip():
                self.log("\n--- Buffered stdout/stderr ---\n")
                self.log(buf + "\n")
        except Exception:
            pass

    def run_calculation(self):

        """1. PREPARE: Validates inputs and starts the background thread."""
        # Clear previous output
        self.txt_output.config(state='normal')
        self.txt_output.delete(1.0, tk.END)
        self.txt_output.config(state='disabled')
        
        # Get and Validate Inputs
        try:
            H_in = float(self.ent_H.get())
            T_in = float(self.ent_T.get())
            d_in = float(self.ent_d.get())
            Uc_in = float(self.ent_Uc.get())
            
            if H_in < 0 or T_in < 0 or d_in < 0:
                raise ValueError("H, T, and d must be positive.")
        except ValueError:
            messagebox.showerror("Input Error", "Please enter valid numeric values.")
            return

        # Disable button to prevent double-clicking while running
        self.btn_calc.config(state='disabled')
        self.log("Running solver... please wait.\n")
        
        # Start the heavy math in a SEPARATE thread
        thread = threading.Thread(target=self._solve_in_background, 
                                  args=(H_in, T_in, d_in, Uc_in))
        thread.daemon = True # Ensures thread dies if we close the window
        thread.start()


    def _solve_in_background(self, H, T, d, Uc):
        """Heavy computations (background thread). All GUI updates are scheduled via after()."""
        import traceback
        try:
            # Case A: No current
            solver0 = FentonStreamFunction(H, T, d, Uc=0.0)
            solver0.solve()

            # Case B: With ambient current (Eulerian)
            solverC = FentonStreamFunction(H, T, d, Uc=Uc)
            solverC.solve()

            # Schedule UI update on the main thread
            self.after(0, self._finalize_output, solver0, solverC, H, T, d, Uc)

        except Exception:
            self.after(0, self._finalize_error, "Solver failed (exception).", traceback.format_exc())


    def _finalize_output(self, solver0, solverC, H_in, T_in, d_in, Uc_in):
        """3. DISPLAY: Prints results (Runs back on Main Thread)."""

        # Clear the "Running..." text
        self.txt_output.config(state='normal')
        self.txt_output.delete(1.0, tk.END)
        self.txt_output.config(state='disabled')

        # Open file logging (UTF-8 to preserve scientific symbols and avoid truncation)
        try:
            f_log = open("output.txt", "w", encoding="utf-8", newline="\n")
        except Exception:
            f_log = None

        def log_wrapper(s: str):
            self.log(s, f_log)


        # ---------------------------- numerical sanity checks ----------------------------
        # Do not attempt to format/report invalid solutions. Instead, surface the reason
        # and return control to the user.
        def _solver_status(slv, label):
            if not getattr(slv, "converged", False):
                msg = getattr(slv, "last_error", "") or "Did not converge."
                return f"[{label}] {msg}"
            # Also guard against silent NaNs/Infs in key outputs
            for attr in ("L", "k", "c"):
                v = getattr(slv, attr, None)
                try:
                    if v is None or (not np.isfinite(float(v))):
                        return f"[{label}] Non-finite result in {attr}."
                except Exception:
                    return f"[{label}] Invalid result in {attr}."
            return ""

        issues = []
        s0 = _solver_status(solver0, "No current")
        if s0:
            issues.append(s0)
        if float(Uc_in) != 0.0:
            sC = _solver_status(solverC, "With current")
            if sC:
                issues.append(sC)

        if issues:
            # Close file logging before returning
            if f_log:
                try:
                    f_log.close()
                except Exception:
                    pass
            self._finalize_error("Numerical failure / non-convergence.", "\n".join(issues))
            return

        # ---------------------------- formatting helpers ----------------------------

        W = 107  # consistent width for ALL boxes/tables in the report

        def hline(ch: str = "-"):
            log_wrapper("+" + (ch * (W - 2)) + "+\n")

        def box_title(title: str):
            hline("-")
            log_wrapper("|" + f"{title}".center(W - 2) + "|\n")
            hline("-")

        def box_text(text: str):
            s = (text or "").replace("\n", " ").rstrip()
            if len(s) > (W - 2):
                s = s[: (W - 5)] + "..."
            log_wrapper("|" + s.ljust(W - 2) + "|\n")

        def fmt_float(v, w):
            """Fit a float into width w (right-aligned), preferring fixed decimals."""
            if v is None:
                return "-".rjust(w)
            try:
                v = float(v)
            except Exception:
                return str(v)[:w].rjust(w)

            # Try fixed-point with decreasing decimals
            for dec in (5, 4, 3, 2, 1, 0):
                s = f"{v:.{dec}f}"
                if len(s) <= w:
                    return s.rjust(w)

            # Fallback to scientific
            for sig in (6, 5, 4, 3):
                s = f"{v:.{sig}e}"
                if len(s) <= w:
                    return s.rjust(w)

            return f"{v:.2e}"[:w].rjust(w)

        def fmt_cell(v, w, align="right"):
            if v is None:
                s = "-"
            elif isinstance(v, str):
                s = v
            else:
                return fmt_float(v, w)
            s = s.replace("\n", " ").strip()
            if len(s) > w:
                s = s[:w]
            if align == "left":
                return s.ljust(w)
            if align == "center":
                return s.center(w)
            return s.rjust(w)

        def table_sep(col_w):
            # internal separator line (keeps box width W)
            parts = []
            for w in col_w:
                parts.append("-" * (w + 2))
            log_wrapper("|" + "+".join(parts) + "|\n")

        def print_table(title, headers, col_w, aligns, rows):
            """Print a boxed table with wrapping and optional section-rows."""
            box_title(title)
            # Header
            line = "|"
            for h, w, al in zip(headers, col_w, aligns):
                line += " " + fmt_cell(h, w, align=al) + " |"
            log_wrapper(line + "\n")
            table_sep(col_w)

            for r in rows:
                # Section header row (spans the full width)
                if isinstance(r, dict) and r.get("_section"):
                    sec = f" {r['_section']} "
                    log_wrapper("|" + sec.center(W - 2, "-") + "|\n")
                    continue

                # Keep numeric cells numeric so we can format them consistently
                # and *prevent* textwrap from splitting long float representations
                # across multiple lines (which produced the “extra digits” lines).
                cells = list(r)

                wrapped = []
                for c, w in zip(cells, col_w):
                    # Numeric cells: do not wrap; formatting happens in fmt_cell/fmt_float.
                    if c is None:
                        wrapped.append(["-"])
                    elif isinstance(c, (int, float, np.integer, np.floating)) and not isinstance(c, bool):
                        wrapped.append([c])
                    else:
                        s = str(c).replace("\n", " ").rstrip()
                        wrapped.append(
                            textwrap.wrap(
                                s,
                                width=w,
                                break_long_words=False,
                                break_on_hyphens=False
                            ) or [""]
                        )

                nlines = max(len(x) for x in wrapped)
                for i in range(nlines):
                    line = "|"
                    for (wlines, w, al) in zip(wrapped, col_w, aligns):
                        seg = wlines[i] if i < len(wlines) else ""
                        line += " " + fmt_cell(seg, w, align=al) + " |"
                    log_wrapper(line + "\n")

            hline("-")
            log_wrapper("\n")

        has_current = (float(Uc_in) != 0.0)

        # ------------------------------- report header ------------------------------
        box_title("NONLINEAR WAVE HYDRODYNAMICS SOLVER (FENTON)")
        box_text(f"Wave height (H)             : {H_in} m")
        box_text(f"Wave period (τ)             : {T_in} s")
        box_text(f"Water depth (d)             : {d_in} m")
        box_text(f"Eulerian current ū₁         : {Uc_in} m/s (positive with wave propagation)")
        hline("-")
        box_text("Status: Full nonlinear system solved successfully.")
        hline("-")
        log_wrapper("\n")

        # ---------------------------- hydrodynamic summary ---------------------------
        # Column widths for 4-col (Param / NoCurrent / WithCurrent / Unit)
        w_param, w_nc, w_wc, w_unit = 42, 16, 16, 20
        headers = ["PARAMETER", "NO CURRENT", "WITH CURRENT", "UNIT"]
        col_w = [w_param, w_nc, w_wc, w_unit]
        aligns = ["left", "right", "right", "left"]

        def wc(v):
            return v if has_current else "-"

        # Pre-compute a few helpful scalars
        g = float(solver0.g)
        d = float(solver0.d)
        sqrt_gd = float(np.sqrt(g * d))
        sqrt_g_over_d = float(np.sqrt(g / d))

        rows = []
        rows.append({"_section": "INPUTS & REFERENCE SCALES"})
        rows += [
            ("Water depth (d)", solver0.d, wc(solverC.d), "m"),
            ("Wave height (H)", solver0.H_target, wc(solverC.H_target), "m"),
            ("Wave period (τ)", solver0.T_target, wc(solverC.T_target), "s"),
            ("H/d", solver0.H_target / solver0.d, wc(solverC.H_target / solverC.d), "-"),
            ("τ√(g/d)", solver0.T_target * sqrt_g_over_d, wc(solverC.T_target * sqrt_g_over_d), "-"),
        ]

        def report_celerity(slv):
            if slv.T_target > 0.0 and np.isfinite(slv.L):
                return slv.L / slv.T_target
            return slv.c

        c0_report = report_celerity(solver0)
        cC_report = report_celerity(solverC)

        rows.append({"_section": "DISPERSION & PHASE (GEOMETRY)"})
        rows += [
            ("Wavelength (L)", solver0.L, wc(solverC.L), "m"),
            ("Wave number (k)", solver0.k, wc(solverC.k), "rad/m"),
            ("kd", float(solver0.z[1]), wc(float(solverC.z[1])), "-"),
            ("Angular frequency (ω)", 2.0 * np.pi / solver0.T_target, wc(2.0 * np.pi / solverC.T_target), "rad/s"),
            ("Celerity / phase speed (c)", c0_report, wc(cC_report), "m/s"),
            ("c/√(gd)", c0_report / sqrt_gd, wc(cC_report / sqrt_gd), "-"),
            ("Crest elevation (ηc)", solver0.eta_crest, wc(solverC.eta_crest), "m"),
            ("Trough elevation (ηt)", solver0.eta_trough, wc(solverC.eta_trough), "m"),
        ]

        rows.append({"_section": "MEAN FLOWS (FENTON SOLUTION-FLAT)"})
        rows += [
            ("Eulerian current (ū₁)", solver0.EulerianCurrent, wc(solverC.EulerianCurrent), "m/s"),
            ("Stokes current (ū₂)", solver0.StokesCurrent, wc(solverC.StokesCurrent), "m/s"),
            ("Mean fluid speed (Ū)", solver0.MeanFluidSpeed, wc(solverC.MeanFluidSpeed), "m/s"),
        ]

        rows.append({"_section": "FLUXES & BERNOULLI CONSTANTS"})
        rows += [
            ("Wave volume flux (q)", solver0.WaveVolumeFlux_q, wc(solverC.WaveVolumeFlux_q), "m²/s"),
            ("Volume flux (Q)", solver0.VolumeFluxQ, wc(solverC.VolumeFluxQ), "m²/s"),
            ("Bernoulli constant (R)", solver0.BernoulliR, wc(solverC.BernoulliR), "m²/s²"),
            ("Reduced Bernoulli (r = R−g d)", solver0.Bernoulli_r, wc(solverC.Bernoulli_r), "m²/s²"),
        ]

        rows.append({"_section": "INTEGRAL QUANTITIES (PER UNIT CREST WIDTH)"})
        rows += [
            ("Kinetic energy (T)", solver0.KineticEnergy / 1000.0, wc(solverC.KineticEnergy / 1000.0), "kJ/m²"),
            ("Potential energy (V)", solver0.PotentialEnergy / 1000.0, wc(solverC.PotentialEnergy / 1000.0), "kJ/m²"),
            ("Total energy (E = T+V)", solver0.EnergyDensity / 1000.0, wc(solverC.EnergyDensity / 1000.0), "kJ/m²"),
            ("Momentum flux (S)", solver0.MomentumFlux / 1000.0, wc(solverC.MomentumFlux / 1000.0), "kN/m"),
            ("Radiation stress (Sₓₓ)", solver0.Sxx / 1000.0, wc(solverC.Sxx / 1000.0), "kN/m"),
            ("Impulse (I)", solver0.Impulse / 1000.0, wc(solverC.Impulse / 1000.0), "10³ kg/(m·s)"),
            ("Wave power (F)", solver0.Power / 1000.0, wc(solverC.Power / 1000.0), "kW/m"),
            ("Group velocity (C𝗀 = F/E)\u3164", solver0.Cg, wc(solverC.Cg), "m/s"),
        ]

        rows.append({"_section": "KINEMATICS (EXTREMES / BED ORBITAL MOTION)"})
        rows += [
            ("Max surface horiz. vel |u|", solver0.u_surf, wc(solverC.u_surf), "m/s"),
            ("Max bed horiz. vel |u|", solver0.u_bed, wc(solverC.u_bed), "m/s"),
            ("Max horiz. accel |aₓ|", solver0.acc_max, wc(solverC.acc_max), "m/s²"),
            ("Velocity asymmetry |uc|/|ut|", solver0.asymmetry, wc(solverC.asymmetry), "-"),
            ("Mean square bed orbital vel ub²", solver0.MeanSquareBedVelocity, wc(solverC.MeanSquareBedVelocity), "m²/s²"),
            ("Bed orbital RMS velocity ub,rms", np.sqrt(max(0.0, solver0.MeanSquareBedVelocity)), wc(np.sqrt(max(0.0, solverC.MeanSquareBedVelocity))), "m/s"),
        ]

        rows.append({"_section": "NONLINEARITY / BREAKING DIAGNOSTICS"})
        warn0 = "BREAKING" if solver0.is_breaking else "STABLE"
        warnC = ("BREAKING" if solverC.is_breaking else "STABLE") if has_current else "-"
        rows += [
            ("Miche breaking limit (Hmax)", solver0.breaking_limit_miche, wc(solverC.breaking_limit_miche), "m"),
            ("Saturation (H/Hmax)", solver0.breaking_index, wc(solverC.breaking_index), "-"),
            ("Breaking status", warn0, warnC, "-"),
            ("Ursell number (U)", solver0.ursell, wc(solverC.ursell), "-"),
            ("Regime (by d/L)", solver0.regime, wc(solverC.regime), "-"),
        ]

        print_table("CALCULATED HYDRODYNAMIC PARAMETERS", headers, col_w, aligns, rows)

        # ------------------------ SOLUTION-FLAT tables (exact set) -------------------
        def print_solution_flat(slv, title):
            # 5-column widths (sum widths + 16 separators/spaces = W)
            w_idx, w_name, w_val, w_adim, w_adval = 2, 37, 13, 25, 14
            headers = ["#", "PARAMETER", "value", "adim param", "adim value"]
            col_w = [w_idx, w_name, w_val, w_adim, w_adval]
            aligns = ["right", "left", "right", "left", "right"]

            g = float(slv.g)
            d = float(slv.d)
            H = float(slv.H_target)
            T = float(slv.T_target)
            L = float(slv.L)
            c = float(L / T) if T > 0.0 and np.isfinite(L) else float(slv.c)
            sqrt_gd = float(np.sqrt(g * d))
            sqrt_gd3 = float(np.sqrt(g * (d ** 3)))

            def kJ(J): return float(J) / 1000.0
            def kN(N): return float(N) / 1000.0
            def kW(W_): return float(W_) / 1000.0

            # Exact 19-line set (Fenton Solution-Flat semantics & scalings).
            rows = [
                (1,  "Water depth",                          f"{d:.5f}",                       "d/d = 1",                         f"{1.0:.5f}"),
                (2,  "Wave length",                          f"{L:.5f}",                       "λ/d",                             f"{(L/d):.5f}"),
                (3,  "Wave height",                          f"{H:.5f}",                       "H/d",                             f"{(H/d):.5f}"),
                (4,  "Wave period",                          f"{T:.5f}",                       "τ√(g/d)",                         f"{(T*np.sqrt(g/d)):.5f}"),
                (5,  "Wave speed",                           f"{c:.5f}",                       "c/√(gd)",                         f"{(c/sqrt_gd):.5f}"),
                (6,  "Eulerian current",                     f"{slv.EulerianCurrent:.5f}",     "ū₁/√(gd)",                        f"{(slv.EulerianCurrent/sqrt_gd):.5f}"),
                (7,  "Stokes current",                       f"{slv.StokesCurrent:.5f}",       "ū₂/√(gd)",                        f"{(slv.StokesCurrent/sqrt_gd):.5f}"),
                (8,  "Mean fluid speed",                     f"{slv.MeanFluidSpeed:.5f}",      "Ū/√(gd)",                         f"{(slv.MeanFluidSpeed/sqrt_gd):.5f}"),
                (9,  "Wave volume flux, q = Ū d − Q",        f"{slv.WaveVolumeFlux_q:.5f}",    "q/√(gd³)",                        f"{(slv.WaveVolumeFlux_q/sqrt_gd3):.5f}"),
                (10, "Bernoulli constant, r = R − gd",       f"{slv.Bernoulli_r:.5f}",         "r/gd",                            f"{(slv.Bernoulli_r/(g*d)):.5f}"),
                (11, "Volume flux",                          f"{slv.VolumeFluxQ:.5f}",         "Q/√(gd³)",                        f"{(slv.VolumeFluxQ/sqrt_gd3):.5f}"),
                (12, "Bernoulli constant",                   f"{slv.BernoulliR:.5f}",          "R/gd",                            f"{(slv.BernoulliR/(g*d)):.5f}"),
                (13, "Momentum flux",                        f"{kN(slv.MomentumFlux):.5f}",    "S/ρgd²",                          f"{(slv.MomentumFluxDepth):.5f}"),
                (14, "Impulse",                              f"{(slv.Impulse/1000.0):.5f}",    "I/(ρ√(gd³))",                     f"{(slv.I_depth):.5f}"),
                (15, "Kinetic energy",                       f"{kJ(slv.KineticEnergy):.5f}",   "T/ρgd²",                          f"{(slv.KE_depth):.5f}"),
                (16, "Potential energy",                     f"{kJ(slv.PotentialEnergy):.5f}", "V/ρgd²",                          f"{(slv.PE_depth):.5f}"),
                (17, "Mean square of bed velocity",          f"{slv.MeanSquareBedVelocity:.5f}","ub²/gd",                         f"{(slv.MeanSquareBedVelocity/(g*d)):.5f}"),
                (18, "Radiation stress",                     f"{kN(slv.Sxx):.5f}",             "S_xx/ρgd²",                       f"{(slv.Sxx_depth):.5f}"),
                (19, "Wave power",                           f"{kW(slv.Power):.5f}",           "F/(ρg³ᐟ²d⁵ᐟ²)\u3164\u3164",           f"{(slv.F_depth):.5f}"),
            ]

            # Convert rows to generic printer rows (already formatted strings for numeric cols).
            out_rows = [(a, b, c, d_, e) for (a, b, c, d_, e) in rows]
            print_table(title, headers, col_w, aligns, out_rows)

        print_solution_flat(solver0, "SOLUTION.RES (NO CURRENT)")
        if has_current:
            print_solution_flat(solverC, "SOLUTION.RES (WITH CURRENT)")

        # --------------------------------- glossary --------------------------------
        glossary_headers = ["TERM / SYMBOL", "MEANING", "UNITS / NONDIM"]
        gw1, gw2, gw3 = 14, 64, 19
        gcol_w = [gw1, gw2, gw3]
        galigns = ["left", "left", "left"]

        terms = [
            ("d", "Still-water depth (bed to mean water level). Reference length scale.", "m ; d/d=1"),
            ("H", "Wave height (crest-to-trough).", "m ; H/d"),
            ("τ", "Wave period.", "s ; τ√(g/d)"),
            ("L", "Wavelength (crest-to-crest).", "m ; L/d"),
            ("k", "Wave number, k = 2π/L.", "rad/m"),
            ("ω", "Angular frequency, ω = 2π/τ.", "rad/s"),
            ("c", "Phase speed (celerity).", "m/s ; c/√(gd)"),
            ("ηc, ηt", "Crest and trough elevations relative to still-water level.", "m"),
            ("ū₁", "Eulerian (depth-mean) current; ū₁ = Uc.", "m/s ; ū₁/√(gd)"),
            ("ū₂", "Stokes / mass-transport current from nonlinear solution.", "m/s ; ū₂/√(gd)"),
            ("Ū", "Mean fluid speed (depth-mean).", "m/s ; Ū/√(gd)"),
            ("Q", "Volume flux (depth-integrated).", "m²/s ; Q/√(gd³)"),
            ("q", "Wave volume flux, q = Ū d − Q.", "m²/s ; q/√(gd³)"),
            ("R", "Bernoulli constant.", "m²/s² ; R/(gd)"),
            ("r", "Reduced Bernoulli constant r = R − g d.", "m²/s² ; r/(gd)"),
            ("S", "Momentum flux (moving frame).", "kN/m ; S/(ρgd²)"),
            ("I", "Wave impulse per unit width.", "10³ kg/(m·s) ; I/(ρ√(gd³))"),
            ("T", "Kinetic energy density.", "kJ/m² ; T/(ρgd²)"),
            ("V", "Potential energy density.", "kJ/m² ; V/(ρgd²)"),
            ("E", "Total energy density E = T + V.", "kJ/m²"),
            ("F", "Wave power (energy flux).", "kW/m ; F/(ρg³ᐟ²d⁵ᐟ²)\u3164"),
            ("Cg", "Group velocity, defined here as Cg = F/E.", "m/s"),
            ("Sₓₓ", "Radiation stress component in wave direction.", "kN/m ; Sₓₓ/(ρgd²)"),
            ("ub²", "Mean square *orbital* bed velocity: ub² = <(ub(t) − ū₁)²>. Non-negative by definition; computed by phase averaging.", "m²/s² ; /gd"),
            ("ub,rms", "Root-mean-square orbital bed velocity: ub,rms = √(ub²).", "m/s"),
            ("usurf,max", "Maximum horizontal velocity at free surface (scanned over phase).", "m/s"),
            ("ubed,max", "Maximum horizontal velocity at seabed (scanned over phase).", "m/s"),
            ("a_x,max", "Maximum horizontal acceleration magnitude.", "m/s²"),
            ("Asymmetry", "Velocity asymmetry indicator |uc|/|ut|.", "-"),
            ("Hmax", "Miche breaking limit used as stability diagnostic.", "m"),
            ("Ursell", "Ursell number, a shallow-water nonlinearity measure.", "-"),
            ("Regime", "Depth regime based on d/L (deep/intermediate/shallow).", "-"),
        ]

        print_table("PARAMETER DEFINITIONS & GLOSSARY", glossary_headers, gcol_w, galigns, terms)

        # Close file logging
        if f_log:
            f_log.close()

        # Re-enable button so user can calculate again
        self.btn_calc.config(state='normal')

        # Scroll to top
        self.txt_output.yview_moveto(0)
if __name__ == "__main__":
    app = FentonApp()
    app.mainloop()