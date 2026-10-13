# ==============================================================================
#  ENGINEERING TECHNICAL REFERENCE & THEORETICAL FORMULATION
# ==============================================================================
#  PROGRAM:      Fenton Nonlinear Wave Suite - Python Tkinter GUI solver
#  FILE:         fenton_gui.py
#  VALIDATION:   See accompanying Lean proofs and numerical audit.
#  Formal proofs concern the real-valued model, not all IEEE-754 executions.
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
#  The hyperbolic ratios use exponential expressions without subtracting
#  large, nearly equal hyperbolic terms:
#
#      S_j(Y) = [exp(jY) - exp(-jY-2jkd)] / [1 + exp(-2jkd)],
#      C_j(Y) = [exp(jY) + exp(-jY-2jkd)] / [1 + exp(-2jkd)].
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
#  The requested order is N = 50. The solver reduces the order when the
#  equilibrated Jacobian exceeds its binary64 conditioning limit. The active
#  residual system always has 2N + 10 equations for the order actually used.
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
#  current iterate and solves J dz = -F(z) using analytic derivatives and a
#  column-equilibrated SVD. Failed height steps are retried with smaller increments.
#  The Fourier order is limited by the conditioning of the equilibrated Jacobian;
#  the report states the order used. Residual tolerances are independent of order.
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
#      python -m PyInstaller --noconfirm --clean --onefile --windowed --name fenton_gui.exe fenton_gui.py
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
# Outward-rounded arithmetic for independent finite-system evaluation.
from decimal import Decimal, Context, ROUND_FLOOR, ROUND_CEILING, ROUND_HALF_EVEN, DecimalException
from functools import lru_cache
from math import factorial


class _Interval:
    down = Context(prec=50, rounding=ROUND_FLOOR)
    up = Context(prec=50, rounding=ROUND_CEILING)
    nearest = Context(prec=50, rounding=ROUND_HALF_EVEN)

    def __init__(self, lower, upper=None):
        def exact(value):
            return Decimal.from_float(float(value)) if isinstance(value, (float, np.floating)) else Decimal(value)
        self.lo = exact(lower)
        self.hi = self.lo if upper is None else exact(upper)
        if not self.lo.is_finite() or not self.hi.is_finite() or self.lo > self.hi:
            raise FloatingPointError("Invalid interval endpoints.")

    @staticmethod
    def coerce(value):
        return value if isinstance(value, _Interval) else _Interval(value)

    def __add__(self, other):
        other = self.coerce(other)
        return _Interval(self.down.add(self.lo, other.lo), self.up.add(self.hi, other.hi))

    __radd__ = __add__

    def __neg__(self):
        return _Interval(self.hi.copy_negate(), self.lo.copy_negate())

    def __sub__(self, other):
        return self + (-self.coerce(other))

    def __rsub__(self, other):
        return self.coerce(other) + (-self)

    def __mul__(self, other):
        other = self.coerce(other)
        pairs = [(a, b) for a in (self.lo, self.hi) for b in (other.lo, other.hi)]
        return _Interval(min(self.down.multiply(a, b) for a, b in pairs),
                         max(self.up.multiply(a, b) for a, b in pairs))

    __rmul__ = __mul__

    def __truediv__(self, other):
        other = self.coerce(other)
        if other.lo <= 0 <= other.hi:
            raise FloatingPointError("Interval division includes zero.")
        return self * _Interval(self.down.divide(Decimal(1), other.hi),
                                self.up.divide(Decimal(1), other.lo))

    def __rtruediv__(self, other):
        return self.coerce(other) / self

    def square(self):
        lower = Decimal(0) if self.lo <= 0 <= self.hi else min(self.lo.copy_abs(), self.hi.copy_abs())
        upper = max(self.lo.copy_abs(), self.hi.copy_abs())
        return _Interval(self.down.multiply(lower, lower), self.up.multiply(upper, upper))

    def __pow__(self, exponent):
        if not isinstance(exponent, int) or exponent < 0:
            raise ValueError("Interval powers require a nonnegative integer.")
        result, base = _Interval(1), self
        while exponent:
            if exponent % 2:
                result = result * base
            base = base.square()
            exponent //= 2
        return result

    def exp(self):
        return _Interval(self.nearest.next_minus(self.nearest.exp(self.lo)),
                         self.nearest.next_plus(self.nearest.exp(self.hi)))

    def sqrt(self):
        if self.lo < 0:
            raise FloatingPointError("Negative square-root interval.")
        lower = max(Decimal(0), self.nearest.next_minus(self.nearest.sqrt(self.lo)))
        return _Interval(lower, self.nearest.next_plus(self.nearest.sqrt(self.hi)))

    def magnitude(self):
        return max(self.lo.copy_abs(), self.hi.copy_abs())

    def record(self):
        return [str(self.lo), str(self.hi)]


_PI_INTERVAL = _Interval("3.14159265358979323846", "3.14159265358979323847")


@lru_cache(maxsize=16)
def _collocation_intervals(order):
    """Enclose sin(q*pi/N), cos(q*pi/N) using Taylor remainders on [-pi,pi]."""
    table = []
    for index in range(2*order):
        q = index if index <= order else index-2*order
        x = _PI_INTERVAL*q/order
        sine, cosine = _Interval(0), _Interval(0)
        st, ct = x, _Interval(1)
        square = x.square()
        for m in range(32):
            sine, cosine = sine+st, cosine+ct
            st = -st*square/((2*m+2)*(2*m+3))
            ct = -ct*square/((2*m+1)*(2*m+2))
        radius = _Interval(x.magnitude())
        sin_error = (radius**65/factorial(65)).hi
        cos_error = (radius**64/factorial(64)).hi
        sine = sine+_Interval(sin_error.copy_negate(), sin_error)
        cosine = cosine+_Interval(cos_error.copy_negate(), cos_error)
        table.append((sine, cosine))
    return tuple(table)


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
# Residual evaluation and field reconstruction use scalar collocation loops.
# Numba JIT compiles these loops; the analytic Jacobian uses NumPy arrays.
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

            x = jj * zsurf
            # Numerical safeguard for diverging trial free-surface ordinates
            if x > 60.0 or x < -60.0:
                rhs[1] = np.inf
                return np.inf
            ep = np.exp(x)
            em = np.exp(-x - 2.0 * jj * kd)
            den = 1.0 + np.exp(-2.0 * jj * kd)
            S = (ep - em) / den
            C = (ep + em) / den

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

        # Scaled ratios avoid overflow and catastrophic cancellation at the bed.
        a = j * Y
        b = j * kd
        ep = np.exp(a)
        em = np.exp(-a - 2.0 * b)
        den = 1.0 + np.exp(-2.0 * b)
        C = (ep + em) / den
        S = (ep - em) / den

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
#   - Newton iteration with analytic Jacobian and equilibrated SVD solve
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

    def __init__(self, H, T, d, Uc=0.0, *, fourier_order=None):
        # Dimensional inputs: H, T, d and imposed Eulerian current Uc
        self.H_target = float(H)    # [m]
        self.T_target = float(T)    # [s]
        self.d        = float(d)    # [m]
        self.Uc       = float(Uc)   # [m/s] (Eulerian / lab-frame)

        # Physical and numerical constants
        self.g = G_STD
        self.N = N_FOURIER if fourier_order is None else int(fourier_order)
        self.requested_order = self.N
        self.jacobian_condition = 0.0
        self.continuation_fraction = 0.0
        self.numerical_note = ""
        self.accuracy_verified = False
        self.accuracy_metrics = {}
        self.finite_system_verified = False
        self.finite_system_certificate = {}
        self.solution_status = "not_solved"

        # Solver control: continuation and Newton iteration tolerances
        self.nstep  = 4        # continuation steps in wave height
        self.number = 40       # max Newton iterations per step
        self.crit   = 1.0e-8   # intermediate-step convergence factor (C++: crit)
        self.criter_final = 1.0e-10  # final-step convergence factor
        self.residual_tolerance = 1.0e-8
        self.residual_tolerance_final = 1.0e-10
        self.residual_max = float("inf")

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
            self.z[5] = self.Current * np.sqrt(self.z[1])
            self.z[6] = 0.0
        else:
            self.z[6] = self.Current * np.sqrt(self.z[1])
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
        jcoeff = j * coeff               # shape (n,)

        for m in range(0, n + 1):
            zsurf = z[10 + m]  # zeta_m = k*eta_m at this free-surface node

            x = j * zsurf
            if np.any(x > 60.0) or np.any(x < -60.0):
                raise FloatingPointError("Divergence: exp(j*zsurf) out of safe range.")
            ep = np.exp(x)
            em = np.exp(-x - 2.0 * j * z[1])
            den = 1.0 + np.exp(-2.0 * j * z[1])
            S = (ep - em) / den
            C = (ep + em) / den

            cosnm = self._cos_nm[m]
            sinnm = self._sin_nm[m]

            psi = float(np.sum(coeff * S * cosnm))
            u   = float(np.sum(jcoeff * C * cosnm))
            v   = float(np.sum(jcoeff * S * sinnm))

            rhs[m + 9] = psi - z[8] - z[7] * z[m + 10]
            rhs[n + m + 10] = 0.5 * ((-z[7] + u) ** 2 + v * v) + z[m + 10] - z[9]

        return float(np.dot(rhs[1:num + 1], rhs[1:num + 1]))
    def _jacobian(self):
        """Analytic derivatives of the full collocation residual vector."""
        n, z, j = self.n, self.z, self._j
        A = np.zeros((self.num + 1, self.num + 1), dtype=DTYPE)
        A[1, 1], A[1, 2] = -self.Hoverd, 1.0
        A[2, 2], A[2, 3] = 1.0, -2.0 * self.height * z[3]
        A[3, 3], A[3, 4] = z[4], z[3]
        A[4, 5], A[4, 7], A[4, 4] = 1.0, 1.0, -1.0
        A[5, 1] = z[6] + z[7] - z[4]
        A[5, 6], A[5, 7], A[5, 4], A[5, 8] = z[1], z[1], -z[1], -1.0
        A[6, self.Current_criterion + 4] = 1.0
        A[6, 1] = -0.5 * self.Current / np.sqrt(z[1])
        A[7, 10:n + 11] = 2.0
        A[7, 10] = A[7, n + 10] = 1.0
        A[8, 10], A[8, n + 10], A[8, 2] = 1.0, -1.0, -1.0

        B = z[n + 11:]
        a = z[10:n + 11, None] * j[None, :]
        ep = np.exp(a)
        depth_exp = np.exp(-2.0 * j * z[1])
        em = np.exp(-a - 2.0 * j * z[1])
        den = 1.0 + depth_exp
        S, C = (ep - em) / den, (ep + em) / den
        dS = 2.0 * j * (em + depth_exp * ep) / den ** 2
        dC = 2.0 * j * (depth_exp * ep - em) / den ** 2
        cos, sin = self._cos_nm, self._sin_nm
        u = (C * cos) @ (j * B)
        v = (S * sin) @ (j * B)
        uw = u - z[7]
        stream = slice(9, n + 10)
        dynamic = slice(n + 10, 2 * n + 11)
        modes = slice(n + 11, 2 * n + 11)
        A[stream, 1] = (dS * cos) @ B
        A[stream, 7] = -z[10:n + 11]
        A[stream, 8] = -1.0
        A[stream, modes] = S * cos
        A[dynamic, 1] = uw * ((dC * cos) @ (j * B)) + v * ((dS * sin) @ (j * B))
        A[dynamic, 7] = -uw
        A[dynamic, 9] = -1.0
        A[dynamic, modes] = j * (uw[:, None] * C * cos + v[:, None] * S * sin)
        m = np.arange(n + 1)
        A[9 + m, 10 + m] = uw
        A[n + 10 + m, 10 + m] = 1.0 + uw * ((S * cos) @ (j * j * B)) + v * ((C * sin) @ (j * j * B))
        return A[1:, 1:]

    def _svd_solve(self, A, b):
        """
        Solve the column-equilibrated system via SVD with truncation:
          wmin = wmax * 1e-15

        Numerical safety:
        - Reject NaNs/Infs before entering LAPACK.
        - Fallback to least-squares if SVD fails to converge.
        """
        if (not np.isfinite(A).all()) or (not np.isfinite(b).all()):
            raise FloatingPointError("Non-finite values in Jacobian system (A or b).")

        scale = np.max(np.abs(A), axis=0)
        scale[scale == 0.0] = 1.0
        scaled_A = A / scale
        try:
            U, s, Vt = np.linalg.svd(scaled_A, full_matrices=False)
            smax = np.max(s) if s.size else 0.0
            self.jacobian_condition = float(smax / s[-1]) if s[-1] > smax * 1.0e-300 else float("inf")
            wmin = smax * 1.0e-15

            s_inv = np.zeros_like(s)
            mask = s > wmin
            s_inv[mask] = 1.0 / s[mask]

            return (Vt.T @ (s_inv * (U.T @ b))) / scale

        except np.linalg.LinAlgError:
            # Conservative fallback: least-squares solution of the same linearised system.
            x, *_ = np.linalg.lstsq(scaled_A, b, rcond=1.0e-15)
            return x / scale



    def _newton(self, iter_count):
        """Damped Newton iteration with an analytic, equilibrated Jacobian."""
        n = self.n
        num = self.num

        # baseline residual
        ss0 = float(self._eqns(self.rhs1))
        if not np.isfinite(ss0):
            raise FloatingPointError("Non-finite residual norm at start of Newton step.")

        z0 = self.z.copy()

        A = self._jacobian()
        b = -self.rhs1[1:num + 1]

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
            try:
                ss1 = float(self._eqns(self.rhs2))
            except FloatingPointError:
                # A divergent trial is rejected; it is not a converged step.
                self.z[:] = z0
                alpha *= 0.5
                continue

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

    def _residual_max(self):
        """Re-evaluate every dimensionless equation at the CURRENT full state."""
        objective = self._eqns(self.rhs1)
        active = self.rhs1[1:self.num + 1]
        if not np.isfinite(objective) or not np.isfinite(active).all():
            raise FloatingPointError("Non-finite residual at accepted state.")
        self.residual_max = float(np.max(np.abs(active)))
        return self.residual_max


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

            a = j * Y
            b = j * kd
            ep = np.exp(a)
            em = np.exp(-a - 2.0 * b)
            den = 1.0 + np.exp(-2.0 * b)
            C = (ep + em) / den
            S = (ep - em) / den

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
        """Use the largest Fourier order supported by binary64 conditioning."""
        self.accuracy_verified = False
        self.accuracy_metrics = {}
        self.finite_system_verified = False
        self.finite_system_certificate = {}
        self.solution_status = "not_solved"
        requested = self.requested_order
        order = self.n
        while True:
            self._solve_at_order()
            if self.converged:
                A = self._jacobian()
                scale = np.maximum(np.max(np.abs(A), axis=0), np.finfo(DTYPE).tiny)
                singular = np.linalg.svd(A / scale, compute_uv=False)
                self.jacobian_condition = float(singular[0] / singular[-1]) if singular[-1] > singular[0] * 1.0e-300 else float("inf")
            if self.jacobian_condition <= 1.0e11 or order <= 12:
                break
            order = max(12, int(order * 0.8))
            controls = (self.nstep, self.number, self.crit, self.criter_final,
                        self.residual_tolerance, self.residual_tolerance_final,
                        self.Current_criterion)
            self.__init__(self.H_target, self.T_target, self.d, self.Uc, fourier_order=order)
            (self.nstep, self.number, self.crit, self.criter_final,
             self.residual_tolerance, self.residual_tolerance_final,
             self.Current_criterion) = controls
            self.requested_order = requested
        if self.converged and self.jacobian_condition > 1.0e11:
            self.converged = False
            self.last_error = "The Fourier system is too ill-conditioned for a reliable binary64 solution."
        if order != requested:
            self.numerical_note = (f"Fourier order used: {order} (requested {requested}); "
                                   "limited by numerical conditioning.")

        if self.converged:
            try:
                self._verify_accuracy()
            except (FloatingPointError, ValueError, DecimalException, np.linalg.LinAlgError) as exc:
                self.converged = False
                self.accuracy_verified = False
                self.solution_status = "accuracy_failure"
                self.last_error = "Accuracy verification failed: " + str(exc)
        elif self.solution_status == "not_solved":
            self.solution_status = "nonconvergence"

    def _solve_at_order(self):
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
        self.residual_max = float("inf")

        # Basic input screening (physical requirements)
        if (not np.isfinite([self.H_target, self.T_target, self.d, self.Uc]).all()
                or self.H_target <= 0.0 or self.T_target <= 0.0 or self.d <= 0.0):
            self.last_error = "Invalid inputs: H, T, d must be finite and > 0; Uc must be finite."
            return
        if self.Current_criterion not in (1, 2):
            self.last_error = "Current criterion must be 1 (Eulerian) or 2 (Stokes)."
            return

        old_err = np.geterr()
        try:
            # Make numerical faults explicit. Underflow is benign in this context.
            np.seterr(over="raise", invalid="raise", divide="raise", under="ignore")

            # Height continuation with rollback to the last converged state.
            fraction = 0.0
            increment = 1.0 / max(1, self.nstep)
            previous_fraction = 0.0
            previous_state = None
            accepted_state = None
            attempts = 0
            while fraction < 1.0:
                attempts += 1
                target = min(1.0, fraction + increment)
                self.height = target * self.Height
                self.Hoverd = target * self.MaxH
                if accepted_state is None:
                    self._init_linear()
                else:
                    self.z[:] = accepted_state
                    if previous_state is not None:
                        ratio = (target - fraction) / (fraction - previous_fraction)
                        self.z[:] += ratio * (accepted_state - previous_state)
                    if self.z[1] <= 0.0 or not np.isfinite(self.z).all():
                        self.z[:] = accepted_state

                step_converged = False
                residual_limit = (self.residual_tolerance_final if target == 1.0
                                  else self.residual_tolerance)
                criter = self.criter_final if target == 1.0 else self.crit
                for it in range(1, self.number + 1):
                    try:
                        err = self._newton(it)
                        residual = self._residual_max()
                    except (FloatingPointError, np.linalg.LinAlgError):
                        break
                    if ((it > 1) and err < criter * abs(self.z[1])
                            and residual <= residual_limit):
                        step_converged = True
                        break
                    if err == 0.0 and residual > residual_limit:
                        break

                if not step_converged:
                    increment *= 0.5
                    if increment < 1.0 / 4096.0 or attempts >= 80:
                        reached = fraction * self.H_target
                        self.last_error = (
                            f"No converged solution at H={self.H_target:g} m, "
                            f"T={self.T_target:g} s, Uc={self.Uc:g} m/s. "
                            f"Height continuation reached H={reached:.6g} m; "
                            f"max|F|={self.residual_max:.3e}. "
                            "The requested wave may be near or beyond the limiting wave height; "
                            "reduce H or increase T."
                        )
                        return
                    continue

                previous_fraction, previous_state = fraction, accepted_state
                fraction, accepted_state = target, self.z.copy()
                self.continuation_fraction = fraction
                increment = min(increment * 1.5, 1.0 / max(1, self.nstep))
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

            # Reimposing exact reporting identities must not invalidate the solution.
            if self._residual_max() > self.residual_tolerance_final:
                self.last_error = ("Post-processing residual check failed: "
                                   f"max|F|={self.residual_max:.3e}.")
                return

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

            self.u_bed_crest, _, _ = self.get_kinematics(0.0, 0.0)
            self.u_surf_crest, _, _ = self.get_kinematics(self.d + self.eta_crest, 0.0)
            u_trough, _, _ = self.get_kinematics(self.d + self.eta_trough, np.pi)
            self.asymmetry = abs(self.u_surf_crest / u_trough) if abs(u_trough) > 0 else 0.0
            self._compute_extrema()
            cf_est = 0.005
            self.tau_bed = 0.5 * RHO * cf_est * self.u_bed ** 2
            orbital_amplitude = self._phase_maximum("bed_orbital")[0]
            self.ExcursionBed = orbital_amplitude * self.T_target / (2.0 * np.pi)

            # The success flag is granted only to a finite, residual-checked result.
            required_outputs = (self.L, self.k, self.c, self.Power, self.Sxx,
                                self.EnergyDensity, self.Cg, self.Impulse,
                                self.PotentialEnergy, self.KineticEnergy,
                                self.MeanSquareBedVelocity, self.MomentumFlux,
                                self.u_bed, self.u_surf, self.acc_max, self.w_max)
            if not np.isfinite(required_outputs).all():
                raise FloatingPointError("Non-finite dimensional output.")
            if self._residual_max() > self.residual_tolerance_final:
                raise FloatingPointError("Final residual exceeds tolerance.")
            self.converged = True

        except FloatingPointError as e:
            self.last_error = f"Floating point failure: {e}"
            self.converged = False
        except Exception as e:
            self.last_error = f"Solver error: {e}"
            self.converged = False
        finally:
            np.seterr(**old_err)

    def _phase_quantity(self, phase, quantity):
        if quantity in ("bed", "bed_orbital"):
            value = self.get_kinematics(0.0, phase)[0]
            if quantity == "bed_orbital":
                value -= self.EulerianCurrent
        else:
            top = self.d * (1.0 + self._surface_keta(phase) / self.z[1])
            values = self.get_kinematics(top, phase)
            value = values[{"surface": 0, "vertical": 1, "acceleration": 2}[quantity]]
        if not np.isfinite(value):
            raise FloatingPointError("Non-finite phase evaluation.")
        return abs(float(value))

    def _phase_maximum(self, quantity):
        """Bracket phase peaks and refine each bracket by golden-section search."""
        phases = np.linspace(0.0, np.pi, max(257, 8 * self.n + 1))
        values = np.array([self._phase_quantity(x, quantity) for x in phases])
        best_index = int(np.argmax(values))
        best, phase = float(values[best_index]), float(phases[best_index])
        golden = (np.sqrt(5.0) - 1.0) / 2.0
        for i in range(1, len(phases) - 1):
            if not (values[i] >= values[i-1] and values[i] >= values[i+1]
                    and (values[i] > values[i-1] or values[i] > values[i+1])):
                continue
            left, right = float(phases[i-1]), float(phases[i+1])
            x1, x2 = right - golden * (right-left), left + golden * (right-left)
            f1, f2 = self._phase_quantity(x1, quantity), self._phase_quantity(x2, quantity)
            for _ in range(64):
                if right - left < 1e-12:
                    break
                if f1 < f2:
                    left, x1, f1 = x1, x2, f2
                    x2 = left + golden * (right-left)
                    f2 = self._phase_quantity(x2, quantity)
                else:
                    right, x2, f2 = x2, x1, f1
                    x1 = right - golden * (right-left)
                    f1 = self._phase_quantity(x1, quantity)
                if f1 > best:
                    best, phase = f1, x1
                if f2 > best:
                    best, phase = f2, x2
        return float(best), float(phase)

    def _compute_extrema(self):
        self.extrema_phases = {}
        for quantity, attribute in (("bed", "u_bed"), ("surface", "u_surf"),
                                    ("vertical", "w_max"), ("acceleration", "acc_max")):
            value, phase = self._phase_maximum(quantity)
            setattr(self, attribute, value)
            self.extrema_phases[quantity] = phase

    def _defining_integrals(self, phase_count, vertical_count):
        """Integrate the dimensional velocity and pressure definitions."""
        phases = 2.0 * np.pi * np.arange(phase_count) / phase_count
        j = self._j
        coefficients = self.Y[1:self.n+1].copy()
        coefficients[-1] *= 0.5
        eta = (np.cos(phases[:, None] * j) @ coefficients) / self.k
        nodes, weights = np.polynomial.legendre.leggauss(vertical_count)
        thickness = (self.d + eta)[:, None] / 2.0
        if np.any(thickness <= 0.0):
            raise FloatingPointError("The interpolated surface intersects the bed.")
        y = -self.d + thickness * (nodes[None, :] + 1.0)
        a = self.k * y[:, :, None] * j
        b = self.z[1] * j
        denominator = 1.0 + np.exp(-2.0*b)
        positive, reflected = np.exp(a), np.exp(-a-2.0*b)
        C, S = (positive+reflected)/denominator, (positive-reflected)/denominator
        amplitude = j * self.B[1:self.n+1] * np.sqrt(self.g/self.k)
        u = self.EulerianCurrent + np.sum(C*amplitude*np.cos(phases[:, None, None]*j), axis=2)
        v = np.sum(S*amplitude*np.sin(phases[:, None, None]*j), axis=2)
        pressure = RHO * (self.Bernoulli_r-self.g*y-0.5*((u-self.c)**2+v*v))
        def average_integral(values):
            return float(np.mean(np.sum(values*weights[None, :]*thickness, axis=1)))
        return {
            "Impulse": average_integral(RHO*u),
            "KineticEnergy": average_integral(0.5*RHO*(u*u+v*v)),
            "PotentialEnergy": float(0.5*RHO*self.g*np.mean(eta*eta)),
            "Sxx": average_integral(pressure+RHO*u*u)-0.5*RHO*self.g*self.d**2,
            "Power": average_integral((pressure+RHO*self.g*y+0.5*RHO*(u*u+v*v))*u),
        }

    def _verify_finite_system(self):
        """Enclose the finite equations at the stored state and selected output formulas."""
        I = _Interval
        z = list(map(I, self.z))
        n = self.n
        H, T, depth, current, gravity, density = map(I, (self.H_target, self.T_target, self.d, self.Uc, self.g, RHO))
        F = [I(0) for _ in range(2*n+10)]
        F[0] = z[2]-z[1]*H/depth
        F[1] = z[2]-H/(gravity*T.square())*z[3].square()
        F[2] = z[4]*z[3]-2*_PI_INTERVAL
        F[3] = z[5]+z[7]-z[4]
        F[4] = z[1]*(z[6]+z[7]-z[4])-z[8]
        F[5] = z[self.Current_criterion+4]-current/(gravity*depth).sqrt()*z[1].sqrt()
        F[6] = z[10]+z[n+10]+2*sum(z[11:n+10], I(0))
        F[7] = z[10]-z[n+10]-z[2]
        table = _collocation_intervals(n)
        denominators = [1+(-2*j*z[1]).exp() for j in range(1, n+1)]
        for m in range(n+1):
            ordinate = z[m+10]
            psi, u, v = I(0), I(0), I(0)
            for j in range(1, n+1):
                ep = (j*ordinate).exp()
                em = (-j*ordinate-2*j*z[1]).exp()
                S, C = (ep-em)/denominators[j-1], (ep+em)/denominators[j-1]
                sine, cosine = table[(m*j) % (2*n)]
                coefficient = z[n+10+j]
                psi += coefficient*S*cosine
                u += j*coefficient*C*cosine
                v += j*coefficient*S*sine
            F[8+m] = psi-z[8]-z[7]*ordinate
            F[n+9+m] = ((u-z[7]).square()+v.square())/2+ordinate-z[9]

        wave_number = z[1]/depth
        speed_scale = (gravity/wave_number).sqrt()
        energy_scale = density*gravity/wave_number.square()
        pulse = z[8]+z[1]*z[5]
        flux_scaled = z[7]*z[1]-z[8]
        kinetic = (z[4]*pulse-z[5]*flux_scaled)/2
        coefficients = list(map(I, self.Y))
        potential = sum((coefficients[j].square()/4 for j in range(1,n)), I(0))+coefficients[n].square()/16
        beta = 2*z[9]-z[4].square()+2*z[4]*z[5]
        stress = 4*kinetic-3*potential+beta*z[1]-2*z[5]*pulse
        power = z[4]*(3*kinetic-2*potential)+beta*(pulse+z[4]*z[1])/2-2*z[4]*z[5]*pulse
        bed_sum = I(0)
        for j in range(1,n+1):
            inverse_cosh = 2*(-j*z[1]).exp()/denominators[j-1]
            bed_sum += (j*z[n+10+j]*inverse_cosh).square()/2
        outputs = {
            "L": 2*_PI_INTERVAL/wave_number,
            "c": z[4]*speed_scale,
            "Impulse": density*speed_scale/wave_number*pulse,
            "KineticEnergy": energy_scale*kinetic,
            "PotentialEnergy": energy_scale*potential,
            "EnergyDensity": energy_scale*(kinetic+potential),
            "Sxx": energy_scale*stress,
            "Power": energy_scale*speed_scale*power,
            "MeanSquareBedVelocity": gravity/wave_number*bed_sum,
        }
        residual_upper = max(value.magnitude() for value in F)
        tolerance = min(Decimal("1e-10"), Decimal.from_float(float(self.residual_tolerance_final)))
        checks = {}
        valid = residual_upper <= tolerance
        for name, interval in outputs.items():
            reported = I(getattr(self, name))
            error = (reported-interval).magnitude()
            minimum = Decimal(0) if interval.lo <= 0 <= interval.hi else min(interval.lo.copy_abs(), interval.hi.copy_abs())
            budget = I.down.add(Decimal("1e-7"), I.down.multiply(Decimal("1e-8"), minimum))
            checks[name] = {"interval": interval.record(), "reported": str(reported.lo),
                            "error_upper": str(error), "budget_lower": str(budget)}
            valid = valid and error <= budget
        self.finite_system_certificate = {
            "arithmetic": "50-digit outward-rounded Decimal intervals",
            "scope": "Residual at the stored finite Fourier state; selected finite-formula output enclosures",
            "root_existence_proved": False,
            "order": n, "equation_count": len(F),
            "residual_intervals": [value.record() for value in F],
            "residual_upper": str(residual_upper), "residual_tolerance": str(tolerance),
            "outputs": checks, "passed": bool(valid),
        }
        self.finite_system_verified = bool(valid)
        return bool(valid)

    def _verify_accuracy(self):
        """Verify numerical budgets before exposing a successful solution."""
        self.accuracy_verified = False
        self.finite_system_verified = False
        self.finite_system_certificate = {}
        failures = []
        metrics = {}
        phases = np.arange(max(2048, 32*self.n)) * (2.0*np.pi/max(2048, 32*self.n))
        j = self._j
        coefficients = self.Y[1:self.n+1].copy()
        coefficients[-1] *= 0.5
        eta = np.cos(phases[:, None]*j) @ coefficients
        a, b = eta[:, None]*j, self.z[1]*j
        denominator = 1.0+np.exp(-2.0*b)
        positive, reflected = np.exp(a), np.exp(-a-2.0*b)
        S, C = (positive-reflected)/denominator, (positive+reflected)/denominator
        cosine, sine = np.cos(phases[:, None]*j), np.sin(phases[:, None]*j)
        B = self.B[1:self.n+1]
        psi = (S*cosine) @ B
        u, v = (C*cosine) @ (j*B), (S*sine) @ (j*B)
        boundary = max(float(np.max(abs(psi-self.z[8]-self.z[7]*eta))),
                       float(np.max(abs(0.5*((u-self.z[7])**2+v*v)+eta-self.z[9]))))
        metrics["surface_residual"] = boundary
        if not np.isfinite(boundary) or boundary > 1e-6:
            failures.append("interpolated-surface residual exceeds 1e-6")

        beta = self.MeanSquareBedAbsoluteVelocity
        reference = {
            "Power": self.c*(3*self.KineticEnergy-2*self.PotentialEnergy)
                     +0.5*beta*(self.Impulse+RHO*self.c*self.d)
                     -2*self.c*self.EulerianCurrent*self.Impulse,
            "Sxx": 4*self.KineticEnergy-3*self.PotentialEnergy+RHO*beta*self.d
                   -2*self.EulerianCurrent*self.Impulse,
        }
        for name, expected in reference.items():
            error = abs(getattr(self, name)-expected)
            budget = 1e-7+1e-8*abs(expected)
            metrics[name+"_spectral_ratio"] = error/budget
            if not np.isfinite(error) or error > budget:
                failures.append(name+" spectral-bed comparison exceeds its budget")

        coarse = self._defining_integrals(max(256, 8*self.n), 64)
        fine = self._defining_integrals(max(512, 16*self.n), 96)
        for name, expected in fine.items():
            budget = 1e-7+1e-8*abs(expected)
            difference = abs(fine[name]-coarse[name])
            error = abs(getattr(self, name)-expected)
            metrics[name+"_quadrature_ratio"] = difference/budget
            metrics[name+"_integral_ratio"] = error/budget
            if not np.isfinite([difference, error]).all() or max(error, difference) > budget:
                failures.append(name+" defining-integral comparison exceeds its budget")
        if not self._verify_finite_system():
            failures.append("finite-system interval verification exceeds its budget")
        metrics["finite_residual_ratio"] = float(_Interval.up.divide(
            Decimal(self.finite_system_certificate["residual_upper"]),
            Decimal(self.finite_system_certificate["residual_tolerance"])))
        metrics["finite_output_ratio"] = max(float(_Interval.up.divide(
            Decimal(check["error_upper"]), Decimal(check["budget_lower"])))
            for check in self.finite_system_certificate["outputs"].values())
        self.accuracy_metrics = metrics
        if failures:
            self.converged = False
            self.solution_status = "resolution_failure"
            self.last_error = ("The collocation system converged, but the available Fourier resolution "
                               "does not satisfy the output accuracy budgets: "+"; ".join(failures)+".")
            return
        self.accuracy_verified = True
        self.solution_status = "verified_numerically"

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
        - Evaluated analytically for the finite Fourier series by orthogonality.
        """
        if self.d <= 0.0 or self.T_target <= 0.0:
            return 0.0

        # Exact phase mean for the finite Fourier series (cosine orthogonality).
        # nph is retained for API compatibility; no phase quadrature is needed.
        kd = float(self.z[1])
        if kd <= 0.0:
            return 0.0
        j = np.arange(1, self.n + 1, dtype=DTYPE)
        e = np.exp(-j * kd)
        sech = 2.0 * e / (1.0 + e * e)
        amplitudes = j * self.B[1:self.n + 1] * sech
        return float(0.5 * self.g * self.d / kd * np.dot(amplitudes, amplitudes))


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
        for i in range(1, self.n):
            pe += 0.25 * (self.Y[i] ** 2)
        pe += (self.Y[self.n] ** 2) / 16.0  # reconstructed amplitude is Y_N/2

        # Fenton (1990), online edition, section 7:
        # beta = <u_bed^2>*k/g is ABSOLUTE, including the Eulerian current.
        # It is different from the separately reported orbital mean square.
        ub2_absolute = float(2.0 * self.z[9] - self.z[4] ** 2
                             + 2.0 * self.z[4] * self.z[5])
        q_term = float(self.z[7] * kd - self.z[8])

        sxx = float(4.0 * ke - 3.0 * pe + ub2_absolute * kd
                    - 2.0 * self.z[5] * pulse)
        f = float(self.z[4] * (3.0 * ke - 2.0 * pe)
                  + 0.5 * ub2_absolute * (pulse + self.z[4] * kd)
                  - 2.0 * self.z[4] * self.z[5] * pulse)

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
        self.EulerianCurrent = float(ce_dimless * np.sqrt(self.g * self.d)) # u1 [m/s]
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
        self.MeanSquareBedAbsoluteVelocity = float(self.MeanSquareBedVelocity
                                                  + self.EulerianCurrent ** 2)

        self.Cg = float(self.Power / self.EnergyDensity) if abs(self.EnergyDensity) > 1e-12 else 0.0


# ==============================================================================
#  GUI APPLICATION CLASS
# ==============================================================================

class FentonApp(tk.Tk):
    def __init__(self):
        super().__init__()
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
        self.log(_runtime_diagnostics() + "\n\n")
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

        valid0 = not s0
        validC = float(Uc_in) != 0.0 and not _solver_status(solverC, "With current")
        if issues and not (valid0 or validC):
            # Close file logging before returning
            if f_log:
                try:
                    f_log.close()
                except Exception:
                    pass
            self._finalize_error("Numerical failure / non-convergence.", "\n".join(issues))
            return

        # A converged state supplies shared input scales while failed columns stay empty.
        if not valid0:
            solver0 = solverC
        if not validC:
            solverC = solver0

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
        box_text("Status: Results shown for converged cases only." if issues
                 else "Status: Finite-system and output accuracy checks passed.")
        for valid, slv, label in ((valid0, solver0, "No current"),
                                  (validC, solverC, "With current")):
            if valid:
                box_text(f"{label}: N={slv.n}; max|F|={slv.residual_max:.3e}")
                certificate = slv.finite_system_certificate
                if certificate.get("passed"):
                    upper = float(certificate["residual_upper"])
                    box_text(f"Finite-system interval residual bound: {upper:.3e}; nine output bounds checked.")
                if slv.numerical_note:
                    box_text(slv.numerical_note)
        hline("-")
        for issue in issues:
            log_wrapper(issue + "\n")
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
            ("Energy-flux velocity (F/E)\u3164", solver0.Cg, wc(solverC.Cg), "m/s"),
        ]

        rows.append({"_section": "KINEMATICS (EXTREMES / BED ORBITAL MOTION)"})
        rows += [
            ("Surface max |u| (estimated)", solver0.u_surf, wc(solverC.u_surf), "m/s"),
            ("Bed max |u| (estimated)", solver0.u_bed, wc(solverC.u_bed), "m/s"),
            ("Surface max |aₓ| (estimated)", solver0.acc_max, wc(solverC.acc_max), "m/s²"),
            ("Velocity asymmetry |uc|/|ut|", solver0.asymmetry, wc(solverC.asymmetry), "-"),
            ("Mean square bed orbital vel ub²", solver0.MeanSquareBedVelocity, wc(solverC.MeanSquareBedVelocity), "m²/s²"),
            ("Bed orbital RMS velocity ub,rms", np.sqrt(max(0.0, solver0.MeanSquareBedVelocity)), wc(np.sqrt(max(0.0, solverC.MeanSquareBedVelocity))), "m/s"),
        ]

        rows.append({"_section": "NONLINEARITY / BREAKING DIAGNOSTICS"})
        warn0 = "BREAKING" if solver0.is_breaking else "BELOW MICHE"
        warnC = ("BREAKING" if solverC.is_breaking else "BELOW MICHE") if has_current else "-"
        rows += [
            ("Miche breaking limit (Hmax)", solver0.breaking_limit_miche, wc(solverC.breaking_limit_miche), "m"),
            ("Saturation (H/Hmax)", solver0.breaking_index, wc(solverC.breaking_index), "-"),
            ("Breaking status", warn0, warnC, "-"),
            ("Ursell number (U)", solver0.ursell, wc(solverC.ursell), "-"),
            ("Regime (by d/L)", solver0.regime, wc(solverC.regime), "-"),
        ]

        rows = [row if isinstance(row, dict) else
                (row[0], row[1] if valid0 else None, row[2] if validC else None, row[3])
                for row in rows]
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
                (17, "Mean square of bed velocity (absolute)", f"{slv.MeanSquareBedAbsoluteVelocity:.5f}","ub²/gd", f"{(slv.MeanSquareBedAbsoluteVelocity/(g*d)):.5f}"),
                (18, "Radiation stress",                     f"{kN(slv.Sxx):.5f}",             "S_xx/ρgd²",                       f"{(slv.Sxx_depth):.5f}"),
                (19, "Wave power",                           f"{kW(slv.Power):.5f}",           "F/(ρg³ᐟ²d⁵ᐟ²)\u3164\u3164",           f"{(slv.F_depth):.5f}"),
            ]

            # Convert rows to generic printer rows (already formatted strings for numeric cols).
            out_rows = [(a, b, c, d_, e) for (a, b, c, d_, e) in rows]
            print_table(title, headers, col_w, aligns, out_rows)

        if valid0:
            print_solution_flat(solver0, "SOLUTION.RES (NO CURRENT)")
        if validC:
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
            ("Cg", "Energy-flux / energy-density ratio F/E.", "m/s"),
            ("Sₓₓ", "Radiation stress component in wave direction.", "kN/m ; Sₓₓ/(ρgd²)"),
            ("ub²", "Mean square *orbital* bed velocity: ub² = <(ub(t) − ū₁)²>. Non-negative by definition; computed by phase averaging.", "m²/s² ; /gd"),
            ("ub,rms", "Root-mean-square orbital bed velocity: ub,rms = √(ub²).", "m/s"),
            ("usurf,max", "Estimated maximum absolute surface horizontal velocity over phase.", "m/s"),
            ("ubed,max", "Estimated maximum absolute bed horizontal velocity over phase.", "m/s"),
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
def _runtime_diagnostics():
    """Identify the actual running file, including packaged Windows programs."""
    import hashlib
    path = os.path.abspath(sys.executable if getattr(sys, "frozen", False) else __file__)
    try:
        with open(path, "rb") as source:
            digest = hashlib.file_digest(source, "sha256").hexdigest() if hasattr(hashlib, "file_digest") else hashlib.sha256(source.read()).hexdigest()
    except OSError:
        digest = "unavailable"
    return (f"Running file: {path}\nSHA256: {digest}\n"
            f"Python: {sys.version.split()[0]}; NumPy: {np.__version__}; "
            f"Numba active: {NUMBA_AVAILABLE}")


def _self_test():
    """Exercise default and short-period waves without opening a window."""
    global NUMBA_AVAILABLE
    print(_runtime_diagnostics(), flush=True)
    original_backend = NUMBA_AVAILABLE
    backends = [False, True] if original_backend else [False]
    all_ok = True
    cases = [(3.0, 9.0, 5.0, u) for u in (0.0, 1.0)]
    cases += [(h, t, 5.0, u) for h in (1.0, 2.0)
              for t in (3.0, 4.0) for u in (0.0, 1.0)]
    cases.append((3.0, 4.0, 5.0, 1.0))
    resolution_limited = {(2.0, 3.0, 5.0, 0.0), (3.0, 4.0, 5.0, 1.0)}
    try:
        for use_numba in backends:
            NUMBA_AVAILABLE = use_numba
            for h, t, d, current in cases:
                solver = FentonStreamFunction(h, t, d, current)
                solver.solve()
                if (h, t, d, current) in resolution_limited:
                    ok = (not solver.converged and not solver.accuracy_verified
                          and solver.solution_status == "resolution_failure"
                          and bool(solver.last_error))
                else:
                    ok = (solver.converged and solver.accuracy_verified
                          and np.isfinite(solver.residual_max)
                          and solver.residual_max <= solver.residual_tolerance_final)
                all_ok = all_ok and bool(ok)
                print(f"{'Numba' if use_numba else 'NumPy'}: H={h:g}, T={t:g}, "
                      f"d={d:g}, Uc={current:g}; N={solver.n}; "
                      f"{'PASS' if ok else 'FAIL'}; {solver.solution_status}; "
                      f"max|F|={solver.residual_max:.6e}", flush=True)
                if not ok:
                    print(solver.last_error, flush=True)
    finally:
        NUMBA_AVAILABLE = original_backend
    print("ALL SOLVER TESTS PASSED" if all_ok else "SOLVER TEST FAILED", flush=True)
    return 0 if all_ok else 1


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Fenton Wave Solver ")
    parser.add_argument("--self-test", action="store_true", help="test default and short-period waves without opening a window")
    args = parser.parse_args()
    if args.self_test:
        raise SystemExit(_self_test())
    app = FentonApp()
    app.mainloop()
