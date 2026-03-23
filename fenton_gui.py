# ==============================================================================
#  ENGINEERING TECHNICAL REFERENCE & THEORETICAL FORMULATION
# ==============================================================================
#  PROGRAM:      Nonlinear Wave Hydrodynamics Solver (Fenton's Stream Function)
#  METHOD:       Fourier Approximation Method for Steady Water Waves (N=50)
#  REFERENCE:    Fenton, J.D. (1999). "Numerical methods for nonlinear waves."
#                In P.L.-F. Liu (Ed.), Advances in Coastal and Ocean Engineering
#                (Vol. 5, pp. 241–324). World Scientific: Singapore.
# ==============================================================================
#
#  1. INTRODUCTION & SCOPE
#  -----------------------------------------------------------------------------
#  This software calculates the hydrodynamics of steady, periodic surface gravity
#  waves using high-order Stream Function theory. Unlike Linear (Airy) Theory,
#  which assumes infinitesimal amplitudes, this method retains full nonlinearity
#  in the boundary conditions.
#
#  Implementation Specifics (C++-Parity Port):
#  - Solver: Newton iteration with finite-difference Jacobian and SVD linear solve
#            (Press et al. truncation wmin = wmax*1e-12), matching the provided C++.
#  - Stability: Continuation in wave height (nstep=4) with linear-theory initialisation
#            and extrapolation between steps, matching Fourier.cpp.
#  - Regime: Finite-depth, Period-input formulation as in the original program.
#
#  2. GOVERNING FIELD EQUATIONS
#  -----------------------------------------------------------------------------
#  The fluid is modeled as inviscid, incompressible, and irrotational.
#  The flow is solved in a frame of reference moving with the wave celerity (c),
#  rendering the flow steady.
#
#  A. Field Equation (Laplace):
#     ∇²ψ = ∂²ψ/∂x² + ∂²ψ/∂z² = 0
#     Where ψ(x,z) is the stream function. Velocities are defined as:
#     u =  ∂ψ/∂z   (Horizontal)
#     w = -∂ψ/∂x   (Vertical)
#
#  B. Bottom Boundary Condition (BBC) at z=0:
#     The seabed is impermeable (a streamline).
#     ψ(x, 0) = -Q
#     Where Q is the volume flux per unit width in the moving frame.
#
#  3. FREE SURFACE BOUNDARY CONDITIONS
#  -----------------------------------------------------------------------------
#  The solution is constrained by two nonlinear conditions at the unknown
#  free surface elevation z = η(x):
#
#  A. Kinematic Boundary Condition (KBC):
#     The free surface is a streamline (constant ψ).
#     ψ(x, η) = 0
#
#  B. Dynamic Boundary Condition (DBC - Bernoulli):
#     Pressure is constant (atmospheric) along the surface.
#     1/2 * [ (∂ψ/∂x)² + (∂ψ/∂z)² ] + gη = R
#     Where R is the Bernoulli constant (Total Energy Head).
#
#  4. NUMERICAL SOLUTION (FOURIER ANSATZ)
#  -----------------------------------------------------------------------------
#  The stream function is approximated by a truncated Fourier series of order N
#  (N=50) that analytically satisfies the Field Equation and Bottom BC:
#
#    ψ(x,z) = -(ū + c) z + Σ_{j=1..N} B_j * [sinh(jkz)/cosh(jkd)] * cos(jkx)
#
#  Deep Water Numerical Stability:
#  To prevent floating-point overflow when kd >> 1, the code replaces the
#  hyperbolic ratio with asymptotic exponentials when arguments > 25.0:
#
#    sinh(jkz)/cosh(jkd) ≈ exp(jk(z-d))
#
#  Optimization Vector (State Space):
#  The solver minimizes residuals for the vector
#
#    X = [k, η_0...η_N, B_1...B_N, Q, R].
#
#  IMPORTANT (Overdetermined Residual System):
#  The residual vector dimension is NOT equal to the number of unknowns.
#  For N=50:
#    - Unknowns: n = 1 + (N+1) + N + 2 = 2N + 4 = 104
#    - Residuals: m = 3 + (N+1) + (N+1) = 2(N+1) + 3 = 105
#  Any C++ (or other) port must allocate m=105 and never assume m==n; otherwise
#  out-of-bounds writes can occur and results may become optimizer/flags dependent.
#
#  5. DERIVED PHYSICAL PARAMETERS & OUTPUT DEFINITIONS
#  -----------------------------------------------------------------------------
#  Upon convergence, the software calculates the following engineering parameters
#  derived from the solved Fourier coefficients (B_j).
#
#  A. FUNDAMENTAL WAVE GEOMETRY & PHASE
#  ------------------------------------
#  1. Wavelength (L):
#     Horizontal distance between crests. Solved via dispersion relation.
#     L = c·T = 2π / k
#
#  2. Celerity (c):
#     Phase velocity. c = L / T.
#
#  B. KINEMATICS (VELOCITIES & ACCELERATIONS)
#  ------------------------------------------
#  1. Horizontal Velocity (u):
#     u(x,z) = c - ū + Σ_{j=1..N} jkB_j * [cosh(jkz)/cosh(jkd)] * cos(jkx)
#
#  2. Vertical Velocity (w):
#     w(x,z) = Σ_{j=1..N} jkB_j * [sinh(jkz)/cosh(jkd)] * sin(jkx)
#
#  3. Max Acceleration (a_x):
#     Total derivative (Convective acceleration).
#     a_x = Du/Dt = u * ∂u/∂x + w * ∂u/∂z
#
#     NOTE (Python-parity detail used by the C++ port):
#     The vertical perturbation term uses +sin(j·phase) (not -sin). A sign error
#     distorts the convective term w·∂u/∂z and breaks Max Accel parity.
#
#  4. Velocity Asymmetry:
#     Asymmetry = |u_crest| / |u_trough|
#
#  C. DYNAMICS (INTEGRAL PROPERTIES)
#  ---------------------------------
#  Computed using exact integral invariants (Fenton Eqs 14-16).
#
#  1. Impulse (I):
#     Total wave momentum (kg·m/s).
#     I = ρ(c d - Q)
#
#  2. Energy Density (E):
#     Mean Energy (J/m²).
#     PE = 1/2 ρ g mean(η²)
#     KE = 1/2 (cI - Qρ U_c)
#     E  = PE + KE
#
#  3. Power / Energy Flux (P):
#     Rate of energy transfer (W/m).
#     P = c(3KE - 2PE) + 1/2 mean(u_b²)(I + ρ c d) + 1/2 ρ Q U_c²
#
#     Note on mean(u_b²) (Mean Square Bed Velocity):
#     To avoid deep-water integration errors, this is computed algebraically:
#     mean(u_b²) = 2(R - g d) - c²
#
#  4. Radiation Stress (Sxx):
#     Excess momentum flux (N/m).
#     Sxx = 4KE - 3PE + ρ mean(u_b²) d + 2ρ I U_c
#
#  5. Stokes current (ū₂) (U_drift):
#     U_drift = I / (ρ d)
#
#  D. STABILITY & REGIME CLASSIFICATION
#  ------------------------------------
#  1. Ursell Number (U_r):
#     U_r = H L² / d³ (Values > 26 indicate significant nonlinearity).
#
#  2. Miche Limit (H_max):
#     Theoretical max height before breaking.
#     H_max = 0.142 L tanh(kd)
#
#  3. Saturation (Breaking Index):
#     Saturation = H / H_max
#     - If > 1.0: Wave is BREAKING.
#     - If < 1.0: Wave is STABLE.
#
#  4. Regime:
#     - Shallow:      d/L < 0.05
#     - Intermediate: 0.05 < d/L < 0.5
#     - Deep:         d/L > 0.5
#
# ==============================================================================
#  6. SOFTWARE USAGE & COMPILATION GUIDE
# ==============================================================================
#  
#  A. PREREQUISITES
#  ----------------
#  - Python 3.8 or higher installed on your system.
#  - Basic familiarity with the command line (Terminal/CMD).
#
#  B. RUNNING FROM SOURCE (DEVELOPMENT)
#  ------------------------------------
#  1. Create a clean virtual environment (CRITICAL to avoid library conflicts):
#     > python -m venv build_env
#
#  2. Activate the environment:
#     - Windows: > build_env\Scripts\activate
#     - Linux/Mac: > source build_env/bin/activate
#
#  3. Install the required dependencies (keep the environment clean):
#     (build_env)> python -m pip install --upgrade pip
#     (build_env)> python -m pip install -U "numpy>=1.20"    # NumPy 2.x supported
#
#     Optional but strongly recommended (major speed-up via JIT compilation):
#     (build_env)> python -m pip install -U "numba>=0.57"       # if this fails, skip (NumPy fallback)
#
#     Notes:
#     - If Numba is not installed (or not supported on your OS/Python build),
#       the solver automatically falls back to the pure-NumPy implementation.
#     - On Windows, installing into a virtual environment avoids DLL conflicts
#       with other Python distributions on the system.
#
#  4. Run the GUI script:
#     (build_env)> python fenton_gui.py
#
#  C. COMPILING TO STANDALONE EXECUTABLE (.EXE)
#  --------------------------------------------
#  To create a portable .exe file that runs without Python:
#
#  1. Activate your clean virtual environment (as above).
#
#  2. Install PyInstaller (keep the environment clean!):
#     (build_env)> pip install pyinstaller
#
#  3. Build the executable:
#     (build_env)> pyinstaller --onefile --noconsole fenton_gui.py
#  
#     Flags explanation:
#     --onefile    : Bundles everything into a single .exe file (e.g., 50MB).
#     --noconsole  : Hides the black command window (GUI mode only).
#
#  4. Locate the output:
#     The finished 'fenton_gui.exe' will be in the 'dist' folder.
#
#  5. Cleanup:
#     (build_env)> deactivate
#
# ==============================================================================
#  BIBLIOGRAPHY
# ==============================================================================
#
#  1.  Fenton, J.D. (1999). "Numerical methods for nonlinear waves." 
#      In P.L.-F. Liu (Ed.), Advances in Coastal and Ocean Engineering (Vol. 5, 
#      pp. 241–324). World Scientific: Singapore.
#      [Primary Source: Comprehensive review of fully-nonlinear methods including 
#      Fourier approximation, Boundary Integral Equation (BIE) methods, and 
#      Local Polynomial Approximation].
#      URL: https://johndfenton.com/Papers/Fenton99Liu-Numerical-methods-for-nonlinear-waves.pdf
#
#  2.  Fenton, J.D. (1988). "The numerical solution of steady water wave problems."
#      Computers & Geosciences, 14(3), 357–368.
#      [The core algorithm for high-accuracy Stream Function Theory].
#      URL: https://doi.org/10.1016/0098-3004(88)90066-0
#
#  3.  Fenton, J.D. (1985). "A fifth-order Stokes theory for steady waves."
#      Journal of Waterway, Port, Coastal, and Ocean Engineering, 111(2), 216–234.
#      [Standard analytical theory for deep/intermediate water pile design].
#      URL: https://doi.org/10.1061/(ASCE)0733-950X(1985)111:2(216)
#
#  4.  Fenton, J.D. (1978). "Wave forces on vertical bodies of revolution."
#      Journal of Fluid Mechanics, 85(2), 241–255.
#      [Foundational diffraction theory for large diameter piles].
#      URL: https://johndfenton.com/Papers/Fenton78-Waves-on-bodies-of-revolution.pdf
#
#  5.  Fenton, J.D. (1990). "Nonlinear wave theories." In B. Le Méhauté & 
#      D.M. Hanes (Eds.), The Sea: Ocean Engineering Science (Vol. 9, Part A).
#      John Wiley & Sons.
#      [Comprehensive guide for selecting wave theories: Stokes vs Cnoidal vs Stream].
#      URL: https://www.johndfenton.com/Papers/Fenton90b-Nonlinear-wave-theories.pdf
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
#  NUMBA ACCELERATION LAYER (OPTIONAL, ZERO-ALGORITHM-CHANGE)
# ==============================================================================
#
# The solver is dominated by repeated evaluations of the nonlinear residual system
# Eqns() inside a finite-difference Jacobian (≈ O(num^2) residual calls per Newton
# step). Pure NumPy vectorisation creates many temporary arrays per collocation
# node; Numba JIT removes those temporaries by compiling tight scalar loops.
#
# IMPORTANT:
# - fastmath is NOT used (to avoid changing floating-point semantics).
# - These kernels implement the exact same algebra as the Python reference path.
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
#  GLOBAL CONSTANTS & CONFIGURATION
# ==============================================================================

# Physical Constants (Matched exactly to C++ Phys namespace)
G_STD = 9.80665         # Standard Gravity [m/s^2]
RHO   = 1025.0          # Density of Seawater [kg/m^3]

# Numerical Configuration
DTYPE = np.float64      # Precision for floating point arithmetic
N_FOURIER = 50          # Order of Fourier Series (N). 50 ensures convergence 
                        # even for highly nonlinear near-breaking waves.
N_NUMBERS = 8           # formatting precision for output text

# Suppress optimization warnings (e.g., initial Jacobian singular matrix)
# that occur normally during the first iterations of the solver.
warnings.filterwarnings('ignore')

@njit(cache=True)
def _eqns_numba(z, rhs, coeff, Tanh, cos_nm, sin_nm, n, num, Hoverd, height, Current, Current_criterion):
    """Numba implementation of Eqns() (finite depth, Period case) with 1-based indexing."""
    pi = np.pi

    # Eqn 1
    rhs[1] = z[2] - z[1] * Hoverd

    # Eqn 2 (Period case)
    rhs[2] = z[2] - height * z[3] * z[3]

    # Eqn 3
    rhs[3] = z[4] * z[3] - 2.0 * pi

    # Eqn 4
    rhs[4] = z[5] + z[7] - z[4]

    # Eqn 5
    rhs[5] = z[1] * (z[6] + z[7] - z[4]) - z[8]

    # coeff and tanh tables
    kd = z[1]
    for i in range(1, n + 1):
        coeff[i] = z[n + i + 10]
        Tanh[i] = np.tanh(i * kd)

    # Eqn 6 (finite depth; correction uses sqrt(z[1]))
    rhs[6] = z[Current_criterion + 4] - Current * np.sqrt(kd)

    # Eqn 7 (mean free surface level; scaling constant irrelevant)
    rhs7 = z[10] + z[n + 10]
    for i in range(1, n):
        rhs7 += 2.0 * z[10 + i]
    rhs[7] = rhs7

    # Eqn 8 (wave height definition)
    rhs[8] = z[10] - z[n + 10] - z[2]

    # Eqns 9..(n+9) and (n+10)..(2n+10): free-surface BCs
    for m in range(0, n + 1):
        zsurf = z[10 + m]  # k(eta-d) at this node

        psi = 0.0
        u = 0.0
        v = 0.0

        for jj in range(1, n + 1):
            cj = coeff[jj]
            tj = Tanh[jj]

            x = jj * zsurf
            # Prevent overflow if the iteration diverges (numerical safeguard)
            if x > 60.0 or x < -60.0:
                rhs[1] = np.inf
                return np.inf
            e = np.exp(x)
            inv_e = 1.0 / e
            sinhkd = 0.5 * (e - inv_e)
            coshkd = 0.5 * (e + inv_e)

            # Hyperbolic rewrite (A-8): C = cosh + sinh*tanh(jkd), S = sinh + cosh*tanh(jkd)
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

    # Sum of squares
    ss = 0.0
    for i in range(1, num + 1):
        ss += rhs[i] * rhs[i]
    return ss


@njit(cache=True)
def _compute_Y_and_B_numba(z, B, Y, cosa, n):
    """Numba implementation of the post-convergence Fourier transform block."""
    for i in range(0, Y.shape[0]):
        Y[i] = 0.0

    two_n = 2 * n

    for j in range(1, n + 1):
        B[j] = z[j + n + 10]

        # s = 0.5 * (z[10] + z[n+10] * ((-1.0)**j))
        sign = 1.0
        if (j % 2) == 1:
            sign = -1.0
        s = 0.5 * (z[10] + z[n + 10] * sign)

        for m in range(1, n):
            s += z[10 + m] * cosa[(m * j) % two_n]

        Y[j] = 2.0 * s / n


@njit(cache=True)
def _surface_keta_numba(Y, n, X):
    """Numba implementation of Surface(X): returns k(eta-d) at phase X."""
    kEta = 0.0
    for j in range(1, n):
        kEta += Y[j] * np.cos(j * X)
    kEta += 0.5 * Y[n] * np.cos(n * X)
    return kEta


@njit(cache=True)
def _point_numba(X, Y, z, Tanh, B, n):
    """Numba implementation of the finite-depth Point(X,Y) kernel."""
    kd = z[1]

    # depth-scaled dimensionless bulk values
    c  = z[4] / np.sqrt(kd)
    ce = z[5] / np.sqrt(kd)
    R  = 1.0 + z[9] / kd

    # local variables in wave scaling
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

    # convert to depth scaling (see C++ comments)
    inv_kd_sqrt = 1.0 / np.sqrt(kd)
    inv_kd_32   = 1.0 / (kd ** 1.5)

    u  *= inv_kd_sqrt
    v  *= inv_kd_sqrt
    ux *= np.sqrt(kd)
    vx *= np.sqrt(kd)

    # add Euler current to u
    u = ce + u

    # time derivatives (steady in moving frame)
    ut = -c * ux
    vt = -c * vx
    uy = vx
    vy = -ux

    dudt = ut + u * ux + v * uy

    return float(u), float(v), float(dudt)


# ==============================================================================
#  CORE SOLVER CLASS (C++-PARITY PORT)
# ==============================================================================
#
# This class implements a direct 1:1 port of the original C++ Fourier/Stream-
# Function solver (Subroutines.cpp/Fourier.cpp) for the finite-depth, Period
# case, including:
#   - Unknown vector z[1..2N+10] (same ordering as C++)
#   - Equation system Eqns() (same algebra, same hyperbolic rewrite)
#   - Newton() with finite-difference Jacobian and SVD linear solve
#   - Height stepping (continuation) identical to Fourier.cpp
#
# Notes
# -----
# - The GUI (Tkinter) is intentionally kept unchanged; only the solver core and
#   derived-quantity post-processing are rewritten for parity.
# - This implementation keeps the original 1-based indexing convention for z[]
#   and associated arrays to reduce porting risk.
#

class FentonStreamFunction:
    """
    Fenton steady-wave solver using the Fourier approximation / stream-function
    method, implemented as a direct port of the provided C++ reference code.

    Public API is preserved to keep the GUI unchanged.
    """

    # --------------------------- construction ---------------------------------

    def __init__(self, H, T, d, Uc=0.0):
        # Inputs (physical)
        self.H_target = float(H)    # [m]
        self.T_target = float(T)    # [s]
        self.d        = float(d)    # [m]
        self.Uc       = float(Uc)   # [m/s] (Eulerian / lab-frame)

        # Constants
        self.g = G_STD
        self.N = N_FOURIER

        # Solver control (C++-style defaults)
        self.nstep  = 4        # continuation steps in wave height
        self.number = 40       # max Newton iterations per step
        self.crit   = 1.0e-8   # intermediate-step convergence factor (C++: crit)
        self.criter_final = 1.0e-10  # final-step convergence factor

        # Problem mode (matches GUI: finite depth, Period input)
        self.Depth = "Finite"
        self.Case  = "Period"

        # Current criterion: 1=Eulerian, 2=Stokes (GUI input is Eulerian current)
        self.Current_criterion = 1

        # Derived input non-dimensional groups (C++ Read_data equivalents)
        # MaxH == H/d;  T_nd == T*sqrt(g/d); Height == (H/d)/(T_nd^2) == H/(g T^2)
        self.MaxH = self.H_target / self.d if self.d > 0 else 0.0
        self.T_nd = self.T_target * np.sqrt(self.g / self.d) if self.d > 0 else 0.0
        self.Height = (self.MaxH / (self.T_nd * self.T_nd)) if self.T_nd > 0 else 0.0

        # Current input in C++ is dimensionless w.r.t. sqrt(g d) for finite depth
        self.Current = self.Uc / np.sqrt(self.g * self.d) if self.d > 0 else 0.0

        # ------------------------- outputs (public) ----------------------------
        self.k = 0.0
        self.L = 0.0
        self.c = 0.0
        self.converged = False

        # Human-readable failure reason (used by GUI when convergence fails)
        self.last_error = ""

        # Robustness for large ambient currents: allow more continuation steps and
        # Newton iterations. This does NOT change equations; it only increases the
        # solver budget to avoid premature failure when |Uc| is high.
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

        # ------------------------- internal C++ arrays -------------------------
        self.n = int(self.N)
        self.num = 2 * self.n + 10

        # 1-based vectors (index 0 unused)
        self.z    = np.zeros(self.num + 1, dtype=DTYPE)
        self.rhs1 = np.zeros(self.num + 1, dtype=DTYPE)
        self.rhs2 = np.zeros(self.num + 1, dtype=DTYPE)
        self.coeff = np.zeros(self.n + 1, dtype=DTYPE)    # coeff[1..n]
        self.Tanh  = np.zeros(self.n + 1, dtype=DTYPE)    # Tanh[1..n]
        self.B     = np.zeros(self.n + 1, dtype=DTYPE)    # B[1..n]
        self.Y     = np.zeros(self.num + 1, dtype=DTYPE)  # Y[0..n] used; keep size

        # Precomputed trig tables as in init()
        self.cosa = np.zeros(2 * self.n + 1, dtype=DTYPE)  # [0..2n]
        self.sina = np.zeros(2 * self.n + 1, dtype=DTYPE)

        # Precompute constant trig tables and collocation lookup tables (C++ init())
        k_idx = np.arange(0, 2 * self.n + 1, dtype=DTYPE)
        self.cosa[:] = np.cos(k_idx * np.pi / self.n)
        self.sina[:] = np.sin(k_idx * np.pi / self.n)

        self._j = np.arange(1, self.n + 1, dtype=DTYPE)
        self._j_int = np.arange(1, self.n + 1, dtype=np.int64)
        self._nm_map = (np.arange(0, self.n + 1, dtype=np.int64)[:, None] * self._j_int[None, :]) % (2 * self.n)
        self._cos_nm = self.cosa[self._nm_map]  # shape (n+1, n)
        self._sin_nm = self.sina[self._nm_map]

        # Extrapolation storage sol[i][1..2]
        self.sol = np.zeros((self.num + 1, 3), dtype=DTYPE)

        # Run-time step variables (C++ globals)
        self.height = 0.0    # stepped 'height' (dimensionless)
        self.Hoverd = 0.0    # stepped H/d

    # --------------------------- C++ port helpers -----------------------------

    def _init_linear(self):
        """
        Port of C++ init() for finite-depth, Period case (with current criterion).
        Produces an initial state in z[1..num] for the first height step.
        """
        n = self.n
        pi = np.pi

        # For finite depth
        sigma = 2.0 * pi * np.sqrt(self.height / self.Hoverd) if self.Hoverd > 0 else 0.0

        # Fenton & McKee (1990) approximation used in the C++ (commented alternatives omitted)
        if sigma > 0:
            self.z[1] = (sigma * sigma) / (np.tanh(sigma ** 1.5) ** (2.0 / 3.0))
        else:
            # very small waves / degenerate: start with something benign
            self.z[1] = 2.0 * pi * max(self.height, 1e-12) / max(self.Hoverd, 1e-12)

        self.z[2] = self.z[1] * self.Hoverd
        self.z[4] = np.sqrt(np.tanh(self.z[1]))
        self.z[3] = 2.0 * pi / self.z[4]

        # Current initialisation (finite)
        if self.Current_criterion == 1:
            self.z[5] = self.Current * np.sqrt(self.z[2])
            self.z[6] = 0.0
        else:
            self.z[6] = self.Current * np.sqrt(self.z[2])
            self.z[5] = 0.0

        self.z[7] = self.z[4]
        self.z[8] = 0.0
        self.z[9] = 0.5 * self.z[7] * self.z[7]

        # Initial surface elevation nodes and Fourier coefficients (B_j)
        self.z[10] = 0.5 * self.z[2]
        for i in range(1, n + 1):
            self.z[n + i + 10] = 0.0
            self.z[i + 10] = 0.5 * self.z[2] * self.cosa[i]

        self.z[n + 11] = 0.5 * self.z[2] / self.z[7]

        # store sol[] for extrapolation (C++ sets sol[10..] to zero for very first)
        for i in range(1, 10):
            self.sol[i, 1] = self.z[i]
        for i in range(10, self.num + 1):
            self.sol[i, 1] = 0.0

    def _eqns(self, rhs_out):
        """
        Port of C++ Eqns(double *rhs). Fills rhs_out[1..num] and returns sum(rhs^2).
        Finite-depth branch only (GUI mode).
        """
        # Numba-accelerated kernel (drops back to pure NumPy if Numba is unavailable).
        #
        # Robustness note:
        # If the JIT path produces NaN/Inf (usually due to a diverging Newton iterate),
        # fall back to the pure-NumPy path below (same algebra, clearer exceptions).
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

        # Eqn 1
        rhs[1] = z[2] - z[1] * self.Hoverd

        # Eqn 2 (Period case)
        rhs[2] = z[2] - self.height * z[3] * z[3]

        # Eqn 3
        rhs[3] = z[4] * z[3] - 2.0 * pi

        # Eqn 4
        rhs[4] = z[5] + z[7] - z[4]

        # Eqn 5
        rhs[5] = z[1] * (z[6] + z[7] - z[4]) - z[8]

        # coeff and tanh tables
        for i in range(1, n + 1):
            self.coeff[i] = z[n + i + 10]
            self.Tanh[i] = np.tanh(i * z[1])

        # Eqn 6 (finite depth; correction uses sqrt(z[1]))
        rhs[6] = z[self.Current_criterion + 4] - self.Current * np.sqrt(z[1])

        # Eqn 7 (mean free surface level; scaling constant irrelevant)
        rhs[7] = z[10] + z[n + 10]
        for i in range(1, n):
            rhs[7] += 2.0 * z[10 + i]

        # Eqn 8 (wave height definition)
        rhs[8] = z[10] - z[n + 10] - z[2]

        # Eqns 9..(n+9) and (n+10)..(2n+10): free-surface BCs
        j = self._j                      # shape (n,)
        coeff = self.coeff[1:n + 1]      # shape (n,)
        tanh = self.Tanh[1:n + 1]        # shape (n,)
        jcoeff = j * coeff               # shape (n,)

        for m in range(0, n + 1):
            zsurf = z[10 + m]  # k(eta-d) at this node

            x = j * zsurf
            if np.any(x > 60.0) or np.any(x < -60.0):
                raise FloatingPointError("Divergence: exp(j*zsurf) out of safe range.")
            e = np.exp(x)
            inv_e = 1.0 / e
            sinhkd = 0.5 * (e - inv_e)
            coshkd = 0.5 * (e + inv_e)

            # Hyperbolic rewrite (A-8): C = cosh + sinh*tanh(jkd), S = sinh + cosh*tanh(jkd)
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
        Port of the "slow Fourier transform" block in Fourier.cpp after convergence.
        Produces B[1..n] and Y[0..n] from final z[].
        """
        # Numba-accelerated kernel (drops back to pure NumPy if Numba is unavailable).
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
        Port of C++ Surface(double x): returns k(eta-d) at phase X (0..pi).
        """
        # Numba-accelerated kernel (drops back to pure NumPy if Numba is unavailable).
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
        Port of the finite-depth branch of C++ Point(X,Y). Returns:
          u_dimless (w.r.t sqrt(g d))
          v_dimless (w.r.t sqrt(g d))
          dudt_dimless (w.r.t g)
        Input:
          X : phase in radians (k x)
          Y : vertical coordinate in wave scaling (k(z-d))
        """
        # Numba-accelerated kernel (drops back to pure NumPy if Numba is unavailable).
        if NUMBA_AVAILABLE:
            return _point_numba(float(X), float(Y), self.z, self.Tanh, self.B, self.n)

        n = self.n
        kd = float(self.z[1])

        # depth-scaled dimensionless bulk values
        c  = float(self.z[4] / np.sqrt(kd))
        ce = float(self.z[5] / np.sqrt(kd))
        R  = float(1.0 + self.z[9] / kd)

        # local variables in wave scaling
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

        # convert to depth scaling (see C++ comments)
        inv_kd_sqrt = 1.0 / np.sqrt(kd)
        inv_kd_32   = 1.0 / (kd ** 1.5)

        u  *= inv_kd_sqrt
        v  *= inv_kd_sqrt
        ux *= np.sqrt(kd)
        vx *= np.sqrt(kd)

        # add Euler current to u
        u = ce + u

        # time derivatives (steady in moving frame)
        ut = -c * ux
        vt = -c * vx
        uy = vx
        vy = -ux

        dudt = ut + u * ux + v * uy
        # dvdt is available if needed:
        # dvdt = vt + u * vx + v * vy

        # Bernoulli/pressure are not required by GUI, but kept for completeness:
        # y = 1.0 + Y / kd
        # Pressure = R - y - 0.5 * (((u - c) ** 2) + v * v)

        return float(u), float(v), float(dudt)

    # ----------------------------- public methods -----------------------------

    def solve(self):
        """
        Solve the steady nonlinear wave problem using Fenton's Fourier / stream-function method.

        The implementation follows the reference C++ structure:
        - Continuation in wave height (nstep)
        - Newton iterations with finite-difference Jacobian
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
                    # IMPORTANT: Update continuation storage BEFORE the convergence break.
                    # Otherwise, sol[:,2] may remain at a previous iterate, and the next
                    # continuation-step extrapolation can start from an invalid state.
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

                # update Y and B for this step (C++ does this each step)
                self._compute_Y_and_B()

            # ------------------------- dimensional post-process --------------------

            kd = float(self.z[1])
            if (not np.isfinite(kd)) or (kd <= 0.0):
                raise FloatingPointError("Invalid wavenumber (kd).")

            k_phys = kd / self.d
            L_phys = 2.0 * np.pi / k_phys
            c_dimless = float(self.z[4] / np.sqrt(kd))   # c / sqrt(g d)
            c_phys = c_dimless * np.sqrt(self.g * self.d)

            if (not np.isfinite(L_phys)) or (L_phys <= 0.0):
                raise FloatingPointError("Invalid wavelength.")
            if not np.isfinite(c_phys):
                raise FloatingPointError("Invalid celerity.")

            # surface nodes correspond to m*pi/n (half wave: crest->trough)
            eta_nodes = np.zeros(self.n + 1, dtype=DTYPE)
            for m in range(0, self.n + 1):
                kEta = float(self.z[10 + m])             # k(eta-d) at node
                eta_nodes[m] = self.d * (1.0 + kEta / kd)

            self.eta_nodes = eta_nodes
            self.k = float(k_phys)
            self.L = float(L_phys)
            self.c = float(c_phys)

            # store Bj as 0-based array for external use
            self.Bj = self.B[1:self.n + 1].copy()

            # crest/trough elevations relative to SWL
            self.eta_crest = float(self.eta_nodes[0] - self.d)
            self.eta_trough = float(self.eta_nodes[-1] - self.d)

            # non-dimensional descriptors
            self.steepness = self.H_target / self.L
            self.rel_depth = self.d / self.L
            self.ursell = (self.H_target * self.L * self.L) / (self.d ** 3)

            # regime classification (engineering convenience)
            if self.rel_depth < 0.05:
                self.regime = "Shallow"
            elif self.rel_depth < 0.5:
                self.regime = "Intermediate"
            else:
                self.regime = "Deep"

            # Miche breaking limit (keep GUI label/behaviour)
            self.breaking_limit_miche = float(0.142 * self.L * np.tanh(self.k * self.d))
            self.breaking_index = float(self.H_target / self.breaking_limit_miche) if self.breaking_limit_miche > 0 else 0.0
            self.is_breaking = bool(self.breaking_limit_miche > 0 and self.H_target > self.breaking_limit_miche)

            # Integral properties (from C++ invariants, then dimensionalised)
            self._calc_integral_props_cpp()

            # Kinematics summary (crest/trough surface and bed under crest)
            self.u_bed, _, _ = self.get_kinematics(0.0, 0.0)

            # Quadratic bed shear estimate kept as-is (engineering heuristic)
            cf_est = 0.005
            self.tau_bed = 0.5 * RHO * cf_est * (self.u_bed ** 2)

            self.ExcursionBed = abs(self.u_bed) * self.T_target / (2.0 * np.pi)

            # crest/trough surface velocities
            self.u_surf, _, _ = self.get_kinematics(self.d + self.eta_crest, 0.0)
            u_trough, _, _ = self.get_kinematics(self.d + self.eta_trough, np.pi)
            self.asymmetry = abs(self.u_surf / u_trough) if abs(u_trough) > 0 else 0.0

            # scan phases for max vertical velocity and horizontal acceleration on the surface
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

        # Sample one full cycle in phase X = kx (0..2π). For a steady wave,
        # spatial averaging over one wavelength equals temporal averaging at a point.
        phases = np.linspace(0.0, 2.0 * np.pi, int(max(36, nph)), endpoint=False)

        ub2 = 0.0
        for ph in phases:
            u_abs, _, _ = self.get_kinematics(z_bed=0.0, phase=float(ph))  # bed: z_bed=0
            u_orb = u_abs - float(self.Uc)  # remove imposed Eulerian current
            ub2 += u_orb * u_orb

        return float(ub2 / len(phases))


    # ------------------------ integral properties (C++ parity) ----------------

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
            # Pressure scaled by (ρ g d) (same expression as C++/Point())
            P = R - y - 0.5 * (urel * urel + v_nd * v_nd)
            integ[idx] = P + urel * urel

        return float(_np_trapz(integ, ys))

    def _calc_integral_props_cpp(self):
        """
        Compute integral quantities using the same invariants as the C++ Output().
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

        # C++ invariants in wave-number scaling
        pulse = float(self.z[8] + kd * self.z[5])
        ke = 0.5 * (self.z[4] * pulse - self.z[5] * Q_dimless * (kd ** 1.5))

        pe = 0.0
        for i in range(1, self.n + 1):
            pe += 0.25 * (self.Y[i] ** 2)

        ub2 = float(2.0 * self.z[9] - self.z[4] * self.z[4])
        q_term = float(self.z[7] * kd - self.z[8])

        sxx = float(4.0 * ke - 3.0 * pe + ub2 * kd + 2.0 * self.z[5] * q_term)
        f = float(self.z[4] * (3.0 * ke - 2.0 * pe) + 0.5 * ub2 * (pulse + self.z[4] * kd) + self.z[4] * self.z[5] * q_term)

        # Convert to depth-scaled dimensionless (as in C++ second column)
        E_depth = float((ke + pe) / (kd ** 2))
        KE_depth = float(ke / (kd ** 2))
        PE_depth = float(pe / (kd ** 2))

        # Store depth-scaled invariants (used for Solution-Flat style reporting)
        self.E_depth = E_depth
        self.KE_depth = KE_depth
        self.PE_depth = PE_depth
        Sxx_depth = float(sxx / (kd ** 2))
        F_depth = float(f / (kd ** 2.5))
        I_depth = float(pulse / (kd ** 1.5))

        self.Sxx_depth = Sxx_depth
        self.F_depth = F_depth
        self.I_depth = I_depth

        # Dimensionalise to GUI units
        self.EnergyDensity = float(RHO * self.g * (self.d ** 2) * E_depth)         # [J/m^2]
        self.Sxx = float(RHO * self.g * (self.d ** 2) * Sxx_depth)                 # [N/m]
        self.Power = float(RHO * (self.g ** 1.5) * (self.d ** 2.5) * F_depth)      # [W/m]
        # Momentum flux (Solution-Flat row 13) in moving frame
        self.MomentumFluxDepth = self._momentum_flux_S_depth(phase=0.0, npts=1200)  # S/(ρ g d²)
        self.MomentumFlux = float(RHO * self.g * (self.d ** 2) * self.MomentumFluxDepth)             # [N/m]
        self.Impulse = float(RHO * np.sqrt(self.g * (self.d ** 3)) * I_depth)  # [kg/(m·s)] per unit crest width

        self.BernoulliR = float(R_dimless * self.g * self.d)                       # [m^2/s^2] head* g? (consistent scalar)
        self.MassTransport = float(cs_dimless * np.sqrt(self.g * self.d))          # [m/s] (Stokes current)

        # Convenience values for Solution-Flat style reporting
        self.EulerianCurrent = float(self.Uc)                              # u1 [m/s]
        self.StokesCurrent = float(self.MassTransport)                      # u2 [m/s]
        self.MeanFluidSpeed = float(ubar_dimless * np.sqrt(self.g * self.d))  # Ū [m/s]
        self.VolumeFluxQ = float(Q_dimless * np.sqrt(self.g * (self.d ** 3))) # Q [m^2/s]
        self.WaveVolumeFlux_q = float(self.MeanFluidSpeed * self.d - self.VolumeFluxQ)  # q [m^2/s]
        self.BernoulliR_dimless = float(R_dimless)                          # R/(g d)
        self.Bernoulli_r = float((R_dimless - 1.0) * self.g * self.d)        # r = R - g d [m^2/s^2]
        self.KineticEnergy = float(RHO * self.g * (self.d ** 2) * self.KE_depth)     # [J/m^2]
        self.PotentialEnergy = float(RHO * self.g * (self.d ** 2) * self.PE_depth)  # [J/m^2]
        # Mean square of bed orbital velocity (Solution-Flat row 17):
        # u_b^2 = < (u_b(t) - ū₁)^2 >  (RMS-orbital-velocity definition; non-negative)
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

        rows.append({"_section": "DISPERSION & PHASE (GEOMETRY)"})
        rows += [
            ("Wavelength (L)", solver0.L, wc(solverC.L), "m"),
            ("Wave number (k)", solver0.k, wc(solverC.k), "rad/m"),
            ("kd", float(solver0.z[1]), wc(float(solverC.z[1])), "-"),
            ("Angular frequency (ω)", 2.0 * np.pi / solver0.T_target, wc(2.0 * np.pi / solverC.T_target), "rad/s"),
            ("Celerity / phase speed (c)", solver0.c, wc(solverC.c), "m/s"),
            ("c/√(gd)", solver0.c / sqrt_gd, wc(solverC.c / sqrt_gd), "-"),
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
            c = float(slv.c)
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