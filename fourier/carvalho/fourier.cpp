/*
================================================================================
fourier.cpp — Fenton-style reference solver for steady nonlinear water waves
================================================================================

PROGRAM ROLE IN THE SUITE
-------------------------
This file is a standalone Fourier / stream-function solver for steady,
two-dimensional, periodic finite-amplitude gravity waves over a horizontal bed,
with optional collinear current and either wavelength- or period-based closure.

The central engineering calculation is the nonlinear mapping between wave height,
period or wavelength, water depth, selected current convention and the resulting
consistent wave field.  In Fenton's formulation, a steady wave train is defined
by the coupled solution for wavelength, celerity, free surface, Fourier
coefficients, wave-frame flux, current variables, Bernoulli constant and integral
quantities.  Linear Airy dispersion is used only for initialization and checking;
it is not the final finite-amplitude theory.

PHYSICAL MODEL
--------------
The model is the classical Rienecker-Fenton / Fenton steady-wave problem:

  • homogeneous inviscid fluid,
  • incompressible and irrotational motion,
  • two-dimensional periodic wave train,
  • horizontal impermeable bed in finite depth,
  • atmospheric pressure on the free surface,
  • no breaking, viscosity, turbulence model or surface tension.

The numerical method is not a Stokes or cnoidal perturbation expansion.  The
Fourier order N is a collocation/spectral resolution parameter: increasing N
increases numerical resolution of the same nonlinear free-boundary problem.

COORDINATE AND DATUM CONVENTION
-------------------------------
Fenton's formulation uses the theoretical datum with still-water level at y=0
and the bed at y=-d.  This legacy source uses the equivalent bed-based physical coordinate
for output and profile sampling:

    bed:              y = 0
    mean water level: y = d

Internally it immediately shifts to a mean-level coordinate:

    y' = y - d,        X = k(x-ct),        Y = k y'

so that the finite-depth bed is Y=-kd and the mean water level is Y=0.  This is
only a vertical datum shift; it is not a different wave theory.

STREAM-FUNCTION FORMULATION
---------------------------
In the wave frame the free surface is steady.  The nondimensional stream function
Ψ satisfies Laplace's equation, the bed condition and the nonlinear free-surface
conditions:

    Ψ_XX + Ψ_YY = 0
    Ψ(X,-kd) = 0
    Ψ(X,ζ(X)) = constant
    1/2(Ψ_X² + Ψ_Y²) + ζ(X) = r*

The Fourier basis is chosen so that the field equation and bed/decay condition
are satisfied analytically.  Newton iteration then solves the nonlinear free-
surface boundary conditions and the global wave-current constraints.

FINITE-DEPTH FOURIER BASIS
--------------------------
The finite-depth representation uses the stable hyperbolic combination:

    S_j(Y) = sinh(jY) + cosh(jY) tanh(jkd)
    C_j(Y) = cosh(jY) + sinh(jY) tanh(jkd)

These are equivalent to sinh[j(Y+kd)]/cosh(jkd) and
cosh[j(Y+kd)]/cosh(jkd), but are more convenient in the shifted coordinate
Y=k(y-d).  In the deep-water limit tanh(jkd) -> 1 and both functions approach
exp(jY).

STATE VECTOR AND RESIDUAL SYSTEM
--------------------------------
The active legacy state vector is the Fenton-style 1-based vector:

    z[1]   = kd
    z[2]   = kH
    z[3]   = T sqrt(gk)
    z[4]   = c sqrt(k/g)
    z[5]   = Eulerian current u1 sqrt(k/g)
    z[6]   = Stokes / mass-transport current u2 sqrt(k/g)
    z[7]   = mean wave-frame speed Ubar sqrt(k/g)
    z[8]   = q sqrt(k^3/g), with q = Ubar*d - Q
    z[9]   = Bernoulli offset r k/g
    z[10..N+10]       = free-surface nodal ordinates ζ_m = kη_m
    z[N+11..2N+10]    = Fourier coefficients B_j

For Fourier order N the half-wave grid has N+1 collocation nodes
X_m=mπ/N, m=0..N.  The complete nonlinear system has

    8 + (N+1) + (N+1) = 2N + 10

residual equations.  The first eight impose global height, period/wavelength,
celerity, current, flux and datum constraints.  The next N+1 equations impose
the free-surface streamline condition.  The final N+1 equations impose
Bernoulli's equation on the free surface.

CURRENT CONVENTIONS
-------------------
The code distinguishes Fenton's two current definitions:

    u1 = c - Ubar              Eulerian current
    u2 = c - Q/d               Stokes / mass-transport current
    q  = Ubar*d - Q            wave-volume-flux variable

The input current criterion selects whether u1 or u2 is prescribed.  The period-
specified problem is therefore closed by both the observed period and the selected
current convention.

OUTPUT FILES
------------
The program writes the Fenton-style output files:

  • solution.res   — integral quantities, dimensionless groups and Fourier coefficients
  • surface.res    — free-surface ordinates, pressure and Bernoulli checks
  • flowfield.res  — velocity, acceleration and pressure profiles

The legacy redundant flat table is intentionally not written by this build.

INPUT FILES AND EMBEDDED CONTROLS
---------------------------------
The optional input file is data.dat.  If it is absent, the command-line interface
accepts finite-depth or deep-water cases with period or wavelength closure.  The historical Points.dat and
Convergence.dat controls are embedded as constants in this source.

COMPILATION
-----------
Windows / MSYS2 / MinGW:

    g++ -std=c++20 fourier.cpp -o fourier.exe -O3 -fno-fast-math -ffp-contract=off -fopenmp -static-libgcc -static-libstdc++

Linux / macOS:

    g++ -std=c++20 -Wall -Wextra -pedantic fourier.cpp -O3 -fno-fast-math -ffp-contract=off -fopenmp -o fourier

ATTRIBUTION
-----------
The hydrodynamic method follows J. D. Fenton's stream-function / Fourier
formulation and the Rienecker-Fenton numerical-wave tradition.  The dense SVD
linear solve follows the Numerical Recipes style used by the historical source
set.
================================================================================
*/

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cstddef>
#include <algorithm>

#include <filesystem>
#include <memory>
#include <charconv>
#include <string_view>
#include <string>
#include <system_error>
#include <vector>
#include <limits>
#include <map>
#include <stdexcept>
#include <cfenv>
#include <utility>

#if defined(_OPENMP)
  #include <omp.h>
#endif


#if defined(_WIN32)
  // Windows: robust executable path discovery (GetModuleFileNameW)
  #ifndef NOMINMAX
    #define NOMINMAX
  #endif
  #ifndef WIN32_LEAN_AND_MEAN
    #define WIN32_LEAN_AND_MEAN
  #endif
  #include <windows.h>
#endif

#if defined(__linux__)
  // Linux: robust executable path discovery (/proc/self/exe)
  #include <unistd.h>
  #include <limits.h>
#endif


namespace Phys {
    using Real = double;
    constexpr Real PI    = 3.141592653589793238462643383279502884;
}
using Phys::Real;

// Enclosures use IEEE binary64 arithmetic in round-to-nearest mode.
// Compile without fast-math or floating-point reassociation.
#if defined(__FAST_MATH__)
#error "Finite-system interval verification requires compilation without fast-math."
#endif
namespace Verified {
struct Interval {
    Real lo, hi;
    Interval(Real value = 0.0) : Interval(value, value) {}
    Interval(Real lower, Real upper) : lo(lower), hi(upper) {
        if (!std::isfinite(lo) || !std::isfinite(hi) || lo > hi)
            throw std::runtime_error("Invalid interval endpoints.");
    }
    static Real down(Real x) { return std::nextafter(x, -std::numeric_limits<Real>::infinity()); }
    static Real up(Real x) { return std::nextafter(x, std::numeric_limits<Real>::infinity()); }
    friend Interval operator+(const Interval& a, const Interval& b) {
        return {down(a.lo+b.lo), up(a.hi+b.hi)};
    }
    friend Interval operator-(const Interval& a) { return {-a.hi, -a.lo}; }
    friend Interval operator-(const Interval& a, const Interval& b) { return a+(-b); }
    friend Interval operator*(const Interval& a, const Interval& b) {
        const Real products[] = {a.lo*b.lo,a.lo*b.hi,a.hi*b.lo,a.hi*b.hi};
        for (Real p : products) if (!std::isfinite(p)) throw std::runtime_error("Interval product overflow.");
        return {down(*std::min_element(products,products+4)), up(*std::max_element(products,products+4))};
    }
    friend Interval operator/(const Interval& a, const Interval& b) {
        if (b.lo <= 0.0 && b.hi >= 0.0) throw std::runtime_error("Interval divisor contains zero.");
        return a*Interval(down(1.0/b.hi),up(1.0/b.lo));
    }
    Interval square() const {
        const Real lower = lo <= 0.0 && hi >= 0.0 ? 0.0 : std::min(std::abs(lo),std::abs(hi));
        const Real upper = magnitude();
        return {std::max(0.0,down(lower*lower)),up(upper*upper)};
    }
    Interval power(unsigned n) const {
        Interval value(1.0), base(*this);
        while (n) {
            if (n & 1U) value=value*base;
            n >>= 1U;
            if (n) base=base.square();
        }
        return value;
    }
    Real magnitude() const { return std::max(std::abs(lo),std::abs(hi)); }
    static Interval exponential_point(Real x) {
        if (x==0.0) return Interval(1.0);
        Real reduced=std::abs(x);
        unsigned halvings=0;
        while (reduced>0.125) { reduced*=0.5; ++halvings; }
        const Interval r(std::max(0.0,down(reduced)),up(reduced));
        Interval sum(1.0), term(1.0);
        for (unsigned j=1;j<=32;++j) { term=term*r/Interval(j); sum=sum+term; }
        // For 0 <= r <= 1/4, exp(r) < 2 bounds the Taylor remainder.
        Interval remainder=r.power(33)*Interval(2.0);
        for (unsigned j=2;j<=33;++j) remainder=remainder/Interval(j);
        sum=Interval(sum.lo,up(sum.hi+remainder.hi));
        if (x<0.0) sum=Interval(1.0)/sum;
        for (unsigned j=0;j<halvings;++j) sum=sum.square();
        return sum;
    }
    Interval exp() const {
        return {exponential_point(lo).lo,exponential_point(hi).hi};
    }
    Interval sqrt() const {
        if (lo<0.0) throw std::runtime_error("Negative square-root interval.");
        if (hi==0.0) return Interval(0.0);
        Real lower=lo==0.0 ? 0.0 : std::max(0.0,down(std::sqrt(lo)));
        Real upper=up(std::sqrt(hi));
        // Squared endpoint tests enclose the root independently of the initial estimate.
        for (int i=0;i<16 && lower>0.0 && up(lower*lower)>lo;++i) lower=std::max(0.0,down(lower));
        for (int i=0;i<16 && down(upper*upper)<hi;++i) upper=up(upper);
        if ((lower>0.0 && up(lower*lower)>lo) || down(upper*upper)<hi)
            throw std::runtime_error("Square-root enclosure failed.");
        return {lower,upper};
    }
};

inline Interval pi() { return {Interval::down(Phys::PI),Interval::up(Phys::PI)}; }

inline std::vector<std::pair<Interval,Interval>> collocation(int n) {
    std::vector<std::pair<Interval,Interval>> table;
    table.reserve((size_t)2*n);
    for (int index=0;index<2*n;++index) {
        const int q=index<=n ? index : index-2*n;
        const Interval x=pi()*Interval(q)/Interval(n), square=x.square();
        Interval sine(0.0),cosine(0.0),st=x,ct(1.0);
        for (int m=0;m<32;++m) {
            sine=sine+st;cosine=cosine+ct;
            st=-st*square/Interval((2*m+2)*(2*m+3));
            ct=-ct*square/Interval((2*m+1)*(2*m+2));
        }
        Interval sin_error=Interval(x.magnitude()).power(65);
        Interval cos_error=Interval(x.magnitude()).power(64);
        for (int j=2;j<=65;++j) sin_error=sin_error/Interval(j);
        for (int j=2;j<=64;++j) cos_error=cos_error/Interval(j);
        table.emplace_back(sine+Interval(-sin_error.hi,sin_error.hi),
                           cosine+Interval(-cos_error.hi,cos_error.hi));
    }
    return table;
}

struct OutputBound {
    Interval enclosure;
    Real reported=0.0, error_upper=0.0, budget_lower=0.0;
};
struct Certificate {
    std::string arithmetic="Outward-rounded IEEE binary64 intervals";
    std::string scope="Finite-state residual and nine selected finite-formula outputs";
    bool root_existence_proved=false, passed=false;
    int order=0, equation_count=0;
    std::vector<Interval> residual_intervals;
    std::vector<Real> state, surface_coefficients;
    Real residual_upper=0.0, residual_tolerance=1e-10;
    std::map<std::string,OutputBound> outputs;
};
} // namespace Verified

namespace fs = std::filesystem;


/*
================================================================================
TECHNICAL ROADMAP FOR THIS TRANSLATION UNIT
================================================================================

This file is intentionally "monolithic": it embeds what historically were several
source modules (Allocation, Solve, Subroutines, Inout, main) into a single build
unit to simplify compilation and reproducibility.

The numerical core uses a Fourier collocation system with adaptive continuation
and numerical acceptance checks before writing result files.

-------------------------------------------------------------------------------
1) COORDINATE SYSTEMS AND FRAMES OF REFERENCE (VERY IMPORTANT)
-------------------------------------------------------------------------------

Two coordinate systems are used throughout, and confusing them is the most common
source of mistakes when extending or validating the solver:

(1) Physical coordinates (x, y)
    - x : horizontal coordinate, positive in the direction of wave propagation
    - y : vertical coordinate, positive upward
    - For finite depth:
        bed is at y = 0
        mean free surface is at y = d  (mean depth)
    - The wave is periodic with wavelength λ and period τ (T in output).

(2) Non-dimensional "solver" coordinates (X, Y)
    - Define the travelling horizontal coordinate ξ := x - c t.
    - X := k ξ   where k = 2π/λ is the wavenumber.
      This makes one wavelength correspond to X ∈ [0, 2π].
    - For *finite depth*, the vertical coordinate used internally by Point() is:
        Y := k (y - d)   so that:
          mean free surface: y = d  -> Y = 0
          bed:              y = 0  -> Y = -kd
      The variable kd appears explicitly as an unknown in the period-
      specified finite-depth problem.

    - For *deep water* the basis functions use exp(j Y), representing the
      kd → ∞ limit of the finite-depth cosh/sinh basis. In deep water there is
      no physical bed; any "vertical extent" used in output profiles is purely
      a plotting choice and does not represent a boundary condition.

Frames:
    - The governing equations and Fourier expansion are formulated in a frame
      moving with the wave at c (steady wave in that frame).
    - Output kinematics (u, v, accelerations) are finally reported in the
      *stationary* frame using the steady-wave identity ∂/∂t = -c ∂/∂x.

-------------------------------------------------------------------------------
2) NON-DIMENSIONALISATION (WHY TWO SETS OF NUMBERS APPEAR)
-------------------------------------------------------------------------------

Fenton's stream-function formulation is often written in a k-based scaling
(think "use g and k to scale everything") because Fourier series are naturally
expressed in terms of X = k ξ (the non-dimensional travelling-phase coordinate).

For finite depth, engineers also want depth-based scaling (g and d), because
design inputs are typically H/d and T√(g/d). Therefore, this code reports many
quantities twice:

  (A) k-based scaling (native to the solver):
      lengths  ~ 1/k
      speeds   ~ √(g/k)
      times    ~ √(1/(gk))

  (B) d-based scaling (more common in coastal engineering):
      lengths  ~ d
      speeds   ~ √(g d)
      times    ~ √(d/g)

Conversion between them introduces factors of kd = k d.
Example:
    u(gd) = u(gk) / √(kd)          (because √(g/k) = √(g d) / √(k d))

You will see these √(kd) and kd^(3/2) factors repeatedly. They are not arbitrary:
they enforce dimensional consistency between the two nondimensional systems.

-------------------------------------------------------------------------------
3) UNKNOWN VECTOR z[] AND RESIDUAL VECTOR rhs[] (LEGACY INDEXING)
-------------------------------------------------------------------------------

This solver uses Numerical-Recipes style 1-based arrays. The active state vector
is the Fenton-style vector z[1..2N+10].  Index 0 is unused.

A convenient mental model is:

  z = [ global scalars | surface ordinates | Fourier coefficients ]

The physical correspondence is:

  z[1]   : kd                         dimensionless water depth
  z[2]   : kH                         dimensionless wave height
  z[3]   : T√(gk)                     k-based dimensionless apparent period
  z[4]   : c√(k/g)                    k-based wave celerity
  z[5]   : u1√(k/g)                   Eulerian current
  z[6]   : u2√(k/g)                   Stokes / mass-transport current
  z[7]   : Ubar√(k/g)                 mean wave-frame speed
  z[8]   : q√(k^3/g), q=Ubar*d-Q      wave-volume-flux variable
  z[9]   : r k/g                      Bernoulli offset

  z[10 + m],  m=0..N:
      ζ_m = kη_m, free-surface nodal ordinate at X_m=mπ/N.

  z[n + 10 + j], j=1..N:
      B_j, stream-function Fourier coefficient.

The residual vector rhs[] has the required size and ordering:

  rhs[1..8]
      global height, period/wavelength, celerity, current, flux and datum
      constraints;

  rhs[9..9+N]
      free-surface streamline / kinematic collocation residuals;

  rhs[n+10..n+10+N]
      free-surface Bernoulli / dynamic collocation residuals.

For Fourier order N the system size is:

  8 + (N+1) + (N+1) = 2N + 10.

The Newton solver drives every component of rhs[] to zero.

-------------------------------------------------------------------------------
4) COLLOCATION GRIDS AND TRIGONOMETRIC TABLES (cosa[], sina[])
-------------------------------------------------------------------------------

For efficiency (and to mirror the original code), the program precomputes:

  cosa[p] = cos(p π / N)    for p=0..2N
  sina[p] = sin(p π / N)    for p=0..2N

At collocation point m and Fourier mode j, the phase is:
  X_m = m π / N
  j X_m = (m j) π / N

The index nm = (m*j) % (2N) maps directly into these precomputed tables, exploiting
periodicity of cos/sin over 2π.

-------------------------------------------------------------------------------
5) LINEAR ALGEBRA AND ROBUSTNESS CHOICES
-------------------------------------------------------------------------------

The Newton step solves a dense system with an analytic Jacobian.
Near limiting steep waves, the Jacobian can be ill-conditioned; therefore the
solver uses an SVD-based least-squares solve (Numerical Recipes dsvdcmp/dsvbksb).

The singular-value cutoff wmin = wmax * 1e-15 is a pragmatic stabilisation:
small singular values are treated as zero to suppress amplification of numerical
noise in nearly-singular directions.

-------------------------------------------------------------------------------
6) EXTENDING THIS CODE SAFELY
-------------------------------------------------------------------------------

If you modify/extend this file later:

  - Treat z[] and rhs[] indices as an API. Any reindexing will silently break the
    solver unless *all* dependent formulas are updated consistently.

  - Keep the "Is_finite / Is_deep" branches semantically aligned: in finite depth
    the vertical basis is sinh/cosh with tanh(j kd) terms; in deep water the basis
    is exp(j Y).

  - Do not replace the NR-style allocation with std::vector without accounting
    for 1-based indexing and the pointer offsets in dvector/dmatrix.

================================================================================
*/

/*
================================================================================
FENTON STEADY NONLINEAR WAVE THEORY (STREAM-FUNCTION / FOURIER-COLLOCATION)
================================================================================

This section maps the Fenton formulation to the variables and equations in this
single-file solver.  The code body uses historical variable names and 1-based
arrays while solving the same nonlinear free-boundary problem.

------------------------------------------------------------------------------
A) PHYSICAL MODEL AND TRAVELLING-WAVE FORMULATION
------------------------------------------------------------------------------

The physical model is a steady progressive gravity wave of permanent form in a
homogeneous inviscid fluid.  The flow is incompressible and irrotational, so a
velocity potential and a stream function exist and both satisfy Laplace's
equation in the fluid domain.

The fixed-frame wave has period T, wavelength L, wavenumber k=2π/L and celerity
c=L/T.  The travelling coordinate is ξ=x-ct.  In the wave frame the free surface
is steady and all unknowns are periodic in X=kξ.

Fenton's formulation uses y=0 at still-water level and y=-d at the bed. This
source uses a bed-based physical coordinate y for output, then shifts to y'=y-d.
Therefore:

    Fenton datum:      bed=-d, mean level=0
    this source output: bed=0,  mean level=d
    internal solver:    Y=k(y-d), so bed=-kd and mean level=0

Only the datum changes.  The nondimensional ordinate stored in z[10+m] is the
same ζ_m=kη_m used in the governing equations.

------------------------------------------------------------------------------
B) STREAM FUNCTION, VELOCITIES AND BOUNDARY CONDITIONS
------------------------------------------------------------------------------

In k-based nondimensional wave-frame variables the stream function Ψ satisfies:

    Ψ_XX + Ψ_YY = 0

and the wave-frame velocities are:

    Uhat =  Ψ_Y
    Vhat = -Ψ_X

The finite-depth bed is a streamline:

    Ψ(X,-kd) = 0

The free surface is also a streamline, and Bernoulli's equation is constant along
it:

    Ψ(X,ζ(X)) = constant
    1/2(Ψ_X² + Ψ_Y²) + ζ(X) = r*

The constant-streamline form used by the code is algebraically rearranged through
z[7] and z[8].  At each collocation node it is assembled as:

    ψ_Fourier(X_m,ζ_m) - z[8] - z[7] ζ_m = 0

where z[8]=q√(k^3/g) and q=Ubar*d-Q.  This is the free-surface residual
R_m^(ψ)=0 written in the legacy variable layout.

The Bernoulli residual is assembled as:

    1/2[(-z[7]+u_m)^2 + v_m^2] + ζ_m - z[9] = 0

where u_m and v_m are the Fourier perturbation velocity sums at the free surface.

------------------------------------------------------------------------------
C) FOURIER BASIS
------------------------------------------------------------------------------

The finite-depth Fourier representation satisfies Laplace's equation and the bed
condition before the nonlinear equations are assembled.  In the shifted coordinate
Y=k(y-d), the stable basis functions are:

    S_j(Y) = sinh(jY) + cosh(jY)tanh(jkd)
    C_j(Y) = cosh(jY) + sinh(jY)tanh(jkd)

The code evaluates ψ, Ψ_Y and -Ψ_X from these S_j and C_j combinations.  In deep
water tanh(jkd)->1, so S_j and C_j reduce to exp(jY), matching the formulation's
kd->∞ limit.

------------------------------------------------------------------------------
D) CURRENT DEFINITIONS
------------------------------------------------------------------------------

The current convention is not a post-processing correction.  It closes a different
nonlinear problem.  The two Fenton current definitions are:

    u1 = c - Ubar        Eulerian current
    u2 = c - Q/d         Stokes / mass-transport current
    q  = Ubar*d - Q      wave-volume-flux variable

The corresponding residuals in this source are:

    rhs[4] = z[5] + z[7] - z[4]
    rhs[5] = z[1]*(z[6] + z[7] - z[4]) - z[8]
    rhs[6] = selected-current variable - specified current

The selected current criterion determines whether z[5] or z[6] is prescribed by
the input current value.

------------------------------------------------------------------------------
E) COLLOCATION AND NONLINEAR SYSTEM SIZE
------------------------------------------------------------------------------

The half-wave collocation grid is:

    X_m = mπ/N,    m=0,1,...,N

Symmetry about crest and trough makes this half-wave grid sufficient.  For order
N the unknowns are:

    9 global scalars + (N+1) surface ordinates + N Fourier coefficients
    = 2N + 10 active entries.

The residuals are:

    8 global equations + (N+1) streamline equations + (N+1) Bernoulli equations
    = 2N + 10 equations.

Thus the implemented Fenton system is square and uses the same active state-vector
meaning used by the governing equations.

------------------------------------------------------------------------------
F) POST-PROCESSING
------------------------------------------------------------------------------

After Newton convergence the code transforms the nodal ζ_m values into cosine
surface coefficients Y[j] for smooth evaluation.  Point() then evaluates local
stationary-frame velocity, pressure and acceleration by combining the solved
Fourier coefficients with the steady-wave identity ∂/∂t=-c∂/∂x.

Integral quantities are reported in both k-based and depth-based scalings where
finite depth is available, matching the dual nondimensional interpretation in the
the formulation and in Fenton's SOLUTION.RES-style output.

================================================================================
*/

/*
------------------------------------------------------------------------------
PORTABILITY + MODERN C++ RESOURCE MANAGEMENT (C++20 / Core Guidelines)
------------------------------------------------------------------------------
The original upstream codebase relied heavily on:
  - global raw FILE* handles, opened/closed manually,
  - Numerical Recipes (NR) 1-based dynamic allocation via malloc/free,
  - optional platform-specific “press any key” pauses (getch()).

This single-file version removes any platform-specific
console dependencies and introduces RAII ownership for:
  - FILE* handles (unique_ptr with fclose deleter),
  - the long-lived NR-allocated vectors/matrices used by the solver.

The numerical solve uses an analytic Jacobian, column equilibration, damped
Newton steps and adaptive height continuation. Fourier order is reduced when
the equilibrated Jacobian condition exceeds 1e11. Every accepted solution must
satisfy a maximum absolute residual of 1e-10 and a sampled surface residual of
1e-6. Finite-depth integral outputs are independently checked by quadrature in
k-based units and outward-rounded intervals in depth-based units. In each
scaling, the absolute-plus-relative output budget is 1e-7 + 1e-8*abs(reference).
These checks do not constitute a Lean proof of program execution, root existence
or error bounds for the continuous free-boundary problem.
------------------------------------------------------------------------------
*/

struct FileCloser {
  void operator()(std::FILE* f) const noexcept {
    if (f) std::fclose(f);
  }
};

using FilePtr = std::unique_ptr<std::FILE, FileCloser>;

// -----------------------------------------------------------------------------
// "Headers.h" content (flattened)
// -----------------------------------------------------------------------------

// Text helpers for depth/case dispatch and numerical output.
static constexpr double pi = 3.14159265358979324;
#define iff(x, y)        if (std::strcmp((x), #y) == 0)
#define HI               " % 23.15e"
#define LO               "\t%8.4f"
#define Is_deep          if (std::strcmp(Depth, "Deep") == 0)
#define Is_finite        if (std::strcmp(Depth, "Finite") == 0)

/*
------------------------------------------------------------------------------
GLOBAL PROGRAM STATE (WHY SO MANY `static` VARIABLES EXIST)
------------------------------------------------------------------------------
This solver is a direct descendant of a procedural FORTRAN/C style program.
For strict fidelity to Fenton-style listings, most state is kept in file-scope
variables rather than being passed explicitly as function arguments.

Advantages (in this legacy context):
  - avoids large argument lists (z, rhs, trig tables, coefficients) everywhere,
  - preserves the original indexing and variable naming used in publications,
  - simplifies the numerical kernels (Eqns, Point, Output) which are performance-
    sensitive and were historically written with globals.

Costs:
  - global state is confined to one problem case; parallel residual evaluations use independent state copies,
  - "hidden" dependencies between routines (documented in comments where possible).

Conventions:
  - Arrays allocated with dvector/dmatrix are 1-based (Numerical Recipes style),
    except where explicitly allocated from 0 for convenience (e.g., Y[0..num]).
  - All variables are `static` to keep internal linkage inside this translation
    unit, reducing the chance of symbol collisions if you embed this into a larger
    project.
------------------------------------------------------------------------------
*/

// Global I/O streams
static std::FILE *monitor = nullptr;
static std::FILE *Input1 = nullptr;
static std::FILE *Elevation = nullptr;
static std::FILE *Flowfield = nullptr;
static std::FILE *Solution = nullptr;

// Global text/buffers (sizes chosen to match fgets(...,400,...) usage)
static char Title[400] = {0};
static char dummy[400] = {0};
static char Case[20] = {0};
static char Currentname[10] = {0};
static char Current1[10] = "Euler";
static char Current2[10] = "Stokes";
static char Depth[100] = {0};
static char Method[100] = {0};

// -----------------------------------------------------------------------------
// Embedded convergence parameters (formerly read from Convergence.dat)
// -----------------------------------------------------------------------------
// These are the exact values shipped in the provided Convergence.dat:
//   50     maximum iterations per height step
//   1.e-9  convergence tolerance
//
// Convergence.dat is intentionally NOT read at runtime.
static constexpr int    kConvergenceMaxIterations = 50;
static constexpr double kConvergenceCriterion     = 1.0e-9;


// -----------------------------------------------------------------------------
// Embedded output sampling controls (formerly read from Points.dat)
// -----------------------------------------------------------------------------
// Points.dat contained:
//
//   50  Number of points on free surface (clustered near crest)
//   18  Number of velocity profiles over half a wavelength
//   10  Number of vertical points in each profile
//
// Points.dat is no longer required at runtime.
static constexpr int kSurfacePointsDefault = 50;
static constexpr int kNprofilesDefault     = 18;
static constexpr int kProfilePointsDefault = 10;


// Diagnostic (present in the provided sources)
static char Diagname[30] = {0};
static char Theory[10] = {0};

// Global ints
static int Current_criterion = 0;
static int n = 0;
static int Nprofiles = 0;
static int nstep = 0;
static int num = 0;
static int number = 0;
static int Points = 0;
static int Surface_points = 0;

// Global arrays
static double **sol = nullptr;
static double *B = nullptr;
static double *cosa = nullptr;
static double *rhs1 = nullptr;
static double *sina = nullptr;
static double *Tanh = nullptr;
static double *Y = nullptr;
static double *z = nullptr;

// Global doubles
/*
------------------------------------------------------------------------------
SCALAR VARIABLES: PHYSICAL MEANING + UNIT SYSTEM
------------------------------------------------------------------------------
Most scalar doubles below correspond to quantities appearing in standard steady
wave theory. Their *numerical values* depend on the current depth regime:

Finite depth (Depth="Finite"):
  - Solver works primarily in k-based non-dimensional variables:
      * lengths  scaled by 1/k
      * speeds   scaled by √(g/k)
      * energies scaled by g/k², etc.
  - Many reported values are also converted to d-based scaling.

Deep water (Depth="Deep"):
  - kd is not a physical unknown (there is no finite bed); the code uses a
    normalisation z[1]=1 and exp(jY) basis functions.

Key scalars used across routines:

  * c, T, L, H
      Physical wave celerity, period, wavelength, height (in depth-based scaling
      for finite depth). They are derived from z[] after convergence.

  * ce, cs, ubar
      Mean currents and mean flow in wave frame:
        ce = Eulerian mean current ū1
        cs = Stokes/mass-transport mean current ū2
        ubar = mean flow in moving frame U

  * Q, q
      Flux-related quantities:
        Q  = depth-integrated discharge in the moving frame
        q  = Ubar*d - Q, the wave-volume-flux variable used in z[8]

  * R, r
      Bernoulli constants (depth-based and k-based variants)

  * u, v, ux, uy, vx, vy, ut, vt, dphidt, dudt, dvdt, Pressure, Bernoulli_check
      Local kinematics and diagnostics produced by Point():
        u, v          velocities (stationary frame)
        ux, uy, ...   spatial derivatives
        ut, vt        time derivatives via steady-wave identity
        dudt, dvdt    material accelerations
        dphidt        ∂φ/∂t (stationary frame)
        Pressure      dynamic pressure from Bernoulli rearrangement
        Bernoulli_check residual of Bernoulli constant (should be ~0)

  * ke, pe, sxx, f, pulse
      Integral quantities:
        ke    kinetic energy-like integral
        pe    potential energy-like integral
        sxx   radiation stress component
        f     wave power
        pulse impulse-like quantity

  * ub2
      Mean square *orbital* bed velocity (as implemented in MeanSquareBedOrbitalVelocity)

Many of these are written exactly as in the original codebase, including the
somewhat terse variable names. Comments in Output() provide the exact formulas
used for each reported integral quantity.
------------------------------------------------------------------------------
*/

static double Bernoulli_check = 0.0;
static double c = 0.0;
static double ce = 0.0;
static double crit = 0.0;
static constexpr double kFinalResidualTolerance = 1.0e-10;
static constexpr double kConditionLimit = 1.0e11;
static double jacobian_condition = std::numeric_limits<double>::infinity();
static double final_residual = std::numeric_limits<double>::infinity();
static double spectral_defect = std::numeric_limits<double>::infinity();
static int requested_order = 0;
static Verified::Certificate finite_system_certificate;
static double integral_error_ratio = 0.0;
static double spectral_error_ratio = 0.0;
static double cs = 0.0;
static double Current = 0.0;
static double dphidt = 0.0;
static double dudt = 0.0;
static double dvdt = 0.0;
static double f = 0.0;
static double H = 0.0;
static double Height = 0.0;
static double Highest = 0.0;
static double Hoverd = 0.0;
static double height = 0.0;
static double kd = 0.0;
static double ke = 0.0;
static double L = 0.0;
static double MaxH = 0.0;
static double pe = 0.0;
static double phi = 0.0;
static double Pressure = 0.0;
static double psi = 0.0;
static double pulse = 0.0;
static double Q = 0.0;
static double q = 0.0;
static double R = 0.0;
static double r = 0.0;
static double s = 0.0;
static double sxx = 0.0;
static double T = 0.0;
static double u = 0.0;
static double ub2 = 0.0;
static double ubar = 0.0;
static double ut = 0.0;
static double ux = 0.0;
static double uy = 0.0;
static double v = 0.0;
static double vt = 0.0;
static double vx = 0.0;
static double vy = 0.0;

// -----------------------------------------------------------------------------
// Numerical Recipes-style allocation (Util.cpp + Allocation.h)
// -----------------------------------------------------------------------------

/*
------------------------------------------------------------------------------
NUMERICAL RECIPES (NR) 1-BASED ALLOCATION HELPERS
------------------------------------------------------------------------------
The classic Numerical Recipes codebase allocates arrays with *user-selected* lower
and upper bounds [nl..nh] and then returns a pointer that is intentionally shifted
so that indexing starts at nl rather than 0.

Example:
    double* v = dvector(1, n);
    v[1] ... v[n] are valid.

Implementation detail:
    The allocated block contains (nh-nl+1+NR_END) elements. The returned pointer
    is offset by (-nl + NR_END). This means:
      - `v` is *not* the original malloc pointer,
      - you must free with free_dvector(v, nl, nh), which re-applies the inverse
        offset before calling free().

Why keep this pattern?
  - It matches the original Fenton-style program listings.
  - It avoids rewriting every loop boundary and index mapping in the solver.
  - It preserves exact numerical behaviour (ordering of operations, etc.).

Safety notes:
  - These routines perform no null checks. If allocation fails, behaviour is
    undefined. In practice, problem sizes are small (N ~ O(10^2)), so this is
    acceptable for the intended use.
  - Because pointer arithmetic is involved, *never* pass these pointers to
    standard containers assuming they start at element 0.

RAII wrappers used in main() ensure these allocations are still released safely.
------------------------------------------------------------------------------
*/

#define NR_END 1
#define FREE_ARG char*
static double *dvector(long nl, long nh)
{
  double *v = (double*)std::malloc((std::size_t)((nh - nl + 1 + NR_END) * sizeof(double)));
  return v - nl + NR_END;
}

static double **dmatrix(long nrl, long nrh, long ncl, long nch)
{
  long i;
  const long nrow = nrh - nrl + 1;
  const long ncol = nch - ncl + 1;
  double **m = (double**)std::malloc((std::size_t)((nrow + NR_END) * sizeof(double*)));
  m += NR_END;
  m -= nrl;

  m[nrl] = (double*)std::malloc((std::size_t)((nrow * ncol + NR_END) * sizeof(double)));
  m[nrl] += NR_END;
  m[nrl] -= ncl;

  for (i = nrl + 1; i <= nrh; i++) m[i] = m[i - 1] + ncol;
  return m;
}

static void free_dvector(double *v_, long nl, long /*nh*/)
{
  std::free((FREE_ARG)(v_ + nl - NR_END));
}

static void free_dmatrix(double **m_, long nrl, long /*nrh*/, long ncl, long /*nch*/)
{
  std::free((FREE_ARG)(m_[nrl] + ncl - NR_END));
  std::free((FREE_ARG)(m_ + nrl - NR_END));
}

// -----------------------------------------------------------------------------
// dpythag / SVD (Dpythag.cpp, Dsvdcmp.cpp, Dsvbksb.cpp)
// -----------------------------------------------------------------------------

/*
------------------------------------------------------------------------------
SINGULAR VALUE DECOMPOSITION (SVD) KERNEL (NUMERICAL RECIPES STYLE)
------------------------------------------------------------------------------
The routines dsvdcmp() and dsvbksb() are a near-direct transcription of the
Numerical Recipes in C SVD implementation.

High-level purpose in this program:
  - Newton's method requires solving J Δz = -rhs.
  - J is a dense analytic Jacobian.
  - For steep waves and/or large N, J can become ill-conditioned.
  - SVD decomposes J into U W V^T, allowing stable least-squares solutions and
    controlled handling of near-singular directions (via small singular values).

Algorithm sketch:
  1) Householder reduction of A to bidiagonal form.
  2) Accumulation of right-hand and left-hand transformations (V and U).
  3) QR iteration on the bidiagonal matrix to obtain singular values (w) and
     final orthogonal factors (U, V).

Important: indexing
  - These routines expect 1-based arrays: a[1..m][1..n], w[1..n], v[1..n][1..n].
  - The internal rv1 work vector is also 1-based.

Modification policy:
  - Treat this code as "library" code. It is not wave-physics-specific.
  - Changes here can subtly alter convergence properties; do not modify unless you
    have strong numerical linear algebra reasons and regression tests.

Known quirk (legacy):
  - The original NR code uses 30 QR iterations as a default; the warning message
    says "30 dsvdcmp iterations" even though the loop bound is 30. We preserve the
    behaviour and the message text.
------------------------------------------------------------------------------
*/

static double dsqrarg;
#define DSQR(a) ((dsqrarg = (a)) == 0.0 ? 0.0 : dsqrarg * dsqrarg)

static double dpythag(double a, double b)
{
  const double absa = std::fabs(a);
  const double absb = std::fabs(b);
  if (absa > absb) return absa * std::sqrt(1.0 + DSQR(absb / absa));
  return (absb == 0.0 ? 0.0 : absb * std::sqrt(1.0 + DSQR(absa / absb)));
}

static double dmaxarg1, dmaxarg2;
#define DMAX(a, b) (dmaxarg1 = (a), dmaxarg2 = (b), (dmaxarg1) > (dmaxarg2) ? (dmaxarg1) : (dmaxarg2))
static int iminarg1, iminarg2;
#define IMIN(a, b) (iminarg1 = (a), iminarg2 = (b), (iminarg1) < (iminarg2) ? (iminarg1) : (iminarg2))
#define SIGN(a, b) ((b) >= 0.0 ? std::fabs(a) : -std::fabs(a))

static bool dsvdcmp(double **a, int m, int n, double w[], double **v)
{
  int flag = 0, i = 0, its = 0, j = 0, jj = 0, k = 0, l = 0, nm = 0;
  double anorm, c_, f_, g_, h_, s_, scale, x_, y_, z_;
  double *rv1;

  rv1 = dvector(1, n);
  g_ = scale = anorm = 0.0;
  for (i = 1; i <= n; i++) {
    l = i + 1;
    rv1[i] = scale * g_;
    g_ = s_ = scale = 0.0;
    if (i <= m) {
      for (k = i; k <= m; k++) scale += std::fabs(a[k][i]);
      if (scale) {
        for (k = i; k <= m; k++) {
          a[k][i] /= scale;
          s_ += a[k][i] * a[k][i];
        }
        f_ = a[i][i];
        g_ = -SIGN(std::sqrt(s_), f_);
        h_ = f_ * g_ - s_;
        a[i][i] = f_ - g_;
        for (j = l; j <= n; j++) {
          for (s_ = 0.0, k = i; k <= m; k++) s_ += a[k][i] * a[k][j];
          f_ = s_ / h_;
          for (k = i; k <= m; k++) a[k][j] += f_ * a[k][i];
        }
        for (k = i; k <= m; k++) a[k][i] *= scale;
      }
    }
    w[i] = scale * g_;
    g_ = s_ = scale = 0.0;
    if (i <= m && i != n) {
      for (k = l; k <= n; k++) scale += std::fabs(a[i][k]);
      if (scale) {
        for (k = l; k <= n; k++) {
          a[i][k] /= scale;
          s_ += a[i][k] * a[i][k];
        }
        f_ = a[i][l];
        g_ = -SIGN(std::sqrt(s_), f_);
        h_ = f_ * g_ - s_;
        a[i][l] = f_ - g_;
        for (k = l; k <= n; k++) rv1[k] = a[i][k] / h_;
        for (j = l; j <= m; j++) {
          for (s_ = 0.0, k = l; k <= n; k++) s_ += a[j][k] * a[i][k];
          for (k = l; k <= n; k++) a[j][k] += s_ * rv1[k];
        }
        for (k = l; k <= n; k++) a[i][k] *= scale;
      }
    }
    anorm = DMAX(anorm, (std::fabs(w[i]) + std::fabs(rv1[i])));
  }
  for (i = n; i >= 1; i--) {
    if (i < n) {
      if (g_) {
        for (j = l; j <= n; j++) v[j][i] = (a[i][j] / a[i][l]) / g_;
        for (j = l; j <= n; j++) {
          for (s_ = 0.0, k = l; k <= n; k++) s_ += a[i][k] * v[k][j];
          for (k = l; k <= n; k++) v[k][j] += s_ * v[k][i];
        }
      }
      for (j = l; j <= n; j++) v[i][j] = v[j][i] = 0.0;
    }
    v[i][i] = 1.0;
    g_ = rv1[i];
    l = i;
  }
  for (i = IMIN(m, n); i >= 1; i--) {
    l = i + 1;
    g_ = w[i];
    for (j = l; j <= n; j++) a[i][j] = 0.0;
    if (g_) {
      g_ = 1.0 / g_;
      for (j = l; j <= n; j++) {
        for (s_ = 0.0, k = l; k <= m; k++) s_ += a[k][i] * a[k][j];
        f_ = (s_ / a[i][i]) * g_;
        for (k = i; k <= m; k++) a[k][j] += f_ * a[k][i];
      }
      for (j = i; j <= m; j++) a[j][i] *= g_;
    } else {
      for (j = i; j <= m; j++) a[j][i] = 0.0;
    }
    ++a[i][i];
  }
  for (k = n; k >= 1; k--) {
    for (its = 1; its <= 100; its++) {
      flag = 1;
      for (l = k; l >= 1; l--) {
        nm = l - 1;
        if ((double)(std::fabs(rv1[l]) + anorm) == anorm) {
          flag = 0;
          break;
        }
        if ((double)(std::fabs(w[nm]) + anorm) == anorm) break;
      }
      if (flag) {
        c_ = 0.0;
        s_ = 1.0;
        for (i = l; i <= k; i++) {
          f_ = s_ * rv1[i];
          rv1[i] = c_ * rv1[i];
          if ((double)(std::fabs(f_) + anorm) == anorm) break;
          g_ = w[i];
          h_ = dpythag(f_, g_);
          w[i] = h_;
          h_ = 1.0 / h_;
          c_ = g_ * h_;
          s_ = -f_ * h_;
          for (j = 1; j <= m; j++) {
            y_ = a[j][nm];
            z_ = a[j][i];
            a[j][nm] = y_ * c_ + z_ * s_;
            a[j][i] = z_ * c_ - y_ * s_;
          }
        }
      }
      z_ = w[k];
      if (l == k) {
        if (z_ < 0.0) {
          w[k] = -z_;
          for (j = 1; j <= n; j++) v[j][k] = -v[j][k];
        }
        break;
      }
      if (its == 100) {
        free_dvector(rv1, 1, n);
        return false;
      }
      x_ = w[l];
      nm = k - 1;
      y_ = w[nm];
      g_ = rv1[nm];
      h_ = rv1[k];
      f_ = ((y_ - z_) * (y_ + z_) + (g_ - h_) * (g_ + h_)) / (2.0 * h_ * y_);
      g_ = dpythag(f_, 1.0);
      f_ = ((x_ - z_) * (x_ + z_) + h_ * ((y_ / (f_ + SIGN(g_, f_))) - h_)) / x_;
      c_ = s_ = 1.0;
      for (j = l; j <= nm; j++) {
        i = j + 1;
        g_ = rv1[i];
        y_ = w[i];
        h_ = s_ * g_;
        g_ = c_ * g_;
        z_ = dpythag(f_, h_);
        rv1[j] = z_;
        c_ = f_ / z_;
        s_ = h_ / z_;
        f_ = x_ * c_ + g_ * s_;
        g_ = g_ * c_ - x_ * s_;
        h_ = y_ * s_;
        y_ *= c_;
        for (jj = 1; jj <= n; jj++) {
          x_ = v[jj][j];
          z_ = v[jj][i];
          v[jj][j] = x_ * c_ + z_ * s_;
          v[jj][i] = z_ * c_ - x_ * s_;
        }
        z_ = dpythag(f_, h_);
        w[j] = z_;
        if (z_) {
          z_ = 1.0 / z_;
          c_ = f_ * z_;
          s_ = h_ * z_;
        }
        f_ = c_ * g_ + s_ * y_;
        x_ = c_ * y_ - s_ * g_;
        for (jj = 1; jj <= m; jj++) {
          y_ = a[jj][j];
          z_ = a[jj][i];
          a[jj][j] = y_ * c_ + z_ * s_;
          a[jj][i] = z_ * c_ - y_ * s_;
        }
      }
      rv1[l] = 0.0;
      rv1[k] = f_;
      w[k] = x_;
    }
  }
  free_dvector(rv1, 1, n);
  return true;
}

static void dsvbksb(double **u, double w[], double **v, int m, int n, double b[], double x[])
{
  int jj, j, i;
  double s_, *tmp;

  tmp = dvector(1, n);
  for (j = 1; j <= n; j++) {
    s_ = 0.0;
    if (w[j]) {
      for (i = 1; i <= m; i++) s_ += u[i][j] * b[i];
      s_ /= w[j];
    }
    tmp[j] = s_;
  }
  for (j = 1; j <= n; j++) {
    s_ = 0.0;
    for (jj = 1; jj <= n; jj++) s_ += v[j][jj] * tmp[jj];
    x[j] = s_;
  }
  free_dvector(tmp, 1, n);
}

// -----------------------------------------------------------------------------
// Solve.cpp (SVD wrapper)
// -----------------------------------------------------------------------------

/**
 * Solve a linear least-squares system using SVD (Numerical Recipes style).
 *
 * This is used inside the Newton iteration to solve:
 *     J Δz = -rhs
 * where J is the analytic Jacobian matrix.
 *
 * SVD-based solving is more robust than naive Gaussian elimination for
 * potentially ill-conditioned spectral/collocation systems near limiting waves.
 */
static void Solve(double **a, double *b, int m, int n_, double *solution, int MP, int NP)
{
  int i;
  double *w, **v;
  double wmax, wmin;

  (void)MP;

  w = dvector(1, NP);
  v = dmatrix(1, NP, 1, NP);

  std::vector<double> scales(static_cast<std::size_t>(n_ + 1), 1.0);
  for (int j = 1; j <= n_; ++j) {
    double maximum = 0.0;
    for (int row = 1; row <= m; ++row)
      maximum = std::max(maximum, std::abs(a[row][j]));
    if (maximum > 0.0) scales[j] = maximum;
    for (int row = 1; row <= m; ++row) a[row][j] /= scales[j];
  }
  // Perform decomposition
  if (!dsvdcmp(a, m, n_, w, v)) {
    jacobian_condition = std::numeric_limits<double>::infinity();
    for (int j = 1; j <= n_; ++j) solution[j] = std::numeric_limits<double>::quiet_NaN();
    free_dmatrix(v, 1, NP, 1, NP);
    free_dvector(w, 1, NP);
    return;
  }

  /*
    SINGULAR VALUE FILTERING
    ------------------------
    After decomposition, w[i] are the singular values. Very small singular values
    correspond to directions in parameter space where the Jacobian provides little
    information (nearly dependent equations / near-singular J).

    If we were to divide by extremely small w[i] during back-substitution, the
    computed correction Δz could be dominated by numerical noise, leading to:
      - erratic Newton steps,
      - non-physical oscillations in the free-surface degrees of freedom,
      - failure to converge near limiting waves.

    The conventional NR remedy is to zero singular values below a relative cutoff.
    Here we use:
        wmin = wmax * 1e-15
    This preserves the well-resolved modes while damping the ill-conditioned ones.

    Note: This is not a "tunable physics parameter"; it is a numerical-stability
    safeguard. Adjust only with regression testing.
  */

  // Set up: see p65 of Press et al.
  wmax = 0.0;
  for (i = 1; i <= n_; i++) if (w[i] >= wmax) wmax = w[i];
  double smallest = std::numeric_limits<double>::infinity();
  for (i = 1; i <= n_; ++i) smallest = std::min(smallest, w[i]);
  jacobian_condition = smallest > 0.0 ? wmax / smallest
                                      : std::numeric_limits<double>::infinity();
  wmin = wmax * 1.e-15;
  for (i = 1; i <= n_; i++) if (w[i] <= wmin) w[i] = 0.0;

  // Back substitute
  dsvbksb(a, w, v, m, n_, b, solution);
  for (int j = 1; j <= n_; ++j) solution[j] /= scales[j];

  free_dmatrix(v, 1, NP, 1, NP);
  free_dvector(w, 1, NP);
}

// -----------------------------------------------------------------------------
// Subroutines.cpp (init, Eqns, Newton)
// -----------------------------------------------------------------------------

/**
 * Initialise the unknown vector z[] for the first (smallest-height) step.
 *
 * This builds a *linear* or weakly-nonlinear starting guess for the Newton solve:
 *   - sets the basic dispersion/period relations,
 *   - sets initial currents depending on the selected criterion,
 *   - initialises the free-surface collocation elevations z[10+m] over m=0..N
 *     as a cosine profile consistent with small-amplitude theory,
 *   - initialises the Fourier coefficients for the stream-function expansion.
 *
 * The overall collocation/Newton structure follows the method described in
 * Fenton (1988) and the formulation in Fenton (1999), Section 3.1.2.
 */
static void init()
{
  int i;
  double sigma;

  Is_finite {
    /*
      FINITE-DEPTH LINEAR START (continuation step 1)
      ---------------------------------------------
      We require an initial guess for the dimensionless mean depth kd.

      Linear dispersion for finite depth is (Fenton 1999, derived from Eq. 3.12 +
      linear ω² = gk tanh(kd)):

          ω² = g k tanh(kd)
          c = ω/k = sqrt(g/k * tanh(kd))

      Hence the non-dimensional wave celerity in k-scaling is:

          c*sqrt(k/g) = sqrt(tanh(kd))

      This code stores that as z[4] (see below) and sets z[3] such that
      z[3] = 2π / z[4], which corresponds to a period-like parameter in the same
      k-based scaling.

      For the Period-specified case (Case="Period"), kd is not known a priori.
      The original driver uses a smooth approximation (legacy, retained verbatim)
      to obtain a reasonable starting kd before Newton’s method enforces the full
      nonlinear system (Fenton 1988; Fenton 1999 §3.1.2).
    */
    iff(Case, Period) {
      sigma = 2. * pi * std::sqrt(height / Hoverd);

      // Period-input kd initialization only.
      // The final wavelength is not obtained from linear Airy theory; Newton's
      // method solves the full Fenton residual system.  This approximation only
      // supplies a finite-depth starting value for kd from the depth-based
      // frequency parameter alpha=sigma^2, where sigma=2π/[T√(g/d)].

      const double alpha = sigma * sigma;
      if (alpha <= 0.0) {
        z[1] = 0.0;
      } else {
        const double exp_arg = alpha * std::log(6.0 / 5.0);
        const double scale = std::exp(std::min(exp_arg, 20.0));
        const double tanh_arg = scale * std::sqrt(alpha);

        z[1] = alpha / std::tanh(tanh_arg);
      }

    } else {
      z[1] = 2. * pi * height / Hoverd;
    }

    z[2] = z[1] * Hoverd;
    z[4] = std::sqrt(std::tanh(z[1]));
    z[3] = 2. * pi / z[4];
  }

  Is_deep {
    z[1] = 1.;
    iff(Case, Period)
      z[2] = 4 * pi * pi * height;
    else
      z[2] = 2 * pi * height;
    z[3] = 2 * pi;
    z[4] = 1.;
  }

  z[7] = z[4];
  Is_finite {
    z[5] = z[6] = Current * std::sqrt(z[1]);
    z[4] = z[7] + z[5];
    z[3] = 2.0 * pi / z[4];
  }
  Is_deep { z[5] = z[6] = 0.0; }
  z[8] = 0.;
  z[9] = 0.5 * z[7] * z[7];
  /*
    TRIGONOMETRIC TABLES FOR COLLOCATION
    -----------------------------------
    The collocation points are located at:
        X_m = m π / N   for m = 0..N
    so that half a wavelength (crest to trough) is covered. The solver needs
    repeated access to cos(j X_m) and sin(j X_m) for modes j=1..N.

    Instead of recomputing cos/sin in the inner loops, we tabulate:
        cosa[p] = cos(p π / N)
        sina[p] = sin(p π / N)
    for p = 0..2N, because cos/sin are 2π-periodic and jX_m spans [0, 2π].

    Index mapping:
        p = (m*j) mod (2N)
    which yields:
        cos(j X_m) = cosa[p]
        sin(j X_m) = sina[p]
  */
  cosa[0] = 1.;
  sina[0] = 0.;
  /*
    INITIAL FREE-SURFACE SHAPE (COLLOCATED kη VALUES)
    -------------------------------------------------
    The solver stores the free surface as discrete values kη_m at the collocation
    phases X_m. A small-amplitude wave has approximately:
        η(ξ) ≈ (H/2) cos(kξ)   (ξ = x - c t)
    therefore:
        kη(X) ≈ (kH/2) cos(X) = (z[2]/2) cos(X)

    We initialise:
        z[10 + m] = (z[2]/2) cos(X_m)
    with X_m = m π / N. This provides a smooth starting shape for Newton.

    Note the special handling:
      - z[10]      corresponds to m=0 (crest, cos=1)
      - z[n+10]    corresponds to m=N (trough, cos=cos(π)=-1)
      - intermediate points follow the cosine profile.

    The Fourier coefficients B_j (stream-function expansion) start at zero; the
    nonlinear solve will populate them.
  */
  z[10] = 0.5 * z[2];
  for (i = 1; i <= n; i++) {
    cosa[i] = std::cos(i * pi / n);
    cosa[i + n] = std::cos((i + n) * pi / n);
    sina[i] = std::sin(i * pi / n);
    sina[i + n] = std::sin((i + n) * pi / n);
    z[n + i + 10] = 0.;
    z[i + 10] = 0.5 * z[2] * cosa[i];
  }
  z[n + 11] = 0.5 * z[2] / z[7];

  // Initial solution (zero height for extrapolation)
  for (i = 1; i <= 9; i++) sol[i][1] = z[i];
  for (i = 10; i <= num; i++) sol[i][1] = 0.;
}

/**
 * Assemble the nonlinear residual vector rhs[] for the collocation system.
 *
 * Unknown vector layout (legacy indexing, matching Fenton-style program listings):
 * z[1]   : kd        finite depth kd, in deep water a dummy normalisation (=1)
 * z[2]   : kH        wave height in k-based scaling
 * z[3]   : T√(gk)    k-based dimensionless apparent period
 * z[4]   : c √(k/g)  wave speed parameter in k-based scaling
 * z[5]   : u1 √(k/g) Eulerian mean current in k-based scaling
 * z[6]   : u2 √(k/g) Stokes/mass-transport current in k-based scaling
 * z[7]   : U √(k/g)  mean fluid speed in the moving frame (k-based)
 * z[8]   : q√(k^3/g) wave-volume-flux variable, q = Ubar*d - Q
 * z[9]   : r k/g     Bernoulli offset in k-based scaling
 *   z[10+m], m=0..N  = kη_m surface elevations at collocation points
 *   z[n+10+j], j=1..N = B_j Fourier coefficients in ψ expansion (Fenton 1999 Eq. 3.5)
 *
 * The core of this routine enforces the two free-surface residual
 * blocks at N+1 collocation points over half a wave (crest-to-trough):
 *
 *   (1) ψ_Fourier(X_m,ζ_m) - z[8] - z[7]ζ_m = 0
 *   (2) ½[(-z[7]+u_m)^2 + v_m^2] + ζ_m - z[9] = 0
 *
 * A Newton solver (see Newton()) drives rhs[] → 0.
 *
 * Returns: Sum of squares of residual components (useful as a diagnostic).
 */
static void VerticalBasis(int mode, double eta, double depth,
                          double& stream, double& velocity)
{
  const double forward = std::exp(mode * eta);
  if (std::strcmp(Depth, "Finite") != 0) {
    stream = velocity = forward;
    return;
  }
  const double reflected = std::exp(-mode * (2.0 * depth + eta));
  const double denominator = 1.0 + std::exp(-2.0 * mode * depth);
  stream = (forward - reflected) / denominator;
  velocity = (forward + reflected) / denominator;
}

static double EqnsAt(const double* state, double* rhs, double* tanh_values)
{
  const bool finite_depth = std::strcmp(Depth, "Finite") == 0;
  const bool wavelength_case = std::strcmp(Case, "Wavelength") == 0;

  rhs[1] = finite_depth ? state[2] - state[1] * Hoverd : state[1] - 1.0;
  rhs[2] = wavelength_case
             ? state[2] - 2.0 * pi * height
             : state[2] - height * state[3] * state[3];
  rhs[3] = state[4] * state[3] - pi - pi;
  rhs[4] = state[5] + state[7] - state[4];
  rhs[5] = state[1] * (state[6] + state[7] - state[4]) - state[8];

  if (finite_depth) {
    for (int j = 1; j <= n; ++j) {
      tanh_values[j] = std::tanh(static_cast<double>(j) * state[1]);
    }
    rhs[6] = state[Current_criterion + 4] - Current * std::sqrt(state[1]);
  } else {
    rhs[6] = wavelength_case
               ? state[Current_criterion + 4] - Current
               : state[Current_criterion + 4] - Current * state[3];
  }

  rhs[7] = state[10] + state[n + 10];
  for (int i = 1; i <= n - 1; ++i) {
    rhs[7] = rhs[7] + state[10 + i] + state[10 + i];
  }
  rhs[8] = state[10] - state[n + 10] - state[2];

  /*
    Each collocation point is independent once the state vector and finite-depth
    tanh factors are fixed.  The loop therefore writes two unique residual entries
    per iteration. The analytic Jacobian uses the same vertical basis and
    collocation phases as this residual evaluation.
  */
  for (int m = 0; m <= n; ++m) {
    double psi_sum = 0.0;
    double u_sum = 0.0;
    double v_sum = 0.0;
    const double eta = state[10 + m];

    for (int j = 1; j <= n; ++j) {
      const int nm = (m * j) % (2 * n);
      const double coefficient = state[n + j + 10];
      double stream_basis, velocity_basis;
      VerticalBasis(j, eta, state[1], stream_basis, velocity_basis);
      const double weighted_mode = static_cast<double>(j) * coefficient;
      psi_sum += coefficient * stream_basis * cosa[nm];
      u_sum += weighted_mode * velocity_basis * cosa[nm];
      v_sum += weighted_mode * stream_basis * sina[nm];
    }

    rhs[m + 9] = psi_sum - state[8] - state[7] * eta;
    rhs[n + m + 10] = 0.5 * (std::pow((-state[7] + u_sum), 2.0) + v_sum * v_sum)
                      + eta - state[9];
  }

  double residual_sum = 0.0;
  for (int j = 1; j <= num; ++j) residual_sum += rhs[j] * rhs[j];
  return residual_sum;
}

static double Eqns(double* rhs)
{
  return EqnsAt(z, rhs, Tanh);
}

// Analytic derivatives of all global constraints and surface equations.
static void Jacobian(double** jac)
{
  for (int row = 1; row <= num; ++row)
    for (int col = 1; col <= num; ++col) jac[row][col] = 0.0;
  const bool finite = std::strcmp(Depth, "Finite") == 0;
  const bool wavelength = std::strcmp(Case, "Wavelength") == 0;
  if (finite) { jac[1][1] = -Hoverd; jac[1][2] = 1.0; }
  else jac[1][1] = 1.0;
  jac[2][2] = 1.0;
  if (!wavelength) jac[2][3] = -2.0 * height * z[3];
  jac[3][3] = z[4]; jac[3][4] = z[3];
  jac[4][5] = jac[4][7] = 1.0; jac[4][4] = -1.0;
  jac[5][1] = z[6] + z[7] - z[4];
  jac[5][6] = jac[5][7] = z[1]; jac[5][4] = -z[1]; jac[5][8] = -1.0;
  jac[6][Current_criterion + 4] = 1.0;
  if (finite) jac[6][1] = -0.5 * Current / std::sqrt(z[1]);
  else if (!wavelength) jac[6][3] = -Current;
  jac[7][10] = jac[7][n + 10] = 1.0;
  for (int m = 1; m < n; ++m) jac[7][10 + m] = 2.0;
  jac[8][10] = 1.0; jac[8][n + 10] = jac[8][2] = -1.0;
#if defined(_OPENMP)
#pragma omp parallel for schedule(static) if(n >= 24)
#endif
  for (int m = 0; m <= n; ++m) {
    const double eta = z[10 + m];
    double us = 0.0, vs = 0.0, uy_sum = 0.0, vy_sum = 0.0;
    double pd = 0.0, ud = 0.0, vd = 0.0;
    for (int j = 1; j <= n; ++j) {
      const int phase = (m * j) % (2 * n);
      double stream, velocity;
      VerticalBasis(j, eta, z[1], stream, velocity);
      const double coeff = z[n + 10 + j], co = cosa[phase], si = sina[phase];
      us += j * coeff * velocity * co;
      vs += j * coeff * stream * si;
      uy_sum += j * j * coeff * stream * co;
      vy_sum += j * j * coeff * velocity * si;
      if (finite) {
        const double e = std::exp(-2.0 * j * z[1]);
        const double forward = std::exp(j * eta);
        const double reflected = std::exp(-j * (2.0 * z[1] + eta));
        const double factor = 2.0 * j / ((1.0 + e) * (1.0 + e));
        const double sd = factor * (reflected + e * forward);
        const double cd = factor * (e * forward - reflected);
        pd += coeff * sd * co;
        ud += j * coeff * cd * co;
        vd += j * coeff * sd * si;
      }
      jac[m + 9][n + 10 + j] = stream * co;
    }
    const double relative_u = us - z[7];
    jac[m + 9][1] = pd;
    jac[m + 9][7] = -eta;
    jac[m + 9][8] = -1.0;
    jac[m + 9][10 + m] = relative_u;
    const int row = n + m + 10;
    jac[row][1] = relative_u * ud + vs * vd;
    jac[row][7] = -relative_u;
    jac[row][9] = -1.0;
    jac[row][10 + m] = 1.0 + relative_u * uy_sum + vs * vy_sum;
    for (int j = 1; j <= n; ++j) {
      const int phase = (m * j) % (2 * n);
      double stream, velocity;
      VerticalBasis(j, eta, z[1], stream, velocity);
      jac[row][n + 10 + j] = j * (relative_u * velocity * cosa[phase] + vs * stream * sina[phase]);
    }
  }
}

static bool Admissible()
{
  for (int i = 1; i <= num; ++i) if (!std::isfinite(z[i])) return false;
  if (z[1] <= 0.0 || z[2] <= 0.0 || z[3] <= 0.0 || z[4] <= 0.0 || z[7] <= 0.0) return false;
  for (int m = 0; m <= n; ++m) {
    if (std::strcmp(Depth, "Finite") == 0 && z[1] + z[10 + m] <= 0.0) return false;
    if (n * std::abs(z[10 + m]) > 650.0) return false;
  }
  return true;
}

static double MaxResidual(const double* residual)
{
  double maximum = 0.0;
  for (int i = 1; i <= num; ++i) {
    if (!std::isfinite(residual[i])) return std::numeric_limits<double>::infinity();
    maximum = std::max(maximum, std::abs(residual[i]));
  }
  return maximum;
}

// Backtracking accepts only admissible states with a decreasing residual norm.
static bool Newton()
{
  if (!Admissible()) return false;
  double* rhs = dvector(1, num);
  double* correction = dvector(1, num);
  double** jacobian = dmatrix(1, num, 1, num);
  std::vector<double> base(static_cast<std::size_t>(num + 1), 0.0);
  for (int i = 1; i <= num; ++i) base[i] = z[i];
  const double initial_merit = Eqns(rhs1);
  for (int i = 1; i <= num; ++i) rhs[i] = -rhs1[i];
  Jacobian(jacobian);
  bool finite = std::isfinite(initial_merit);
  for (int i = 1; i <= num; ++i)
    for (int j = 1; j <= num; ++j) finite = finite && std::isfinite(jacobian[i][j]);
  bool accepted = false;
  if (finite) {
    Solve(jacobian, rhs, num, num, correction, num, num);
    for (double damping = 1.0; damping >= 1.0 / 1048576.0; damping *= 0.5) {
      for (int i = 1; i <= num; ++i) z[i] = base[i] + damping * correction[i];
      if (!Admissible()) continue;
      const double merit = Eqns(rhs1);
      if (std::isfinite(merit) &&
          (merit < initial_merit * (1.0 - 1.0e-4 * damping) || MaxResidual(rhs1) <= kFinalResidualTolerance)) {
        accepted = true;
        break;
      }
    }
  }
  if (!accepted) for (int i = 1; i <= num; ++i) z[i] = base[i];
  free_dvector(rhs, 1, num);
  free_dvector(correction, 1, num);
  free_dmatrix(jacobian, 1, num, 1, num);
  return accepted;
}

static void RefreshCondition()
{
  jacobian_condition = std::numeric_limits<double>::infinity();
  if (!Admissible()) return;
  double** jac = dmatrix(1, num, 1, num);
  double* zero = dvector(1, num);
  double* correction = dvector(1, num);
  Jacobian(jac);
  bool finite = true;
  for (int i = 1; i <= num; ++i) {
    zero[i] = 0.0;
    for (int j = 1; j <= num; ++j) finite = finite && std::isfinite(jac[i][j]);
  }
  if (finite) Solve(jac, zero, num, num, correction, num, num);
  free_dmatrix(jac, 1, num, 1, num);
  free_dvector(zero, 1, num);
  free_dvector(correction, 1, num);
}

static bool ContinueToTarget()
{
  double fraction = 0.0, previous_fraction = 0.0;
  double increment = std::min(0.25, 1.0 / nstep);
  const double minimum_increment = 1.0 / 4096.0;
  std::vector<double> accepted(static_cast<std::size_t>(num + 1), 0.0), previous = accepted;
  bool have_previous = false;
  final_residual = std::numeric_limits<double>::infinity();
  for (int trial = 0; trial < 160 && fraction < 1.0; ++trial) {
    const double target = std::min(1.0, fraction + increment);
    height = target * Height;
    Hoverd = target * MaxH;
    if (fraction == 0.0) init();
    else {
      const double ratio = have_previous ? (target - fraction) / (fraction - previous_fraction) : 0.0;
      for (int i = 1; i <= num; ++i) z[i] = accepted[i] + ratio * (accepted[i] - previous[i]);
      if (!Admissible()) for (int i = 1; i <= num; ++i) z[i] = accepted[i];
    }
    bool converged = false;
    int iterations = 0;
    const double tolerance = target == 1.0 ? kFinalResidualTolerance : std::min(crit, 1.0e-10);
    for (; iterations <= number; ++iterations) {
      if (!Admissible()) break;
      Eqns(rhs1);
      final_residual = MaxResidual(rhs1);
      if (final_residual <= tolerance) { converged = true; break; }
      if (iterations == number || !Newton()) break;
    }
    std::fprintf(monitor, "Height fraction %.8f, iterations %d, max|F| %.3e: %s\n",
                 target, iterations, final_residual, converged ? "PASS" : "retry");
    if (converged) {
      previous = accepted;
      previous_fraction = fraction;
      have_previous = fraction > 0.0;
      for (int i = 1; i <= num; ++i) accepted[i] = z[i];
      fraction = target;
      increment = std::min(0.25, increment * (iterations <= 8 ? 1.5 : 1.0));
    } else {
      if (fraction > 0.0) {
        for (int i = 1; i <= num; ++i) z[i] = accepted[i];
        height = fraction * Height;
        Hoverd = fraction * MaxH;
      }
      increment *= 0.5;
      if (increment < minimum_increment) return false;
    }
  }
  if (fraction != 1.0 || !Admissible()) return false;
  height = Height;
  Hoverd = MaxH;
  Eqns(rhs1);
  final_residual = MaxResidual(rhs1);
  return final_residual <= kFinalResidualTolerance;
}

// Independent free-surface boundary residuals between collocation nodes.
static double SpectralDefect()
{
  double maximum = 0.0;
  const int samples = std::max(2048, 32 * n);
  for (int m = 0; m < samples; ++m) {
    const double phase = pi * (m + 0.5) / samples;
    double eta = 0.5 * Y[n] * std::cos(n * phase);
    for (int j = 1; j < n; ++j) eta += Y[j] * std::cos(j * phase);
    if (std::strcmp(Depth, "Finite") == 0 && z[1] + eta <= 0.0)
      return std::numeric_limits<double>::infinity();
    double ps = 0.0, us = 0.0, vs = 0.0;
    for (int j = 1; j <= n; ++j) {
      double stream, velocity;
      VerticalBasis(j, eta, z[1], stream, velocity);
      ps += B[j] * stream * std::cos(j * phase);
      us += j * B[j] * velocity * std::cos(j * phase);
      vs += j * B[j] * stream * std::sin(j * phase);
    }
    const double kin = ps - z[8] - z[7] * eta;
    const double dyn = 0.5 * ((us - z[7]) * (us - z[7]) + vs * vs) + eta - z[9];
    if (!std::isfinite(kin) || !std::isfinite(dyn)) return std::numeric_limits<double>::infinity();
    maximum = std::max({maximum, std::abs(kin), std::abs(dyn)});
  }
  return maximum;
}

// -----------------------------------------------------------------------------
// Inout.cpp (I/O, reporting, kinematics)
// -----------------------------------------------------------------------------

static void Title_block(std::FILE* file);
static void Title_block_deep(std::FILE* file);

static void Input_Data_block(std::FILE* file)
{
  std::fprintf(file, "# %s", Title);
  std::fprintf(file, "\n\n# Printing input data here to check");
  std::fprintf(file, "\n\n# Height/Depth:%6.8f", MaxH);
  iff(Case, Wavelength) {
    std::fprintf(file, "\n# Length/Depth:%7.8f", L);
  }
  iff(Case, Period) {
    std::fprintf(file, "\n# Dimensionless Period T*sqrt(g/d):%8.8f", T);
  }
  std::fprintf(file, "\n# Current criterion: %s,  Dimensionless value:%6.8lf", Currentname, Current);

  if (std::strcmp(Theory, "Stokes") == 0) {
    if (n <= 5) std::sprintf(Method, "\n# Solution by %d-order Stokes theory", n);
    else {
      n = 5;
      std::sprintf(Method, "\n# Solution by %d-order Stokes theory", n);
      std::printf("\n\n# (A value of N > 5 has been specified for the Stokes theory.");
      std::printf("\n# I do not have a theory for that. The program has set N = 5)");
    }
  }
  if (std::strcmp(Theory, "Fourier") == 0)
    std::sprintf(Method, "\n# Solution by %d-term Fourier series", n);

  std::fprintf(file, "\n\n%s\n", Method);
}

/**
 * Read one problem case from data.dat and associated control files.
 *
 * data.dat may contain multiple cases; a title line equal to "FINISH" ends the run.
 *
 * This routine:
 *   - reads wave height and chooses finite vs deep regime (legacy sign convention),
 *   - reads Case: "Wavelength" or "Period",
 *   - reads current criterion and value,
 *   - reads N and nstep (spectral order and height stepping),
 *   - reads embedded convergence controls (formerly Convergence.dat),
 *   - sets embedded output sampling controls (Points.dat no longer required).
 *
 * Returns:
 *   1 if a case was read, 0 if FINISH was encountered.
 */
static int Read_data(void)
{
  if (!std::fgets(Title, sizeof(Title), Input1)) return 0;
  Title[std::strcspn(Title, "\r\n")] = '\0';
  if (std::strcmp(Title, "FINISH") == 0) return 0;
  auto read_double = [](double& value) {
    return std::fgets(dummy, sizeof(dummy), Input1) &&
           std::sscanf(dummy, "%lf", &value) == 1;
  };
  auto read_int = [](int& value) {
    return std::fgets(dummy, sizeof(dummy), Input1) &&
           std::sscanf(dummy, "%d", &value) == 1;
  };
  double value = 0.0;
  if (!read_double(MaxH) || !std::fgets(dummy, sizeof(dummy), Input1) ||
      std::sscanf(dummy, "%19s", Case) != 1 || !read_double(value) ||
      !read_int(Current_criterion) || !read_double(Current) ||
      !read_int(n) || !read_int(nstep)) return 0;
  if (MaxH >= 0.0) {
    std::strcpy(Depth, "Finite");
    if (!std::isfinite(value) || value <= 0.0) return 0;
    if (std::strcmp(Case, "Wavelength") == 0) { L = value; Height = MaxH / L; }
    else if (std::strcmp(Case, "Period") == 0) { T = value; Height = MaxH / (T * T); }
    else return 0;
  } else {
    std::strcpy(Depth, "Deep");
    Height = -MaxH;
  }
  if (Current_criterion == 1) std::strcpy(Currentname, Current1);
  if (Current_criterion == 2) std::strcpy(Currentname, Current2);
  number = kConvergenceMaxIterations;
  crit = kConvergenceCriterion;
  Surface_points = kSurfacePointsDefault;
  Nprofiles = kNprofilesDefault;
  Points = kProfilePointsDefault;
  Input_Data_block(monitor);
  return 1;
}

static void Title_block(std::FILE* file)
{
  // Highest-wave reference curve:
  //   Hm/d as a function of λ/d is fitted by a rational approximation to
  //   Williams (1981) highest-wave computations, presented as Eq. (32) in
  //   Fenton (1990) “Nonlinear wave theories” (see the bibliography in the
  //   header comment). This is used here only as a convenient benchmark in
  //   the printed headers (percentage of highest possible for this λ/d).
  L = 2 * pi / z[1];
  Highest = (0.0077829 * L * L * L + 0.0095721 * L * L + 0.141063 * L)
          / (0.0093407 * L * L * L + 0.0317567 * L * L + 0.078834 * L + 1);
  std::fprintf(file, "# %s", Title);
  std::fprintf(file, "\n\n%s\n", Method);
  std::fprintf(file, "\n# Height/Depth:%6.3f, %3.0lf%% of the maximum of H/d =%6.3f for this length:",
               z[2] / z[1], z[2] / z[1] / Highest * 100., Highest);
  std::fprintf(file, "\n# Length/Depth:%7.2f", 2 * pi / z[1]);
  std::fprintf(file, "\n# Dimensionless Period T*sqrt(g/d):%7.2f", z[3] / std::sqrt(z[1]));
  std::fprintf(file, "\n# Current criterion: %s,  Dimensionless value:%6.3lf\n", Currentname, Current);
}

static void Title_block_deep(std::FILE* file)
{
  Highest = 0.141063;
  std::fprintf(file, "# %s", Title);
  std::fprintf(file, "\n\n%s\n", Method);
  std::fprintf(file, "\n# Height/Length:%6.5f, %3.0lf%% of the maximum of H/L =%6.5f",
               z[2] / 2 / pi, (z[2] / 2 / pi) / Highest * 100., Highest);
  std::fprintf(file, "\n# Dimensionless Period T*sqrt(g/L):%7.5f", z[3] / std::sqrt(2 * pi));
  std::fprintf(file, "\n# Current criterion: %s,  Dimensionless value:%6.3lf\n", Currentname, Current);
}

static void Results(const char *Description, double x, double y)
{
  /*
    RESULTS LINE PRINTER (solution.res)

    The original multi-file program wrote two outputs:
      - solution.res       : human-readable table (this file)
      - Solution-Flat.res  : redundant “flat” tabular mirror

    The legacy redundant tabular file “Solution-Flat.res” is not produced by this build. This routine now writes only to solution.res, preserving formatting and numerical values.
  */
  std::fprintf(Solution, "\n%s" LO, Description, x);
  Is_finite {
    std::fprintf(Solution, "" LO, y);
  }
}

/**
 * Evaluate the free-surface elevation in k-based non-dimensional form: kη(X).
 *
 * After convergence, the collocation elevations z[10+m] are post-processed into
 * a cosine-series representation:
 *     kη(X) ≈ Σ_{j=1}^{N} E_j cos(j X)
 * which is used for smooth evaluation and output sampling.
 *
 * This matches the Fourier representation used in stream-function/Fourier methods
 * (see Fenton 1999 discussion leading to Eq. 3.7–3.8, where η(X) is unknown and
 * discretised/collocated).
 */
static double Surface(double x)
{
  int j;
  static double kEta;

  kEta = 0.;

  /*
    COSINE SERIES EVALUATION DETAILS
    -------------------------------
    The coefficients Y[j] are produced in main() by a discrete cosine transform-
    like reconstruction from the collocation values kη_m.

    The series is written in the "half-range" cosine form consistent with the
    collocation grid on [0, π]. For that reason, the highest mode (j=N) carries
    a 1/2 factor in the standard DCT-I normalisation.

    This is why the final term is:
        0.5 * Y[N] * cos(N X)
    while modes 1..N-1 have full weight.
  */
  for (j = 1; j < n; j++) kEta += Y[j] * std::cos(j * x);
  kEta += 0.5 * Y[n] * std::cos(n * x);
  return kEta;
}

/**
 * Evaluate kinematics and pressure at a single point (X, Y_).
 *
 * Inputs:
 *   X   : phase variable in the travelling-wave formulation; nondimensional X = k ξ with ξ = x - c t.
 *   Y_  : vertical coordinate in the solver’s k-based system.
 *         For finite depth, Point() expects Y_ = k (y - d), where y is physical
 *         elevation above the bed (y=0 at bed, y=d at mean free surface).
 *
 * Outputs (written into global variables for legacy compatibility):
 *   u, v              : horizontal/vertical velocities in the stationary frame
 *   ux, uy, vx, vy    : spatial derivatives used for accelerations
 *   dphidt            : ∂φ/∂t in the stationary frame (via steady-wave relations)
 *   Pressure          : dynamic pressure check from Bernoulli
 *   Bernoulli_check   : residual of Bernoulli constant (diagnostic)
 *
 * Theory:
 *   The stream-function Fourier series representation corresponds to Fenton (1999)
 *   Eq. (3.5) (finite depth) and Eq. (3.6) (deep limit), differentiated to obtain
 *   velocities (U,V) in the moving frame and then shifted to the stationary frame.
 */
static void Point(double X, double Y_)
{
  int j;
  double C = 0.0, S = 0.0, Sin = 0.0, Cos = 0.0, y = 0.0;

  u = v = ux = vx = phi = psi = 0.;

  /*
    VERTICAL COORDINATE NORMALISATION (FINITE DEPTH)
    ------------------------------------------------
    In finite depth we often need the *depth-based* non-dimensional elevation
    y/d (with y=0 at bed, y=d at mean free surface) for Bernoulli terms.

    Recall:
        Y_ = k (y - d)
        kd = k d

    Therefore:
        y/d = 1 + (Y_ / (k d)) = 1 + Y_/kd

    This `y` variable is thus y/d (dimensionless by depth).
  */
  Is_finite y = 1. + Y_ / kd;

  /*
    FOURIER EVALUATION OF POTENTIAL/STREAMFUNCTION AND DERIVATIVES
    --------------------------------------------------------------
    The converged solution provides coefficients B[j] for harmonic modes j=1..N.

    For each mode:
      - Horizontal dependence:
          sin(jX) and cos(jX)
      - Vertical dependence (finite depth):
          cosh(jY) + sinh(jY) tanh(jkd)   and
          sinh(jY) + cosh(jY) tanh(jkd)
        These combinations arise from satisfying the bottom boundary condition
        exactly while still satisfying Laplace's equation for each harmonic.

      - Vertical dependence (deep water):
          exp(jY) (the kd→∞ limit)

    Quantities accumulated:
      phi, psi:
        "potential" and streamfunction-like series values in the program's
        internal nondimensionalisation.
      u, v:
        derivatives used to obtain velocities (still in the moving-frame/k-scaling
        until converted below).
      ux, vx:
        spatial derivatives needed for accelerations and Bernoulli checks.

    Note: The naming here is inherited. The relationships to standard φ/ψ are
    described in Fenton (1999); the key point is that differentiation of the
    harmonic basis yields these velocity components consistently.
  */
  for (j = 1; j <= n; j++) {
    Cos = std::cos(j * X);
    Sin = std::sin(j * X);
    VerticalBasis(j, Y_, kd, S, C);
    phi += B[j] * C * Sin;
    psi += B[j] * S * Cos;
    u += j * B[j] * C * Cos;
    v += j * B[j] * S * Sin;
    ux += -j * j * B[j] * C * Sin;
    vx += j * j * B[j] * S * Cos;
  }

  Is_finite {
    /*
      SCALE CONVERSION: k-BASED -> DEPTH-BASED (FINITE DEPTH)
      -------------------------------------------------------
      Up to this point, the harmonic sums have been evaluated in the solver's
      native k-based nondimensionalisation (g & k).

      For reporting (and for the Bernoulli terms that use y/d), we convert to
      the depth-based nondimensionalisation (g & d).

      Using kd = k d:

        φ_(gd)   = φ_(gk) / (kd)^(3/2)
        ψ_(gd)   = ψ_(gk) / (kd)^(3/2)
        u_(gd)   = u_(gk) / √(kd)
        v_(gd)   = v_(gk) / √(kd)

      Spatial derivatives transform accordingly; see the factors below.

      After converting, we add the Eulerian mean current ce (depth-based) to get
      stationary-frame absolute velocities, and we form ∂φ/∂t and accelerations
      using the steady-wave identity ∂/∂t = -c ∂/∂x.
    */
    phi /= std::pow(kd, 1.5);
    psi /= std::pow(kd, 1.5);
    u /= std::pow(kd, 0.5);
    v /= std::pow(kd, 0.5);
    ux *= std::pow(kd, 0.5);
    vx *= std::pow(kd, 0.5);
    u = ce + u;
    phi = ce * X + phi;
    psi = ce * y + psi;
    /*
      TIME DERIVATIVE OF POTENTIAL (STEADY WAVE RELATION)
      ---------------------------------------------------
      In a steady travelling wave, fields depend on (x - c t). Therefore:
          ∂/∂t = -c ∂/∂x
      and since u = ∂φ/∂x, we have:
          ∂φ/∂t = -c u

      This is the quantity labelled dphidt in the output tables.
    */
    dphidt = -c * u;

    ut = -c * ux;
    vt = -c * vx;
    uy = vx;
    vy = -ux;
    /*
      MATERIAL (LAGRANGIAN) ACCELERATIONS
      ----------------------------------
      The local (Eulerian) time derivatives are:
          ut = ∂u/∂t
          vt = ∂v/∂t
      and the convective terms are (u·∇)u, (u·∇)v.

      In 2D:
          Du/Dt = ut + u ux + v uy
          Dv/Dt = vt + u vx + v vy

      These appear in Morison-type force calculations and are written to
      flowfield.res.
    */
    dudt = ut + u * ux + v * uy;
    dvdt = vt + u * vx + v * vy;
    /*
      PRESSURE AND BERNOULLI CONSISTENCY CHECK
      ---------------------------------------
      Dynamic free-surface condition enforces a Bernoulli constant R along the
      free surface. Away from the free surface, the same Bernoulli expression
      provides a convenient diagnostic:

          ∂φ/∂t + p/ρ + y + ½(u²+v²) = constant

      The program reports:

        Pressure:
          rearranged pressure term (p/ρ) from the Bernoulli constant expression,
          using the moving-frame velocity (u-c) where appropriate.

        Bernoulli_check:
          residual of the Bernoulli equation after substituting all computed terms.
          Ideally ~ 0 everywhere; deviations indicate truncation error or incomplete
          convergence.

      Note: y is y/d here (dimensionless by depth), consistent with the conversion
      applied above.
    */
    Pressure = R - y - 0.5 * ((u - c) * (u - c) + v * v);
    Bernoulli_check = dphidt + Pressure + y + 0.5 * (u * u + v * v) - (R - 0.5 * c * c);
  }

  Is_deep {
    u = z[5] + u;
    phi = z[5] * X + phi;
    dphidt = -z[4] * u;

    ut = -z[4] * ux;
    vt = -z[4] * vx;
    uy = vx;
    vy = -ux;
    dudt = ut + u * ux + v * uy;
    dvdt = vt + u * vx + v * vy;
    Pressure = z[9] - Y_ - 0.5 * ((u - z[4]) * (u - z[4]) + v * v);
    Bernoulli_check = dphidt + Pressure + Y_ + 0.5 * (u * u + v * v) - (z[9] - 0.5 * z[4] * z[4]);
  }
}

struct KinematicOutputs {
  double surface_velocity = 0.0, bed_velocity = 0.0;
  double surface_acceleration = 0.0, bed_acceleration = 0.0;
  double asymmetry = 0.0;
};
static KinematicOutputs kinematics;

static double PhaseQuantity(double phase, bool bed, bool acceleration)
{
  const double ordinate = bed ? -kd : Surface(phase);
  if (std::strcmp(Depth,"Finite")==0 && !bed && kd+ordinate<=0.0)
    throw std::runtime_error("Interpolated free surface intersects the bed.");
  Point(phase,ordinate);
  const double value = acceleration ? dudt : u;
  if (!std::isfinite(value)) throw std::runtime_error("Non-finite kinematic phase evaluation.");
  return std::abs(value);
}

  static double PhaseMaximum(bool bed, bool acceleration) {
    const int count = std::max(257, 8 * n + 1);
    std::vector<double> phases(count), values(count);
    double best = 0.0;
    for (int i = 0; i < count; ++i) {
      phases[i] = pi * i / (count - 1);
      values[i] = PhaseQuantity(phases[i], bed, acceleration);
      best = std::max(best, values[i]);
    }
    const double golden = (std::sqrt(5.0) - 1.0) / 2.0;
    for (int i = 1; i < count - 1; ++i) {
      if (!(values[i] >= values[i-1] && values[i] >= values[i+1] &&
            (values[i] > values[i-1] || values[i] > values[i+1]))) continue;
      double left = phases[i-1], right = phases[i+1];
      double x1 = right - golden*(right-left), x2 = left + golden*(right-left);
      double f1 = PhaseQuantity(x1, bed, acceleration);
      double f2 = PhaseQuantity(x2, bed, acceleration);
      best = std::max({best, f1, f2});
      for (int iteration = 0; iteration < 64 && right-left >= 1e-12; ++iteration) {
        if (f1 < f2) {
          left = x1; x1 = x2; f1 = f2;
          x2 = left + golden*(right-left);
          f2 = PhaseQuantity(x2, bed, acceleration);
        } else {
          right = x2; x2 = x1; f2 = f1;
          x1 = right - golden*(right-left);
          f1 = PhaseQuantity(x1, bed, acceleration);
        }
        best = std::max({best, f1, f2});
      }
    }
    return best;
  }


static void CalculateKinematicOutputs()
{
  kinematics = KinematicOutputs{};
  kinematics.surface_velocity = PhaseMaximum(false,false);
  kinematics.surface_acceleration = PhaseMaximum(false,true);
  if (std::strcmp(Depth,"Finite")==0) {
    kinematics.bed_velocity = PhaseMaximum(true,false);
    kinematics.bed_acceleration = PhaseMaximum(true,true);
    if (!std::isfinite(ub2) || ub2<0.0) throw std::runtime_error("Invalid bed orbital velocity variance.");
  }
  Point(0.0,Surface(0.0));
  const double crest_velocity = u;
  Point(pi,Surface(pi));
  if (!std::isfinite(crest_velocity) || !std::isfinite(u))
    throw std::runtime_error("Non-finite surface velocity for asymmetry.");
  kinematics.asymmetry = std::abs(u)>0.0 ? std::abs(crest_velocity/u) : 0.0;
  if (!std::isfinite(kinematics.asymmetry)) throw std::runtime_error("Non-finite velocity asymmetry.");
}

static void WriteKinematicOutputs(std::FILE* file)
{
  const bool finite = std::strcmp(Depth,"Finite")==0;
  const char* velocity_scale = finite ? "u/sqrt(g*d)" : "u/sqrt(g/k)";
  std::fprintf(file,"\n\n");
  std::fprintf(file,"# Material accelerations are divided by g. Bed quantities require finite depth.\n");
  auto row = [file](const char* label, double value, const char* scaling, bool available=true) {
    if (available) std::fprintf(file,"# %-38s % .16e  %s\n",label,value,scaling);
    else std::fprintf(file,"# %-38s N/A\n",label);
  };
  row("Surface max |u|",kinematics.surface_velocity,velocity_scale);
  row("Bed max |u|",kinematics.bed_velocity,velocity_scale,finite);
  row("Surface max |aₓ|",kinematics.surface_acceleration,"a_x/g");
  row("Bed max |aₓ|",kinematics.bed_acceleration,"a_x/g",finite);
  row("Velocity asymmetry |uc|/|ut|",kinematics.asymmetry,"-");
  row("Mean square bed orbital vel ub²",finite ? ub2/kd : 0.0,"ub²/(g*d)",finite);
  row("Bed orbital RMS velocity ub,rms",finite ? std::sqrt(ub2/kd) : 0.0,"ub,rms/sqrt(g*d)",finite);
}

static double MeanSquareBedOrbitalVelocity()
{
  if (std::strcmp(Depth, "Finite") != 0) return 0.0;
  double result = 0.0;
  for (int j = 1; j <= n; ++j) {
    const double inverse_cosh = 2.0 * std::exp(-j*z[1]) / (1.0 + std::exp(-2.0*j*z[1]));
    const double amplitude = j*B[j]*inverse_cosh;
    result += 0.5*amplitude*amplitude;
  }
  return result;
}

static void CalculateIntegralOutputs()
{
  kd = z[1];
  L = 2.0*pi/kd; H = z[2]/kd; T = z[3]/std::sqrt(kd);
  c = z[4]/std::sqrt(kd); ce = z[5]/std::sqrt(kd);
  cs = z[6]/std::sqrt(kd); ubar = z[7]/std::sqrt(kd);
  if (std::strcmp(Depth, "Finite") != 0) return;
  q = z[7]*kd-z[8]; Q = q/std::pow(kd,1.5);
  r = z[9]+kd; R = r/kd;
  pulse = z[8]+kd*z[5];
  ke = 0.5*(z[4]*pulse-z[5]*q);
  pe = Y[n]*Y[n]/16.0;
  for (int j=1;j<n;++j) pe += 0.25*Y[j]*Y[j];
  const double beta = 2.0*z[9]-z[4]*z[4]+2.0*z[4]*z[5];
  sxx = 4.0*ke-3.0*pe+beta*kd-2.0*z[5]*pulse;
  f = z[4]*(3.0*ke-2.0*pe)+0.5*beta*(pulse+z[4]*kd)-2.0*z[4]*z[5]*pulse;
  ub2 = MeanSquareBedOrbitalVelocity();
  s = sxx-2.0*z[4]*pulse+(z[4]*z[4]+0.5*kd)*kd;
}

static std::map<std::string,double> SelectedOutputs()
{
  // Values use g and k for nondimensionalisation; density is factored out.
  return {{"L",2.0*pi},{"c",z[4]},{"Impulse",pulse},{"KineticEnergy",ke},
          {"PotentialEnergy",pe},{"EnergyDensity",ke+pe},{"Sxx",sxx},{"Power",f},
          {"MeanSquareBedVelocity",ub2}};
}

static bool VerifyFiniteSystem() {
        using I=Verified::Interval;

        finite_system_certificate=Verified::Certificate{};
        volatile Real smallest=std::numeric_limits<Real>::denorm_min();
        if (!std::numeric_limits<Real>::is_iec559 || std::fegetround()!=FE_TONEAREST ||
            smallest*2.0!=std::nextafter(std::numeric_limits<Real>::denorm_min(),1.0))
            throw std::runtime_error("Interval verification requires IEEE binary64, round-to-nearest and gradual underflow.");
        std::vector<I> state;state.reserve(static_cast<size_t>(num+1));
        state.emplace_back(0.0);
        for (int i=1;i<=num;++i) state.emplace_back(z[i]);
        const bool finite = std::strcmp(Depth,"Finite")==0;
        const bool wavelength = std::strcmp(Case,"Wavelength")==0;
        std::vector<I> F((size_t)(2*n+10));
        F[0]=finite ? state[2]-state[1]*I(MaxH) : state[1]-I(1);
        F[1]=wavelength ? state[2]-I(2)*Verified::pi()*I(Height)
                        : state[2]-I(Height)*state[3].square();
        F[2]=state[4]*state[3]-I(2)*Verified::pi();
        F[3]=state[5]+state[7]-state[4];
        F[4]=state[1]*(state[6]+state[7]-state[4])-state[8];
        F[5]=state[(size_t)Current_criterion+4]-I(Current)*(finite ? state[1].sqrt()
                           : (wavelength ? I(1) : state[3]));
        F[6]=state[10]+state[(size_t)n+10];
        for (int j=11;j<n+10;++j) F[6]=F[6]+I(2)*state[(size_t)j];
        F[7]=state[10]-state[(size_t)n+10]-state[2];
        const auto table=Verified::collocation(n);
        std::vector<I> denominators;denominators.reserve((size_t)n);
        for (int j=1;j<=n;++j) denominators.push_back(I(1)+(-I(2*j)*state[1]).exp());
        for (int m=0;m<=n;++m) {
            const I ordinate=state[(size_t)m+10];
            I psi,u,v;
            for (int j=1;j<=n;++j) {
                const I ep=(I(j)*ordinate).exp(),em=(-I(j)*ordinate-I(2*j)*state[1]).exp();
                const I S=finite ? (ep-em)/denominators[(size_t)j-1] : ep;
                const I C=finite ? (ep+em)/denominators[(size_t)j-1] : ep;
                const auto& trig=table[(size_t)((m*j)%(2*n))];
                const I coefficient=state[(size_t)n+10+j];
                psi=psi+coefficient*S*trig.second;
                u=u+I(j)*coefficient*C*trig.second;
                v=v+I(j)*coefficient*S*trig.first;
            }
            F[(size_t)8+m]=psi-state[8]-state[7]*ordinate;
            F[(size_t)n+9+m]=((u-state[7]).square()+v.square())/I(2)+ordinate-state[9];
        }
        const I wave_number=state[1],speed_scale=(I(1)/wave_number).sqrt();
        const I energy_scale=I(1)/wave_number.square(),density(1),gravity(1);
        const I pulse=state[8]+state[1]*state[5],flux_scaled=state[7]*state[1]-state[8];
        const I kinetic=(state[4]*pulse-state[5]*flux_scaled)/I(2);
        I potential;
        for (int j=1;j<n;++j) potential=potential+I(Y[(size_t)j]).square()/I(4);
        potential=potential+I(Y[(size_t)n]).square()/I(16);
        const I beta=I(2)*state[9]-state[4].square()+I(2)*state[4]*state[5];
        const I stress=I(4)*kinetic-I(3)*potential+beta*state[1]-I(2)*state[5]*pulse;
        const I power=state[4]*(I(3)*kinetic-I(2)*potential)+beta*(pulse+state[4]*state[1])/I(2)
                     -I(2)*state[4]*state[5]*pulse;
        I bed_sum;
        for (int j=1;j<=n;++j) {
            const I inverse_cosh=I(2)*(-I(j)*state[1]).exp()/denominators[(size_t)j-1];
            bed_sum=bed_sum+(I(j)*state[(size_t)n+10+j]*inverse_cosh).square()/I(2);
        }
        const std::map<std::string,I> outputs={
            {"L",I(2)*Verified::pi()/wave_number},{"c",state[4]*speed_scale},
            {"Impulse",density*speed_scale/wave_number*pulse},{"KineticEnergy",energy_scale*kinetic},
            {"PotentialEnergy",energy_scale*potential},{"EnergyDensity",energy_scale*(kinetic+potential)},
            {"Sxx",energy_scale*stress},{"Power",energy_scale*speed_scale*power},
            {"MeanSquareBedVelocity",gravity/wave_number*bed_sum}
        };
        const std::map<std::string,Real> reported={
            {"L",L},{"c",c},{"Impulse",::pulse/std::pow(kd,1.5)},
            {"KineticEnergy",ke/(kd*kd)},{"PotentialEnergy",pe/(kd*kd)},
            {"EnergyDensity",(ke+pe)/(kd*kd)},{"Sxx",sxx/(kd*kd)},
            {"Power",f/std::pow(kd,2.5)},{"MeanSquareBedVelocity",ub2/kd}
        };
        auto& certificate=finite_system_certificate;
        certificate.state.assign(1,0.0);
        certificate.state.insert(certificate.state.end(),z+1,z+num+1);
        certificate.surface_coefficients.assign(Y,Y+n+1);
        certificate.scope=finite ? "Collocation residuals and nine depth-scaled outputs"
                                 : "Deep-water collocation residuals; finite-depth integral outputs excluded";
        certificate.order=n;certificate.equation_count=(int)F.size();certificate.residual_intervals=F;
        // The lower adjacent constants are conservative decimal tolerance bounds.
        certificate.residual_tolerance=I::down(1e-10);
        for (const auto& residual : F) certificate.residual_upper=std::max(certificate.residual_upper,residual.magnitude());
        bool valid=certificate.residual_upper<=certificate.residual_tolerance;
        if (finite) for (const auto& entry : outputs) {
            const auto& interval=entry.second;
            const Real minimum=interval.lo<=0.0 && interval.hi>=0.0 ? 0.0 : std::min(std::abs(interval.lo),std::abs(interval.hi));
            Verified::OutputBound check;
            check.enclosure=interval;check.reported=reported.at(entry.first);
            check.error_upper=(I(check.reported)-interval).magnitude();
            check.budget_lower=I::down(I::down(1e-7)+I::down(I::down(1e-8)*minimum));
            valid=valid && check.error_upper<=check.budget_lower;
            certificate.outputs.emplace(entry.first,check);
        }
        certificate.passed=valid;
        return valid;
    }

    static std::pair<std::vector<Real>, std::vector<Real>> gauss_legendre(int count) {
        std::vector<Real> nodes((size_t)count), weights((size_t)count);
        for (int i = 0; i < (count+1)/2; ++i) {
            Real root = std::cos(Phys::PI*(i+0.75)/(count+0.5)), derivative = 0.0;
            bool ready = false;
            for (int iteration = 0; iteration < 100; ++iteration) {
                Real p0 = 1.0, p1 = root;
                for (int j = 2; j <= count; ++j) {
                    const Real next = ((2*j-1)*root*p1-(j-1)*p0)/j;
                    p0 = p1; p1 = next;
                }
                derivative = count*(root*p1-p0)/(root*root-1.0);
                const Real delta = p1/derivative;
                if (std::abs(delta) <= 4*std::numeric_limits<Real>::epsilon()) { ready = true; break; }
                root -= delta;
            }
            if (!ready) throw std::runtime_error("Gauss-Legendre node iteration did not converge.");
            const Real weight = 2.0/((1.0-root*root)*derivative*derivative);
            nodes[(size_t)i] = -root; nodes[(size_t)(count-1-i)] = root;
            weights[(size_t)i] = weights[(size_t)(count-1-i)] = weight;
        }
        return {nodes, weights};
    }


static std::map<std::string, Real> DefiningIntegrals(int phase_count, int vertical_count) {
        const auto quadrature = gauss_legendre(vertical_count);
        Real impulse = 0.0, kinetic = 0.0, potential = 0.0, stress = 0.0, power = 0.0;
        std::vector<Real> cosine((size_t)n+1), sine((size_t)n+1);
        for (int i = 0; i < phase_count; ++i) {
            const Real phase = 2.0*Phys::PI*i/phase_count;
            const Real eta = Surface(phase), thickness = (z[1]+eta)/2.0;
            if (!(thickness > 0.0) || !std::isfinite(thickness))
                throw std::runtime_error("The interpolated surface intersects the bed or is non-finite.");
            potential += 0.5*eta*eta/phase_count;
            for (int j = 1; j <= n; ++j) {
                cosine[(size_t)j] = std::cos(j*phase);
                sine[(size_t)j] = std::sin(j*phase);
            }
            for (int m = 0; m < vertical_count; ++m) {
                const Real y = -z[1]+thickness*(quadrature.first[(size_t)m]+1.0);
                Real u = z[5], v = 0.0;
                for (int j = 1; j <= n; ++j) {
                    const Real a = j*y, b = j*z[1];
                    const Real denominator = 1.0+std::exp(-2*b);
                    const Real positive = std::exp(a), reflected = std::exp(-a-2*b);
                    const Real amplitude = j*B[(size_t)j];
                    u += (positive+reflected)/denominator*amplitude*cosine[(size_t)j];
                    v += (positive-reflected)/denominator*amplitude*sine[(size_t)j];
                }
                const Real pressure = z[9]-y-0.5*((u-z[4])*(u-z[4])+v*v);
                const Real weight = quadrature.second[(size_t)m]*thickness/phase_count;
                impulse += weight*u;
                kinetic += weight*0.5*(u*u+v*v);
                stress += weight*(pressure+u*u);
                power += weight*(pressure+y+0.5*(u*u+v*v))*u;
            }
        }
        return {{"Impulse", impulse}, {"KineticEnergy", kinetic}, {"PotentialEnergy", potential},
                {"Sxx", stress-0.5*z[1]*z[1]}, {"Power", power}};
    }

static double ComparisonRatio(double value, double reference)
{
  if (!std::isfinite(value) || !std::isfinite(reference))
    return std::numeric_limits<double>::infinity();
  return std::abs(value-reference)/(1e-7+1e-8*std::abs(reference));
}

static bool VerifyAccuracy()
{
  spectral_defect = SpectralDefect();
  integral_error_ratio = spectral_error_ratio = 0.0;
  const bool interval_pass = VerifyFiniteSystem();
  if (!std::isfinite(spectral_defect) || spectral_defect > 1e-6 || !interval_pass) return false;
  if (std::strcmp(Depth,"Finite") != 0) return true;
  const auto reported = SelectedOutputs();
  const double bed_absolute = ub2+z[5]*z[5];
  const double spectral_stress = 4.0*ke-3.0*pe+bed_absolute*kd-2.0*z[5]*pulse;
  const double spectral_power = z[4]*(3.0*ke-2.0*pe)
      +0.5*bed_absolute*(pulse+z[4]*kd)-2.0*z[4]*z[5]*pulse;
  spectral_error_ratio = std::max(ComparisonRatio(sxx,spectral_stress),ComparisonRatio(f,spectral_power));
  const auto coarse = DefiningIntegrals(std::max(256,8*n),64);
  const auto fine = DefiningIntegrals(std::max(512,16*n),96);
  for (const auto& entry : fine) {
    integral_error_ratio = std::max({integral_error_ratio,
        ComparisonRatio(coarse.at(entry.first),entry.second),
        ComparisonRatio(reported.at(entry.first),entry.second)});
  }
  return spectral_error_ratio<=1.0 && integral_error_ratio<=1.0;
}

/**
 * Write all output files for the converged solution:
 *   - integral quantities (solution.res),
 *   - Fourier coefficients,
 *   - free-surface coordinates (surface.res),
 *   - flowfield profiles (flowfield.res).
 *
 * Integral quantities follow the notation in Fenton (1988) and related summaries
 * in Fenton (1990, 1999). The program reports them in two non-dimensional systems
 * for finite depth:
 *   (i)  g & wavenumber k (native solver scaling),
 *   (ii) g & mean depth d (converted scaling).
 */
static void Output(void)
{
  int i, I;
  double X = 0.0, eta = 0.0, y = 0.0;
  kd = z[1];

  std::fprintf(monitor, "\n\n# Solution summary:\n\n");
  Is_finite Title_block(monitor);
  Is_deep Title_block_deep(monitor);

  Is_finite Title_block(Solution);
  Is_deep Title_block_deep(Solution);


  std::fprintf(monitor, "\n# Requested Fourier order: %d; active order: %d", requested_order, n);
  std::fprintf(monitor, "\n# Maximum absolute residual: %.16e", final_residual);
  std::fprintf(monitor, "\n# Equilibrated Jacobian condition: %.16e", jacobian_condition);
  std::fprintf(monitor, "\n# Between-node boundary defect: %.16e", spectral_defect);
  std::fprintf(monitor, "\n# Interval residual upper bound: %.16e", finite_system_certificate.residual_upper);
  std::fprintf(monitor, "\n# Verification scope: %s", finite_system_certificate.scope.c_str());
  Is_finite std::fprintf(monitor, "\n# Output budget in each nondimensional scaling: 1e-7 + 1e-8 * abs(reference)");
  Is_finite std::fprintf(monitor, "\n# Spectral/integral comparison-to-budget ratios: %.8e / %.8e",
                        spectral_error_ratio, integral_error_ratio);
  Is_finite std::fprintf(Solution, "\n# Stokes-Ursell number %7.3f", 0.5 * z[2] / std::pow(z[1], 3));
  std::fprintf(Solution, "\n\n# Integral quantities - notation from Fenton (1988)");
  Is_finite {
    std::fprintf(Solution, "\n# (1) Quantity, (2) symbol, solution non-dimensionalised by (3) g & wavenumber, and (4) g & mean depth\n");
    Results("# Water depth                        (d)", z[1], 1.);
  }
  Is_deep {
    std::fprintf(Solution, "\n# (1) Quantity, (2) symbol, solution non-dimensionalised by (3) g & wavenumber\n");
  }
  Results("# Wave length                   (lambda)", 2 * pi, L);
  Results("# Wave height                        (H)", z[2], H);
  Results("# Wave period                      (tau)", z[3], T);
  Results("# Wave speed                         (c)", z[4], c);
  Results("# Eulerian current                 (u1_)", z[5], ce);
  Results("# Stokes current                   (u2_)", z[6], cs);
  Results("# Mean fluid speed in frame of wave (U_)", z[7], ubar);
  Results("# Volume flux due to waves           (q)", z[8], z[8] / std::pow(kd, 1.5));
  Results("# Bernoulli constant                 (r)", z[9], z[9] / kd);

  Is_finite {
    Results("# Volume flux                        (Q)", Q * std::pow(kd, 1.5), Q);
    Results("# Bernoulli constant                 (R)", R * kd, R);
    Results("# Momentum flux                      (S)", s, s / kd / kd);
    Results("# Impulse                            (I)", pulse, pulse / std::pow(kd, 1.5));
    Results("# Kinetic energy                     (T)", ke, ke / kd / kd);
    Results("# Potential energy                   (V)", pe, pe / kd / kd);
    Results("# Mean square of bed velocity     (ub2_)", ub2, ub2 / kd);
    Results("# Radiation stress                 (Sxx)", sxx, sxx / kd / kd);
    Results("# Wave power                         (F)", f, f / std::pow(kd, 2.5));
  }

  std::fprintf(Solution, "\n\n# Dimensionless coefficients in Fourier series");
  std::fprintf(Solution, "\n# Potential/Streamfn\tSurface elevations");
  std::fprintf(Solution, "\n# j, B[j], & E[j], j=1..N\n");
  for (i = 1; i <= n; i++) {
    std::fprintf(Solution, "\n%2d\t%15.7e\t%15.7e", i, B[i], Y[i]);
  }
  std::fprintf(Solution, "\n\n");

  WriteKinematicOutputs(monitor);

  // Surface coordinates
  std::fprintf(Elevation, "# %s\n", Title);
  std::fprintf(Elevation, "%s\n", Method);
  std::fprintf(Elevation, "\n# Surface of wave - trough-crest-trough,");
  std::fprintf(Elevation, " note quadratic point spacing clustered around crest");
  Is_finite {
    std::fprintf(Elevation, "\n# Non-dimensionalised with respect to depth");
    std::fprintf(Elevation, "\n# X/d, eta/d, & check of surface pressure\n");
    std::fprintf(Elevation, "\n0.\t0.\t0. # Dummy point to scale plot\n");
  }
  Is_deep {
    std::fprintf(Elevation, "\n# Non-dimensionalised with respect to wavenumber");
    std::fprintf(Elevation, "\n# kX, k eta, & check of surface pressure\n");
  }

  for (i = -Surface_points / 2; i <= Surface_points / 2; i++) {
    X = 4 * pi * (i * std::fabs((double)i) / Surface_points / Surface_points);
    eta = Surface(X);
    Point(X, eta);
    Is_finite std::fprintf(Elevation, "\n%8.4lf\t%7.4f\t%7.0e", X / kd, 1 + eta / kd, Pressure);
    Is_deep std::fprintf(Elevation, "\n%8.4lf\t%7.4f\t%7.0e", X, eta, Pressure);
  }
  std::fprintf(Elevation, "\n\n");

  // Velocity and acceleration profiles
  std::fprintf(Flowfield, "# %s\n", Title);
  std::fprintf(Flowfield, "%s\n", Method);
  std::fprintf(Flowfield, "\n# Velocity and acceleration profiles and Bernoulli checks\n");
  std::fprintf(Flowfield, "\n# All quantities are dimensionless with respect to g and/or ");
  Is_deep std::fprintf(Flowfield, "k\n");
  Is_finite std::fprintf(Flowfield, "d\n");
  std::fprintf(Flowfield, "\n#*******************************************************************************");
  Is_finite {
    std::fprintf(Flowfield, "\n# y        u       v    dphi/dt   du/dt   dv/dt  du/dx   du/dy Bernoulli check  ");
    std::fprintf(Flowfield, "\n# -     -------------   -------  ------   -----  ------------- ---------------  ");
    std::fprintf(Flowfield, "\n# d        sqrt(gd)       gd        g       g       sqrt(g/d)        gd         ");
  }
  Is_deep {
    std::fprintf(Flowfield, "\n# ky       u       v    dphi/dt   du/dt   dv/dt  du/dx   du/dy Bernoulli check  ");
    std::fprintf(Flowfield, "\n#       -------------   -------  ------   -----  ------------- ---------------  ");
    std::fprintf(Flowfield, "\n#         sqrt(g/k)       g/k       g       g       sqrt(gk)        g/k         ");
  }
  std::fprintf(Flowfield, "\n#*******************************************************************************");

  std::fprintf(Flowfield, "\n# Note that increasing X/d and 'Phase' here describes half of a wave for");
  std::fprintf(Flowfield, "\n# X/d >= 0. In a physical problem, where we might have a constant x, because");
  std::fprintf(Flowfield, "\n# the phase X = x - c * t, then as time t increases, X becomes increasingly");
  std::fprintf(Flowfield, "\n# negative and not positive as passing down the page here implies.");

  for (I = 0; I <= Nprofiles; ++I) {
    X = pi * I / (Nprofiles);
    eta = Surface(X);
    Is_finite std::fprintf(Flowfield, "\n\n# X/d = %8.7f, Phase = %6.1f\n", X / kd, X * 180. / pi);
    Is_deep std::fprintf(Flowfield, "\n\n# Phase = %6.1f\n", X * 180. / pi);

    for (i = 0; i <= Points; ++i) {
      Is_finite {
        y = (i) * (1. + eta / kd) / (Points);
        Point(X, kd * (y - 1.));
      }
      Is_deep {
        y = -pi + (double)i / Points * (eta + pi);
        Point(X, y);
      }
      std::fprintf(Flowfield, "% 10.6f % 10.6f % 10.6f % 10.6f % 10.6f % 10.6f % 10.6f % 10.6f % 10.6f",
                   y, u, v, dphidt, ut, vt, ux, uy, Bernoulli_check);
    }
  }
  std::fprintf(Flowfield, "\n\n");

  // Diagnostic output for Um (as in original)
  double psi0, psi1, Um_, dd;

  Point(0., Surface(0.));
  psi1 = psi;
  Point(0., -kd);
  psi0 = psi;

  dd = 1. + Surface(0.) / kd;
  Um_ = (psi1 - psi0) / dd;
  std::printf("\nUm %g", Um_);
}


// -----------------------------------------------------------------------------
// main (Fourier.cpp)
// -----------------------------------------------------------------------------

int main(int argc, char** argv)
{
  /*
    =============================================================================
    PROGRAM ENTRY POINT (C++20)
    =============================================================================

    Execution modes
    --------------
    1) File mode (preferred, matches original Fenton-style driver):
       If a file named **data.dat** exists in the current working directory,
       this program reads one case from that file (until the first "FINISH" title).

       This preserves the legacy input semantics used in the original multi-file
       program and in Fenton-style example drivers (see Fenton, 1988; Fenton, 1999,
       §3.1.2 “Fourier approximation methods”).

    2) Command-line (CLI) mode:
       If **data.dat** does not exist, the program expects the problem definition
       to be provided on the command line, in a format aligned with the *data.dat*
       ordering (see the attached sample data.dat):

           fourier  H_over_d  Case  Value  CurrentCriterion  [CurrentValue]  N  Steps

       where (Instructions.pdf Table 6-1):
         - H_over_d         = relative wave height input.
                               * Finite depth:  H_over_d = H/d  (H_over_d > 0)
                               * Deep water  :  H_over_d < 0 is interpreted as -H/L
         - Case             = "Wavelength" or "Period"
         - Value            = λ/d (if Case="Wavelength") OR τ*sqrt(g/d) (if Case="Period")
         - CurrentCriterion = 1 (Eulerian mean current ū1) or 2 (Stokes/mass-transport current ū2)
         - CurrentValue     = OPTIONAL, nondimensional current magnitude:
                               CurrentValue = current / sqrt(g d)
                             If omitted, defaults to 0.0.
                             IMPORTANT: if provided, it is the argument immediately
                             after CurrentCriterion, matching the data.dat layout.
         - N                = number of Fourier components (spectral order)
         - Steps            = number of height steps used to reach the final H/d (robust stepping)

    Notes on “current”
    ------------------
    This implementation follows Fenton’s distinction between:
      - Eulerian mean current ū1 (Fenton 1999 Eq. 3.13, then Eq. 3.15 for τ-known cases)
      - Stokes/mass-transport current ū2 (Fenton 1999 Eq. 3.14, then Eq. 3.16)

    In the (finite-depth, τ-known) problem family, the relevant nondimensional input
    is ū1 / sqrt(g d) or ū2 / sqrt(g d), consistent with Fenton (1999), Eq. (3.15–3.16).
  */

  
auto print_usage = [&](std::FILE* out) {
  const char* exe = (argc > 0 && argv && argv[0]) ? argv[0] : "fourier";
  std::fprintf(out,
    "Usage:\n"
    "  %s\n"
    "      Reads data.dat if present (searched in current directory, then next to the executable).\n"
    "  %s H_over_d Case Value CurrentCriterion [CurrentValue] N Steps\n\n"
    "CLI arguments (used only if data.dat is not found):\n"
    "    H_over_d         = relative wave height.\n"
    "                      Finite depth:  H_over_d = H/d  (must be > 0)\n"
    "                      Deep water  :  H_over_d < 0 means H_over_d = -H/L (Instructions.pdf Table 6-1)\n"
    "    Case             = \"Wavelength\" or \"Period\"  (Instructions.pdf Table 6-1)\n"
    "    Value            = if Case=Wavelength: L/d\n"
    "                      if Case=Period    : T*sqrt(g/d)  (dimensionless period, Instructions.pdf §6.1.3)\n"
    "    CurrentCriterion = 1 for Euler (ū1), or 2 for Stokes/mass-transport (ū2) (Instructions.pdf §6.1.4)\n"
    "    CurrentValue     = OPTIONAL, dimensionless current magnitude for the chosen criterion\n"
    "                        (defaults to 0.0 if omitted)\n"
    "                        CurrentValue = current / sqrt(g d)\n"
    "                        IMPORTANT: must appear immediately after CurrentCriterion\n"
    "    N                = requested Fourier order (1..512); reduced if ill-conditioned\n"
    "    Steps            = initial height subdivision (1..4096); continuation is adaptive\n\n"
    "Output files (lowercase): solution.res, surface.res, flowfield.res\n",
    exe, exe
  );
};

  auto parse_int = [](std::string_view s, int& out) -> bool {
    // Core Guidelines: prefer from_chars for locale-independent fast parsing.
    const char* b = s.data();
    const char* e = s.data() + s.size();
    auto [p, ec] = std::from_chars(b, e, out);
    return (ec == std::errc{} && p == e);
  };

  auto parse_double = [](std::string_view s, double& out) -> bool {
    // from_chars for floating is C++20 but not uniformly implemented on older libstdc++.
    // For maximum cross-platform compatibility, fall back to strtod while validating.
    char* end = nullptr;
    const std::string tmp(s);
    out = std::strtod(tmp.c_str(), &end);
    return end && end != tmp.c_str() && *end == '\0' && std::isfinite(out);
  };

  monitor = stdout;
  std::strcpy(Theory, "Fourier");         // This program implements the Fourier/stream-function method.
  std::strcpy(Diagname, "catalogue.res"); // retained legacy variable (not used for file output here).

  // ---------------------------------------------------------------------------
  // Determine execution mode (file vs CLI)
  // ---------------------------------------------------------------------------
  // Locate data.dat.
  //
  // Per Fenton's user manual (Instructions.pdf, §6 "Input data"): the controlling files
  // should be in the same directory as FOURIER.EXE. On Windows it is common to
  // start a program by double-clicking the executable; in that case the *current
  // working directory* may not be the executable directory. To make the program
  // behave predictably across Windows and Linux, we search for data.dat in:
  //   (1) the current working directory, then
  //   (2) the executable directory (best-effort).
  //
  // The directory where data.dat is found is also used as the output directory,
  // so the produced *.res files appear next to the controlling input file.
  auto find_data_file = [&]() -> fs::path {
    const fs::path cwd = fs::current_path();
    const fs::path in_cwd = cwd / "data.dat";
    if (fs::exists(in_cwd)) return in_cwd;

    
    // Best-effort executable directory discovery (cross-platform):
    //   - Windows: GetModuleFileNameW
    //   - Linux  : /proc/self/exe
    //   - Fallback: argv[0] resolved against current directory
    fs::path exe_dir = cwd;

#if defined(_WIN32)
    {
      // Core Guidelines: avoid raw owning pointers; use std::wstring as an owning buffer.
      std::wstring wbuf(32768, L'\0');
      const DWORD len = ::GetModuleFileNameW(nullptr, wbuf.data(), static_cast<DWORD>(wbuf.size()));
      if (len > 0 && len < wbuf.size()) {
        wbuf.resize(len);
        const fs::path exe_path = fs::path(wbuf);
        if (!exe_path.empty()) exe_dir = exe_path.parent_path();
      }
    }
#elif defined(__linux__)
    {
      // readlink does not null-terminate; we handle that explicitly.
      char buf[PATH_MAX] = {0};
      const ssize_t len = ::readlink("/proc/self/exe", buf, sizeof(buf) - 1);
      if (len > 0) {
        buf[len] = '\0';
        const fs::path exe_path = fs::path(buf);
        if (!exe_path.empty()) exe_dir = exe_path.parent_path();
      }
    }
#endif

    // Fallback for platforms without a reliable self-exe path (or if the above fails).
    if (exe_dir.empty() || exe_dir == cwd) {
      if (argc > 0 && argv && argv[0] && std::strlen(argv[0]) > 0) {
        fs::path exe_path = fs::path(argv[0]);
        if (exe_path.is_relative()) exe_path = cwd / exe_path;

        std::error_code ec{};
        fs::path canon = fs::weakly_canonical(exe_path, ec);
        exe_dir = (!ec ? canon.parent_path() : exe_path.parent_path());
        if (exe_dir.empty()) exe_dir = cwd;
      }
    }

    const fs::path in_exe = exe_dir / "data.dat";
    if (fs::exists(in_exe)) return in_exe;

    return {};
  };

  const fs::path data_path = find_data_file();
  const bool have_data_file = !data_path.empty();

  // Default I/O directory:
  //   - file mode: directory containing data.dat
  //   - CLI mode : current working directory
  fs::path io_dir = fs::current_path();
  if (have_data_file) io_dir = data_path.parent_path();

std::fprintf(stdout, "Working directory: %s\n", fs::current_path().string().c_str());
if (have_data_file) {
  std::fprintf(stdout, "Using input file:  %s\n", fs::absolute(data_path).string().c_str());
  std::fprintf(stdout, "I/O directory:     %s\n", fs::absolute(io_dir).string().c_str());
} else {
  std::fprintf(stdout, "data.dat not found; entering CLI mode.\n");
  std::fprintf(stdout, "I/O directory:     %s\n", fs::absolute(io_dir).string().c_str());
}



  if (have_data_file) {
    // --------------------------
    // File mode: read data.dat
    // --------------------------
    FilePtr input{std::fopen(data_path.string().c_str(), "r")};
    if (!input) {
      std::perror(data_path.string().c_str());
      return 1;
    }
    Input1 = input.get();

    if (!Read_data()) {
      std::fprintf(stderr, "No valid case found in data.dat (incomplete input or FINISH).\n");
      return 1;
    }
  } else {
    // --------------------------
    // CLI mode: parse arguments
    // --------------------------

// The CLI mirrors the *data.dat* column order from Instructions.pdf (Table 6-1),
// excluding the title line and the terminating "FINISH" sentinel:
//
//   H_over_d
//   Case ("Wavelength" or "Period")
//   Value (L/d or T*sqrt(g/d) respectively)
//   CurrentCriterion (1 Euler, 2 Stokes)
//   CurrentValue (optional; default 0.0)   <-- must appear immediately after criterion
//   N
//   Steps
//
// Therefore:
//   no CurrentValue:  fourier H Case Value Crit N Steps        (argc == 7)
//   with CurrentValue: fourier H Case Value Crit Cur N Steps   (argc == 8)
//
// We also support "-h/--help" explicitly.

if (argc == 2 && (std::strcmp(argv[1], "-h") == 0 || std::strcmp(argv[1], "--help") == 0)) {
  print_usage(stdout);
  return 0;
}

if (!(argc == 7 || argc == 8)) {
  std::fprintf(stderr, "data.dat not found (searched current directory and executable directory).\n\n");
  print_usage(stderr);
  return 2;
}

double H_in   = 0.0;
double value  = 0.0;
int    crit_id = 0;
double cur_nd  = 0.0;
int    N       = 0;
int    steps   = 0;

if (!parse_double(argv[1], H_in)) {
  std::fprintf(stderr, "Invalid H_over_d argument.\n\n");
  print_usage(stderr);
  return 2;
}

// Case is a string token (exactly "Wavelength" or "Period" in the supplied tools).
const std::string case_arg = argv[2];

if (!parse_double(argv[3], value) || !parse_int(argv[4], crit_id)) {
  std::fprintf(stderr, "Invalid Value or CurrentCriterion argument.\n\n");
  print_usage(stderr);
  return 2;
}

const bool have_cur = (argc == 8);
const int  idx_N    = have_cur ? 6 : 5;
const int  idx_S    = have_cur ? 7 : 6;

if (have_cur) {
  if (!parse_double(argv[5], cur_nd)) {
    std::fprintf(stderr, "Invalid CurrentValue argument.\n\n");
    print_usage(stderr);
    return 2;
  }
}

if (!parse_int(argv[idx_N], N) || !parse_int(argv[idx_S], steps)) {
  std::fprintf(stderr, "Invalid N or Steps argument.\n\n");
  print_usage(stderr);
  return 2;
}

// Range checks consistent with Instructions.pdf Table 6-1 semantics.
if (!(crit_id == 1 || crit_id == 2) || N <= 0 || steps <= 0) {
  std::fprintf(stderr, "Arguments out of range (CurrentCriterion must be 1 or 2; N and Steps must be > 0).\n\n");
  print_usage(stderr);
  return 2;
}

if (value <= 0.0) {
  std::fprintf(stderr, "Value must be > 0 (L/d or T*sqrt(g/d)).\n\n");
  print_usage(stderr);
  return 2;
}

// Configure solver state to match the data.dat semantics, preserving all legacy branches
// and formulas in Read_data().

    // Configure solver state to match the *data.dat* “Period” finite-depth case.

std::strncpy(Title, "COMMAND LINE CASE", sizeof(Title) - 1);
Title[sizeof(Title) - 1] = '\0';

// H_over_d sign convention (Instructions.pdf Table 6-1):
//   MaxH >= 0  -> finite depth, MaxH = H/d
//   MaxH <  0  -> deep water,  MaxH = -H/L (and the length lines are ignored by the legacy code)
MaxH = H_in;

if (MaxH >= 0.0) {
  std::strcpy(Depth, "Finite");

  // Case selection is preserved exactly as in the file driver.
  if (case_arg == "Wavelength" || case_arg == "wavelength") {
    std::strcpy(Case, "Wavelength");
    L = value;                 // L/d
    Height = MaxH / L;         // legacy internal height parameter for finite depth wavelength-known case
  } else if (case_arg == "Period" || case_arg == "period") {
    std::strcpy(Case, "Period");
    T = value;                 // T*sqrt(g/d)
    Height = MaxH / (T * T);   // legacy internal height parameter for finite depth period-known case
  } else {
    std::fprintf(stderr, "Case must be \\\"Wavelength\\\" or \\\"Period\\\" (got: %s).\n\n", case_arg.c_str());
    print_usage(stderr);
    return 2;
  }
} else {
  std::strcpy(Depth, "Deep");
  // In the deep-water branch of Read_data(), the Case and Value lines are read but skipped.
  // We keep the same interface here for correspondence with Instructions.pdf Table 6-1.
  if (case_arg == "Wavelength" || case_arg == "wavelength") std::strcpy(Case, "Wavelength");
  else if (case_arg == "Period" || case_arg == "period")    std::strcpy(Case, "Period");
  else {
    std::fprintf(stderr, "Case must be \\\"Wavelength\\\" or \\\"Period\\\" (got: %s).\n\n", case_arg.c_str());
    print_usage(stderr);
    return 2;
  }

  Height = -MaxH;              // legacy: Height = H/L for deep water (MaxH stores -H/L)
}

Current_criterion = crit_id;
Current           = have_cur ? cur_nd : 0.0;

if (Current_criterion == 1) std::strcpy(Currentname, Current1);
if (Current_criterion == 2) std::strcpy(Currentname, Current2);

n     = N;
nstep = steps;

    // Convergence controls (embedded constants; formerly read from Convergence.dat)
    number = kConvergenceMaxIterations;
    crit   = kConvergenceCriterion;

    // Output sampling controls (embedded constants; formerly read from Points.dat)
    Surface_points = kSurfacePointsDefault;
    Nprofiles      = kNprofilesDefault;
    Points         = kProfilePointsDefault;

    // Print a faithful echo of the “input block” to stdout, as in file mode.
    Input_Data_block(monitor);
  }

  if (!std::isfinite(MaxH) || MaxH == 0.0 || !std::isfinite(Height) || Height <= 0.0 ||
      !std::isfinite(Current) || n < 1 || n > 512 || nstep < 1 || nstep > 4096 ||
      (Current_criterion != 1 && Current_criterion != 2) ||
      (std::strcmp(Case, "Period") != 0 && std::strcmp(Case, "Wavelength") != 0)) {
    std::fprintf(stderr, "Invalid input: finite nonzero height, positive period/wavelength, "
                         "criterion 1 or 2, 1 <= N <= 512 and 1 <= Steps <= 4096 are required.\n");
    return 2;
  }

  // ---------------------------------------------------------------------------
  // Allocate solver work arrays (NR-style 1-based indexing, but RAII-managed)
  // ---------------------------------------------------------------------------
  num = 2 * n + 10;

#if defined(_OPENMP)
  // Respect an explicit OMP_NUM_THREADS setting.  Otherwise use a conservative
  // workload-dependent limit to avoid excessive OpenMP overhead on high-core-count
  // systems for the relatively small dense systems used by this solver.
  if (std::getenv("OMP_NUM_THREADS") == nullptr) {
    const int useful_threads = std::max(1, std::min(8, num / 24));
    omp_set_num_threads(std::min(useful_threads, omp_get_num_procs()));
  }
  omp_set_dynamic(0);
#endif


  // RAII owners for long-lived allocations. Raw pointers are retained for
  // compatibility with the original indexing and formula structure.
  struct DVectorDeleter {
    long nl{};
    long nh{};
    void operator()(double* p_) const noexcept {
      if (p_) free_dvector(p_, nl, nh);
    }
  };
  struct DMatrixDeleter {
    long nrl{}, nrh{}, ncl{}, nch{};
    void operator()(double** p_) const noexcept {
      if (p_) free_dmatrix(p_, nrl, nrh, ncl, nch);
    }
  };

  using DVecPtr = std::unique_ptr<double, DVectorDeleter>;
  using DMatPtr = std::unique_ptr<double*, DMatrixDeleter>;

  auto make_dvector_owner = [](long nl, long nh) -> DVecPtr {
    return DVecPtr(dvector(nl, nh), DVectorDeleter{nl, nh});
  };
  auto make_dmatrix_owner = [](long nrl, long nrh, long ncl, long nch) -> DMatPtr {
    return DMatPtr(dmatrix(nrl, nrh, ncl, nch), DMatrixDeleter{nrl, nrh, ncl, nch});
  };

  auto Y_owner     = make_dvector_owner(0, num);         Y = Y_owner.get();
  auto z_owner     = make_dvector_owner(1, num);         z = z_owner.get();
  auto rhs1_owner  = make_dvector_owner(1, num);         rhs1 = rhs1_owner.get();
  auto cosa_owner  = make_dvector_owner(0, 2 * n);       cosa = cosa_owner.get();
  auto sina_owner  = make_dvector_owner(0, 2 * n);       sina = sina_owner.get();
  auto sol_owner   = make_dmatrix_owner(0, num, 1, 2);   sol = sol_owner.get();
  auto B_owner     = make_dvector_owner(1, n);           B = B_owner.get();
  auto Tanh_owner  = make_dvector_owner(1, n);           Tanh = Tanh_owner.get();

  requested_order = n;
  const int minimum_order = std::min(n, 12);
  bool solved = false;
  for (;;) {
    num = 2 * n + 10;
    solved = ContinueToTarget();
    RefreshCondition();
    std::fprintf(monitor, "\nFourier order %d: residual %.3e, condition %.3e\n",
                 n, final_residual, jacobian_condition);
    if (jacobian_condition <= kConditionLimit || n <= minimum_order) break;
    n = std::max(minimum_order, static_cast<int>(0.8 * n));
    std::fprintf(monitor, "Retrying with Fourier order %d.\n", n);
  }
  if (!solved || !std::isfinite(jacobian_condition) || jacobian_condition > kConditionLimit) {
    std::fprintf(stderr, "Numerical failure: no verified solution at the requested height. "
                         "No result files were written.\n");
    return 3;
  }
  std::snprintf(Method, sizeof(Method), "# Solution by %d-term Fourier series", n);
  Y[0] = 0.0;
  for (int j = 1; j <= n; ++j) {
    B[j] = z[j + n + 10];
    double sum = 0.5 * (z[10] + z[n + 10] * (j % 2 ? -1.0 : 1.0));
    for (int m = 1; m < n; ++m) sum += z[10 + m] * cosa[(m * j) % (2 * n)];
    Y[j] = 2.0 * sum / n;
  }
  CalculateIntegralOutputs();
  try {
    if (!VerifyAccuracy()) {
      std::fprintf(stderr,"Numerical accuracy checks failed: surface residual %.3e, interval residual %.3e, "
                          "spectral comparison/budget %.3e, integral comparison/budget %.3e. "
                          "No result files were written for this run.\n",
                          spectral_defect,finite_system_certificate.residual_upper,
                          spectral_error_ratio,integral_error_ratio);
      return 4;
    }
    CalculateKinematicOutputs();
  } catch (const std::exception& error) {
    std::fprintf(stderr,"Numerical verification failed: %s. No result files were written for this run.\n",error.what());
    return 4;
  }

// To make it obvious to users *where* files are written (especially on Windows,
// where launching an executable from Explorer may yield a surprising working
// directory), we print absolute output paths here.
const fs::path out_solution  = fs::absolute(io_dir / "solution.res");
const fs::path out_surface   = fs::absolute(io_dir / "surface.res");
const fs::path out_flowfield = fs::absolute(io_dir / "flowfield.res");

std::fprintf(stdout, "\nWriting output files to:\n");
std::fprintf(stdout, "  %s\n", out_solution.string().c_str());
std::fprintf(stdout, "  %s\n", out_surface.string().c_str());
std::fprintf(stdout, "  %s\n", out_flowfield.string().c_str());
  FilePtr solution_fp{std::fopen((io_dir / "solution.res").string().c_str(), "w")};
  FilePtr elevation_fp{std::fopen((io_dir / "surface.res").string().c_str(), "w")};
  FilePtr flowfield_fp{std::fopen((io_dir / "flowfield.res").string().c_str(), "w")};

  if (!solution_fp || !elevation_fp || !flowfield_fp) {
    std::perror("output");
    return 1;
  }

  Solution  = solution_fp.get();
  Elevation = elevation_fp.get();
  Flowfield = flowfield_fp.get();

  Output();

  auto close_result = [](FilePtr& file, const char* name) {
    bool failed = std::ferror(file.get()) != 0;
    if (std::fclose(file.release()) != 0) failed = true;
    if (failed) std::fprintf(stderr, "Failed to finish writing %s.\n", name);
    return !failed;
  };
  const bool solution_closed = close_result(solution_fp, "solution.res");
  const bool elevation_closed = close_result(elevation_fp, "surface.res");
  const bool flowfield_closed = close_result(flowfield_fp, "flowfield.res");
  Solution = Elevation = Flowfield = nullptr;
  if (!solution_closed || !elevation_closed || !flowfield_closed) return 1;

  std::printf("\nFinished\n");
  std::fflush(stdout);
  std::fflush(stderr);
#if defined(_WIN32)
  // Release a console used exclusively by this process.
  DWORD console_processes[2] = {};
  const DWORD console_count = ::GetConsoleProcessList(console_processes, 2);
  if (console_count == 1 && console_processes[0] == ::GetCurrentProcessId())
    ::FreeConsole();
#endif
  return 0;
}
