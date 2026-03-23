/*
================================================================================
fourier.cpp  —  Single-file implementation of J.D. Fenton-style steady-wave
Fourier/Newton “Stream Function” solver (amalgamated / modernised build unit)
Original source code - Fourier.zip, version 23 July 2015.
See URL: https://johndfenton.com/Steady-waves/Fourier.html
================================================================================

OVERVIEW
--------
This program computes *steady, two-dimensional, periodic gravity waves* propagating
without change of form over either:

  • water of *finite* mean depth (mean depth d), or
  • *deep water* (infinite depth limit).

The numerical method is the classical *Fourier approximation / stream-function*
approach combined with *collocation* on the free surface and a *Newton iteration*
for the resulting nonlinear algebraic system.

The implementation you are reading is an **amalgamation** (single translation unit)
of the original multi-file C++ port that in turn follows the structure and notation
of J.D. Fenton’s published algorithms and example program listings, principally:

  • Fenton (1999) — Review of fully-nonlinear wave computation methods (Fourier,
                    BIE, polynomial approximation).
  • Fenton (1990) — Survey and selection guidance: Stokes vs cnoidal vs stream-
                    function, etc.
  • Fenton (1988) — Core stream-function/Fourier algorithm for steady, periodic
                    nonlinear waves.

The code computes both:
  (A) *global / integral* wave properties (c, T, Q, energies, radiation stress, …), and
  (B) *local* kinematics and pressures (u, v, ∂φ/∂t, accelerations, Bernoulli check)
      at arbitrary phase and elevation points, and outputs those to files.

KEY CAPABILITIES
----------------
1) Problem types (input “Case”):
   • Wavelength-specified cases (“Wavelength”): solve for a wave with given λ/d and H/d.
   • Period-specified cases (“Period”): solve for a wave with given T√(g/d) and H/d.
     This case requires careful treatment of the current definition (Fenton, 1999
     discussion in the Fourier/Newton method section).

2) Water depth regimes:
   • Finite depth: kd = k d is solved as part of the system; outputs are provided in
     two consistent non-dimensionalisations:
       (i) by g and wavenumber k, and
       (ii) by g and mean depth d.
   • Deep water: the kd → ∞ limit is handled using the exponential vertical basis
     consistent with Fenton (1999) Eq. (3.6).

3) Currents (two alternative criteria, as in Fenton’s formulations):
   • “Euler” current criterion: specify the Eulerian mean current ū₁ (often denoted Uc).
   • “Stokes” current criterion: specify the Stokes mean mass-transport current ū₂.
   The program solves the wave on a current consistent with the chosen criterion.

4) Robust convergence for steep waves:
   The program can solve a sequence of increasing heights (“height stepping”)
   and extrapolate unknowns between steps. This is a standard robustness strategy
   noted in Fenton (1988, 1999) for near-limiting waves and/or long waves.

INPUTS AND OUTPUTS
------------------
This build has **one optional input file**:

  • data.dat
      The main run-control file. The attached sample shows the intended format:

        Title line             (string)
        H/d                    (positive => finite depth; negative legacy flag 
                                => deep water)
        Case                   ("Wavelength" or "Period")
        λ/d  OR  T*sqrt(g/d)   (depending on Case)
        CurrentCriterion       (1 = Eulerian ū1, 2 = Stokes/mass-transport ū2)
        CurrentValue           (nondimensional: current / sqrt(g d))
        N                      (number of Fourier components)
        Steps                  (number of height steps)

If **data.dat is absent**, you may run the program in **CLI mode** (see USAGE section
and the runtime Usage() printout). In CLI mode, the program solves the *finite-depth,
Period-specified* family (Case="Period") because the inputs provided are H/d and
T*sqrt(g/d) (Fenton 1999, discussion around Eq. 3.12–3.16).

IMPORTANT CHANGE: Points.dat EMBEDDED
-------------------------------------
The original program used a Points.dat file to control the output sampling density
(surface point count, number of profiles, points per profile).
**Points.dat is no longer required**: its supplied contents are embedded as constants.

IMPORTANT CHANGE: Convergence.dat EMBEDDED
------------------------------------------
The original program read a small Convergence.dat file containing:
  - maximum Newton iterations per height step
  - convergence tolerance (criterion)

**Convergence.dat is no longer required**. The values that were previously read from 
that file are now *embedded in this source* as constants:

    maximum iterations: 50
    convergence criterion: 1e-9

If you need different convergence behaviour you must edit those constants and
recompile.

FILES WRITTEN
-------------
The program writes the following result files:

  • solution.res   — summary of integral quantities, derived constants, and Fourier coefficients
  • surface.res    — free-surface coordinates (trough–crest–trough) and pressure / Bernoulli checks
  • flowfield.res  — phase-by-phase vertical profiles of velocity, accelerations, and Bernoulli checks

Note: The legacy redundant tabular file “Solution-Flat.res” has been intentionally removed.

NUMERICAL METHOD SUMMARY (WITH KEY FENTON REFERENCES)
-----------------------------------------------------
Governing equations in the travelling-wave formulation (consistent notation):

  • Stationary physical coordinates: (x, y), with bed at y = 0 and mean depth d
    (mean free surface at y = d).

  • Mean-level vertical coordinate: y' := y - d
      mean free surface: y' = 0
      bed:              y' = -d

  • Travelling horizontal coordinate: ξ := x - c t

  • Non-dimensional solver coordinates (used throughout the nonlinear system):
      X := k ξ
      Y := k y'
      η(X) := k η_phys(ξ)
    where k := 2π/λ and kd := k d.
    One wavelength corresponds to X ∈ [0, 2π]. In finite depth, the bed is Y = -kd.

Stream-function formulation in the travelling frame:

  • For 2D incompressible flow, a stream function ψ(X,Y) exists such that the
    *non-dimensional* travelling-frame velocities (Û, V̂) are:
        Û =  ∂ψ/∂Y,    V̂ = -∂ψ/∂X        (Fenton, 1999, Eq. 3.1)
    with the usual k-based scaling: velocities ~ √(g/k), ψ ~ √(g/k³).

  • Laplace equation in the fluid:
        ∂²ψ/∂X² + ∂²ψ/∂Y² = 0            (Fenton, 1999, Eq. 3.1)

  • Boundary conditions (finite depth):
        ψ(X, -kd)     = 0                (Fenton, 1999, Eq. 3.2)
        ψ(X,  η(X))   = -q               (Fenton, 1999, Eq. 3.3)
        ½(ψ_X²+ψ_Y²) + η(X) = R          (Fenton, 1999, Eq. 3.4)

    Deep water is recovered by the kd → ∞ limit, together with decay as Y → -∞.

Fourier/stream-function representation:
  ψ is expanded in a truncated Fourier series that satisfies Laplace and the
  bottom/decay condition identically, leaving only the nonlinear free-surface
  conditions to determine the coefficients (Fenton, 1999, Eq. 3.5; deep-water
  limit Eq. 3.6).

Collocation and nonlinear solve:
  The free-surface conditions are enforced at N+1 collocation points over half a
  wave (crest-to-trough symmetry), producing a nonlinear system. This program
  solves that system by Newton’s method, with the Jacobian obtained by numerical
  differencing (Fenton, 1988; see also Fenton, 1999 §3.1.2).

MEAN SQUARE BED ORBITAL VELOCITY (ub2)
--------------------------------------
This program reports a quantity labelled:

    “Mean square of bed velocity (ub2)”

Here, ub2 is computed by numerical phase-averaging of the near-bed *orbital*
horizontal velocity, taken relative to the Eulerian mean current ū₁:

    ub2 = ⟨ ( u_bed(t) - ū₁ )² ⟩

For a steady travelling wave, temporal averaging over one period is equivalent to
spatial averaging over one wavelength; this implementation samples phase uniformly
in X = k ξ over [0, 2π).

Implementation details:
  • Sample phase X uniformly over [0, 2π).
  • At each phase, evaluate the stationary-frame bed velocity u_bed at y = 0.
  • Subtract ū₁ (Eulerian mean current) to obtain the orbital component.
  • Square and average over the samples.

Note:
  Some legacy integral expressions in the original program use an internal
  auxiliary quantity derived from Bernoulli constants. That internal quantity is
  preserved for the Sxx/F-related formulas; only the *reported* ub2 quantity is
  defined as the physically standard mean-square orbital velocity.

COMPILATION
-----------
This is a single C++20 source file.

  Windows (MSYS2/MinGW):
      g++ -std=c++20 fourier.cpp -o fourier.exe -O2 -static -static-libgcc -static-libstdc++

  Linux/macOS:
      g++ -std=c++20 -Wall -Wextra -pedantic fourier.cpp -O2 -o fourier -static-libstdc++ -static-libgcc

USAGE
-----
Place data.dat in the working directory (or run in CLI mode), then run:

    fourier.exe

The data.dat format supports multiple cases terminated by a title line "FINISH"; this build reads the first case only (single-run driver).

LICENSE / ATTRIBUTION NOTES
----------------------------
• The hydrodynamic method and most of the mathematical formulation follows the
  published work of J.D. Fenton and earlier stream-function/Fourier approaches cited
  within those works (Chappelear, Dean, Chaplin, Rienecker & Fenton, etc.).
• The linear algebra helper (SVD-based solve) is consistent with the “Numerical Recipes”
  style routines (dsvdcmp/dsvbksb) included in the supplied source set.

BIBLIOGRAPHY
-------------

  1.  Fenton, J. D. (1999). "Numerical methods for nonlinear waves."
      in P. L.-F. Liu, ed., Advances in Coastal and Ocean Engineering, Vol. 5,
      World Scientific, Singapore, pp. 241-324.
      [Relevance: Review of fully-nonlinear wave computation methods (Fourier,
      BIE, polynomial approximation).]
      URL: https://johndfenton.com/Papers/Fenton99Liu-Numerical-methods-for-nonlinear-waves.pdf

  2.  Fenton, J. D. (1999). "The cnoidal theory of water waves."
      in J. B. Herbich, ed., Developments in Offshore Engineering, Gulf
      Publishing, Houston, pp. 275-337.
      [Relevance: Cnoidal-wave theory (finite-depth, long waves) and numerical
      formulations.]
      URL: https://johndfenton.com/Papers/Fenton99Cnoidal-The-cnoidal-theory-of-water-waves.pdf

  3.  Fenton, J. D. & Kennedy, A. B. (1996). "Fast methods for computing the
      shoaling of nonlinear waves."
      in Proc. 25th Int. Conf. Coastal Engng, Vol. 1, Orlando, pp. 1130-1143.
      [Relevance: Nonlinear wave propagation/shoaling over varying bathymetry.]
      URL: https://johndfenton.com/Papers/Fenton96%2BKennedy-Fast-methods-for-computing-the-shoaling-of-nonlinear-waves.pdf

  4.  Fenton, J. D. (1995). "A numerical cnoidal theory for steady water waves."
      in Proc. 12th Australasian Coastal and Ocean Engng Conference, Melbourne,
      pp. 175-180.
      [Relevance: Cnoidal-wave theory (finite-depth, long waves) and numerical
      formulations.]
      URL: https://johndfenton.com/Papers/Fenton95-A-numerical-cnoidal-theory-for-steady-water-waves.pdf

  5.  Townsend, M. & Fenton, J. D. (1995). "Numerical comparisons of wave
      analysis methods."
      in Proc. 12th Australasian Coastal and Ocean Engng Conference, Melbourne,
      pp. 169-173.
      [Relevance: Wave analysis/inversion from pressure measurements;
      conditioning and method comparison.]
      URL: https://johndfenton.com/Papers/Townsend95%2BFenton-Numerical-comparisons-of-wave-analysis-methods.pdf

  6.  Kennedy, A. B. & Fenton, J. D. (1995). "Simulation of the propagation of
      surface gravity waves using local polynomial approximation."
      in Proc. 12th Australasian Coastal and Ocean Engng Conference, Melbourne,
      pp. 287-292.
      [Relevance: Nonlinear wave propagation/shoaling over varying bathymetry.]
      URL: https://johndfenton.com/Papers/Kennedy95%2BFenton-Simulation-of-the-propagation-of-surface-gravity-waves-using-local-polynomial-approximation.pdf

  7.  Fenton, J. D. (1993). "Simulating wave shoaling with boundary integral
      equations."
      in Proc. 11th Australasian Conference on Coastal and Ocean Engng,
      Townsville, pp. 71-76.
      [Relevance: Boundary-integral formulation for nonlinear wave
      transformation (with singularity subtraction).]
      URL: https://johndfenton.com/Papers/Fenton93-Simulating-wave-shoaling-with-boundary-integral-equations.pdf

  8.  Fenton, J. D. (1990). "Nonlinear wave theories."
      in B. Le Méhauté & D. M. Hanes, eds, The Sea - Ocean Engineering Science,
      Part A, Vol. 9, Wiley, New York, pp. 3-25.
      [Relevance: Survey and selection guidance: Stokes vs cnoidal vs
      stream-function, etc.]
      URL: https://johndfenton.com/Papers/Fenton90b-Nonlinear-wave-theories.pdf

  9.  Drennan, W. M., Fenton, J. D. & Donelan, M. A. (1990). "Numerical
      simulation of nonlinear wave groups."
      in Proc. 11th Ann. Conf. Canadian Applied Math. Soc., Halifax.
      [Relevance: Relevant to nonlinear/free-surface wave modelling.]
      URL: https://johndfenton.com/Papers/Drennan90-Numerical-simulation-of-nonlinear-wave-groups.pdf

 10.  Fenton, J. D. & McKee, W. D. (1990). "On calculating the lengths of water
      waves."
      Coastal Engineering 14, 499-513.
      [Relevance: Wavelength determination for nonlinear waves in finite depth.]
      URL: https://johndfenton.com/Papers/Fenton90c%2BMcKee-On-calculating-the-lengths-of-water-waves.pdf

 11.  Fenton, J. D. (1988). "The numerical solution of steady water wave
      problems."
      Computers and Geosciences 14, 357-368.
      [Relevance: Core stream-function/Fourier algorithm for steady, periodic
      nonlinear waves.]
      URL: https://johndfenton.com/Papers/Fenton88-The-numerical-solution-of-steady-water-wave-problems.pdf

 12.  Fenton, J. D. (1986). "Polynomial approximation and water waves."
      in Proc. 20th Int. Conf. Coastal Engng, Vol. 1, Taipei, pp. 193-207.
      [Relevance: Relevant to nonlinear/free-surface wave modelling.]
      URL: https://johndfenton.com/Papers/Fenton86-Polynomial-approximation-and-water-waves.pdf

 13.  Fenton, J. D. (1985). "A fifth-order Stokes theory for steady waves."
      J. Waterway Port Coastal and Ocean Engng 111 , 216-234.
      [Relevance: Closed-form Stokes expansion (5th order) for steady waves in
      finite depth.]
      URL: https://johndfenton.com/Papers/Fenton85d-A-fifth-order-Stokes-theory-for-steady-waves.pdf

 14.  Fenton, J. D. (1983). "On the application of steady wave theories."
      in Proc. 6th Australasian Conf. Coastal and Ocean Engng, Christchurch,
      pp. 65-70.
      [Relevance: Guidance on applicability/limits of steady-wave theories in
      engineering.]
      URL: https://johndfenton.com/Papers/Fenton83-On-the-application-of-steady-wave-theories.pdf

 15.  Fenton, J. D. & Rienecker, M. M. (1982). "A Fourier method for solving
      nonlinear water wave problems."
      J. Fluid Mechanics 118, 411-443.
      [Relevance: Fourier-series method for solving fully nonlinear steady water
      waves.]
      URL: https://johndfenton.com/Papers/Fenton82c%2BRienecker-A-Fourier-method-for-solving-nonlinear-water-wave-problems.pdf

 16.  Schwartz, L. W. & Fenton, J. D. (1982). "Strongly-nonlinear waves."
      in M. Van Dyke, J. V. Wehausen & J. L. Lumley, eds, Ann. Rev. Fluid Mech.
      14 , 39-60.
      [Relevance: Fundamental properties/approximations for strongly nonlinear
      wave motion.]
      URL: https://johndfenton.com/Papers/Schwartz82-Strongly-nonlinear-waves.pdf

 17.  Rienecker, M. M. & Fenton, J. D. (1981). "A Fourier approximation method
      for steady water waves."
      J. Fluid Mechanics 104, 119-137.
      [Relevance: Fourier-series method for solving fully nonlinear steady water
      waves.]
      URL: https://johndfenton.com/Papers/Rienecker81%2BFenton-A-Fourier-approximation-method-for-steady-water-waves.pdf

 18.  Fenton, J. D. & Rienecker, M. M. (1980). "Accurate numerical solutions for
      nonlinear waves."
      in Proc. 17th Int. Conf. Coastal Engng, Sydney, pp. 50-69.
      [Relevance: Benchmark accurate numerical solutions for nonlinear wave
      profiles/kinematics.]
      URL: https://johndfenton.com/Papers/Fenton80%2BRienecker-Accurate-numerical-solutions-for-nonlinear-waves.pdf

 19.  Fenton, J. D. (1979). "A high-order cnoidal wave theory."
      J. Fluid Mechanics 94, 129-161.
      [Relevance: Cnoidal-wave theory (finite-depth, long waves) and numerical
      formulations.]
      URL: https://johndfenton.com/Papers/Fenton79-A-high-order-cnoidal-wave-theory.pdf

 20.  Fenton, J. D. & Mills, D. A. (1976). "Shoaling waves: numerical solution of
      exact equations."
      in D. G. Provis & R. Radok, eds, Proc. IUTAM Symposium on Waves on
      Water of Variable Depth, Canberra, Springer-Verlag, pp. 93-100.
      [Relevance: Nonlinear wave propagation/shoaling over varying bathymetry.]
      URL: https://johndfenton.com/Papers/Fenton76%2BMills-Shoaling-waves-Numerical-solution-of-exact-equations.pdf

 21.  Fenton, J. D. (1972). "A ninth-order solution for the solitary wave."
      J. Fluid Mechanics 53, 257-271.
      [Relevance: High-order analytic/series solution for solitary waves.]
      URL: https://johndfenton.com/Papers/Fenton72-A-ninth-order-solution-for-the-solitary-wave.pdf

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


namespace fs = std::filesystem;


/*
================================================================================
TECHNICAL ROADMAP FOR THIS TRANSLATION UNIT
================================================================================

This file is intentionally "monolithic": it embeds what historically were several
source modules (Allocation, Solve, Subroutines, Inout, main) into a single build
unit to simplify compilation and reproducibility.

The numerical heart of the program is unchanged from the legacy Fourier/stream-
function driver: *only comments, I/O robustness, and lifetime management were
modernised*.

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

This solver uses Numerical-Recipes style 1-based arrays. The unknown vector z[]
therefore starts at index 1.

Scalar unknowns (indices 1..9) represent global wave parameters.
Vector unknowns represent:
  - collocated free-surface elevations (kη at N+1 points)
  - Fourier coefficients for the streamfunction/potential series

A convenient mental model is:

  z = [ global scalars | surface elevations | Fourier coefficients ]

More precisely (see Eqns() for enforcement):

  z[1]   : kd        finite depth kd, in deep water a dummy normalisation (=1)
  z[2]   : kH        wave height in k-based scaling
  z[3]   : kc / ω    period parameter used by the legacy code (ties to τ and dispersion)
  z[4]   : c √(k/g)  wave speed parameter in k-based scaling
  z[5]   : u1 √(k/g) Eulerian mean current in k-based scaling
  z[6]   : u2 √(k/g) Stokes/mass-transport current in k-based scaling
  z[7]   : U √(k/g)  mean fluid speed in the moving frame (k-based)
  z[8]   : (Ū*kd-q)  auxiliary streamline-constant combination (k-based)
  z[9]   : r         Bernoulli constant in k-based scaling

  z[10 + m],  m=0..N  : collocated kη_m values at free-surface collocation points
                        over half a wave (crest to trough). The collocation phases
                        are X_m = m π / N (implied by the cosine/sine tables).

  z[n + 10 + j], j=1..N : B_j coefficients in the Fourier/stream-function series.

Residual vector rhs[] is organised in the same overall order:
  rhs[1..8]            : global closure equations (height, period/length, current, mean level)
  rhs[9..9+N]          : free-surface streamfunction (kinematic) BC at each collocation point
  rhs[n+10..n+10+N]    : free-surface Bernoulli (dynamic) BC at each collocation point

The Newton solver drives rhs → 0.

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

The Newton step solves a dense system with a numerically-differenced Jacobian.
Near limiting steep waves, the Jacobian can be ill-conditioned; therefore the
solver uses an SVD-based least-squares solve (Numerical Recipes dsvdcmp/dsvbksb).

The singular-value cutoff wmin = wmax * 1e-12 is a pragmatic stabilisation:
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
FENTON STEADY NONLINEAR WAVE THEORY (STREAM-FUNCTION / FOURIER-APPROXIMATION)
================================================================================

This section documents the hydrodynamic model and the mapping from the published
Fenton formulation to the variables and nonlinear algebraic system solved here.

The goal is to make the "wave theory" content explicit enough that:
  • the governing equations can be audited line-by-line against the implementation,
  • sign conventions and reference frames are unambiguous,
  • future maintenance does not accidentally corrupt the mathematical model.

------------------------------------------------------------------------------
A) PHYSICAL MODEL AND SCOPE
------------------------------------------------------------------------------

The wave is modelled as a *two-dimensional*, *inviscid*, *incompressible*,
*irrotational* gravity wave of permanent form, propagating over a horizontal bed.

Assumptions:
  1. Two-dimensional flow (no spanwise variation).
  2. Incompressibility: ∇·u = 0.
  3. Irrotationality: ∇×u = 0, so a velocity potential φ exists (u = ∇φ).
  4. No viscosity, no turbulence closure, no wave breaking.
  5. Free surface is at atmospheric pressure; surface tension is neglected.
  6. Strict periodicity with wavelength λ and period τ.

Under these assumptions, the "steady travelling wave" problem is a classical
nonlinear boundary-value problem that possesses a two-parameter family of
solutions. Practical input cases fix two independent parameters among:
  { H/d, λ/d, τ√(g/d), current specification }.

------------------------------------------------------------------------------
B) FRAMES OF REFERENCE AND SOLVER COORDINATES
------------------------------------------------------------------------------

Three coordinate descriptions are used. Keeping them distinct is essential for
sign conventions, current definitions, and interpreting output:

  (x, y)   Stationary (laboratory/Earth) coordinates
           x positive in the propagation direction; y positive upward.
           Finite depth convention in this program: bed at y = 0, mean depth d
           (mean free surface at y = d).

  (ξ, y')  Travelling / mean-level coordinates
           ξ := x - c t      (steady in ξ for a permanent-form wave)
           y' := y - d       (mean free surface at y' = 0; bed at y' = -d)

  (X, Y)   Non-dimensional travelling coordinates used by the algebraic system
           X := k ξ
           Y := k y'
           with k := 2π/λ and kd := k d.

Unless explicitly stated otherwise, this theory section uses (X, Y) in the last,
non-dimensional sense above.

Velocities:
  Stationary-frame velocities are (u, v).
  Travelling-frame velocities are (U, V) with:
      U = u - c
      V = v

After convergence, stationary-frame time-derivatives are obtained using the
steady-wave identity in the travelling coordinate:
    ∂/∂t = -c ∂/∂x  (equivalently, ∂/∂t = -c ∂/∂ξ).

------------------------------------------------------------------------------
C) STREAM FUNCTION AND POTENTIAL (2D IRROTATIONAL FLOW)
------------------------------------------------------------------------------

In 2D incompressible flow a stream function ψ exists such that, in the k-based
non-dimensionalisation (X,Y) and with velocities scaled by √(g/k):

    Û =  ∂ψ/∂Y
    V̂ = -∂ψ/∂X

Irrotationality implies ψ is harmonic:
    ∂²ψ/∂X² + ∂²ψ/∂Y² = 0

A velocity potential φ exists as well:
    Û = ∂φ/∂X
    V̂ = ∂φ/∂Y

In 2D potential flow, φ and ψ are harmonic conjugates (Cauchy–Riemann):
    φ_X =  ψ_Y
    φ_Y = -ψ_X

Implementation note:
  The Fourier representation exploits this conjugacy: ψ is represented with cosine
  modes in X, and the corresponding potential uses sine modes (and vice-versa).

------------------------------------------------------------------------------
D) BOUNDARY CONDITIONS (STEADY TRAVELLING WAVE)
------------------------------------------------------------------------------

Let the free surface be Y = η(X) in the non-dimensional travelling coordinates.

Finite depth:
  Bed is at Y = -kd and is impermeable. In stream-function form this is equivalent
  to the bed being a streamline; the constant is chosen as:
      ψ(X, -kd) = 0

Deep water:
  There is no bed; instead the Fourier representation is required to decay as
  Y → -∞ (the kd → ∞ limit of the finite-depth basis).

Free-surface kinematic condition (streamline condition):
  In the travelling frame the free surface is steady and must coincide with a
  streamline:
      ψ(X, η(X)) = constant = -q

Free-surface dynamic condition (Bernoulli, p = 0 on the surface):
  In the same k-based non-dimensionalisation, Bernoulli reduces to:
      ½( ψ_X² + ψ_Y² ) + η(X) = R

Here q and R are global constants to be solved for (their dimensional counterparts
are recovered by the k-based scaling).

The nonlinear solve finds η(X), the Fourier coefficients, and the global constants
that satisfy these conditions (together with the selected current criterion).

------------------------------------------------------------------------------
E) FOURIER-APPROXIMATION / STREAM-FUNCTION REPRESENTATION
------------------------------------------------------------------------------

Fenton's Fourier approximation constructs a spectral representation that:
  • satisfies Laplace's equation exactly,
  • satisfies the bottom condition exactly,
  • leaves only the nonlinear free-surface conditions to be enforced by collocation.

In the common bed-based vertical coordinate (bed at y=0, mean depth d), a typical
finite-depth representation is:

  ψ(ξ,y) = -Ū y + Σ{j=1..N} B_j * sinh(j k y) / cosh(j k d) * cos(j k ξ)

  where ξ := x - c t is the travelling horizontal coordinate. In non-dimensional
  variables X = k ξ, the phase dependence is simply cos(jX).

where:
  • k = 2π/λ is the wavenumber,
  • Ū is the mean horizontal speed on any horizontal line in the travelling frame,
  • B_j are Fourier coefficients,
  • N is the truncation order.

This form automatically enforces ψ=0 on the bed (y=0) because sinh(0)=0.

---- Vertical coordinate shift used here --------------------------------------

This code works with the shifted solver coordinate (mean-surface origin):
  Y  = k (y - d)

so:
  mean surface: y=d   -> Y  = 0
  bed:          y=0   -> Y  = -kd, where kd = k d

Using the identity:
  sinh(jk y)/cosh(jk d) = sinh(j(Y+kd))/cosh(jkd)
                        = sinh(jY) + cosh(jY) tanh(jkd)

the finite-depth basis becomes exactly the combination used in Eqns() and Point():

  S_j(Y) = sinh(jY) + cosh(jY) tanh(jkd)

Deep-water limit:
  as kd → ∞, tanh(jkd) → 1 and S_j(Y) → exp(jY), giving the exponential basis
  used in the "Deep" branches.

---- Separation of the uniform-flow term -------------------------------------

The uniform-flow term (-Ū y) is not embedded inside the Fourier sums in this code.
Instead, its contributions appear explicitly:
  • in the streamline boundary condition (through z[7] and z[8]),
  • in the Bernoulli kinetic term (through the subtraction -z[7]).

This keeps the harmonic sums focused on the non-uniform wave-induced structure.

------------------------------------------------------------------------------
F) CURRENT DEFINITIONS: EULERIAN ū₁ VS MASS-TRANSPORT ū₂
------------------------------------------------------------------------------

With current, "mean current" and "wave speed" depend on how averaging and
reference frames are defined.

Depth-averaged mean velocities in the *stationary* frame:
  ū₁ : Eulerian mean current (time mean at a point, then depth-averaged)
  ū₂ : mass-transport (Stokes) mean current tied to net discharge

Travelling-frame mean speed Ū is related to ū₁ by:
  ū₁ = c - Ū

Net discharge Q relates to ū₂ by:
  ū₂ = c - Q/d

These relationships are enforced by the global residual equations:
  rhs[4] : z[5] + z[7] - z[4] = 0
  rhs[5] : z[1]*(z[6] + z[7] - z[4]) - z[8] = 0

Interpretation of variables in k-based scaling:
  z[4] = c √(k/g)
  z[5] = ū₁ √(k/g)
  z[6] = ū₂ √(k/g)
  z[7] = Ū √(k/g)

The period-specified family is closed by prescribing ū₁/√(gd) or ū₂/√(gd)
(CurrentCriterion selects which). rhs[6] implements that closure.

------------------------------------------------------------------------------
G) COLLOCATION SYSTEM AND WHY ONLY HALF A WAVE IS SOLVED
------------------------------------------------------------------------------

The solver exploits the symmetry of a periodic wave about crest and trough.
It solves on X ∈ [0, π] (crest to trough), using collocation nodes:

  X_m = m π / N,   m = 0..N

This is sufficient because the fields repeat with period 2π in X, and the chosen
Fourier structure (cosine/sine) enforces the correct parity of η, ψ, and φ.

Unknowns:
  • global scalars (kd, kH, period/celerity, current/discharge constants, R),
  • N+1 nodal surface elevations kη_m,
  • N Fourier coefficients B_j.

Equations:
  • 2(N+1) collocation equations (streamline + Bernoulli on the surface),
  • additional global closure equations (height, mean level, current, scalings),
giving a square nonlinear system of size num = 2N + 10.

------------------------------------------------------------------------------
H) HOW THE IMPLEMENTATION USES THE THEORY
------------------------------------------------------------------------------

At each collocation point:
  • ψ_Fourier(X_m, kη_m) is evaluated from the truncated series.
  • ũ = ∂ψ_Fourier/∂Y and ṽ = -∂ψ_Fourier/∂X are evaluated.

These represent deviations from the Eulerian current. The travelling-frame
relative velocity is:
  u - c = (ū₁ + ũ) - (ū₁ + Ū) = ũ - Ū
so the Bernoulli kinetic term uses (ũ - Ū)² + ṽ², exactly as coded via
(-z[7] + u)² + v².

The streamline condition on the free surface in shifted coordinates can be written:
  ψ_Fourier(X, kη) = Ū (kd + kη) - q
where q is a discharge constant in the solver’s k-based scaling. The legacy scalar
z[8] stores the combination (Ū*kd - q), so the residual is assembled as:
  psi - z[8] - z[7]*(kη) = 0.

After convergence:
  • the nodal kη_m values are converted to cosine-series coefficients for smooth
    evaluation of η(X),
  • Point() evaluates kinematics/pressure throughout the field, converting between
    k-based and depth-based nondimensionalisations and applying steady-wave
    relations to obtain accelerations.

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

IMPORTANT: We *do not* change any numerical formulas or iteration logic.
Only lifetime/ownership and cross-platform I/O robustness is improved.
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

/*
------------------------------------------------------------------------------
LEGACY MACROS (FENTON-STYLE DRIVER COMPATIBILITY)
------------------------------------------------------------------------------
The original multi-file program used a small header ("Headers.h") with many
preprocessor macros to keep the C code compact:

  - `Skip` / `skip` : consume the remainder of the current input line.
  - `Read` / `read` : scanf a value then Skip to end-of-line (two spellings kept).
  - `iff(x, Token)` : string-compare convenience that reads like a mini language.
  - `Is_finite` / `Is_deep` : branch on the chosen depth regime.
  - `HI` / `LO` : legacy printf field-width specifiers used in output tables.

These macros are retained verbatim to ensure:
  1) the *token stream* of numerical formulas remains identical to the legacy code,
  2) input parsing behaves exactly like the historical driver,
  3) the z[]/rhs[] indexing semantics are preserved.

Important cautions:
  - Macros do not respect scope; they can hide side effects.
  - `Readtext` assumes the line ends with '\n' and overwrites it with '\0'.
  - `Skip` hardcodes a 400-byte buffer size (`dummy[400]`), matching the legacy.

Do not "modernise" these macros unless you are prepared to revalidate the solver
numerically against known benchmark cases (e.g., Fenton's published tables).
------------------------------------------------------------------------------
*/

#define Skip             std::fgets(dummy, 400, Input1)
#define skip             Skip
#define read(x, y)       std::fscanf(Input1, "%" #y, &(x)); skip
#define pi               3.14159265358979324
#define iff(x, y)        if (std::strcmp((x), #y) == 0)
#define HI               " % 13.8f"
#define LO               " % 13.8f"
#define Readtext(x)      do { std::fgets((x), 400, Input1); (x)[std::strlen((x)) - 1] = '\0'; } while (0)
#define Read(x, y)       do { std::fscanf(Input1, "%" #y, &(x)); Skip; } while (0)
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
  - thread-unsafety (not relevant here: the solver is single-case, single-thread),
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
static int ns = 0;
static int nstep = 0;
static int num = 0;
static int number = 0;
static int Points = 0;
static int Surface_points = 0;

// Global arrays
static double **sol = nullptr;
static double *B = nullptr;
static double *coeff = nullptr;
static double *cosa = nullptr;
static double *rhs1 = nullptr;
static double *rhs2 = nullptr;
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
        Q  = depth-integrated discharge in the moving frame (depth-based)
        q  = auxiliary/legacy flux variable used by the original formulation

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
static double criter = 0.0;
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
static double sum = 0.0;
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
static float *vector(long nl, long nh)
{
  float *v = (float*)std::malloc((std::size_t)((nh - nl + 1 + NR_END) * sizeof(float)));
  return v - nl + NR_END;
}

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

static void free_ivector(int *v_, long nl, long /*nh*/)
{
  std::free((FREE_ARG)(v_ + nl - NR_END));
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
  - J is a dense Jacobian built by finite differences.
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

static void dsvdcmp(double **a, int m, int n, double w[], double **v)
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
    for (its = 1; its <= 30; its++) {
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
      if (its == 50) {
        std::fprintf(stderr, "no convergence in 30 dsvdcmp iterations");
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
 * where J is the numerically-differenced Jacobian matrix.
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

  // Perform decomposition
  dsvdcmp(a, m, n_, w, v);

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
        wmin = wmax * 1e-12
    This preserves the well-resolved modes while damping the ill-conditioned ones.

    Note: This is not a "tunable physics parameter"; it is a numerical-stability
    safeguard. Adjust only with regression testing.
  */

  // Set up: see p65 of Press et al.
  wmax = 0.0;
  for (i = 1; i <= n_; i++) if (w[i] >= wmax) wmax = w[i];
  wmin = wmax * 1.e-12;
  for (i = 1; i <= n_; i++) if (w[i] <= wmin) w[i] = 0.0;

  // Back substitute
  dsvbksb(a, w, v, m, n_, b, solution);

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

      // Carvalho (2006) GEP-evolved approximation:
      //   kd ≈ α / tanh( (6⁄5)ᵅ · α¹ᐟ² )
      //
      // In this code:
      //   α = σ² = (ω² d / g)  and  α = (kd) tanh(kd)
      //   σ = 2π / τ, with τ = T√(g/d) (dimensionless period input)
      //   α = σ²

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
      z[2] = 4 * pi * pi * Height;
    else
      z[2] = 2 * pi * Height;
    z[3] = 2 * pi;
    z[4] = 1.;
  }

  Is_finite {
    /*
      INITIAL CURRENT GUESS
      ---------------------
      The nonlinear steady-wave family depends on the specification of “current”.
      Fenton (1999) distinguishes:
        - Eulerian mean current ū1 = c - Ū  (Eq. 3.13)
        - Stokes / mass-transport current ū2 = c - Q/d (Eq. 3.14)

      Here, Current_criterion selects whether the *input* CurrentValue corresponds
      to ū1 (Euler) or ū2 (Stokes). The exact enforcement of the chosen criterion
      is performed inside Eqns() (see rhs[6] and Fenton 1999 Eq. 3.15–3.16).

      The lines below are ONLY a starting guess for Newton’s method; the subsequent
      iterations will correct these values to satisfy the full system.
    */
    if (Current_criterion == 1) {
      z[5] = Current * std::sqrt(z[2]);
      z[6] = 0.;
    } else {
      z[6] = Current * std::sqrt(z[2]);
      z[5] = 0.;
    }
  }

  Is_deep {
    z[5] = 0.;
    z[6] = 0.;
  }

  /*
    INITIAL GUESSES FOR MEAN FLOW / FLUX / BERNOULLI CONSTANT
    --------------------------------------------------------

    The initial values below are a *small-amplitude* consistent guess:
      - z[7] (U) is set equal to z[4] so that the mean flow in the wave frame
        matches the phase speed in the linear limit. This makes the moving-frame
        velocities small initially, stabilising the first Newton steps.

      - z[8] (q) is initialised to zero. q is a legacy auxiliary variable that
        becomes nonzero as the nonlinear solution develops and mass transport /
        set-down effects appear.

      - z[9] (r) is seeded as ½ U^2, consistent with a Bernoulli constant for
        a uniform flow state. Newton's method will then adjust it to satisfy the
        nonlinear dynamic condition on the free surface.

    These seeds matter only for convergence speed; the final converged solution
    is determined by Eqns() constraints.
  */
  z[7] = z[4];
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
 * z[3]   : kc / ω    period parameter used by the legacy code (ties to τ and dispersion)
 * z[4]   : c √(k/g)  wave speed parameter in k-based scaling
 * z[5]   : u1 √(k/g) Eulerian mean current in k-based scaling
 * z[6]   : u2 √(k/g) Stokes/mass-transport current in k-based scaling
 * z[7]   : U √(k/g)  mean fluid speed in the moving frame (k-based)
 * z[8]   : (Ū*kd-q)  auxiliary streamline-constant combination (k-based)
 * z[9]   : r         Bernoulli constant in k-based scaling
 *   z[10+m], m=0..N  = kη_m surface elevations at collocation points
 *   z[n+10+j], j=1..N = B_j Fourier coefficients in ψ expansion (Fenton 1999 Eq. 3.5)
 *
 * The core of this routine enforces the two free-surface boundary conditions
 * at N+1 collocation points over half a wave (crest-to-trough):
 *
 *   (1) ψ(X, η(X)) = -Q          (Fenton 1999 Eq. 3.3 / 3.7)
 *   (2) ½(ψ_X² + ψ_Y²) + η = R   (Fenton 1999 Eq. 3.4 / 3.8)
 *
 * A Newton solver (see Newton()) drives rhs[] → 0.
 *
 * Returns: Sum of squares of residual components (useful as a diagnostic).
 */
static double Eqns(double *rhs)
{

  /*
    ---------------------------------------------------------------------------
    rhs[] INDEX MAP (QUICK REFERENCE)
    ---------------------------------------------------------------------------

    For N Fourier modes (variable `n`), the total system size is:
        num = 2*N + 10

    Indices:

      rhs[1]  : (finite depth) kH - kd*(H/d) = 0
      rhs[2]  : ties internal continuation height variable to kH (depends on Case)
      rhs[3]  : z[4]*z[3] - 2π = 0  (period/celerity scaling identity)
      rhs[4]  : Euler current definition: u1 + U - c = 0
      rhs[5]  : Stokes current definition (mass transport / discharge relation)
      rhs[6]  : user-specified current magnitude for chosen criterion
      rhs[7]  : mean water level condition (mean of kη_m over wavelength is zero)
      rhs[8]  : wave height condition: kη_crest - kη_trough - kH = 0

      For each collocation point m = 0..N:

        rhs[9 + m]          : streamfunction BC on free surface (kinematic):
                              ψ(X_m, η_m) + Q = 0  (in program variables)

        rhs[(N+10) + m]     : Bernoulli BC on free surface (dynamic):
                              ½(U^2+V^2) + η_m - r = 0  (in program variables)

    This organisation matches the original driver and is relied upon by Newton()
    when forming the Jacobian and updating z[].

    ---------------------------------------------------------------------------
  */
  int i, j, m, nm;
  double coshkd, e, sinhkd, Sum, u_, v_;

  /*
    ---------------------------------------------------------------------------
    “GLOBAL” (NON-COLLOCATION) CONSTRAINT EQUATIONS
    ---------------------------------------------------------------------------

    These residual components implement the additional equations needed to close
    the steady-wave system beyond the 2(N+1) free-surface collocation equations.

    Fenton (1999) counts 2N+5 unknowns for the canonical formulation (Eq. 3.7–3.11
    plus Eq. 3.15/3.16 variants). This legacy driver uses a slightly expanded
    unknown set (2N+10 entries in z[]) for convenience and historical reasons,
    but the constraints map directly to the same physical requirements:

      - known H/d (relative height),
      - known λ/d (wavelength case) OR known τ*sqrt(g/d) and current (period case),
      - chosen definition of mean current (Eulerian ū1 or Stokes ū2),
      - consistency relations between the various nondimensional groups.

    See Fenton (1999) §3.1.2 and equations (3.9)–(3.16).
  */

  // rhs[1]: enforce kH = kd*(H/d) in finite depth.
  // This is a direct consequence of the definitions:
  //   kd = k d,   H/d = Hoverd,   kH = z[2]  =>  z[2] - z[1]*Hoverd = 0.
  // Compare with the wave-height specification logic around Fenton (1999) Eq. (3.10),
  // where a known H/d is used while kd is unknown in τ-specified problems.
  Is_deep   rhs[1] = 0.0;
  Is_finite rhs[1] = z[2] - z[1] * Hoverd;

  // rhs[2]: connects the “height parameter” used internally by the legacy code to
  // the user-specified wave scale, depending on whether λ or τ is specified.
  //
  // - Wavelength case: k = 2π/λ is known once λ/d is specified, so kH can be related
  //   directly to H/λ (cf. Fenton 1999 Eq. 3.11 for kd when λ/d known).
  //
  // - Period case: τ is specified instead of λ; the system uses z[3] and z[4] to
  //   represent period/celerity consistency, and this equation ties the chosen
  //   “height” continuation variable to z[3] (legacy driver retained exactly).
  iff(Case, Wavelength)
    rhs[2] = z[2] - 2.0 * pi * height;
  else
    rhs[2] = z[2] - height * z[3] * z[3];

  // rhs[3]: period/celerity consistency (linear scaling identity preserved).
  // In the linear limit, τ*sqrt(gk) = 2π / sqrt(tanh(kd)); the code stores the
  // corresponding factors as z[3] and z[4] and enforces: z[4]*z[3] = 2π.
  rhs[3] = z[4] * z[3] - pi - pi;

  // rhs[4]: Eulerian current definition ū1 = c - Ū  (Fenton 1999 Eq. 3.13).
  // Rearranged to residual form: ū1 + Ū - c = 0.
  rhs[4] = z[5] + z[7] - z[4];

  // rhs[5]: Stokes/mass-transport current definition ū2 = c - Q/d (Fenton 1999 Eq. 3.14),
  // expressed in the program’s nondimensional variables.
  rhs[5] = z[1] * (z[6] + z[7] - z[4]) - z[8];

  for (i = 1; i <= n; i++) {
    coeff[i] = z[n + i + 10];
    Is_finite Tanh[i] = std::tanh(i * z[1]);
  }

  // rhs[6]: enforce the USER-SPECIFIED current magnitude for the chosen criterion.
  // In the τ-specified (period) problems, Fenton (1999) shows that the equations
  // may be closed by specifying τ*sqrt(g/d) and one of ū1/sqrt(gd) or ū2/sqrt(gd)
  // (Eq. 3.15 or Eq. 3.16).
  //
  // In this implementation, Current is read/stored as that *depth-based*
  // nondimensional value (current / sqrt(g d)). The solver’s native scaling is
  // wavenumber-based; conversion introduces sqrt(kd)=sqrt(z[1]).
  Is_finite rhs[6] = z[Current_criterion + 4] - Current * std::sqrt(z[1]);

  Is_deep {
    iff(Case, Wavelength)
      rhs[6] = z[Current_criterion + 4] - Current;
    else
      rhs[6] = z[Current_criterion + 4] - Current * z[3];
  }

  // rhs[7]: mean water level / mean depth condition.
  // The Fenton (1999) formulation includes an equation for mean depth using a
  // trapezoidal rule over the collocation points (Eq. 3.9). This legacy driver
  // uses an equivalent condition expressed directly in terms of the collocated
  // surface elevations kη_m: the (periodic) trapezoidal sum of kη_m is enforced
  // to be zero, i.e. η is measured relative to the mean level.
  rhs[7] = z[10] + z[n + 10];
  for (i = 1; i <= n - 1; i++) rhs[7] = rhs[7] + z[10 + i] + z[10 + i];

  // rhs[8]: wave height condition in k-scaling: kη_crest - kη_trough = kH.
  // This matches the wave-height specification equation (Fenton 1999 Eq. 3.10)
  // when expressed in terms of kH.
  rhs[8] = z[10] - z[n + 10] - z[2];

  /*
    ---------------------------------------------------------------------------
    COLLOCATION LOOP OVER FREE SURFACE POINTS
    ---------------------------------------------------------------------------

    At each collocation index m (0..N), the free surface elevation is represented
    by the unknown value:
        kη_m = z[10 + m]

    We evaluate the truncated Fourier representation of ψ and its derivatives at:
        X = X_m = m π / N
        Y = kη_m   (because in k-based scaling the mean free surface is Y=0)

    The basis functions depend on depth regime:

      Finite depth:
        The vertical structure that satisfies Laplace's equation and the bed
        boundary condition is expressed using combinations of sinh/cosh and
        tanh(j kd). In Fenton (1999), this corresponds to Eq. (3.5).

      Deep water:
        The kd → ∞ limit reduces the basis to exp(j Y) (Eq. (3.6)).

    The sums below compute:
        ψ(X_m, η_m),  U = ∂ψ/∂Y,  V = -∂ψ/∂X
    up to the truncation order N, then enforce the two surface boundary
    conditions by setting the corresponding residual entries.
    ---------------------------------------------------------------------------
  */
  for (m = 0; m <= n; m++) {
    psi = 0.;
    u_ = 0.;
    v_ = 0.;
    for (j = 1; j <= n; j++) {
      /*
        TRIG TABLE INDEX nm
        -------------------
        X_m = m π / N. Therefore j X_m = (m*j) π / N.
        Because cos/sin are 2π-periodic, we can reduce the phase modulo 2π.
        The tables cosa[]/sina[] are stored as cos(p π / N), sin(p π / N) for
        p = 0..2N, so the appropriate index is:

            nm = (m*j) mod (2N)   where 2N == (n+n) in this code.
      */
      nm = (m * j) % (n + n);
      Is_finite {
        /*
          Vertical basis evaluation at the free surface:
            Y = kη_m = z[10+m]

          We need sinh(jY) and cosh(jY). Using exponentials improves speed and
          was common in the original code:

            sinh(jY) = (e^{jY} - e^{-jY})/2
            cosh(jY) = (e^{jY} + e^{-jY})/2

          Here:
            e = exp(jY)
            1/e = exp(-jY)

          Note: For steep waves, Y near crest can be positive, making exp(jY)
          large for high j. The truncation order N must therefore be chosen with
          care for numerical stability.
        */
        e = std::exp(j * (z[10 + m]));
        sinhkd = 0.5 * (e - 1. / e);
        coshkd = 0.5 * (e + 1. / e);
        psi = psi + coeff[j] * (sinhkd + coshkd * Tanh[j]) * cosa[nm];
        u_ = u_ + j * coeff[j] * (coshkd + sinhkd * Tanh[j]) * cosa[nm];
        v_ = v_ + j * coeff[j] * (sinhkd + coshkd * Tanh[j]) * sina[nm];
      }
      Is_deep {
        e = std::exp(j * (z[10 + m]));
        psi = psi + coeff[j] * e * cosa[nm];
        u_ = u_ + j * coeff[j] * e * cosa[nm];
        v_ = v_ + j * coeff[j] * e * sina[nm];
      }
    }
    // Collocation equation #1 (kinematic / streamline condition on the free surface)
    //
    // In the travelling frame the free surface is a streamline, so ψ is constant there.
    // With the solver vertical coordinate Y = k(y - d) (mean surface at Y=0, bed at Y=-kd),
    // the Fourier sum `psi` represents ψ_Fourier(X, η), where η = kη_phys and X = k(x - c t).
    // The free-surface streamline condition may be expressed in the algebraic form:
    //
    //   psi(X, kη) = Ū * (kd + kη) - q
    //
    // where Ū is the travelling-frame mean speed (z[7]) and q is a discharge constant.
    // The scalar z[8] stores Ū*kd - q, so the residual is assembled as:
    //   psi - z[8] - z[7]*(kη) = 0.
    rhs[m + 9] = psi - z[8] - z[7] * z[m + 10];

    // Collocation equation #2 (dynamic / Bernoulli condition on free surface)
    //   1/2(ψ_X^2 + ψ_Y^2) + η = R   (k-based non-dimensional Bernoulli)
    // In dimensionless form this is Eq. (3.8) of Fenton (1999).
    // The velocities in the moving frame enter via (U,V) = (ψ_Y, -ψ_X) (Eq. 3.1),
    // and the squared speed term is assembled here as (Ū - U)^2 + V^2.
    rhs[n + m + 10] = 0.5 * (std::pow((-z[7] + u_), 2.) + v_ * v_) + z[m + 10] - z[9];
  }

  for (j = 1, Sum = 0.; j <= num; j++) Sum += rhs[j] * rhs[j];
  return Sum;
}

/**
 * Perform one Newton iteration on the nonlinear system defined by Eqns().
 *
 * This implementation follows the strategy described in Fenton (1988):
 *   - Evaluate rhs at the current state z[]
 *   - Build a numerical Jacobian by forward differencing each component z[i]
 *   - Solve the linear system J Δz = -rhs for the correction Δz
 *     using an SVD-based solver for robustness
 *   - Update z ← z + Δz
 *
 * The return value is the mean absolute correction applied to the free-surface
 * unknowns (z[10..10+N]), which is used as a convergence indicator.
 */
static double Newton(int count)
{
  double **a = nullptr;
  double *rhs = nullptr;
  double *x = nullptr;
  double h, sum_;

  int i, j;

  Eqns(rhs1);

  /*
    JACOBIAN CONSTRUCTION BY FINITE DIFFERENCES
    ------------------------------------------
    The Jacobian matrix J = ∂rhs/∂z is not derived analytically here. Instead, we
    approximate each column i by a forward difference:

        J[:, i] ≈ ( rhs(z + h e_i) - rhs(z) ) / h

    where e_i is the i-th unit vector.

    Step size selection:
      h = 0.01 * z[i]        (1% relative perturbation)
      but if |z[i]| < 1e-4 then h = 1e-5 (absolute floor)

    Rationale:
      - A relative step captures scale differences across parameters.
      - The absolute floor avoids h=0 when z[i] is (near) zero.
      - Forward differences are used because they require only one extra rhs
        evaluation per column (cheaper than central differences), which is
        acceptable given the moderate system sizes typical of this solver.

    Trade-offs:
      - Too small h -> subtraction cancellation (roundoff dominates)
      - Too large h -> poor linearisation (truncation error dominates)

    The chosen heuristic matches the legacy code and has proven robust across the
    intended operating range.
  */

  if (count >= 1) {
    ++count;
    rhs = dvector(1, num);
    x = dvector(1, num);
    a = dmatrix(1, num, 1, num);
  }

  for (i = 1; i <= num; i++) {
    h = 0.01 * z[i];
    if (std::fabs(z[i]) < 1.e-4) h = 1.e-5;
    z[i] = z[i] + h;
    Eqns(rhs2);
    z[i] = z[i] - h;
    rhs[i] = -rhs1[i];
    for (j = 1; j <= num; j++) a[j][i] = (rhs2[j] - rhs1[j]) / h;
  }

  Solve(a, rhs, (int)num, (int)num, x, (int)num, (int)num);

  for (i = 1; i <= num; i++) z[i] += x[i];

  /*
    CONVERGENCE METRIC USED BY THE LEGACY DRIVER
    --------------------------------------------
    The original program monitored convergence primarily via the mean absolute
    correction to the *free-surface* degrees of freedom (kη_m). Those are stored
    in z[10..10+N].

    This is a sensible choice because:
      - surface elevations are directly tied to the nonlinear boundary conditions,
      - they typically exhibit the largest sensitivity in steep-wave regimes,
      - they provide an intuitive "physical" measure of how much the wave shape
        is still changing.

    We compute:
        error = (1/N) Σ_{m=0..N} |Δ(kη_m)|

    and compare it against the tolerance criterion (crit or tighter on final step).
  */
  for (sum_ = 0., i = 10; i <= n + 10; i++) sum_ += std::fabs(x[i]);
  sum_ /= n;

  free_dvector(rhs, 1, num);
  free_dvector(x, 1, num);
  free_dmatrix(a, 1, num, 1, num);

  return sum_;
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

  std::fprintf(file, "\n%s\n", Method);
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
  Readtext(Title);
  iff(Title, FINISH) return 0;
  Read(MaxH, lf);

  /*
    DEPTH-REGIME SELECTION VIA LEGACY SIGN CONVENTION
    -------------------------------------------------
    The historical input format uses the sign of the "Height/Depth" line to
    select finite vs deep water:

      MaxH >= 0  => finite depth
                   MaxH is interpreted as H/d

      MaxH <  0  => deep water
                   MaxH is interpreted as -H/L  (note the minus sign)

    This convention is preserved for compatibility with existing data.dat files
    and with the documented interface in Instructions.pdf (Table 6-1).

    Internally, the variable `Height` is a continuation/control parameter used by
    init() and Eqns() to seed the solver. Its meaning depends on the case:

      Finite depth, Case=Wavelength:
        Height = (H/d) / (L/d) = H/L

      Finite depth, Case=Period:
        Height = (H/d) / (T√(g/d))^2  (legacy scaling used to relate kH to τ)

      Deep water:
        Height = H/L (because MaxH stores -H/L)

    This is why both `Height` and `Hoverd` appear in the code: Hoverd is the
    physically meaningful H/d at each continuation step, while Height is the
    legacy "driver variable" that keeps the original formulas intact.
  */

  if (MaxH >= 0.) {
    std::strcpy(Depth, "Finite");
    std::fscanf(Input1, "%s", Case);
    Skip;
    iff(Case, Wavelength) {
      Read(L, lf);
      Height = MaxH / L;
    }
    iff(Case, Period) {
      Read(T, lf);
      Height = MaxH / (T * T);
    }
  }

  if (MaxH < 0.) {
    std::strcpy(Depth, "Deep");
    std::fscanf(Input1, "%s", Case);
    Skip;
    Skip;
    Height = -MaxH;
  }

  Read(Current_criterion, d);
  Read(Current, lf);
  if (Current_criterion == 1) std::strcpy(Currentname, Current1);
  if (Current_criterion == 2) std::strcpy(Currentname, Current2);

  Read(n, d);
  Read(nstep, d);

  Input_Data_block(monitor);

  if (std::strcmp(Theory, "Stokes") == 0) {
    iff(Case, Wavelength)
      if (L > 10.) {
        std::printf("\nThe dimensionless wavelength is greater than 10.");
        std::printf("\nStokes theory should not be applied. Exiting.");
        std::exit(1);
      }
    iff(Case, Period)
      if (T > 10.) {
        std::printf("\nThe dimensionless period is greater than 10.");
        std::printf("\nStokes theory should not be applied. Exiting.");
        std::exit(1);
      }
  }

  // Convergence criteria (embedded defaults; Convergence.dat is not required)
  // Values correspond to the original Convergence.dat shipped with this codebase.
  number = kConvergenceMaxIterations;
  crit   = kConvergenceCriterion;


  // Output sampling resolution (embedded; formerly read from Points.dat)
  Surface_points = kSurfacePointsDefault;
  Nprofiles      = kNprofilesDefault;
  Points         = kProfilePointsDefault;


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
  std::fprintf(file, "\n%s\n", Method);
  std::fprintf(file, "\n# Height/Depth:%6.5f, %3.0lf%% of the maximum of H/d =%6.5f for this length:",
               z[2] / z[1], z[2] / z[1] / Highest * 100., Highest);
  std::fprintf(file, "\n# Length/Depth:%7.5f", 2 * pi / z[1]);
  std::fprintf(file, "\n# Dimensionless Period T*sqrt(g/d):%7.5f", z[3] / std::sqrt(z[1]));
  std::fprintf(file, "\n# Current criterion: %s,  Dimensionless value:%6.5lf\n", Currentname, Current);
}

static void Title_block_deep(std::FILE* file)
{
  Highest = 0.141063;
  std::fprintf(file, "# %s", Title);
  std::fprintf(file, "\n%s\n", Method);
  std::fprintf(file, "\n# Height/Length:%6.5f, %3.0lf%% of the maximum of H/L =%6.5f",
               z[2] / 2 / pi, (z[2] / 2 / pi) / Highest * 100., Highest);
  std::fprintf(file, "\n# Dimensionless Period T*sqrt(g/L):%7.5f", z[3] / std::sqrt(2 * pi));
  std::fprintf(file, "\n# Current criterion: %s,  Dimensionless value:%6.5lf\n", Currentname, Current);
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
  double coshdelta, sinhdelta;

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
  y = 1. + Y_ / kd;

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
    Is_finite {
      coshdelta = std::cosh(j * Y_);
      sinhdelta = std::sinh(j * Y_);
      C = coshdelta + sinhdelta * Tanh[j];
      S = sinhdelta + coshdelta * Tanh[j];
    }
    Is_deep C = S = std::exp(j * Y_);
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

static double MeanSquareBedOrbitalVelocity(int nph = 720)
{
  // Only defined/used for finite depth in this legacy program.
  // In deep water there is no physical bed at finite depth to evaluate.
  if (std::strcmp(Depth, "Finite") != 0) return 0.0;
  if (kd <= 0.0) return 0.0;

  // Definition used here:
  //   ub² = < (ub - ū₁)² >
  // with phase-averaging over one wave period / wavelength.
  const int N = std::max(36, nph);

  const double twoPi = 2.0 * pi;

  // In finite depth, Point() expects Y_ = k(y - d). At the bed y = 0 => Y_ = -k d = -kd.
  const double Y_bed = -kd;

  double sumsq = 0.0;

  for (int i = 0; i < N; ++i) {
    const double X = twoPi * (static_cast<double>(i) / static_cast<double>(N));

    // Evaluates u, v, pressures, etc. at (X, Y_bed). Results are stored in globals.
    Point(X, Y_bed);

    // Point() returns velocities in the depth-based non-dimensionalisation (scaled by √(g d)).
    // Convert to the solver's native k-based non-dimensionalisation (scaled by √(g/k)):
    const double u_abs_gk = u * std::sqrt(kd);

    // Remove the Eulerian mean current u1 (k-based variable is z[5]) to obtain orbital velocity.
    const double u_orb_gk = u_abs_gk - z[5];

    sumsq += u_orb_gk * u_orb_gk;
  }

  return sumsq / static_cast<double>(N);
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

  std::fprintf(monitor, "\n\n# Solution summary:\n\n");
  Is_finite Title_block(monitor);
  Is_deep Title_block_deep(monitor);

  Is_finite Title_block(Solution);
  Is_deep Title_block_deep(Solution);

  Is_finite {
    /*
      DERIVING "ENGINEERING" VARIABLES FROM THE SOLVER STATE (FINITE DEPTH)
      --------------------------------------------------------------------
      After Newton convergence at the final continuation step, z[] contains the
      solution in the solver's k-based nondimensionalisation.

      We now compute commonly used dimensional groups (still non-dimensional but
      in depth-based scaling) for reporting:

        kd   = k d
        L/d  = 2π / kd
        H/d  = (kH) / (kd)
        T√(g/d) = z[3] / √(kd)
        c/√(g d) = z[4] / √(kd)

      and corresponding currents:
        ū1/√(g d) = z[5] / √(kd)
        ū2/√(g d) = z[6] / √(kd)
        U/√(g d)  = z[7] / √(kd)

      Flux and Bernoulli constants also require kd^(3/2) factors because of the
      scaling of streamfunction/potential with √(g/k³).
    */
    kd = z[1];
    L = 2 * pi / z[1];
    H = z[2] / z[1];
    T = z[3] / std::sqrt(z[1]);
    c = z[4] / std::sqrt(z[1]);
    ce = z[5] / std::sqrt(z[1]);
    cs = z[6] / std::sqrt(z[1]);
    ubar = z[7] / std::sqrt(z[1]);
    /*
      DISCHARGE (VOLUME FLUX) IN THE MOVING FRAME
      ------------------------------------------
      In Fenton's notation, Q is the depth-integrated volume flux per unit width
      in the frame moving with the wave. In the solver variables, z[8] represents
      a k-based flux-related term, so the conversion to depth-based scaling
      introduces kd^(3/2):

          Q = U - q / (kd)^(3/2)

      where U is the mean flow in the wave frame (depth-based).
    */
    Q = ubar - z[8] / std::pow(kd, 1.5);
    R = 1 + z[9] / z[1];

    pulse = z[8] + z[1] * z[5];
    ke = 0.5 * (z[4] * pulse - z[5] * Q * std::pow(kd, 1.5));

    pe = 0;
    for (i = 1; i <= n; ++i) pe += 0.25 * std::pow(Y[i], 2);

    // -------------------------------------------------------------------------
    // Integral-wave auxiliary quantity used in Fenton-style expressions:
    //   ub2_internal = 2 r - c^2 (in the program's k-based non-dimensionalisation)
    //
    // This value appears inside the original formulas for radiation stress (Sxx)
    // and wave power (F). We preserve that usage unchanged to keep those formulas
    // identical to the supplied codebase.
    // -------------------------------------------------------------------------
    const double ub2_internal = 2. * z[9] - z[4] * z[4];

    sxx = 4. * ke - 3. * pe + ub2_internal * z[1] + 2. * z[5] * (z[7] * z[1] - z[8]);
    f = z[4] * (3. * ke - 2. * pe) + 0.5 * ub2_internal * (pulse + z[4] * z[1]) + z[4] * z[5] * (z[7] * z[1] - z[8]);

    // -------------------------------------------------------------------------
    // Reported mean square *orbital* bed velocity ub2:
    //   ub2 = < (ub - u1)^2 >
    // where u1 is the Eulerian current and <·> denotes phase averaging.
    // -------------------------------------------------------------------------
    ub2 = MeanSquareBedOrbitalVelocity();
    q = z[7] * z[1] - z[8];
    r = z[9] + z[1];
    s = sxx - 2. * z[4] * pulse + (z[4] * z[4] + 0.5 * z[1]) * z[1];
  }

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
  Results("# Eulerian current                  (u1)", z[5], ce);
  Results("# Stokes current                    (u2)", z[6], cs);
  Results("# Mean fluid speed in frame of wave  (U)", z[7], ubar);
  Results("# Volume flux due to waves           (q)", z[8], z[8] / std::pow(kd, 1.5));
  Results("# Bernoulli constant                 (r)", z[9], z[9] / kd);

  Is_finite {
    Results("# Volume flux                        (Q)", Q * std::pow(kd, 1.5), Q);
    Results("# Bernoulli constant                 (R)", R * kd, R);
    Results("# Momentum flux                      (S)", s, s / kd / kd);
    Results("# Impulse                            (I)", pulse, pulse / std::pow(kd, 1.5));
    Results("# Kinetic energy                     (T)", ke, ke / kd / kd);
    Results("# Potential energy                   (V)", pe, pe / kd / kd);
    Results("# Mean square of bed velocity      (ub2)", ub2, ub2 / kd);
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
    "    N                = number of Fourier components (Instructions.pdf §6.1.5)\n"
    "    Steps            = number of height steps to reach the final H/d (Instructions.pdf §6.1.6)\n\n"
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
    return end && *end == '\0';
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
      std::fprintf(stderr, "No case found in data.dat (FINISH encountered).\n");
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

  // ---------------------------------------------------------------------------
  // Allocate solver work arrays (NR-style 1-based indexing, but RAII-managed)
  // ---------------------------------------------------------------------------
  num = 2 * n + 10;

  const double dhe = Height / nstep;   // increment in “height” variable used by the legacy algorithm
  const double dho = MaxH   / nstep;   // increment in H/d used for step-wise continuation

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
  auto rhs2_owner  = make_dvector_owner(1, num);         rhs2 = rhs2_owner.get();
  auto coeff_owner = make_dvector_owner(0, n);           coeff = coeff_owner.get();
  auto cosa_owner  = make_dvector_owner(0, 2 * n);       cosa = cosa_owner.get();
  auto sina_owner  = make_dvector_owner(0, 2 * n);       sina = sina_owner.get();
  auto sol_owner   = make_dmatrix_owner(0, num, 1, 2);   sol = sol_owner.get();
  auto B_owner     = make_dvector_owner(1, n);           B = B_owner.get();
  auto Tanh_owner  = make_dvector_owner(1, n);           Tanh = Tanh_owner.get();

  // ---------------------------------------------------------------------------
  // Continuation in wave height (Fenton 1988 continuation/robustness strategy)
  // ---------------------------------------------------------------------------
  /*
    HEIGHT STEPPING / CONTINUATION STRATEGY
    --------------------------------------
    Solving directly at the target H/d can fail when the wave is steep because
    Newton's method is local: it requires a sufficiently good initial guess.

    The standard remedy in steady-wave solvers is *continuation*:
      1) Start from a small-amplitude wave where linear theory provides a good
         initial guess (init()).
      2) Increase the target height gradually in `nstep` increments.
      3) At each increment, use the previous converged solution as a predictor
         for the next one. This code uses a simple secant predictor:
             z_pred = 2*z_prev - z_prevprev
         which approximates the tangent direction along the solution branch.

    Variables:
      - `height` : legacy internal continuation parameter used in Eqns() branches
      - `Hoverd` : physically meaningful H/d at the current step

    On the *final* step (ns == nstep), the solver tightens the tolerance from
    crit (typically 1e-9) to 1e-10, matching the legacy program behaviour.
  */
  for (ns = 1; ns <= nstep; ns++) {
    height = ns * dhe;   // legacy “height parameter” used internally by init/Eqns branches
    Hoverd = ns * dho;   // H/d at this continuation step

    std::fprintf(monitor, "\n\nHeight step %2d of %2d\n", ns, nstep);

    // Step 1: small-amplitude initial guess (linear / weakly nonlinear start)
    if (ns <= 1) {
      init();
    } else {
      // Step >1: extrapolation from two previous solutions (simple secant continuation)
      for (int i = 1; i <= num; i++) z[i] = 2.0 * sol[i][2] - sol[i][1];
    }

    // Newton iterations for this step
    for (int iter = 1; iter <= number; iter++) {
      std::fprintf(monitor, "\nIteration%3d:", iter);

      const double error = Newton(iter);

      std::fprintf(stdout, " Mean of corrections to free surface: %8.1e", error);

      // Tighten the final step (legacy behavior preserved)
      criter = (ns == nstep) ? 1.0e-10 : crit;

      if ((error < criter * std::fabs(z[1])) && iter > 1) break;

      if (iter == number) {
        std::fprintf(stdout,
          "\nNote that the program still had not converged to the degree specified\n"
        );
      }

      // Extrapolation bookkeeping: retain last two solutions for the secant predictor
      if (ns == 1) {
        for (int i = 1; i <= num; i++) sol[i][2] = z[i];
      } else {
        for (int i = 1; i <= num; i++) {
          sol[i][1] = sol[i][2];
          sol[i][2] = z[i];
        }
      }
    }

    // -------------------------------------------------------------------------
    // POST-PROCESS: RECOVER COSINE-SERIES COEFFICIENTS FOR kη(X)
    // -------------------------------------------------------------------------
    Y[0] = 0.0;
    for (int j = 1; j <= n; ++j) {
      B[j] = z[j + n + 10];

      // Discrete cosine transform-like reconstruction of surface elevation Fourier
      // coefficients E_j (stored here in Y[j]) from the collocation values kη_m.
      //
      // This is a “slow Fourier transform” consistent with Fenton-style post-processing:
      // it forms the cosine-series that is later evaluated by Surface(X).
      sum = 0.5 * (z[10] + z[n + 10] * std::pow(-1.0, static_cast<double>(j)));
      for (int m = 1; m <= n - 1; ++m) {
        sum += z[10 + m] * cosa[(m * j) % (n + n)];
      }
      Y[j] = 2.0 * sum / n;
    }
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

  std::printf("\nFinished\n");
  return 0;
}