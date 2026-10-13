/* ==============================================================================
 *  ENGINEERING TECHNICAL REFERENCE & THEORETICAL FORMULATION
 * ==============================================================================
 *  PROGRAM:      Fenton Nonlinear Wave Suite - native Windows GUI solver
 *  FILE:         fenton_gui.cpp
 *  METHOD:       Stream-function / Fourier collocation method for steady,
 *                two-dimensional, finite-amplitude gravity waves.
 *  REFERENCE:    Rienecker-Fenton / Fenton numerical-wave formulation,
 *                as implemented in this source file.
 * ==============================================================================
 *
 *  1. SCOPE OF THIS IMPLEMENTATION
 *  -----------------------------------------------------------------------------
 *  This executable computes steady, periodic, finite-amplitude gravity waves in
 *  finite water depth, with optional collinear current.  The calculation follows
 *  the implemented reference formulation: the Fourier basis satisfies Laplace's
 *  equation and the bed condition analytically, while Newton iteration solves the
 *  nonlinear free-surface streamline condition, Bernoulli condition, wave-height
 *  constraint, period relation, current closure and global flux identities.
 *
 *  The central engineering mapping is
 *
 *      L = F(H, T, d, U_c)
 *
 *  where H is crest-to-trough wave height, T is the apparent fixed-frame period,
 *  d is still-water depth, U_c is the imposed collinear current under the selected
 *  current convention, and L is the dynamically consistent wavelength.  In the
 *  Fenton formulation a steady wave train is fundamentally defined by H, d and L;
 *  when T is specified instead of L, the current or wave-speed convention is part
 *  of the closure because the observed period is Doppler-shifted by current.
 *
 *  2. COORDINATES, NONDIMENSIONAL VARIABLES AND VELOCITIES
 *  -----------------------------------------------------------------------------
 *  The moving wave-frame phase coordinate and vertical coordinate are
 *
 *      X = k(x - c t),        Y = k y,        k = 2*pi/L,
 *
 *  with y = 0 at still-water level and y = -d at the bed.  One wavelength is
 *  X in [0, 2*pi], and the bed is Y = -k d.  The free surface is represented by
 *
 *      Y = zeta(X) = k eta(X).
 *
 *  The nondimensional stream function is
 *
 *      Psi = psi * sqrt(k^3/g),
 *
 *  and the wave-frame velocities are
 *
 *      U_hat = U * sqrt(k/g) = dPsi/dY,
 *      V_hat = V * sqrt(k/g) = -dPsi/dX.
 *
 *  3. STREAM-FUNCTION / FOURIER REPRESENTATION
 *  -----------------------------------------------------------------------------
 *  The finite-depth Fourier representation used by the residual equations is
 *
 *      Psi(X,Y) = -Ubar*(Y + kd)
 *                 + sum_{j=1..N} B_j sinh[j(Y+kd)]/cosh(jkd) cos(jX),
 *
 *  where Ubar* = Ubar sqrt(k/g).  At the bed, Y = -kd, both the linear term and
 *  all hyperbolic sine terms vanish, so the impermeable-bed condition is already
 *  satisfied before the nonlinear equations are assembled.
 *
 *  The hyperbolic ratios use exponential expressions without subtracting
 *  large, nearly equal hyperbolic terms:
 *
 *      S_j(Y) = [exp(jY) - exp(-jY-2jkd)] / [1 + exp(-2jkd)],
 *      C_j(Y) = [exp(jY) + exp(-jY-2jkd)] / [1 + exp(-2jkd)].
 *
 *  For large jkd, tanh(jkd) tends to 1 and both functions tend to exp(jY).
 *
 *  4. STATE VECTOR z[i]
 *  -----------------------------------------------------------------------------
 *  This file implements the Fenton-style 1-based state vector.  Index 0
 *  is intentionally unused so that the residual equations map directly to the
 *  documented z[i] notation:
 *
 *      z[1]  = kd
 *      z[2]  = kH
 *      z[3]  = T sqrt(gk)
 *      z[4]  = c sqrt(k/g)
 *      z[5]  = Eulerian current variable ubar_1 sqrt(k/g)
 *      z[6]  = Stokes / mass-transport current variable ubar_2 sqrt(k/g)
 *      z[7]  = mean wave-frame velocity Ubar sqrt(k/g)
 *      z[8]  = q sqrt(k^3/g), with q = Ubar d - Q
 *      z[9]  = r k/g, with r = R - g d
 *      z[10] ... z[N+10]       = free-surface ordinates zeta_m = k eta_m
 *      z[N+11] ... z[2N+10]    = Fourier coefficients B_j
 *
 *  Therefore, for Fourier order N, the active vector length is
 *
 *      num = 2N + 10.
 *
 *  The requested order is N = 50. The solver reduces the order when the
 *  equilibrated Jacobian exceeds its binary64 conditioning limit. The active
 *  residual system has 2N + 10 equations for the order actually used.
 *
 *  5. RESIDUAL SYSTEM IMPLEMENTED IN eqns()
 *  -----------------------------------------------------------------------------
 *  At convergence every component of F(z) is zero.  The first eight residuals are
 *  the global scalar constraints; the next N+1 enforce the free-surface streamline
 *  condition; the final N+1 enforce Bernoulli's equation at the same collocation
 *  nodes X_m = m*pi/N, m = 0..N.
 *
 *      r1 = z2 - z1(H/d)
 *      r2 = z2 - Hs z3^2
 *      r3 = z4 z3 - 2*pi
 *      r4 = z5 + z7 - z4
 *      r5 = z1(z6 + z7 - z4) - z8
 *      r6 = z[c+4] - U_c sqrt(z1)
 *      r7 = z10 + z[N+10] + 2 sum_{i=1..N-1} z[10+i]
 *      r8 = z10 - z[N+10] - z2
 *
 *  For each free-surface node m = 0..N,
 *
 *      r[9+m] = psi_m - z8 - z7 z[10+m]
 *
 *  and
 *
 *      r[N+10+m] = 0.5*((-z7 + u_m)^2 + v_m^2) + z[10+m] - z9.
 *
 *  These equations solve the wave form, wave speed, current variables, mean-flow
 *  flux, Bernoulli offset and Fourier coefficients as one coupled nonlinear
 *  free-boundary problem.  They are not empirical fitting equations and are not a
 *  post-processing correction to linear Airy theory.
 *
 *  6. CURRENT CONVENTIONS
 *  -----------------------------------------------------------------------------
 *  The solver keeps the current variables separate, as required by the Fenton formulation:
 *
 *      ubar_1 = c - Ubar,
 *      ubar_2 = c - Q/d,
 *      q      = Ubar d - Q.
 *
 *  The GUI uses the Eulerian-current criterion by default, but the state vector and
 *  residual equation r6 preserve the Fenton current-selector convention.
 *
 *  7. NUMERICAL SOLUTION STRATEGY
 *  -----------------------------------------------------------------------------
 *  Wave height is introduced by continuation from a near-linear wave to the target
 *  value.  At each continuation step, Newton iteration linearizes F(z) about the
 *  current iterate and solves J dz = -F(z) using analytic derivatives and a
 *  column-equilibrated one-sided Jacobi SVD with relative cutoff 1e-15.
 *  Failed height steps are retried with smaller increments. The Fourier order
 *  is reduced when the equilibrated Jacobian condition exceeds 1e11.
 *  The report states the order used; residual tolerances are independent of order.
 *
 *  The linear Airy/Fenton-McKee estimates are used only as initialization and
 *  checking aids.  The final result is the converged nonlinear stream-function /
 *  Fourier collocation solution.
 *
 *  8. PHYSICAL AND NUMERICAL ADMISSIBILITY
 *  -----------------------------------------------------------------------------
 *  H is the height of one deterministic regular wave. Numerical acceptance
 *  requires finite values, a small surface correction and max_i |r_i| <= 1e-10
 *  at the final height. Intermediate steps use a residual tolerance of 1e-8.
 *  The final residual is evaluated again after dimensional post-processing.
 *  Accepted outputs also require interpolated-surface residual <= 1e-6 and
 *  spectral/integral agreement within 1e-7 + 1e-8*abs(reference).
 *  A separate interval evaluator checks all finite equations and nine outputs.
 *  IEEE binary64 enclosures use bounded Taylor series and outward rounding.
 *  Phase extrema and quadrature refinement are numerical estimates.
 *  These checks do not constitute a formal proof of all floating-point executions.
 *  Miche and Fenton limiting-height ratios and the spectral tail are reported
 *  as diagnostics; they are not a proof of the physical wave's stability.
 *
 *  9. OUTPUT INTERPRETATION
 *  -----------------------------------------------------------------------------
 *  The report follows the distinction between k-based and d-based
 *  nondimensional quantities: kd, kH, T sqrt(gk), c sqrt(k/g), current variables,
 *  fluxes, Bernoulli constants, impulse, energy, momentum flux, radiation stress
 *  and wave power are reported alongside their engineering depth-scaled forms.
 *
 *  The same calculation is run for no-current and, when requested, with-current
 *  cases so the output can be compared under the selected current convention.
 *
 *  10. BUILDING FROM SOURCE
 *  -----------------------------------------------------------------------------
 *  Windows / MinGW-w64 release-like GUI build:
 *
 *      g++ fenton_gui.cpp -o fenton_gui.exe -O3 -std=c++20 ^
 *          -march=native -flto=auto -fopenmp -static -static-libgcc ^
 *          -static-libstdc++ -mwindows -pthread ^
 *          -lgdi32 -luser32 -lkernel32 -lcomctl32
 *
 *  Keep the same compiler family, CPU family and optimization flags when bitwise
 *  reproducibility is required, because Newton acceptance and SVD
 *  truncation can be sensitive to floating-point differences.
 * ==============================================================================
 */

#define _USE_MATH_DEFINES

#ifndef WIN32_LEAN_AND_MEAN
#define WIN32_LEAN_AND_MEAN
#endif

#ifndef NOMINMAX
#define NOMINMAX 1
#endif

#ifndef FENTON_HEADLESS
#include <windows.h>
#include <commctrl.h>
#endif

// --- nanosleep64 stub (MSYS2/MinGW GCC 15 + -static workaround) -----------------
#if defined(_WIN32) && defined(__MINGW32__) && !defined(FENTON_HEADLESS)

#include <time.h>    // struct _timespec64
#include <errno.h>
#include <limits.h>  // LLONG_MAX

extern "C" __attribute__((used))
int __cdecl nanosleep64(const struct _timespec64* request, struct _timespec64* remain)
{
    if (!request) { errno = EINVAL; return -1; }

    const long long sec  = (long long)request->tv_sec;
    const long      nsec = (long)request->tv_nsec;

    if (sec < 0 || nsec < 0 || nsec >= 1000000000L) { errno = EINVAL; return -1; }

    // Convert to 100ns units for WaitableTimer (ceil so we don't undersleep)
    const long long max_sec = LLONG_MAX / 10000000LL;
    const long long sec_c   = (sec > max_sec) ? max_sec : sec;

    long long total_100ns = sec_c * 10000000LL + ((long long)nsec + 99LL) / 100LL;
    if (total_100ns <= 0) total_100ns = 1;

    HANDLE t = CreateWaitableTimerW(NULL, TRUE, NULL);
    if (!t) {
        // Fallback: Sleep (ms resolution)
        unsigned long long total_ms =
            (unsigned long long)sec_c * 1000ULL +
            (unsigned long long)((nsec + 999999L) / 1000000L);

        DWORD ms = (total_ms > 0xFFFFFFFEULL) ? 0xFFFFFFFEUL : (DWORD)total_ms;
        Sleep(ms);
    } else {
        LARGE_INTEGER due;
        due.QuadPart = -total_100ns; // relative
        if (SetWaitableTimer(t, &due, 0, NULL, NULL, FALSE))
            WaitForSingleObject(t, INFINITE);
        CloseHandle(t);
    }

    if (remain) {
        remain->tv_sec  = 0;
        remain->tv_nsec = 0;
    }
    return 0;
}

#endif // defined(_WIN32) && defined(__MINGW32__)
// -------------------------------------------------------------------------------

// Standard C++ support for numerical solution, reporting and Win32 GUI plumbing
#include <algorithm>
#include <atomic>
#include <cmath>
#include <cfenv>
#include <cstdint>
#include <iomanip>
#include <limits>
#include <memory>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>
#include <fstream>

#ifdef _OPENMP
#include <omp.h>
#endif

// ==============================================================================
//  OPENMP EXECUTION CONTROL
// ==============================================================================

[[maybe_unused]] static int openmp_thread_count() noexcept {
#ifdef _OPENMP
    const int available = std::max(1, omp_get_num_procs());
    const int requested = std::max(1, omp_get_max_threads());
    return std::max(1, std::min(available, requested));
#else
    return 1;
#endif
}


// ==============================================================================
//  PHYSICAL CONSTANTS (dimensional convention)
// ==============================================================================

namespace Phys {
    using Real = double;
    constexpr Real PI    = 3.141592653589793238462643383279502884;
    constexpr Real RHO   = 1025.0;
    constexpr Real G_STD = 9.80665;
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

// ==============================================================================
//  SMALL LINEAR ALGEBRA (dense stabilized solves for num = 2N + 10)
// ==============================================================================

namespace LinAlg {

// Column-equilibrated Jacobi SVD: W = (A/D) V has column norms sigma_j.
static std::vector<Real> svd_solve(const std::vector<Real>& matrix,
                                  const std::vector<Real>& rhs, int n, Real& condition) {
    std::vector<Real> W((size_t)n*n), V((size_t)n*n, 0.0), scale((size_t)n, 1.0);
    for (int j = 0; j < n; ++j) {
        Real maximum = 0.0;
        for (int i = 0; i < n; ++i) maximum = std::max(maximum, std::abs(matrix[(size_t)i*n+j]));
        if (maximum > 0.0) scale[(size_t)j] = maximum;
    }
    for (int j = 0; j < n; ++j) {
        V[(size_t)j*n+j] = 1.0;
        for (int i = 0; i < n; ++i) {
            const Real value = matrix[(size_t)i*n+j];
            if (!std::isfinite(value)) throw std::runtime_error("Non-finite Jacobian.");
            W[(size_t)j*n+i] = value / scale[(size_t)j];
        }
        if (!std::isfinite(rhs[(size_t)j])) throw std::runtime_error("Non-finite linear right-hand side.");
    }
    const Real tolerance = 32.0 * std::numeric_limits<Real>::epsilon();
    bool orthogonal = false;
    for (int sweep = 0; sweep < 100; ++sweep) {
        orthogonal = true;
        for (int p = 0; p < n; ++p) for (int q = p+1; q < n; ++q) {
            Real pp = 0.0, qq = 0.0, pq = 0.0;
            for (int i = 0; i < n; ++i) {
                const Real a = W[(size_t)p*n+i], b = W[(size_t)q*n+i];
                pp += a*a; qq += b*b; pq += a*b;
            }
            if (!(pp > 0.0) || !(qq > 0.0) ||
                std::abs(pq) <= tolerance * std::sqrt(pp) * std::sqrt(qq)) continue;
            orthogonal = false;
            const Real tau = (qq-pp)/(2.0*pq);
            const Real t = std::copysign(1.0,tau)/(std::abs(tau)+std::hypot(1.0,tau));
            const Real c = 1.0/std::hypot(1.0,t), sn = t*c;
            for (int i = 0; i < n; ++i) {
                const Real a = W[(size_t)p*n+i], b = W[(size_t)q*n+i];
                W[(size_t)p*n+i] = c*a-sn*b;
                W[(size_t)q*n+i] = sn*a+c*b;
                const Real vp = V[(size_t)p*n+i], vq = V[(size_t)q*n+i];
                V[(size_t)p*n+i] = c*vp-sn*vq;
                V[(size_t)q*n+i] = sn*vp+c*vq;
            }
        }
        if (orthogonal) break;
    }
    if (!orthogonal) throw std::runtime_error("Jacobi SVD did not converge.");
    std::vector<Real> sigma((size_t)n, 0.0), x((size_t)n, 0.0);
    Real largest = 0.0;
    for (int j = 0; j < n; ++j) {
        for (int i = 0; i < n; ++i) sigma[(size_t)j] = std::hypot(sigma[(size_t)j], W[(size_t)j*n+i]);
        largest = std::max(largest, sigma[(size_t)j]);
    }
    if (!std::isfinite(largest)) throw std::runtime_error("Non-finite singular value.");
    const Real smallest = *std::min_element(sigma.begin(), sigma.end());
    condition = smallest > largest*1e-300 ? largest/smallest : std::numeric_limits<Real>::infinity();
    for (int j = 0; j < n; ++j) {
        const Real sj = sigma[(size_t)j];
        if (!(sj > largest*1e-15)) continue;
        Real projection = 0.0;
        for (int i = 0; i < n; ++i) projection += (W[(size_t)j*n+i]/sj)*rhs[(size_t)i];
        for (int i = 0; i < n; ++i) x[(size_t)i] += V[(size_t)j*n+i]*(projection/sj);
    }
    for (int j = 0; j < n; ++j) x[(size_t)j] /= scale[(size_t)j];
    return x;
}

} // namespace LinAlg

// ==============================================================================
//  FENTON STREAM-FUNCTION SOLVER (z[i] formulation)
// ==============================================================================


class FentonStreamFunction {
public:
    // ----------------------- physical input variables -------------------------
    Real H_target;   // [m] wave height H, crest-to-trough
    Real T_target;   // [s] apparent fixed-frame period T
    Real d;          // [m] still-water depth d
    Real Uc;         // [m/s] imposed collinear current U_c, GUI Eulerian convention
    Real g;          // [m/s^2] gravitational acceleration g
    int  N;          // Fourier order / half-wave interval count N (requested N=50)

    // ----------------------- continuation / Newton controls --------------------
    int  nstep;          // continuation steps in H/(gT^2) toward the target wave height
    int  number;         // maximum Newton iterations per continuation step
    Real crit;           // intermediate-step convergence threshold on surface correction
    Real criter_final;   // final-step convergence threshold on surface correction

    // Current selector used in residual r6: 1 = Eulerian current, 2 = Stokes/mass-transport current.
    int Current_criterion;

    // Derived nondimensional inputs
    Real MaxH;     // H/d, prescribed relative wave height
    Real T_nd;     // T sqrt(g/d), depth-based dimensionless period
    Real Height;   // H/(gT^2), period-form height parameter H_s
    Real Current;  // U_c/sqrt(gd), depth-based imposed-current Froude variable

    // ----------------------- dimensional and report outputs --------------------
    Real k;        // [rad/m] wave number, k = 2*pi/L
    Real L;        // [m] dynamically consistent wavelength L
    Real c;        // [m/s] fixed-frame wave celerity, c = L/T

    std::vector<Real> eta_nodes; // size N+1, free-surface node elevations as bed-origin coordinates [m]
    std::vector<Real> Bj;        // size N, Fourier stream-function coefficients B_1..B_N

    Real eta_crest;    // [m] relative to SWL
    Real eta_trough;   // [m] relative to SWL

    Real steepness;
    Real rel_depth;
    Real ursell;
    std::string regime;

    // Highest-wave and numerical-admissibility diagnostics.
    Real breaking_limit_miche = 0.0;
    Real breaking_index = 0.0;
    bool is_breaking = false;
    Real highest_wave_limit_fenton;     // [m] theoretical highest steady wave from Fenton (1990), Eq. (32)
    Real highest_wave_ratio;            // H/Hm
    Real fourier_accuracy_limit;         // [m] 98% of Hm, limiting-height diagnostic
    Real surface_tail_coefficient;       // |E_N|, spectral-resolution indicator
    Real max_residual;                   // ||F(z)||_infinity after final convergence
    Real crest_stagnation_margin;        // [m/s] c-u at the crest; must remain positive below the limiting wave
    bool near_limiting_wave;
    bool practical_breaking_risk;

    // Integral properties reported in the SOLUTION.RES convention
    Real EulerianCurrent;
    Real StokesCurrent;
    Real MeanFluidSpeed;

    Real WaveVolumeFlux_q;
    Real VolumeFluxQ;

    Real BernoulliR;
    Real Bernoulli_r;

    Real KineticEnergy;
    Real PotentialEnergy;
    Real EnergyDensity;

    Real MomentumFlux;
    Real MomentumFluxDepth;

    Real Sxx;
    Real Sxx_depth;

    Real Impulse;
    Real I_depth;

    Real Power;
    Real F_depth;

    Real Cg;

    // Depth-scaled invariants used in the dimensionless-output table
    Real E_depth;
    Real KE_depth;
    Real PE_depth;

    // Stokes / mass-transport current derived from the nonlinear flux solution
    Real MassTransport;

    // Bed and free-surface kinematic diagnostics derived after convergence
    Real MeanSquareBedAbsoluteVelocity = 0.0;
    Real MeanSquareBedVelocity;  // ub^2 [m^2/s^2], phase-mean orbital bed velocity square
    Real u_bed;
    Real u_surf;
    Real u_bed_crest = 0.0;
    Real u_surf_crest = 0.0;
    std::map<std::string, Real> extrema_phases;
    Real acc_max;
    Real acc_bed_max;
    Real w_max;
    Real asymmetry;

    // Convenience quantities for secondary engineering diagnostics
    Real tau_bed;
    Real ExcursionBed;

    int requested_order;
    Real jacobian_condition = 0.0;
    Real continuation_fraction = 0.0;
    std::string numerical_note;
    bool converged;
    bool accuracy_verified = false;
    bool finite_system_verified = false;
    Verified::Certificate finite_system_certificate;
    std::map<std::string, Real> accuracy_metrics;
    std::string solution_status = "not_solved";
    std::string last_error;

    // -------------------------- construction --------------------------------
    explicit FentonStreamFunction(Real H, Real T, Real depth, Real current = 0.0, int fourier_order = 50)
        : H_target(H), T_target(T), d(depth), Uc(current),
          g(Phys::G_STD), N(std::max(1, fourier_order)),
          nstep(4), number(40), crit(1e-8), criter_final(1e-10),
          Current_criterion(1),
          MaxH(0.0), T_nd(0.0), Height(0.0), Current(0.0),
          k(0.0), L(0.0), c(0.0),
          eta_nodes((size_t)N + 1, depth),
          Bj((size_t)N, 0.0),
          eta_crest(0.0), eta_trough(0.0),
          steepness(0.0), rel_depth(0.0), ursell(0.0), regime(""),
          highest_wave_limit_fenton(0.0), highest_wave_ratio(0.0), fourier_accuracy_limit(0.0),
          surface_tail_coefficient(0.0), max_residual(0.0), crest_stagnation_margin(0.0),
          near_limiting_wave(false), practical_breaking_risk(false),
          EulerianCurrent(0.0), StokesCurrent(0.0), MeanFluidSpeed(0.0),
          WaveVolumeFlux_q(0.0), VolumeFluxQ(0.0),
          BernoulliR(0.0), Bernoulli_r(0.0),
          KineticEnergy(0.0), PotentialEnergy(0.0), EnergyDensity(0.0),
          MomentumFlux(0.0), MomentumFluxDepth(0.0),
          Sxx(0.0), Sxx_depth(0.0),
          Impulse(0.0), I_depth(0.0),
          Power(0.0), F_depth(0.0),
          Cg(0.0),
          MeanSquareBedVelocity(0.0),
          u_bed(0.0), u_surf(0.0), acc_max(0.0), acc_bed_max(0.0), w_max(0.0), asymmetry(0.0),
          tau_bed(0.0), ExcursionBed(0.0),
          requested_order(N), converged(false), last_error("")
    {
        // Input screening for the physical problem: H, T and d must be positive.
        if (d > 0.0) {
            MaxH   = H_target / d;
            T_nd   = T_target * std::sqrt(g / d);
            Height = (T_nd > 0.0) ? (MaxH / (T_nd * T_nd)) : 0.0;
            Current = Uc / std::sqrt(g * d);
        }

        // Large imposed-current cases may require more continuation/Newton work.
        if (std::abs(Current) >= 1.0) {
            nstep  = std::max(nstep, 8);
            number = std::max(number, 80);
        }

        // Internal arrays follow the Fenton 1-based z[1..num] state-vector layout.
        n = N;
        num = 2 * n + 10;

        z.assign((size_t)num + 1, 0.0);
        rhs1.assign((size_t)num + 1, 0.0);
        rhs2.assign((size_t)num + 1, 0.0);

        Tanh.assign((size_t)n + 1, 0.0);
        B.assign((size_t)n + 1, 0.0);
        Y.assign((size_t)num + 1, 0.0);

        cosa.assign((size_t)(2 * n + 1), 0.0);
        sina.assign((size_t)(2 * n + 1), 0.0);

        cos_nm.assign((size_t)(n + 1) * (size_t)n, 0.0);
        sin_nm.assign((size_t)(n + 1) * (size_t)n, 0.0);

        sol.assign((size_t)(num + 1) * 3, 0.0);

        init_trig_tables();
    }

    // ------------------------------- API ------------------------------------
    void solve() {
        accuracy_verified = false;
        finite_system_verified = false;
        finite_system_certificate = Verified::Certificate{};
        accuracy_metrics.clear();
        solution_status = "not_solved";
        const int requested = requested_order;
        int order = N;
        while (true) {
            solve_at_order();
            if (converged) {
                try {
                    LinAlg::svd_solve(jacobian(), std::vector<Real>((size_t)num, 0.0), num, jacobian_condition);
                } catch (const std::exception& e) {
                    converged = false;
                    last_error = e.what();
                    jacobian_condition = std::numeric_limits<Real>::infinity();
                }
            }
            if (jacobian_condition <= 1e11 || order <= 12) break;
            order = std::max(12, (int)(order * 0.8));
            FentonStreamFunction next(H_target, T_target, d, Uc, order);
            next.nstep = nstep; next.number = number;
            next.crit = crit; next.criter_final = criter_final;
            next.Current_criterion = Current_criterion;
            next.requested_order = requested;
            *this = std::move(next);
        }
        if (converged && jacobian_condition > 1e11) {
            converged = false;
            last_error = "The Fourier system is too ill-conditioned for a reliable binary64 solution.";
        }
        if (order != requested) {
            numerical_note = "Fourier order used: " + std::to_string(order) +
                " (requested " + std::to_string(requested) + "); limited by numerical conditioning.";
        }
        if (converged) {
            try {
                verify_accuracy();
            } catch (const std::exception& e) {
                converged = false;
                accuracy_verified = false;
                solution_status = "accuracy_failure";
                last_error = std::string("Accuracy verification failed: ") + e.what();
            }
        } else {
            solution_status = "nonconvergence";
        }
    }

    void solve_at_order() {
        converged = false;
        last_error.clear();

        if (!std::isfinite(H_target) || !std::isfinite(T_target) ||
            !std::isfinite(d) || !std::isfinite(Uc) ||
            !(H_target > 0.0) || !(T_target > 0.0) || !(d > 0.0)) {
            last_error = "Invalid inputs: H, T, and d must be finite and > 0; current must be finite.";
            return;
        }

        max_residual = std::numeric_limits<Real>::infinity();
        if (Current_criterion != 1 && Current_criterion != 2) {
            last_error = "Current criterion must be 1 (Eulerian) or 2 (Stokes).";
            return;
        }

        try {
            solve_internal();
            converged = true;
        } catch (const std::exception& e) {
            converged = false;
            last_error = std::string("Solver error: ") + e.what();
        } catch (...) {
            converged = false;
            last_error = "Solver error: unknown exception.";
        }
    }

    // GUI-facing kinematics: fixed-frame u, w and a_x at bed-origin z [m] and phase X=kx [rad].
    void get_kinematics(Real z_bed, Real phase, Real& u_abs, Real& w_abs, Real& ax) const {
        const Real kd = z[1];
        if (!(kd > 0.0) || !(d > 0.0)) {
            u_abs = w_abs = ax = 0.0;
            return;
        }

        const Real k_phys = kd / d;
        const Real X = phase;
        const Real Yloc = k_phys * (z_bed - d); // Y = k y, with bed-origin z converted to y = z-d

        Real u_nd = 0.0, v_nd = 0.0, dudt_nd = 0.0;
        point(X, Yloc, u_nd, v_nd, dudt_nd);

        const Real scale_v = std::sqrt(g * d);
        u_abs = u_nd * scale_v;
        w_abs = v_nd * scale_v;
        ax    = dudt_nd * g;
    }

    // Accessor used for SOLUTION.RES-style reporting: dimensionless wavenumber kd.
    Real kd_dimless() const { return (z.size() > 1) ? z[1] : 0.0; }

    // Numerical comparisons are acceptance checks, not formal floating-point proofs.
    void verify_accuracy() {
        accuracy_verified = false;
        finite_system_verified = false;
        finite_system_certificate = Verified::Certificate{};
        accuracy_metrics.clear();
        std::vector<std::string> failures;
        const int phase_count = std::max(2048, 32*n);
        Real boundary = 0.0;
        for (int i = 0; i < phase_count; ++i) {
            const Real phase = 2.0*Phys::PI*i/phase_count;
            const Real eta = surface_keta(phase);
            Real psi = 0.0, u = 0.0, v = 0.0;
            for (int j = 1; j <= n; ++j) {
                const Real a = j*eta, b = j*z[1];
                const Real denominator = 1.0+std::exp(-2.0*b);
                const Real positive = std::exp(a), reflected = std::exp(-a-2.0*b);
                const Real S = (positive-reflected)/denominator;
                const Real C = (positive+reflected)/denominator;
                psi += S*std::cos(j*phase)*B[(size_t)j];
                u += C*std::cos(j*phase)*j*B[(size_t)j];
                v += S*std::sin(j*phase)*j*B[(size_t)j];
            }
            const Real streamline = psi-z[8]-z[7]*eta;
            const Real bernoulli = 0.5*((u-z[7])*(u-z[7])+v*v)+eta-z[9];
            if (!std::isfinite(streamline) || !std::isfinite(bernoulli))
                throw std::runtime_error("Non-finite interpolated-surface residual.");
            boundary = std::max(boundary, std::max(std::abs(streamline), std::abs(bernoulli)));
        }
        accuracy_metrics["surface_residual"] = boundary;
        if (boundary > 1e-6) failures.push_back("interpolated-surface residual exceeds 1e-6");

        const Real beta = MeanSquareBedAbsoluteVelocity;
        const std::map<std::string, Real> spectral = {
            {"Power", c*(3*KineticEnergy-2*PotentialEnergy)
                +0.5*beta*(Impulse+Phys::RHO*c*d)-2*c*EulerianCurrent*Impulse},
            {"Sxx", 4*KineticEnergy-3*PotentialEnergy+Phys::RHO*beta*d-2*EulerianCurrent*Impulse}
        };
        const std::map<std::string, Real> values = {
            {"Power", Power}, {"Sxx", Sxx}, {"Impulse", Impulse},
            {"KineticEnergy", KineticEnergy}, {"PotentialEnergy", PotentialEnergy}
        };
        for (const auto& entry : spectral) {
            const Real error = std::abs(values.at(entry.first)-entry.second);
            const Real budget = 1e-7+1e-8*std::abs(entry.second);
            accuracy_metrics[entry.first+"_spectral_ratio"] = error/budget;
            if (!std::isfinite(error) || error > budget)
                failures.push_back(entry.first+" spectral-bed comparison exceeds its budget");
        }
        const auto coarse = defining_integrals(std::max(256, 8*n), 64);
        const auto fine = defining_integrals(std::max(512, 16*n), 96);
        for (const auto& entry : fine) {
            const Real budget = 1e-7+1e-8*std::abs(entry.second);
            const Real difference = std::abs(entry.second-coarse.at(entry.first));
            const Real error = std::abs(values.at(entry.first)-entry.second);
            accuracy_metrics[entry.first+"_quadrature_ratio"] = difference/budget;
            accuracy_metrics[entry.first+"_integral_ratio"] = error/budget;
            if (!std::isfinite(difference) || !std::isfinite(error) || std::max(error,difference) > budget)
                failures.push_back(entry.first+" defining-integral comparison exceeds its budget");
        }
        if (!verify_finite_system()) failures.push_back("finite-system interval verification exceeds its budget");
        accuracy_metrics["finite_residual_ratio"] = finite_system_certificate.residual_upper/finite_system_certificate.residual_tolerance;
        Real output_ratio=0.0;
        for (const auto& entry : finite_system_certificate.outputs)
            output_ratio=std::max(output_ratio,entry.second.error_upper/entry.second.budget_lower);
        accuracy_metrics["finite_output_ratio"]=output_ratio;
        if (!failures.empty()) {
            converged = false;
            solution_status = "resolution_failure";
            last_error = "The collocation system converged, but the available Fourier resolution "
                         "does not satisfy the output accuracy budgets: ";
            for (size_t i = 0; i < failures.size(); ++i) {
                if (i) last_error += "; ";
                last_error += failures[i];
            }
            last_error += ".";
            return;
        }
        accuracy_verified = true;
        solution_status = "verified_numerically";
    }

    bool verify_finite_system() {
        using I=Verified::Interval;
        finite_system_verified=false;
        finite_system_certificate=Verified::Certificate{};
        volatile Real smallest=std::numeric_limits<Real>::denorm_min();
        if (!std::numeric_limits<Real>::is_iec559 || std::fegetround()!=FE_TONEAREST ||
            smallest*2.0!=std::nextafter(std::numeric_limits<Real>::denorm_min(),1.0))
            throw std::runtime_error("Interval verification requires IEEE binary64, round-to-nearest and gradual underflow.");
        std::vector<I> state;state.reserve(z.size());
        for (Real value : z) state.emplace_back(value);
        const I H(H_target),T(T_target),depth(d),current(Uc),gravity(g),density(Phys::RHO);
        std::vector<I> F((size_t)(2*n+10));
        F[0]=state[2]-state[1]*H/depth;
        F[1]=state[2]-H/(gravity*T.square())*state[3].square();
        F[2]=state[4]*state[3]-I(2)*Verified::pi();
        F[3]=state[5]+state[7]-state[4];
        F[4]=state[1]*(state[6]+state[7]-state[4])-state[8];
        F[5]=state[(size_t)Current_criterion+4]-current/(gravity*depth).sqrt()*state[1].sqrt();
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
                const I S=(ep-em)/denominators[(size_t)j-1],C=(ep+em)/denominators[(size_t)j-1];
                const auto& trig=table[(size_t)((m*j)%(2*n))];
                const I coefficient=state[(size_t)n+10+j];
                psi=psi+coefficient*S*trig.second;
                u=u+I(j)*coefficient*C*trig.second;
                v=v+I(j)*coefficient*S*trig.first;
            }
            F[(size_t)8+m]=psi-state[8]-state[7]*ordinate;
            F[(size_t)n+9+m]=((u-state[7]).square()+v.square())/I(2)+ordinate-state[9];
        }
        const I wave_number=state[1]/depth,speed_scale=(gravity/wave_number).sqrt();
        const I energy_scale=density*gravity/wave_number.square();
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
            {"L",L},{"c",c},{"Impulse",Impulse},{"KineticEnergy",KineticEnergy},
            {"PotentialEnergy",PotentialEnergy},{"EnergyDensity",EnergyDensity},
            {"Sxx",Sxx},{"Power",Power},{"MeanSquareBedVelocity",MeanSquareBedVelocity}
        };
        auto& certificate=finite_system_certificate;
        certificate.state=z;
        certificate.surface_coefficients.assign(Y.begin(),Y.begin()+n+1);
        certificate.order=n;certificate.equation_count=(int)F.size();certificate.residual_intervals=F;
        // The lower adjacent constants are conservative decimal tolerance bounds.
        certificate.residual_tolerance=I::down(1e-10);
        for (const auto& residual : F) certificate.residual_upper=std::max(certificate.residual_upper,residual.magnitude());
        bool valid=certificate.residual_upper<=certificate.residual_tolerance;
        for (const auto& entry : outputs) {
            const auto& interval=entry.second;
            const Real minimum=interval.lo<=0.0 && interval.hi>=0.0 ? 0.0 : std::min(std::abs(interval.lo),std::abs(interval.hi));
            Verified::OutputBound check;
            check.enclosure=interval;check.reported=reported.at(entry.first);
            check.error_upper=(I(check.reported)-interval).magnitude();
            check.budget_lower=I::down(I::down(1e-7)+I::down(I::down(1e-8)*minimum));
            valid=valid && check.error_upper<=check.budget_lower;
            certificate.outputs.emplace(entry.first,check);
        }
        certificate.passed=valid;finite_system_verified=valid;
        return valid;
    }

private:
    Real phase_quantity(Real phase, const std::string& quantity) const {
        Real u, w, acceleration;
        const bool bed = quantity == "bed" || quantity == "bed_orbital" || quantity == "bed_acceleration";
        const Real top = bed ? 0.0 : d*(1.0+surface_keta(phase)/z[1]);
        get_kinematics(top, phase, u, w, acceleration);
        const Real value = quantity == "bed_orbital" ? u-EulerianCurrent
                         : quantity == "vertical" ? w
                         : (quantity == "acceleration" || quantity == "bed_acceleration") ? acceleration : u;
        if (!std::isfinite(value)) throw std::runtime_error("Non-finite phase evaluation.");
        return std::abs(value);
    }

    std::pair<Real, Real> phase_maximum(const std::string& quantity) const {
        const int count = std::max(257, 8*n+1);
        std::vector<Real> phases((size_t)count), values((size_t)count);
        Real best = -1.0, phase = 0.0;
        for (int i = 0; i < count; ++i) {
            phases[(size_t)i] = Phys::PI*i/(count-1);
            values[(size_t)i] = phase_quantity(phases[(size_t)i], quantity);
            if (values[(size_t)i] > best) { best = values[(size_t)i]; phase = phases[(size_t)i]; }
        }
        const Real golden = (std::sqrt(5.0)-1.0)/2.0;
        for (int i = 1; i < count-1; ++i) {
            if (!(values[(size_t)i] >= values[(size_t)i-1] && values[(size_t)i] >= values[(size_t)i+1]
                && (values[(size_t)i] > values[(size_t)i-1] || values[(size_t)i] > values[(size_t)i+1]))) continue;
            Real left = phases[(size_t)i-1], right = phases[(size_t)i+1];
            Real x1 = right-golden*(right-left), x2 = left+golden*(right-left);
            Real f1 = phase_quantity(x1, quantity), f2 = phase_quantity(x2, quantity);
            for (int iteration = 0; iteration < 64; ++iteration) {
                if (right-left < 1e-12) break;
                if (f1 < f2) {
                    left = x1; x1 = x2; f1 = f2;
                    x2 = left+golden*(right-left); f2 = phase_quantity(x2, quantity);
                } else {
                    right = x2; x2 = x1; f2 = f1;
                    x1 = right-golden*(right-left); f1 = phase_quantity(x1, quantity);
                }
                if (f1 > best) { best = f1; phase = x1; }
                if (f2 > best) { best = f2; phase = x2; }
            }
        }
        return {best, phase};
    }

    void compute_extrema() {
        extrema_phases.clear();
        for (const auto& entry : std::vector<std::pair<std::string, Real*>>{
                {"bed", &u_bed}, {"surface", &u_surf}, {"vertical", &w_max}, {"acceleration", &acc_max}, {"bed_acceleration", &acc_bed_max}}) {
            const auto peak = phase_maximum(entry.first);
            *entry.second = peak.first;
            extrema_phases[entry.first] = peak.second;
        }
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

    std::map<std::string, Real> defining_integrals(int phase_count, int vertical_count) const {
        const auto quadrature = gauss_legendre(vertical_count);
        Real impulse = 0.0, kinetic = 0.0, potential = 0.0, stress = 0.0, power = 0.0;
        std::vector<Real> cosine((size_t)n+1), sine((size_t)n+1);
        for (int i = 0; i < phase_count; ++i) {
            const Real phase = 2.0*Phys::PI*i/phase_count;
            const Real eta = surface_keta(phase)/k, thickness = (d+eta)/2.0;
            if (!(thickness > 0.0) || !std::isfinite(thickness))
                throw std::runtime_error("The interpolated surface intersects the bed or is non-finite.");
            potential += 0.5*Phys::RHO*g*eta*eta/phase_count;
            for (int j = 1; j <= n; ++j) {
                cosine[(size_t)j] = std::cos(j*phase);
                sine[(size_t)j] = std::sin(j*phase);
            }
            for (int m = 0; m < vertical_count; ++m) {
                const Real y = -d+thickness*(quadrature.first[(size_t)m]+1.0);
                Real u = EulerianCurrent, v = 0.0;
                for (int j = 1; j <= n; ++j) {
                    const Real a = j*k*y, b = j*z[1];
                    const Real denominator = 1.0+std::exp(-2*b);
                    const Real positive = std::exp(a), reflected = std::exp(-a-2*b);
                    const Real amplitude = j*B[(size_t)j]*std::sqrt(g/k);
                    u += (positive+reflected)/denominator*amplitude*cosine[(size_t)j];
                    v += (positive-reflected)/denominator*amplitude*sine[(size_t)j];
                }
                const Real pressure = Phys::RHO*(Bernoulli_r-g*y-0.5*((u-c)*(u-c)+v*v));
                const Real weight = quadrature.second[(size_t)m]*thickness/phase_count;
                impulse += weight*Phys::RHO*u;
                kinetic += weight*0.5*Phys::RHO*(u*u+v*v);
                stress += weight*(pressure+Phys::RHO*u*u);
                power += weight*(pressure+Phys::RHO*g*y+0.5*Phys::RHO*(u*u+v*v))*u;
            }
        }
        return {{"Impulse", impulse}, {"KineticEnergy", kinetic}, {"PotentialEnergy", potential},
                {"Sxx", stress-0.5*Phys::RHO*g*d*d}, {"Power", power}};
    }

    int n = 0;
    int num = 0;

    // Fenton 1-based vectors; index 0 is intentionally unused.
    std::vector<Real> z, rhs1, rhs2, Tanh, B, Y;

    // Collocation trigonometric tables for X_m = m*pi/N and harmonics j = 1..N.
    std::vector<Real> cosa, sina;
    std::vector<Real> cos_nm, sin_nm; // (N+1) x N, row-major values cos(jX_m), sin(jX_m)

    // Continuation storage for extrapolating the state vector between height steps.
    std::vector<Real> sol;

    // Current continuation-step variables
    Real height = 0.0; // stepped H_s = H/(gT^2) value used in r2
    Real Hoverd = 0.0; // stepped H/d value used in r1

private:
    static constexpr Real kFourierApplicabilityFraction = 0.98;
    static constexpr Real kFinalResidualTolerance = 1.0e-10;
    static constexpr Real kIntermediateResidualTolerance = 1.0e-8;

    inline Real& SOL(int i, int k) { return sol[(size_t)i * 3 + (size_t)k]; }
    inline Real  SOL(int i, int k) const { return sol[(size_t)i * 3 + (size_t)k]; }

    static Real fenton_highest_hoverd(Real lambda_over_d) {
        if (!(lambda_over_d > 0.0) || !std::isfinite(lambda_over_d)) {
            return std::numeric_limits<Real>::quiet_NaN();
        }

        // Fenton (1990), Eq. (32), evaluated in a form that avoids overflow for
        // very large wavelength/depth ratios.
        const Real x = lambda_over_d;
        if (x <= 1.0) {
            const Real numerator = x * (0.141063 + x * (0.0095721 + 0.0077829 * x));
            const Real denominator = 1.0 + x * (0.0788340 + x * (0.0317567 + 0.0093407 * x));
            return numerator / denominator;
        }

        const Real inv = 1.0 / x;
        const Real inv2 = inv * inv;
        const Real inv3 = inv2 * inv;
        const Real numerator = 0.0077829 + 0.0095721 * inv + 0.141063 * inv2;
        const Real denominator = 0.0093407 + 0.0317567 * inv + 0.0788340 * inv2 + inv3;
        return numerator / denominator;
    }

    Real residual_inf_norm(const std::vector<Real>& residual) const {
        Real value = 0.0;
        for (int i = 1; i <= num; ++i) value = std::max(value, std::abs(residual[(size_t)i]));
        return value;
    }

    void init_trig_tables() {
        // cosa[k] = cos(k*pi/N), k = 0..2N, used by the half-wave collocation grid.
        for (int i = 0; i <= 2 * n; ++i) {
            const Real ang = (Real)i * Phys::PI / (Real)n;
            cosa[(size_t)i] = std::cos(ang);
            sina[(size_t)i] = std::sin(ang);
        }

        // cos_nm[m,j] and sin_nm[m,j] store cos(jX_m), sin(jX_m), X_m = m*pi/N.
        for (int m = 0; m <= n; ++m) {
            for (int j = 1; j <= n; ++j) {
                const int nm = (m * j) % (2 * n);
                cos_nm[(size_t)m * (size_t)n + (size_t)(j - 1)] = cosa[(size_t)nm];
                sin_nm[(size_t)m * (size_t)n + (size_t)(j - 1)] = sina[(size_t)nm];
            }
        }
    }

    // ----------------------------------------------------------------------
    // Linear/Fenton-McKee initial estimate used only to seed the nonlinear residual solve.
    // ----------------------------------------------------------------------
    void init_linear() {
        const Real pi = Phys::PI;

        const Real sigma = (Hoverd > 0.0) ? (2.0 * pi * std::sqrt(height / Hoverd)) : 0.0;

        if (sigma > 0.0) {
            const Real t = std::tanh(std::pow(sigma, 1.5));
            // Fenton-McKee approximation for the initial kd estimate; not the final theory.
            z[1] = (sigma * sigma) / std::pow(t, (Real)(2.0 / 3.0));
        } else {
            z[1] = 2.0 * pi * std::max(height, (Real)1e-12) / std::max(Hoverd, (Real)1e-12);
        }

        z[2] = z[1] * Hoverd;
        z[4] = std::sqrt(std::tanh(z[1]));
        z[3] = 2.0 * pi / z[4];

        // Current-variable initialization for the selected Eulerian/Stokes closure.
        if (Current_criterion == 1) {
            z[5] = Current * std::sqrt(z[1]); // ubar_1 sqrt(k/g), Eulerian current variable
            z[6] = 0.0;                       // ubar_2 sqrt(k/g), Stokes current variable
        } else {
            z[6] = Current * std::sqrt(z[1]);
            z[5] = 0.0;
        }

        z[7] = z[4];       // Ubar sqrt(k/g), mean wave-frame velocity variable
        z[8] = 0.0;        // q sqrt(k^3/g), flux-related variable
        z[9] = 0.5 * z[7] * z[7];

        z[10] = 0.5 * z[2];
        for (int i = 1; i <= n; ++i) {
            z[n + i + 10] = 0.0;                 // Fourier coefficient B_i
            z[i + 10] = 0.5 * z[2] * cosa[(size_t)i];
        }
        z[n + 11] = 0.5 * z[2] / z[7];

        // Store the initial state for continuation extrapolation.
        for (int i = 1; i < 10; ++i) SOL(i, 1) = z[(size_t)i];
        for (int i = 10; i <= num; ++i) SOL(i, 1) = 0.0;
    }

    // ----------------------------------------------------------------------
    // residual vector F(z): fills rhs_out[1..2N+10] and returns sum(F_i^2).
    // ----------------------------------------------------------------------
    Real eqns_for_state(const std::vector<Real>& state,
                        std::vector<Real>& rhs_out,
                        std::vector<Real>& tanh_workspace) const {
        const Real pi = Phys::PI;
        if (rhs_out.size() != (size_t)num + 1) rhs_out.resize((size_t)num + 1);
        std::fill(rhs_out.begin(), rhs_out.end(), 0.0);
        if (tanh_workspace.size() != (size_t)n + 1) tanh_workspace.resize((size_t)n + 1);

        rhs_out[1] = state[2] - state[1] * Hoverd;
        rhs_out[2] = state[2] - height * state[3] * state[3];
        rhs_out[3] = state[4] * state[3] - 2.0 * pi;
        rhs_out[4] = state[5] + state[7] - state[4];
        rhs_out[5] = state[1] * (state[6] + state[7] - state[4]) - state[8];

        for (int i = 1; i <= n; ++i) {
            tanh_workspace[(size_t)i] = std::tanh((Real)i * state[1]);
        }

        rhs_out[6] = state[(size_t)(Current_criterion + 4)] - Current * std::sqrt(state[1]);

        rhs_out[7] = state[10] + state[(size_t)(n + 10)];
        for (int i = 1; i < n; ++i) rhs_out[7] += 2.0 * state[(size_t)(10 + i)];

        rhs_out[8] = state[10] - state[(size_t)(n + 10)] - state[2];

        for (int m = 0; m <= n; ++m) {
            const Real zsurf = state[(size_t)(10 + m)];

            Real psi = 0.0;
            Real u   = 0.0;
            Real v   = 0.0;

            const Real* cosrow = &cos_nm[(size_t)m * (size_t)n];
            const Real* sinrow = &sin_nm[(size_t)m * (size_t)n];

            for (int j = 1; j <= n; ++j) {
                const Real x = (Real)j * zsurf;
                if (x > 60.0 || x < -60.0) {
                    throw std::runtime_error("Divergence: exp(j*zsurf) out of safe range.");
                }

                const Real ep = std::exp(x);
                const Real em = std::exp(-x - 2.0*j*state[1]);
                const Real den = 1.0 + std::exp(-2.0*j*state[1]);
                const Real S = (ep-em)/den;
                const Real C = (ep+em)/den;

                const Real c_nm = cosrow[(size_t)(j - 1)];
                const Real s_nm = sinrow[(size_t)(j - 1)];

                const Real cj = state[(size_t)(n + j + 10)];
                const Real j_cj = (Real)j * cj;

                psi += cj * S * c_nm;
                u   += j_cj * C * c_nm;
                v   += j_cj * S * s_nm;
            }

            rhs_out[(size_t)(m + 9)] = psi - state[8] - state[7] * state[(size_t)(m + 10)];
            rhs_out[(size_t)(n + m + 10)] =
                0.5 * (((-state[7] + u) * (-state[7] + u)) + v * v)
                + state[(size_t)(m + 10)] - state[9];
        }

        Real ss = 0.0;
        for (int i = 1; i <= num; ++i) {
            const Real ri = rhs_out[(size_t)i];
            ss += ri * ri;
        }
        return ss;
    }

    Real eqns(std::vector<Real>& rhs_out) const {
        std::vector<Real> tanh_workspace((size_t)n + 1, 0.0);
        return eqns_for_state(z, rhs_out, tanh_workspace);
    }

    // ----------------------------------------------------------------------
    // SVD solve for the Newton correction with small singular-value truncation.
    // Solves J dz = -F for the square system, num = 2N + 10.
    // ----------------------------------------------------------------------
    std::vector<Real> svd_solve(const std::vector<Real>& A, const std::vector<Real>& b, int size) {
        return LinAlg::svd_solve(A, b, size, jacobian_condition);
    }

    std::vector<Real> jacobian() const {
        std::vector<Real> A((size_t)num*num, 0.0);
        auto at = [&](int r, int c) -> Real& { return A[(size_t)(r-1)*num+c-1]; };
        at(1,1)=-Hoverd; at(1,2)=1.0;
        at(2,2)=1.0; at(2,3)=-2.0*height*z[3];
        at(3,3)=z[4]; at(3,4)=z[3];
        at(4,5)=1.0; at(4,7)=1.0; at(4,4)=-1.0;
        at(5,1)=z[6]+z[7]-z[4];
        at(5,6)=z[1]; at(5,7)=z[1]; at(5,4)=-z[1]; at(5,8)=-1.0;
        at(6,Current_criterion+4)=1.0; at(6,1)=-0.5*Current/std::sqrt(z[1]);
        for (int m=0; m<=n; ++m) at(7,10+m)=(m==0 || m==n) ? 1.0 : 2.0;
        at(8,10)=1.0; at(8,n+10)=-1.0; at(8,2)=-1.0;
        for (int m=0; m<=n; ++m) {
            const int sr=9+m, dr=n+10+m;
            const Real surface=z[(size_t)10+m];
            Real u=0.0, v=0.0, duk=0.0, dvk=0.0, duy=0.0, dvy=0.0;
            std::vector<Real> S((size_t)n), C((size_t)n);
            for (int j=1; j<=n; ++j) {
                const size_t idx=(size_t)m*n+j-1;
                const Real co=cos_nm[idx], si=sin_nm[idx], bj=z[(size_t)n+10+j];
                const Real ep=std::exp(j*surface), q=std::exp(-2.0*j*z[1]);
                const Real em=std::exp(-j*surface-2.0*j*z[1]), den=1.0+q;
                const Real sj=(ep-em)/den, cj=(ep+em)/den;
                const Real ds=2.0*j*(em+q*ep)/(den*den), dc=2.0*j*(q*ep-em)/(den*den);
                S[(size_t)j-1]=sj; C[(size_t)j-1]=cj;
                u+=j*bj*cj*co; v+=j*bj*sj*si;
                at(sr,1)+=bj*ds*co;
                duk+=j*bj*dc*co; dvk+=j*bj*ds*si;
                duy+=j*j*bj*sj*co; dvy+=j*j*bj*cj*si;
                at(sr,n+10+j)=sj*co;
            }
            const Real uw=u-z[7];
            at(sr,7)=-surface; at(sr,8)=-1.0; at(sr,10+m)=uw;
            at(dr,1)=uw*duk+v*dvk;
            at(dr,7)=-uw; at(dr,9)=-1.0; at(dr,10+m)=1.0+uw*duy+v*dvy;
            for (int j=1; j<=n; ++j) {
                const size_t idx=(size_t)m*n+j-1;
                at(dr,n+10+j)=j*(uw*C[(size_t)j-1]*cos_nm[idx]+v*S[(size_t)j-1]*sin_nm[idx]);
            }
        }
        return A;
    }

    // ----------------------------------------------------------------------
    // Newton iteration for the nonlinear residual system F(z) = 0.
    // ----------------------------------------------------------------------
    Real newton_step(int /*iter_count*/) {
        const Real ss0 = eqns(rhs1);
        if (!std::isfinite(ss0)) throw std::runtime_error("Non-finite initial residual.");
        const std::vector<Real> z0 = z;
        const auto A = jacobian();
        std::vector<Real> b((size_t)num);
        for (int i = 1; i <= num; ++i) b[(size_t)i-1] = -rhs1[(size_t)i];
        const auto dx = svd_solve(A, b, num);
        for (Real value : dx) if (!std::isfinite(value)) throw std::runtime_error("Non-finite Newton step.");
        auto best = z0;
        for (Real alpha = 1.0; alpha >= 1e-4; alpha *= 0.5) {
            auto trial = z0;
            bool valid = true;
            for (int i = 1; i <= num; ++i) {
                trial[(size_t)i] += alpha * dx[(size_t)i-1];
                valid = valid && std::isfinite(trial[(size_t)i]);
            }
            if (!valid || !(trial[1] > 0.0)) continue;
            try {
                std::vector<Real> work((size_t)n+1);
                const Real ss = eqns_for_state(trial, rhs2, work);
                if (std::isfinite(ss) && ss <= ss0) { best = std::move(trial); break; }
            } catch (const std::exception&) { continue; }
        }
        Real correction = 0.0;
        for (int i = 10; i <= n+10; ++i) correction += std::abs(best[(size_t)i]-z0[(size_t)i]);
        z = std::move(best);
        return correction/(n+1);
    }

    // ----------------------------------------------------------------------
    // Reconstruct free-surface cosine coefficients Y_j and copy Fourier coefficients B_j.
    // ----------------------------------------------------------------------
    void compute_Y_and_B() {
        std::fill(Y.begin(), Y.end(), 0.0);

        for (int j = 1; j <= n; ++j) {
            B[(size_t)j] = z[(size_t)(j + n + 10)];

            Real s = 0.5 * (z[10] + z[(size_t)(n + 10)] * (((j % 2) == 0) ? 1.0 : -1.0));
            for (int m = 1; m < n; ++m) {
                const int idx = (m * j) % (2 * n);
                s += z[(size_t)(10 + m)] * cosa[(size_t)idx];
            }
            Y[(size_t)j] = 2.0 * s / (Real)n;
        }

        // Refresh tanh(jkd) for kinematic and post-processing evaluations.
        for (int i = 1; i <= n; ++i) Tanh[(size_t)i] = std::tanh((Real)i * z[1]);
    }

    // ----------------------------------------------------------------------
    // Evaluate the free-surface ordinate zeta(X) = k eta(X) from cosine coefficients.
    // ----------------------------------------------------------------------
    Real surface_keta(Real X) const {
        Real kEta = 0.0;
        for (int j = 1; j < n; ++j) kEta += Y[(size_t)j] * std::cos((Real)j * X);
        kEta += 0.5 * Y[(size_t)n] * std::cos((Real)n * X);
        return kEta;
    }

    // ----------------------------------------------------------------------
    // Point kinematics: returns depth-scaled fixed-frame u, v and material acceleration.
    // ----------------------------------------------------------------------
    void point(Real X, Real Yloc, Real& u_out, Real& v_out, Real& dudt_out) const {
        const Real kd = z[1];
        const Real kd_sqrt = std::sqrt(kd);

        const Real c_nd  = z[4] / kd_sqrt; // c/sqrt(gd)
        const Real ce_nd = z[5] / kd_sqrt; // ubar_1/sqrt(gd), Eulerian current variable

        Real u = 0.0, v = 0.0, ux = 0.0, vx = 0.0;

        for (int j = 1; j <= n; ++j) {
            const Real Cos = std::cos((Real)j * X);
            const Real Sin = std::sin((Real)j * X);

            const Real a = (Real)j * Yloc;
            const Real b = (Real)j * kd;
            const Real ep = std::exp(a), em = std::exp(-a-2.0*b);
            const Real denominator = 1.0 + std::exp(-2.0*b);
            const Real C = (ep+em)/denominator;
            const Real S = (ep-em)/denominator;

            const Real Bj_ = B[(size_t)j];

            u  += (Real)j * Bj_ * C * Cos;
            v  += (Real)j * Bj_ * S * Sin;
            ux += -((Real)j * (Real)j) * Bj_ * C * Sin;
            vx +=  ((Real)j * (Real)j) * Bj_ * S * Cos;
        }

        const Real inv_kd_sqrt = 1.0 / kd_sqrt;

        u  *= inv_kd_sqrt;
        v  *= inv_kd_sqrt;
        ux *= kd_sqrt;
        vx *= kd_sqrt;

        u = ce_nd + u;

        const Real ut = -c_nd * ux;
        const Real uy = vx;

        const Real dudt = ut + u * ux + v * uy;

        u_out = u;
        v_out = v;
        dudt_out = dudt;
    }

    // ----------------------------------------------------------------------
    // Mean square orbital bed velocity ub^2 = <(u_b - ubar_1)^2>.
    // ----------------------------------------------------------------------
    Real mean_square_bed_orbital_velocity(int /*nph*/ = 720) {
        if (!(d > 0.0) || !(T_target > 0.0) || !(z[1] > 0.0)) return 0.0;
        Real sum = 0.0;
        for (int j = 1; j <= n; ++j) {
            const Real e = std::exp(-(Real)j*z[1]);
            const Real sech = 2.0*e/(1.0+e*e);
            const Real amplitude = (Real)j*B[(size_t)j]*sech;
            sum += amplitude*amplitude;
        }
        return 0.5*g*d/z[1]*sum;
    }

    // ----------------------------------------------------------------------
    // Depth-scaled momentum flux S/(rho g d^2), evaluated in the moving-frame formulation.
    // ----------------------------------------------------------------------
    Real momentum_flux_S_depth(Real phase = 0.0, int npts = 1200) {
        const Real kd = z[1];
        if (!(kd > 0.0)) return 0.0;

        const Real c_nd = z[4] / std::sqrt(kd);
        const Real R_nd = 1.0 + z[9] / kd;

        const Real X = phase;
        const Real kEta = surface_keta(X);
        const Real eta_over_d = 1.0 + kEta / kd; // free-surface elevation as bed-origin y/d

        const int Np = std::max(50, npts);
        const Real dy = eta_over_d / (Real)(Np - 1);

        Real integ = 0.0;
        for (int i = 0; i < Np; ++i) {
            const Real y = dy * (Real)i;
            const Real Yloc = kd * (y - 1.0); // Y = k y, with y/d measured from still-water level

            Real u_nd, v_nd, dudt_nd;
            point(X, Yloc, u_nd, v_nd, dudt_nd);

            const Real urel = u_nd - c_nd;
            const Real P = R_nd - y - 0.5 * (urel * urel + v_nd * v_nd);
            const Real f = P + urel * urel;

            const Real w = (i == 0 || i == (Np - 1)) ? 0.5 : 1.0;
            integ += w * f;
        }

        return integ * dy;
    }

    // ----------------------------------------------------------------------
    // Integral properties in the SOLUTION.RES nondimensional and dimensional convention.
    // ----------------------------------------------------------------------
    void calc_integral_props_cpp() {
        const Real kd = z[1];
        if (!(kd > 0.0)) {
            Power = EnergyDensity = Sxx = Impulse = Cg = 0.0;
            MassTransport = 0.0;
            BernoulliR = 0.0;
            return;
        }

        const Real kd_sqrt = std::sqrt(kd);

        const Real ce_dimless    = z[5] / kd_sqrt;
        const Real cs_dimless    = z[6] / kd_sqrt;
        const Real ubar_dimless  = z[7] / kd_sqrt;

        const Real kd_32 = kd * kd_sqrt;

        const Real Q_dimless = ubar_dimless - z[8] / kd_32;
        const Real R_dimless = 1.0 + z[9] / kd;

        const Real pulse = z[8] + kd * z[5];
        const Real ke = 0.5 * (z[4] * pulse - z[5] * Q_dimless * kd_32);

        Real pe = 0.0;
        for (int i = 1; i < n; ++i) pe += 0.25 * (Y[(size_t)i] * Y[(size_t)i]);
        pe += Y[(size_t)n] * Y[(size_t)n] / 16.0;

        const Real beta = 2.0*z[9] - z[4]*z[4] + 2.0*z[4]*z[5];
        const Real sxx = 4.0*ke - 3.0*pe + beta*kd - 2.0*z[5]*pulse;
        const Real f = z[4]*(3.0*ke - 2.0*pe) + 0.5*beta*(pulse+z[4]*kd) - 2.0*z[4]*z[5]*pulse;

        const Real kd2 = kd * kd;
        const Real kd25 = kd2 * kd_sqrt;

        const Real E_depth_loc = (ke + pe) / kd2;
        const Real KE_depth_loc = ke / kd2;
        const Real PE_depth_loc = pe / kd2;

        E_depth = E_depth_loc;
        KE_depth = KE_depth_loc;
        PE_depth = PE_depth_loc;

        Sxx_depth = sxx / kd2;
        F_depth = f / kd25;
        I_depth = pulse / kd_32;

        // Convert depth-scaled invariants to dimensional engineering units.
        EnergyDensity = Phys::RHO * g * (d * d) * E_depth;            // [J/m^2]
        Sxx           = Phys::RHO * g * (d * d) * Sxx_depth;          // [N/m]
        Power         = Phys::RHO * std::pow(g, 1.5) * std::pow(d, 2.5) * F_depth; // [W/m]

        MomentumFluxDepth = momentum_flux_S_depth(0.0, 1200);
        MomentumFlux      = Phys::RHO * g * (d * d) * MomentumFluxDepth;           // [N/m]

        Impulse = Phys::RHO * std::sqrt(g * std::pow(d, 3.0)) * I_depth;           // [kg/(m·s)]

        BernoulliR = R_dimless * g * d;
        MassTransport = cs_dimless * std::sqrt(g * d);

        // Convenience values for the output tables and glossary.
        EulerianCurrent = ce_dimless * std::sqrt(g*d);
        StokesCurrent = MassTransport;
        MeanFluidSpeed = ubar_dimless * std::sqrt(g * d);

        VolumeFluxQ = Q_dimless * std::sqrt(g * std::pow(d, 3.0));
        WaveVolumeFlux_q = MeanFluidSpeed * d - VolumeFluxQ;

        Bernoulli_r = (R_dimless - 1.0) * g * d;

        KineticEnergy   = Phys::RHO * g * (d * d) * KE_depth; // [J/m^2]
        PotentialEnergy = Phys::RHO * g * (d * d) * PE_depth; // [J/m^2]

        MeanSquareBedVelocity = mean_square_bed_orbital_velocity(720);
        MeanSquareBedAbsoluteVelocity = MeanSquareBedVelocity + EulerianCurrent*EulerianCurrent;

        Cg = (std::abs(EnergyDensity) > 1e-12) ? (Power / EnergyDensity) : 0.0;
    }

    // ----------------------------------------------------------------------
    // Core continuation, Newton solution of F(z)=0, and dimensional post-processing.
    // ----------------------------------------------------------------------
    void solve_internal() {
        Real fraction=0.0, increment=1.0/std::max(1,nstep), previous_fraction=0.0;
        std::vector<Real> accepted, previous;
        int attempts=0;
        while (fraction < 1.0) {
            ++attempts;
            const Real target=std::min(1.0, fraction+increment);
            height=target*Height; Hoverd=target*MaxH;
            if (accepted.empty()) {
                init_linear();
            } else {
                z=accepted;
                if (!previous.empty()) {
                    const Real ratio=(target-fraction)/(fraction-previous_fraction);
                    for (int i=1; i<=num; ++i) z[(size_t)i]+=ratio*(accepted[(size_t)i]-previous[(size_t)i]);
                }
                if (!(z[1]>0.0) || !std::all_of(z.begin(),z.end(),[](Real v){return std::isfinite(v);})) z=accepted;
            }
            bool step_converged=false;
            const Real residual_limit=target==1.0 ? kFinalResidualTolerance : kIntermediateResidualTolerance;
            const Real correction_limit=target==1.0 ? criter_final : crit;
            for (int it=1; it<=number; ++it) {
                Real err;
                try {
                    err=newton_step(it);
                    if (!std::isfinite(eqns(rhs1))) break;
                    max_residual=residual_inf_norm(rhs1);
                } catch (const std::exception&) { break; }
                if (it>1 && err<correction_limit*std::abs(z[1]) && max_residual<=residual_limit) {
                    step_converged=true; break;
                }
                if (err==0.0 && max_residual>residual_limit) break;
            }
            if (!step_converged) {
                increment*=0.5;
                if (increment<1.0/4096.0 || attempts>=80) {
                    std::ostringstream msg;
                    msg << "No converged solution at H=" << H_target << " m, T=" << T_target
                        << " s, Uc=" << Uc << " m/s. Height continuation reached H="
                        << fraction*H_target << " m; max|F|=" << std::scientific << std::setprecision(3)
                        << max_residual << ". The requested wave may be near or beyond the limiting wave height; reduce H or increase T.";
                    throw std::runtime_error(msg.str());
                }
                continue;
            }
            previous_fraction=fraction; previous=accepted;
            fraction=target; accepted=z; continuation_fraction=fraction;
            increment=std::min(increment*1.5,1.0/std::max(1,nstep));
            compute_Y_and_B();
        }

        // Evaluate the complete nonlinear system and record the spectral tail.
        {
            const Real final_ss = eqns(rhs1);
            if (!std::isfinite(final_ss)) {
                throw std::runtime_error("Final nonlinear residual norm is non-finite.");
            }
            max_residual = residual_inf_norm(rhs1);
            if (max_residual > kFinalResidualTolerance) {
                std::ostringstream msg;
                msg.setf(std::ios::scientific);
                msg << std::setprecision(3)
                    << "False convergence rejected: maximum residual = " << max_residual
                    << ", tolerance = " << kFinalResidualTolerance << ".";
                throw std::runtime_error(msg.str());
            }

            surface_tail_coefficient = std::abs(Y[(size_t)n]);
        }

        // ------------------------- dimensional post-process --------------------
        const Real kd = z[1];
        if (!(kd > 0.0) || !std::isfinite(kd)) throw std::runtime_error("Invalid wavenumber (kd).");

        const Real k_phys = kd / d;
        const Real L_phys = 2.0 * Phys::PI / k_phys;
        const Real c_phys = L_phys / T_target;

        if (!(L_phys > 0.0) || !std::isfinite(L_phys)) throw std::runtime_error("Invalid wavelength.");
        if (!(c_phys > 0.0) || !std::isfinite(c_phys)) throw std::runtime_error("Invalid celerity.");

        // Surface nodes correspond to X_m = m*pi/N over the symmetric half-wave crest-to-trough.
        for (int m = 0; m <= n; ++m) {
            const Real kEta = z[(size_t)(10 + m)]; // zeta_m = k eta_m at node
            eta_nodes[(size_t)m] = d * (1.0 + kEta / kd);
        }

        k = k_phys;
        L = L_phys;
        c = c_phys;

        const Real z3_period = T_target * std::sqrt(g * k_phys);
        const Real z4_period = c_phys * std::sqrt(k_phys / g);
        if (std::isfinite(z3_period) && std::isfinite(z4_period) && z3_period > 0.0 && z4_period > 0.0) {
            z[3] = z3_period;
            z[4] = z4_period;

            if (Current_criterion == 1) {
                z[5] = Uc * std::sqrt(k_phys / g);
                z[7] = z[4] - z[5];
                z[6] = z[8] / z[1] - z[7] + z[4];
            } else {
                z[6] = Uc * std::sqrt(k_phys / g);
                z[7] = z[8] / z[1] - z[6] + z[4];
                z[5] = z[4] - z[7];
            }
        }

        // Any dimensional consistency normalization must leave the complete
        // nonlinear solution unchanged within the same residual tolerance.
        {
            const Real normalized_ss = eqns(rhs1);
            if (!std::isfinite(normalized_ss)) {
                throw std::runtime_error("Non-finite residual after dimensional consistency normalization.");
            }
            max_residual = residual_inf_norm(rhs1);
            if (max_residual > kFinalResidualTolerance) {
                throw std::runtime_error(
                    "Dimensional consistency normalization invalidated the nonlinear solution.");
            }
        }

        // Store B_j as a zero-based C++ array for reporting and kinematics.
        for (int j = 1; j <= n; ++j) Bj[(size_t)(j - 1)] = B[(size_t)j];

        eta_crest  = eta_nodes[0] - d;
        eta_trough = eta_nodes[(size_t)n] - d;

        steepness = H_target / L;
        rel_depth = d / L;
        ursell    = (H_target * L * L) / (d * d * d);

        if (rel_depth < 0.05) regime = "Shallow";
        else if (rel_depth < 0.5) regime = "Intermediate";
        else regime = "Deep";

        // Fenton's wavelength-dependent theoretical highest-wave curve is the
        // steady-wave domain boundary for this steady Fourier formulation.
        const Real highest_hoverd = fenton_highest_hoverd(L / d);
        if (!(highest_hoverd > 0.0) || !std::isfinite(highest_hoverd)) {
            throw std::runtime_error("Could not evaluate the theoretical highest-wave limit.");
        }

        highest_wave_limit_fenton = highest_hoverd * d;
        fourier_accuracy_limit = kFourierApplicabilityFraction * highest_wave_limit_fenton;
        highest_wave_ratio = H_target / highest_wave_limit_fenton;
        near_limiting_wave = highest_wave_ratio >= 0.90;
        practical_breaking_risk = MaxH > 0.55;

        breaking_limit_miche = 0.142 * L * std::tanh(k*d);
        breaking_index = breaking_limit_miche > 0.0 ? H_target/breaking_limit_miche : 0.0;
        is_breaking = breaking_limit_miche > 0.0 && H_target > breaking_limit_miche;

        calc_integral_props_cpp();

        Real w, ax, u_trough;
        get_kinematics(0.0, 0.0, u_bed_crest, w, ax);
        get_kinematics(d + eta_crest, 0.0, u_surf_crest, w, ax);
        get_kinematics(d + eta_trough, Phys::PI, u_trough, w, ax);
        asymmetry = std::abs(u_trough) > 0.0 ? std::abs(u_surf_crest/u_trough) : 0.0;
        crest_stagnation_margin = c-u_surf_crest;
        compute_extrema();
        const Real cf_est = 0.005;
        tau_bed = 0.5*Phys::RHO*cf_est*u_bed*u_bed;
        ExcursionBed = phase_maximum("bed_orbital").first*T_target/(2.0*Phys::PI);

        for (Real value : {L, k, c, Power, Sxx, EnergyDensity, Cg, Impulse,
                           PotentialEnergy, KineticEnergy, MeanSquareBedVelocity,
                           MomentumFlux, u_bed, u_surf, acc_max, acc_bed_max, w_max}) {
            if (!std::isfinite(value)) throw std::runtime_error("Non-finite dimensional output.");
        }
        if (!std::isfinite(eqns(rhs1))) throw std::runtime_error("Non-finite final residual.");
        max_residual = residual_inf_norm(rhs1);
        if (max_residual > kFinalResidualTolerance) throw std::runtime_error("Final residual exceeds tolerance.");
    }
};

// ==============================================================================
//  OUTPUT FORMATTING (SOLUTION.RES-style report)
// ==============================================================================

namespace ReportFmt {

static constexpr int W = 107; // fixed report width, including table borders

// Count UTF-8 codepoints for fixed-width SOLUTION.RES-style tables.
static size_t utf8_len(const std::string& s) {
    size_t n = 0;
    for (unsigned char c : s) {
        if ((c & 0xC0) != 0x80) ++n;
    }
    return n;
}

static std::string utf8_trunc(const std::string& s, size_t max_cp) {
    if (utf8_len(s) <= max_cp) return s;
    size_t cp = 0;
    size_t i = 0;
    for (; i < s.size(); ++i) {
        unsigned char c = (unsigned char)s[i];
        if ((c & 0xC0) != 0x80) {
            if (cp == max_cp) break;
            ++cp;
        }
    }
    return s.substr(0, i);
}

static std::string pad_left(const std::string& s, int w) {
    std::string t = utf8_trunc(s, (size_t)std::max(0, w));
    const int len = (int)utf8_len(t);
    if (len >= w) return t;
    return std::string((size_t)(w - len), ' ') + t;
}

static std::string pad_right(const std::string& s, int w) {
    std::string t = utf8_trunc(s, (size_t)std::max(0, w));
    const int len = (int)utf8_len(t);
    if (len >= w) return t;
    return t + std::string((size_t)(w - len), ' ');
}

static std::string pad_center(const std::string& s, int w, char fill = ' ') {
    std::string t = utf8_trunc(s, (size_t)std::max(0, w));
    const int len = (int)utf8_len(t);
    if (len >= w) return t;
    const int left = (w - len) / 2;
    const int right = w - len - left;
    return std::string((size_t)left, fill) + t + std::string((size_t)right, fill);
}

static std::string py_str_float(double v) {
    if (!std::isfinite(v)) return "nan";
    std::ostringstream ss;
    ss.setf(std::ios::fmtflags(0), std::ios::floatfield);
    ss << std::setprecision(15) << v;
    std::string s = ss.str();
    const bool has_exp = (s.find('e') != std::string::npos) || (s.find('E') != std::string::npos);
    const bool has_dot = (s.find('.') != std::string::npos);
    if (!has_exp && !has_dot) s += ".0";
    // Normalize non-finite numeric text for stable report rendering.
    if (s == "inf" || s == "+inf") return "inf";
    if (s == "-inf") return "-inf";
    return s;
}

static void hline(std::ostringstream& out, char ch = '-') {
    out << "+" << std::string((size_t)(W - 2), ch) << "+\n";
}

static void box_title(std::ostringstream& out, const std::string& title) {
    hline(out, '-');
    out << "|" << pad_center(title, W - 2) << "|\n";
    hline(out, '-');
}

static void box_text(std::ostringstream& out, const std::string& text) {
    std::string s = text;
    for (char& c : s) if (c == '\n') c = ' ';
    while (!s.empty() && (s.back() == '\r' || s.back() == ' ')) s.pop_back();

    if (utf8_len(s) > (size_t)(W - 2)) {
        // python: s[:W-5] + "..."
        s = utf8_trunc(s, (size_t)(W - 5)) + "...";
    }
    out << "|" << pad_right(s, W - 2) << "|\n";
}

static std::string fmt_float(double v, int w) {
    if (std::isnan(v)) return pad_left("nan", w);
    if (!std::isfinite(v)) return pad_left("nan", w);

    // fixed-point, decreasing decimals
    for (int dec : {5, 4, 3, 2, 1, 0}) {
        std::ostringstream ss;
        ss.setf(std::ios::fixed);
        ss << std::setprecision(dec) << v;
        std::string s = ss.str();
        if ((int)utf8_len(s) <= w) return pad_left(s, w);
    }

    // scientific fallback
    for (int sig : {6, 5, 4, 3}) {
        std::ostringstream ss;
        ss.setf(std::ios::scientific);
        ss << std::setprecision(sig) << v;
        std::string s = ss.str();
        if ((int)utf8_len(s) <= w) return pad_left(s, w);
    }

    std::ostringstream ss;
    ss.setf(std::ios::scientific);
    ss << std::setprecision(2) << v;
    std::string s = ss.str();
    s = utf8_trunc(s, (size_t)w);
    return pad_left(s, w);
}

struct Cell {
    enum Kind { NONE, NUM, STR } kind = NONE;
    double num = 0.0;
    std::string str;

    static Cell none() { return Cell(); }
    static Cell numv(double v) { Cell c; c.kind = NUM; c.num = v; return c; }
    static Cell strv(const std::string& s) { Cell c; c.kind = STR; c.str = s; return c; }
};

static std::string fmt_cell(const Cell& c, int w, const std::string& align) {
    std::string s;
    if (c.kind == Cell::NONE) {
        s = "-";
    } else if (c.kind == Cell::NUM) {
        return fmt_float(c.num, w);
    } else {
        s = c.str;
    }

    for (char& ch : s) if (ch == '\n') ch = ' ';
    // strip
    while (!s.empty() && (s.front() == ' ' || s.front() == '\t')) s.erase(s.begin());
    while (!s.empty() && (s.back() == ' '  || s.back() == '\t' || s.back() == '\r')) s.pop_back();

    s = utf8_trunc(s, (size_t)w);

    if (align == "left")   return pad_right(s, w);
    if (align == "center") return pad_center(s, w);
    return pad_left(s, w);
}

static void table_sep(std::ostringstream& out, const std::vector<int>& col_w) {
    out << "|";
    for (size_t i = 0; i < col_w.size(); ++i) {
        out << std::string((size_t)(col_w[i] + 2), '-');
        if (i + 1 < col_w.size()) out << "+";
    }
    out << "|\n";
}

static std::vector<std::string> wrap_text(const std::string& s_in, int width) {
    std::string s = s_in;
    for (char& c : s) if (c == '\n') c = ' ';
    while (!s.empty() && (s.back() == '\r' || s.back() == ' ')) s.pop_back();

    // trivial case
    if ((int)utf8_len(s) <= width) return { s };

    std::istringstream iss(s);
    std::string word;
    std::vector<std::string> lines;
    std::string line;

    auto line_len = [&](const std::string& x) { return (int)utf8_len(x); };

    while (iss >> word) {
        if (line.empty()) {
            line = word;
        } else {
            const int candidate = line_len(line) + 1 + line_len(word);
            if (candidate <= width) {
                line += " " + word;
            } else {
                lines.push_back(line);
                line = word;
            }
        }

        // If a word exceeds width (break_long_words=False), we keep it as-is;
        // fmt_cell will truncate it later to preserve fixed-width table behavior.
        if (!line.empty() && line_len(line) > width) {
            lines.push_back(line);
            line.clear();
        }
    }
    if (!line.empty()) lines.push_back(line);
    if (lines.empty()) lines.push_back("");

    return lines;
}

struct Row {
    bool is_section = false;
    std::string section;
    std::vector<Cell> cells;
};

static void print_table(std::ostringstream& out,
                        const std::string& title,
                        const std::vector<std::string>& headers,
                        const std::vector<int>& col_w,
                        const std::vector<std::string>& aligns,
                        const std::vector<Row>& rows)
{
    box_title(out, title);

    // Table header
    std::string line = "|";
    for (size_t i = 0; i < headers.size(); ++i) {
        Cell c = Cell::strv(headers[i]);
        line += " " + fmt_cell(c, col_w[i], aligns[i]) + " |";
    }
    out << line << "\n";
    table_sep(out, col_w);

    // Table body
    for (const auto& r : rows) {
        if (r.is_section) {
            const std::string sec = " " + r.section + " ";
            out << "|" << pad_center(sec, W - 2, '-') << "|\n";
            continue;
        }

        std::vector<std::vector<std::string>> wrapped;
        wrapped.reserve(r.cells.size());
        size_t nlines = 1;

        for (size_t i = 0; i < r.cells.size(); ++i) {
            const Cell& c = r.cells[i];
            if (c.kind == Cell::NUM || c.kind == Cell::NONE) {
                wrapped.push_back({ "" }); // placeholder; numeric rendered via fmt_cell
            } else {
                auto wl = wrap_text(c.str, col_w[i]);
                nlines = std::max(nlines, wl.size());
                wrapped.push_back(std::move(wl));
            }
        }

        for (size_t li = 0; li < nlines; ++li) {
            std::string l = "|";
            for (size_t ci = 0; ci < r.cells.size(); ++ci) {
                Cell cell = r.cells[ci];
                if (cell.kind == Cell::STR) {
                    cell.str = (li < wrapped[ci].size()) ? wrapped[ci][li] : "";
                } else if (cell.kind == Cell::NUM || cell.kind == Cell::NONE) {
                    // numeric/no-value: only on first line
                    if (li > 0) cell = Cell::strv("");
                }
                l += " " + fmt_cell(cell, col_w[ci], aligns[ci]) + " |";
            }
            out << l << "\n";
        }
    }

    hline(out, '-');
    out << "\n";
}

} // namespace ReportFmt

// ----------------------------------------------------------------------------
// Report generator: no-current and Eulerian-current cases under the convention.
// ----------------------------------------------------------------------------
[[maybe_unused]] static std::string generate_output(double H_in, double T_in, double d_in, double Uc_in) {
    using namespace ReportFmt;

    const bool has_current = (Uc_in != 0.0);
    FentonStreamFunction solver0(H_in, T_in, d_in, 0.0);
    FentonStreamFunction solverC(H_in, T_in, d_in, Uc_in);

    if (has_current) {
#ifdef _OPENMP
        const int total_threads = openmp_thread_count();
        if (total_threads >= 2) {
            const int threads_no_current = (total_threads + 1) / 2;
            const int threads_with_current = total_threads / 2;

#pragma omp parallel sections num_threads(2)
            {
#pragma omp section
                {
                    omp_set_num_threads(threads_no_current);
                    solver0.solve();
                }
#pragma omp section
                {
                    omp_set_num_threads(threads_with_current);
                    solverC.solve();
                }
            }
        } else {
            solver0.solve();
            solverC.solve();
        }
#else
        solver0.solve();
        solverC.solve();
#endif
    } else {
        // A single case receives the complete OpenMP team.
        solver0.solve();
    }

    // Numerical sanity checks for convergence and finite nonlinear-wave outputs.
    auto solver_issue = [&](const FentonStreamFunction& s, const char* label) -> std::string {
        if (!s.converged) {
            if (!s.last_error.empty()) return std::string("[") + label + "] " + s.last_error;
            return std::string("[") + label + "] Did not converge.";
        }
        for (const Real v : { s.L, s.k, s.c }) {
            if (!std::isfinite(v)) return std::string("[") + label + "] Non-finite result.";
        }
        return "";
    };

    std::vector<std::string> issues;
    {
        auto s0 = solver_issue(solver0, "No current");
        if (!s0.empty()) issues.push_back(s0);
        if (has_current) {
            auto sC = solver_issue(solverC, "With current");
            if (!sC.empty()) issues.push_back(sC);
        }
    }

    std::ostringstream out;

    const bool valid0 = solver_issue(solver0, "No current").empty();
    const bool validC = has_current && solver_issue(solverC, "With current").empty();
    if (!valid0 && !validC) {
        box_title(out, "INVALID, OUT-OF-DOMAIN, OR NON-CONVERGED WAVE");
        for (const auto& msg : issues) out << msg << "\n";
        hline(out, '-');
        return out.str();
    }

    // Shared input scales come from a converged state; failed columns stay empty.
    if (!valid0) solver0 = solverC;
    if (!validC) solverC = solver0;

    // ------------------------------- report header ------------------------------
    box_title(out, "NONLINEAR WAVE HYDRODYNAMICS SOLVER (FENTON)");
    box_text(out, std::string("Regular wave height (H)     : ") + py_str_float(H_in) + " m");
    box_text(out, std::string("Wave period (τ)             : ") + py_str_float(T_in) + " s");
    box_text(out, std::string("Water depth (d)             : ") + py_str_float(d_in) + " m");
    box_text(out, std::string("Eulerian current ū₁         : ") + py_str_float(Uc_in) + " m/s (positive with wave propagation)");
    hline(out, '-');
    box_text(out, issues.empty() ? "Status: Finite-system and output accuracy checks passed."
                                 : "Status: Results shown for accepted cases only.");
    auto show_status = [&](bool valid, const FentonStreamFunction& slv, const char* label) {
        if (!valid) return;
        std::ostringstream status;
        status << label << ": N=" << slv.N << "; max|F|=" << std::scientific
               << std::setprecision(3) << slv.max_residual;
        box_text(out, status.str());
        if (slv.finite_system_verified) {
            std::ostringstream bound;
            bound << "Finite-system interval residual bound: " << std::scientific << std::setprecision(3)
                  << slv.finite_system_certificate.residual_upper << "; nine output bounds checked.";
            box_text(out,bound.str());
        }
        if (!slv.numerical_note.empty()) box_text(out, slv.numerical_note);
    };
    show_status(valid0, solver0, "No current");
    show_status(validC, solverC, "With current");
    hline(out, '-');
    for (const auto& msg : issues) out << msg << "\n";
    out << "\n";

    // ---------------------------- hydrodynamic summary ---------------------------
    const Real g = solver0.g;
    const Real d = solver0.d;
    const Real sqrt_gd = std::sqrt(g * d);
    const Real sqrt_g_over_d = std::sqrt(g / d);

    auto wc_num = [&](double v) -> Cell { return has_current ? Cell::numv(v) : Cell::strv("-"); };
    auto wc_str = [&](const std::string& s) -> Cell { return has_current ? Cell::strv(s) : Cell::strv("-"); };
    auto sci_text = [](double v) -> std::string {
        if (!std::isfinite(v)) return "nan";
        std::ostringstream ss;
        ss.setf(std::ios::scientific);
        ss << std::setprecision(3) << v;
        return ss.str();
    };
    auto wc_sci = [&](double v) -> Cell { return has_current ? Cell::strv(sci_text(v)) : Cell::strv("-"); };
    auto report_celerity = [](const FentonStreamFunction& slv) -> Real {
        if (slv.T_target > 0.0 && std::isfinite(slv.L)) return slv.L / slv.T_target;
        return slv.c;
    };

    const Real c0_report = report_celerity(solver0);
    const Real cC_report = report_celerity(solverC);

    const std::vector<std::string> headers = { "PARAMETER", "NO CURRENT", "WITH CURRENT", "UNIT" };
    const std::vector<int> col_w = { 42, 16, 16, 20 };
    const std::vector<std::string> aligns = { "left", "right", "right", "left" };

    std::vector<Row> rows;

    rows.push_back(Row{ true, "INPUTS & REFERENCE SCALES", {} });
    rows.push_back(Row{ false, "", { Cell::strv("Water depth (d)"), Cell::numv(solver0.d), wc_num(solverC.d), Cell::strv("m") } });
    rows.push_back(Row{ false, "", { Cell::strv("Regular wave height (H)"), Cell::numv(solver0.H_target), wc_num(solverC.H_target), Cell::strv("m") } });
    rows.push_back(Row{ false, "", { Cell::strv("Wave period (τ)"), Cell::numv(solver0.T_target), wc_num(solverC.T_target), Cell::strv("s") } });
    rows.push_back(Row{ false, "", { Cell::strv("H/d"), Cell::numv(solver0.H_target / solver0.d), wc_num(solverC.H_target / solverC.d), Cell::strv("-") } });
    rows.push_back(Row{ false, "", { Cell::strv("τ√(g/d)"), Cell::numv(solver0.T_target * sqrt_g_over_d), wc_num(solverC.T_target * sqrt_g_over_d), Cell::strv("-") } });

    rows.push_back(Row{ true, "DISPERSION & PHASE (GEOMETRY)", {} });
    rows.push_back(Row{ false, "", { Cell::strv("Wavelength (L)"), Cell::numv(solver0.L), wc_num(solverC.L), Cell::strv("m") } });
    rows.push_back(Row{ false, "", { Cell::strv("Wave number (k)"), Cell::numv(solver0.k), wc_num(solverC.k), Cell::strv("rad/m") } });
    rows.push_back(Row{ false, "", { Cell::strv("kd"), Cell::numv(solver0.kd_dimless()), wc_num(solverC.kd_dimless()), Cell::strv("-") } });
    rows.push_back(Row{ false, "", { Cell::strv("Angular frequency (ω)"), Cell::numv(2.0 * Phys::PI / solver0.T_target), wc_num(2.0 * Phys::PI / solverC.T_target), Cell::strv("rad/s") } });
    rows.push_back(Row{ false, "", { Cell::strv("Celerity / phase speed (c)"), Cell::numv(c0_report), wc_num(cC_report), Cell::strv("m/s") } });
    rows.push_back(Row{ false, "", { Cell::strv("c/√(gd)"), Cell::numv(c0_report / sqrt_gd), wc_num(cC_report / sqrt_gd), Cell::strv("-") } });
    rows.push_back(Row{ false, "", { Cell::strv("Crest elevation (ηc)"), Cell::numv(solver0.eta_crest), wc_num(solverC.eta_crest), Cell::strv("m") } });
    rows.push_back(Row{ false, "", { Cell::strv("Trough elevation (ηt)"), Cell::numv(solver0.eta_trough), wc_num(solverC.eta_trough), Cell::strv("m") } });

    rows.push_back(Row{ true, "MEAN FLOWS (FENTON SOLUTION-FLAT)", {} });
    rows.push_back(Row{ false, "", { Cell::strv("Eulerian current (ū₁)"), Cell::numv(solver0.EulerianCurrent), wc_num(solverC.EulerianCurrent), Cell::strv("m/s") } });
    rows.push_back(Row{ false, "", { Cell::strv("Stokes current (ū₂)"), Cell::numv(solver0.StokesCurrent), wc_num(solverC.StokesCurrent), Cell::strv("m/s") } });
    rows.push_back(Row{ false, "", { Cell::strv("Mean fluid speed (Ū)"), Cell::numv(solver0.MeanFluidSpeed), wc_num(solverC.MeanFluidSpeed), Cell::strv("m/s") } });

    rows.push_back(Row{ true, "FLUXES & BERNOULLI CONSTANTS", {} });
    rows.push_back(Row{ false, "", { Cell::strv("Wave volume flux (q)"), Cell::numv(solver0.WaveVolumeFlux_q), wc_num(solverC.WaveVolumeFlux_q), Cell::strv("m²/s") } });
    rows.push_back(Row{ false, "", { Cell::strv("Volume flux (Q)"), Cell::numv(solver0.VolumeFluxQ), wc_num(solverC.VolumeFluxQ), Cell::strv("m²/s") } });
    rows.push_back(Row{ false, "", { Cell::strv("Bernoulli constant (R)"), Cell::numv(solver0.BernoulliR), wc_num(solverC.BernoulliR), Cell::strv("m²/s²") } });
    rows.push_back(Row{ false, "", { Cell::strv("Reduced Bernoulli (r = R−g d)"), Cell::numv(solver0.Bernoulli_r), wc_num(solverC.Bernoulli_r), Cell::strv("m²/s²") } });

    rows.push_back(Row{ true, "INTEGRAL QUANTITIES (PER UNIT CREST WIDTH)", {} });
    rows.push_back(Row{ false, "", { Cell::strv("Kinetic energy (T)"), Cell::numv(solver0.KineticEnergy / 1000.0), wc_num(solverC.KineticEnergy / 1000.0), Cell::strv("kJ/m²") } });
    rows.push_back(Row{ false, "", { Cell::strv("Potential energy (V)"), Cell::numv(solver0.PotentialEnergy / 1000.0), wc_num(solverC.PotentialEnergy / 1000.0), Cell::strv("kJ/m²") } });
    rows.push_back(Row{ false, "", { Cell::strv("Total energy (E = T+V)"), Cell::numv(solver0.EnergyDensity / 1000.0), wc_num(solverC.EnergyDensity / 1000.0), Cell::strv("kJ/m²") } });
    rows.push_back(Row{ false, "", { Cell::strv("Momentum flux (S)"), Cell::numv(solver0.MomentumFlux / 1000.0), wc_num(solverC.MomentumFlux / 1000.0), Cell::strv("kN/m") } });
    rows.push_back(Row{ false, "", { Cell::strv("Radiation stress (Sₓₓ)"), Cell::numv(solver0.Sxx / 1000.0), wc_num(solverC.Sxx / 1000.0), Cell::strv("kN/m") } });
    rows.push_back(Row{ false, "", { Cell::strv("Impulse (I)"), Cell::numv(solver0.Impulse / 1000.0), wc_num(solverC.Impulse / 1000.0), Cell::strv("10³ kg/(m·s)") } });
    rows.push_back(Row{ false, "", { Cell::strv("Wave power (F)"), Cell::numv(solver0.Power / 1000.0), wc_num(solverC.Power / 1000.0), Cell::strv("kW/m") } });
    rows.push_back(Row{ false, "", { Cell::strv("Energy-flux velocity (F/E)"), Cell::numv(solver0.Cg), wc_num(solverC.Cg), Cell::strv("m/s") } });

    rows.push_back(Row{ true, "KINEMATICS (EXTREMES / BED ORBITAL MOTION)", {} });
    rows.push_back(Row{ false, "", { Cell::strv("Surface max |u|"), Cell::numv(solver0.u_surf), wc_num(solverC.u_surf), Cell::strv("m/s") } });
    rows.push_back(Row{ false, "", { Cell::strv("Bed max |u|"), Cell::numv(solver0.u_bed), wc_num(solverC.u_bed), Cell::strv("m/s") } });
    rows.push_back(Row{ false, "", { Cell::strv("Surface max |aₓ|"), Cell::numv(solver0.acc_max), wc_num(solverC.acc_max), Cell::strv("m/s²") } });
    rows.push_back(Row{ false, "", { Cell::strv("Bed max |aₓ|"), Cell::numv(solver0.acc_bed_max), wc_num(solverC.acc_bed_max), Cell::strv("m/s²") } });
    rows.push_back(Row{ false, "", { Cell::strv("Velocity asymmetry |uc|/|ut|"), Cell::numv(solver0.asymmetry), wc_num(solverC.asymmetry), Cell::strv("-") } });
    rows.push_back(Row{ false, "", { Cell::strv("Mean square bed orbital vel ub²"), Cell::numv(solver0.MeanSquareBedVelocity), wc_num(solverC.MeanSquareBedVelocity), Cell::strv("m²/s²") } });
    rows.push_back(Row{ false, "", { Cell::strv("Bed orbital RMS velocity ub,rms"), Cell::numv(std::sqrt(std::max(0.0, solver0.MeanSquareBedVelocity))), wc_num(std::sqrt(std::max(0.0, solverC.MeanSquareBedVelocity))), Cell::strv("m/s") } });

    rows.push_back(Row{ true, "NONLINEARITY / BREAKING DIAGNOSTICS", {} });
    rows.push_back(Row{ false, "", { Cell::strv("Miche breaking limit (Hmax)"), Cell::numv(solver0.breaking_limit_miche), wc_num(solverC.breaking_limit_miche), Cell::strv("m") } });
    rows.push_back(Row{ false, "", { Cell::strv("Saturation (H/Hmax)"), Cell::numv(solver0.breaking_index), wc_num(solverC.breaking_index), Cell::strv("-") } });
    rows.push_back(Row{ false, "", { Cell::strv("Breaking status"), Cell::strv(solver0.is_breaking ? "BREAKING" : "BELOW MICHE"), wc_str(solverC.is_breaking ? "BREAKING" : "BELOW MICHE"), Cell::strv("-") } });
    const std::string theoretical0 = solver0.highest_wave_ratio >= 1.0 ? "ABOVE LIMIT" : (solver0.near_limiting_wave ? "NEAR LIMIT" : "BELOW LIMIT");
    const std::string theoreticalC = has_current
        ? (solverC.highest_wave_ratio >= 1.0 ? "ABOVE LIMIT" : (solverC.near_limiting_wave ? "NEAR LIMIT" : "BELOW LIMIT"))
        : "-";
    const std::string practical0 = solver0.practical_breaking_risk ? "ELEVATED RISK" : "H/d <= 0.55";
    const std::string practicalC = has_current
        ? (solverC.practical_breaking_risk ? "ELEVATED RISK" : "H/d <= 0.55")
        : "-";
    rows.push_back(Row{ false, "", { Cell::strv("Fenton highest steady wave (Hm)"), Cell::numv(solver0.highest_wave_limit_fenton), wc_num(solverC.highest_wave_limit_fenton), Cell::strv("m") } });
    rows.push_back(Row{ false, "", { Cell::strv("Limiting-wave ratio (H/Hm)"), Cell::numv(solver0.highest_wave_ratio), wc_num(solverC.highest_wave_ratio), Cell::strv("-") } });
    rows.push_back(Row{ false, "", { Cell::strv("98% Fenton limiting height"), Cell::numv(solver0.fourier_accuracy_limit), wc_num(solverC.fourier_accuracy_limit), Cell::strv("m") } });
    rows.push_back(Row{ false, "", { Cell::strv("Steady-wave domain status"), Cell::strv(theoretical0), wc_str(theoreticalC), Cell::strv("-") } });
    rows.push_back(Row{ false, "", { Cell::strv("Practical depth-breaking screen"), Cell::strv(practical0), wc_str(practicalC), Cell::strv("-") } });
    rows.push_back(Row{ false, "", { Cell::strv("Crest stagnation margin (c−ucrest)"), Cell::numv(solver0.crest_stagnation_margin), wc_num(solverC.crest_stagnation_margin), Cell::strv("m/s") } });
    rows.push_back(Row{ false, "", { Cell::strv("Maximum nonlinear residual"), Cell::strv(sci_text(solver0.max_residual)), wc_sci(solverC.max_residual), Cell::strv("-") } });
    rows.push_back(Row{ false, "", { Cell::strv("Surface spectral tail |E_N|"), Cell::strv(sci_text(solver0.surface_tail_coefficient)), wc_sci(solverC.surface_tail_coefficient), Cell::strv("-") } });
    rows.push_back(Row{ false, "", { Cell::strv("Ursell number (U)"), Cell::numv(solver0.ursell), wc_num(solverC.ursell), Cell::strv("-") } });
    rows.push_back(Row{ false, "", { Cell::strv("Regime (by d/L)"), Cell::strv(solver0.regime), wc_str(solverC.regime), Cell::strv("-") } });

    for (auto& row : rows) {
        if (row.cells.size() == 4) {
            if (!valid0) row.cells[1] = Cell::none();
            if (!validC) row.cells[2] = Cell::none();
        }
    }
    print_table(out, "CALCULATED HYDRODYNAMIC PARAMETERS", headers, col_w, aligns, rows);

    // ------------------------ dimensionless-variable tables ----------------
    auto print_solution_flat = [&](const FentonStreamFunction& slv, const std::string& title) {
        const std::vector<std::string> h = { "#", "PARAMETER", "value", "adim param", "adim value" };
        const std::vector<int> cw = { 2, 37, 13, 25, 14 };
        const std::vector<std::string> al = { "right", "left", "right", "left", "right" };

        const Real g_ = slv.g;
        const Real d_ = slv.d;
        const Real H_ = slv.H_target;
        const Real T_ = slv.T_target;
        const Real L_ = slv.L;
        const Real c_ = L_ / T_;

        const Real sqrt_gd_ = std::sqrt(g_ * d_);
        const Real sqrt_gd3_ = std::sqrt(g_ * (d_ * d_ * d_));

        auto f5 = [&](double v) {
            std::ostringstream ss; ss.setf(std::ios::fixed); ss << std::setprecision(5) << v; return ss.str();
        };

        auto kJ = [&](double J) { return J / 1000.0; };
        auto kN = [&](double N) { return N / 1000.0; };
        auto kW = [&](double W_) { return W_ / 1000.0; };

        std::vector<Row> rr;
        rr.push_back(Row{ false, "", { Cell::numv(1),  Cell::strv("Water depth"),                     Cell::strv(f5(d_)),                    Cell::strv("d/d = 1"),                         Cell::strv(f5(1.0)) } });
        rr.push_back(Row{ false, "", { Cell::numv(2),  Cell::strv("Wave length"),                     Cell::strv(f5(L_)),                    Cell::strv("λ/d"),                             Cell::strv(f5(L_ / d_)) } });
        rr.push_back(Row{ false, "", { Cell::numv(3),  Cell::strv("Wave height"),                     Cell::strv(f5(H_)),                    Cell::strv("H/d"),                             Cell::strv(f5(H_ / d_)) } });
        rr.push_back(Row{ false, "", { Cell::numv(4),  Cell::strv("Wave period"),                     Cell::strv(f5(T_)),                    Cell::strv("τ√(g/d)"),                         Cell::strv(f5(T_ * std::sqrt(g_ / d_))) } });
        rr.push_back(Row{ false, "", { Cell::numv(5),  Cell::strv("Wave speed"),                      Cell::strv(f5(c_)),                    Cell::strv("c/√(gd)"),                         Cell::strv(f5(c_ / sqrt_gd_)) } });
        rr.push_back(Row{ false, "", { Cell::numv(6),  Cell::strv("Eulerian current"),                Cell::strv(f5(slv.EulerianCurrent)),   Cell::strv("ū₁/√(gd)"),                        Cell::strv(f5(slv.EulerianCurrent / sqrt_gd_)) } });
        rr.push_back(Row{ false, "", { Cell::numv(7),  Cell::strv("Stokes current"),                  Cell::strv(f5(slv.StokesCurrent)),     Cell::strv("ū₂/√(gd)"),                        Cell::strv(f5(slv.StokesCurrent / sqrt_gd_)) } });
        rr.push_back(Row{ false, "", { Cell::numv(8),  Cell::strv("Mean fluid speed"),                Cell::strv(f5(slv.MeanFluidSpeed)),    Cell::strv("Ū/√(gd)"),                         Cell::strv(f5(slv.MeanFluidSpeed / sqrt_gd_)) } });
        rr.push_back(Row{ false, "", { Cell::numv(9),  Cell::strv("Wave volume flux, q = Ū d − Q"),   Cell::strv(f5(slv.WaveVolumeFlux_q)),  Cell::strv("q/√(gd³)"),                        Cell::strv(f5(slv.WaveVolumeFlux_q / sqrt_gd3_)) } });
        rr.push_back(Row{ false, "", { Cell::numv(10), Cell::strv("Bernoulli constant, r = R − gd"),  Cell::strv(f5(slv.Bernoulli_r)),       Cell::strv("r/gd"),                            Cell::strv(f5(slv.Bernoulli_r / (g_ * d_))) } });
        rr.push_back(Row{ false, "", { Cell::numv(11), Cell::strv("Volume flux"),                     Cell::strv(f5(slv.VolumeFluxQ)),       Cell::strv("Q/√(gd³)"),                        Cell::strv(f5(slv.VolumeFluxQ / sqrt_gd3_)) } });
        rr.push_back(Row{ false, "", { Cell::numv(12), Cell::strv("Bernoulli constant"),              Cell::strv(f5(slv.BernoulliR)),        Cell::strv("R/gd"),                            Cell::strv(f5(slv.BernoulliR / (g_ * d_))) } });
        rr.push_back(Row{ false, "", { Cell::numv(13), Cell::strv("Momentum flux"),                   Cell::strv(f5(kN(slv.MomentumFlux))),  Cell::strv("S/ρgd²"),                          Cell::strv(f5(slv.MomentumFluxDepth)) } });
        rr.push_back(Row{ false, "", { Cell::numv(14), Cell::strv("Impulse"),                         Cell::strv(f5(slv.Impulse / 1000.0)),  Cell::strv("I/(ρ√(gd³))"),                     Cell::strv(f5(slv.I_depth)) } });
        rr.push_back(Row{ false, "", { Cell::numv(15), Cell::strv("Kinetic energy"),                  Cell::strv(f5(kJ(slv.KineticEnergy))), Cell::strv("T/ρgd²"),                          Cell::strv(f5(slv.KE_depth)) } });
        rr.push_back(Row{ false, "", { Cell::numv(16), Cell::strv("Potential energy"),                Cell::strv(f5(kJ(slv.PotentialEnergy))),Cell::strv("V/ρgd²"),                          Cell::strv(f5(slv.PE_depth)) } });
        rr.push_back(Row{ false, "", { Cell::numv(17), Cell::strv("Mean square bed velocity (absolute)"),     Cell::strv(f5(slv.MeanSquareBedAbsoluteVelocity)),Cell::strv("ub²/gd"),                          Cell::strv(f5(slv.MeanSquareBedAbsoluteVelocity / (g_ * d_))) } });
        rr.push_back(Row{ false, "", { Cell::numv(18), Cell::strv("Radiation stress"),                Cell::strv(f5(kN(slv.Sxx))),           Cell::strv("S_xx/ρgd²"),                       Cell::strv(f5(slv.Sxx_depth)) } });
        rr.push_back(Row{ false, "", { Cell::numv(19), Cell::strv("Wave power"),                      Cell::strv(f5(kW(slv.Power))),         Cell::strv("F/(ρg³ᐟ²d⁵ᐟ²)ㅤㅤ"),               Cell::strv(f5(slv.F_depth)) } });

        print_table(out, title, h, cw, al, rr);
    };

    if (valid0) print_solution_flat(solver0, "SOLUTION.RES (NO CURRENT)");
    if (validC) print_solution_flat(solverC, "SOLUTION.RES (WITH CURRENT)");

    // --------------------------------- glossary --------------------------------
    const std::vector<std::string> gh = { "TERM / SYMBOL", "MEANING", "UNITS / NONDIM" };
    const std::vector<int> gcw = { 14, 64, 19 };
    const std::vector<std::string> gal = { "left", "left", "left" };

    std::vector<Row> terms;
    auto add_term = [&](const std::string& a, const std::string& b, const std::string& c) {
        terms.push_back(Row{ false, "", { Cell::strv(a), Cell::strv(b), Cell::strv(c) } });
    };

    add_term("d", "Still-water depth (bed to mean water level). Reference length scale.", "m ; d/d=1");
    add_term("H", "Regular periodic-wave height (crest-to-trough); not significant wave height Hs.", "m ; H/d");
    add_term("τ", "Wave period.", "s ; τ√(g/d)");
    add_term("L", "Wavelength (crest-to-crest).", "m ; L/d");
    add_term("k", "Wave number, k = 2π/L.", "rad/m");
    add_term("ω", "Angular frequency, ω = 2π/τ.", "rad/s");
    add_term("c", "Phase speed (celerity).", "m/s ; c/√(gd)");
    add_term("ηc, ηt", "Crest and trough elevations relative to still-water level.", "m");
    add_term("ū₁", "Eulerian (depth-mean) current; ū₁ = Uc.", "m/s ; ū₁/√(gd)");
    add_term("ū₂", "Stokes / mass-transport current from nonlinear solution.", "m/s ; ū₂/√(gd)");
    add_term("Ū", "Mean fluid speed (depth-mean).", "m/s ; Ū/√(gd)");
    add_term("Q", "Volume flux (depth-integrated).", "m²/s ; Q/√(gd³)");
    add_term("q", "Wave volume flux, q = Ū d − Q.", "m²/s ; q/√(gd³)");
    add_term("R", "Bernoulli constant.", "m²/s² ; R/(gd)");
    add_term("r", "Reduced Bernoulli constant r = R − g d.", "m²/s² ; r/(gd)");
    add_term("S", "Momentum flux (moving frame).", "kN/m ; S/(ρgd²)");
    add_term("I", "Wave impulse per unit width.", "10³ kg/(m·s) ; I/(ρ√(gd³))");
    add_term("T", "Kinetic energy density.", "kJ/m² ; T/(ρgd²)");
    add_term("V", "Potential energy density.", "kJ/m² ; V/(ρgd²)");
    add_term("E", "Total energy density E = T + V.", "kJ/m²");
    add_term("F", "Wave power (energy flux).", "kW/m ; F/(ρg³ᐟ²d⁵ᐟ²)ㅤ");
    add_term("Cg", "Energy-flux velocity Cg = F/E; energy flux divided by energy density.", "m/s");
    add_term("Sₓₓ", "Radiation stress component in wave direction.", "kN/m ; Sₓₓ/(ρgd²)");
    add_term("ub²", "Mean square *orbital* bed velocity: ub² = <(ub(t) − ū₁)²>. Non-negative by definition; evaluated by Fourier orthogonality.", "m²/s² ; /gd");
    add_term("ub,rms", "Root-mean-square orbital bed velocity: ub,rms = √(ub²).", "m/s");
    add_term("usurf,max", "Maximum horizontal velocity magnitude at the free surface, refined over phase.", "m/s");
    add_term("ubed,max", "Maximum horizontal velocity magnitude at the seabed, refined over phase.", "m/s");
    add_term("a_x,surf,max", "Maximum absolute surface material horizontal acceleration, refined over phase.", "m/s²");
    add_term("a_x,bed,max", "Maximum absolute bed material horizontal acceleration, refined over phase.", "m/s²");
    add_term("Asymmetry", "Velocity asymmetry indicator |uc|/|ut|.", "-");
    add_term("Hm", "Theoretical highest regular steady wave from Fenton (1990), Eq. (32).", "m ; H/Hm");
    add_term("E_N", "Last free-surface Fourier coefficient; spectral-resolution diagnostic.", "-");
    add_term("c−ucrest", "Crest stagnation margin; positive below the limiting wave.", "m/s");
    add_term("Ursell", "Ursell number, a shallow-water nonlinearity measure.", "-");
    add_term("Regime", "Depth regime based on d/L (deep/intermediate/shallow).", "-");

    print_table(out, "PARAMETER DEFINITIONS & GLOSSARY", gh, gcw, gal, terms);

    return out.str();
}
// Convert LF to CRLF for Win32 multiline EDIT control display.
#ifndef FENTON_HEADLESS
static std::string to_windows_newlines(const std::string& s) {
    std::string out;
    out.reserve(s.size() + s.size() / 20);
    for (size_t i = 0; i < s.size(); ++i) {
        if (s[i] == '\n') out += "\r\n";
        else out += s[i];
    }
    return out;
}

// ==============================================================================
//  GUI (Win32): inputs H, T, d and U_c; output displays the reference formulation-style report
// ==============================================================================

#define IDC_EDIT_H      101
#define IDC_EDIT_T      102
#define IDC_EDIT_D      103
#define IDC_EDIT_UC     104
#define IDC_BTN_CALC    110
#define IDC_OUTPUT      111

static HWND  g_hEditH = nullptr;
static HWND  g_hEditT = nullptr;
static HWND  g_hEditD = nullptr;
static HWND  g_hEditUc = nullptr;
static HWND  g_hOutput = nullptr;
static HWND  g_hBtnCalc = nullptr;
static HFONT g_hUIFont = nullptr;
static HFONT g_hMonoFont = nullptr;
static HWND  g_hMain = nullptr;

static std::atomic<bool> g_closing(false);

static constexpr UINT WM_APP_RESULT = WM_APP + 1;

// The formatted report uses a fixed reference formulation/SOLUTION.RES-style line width.
// The output area is sized so that the complete line is visible without clipping.
static constexpr int OUTPUT_COLS = 107;

static bool parse_double_w(const std::wstring& in, double& out_val) {
    std::wstring s = in;
    // Accept decimal comma in numeric input fields.
    std::replace(s.begin(), s.end(), L',', L'.');

    // Trim surrounding whitespace before numeric conversion.
    auto is_ws = [](wchar_t c) { return c == L' ' || c == L'\t' || c == L'\r' || c == L'\n'; };
    while (!s.empty() && is_ws(s.front())) s.erase(s.begin());
    while (!s.empty() && is_ws(s.back())) s.pop_back();

    if (s.empty()) return false;

    try {
        size_t pos = 0;
        out_val = std::stod(s, &pos);
        if (pos != s.size()) return false;
        if (!std::isfinite(out_val)) return false;
        return true;
    } catch (...) {
        return false;
    }
}

static bool read_edit_double(HWND hEdit, double& v) {
    wchar_t buf[128];
    GetWindowTextW(hEdit, buf, 128);
    return parse_double_w(buf, v);
}

struct CalcParams {
    double H, T, d, Uc;
};

static DWORD WINAPI CalcThread(LPVOID lpParam) {
    std::unique_ptr<CalcParams> p((CalcParams*)lpParam);

    // Execute the nonlinear stream-function solve and report generation.
    std::string txt = generate_output(p->H, p->T, p->d, p->Uc);

    // Write output.txt in the same SOLUTION.RES-style report format.
    try {
        std::ofstream f("output.txt", std::ios::binary);
        f.write(txt.data(), (std::streamsize)txt.size());
    } catch (...) {
        // Ignore output-file write failures; the GUI report remains available.
    }

    // Convert report text to Windows newlines before displaying it.
    std::string display_txt = to_windows_newlines(txt);

    // Convert UTF-8 report text to UTF-16 for the Win32 EDIT control.
    const int wlen = MultiByteToWideChar(CP_UTF8, 0, display_txt.c_str(), -1, nullptr, 0);
    std::wstring* wres = new std::wstring((size_t)wlen, L'\0');
    MultiByteToWideChar(CP_UTF8, 0, display_txt.c_str(), -1, &(*wres)[0], wlen);

    // Post the completed report back to the UI thread.
    if (!g_closing.load() && g_hMain && IsWindow(g_hMain)) {
        PostMessageW(g_hMain, WM_APP_RESULT, 0, (LPARAM)wres);
    } else {
        delete wres;
    }
    return 0;
}

static void resize_gui_to_fit_output(HWND hwnd) {
    const int sw = GetSystemMetrics(SM_CXSCREEN);
    const int sh = GetSystemMetrics(SM_CYSCREEN);
    const int desired_client_w = std::min(1600, sw - 80);
    const int desired_client_h = std::min(950, sh - 120);

    RECT rc{ 0, 0, desired_client_w, desired_client_h };
    const DWORD style = (DWORD)GetWindowLongPtrW(hwnd, GWL_STYLE);
    const DWORD exstyle = (DWORD)GetWindowLongPtrW(hwnd, GWL_EXSTYLE);
    AdjustWindowRectEx(&rc, style, FALSE, exstyle);
    const int win_w = rc.right - rc.left;
    const int win_h = rc.bottom - rc.top;

    SetWindowPos(hwnd, nullptr, 0, 0, win_w, win_h, SWP_NOMOVE | SWP_NOZORDER);

    if (g_hOutput && IsWindow(g_hOutput)) {
        SendMessageW(g_hOutput, EM_SETMARGINS, EC_LEFTMARGIN | EC_RIGHTMARGIN, MAKELPARAM(0, 0));
        SetWindowPos(g_hOutput, nullptr,
            300, 20,
            std::max(100, desired_client_w - 320),
            std::max(100, desired_client_h - 40),
            SWP_NOZORDER);
    }
}

static LRESULT CALLBACK WndProc(HWND hwnd, UINT msg, WPARAM wParam, LPARAM lParam) {
    switch (msg) {
    case WM_CREATE: {
        g_hMain = hwnd;

        g_hUIFont = CreateFontW(19, 0, 0, 0, FW_NORMAL, FALSE, FALSE, FALSE,
            DEFAULT_CHARSET, OUT_DEFAULT_PRECIS, CLIP_DEFAULT_PRECIS,
            DEFAULT_QUALITY, DEFAULT_PITCH | FF_SWISS, L"Segoe UI");

        g_hMonoFont = CreateFontW(22, 0, 0, 0, FW_NORMAL, FALSE, FALSE, FALSE,
            DEFAULT_CHARSET, OUT_DEFAULT_PRECIS, CLIP_DEFAULT_PRECIS,
            DEFAULT_QUALITY, FIXED_PITCH | FF_MODERN, L"Consolas");

        int y = 20;

        auto make_label = [&](const wchar_t* t, int yy) {
            HWND h = CreateWindowW(L"STATIC", t, WS_CHILD | WS_VISIBLE, 20, yy, 160, 25, hwnd, nullptr, nullptr, nullptr);
            SendMessageW(h, WM_SETFONT, (WPARAM)g_hUIFont, TRUE);
        };

        auto make_edit = [&](const wchar_t* def, int id, int yy) -> HWND {
            HWND h = CreateWindowW(L"EDIT", def, WS_CHILD | WS_VISIBLE | WS_BORDER,
                190, yy, 90, 25, hwnd, (HMENU)(INT_PTR)id, nullptr, nullptr);
            SendMessageW(h, WM_SETFONT, (WPARAM)g_hUIFont, TRUE);
            return h;
        };

        make_label(L"Regular Wave H (m):", y); g_hEditH = make_edit(L"3.0", IDC_EDIT_H, y); y += 35;
        make_label(L"Wave Period (s):", y); g_hEditT = make_edit(L"9.0", IDC_EDIT_T, y); y += 35;
        make_label(L"Water Depth (m):", y); g_hEditD = make_edit(L"5.0", IDC_EDIT_D, y); y += 35;
        make_label(L"Current (m/s):", y);   g_hEditUc = make_edit(L"1.0", IDC_EDIT_UC, y); y += 45;

        g_hBtnCalc = CreateWindowW(L"BUTTON", L"CALCULATE",
            WS_CHILD | WS_VISIBLE | BS_DEFPUSHBUTTON,
            20, y, 260, 40, hwnd, (HMENU)IDC_BTN_CALC, nullptr, nullptr);
        SendMessageW(g_hBtnCalc, WM_SETFONT, (WPARAM)g_hUIFont, TRUE);

        g_hOutput = CreateWindowW(L"EDIT", L"",
            WS_CHILD | WS_VISIBLE | WS_BORDER | ES_MULTILINE | WS_VSCROLL | ES_READONLY | WS_HSCROLL,
            300, 20, 800, 520, hwnd, (HMENU)IDC_OUTPUT, nullptr, nullptr);
        SendMessageW(g_hOutput, WM_SETFONT, (WPARAM)g_hMonoFont, TRUE);

        // Make the output area wide enough for the fixed-width SOLUTION.RES-style report.
        resize_gui_to_fit_output(hwnd);

        return 0;
    }
    case WM_COMMAND: {
        if (LOWORD(wParam) == IDC_BTN_CALC) {
            double H, T, d, Uc;
            if (!read_edit_double(g_hEditH, H) ||
                !read_edit_double(g_hEditT, T) ||
                !read_edit_double(g_hEditD, d) ||
                !read_edit_double(g_hEditUc, Uc))
            {
                MessageBoxW(hwnd, L"Invalid numeric input.\n\n"
                                  L"Accepted formats: 3.5 or 3,5 (decimal comma supported).",
                            L"Input Error", MB_ICONERROR);
                return 0;
            }

            if (!std::isfinite(H) || !std::isfinite(T) || !std::isfinite(d) ||
                !std::isfinite(Uc) || !(H > 0.0) || !(T > 0.0) || !(d > 0.0)) {
                MessageBoxW(hwnd,
                            L"H, T, and d must be finite and positive; current must be finite.",
                            L"Input Error", MB_ICONERROR);
                return 0;
            }

            EnableWindow(g_hBtnCalc, FALSE);
            SetWindowTextW(g_hBtnCalc, L"SOLVING...");

            CalcParams* p = new CalcParams{ H, T, d, Uc };

            HANDLE hTh = CreateThread(nullptr, 0, CalcThread, p, 0, nullptr);
            if (hTh) CloseHandle(hTh);
            return 0;
        }
        break;
    }
    case WM_APP_RESULT: {
        std::unique_ptr<std::wstring> pStr((std::wstring*)lParam);
        if (g_hOutput && IsWindow(g_hOutput)) {
            SetWindowTextW(g_hOutput, pStr->c_str());
        }
        if (g_hBtnCalc && IsWindow(g_hBtnCalc)) {
            EnableWindow(g_hBtnCalc, TRUE);
            SetWindowTextW(g_hBtnCalc, L"CALCULATE");
        }
        return 0;
    }
    case WM_DESTROY:
        g_closing.store(true);
        PostQuitMessage(0);
        return 0;
    }
    return DefWindowProcW(hwnd, msg, wParam, lParam);
}

int WINAPI WinMain(HINSTANCE hInst, HINSTANCE, LPSTR, int nCmdShow) {
    WNDCLASSEXW wc{};
    wc.cbSize = sizeof(WNDCLASSEXW);
    wc.lpfnWndProc = WndProc;
    wc.hInstance = hInst;
    wc.hIcon = LoadIcon(nullptr, IDI_APPLICATION);
    wc.hCursor = LoadCursor(nullptr, IDC_ARROW);
    wc.hbrBackground = (HBRUSH)(COLOR_WINDOW + 1);
    wc.lpszClassName = L"FentonClassStable";

    if (!RegisterClassExW(&wc)) return 0;

    const DWORD main_style = WS_OVERLAPPED | WS_CAPTION | WS_SYSMENU | WS_MINIMIZEBOX;
    const DWORD main_exstyle = 0;
    const int sw = GetSystemMetrics(SM_CXSCREEN);
    const int sh = GetSystemMetrics(SM_CYSCREEN);
    const int desired_client_w = std::min(1600, sw - 80);
    const int desired_client_h = std::min(950, sh - 120);
    RECT rc{ 0, 0, desired_client_w, desired_client_h };
    AdjustWindowRectEx(&rc, main_style, FALSE, main_exstyle);

    HWND hwnd = CreateWindowExW(
        main_exstyle, L"FentonClassStable", L"Fenton Wave Solver",
        main_style,
        CW_USEDEFAULT, CW_USEDEFAULT, rc.right - rc.left, rc.bottom - rc.top,
        nullptr, nullptr, hInst, nullptr
    );

    if (!hwnd) return 0;

    ShowWindow(hwnd, nCmdShow);
    UpdateWindow(hwnd);

    MSG msg;
    while (GetMessageW(&msg, nullptr, 0, 0)) {
        TranslateMessage(&msg);
        DispatchMessageW(&msg);
    }
    return 0;
}

#endif // FENTON_HEADLESS
