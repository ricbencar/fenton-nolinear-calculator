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
 *  For numerical stability the hyperbolic ratios are evaluated in the stable
 *  form
 *
 *      S_j(Y) = sinh(jY) + cosh(jY) tanh(jkd),
 *      C_j(Y) = cosh(jY) + sinh(jY) tanh(jkd).
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
 *  With the production value N = 50, this gives 110 active 1-based state entries
 *  and 110 nonlinear residual equations, matching the full full residual
 *  system for the stream-function / Fourier collocation solver.
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
 *  current iterate and solves J dz = -F(z).  Numerical finite-difference Jacobians
 *  and singular-value stabilized dense linear algebra are used because high-order
 *  stream-function waves, near-limiting waves and current cases can make the
 *  Jacobian ill-conditioned.
 *
 *  The linear Airy/Fenton-McKee estimates are used only as initialization and
 *  checking aids.  The final result is the converged nonlinear stream-function /
 *  Fourier collocation solution.
 *
 *  8. OUTPUT INTERPRETATION
 *  -----------------------------------------------------------------------------
 *  The report follows the distinction between k-based and d-based
 *  nondimensional quantities: kd, kH, T sqrt(gk), c sqrt(k/g), current variables,
 *  fluxes, Bernoulli constants, impulse, energy, momentum flux, radiation stress
 *  and wave power are reported alongside their engineering depth-scaled forms.
 *
 *  The same calculation is run for no-current and, when requested, with-current
 *  cases so the output can be compared under the selected current convention.
 *
 *  9. BUILDING FROM SOURCE
 *  -----------------------------------------------------------------------------
 *  Windows / MinGW-w64 release-like GUI build:
 *
 *      g++ fenton_gui.cpp -o fenton_gui.exe -O3 -std=c++20 -march=native ^
 *          -flto=auto -fopenmp -static -static-libgcc -static-libstdc++ ^
 *          -mwindows -pthread -lgdi32 -luser32 -lkernel32 -lcomctl32
 *
 *  Keep the same compiler family, CPU family and optimization flags when bitwise
 *  reproducibility is required, because Newton/trust-region acceptance and SVD
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

#include <windows.h>
#include <commctrl.h>

// --- nanosleep64 stub (MSYS2/MinGW GCC 15 + -static workaround) -----------------
#if defined(_WIN32) && defined(__MINGW32__)

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
#include <cstdint>
#include <iomanip>
#include <limits>
#include <memory>
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

static int openmp_thread_count() noexcept {
#ifdef _OPENMP
    const int available = std::max(1, omp_get_num_procs());
    return std::max(1, std::min(available, omp_get_max_threads()));
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

// ==============================================================================
//  SMALL LINEAR ALGEBRA (dense stabilized solves for num = 2N + 10)
// ==============================================================================

namespace LinAlg {

static inline bool is_finite(Real x) {
    return std::isfinite(x);
}

static inline Real dot(const std::vector<Real>& a, const std::vector<Real>& b) {
    Real s = 0.0;
    const size_t n = a.size();
    for (size_t i = 0; i < n; ++i) s += a[i] * b[i];
    return s;
}

static inline Real norm2_sq(const std::vector<Real>& a) {
    return dot(a, a);
}

static inline Real norm2(const std::vector<Real>& a) {
    return std::sqrt(norm2_sq(a));
}

static inline Real norm_inf(const std::vector<Real>& a) {
    Real m = 0.0;
    for (Real v : a) m = std::max(m, std::abs(v));
    return m;
}

// Solve A x = b with Gaussian elimination + partial pivoting.
// A is overwritten (row operations). b is overwritten. Returns false if singular.
static bool solve_linear_system(int n, std::vector<Real>& A, std::vector<Real>& b, std::vector<Real>& x) {
    x = b;

    for (int k = 0; k < n - 1; ++k) {
        Real max_v = 0.0;
        int max_i = k;

        for (int i = k; i < n; ++i) {
            const Real v = std::abs(A[i * n + k]);
            if (v > max_v) { max_v = v; max_i = i; }
        }
        if (max_v < 1e-30) return false;

        if (max_i != k) {
            for (int j = k; j < n; ++j) std::swap(A[k * n + j], A[max_i * n + j]);
            std::swap(x[k], x[max_i]);
        }

        const Real pivot = A[k * n + k];
        const Real inv_pivot = 1.0 / pivot;

        for (int i = k + 1; i < n; ++i) {
            const Real m = A[i * n + k] * inv_pivot;
            A[i * n + k] = 0.0;
            for (int j = k + 1; j < n; ++j) A[i * n + j] -= m * A[k * n + j];
            x[i] -= m * x[k];
        }
    }

    if (std::abs(A[(n - 1) * n + (n - 1)]) < 1e-30) return false;

    for (int i = n - 1; i >= 0; --i) {
        Real sum = 0.0;
        for (int j = i + 1; j < n; ++j) sum += A[i * n + j] * x[j];
        x[i] = (x[i] - sum) / A[i * n + i];
    }
    return true;
}



// Solve min ||A x - b||_2 for a tall (rows >= cols) matrix A using Householder QR.
// A_in is row-major of size rows*cols. b_in length rows.
// Returns false if the problem is rank-deficient / numerically singular.
static bool qr_solve_least_squares(int rows, int cols,
                                  const std::vector<Real>& A_in,
                                  const std::vector<Real>& b_in,
                                  std::vector<Real>& x_out)
{
    if (rows < cols) return false;
    if ((int)A_in.size() != rows * cols) return false;
    if ((int)b_in.size() != rows) return false;

    std::vector<Real> A = A_in; // working copy (row-major)
    std::vector<Real> b = b_in; // will become Q^T b

    std::vector<Real> v; // Householder vector (reused)
    v.reserve((size_t)rows);

    for (int k = 0; k < cols; ++k) {
        // Compute 2-norm of column k from row k..rows-1
        Real sigma = 0.0;
        for (int i = k; i < rows; ++i) {
            const Real a = A[(size_t)i * (size_t)cols + (size_t)k];
            sigma += a * a;
        }
        if (sigma <= 0.0) continue;

        const Real x0 = A[(size_t)k * (size_t)cols + (size_t)k];
        const Real normx = std::sqrt(sigma);

        // alpha = -sign(x0) * ||x||
        const Real alpha = -std::copysign(normx, x0);

        // v = x; v0 = x0 - alpha
        v.assign((size_t)(rows - k), 0.0);
        v[0] = x0 - alpha;
        for (int i = k + 1; i < rows; ++i) {
            v[(size_t)(i - k)] = A[(size_t)i * (size_t)cols + (size_t)k];
        }

        // tau = 2 / (v^T v)
        Real vTv = 0.0;
        for (Real vi : v) vTv += vi * vi;
        if (vTv <= 0.0) continue;
        const Real tau = 2.0 / vTv;

        // Apply to A columns k..cols-1
        for (int j = k; j < cols; ++j) {
            Real dot_v = 0.0;
            for (int i = 0; i < rows - k; ++i) {
                dot_v += v[(size_t)i] * A[(size_t)(k + i) * (size_t)cols + (size_t)j];
            }
            const Real s = tau * dot_v;
            for (int i = 0; i < rows - k; ++i) {
                A[(size_t)(k + i) * (size_t)cols + (size_t)j] -= s * v[(size_t)i];
            }
        }

        // Apply to b
        Real dot_b = 0.0;
        for (int i = 0; i < rows - k; ++i) {
            dot_b += v[(size_t)i] * b[(size_t)(k + i)];
        }
        const Real sb = tau * dot_b;
        for (int i = 0; i < rows - k; ++i) {
            b[(size_t)(k + i)] -= sb * v[(size_t)i];
        }

        // Explicitly set below-diagonal entries to zero (numerical hygiene)
        A[(size_t)k * (size_t)cols + (size_t)k] = alpha;
        for (int i = k + 1; i < rows; ++i) {
            A[(size_t)i * (size_t)cols + (size_t)k] = 0.0;
        }
    }

    // Back-substitution on R (upper triangle in A)
    x_out.assign((size_t)cols, 0.0);
    for (int i = cols - 1; i >= 0; --i) {
        const Real rii = A[(size_t)i * (size_t)cols + (size_t)i];
        if (std::abs(rii) < 1e-30) return false;

        Real sum = 0.0;
        for (int j = i + 1; j < cols; ++j) {
            sum += A[(size_t)i * (size_t)cols + (size_t)j] * x_out[(size_t)j];
        }
        x_out[(size_t)i] = (b[(size_t)i] - sum) / rii;
    }
    return true;
}

// Convenience: compute y = J * x for row-major J (m x n).
static inline void mat_vec_mul(int rows, int cols,
                               const std::vector<Real>& J,
                               const std::vector<Real>& x,
                               std::vector<Real>& y)
{
    y.assign((size_t)rows, 0.0);
    for (int i = 0; i < rows; ++i) {
        Real s = 0.0;
        const size_t off = (size_t)i * (size_t)cols;
        for (int j = 0; j < cols; ++j) s += J[off + (size_t)j] * x[(size_t)j];
        y[(size_t)i] = s;
    }
}

// ------------------------------------------------------------------------------
// Symmetric eigen-decomposition (Jacobi rotations)
// ------------------------------------------------------------------------------
// Dependency-free stabilized eigen-decomposition used by the dense SVD helpers.
// For the state-vector system, n is num = 2N + 10 (110 when N=50).
//
// Input:
//   A_in : symmetric matrix (n x n) in row-major.
// Output:
//   w    : eigenvalues (n)
//   V    : eigenvectors (n x n) in row-major, columns are eigenvectors.
//
// Notes:
// - Cyclic Jacobi with thresholding.
// - Eigenvalues are NOT sorted.
// ------------------------------------------------------------------------------
static bool jacobi_eigen_sym(const std::vector<Real>& A_in, int n,
                             std::vector<Real>& w,
                             std::vector<Real>& V,
                             int max_sweeps = 60)
{
    if ((int)A_in.size() != n * n) return false;

    std::vector<Real> A = A_in;
    w.assign((size_t)n, 0.0);
    V.assign((size_t)n * (size_t)n, 0.0);
    for (int i = 0; i < n; ++i) V[(size_t)i * (size_t)n + (size_t)i] = 1.0;

    auto a = [&](int r, int c) -> Real& { return A[(size_t)r * (size_t)n + (size_t)c]; };
    auto v = [&](int r, int c) -> Real& { return V[(size_t)r * (size_t)n + (size_t)c]; };

    const Real eps = std::numeric_limits<Real>::epsilon();

    for (int sweep = 0; sweep < max_sweeps; ++sweep) {
        Real off = 0.0;
        for (int p = 0; p < n; ++p) {
            for (int q = p + 1; q < n; ++q) off += std::abs(a(p, q));
        }
        if (off <= eps) break;

        for (int p = 0; p < n; ++p) {
            for (int q = p + 1; q < n; ++q) {
                const Real apq = a(p, q);
                if (std::abs(apq) <= eps) continue;

                const Real app = a(p, p);
                const Real aqq = a(q, q);
                const Real tau = (aqq - app) / (2.0 * apq);

                const Real t = std::copysign((Real)1.0, tau) /
                               (std::abs(tau) + std::sqrt((Real)1.0 + tau * tau));
                const Real c = 1.0 / std::sqrt(1.0 + t * t);
                const Real s = t * c;

                a(p, q) = 0.0;
                a(q, p) = 0.0;
                a(p, p) = app - t * apq;
                a(q, q) = aqq + t * apq;

                for (int k = 0; k < n; ++k) {
                    if (k == p || k == q) continue;
                    const Real aik = a(p, k);
                    const Real aqk = a(q, k);
                    a(p, k) = c * aik - s * aqk;
                    a(k, p) = a(p, k);
                    a(q, k) = s * aik + c * aqk;
                    a(k, q) = a(q, k);
                }

                for (int k = 0; k < n; ++k) {
                    const Real vip = v(k, p);
                    const Real viq = v(k, q);
                    v(k, p) = c * vip - s * viq;
                    v(k, q) = s * vip + c * viq;
                }
            }
        }
    }

    for (int i = 0; i < n; ++i) w[(size_t)i] = A[(size_t)i * (size_t)n + (size_t)i];
    return true;
}

static void jt_j(int m, int n, const std::vector<Real>& J, std::vector<Real>& A) {
    A.assign((size_t)n * (size_t)n, 0.0);
    for (int i = 0; i < m; ++i) {
        const size_t row = (size_t)i * (size_t)n;
        for (int a = 0; a < n; ++a) {
            const Real Jia = J[row + (size_t)a];
            if (Jia == 0.0) continue;
            for (int b = a; b < n; ++b) {
                A[(size_t)a * (size_t)n + (size_t)b] += Jia * J[row + (size_t)b];
            }
        }
    }
    for (int a = 0; a < n; ++a) {
        for (int b = a + 1; b < n; ++b) {
            A[(size_t)b * (size_t)n + (size_t)a] = A[(size_t)a * (size_t)n + (size_t)b];
        }
    }
}

static bool svd_via_jtj(int m, int n,
                        const std::vector<Real>& J,
                        std::vector<Real>& s,
                        std::vector<Real>& V)
{
    std::vector<Real> A;
    jt_j(m, n, J, A);

    std::vector<Real> w;
    if (!jacobi_eigen_sym(A, n, w, V, 60)) return false;

    for (Real& ev : w) {
        if (ev < 0.0 && ev > (Real)-1e-14) ev = 0.0;
    }

    std::vector<int> idx((size_t)n);
    for (int i = 0; i < n; ++i) idx[(size_t)i] = i;
    std::sort(idx.begin(), idx.end(), [&](int a, int b) { return w[(size_t)a] > w[(size_t)b]; });

    std::vector<Real> V_sorted((size_t)n * (size_t)n, 0.0);
    s.assign((size_t)n, 0.0);

    for (int col = 0; col < n; ++col) {
        const int src = idx[(size_t)col];
        const Real ev = std::max((Real)0.0, w[(size_t)src]);
        s[(size_t)col] = std::sqrt(ev);
        for (int r = 0; r < n; ++r) {
            V_sorted[(size_t)r * (size_t)n + (size_t)col] = V[(size_t)r * (size_t)n + (size_t)src];
        }
    }

    V.swap(V_sorted);
    return true;
}

// --------------------------------------------------------------------------------------
// One-sided Jacobi SVD for tall/skinny dense matrices (m >= n).
//
// Why this exists:
//   The residual system can become ill-conditioned for high finite-amplitude
//   waves, near-limiting waves and current cases.  Directly orthogonalizing the
//   dense Jacobian columns provides a stable singular-value basis for Newton
//   corrections without forming J^T J as the primary solve path.
//
// Algorithm:
//   We orthogonalize the columns of A = J in-place by applying Jacobi rotations
//   to column pairs (p, q) until off-diagonal correlations are negligible.
//   Accumulated rotations are stored in V (right singular vectors). At the end,
//   the singular values are the column norms of the orthogonalized matrix.
//
// Notes:
//   - U is not explicitly constructed.
//   - The Newton correction only requires singular values and right singular
//     vectors for the dense residual Jacobian.
//   - Output singular values are sorted descending, with V columns permuted.
// --------------------------------------------------------------------------------------
static bool svd_jacobi_onesided(int m, int n,
                                const std::vector<Real>& J,
                                std::vector<Real>& s,
                                std::vector<Real>& V)
{
    if (m < n || m <= 0 || n <= 0) return false;

    // Column-major working arrays make every Jacobi column operation contiguous.
    std::vector<Real> A((size_t)m * (size_t)n, 0.0);
    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < n; ++j) {
            A[(size_t)j * (size_t)m + (size_t)i] =
                J[(size_t)i * (size_t)n + (size_t)j];
        }
    }

    std::vector<Real> Vcol((size_t)n * (size_t)n, 0.0);
    for (int i = 0; i < n; ++i) Vcol[(size_t)i * (size_t)n + (size_t)i] = 1.0;

    auto col_dot = [&](int p, int q) -> Real {
        const Real* ap = A.data() + (size_t)p * (size_t)m;
        const Real* aq = A.data() + (size_t)q * (size_t)m;
        Real sum = 0.0;
#ifdef _OPENMP
#pragma omp simd reduction(+:sum)
#endif
        for (int i = 0; i < m; ++i) sum += ap[(size_t)i] * aq[(size_t)i];
        return sum;
    };

    auto col_norm2 = [&](int p) -> Real {
        const Real* ap = A.data() + (size_t)p * (size_t)m;
        Real sum = 0.0;
        Real compensation = 0.0;
        for (int i = 0; i < m; ++i) {
            const Real product = ap[(size_t)i] * ap[(size_t)i];
            const Real y = product - compensation;
            const Real t = sum + y;
            compensation = (t - sum) - y;
            sum = t;
        }
        return sum;
    };

    const Real eps = std::sqrt(std::numeric_limits<Real>::epsilon());
    const int max_sweeps = 30;
    std::vector<Real> norm2_cache((size_t)n, 0.0);
    for (int j = 0; j < n; ++j) norm2_cache[(size_t)j] = col_norm2(j);

    for (int sweep = 0; sweep < max_sweeps; ++sweep) {
        Real max_corr = 0.0;

        for (int p = 0; p < n; ++p) {
            const Real app = norm2_cache[(size_t)p];
            if (app <= 0.0) continue;

            for (int q = p + 1; q < n; ++q) {
                const Real aqq = norm2_cache[(size_t)q];
                if (aqq <= 0.0) continue;

                const Real apq = col_dot(p, q);
                const Real denom = std::sqrt(app * aqq);
                if (denom <= 0.0) continue;

                const Real corr = std::abs(apq) / denom;
                if (corr > max_corr) max_corr = corr;
                if (corr < 10.0 * eps) continue;

                const Real tau = (aqq - app) / (2.0 * apq);
                const Real t = std::copysign((Real)1.0, tau) /
                               (std::abs(tau) + std::sqrt((Real)1.0 + tau * tau));
                const Real c_rot = 1.0 / std::sqrt(1.0 + t * t);
                const Real s_rot = t * c_rot;

                Real* ap = A.data() + (size_t)p * (size_t)m;
                Real* aq = A.data() + (size_t)q * (size_t)m;
                Real norm_p = 0.0, comp_p = 0.0;
                Real norm_q = 0.0, comp_q = 0.0;
                for (int i = 0; i < m; ++i) {
                    const Real aip = ap[(size_t)i];
                    const Real aiq = aq[(size_t)i];
                    const Real new_p = c_rot * aip - s_rot * aiq;
                    const Real new_q = s_rot * aip + c_rot * aiq;
                    ap[(size_t)i] = new_p;
                    aq[(size_t)i] = new_q;

                    const Real prod_p = new_p * new_p;
                    const Real yp = prod_p - comp_p;
                    const Real tp = norm_p + yp;
                    comp_p = (tp - norm_p) - yp;
                    norm_p = tp;

                    const Real prod_q = new_q * new_q;
                    const Real yq = prod_q - comp_q;
                    const Real tq = norm_q + yq;
                    comp_q = (tq - norm_q) - yq;
                    norm_q = tq;
                }
                norm2_cache[(size_t)p] = norm_p;
                norm2_cache[(size_t)q] = norm_q;

                Real* vp = Vcol.data() + (size_t)p * (size_t)n;
                Real* vq = Vcol.data() + (size_t)q * (size_t)n;
#ifdef _OPENMP
#pragma omp simd
#endif
                for (int i = 0; i < n; ++i) {
                    const Real vip = vp[(size_t)i];
                    const Real viq = vq[(size_t)i];
                    vp[(size_t)i] = c_rot * vip - s_rot * viq;
                    vq[(size_t)i] = s_rot * vip + c_rot * viq;
                }
            }
        }

        if (max_corr < 10.0 * eps) break;
    }

    s.assign((size_t)n, 0.0);
    for (int j = 0; j < n; ++j) {
        s[(size_t)j] = std::sqrt(std::max((Real)0.0, norm2_cache[(size_t)j]));
    }

    std::vector<int> idx((size_t)n);
    for (int i = 0; i < n; ++i) idx[(size_t)i] = i;
    std::sort(idx.begin(), idx.end(), [&](int a, int b) {
        return s[(size_t)a] > s[(size_t)b];
    });

    std::vector<Real> s_sorted((size_t)n, 0.0);
    V.assign((size_t)n * (size_t)n, 0.0);
    for (int col = 0; col < n; ++col) {
        const int src = idx[(size_t)col];
        s_sorted[(size_t)col] = s[(size_t)src];
        const Real* source_col = Vcol.data() + (size_t)src * (size_t)n;
        for (int row = 0; row < n; ++row) {
            V[(size_t)row * (size_t)n + (size_t)col] = source_col[(size_t)row];
        }
    }
    s.swap(s_sorted);
    return true;
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
    int  N;          // Fourier order / half-wave interval count N (production N=50)

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

    Real breaking_limit_miche;
    Real breaking_index;
    bool is_breaking;

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
    Real MeanSquareBedVelocity;  // ub^2 [m^2/s^2], phase-mean orbital bed velocity square
    Real u_bed;
    Real u_surf;
    Real acc_max;
    Real w_max;
    Real asymmetry;

    // Convenience quantities for secondary engineering diagnostics
    Real tau_bed;
    Real ExcursionBed;

    bool converged;
    std::string last_error;

    // -------------------------- construction --------------------------------
    explicit FentonStreamFunction(Real H, Real T, Real depth, Real current = 0.0)
        : H_target(H), T_target(T), d(depth), Uc(current),
          g(Phys::G_STD), N(50),
          nstep(4), number(40), crit(1e-8), criter_final(1e-10),
          Current_criterion(1),
          MaxH(0.0), T_nd(0.0), Height(0.0), Current(0.0),
          k(0.0), L(0.0), c(0.0),
          eta_nodes((size_t)50 + 1, depth),
          Bj((size_t)50, 0.0),
          eta_crest(0.0), eta_trough(0.0),
          steepness(0.0), rel_depth(0.0), ursell(0.0), regime(""),
          breaking_limit_miche(0.0), breaking_index(0.0), is_breaking(false),
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
          u_bed(0.0), u_surf(0.0), acc_max(0.0), w_max(0.0), asymmetry(0.0),
          tau_bed(0.0), ExcursionBed(0.0),
          converged(false), last_error("")
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
        converged = false;
        last_error.clear();

        if (!(H_target > 0.0) || !(T_target > 0.0) || !(d > 0.0)) {
            last_error = "Invalid inputs: H, T, and d must be > 0.";
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

private:
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
    inline Real& SOL(int i, int k) { return sol[(size_t)i * 3 + (size_t)k]; }
    inline Real  SOL(int i, int k) const { return sol[(size_t)i * 3 + (size_t)k]; }

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
    // Linear/Fenton-McKee initial estimate used only to seed the nonlinear nonlinear residual solve.
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
            z[5] = Current * std::sqrt(z[2]); // ubar_1 sqrt(k/g), Eulerian current variable
            z[6] = 0.0;                       // ubar_2 sqrt(k/g), Stokes current variable
        } else {
            z[6] = Current * std::sqrt(z[2]);
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

                const Real e = std::exp(x);
                const Real inv_e = 1.0 / e;
                const Real sinhkd = 0.5 * (e - inv_e);
                const Real coshkd = 0.5 * (e + inv_e);
                const Real tanh_jkd = tanh_workspace[(size_t)j];

                const Real S = sinhkd + coshkd * tanh_jkd;
                const Real C = coshkd + sinhkd * tanh_jkd;

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
    static std::vector<Real> svd_solve(const std::vector<Real>& A, const std::vector<Real>& b, int n) {
        std::vector<Real> s, V;
        if (!LinAlg::svd_jacobi_onesided(n, n, A, s, V)) {
            throw std::runtime_error("SVD failure (jacobi_onesided).");
        }

        Real wmax = 0.0;
        for (Real si : s) wmax = std::max(wmax, si);
        const Real wmin = wmax * (Real)1e-12;

        // Atb = A^T b
        std::vector<Real> Atb((size_t)n, 0.0);
        for (int j = 0; j < n; ++j) {
            Real sum = 0.0;
            for (int i = 0; i < n; ++i) sum += A[(size_t)i * (size_t)n + (size_t)j] * b[(size_t)i];
            Atb[(size_t)j] = sum;
        }

        // tmp = V^T Atb
        std::vector<Real> tmp((size_t)n, 0.0);
        for (int i = 0; i < n; ++i) {
            Real sum = 0.0;
            for (int j = 0; j < n; ++j) sum += V[(size_t)j * (size_t)n + (size_t)i] * Atb[(size_t)j];
            tmp[(size_t)i] = sum;
        }

        // scale_i = tmp_i / s_i^2 (truncated)
        std::vector<Real> scale((size_t)n, 0.0);
        for (int i = 0; i < n; ++i) {
            const Real si = s[(size_t)i];
            if (si > wmin && si > 0.0) scale[(size_t)i] = tmp[(size_t)i] / (si * si);
            else scale[(size_t)i] = 0.0;
        }

        // x = V * scale
        std::vector<Real> x((size_t)n, 0.0);
        for (int j = 0; j < n; ++j) {
            Real sum = 0.0;
            for (int i = 0; i < n; ++i) sum += V[(size_t)j * (size_t)n + (size_t)i] * scale[(size_t)i];
            x[(size_t)j] = sum;
        }
        return x;
    }

    // ----------------------------------------------------------------------
    // Newton iteration for the nonlinear residual system F(z) = 0.
    // ----------------------------------------------------------------------
    Real newton_step(int /*iter_count*/) {
        const Real ss0 = eqns(rhs1);
        if (!std::isfinite(ss0)) {
            throw std::runtime_error("Non-finite residual norm at start of Newton step.");
        }

        const std::vector<Real> z0 = z;

        std::vector<Real> A((size_t)num * (size_t)num, 0.0);
        std::vector<Real> b((size_t)num, 0.0);
        for (int i = 1; i <= num; ++i) b[(size_t)(i - 1)] = -rhs1[(size_t)i];

        std::vector<std::string> jacobian_errors((size_t)num + 1);
        const int jacobian_threads = std::min(openmp_thread_count(), num);

        // Every Jacobian column is independent. Each OpenMP worker owns its
        // perturbed state and residual vector, so no solver state is shared or
        // modified while the dense finite-difference Jacobian is assembled.
#ifdef _OPENMP
#pragma omp parallel num_threads(jacobian_threads) if(num >= 24)
#endif
        {
            std::vector<Real> z_perturbed((size_t)num + 1, 0.0);
            std::vector<Real> rhs_column((size_t)num + 1, 0.0);
            std::vector<Real> tanh_workspace((size_t)n + 1, 0.0);

#ifdef _OPENMP
#pragma omp for schedule(static)
#endif
            for (int i = 1; i <= num; ++i) {
                try {
                    Real h = (Real)0.01 * z0[(size_t)i];
                    if (std::abs(z0[(size_t)i]) < (Real)1e-4) h = (Real)1e-5;
                    if (std::abs(h) > 1.0) h = std::copysign((Real)1.0, h);

                    std::copy(z0.begin(), z0.end(), z_perturbed.begin());
                    z_perturbed[(size_t)i] += h;
                    eqns_for_state(z_perturbed, rhs_column, tanh_workspace);

                    const Real inv_h = 1.0 / h;
                    for (int r = 1; r <= num; ++r) {
                        A[(size_t)(r - 1) * (size_t)num + (size_t)(i - 1)] =
                            (rhs_column[(size_t)r] - rhs1[(size_t)r]) * inv_h;
                    }
                } catch (const std::exception& e) {
                    jacobian_errors[(size_t)i] = e.what();
                } catch (...) {
                    jacobian_errors[(size_t)i] = "Unknown failure during Jacobian evaluation.";
                }
            }
        }

        for (int i = 1; i <= num; ++i) {
            if (!jacobian_errors[(size_t)i].empty()) {
                throw std::runtime_error(jacobian_errors[(size_t)i]);
            }
        }

        std::vector<Real> dx = svd_solve(A, b, num);
        for (Real v : dx) {
            if (!std::isfinite(v)) throw std::runtime_error("Non-finite Newton correction vector (dx).");
        }

        // Damped Newton step: preserve kd > 0 and accept only non-worsening residual norm.
        Real alpha = 1.0;
        Real ss_best = ss0;
        std::vector<Real> z_best = z0;

        while (alpha >= (Real)1e-4) {
            std::vector<Real> z_try = z0;
            for (int i = 1; i <= num; ++i) z_try[(size_t)i] = z0[(size_t)i] + alpha * dx[(size_t)(i - 1)];

            bool ok = (z_try[1] > 0.0);
            if (ok) {
                for (int i = 1; i <= num; ++i) {
                    if (!std::isfinite(z_try[(size_t)i])) { ok = false; break; }
                }
            }
            if (!ok) { alpha *= 0.5; continue; }

            z = z_try;
            const Real ss1 = eqns(rhs2);
            if (std::isfinite(ss1) && (ss1 <= ss_best)) {
                ss_best = ss1;
                z_best = z_try;
                if (ss1 <= ss0) break;
            }
            alpha *= 0.5;
        }

        z = z_best;

        // Convergence monitor: mean absolute correction of free-surface ordinates z[10..N+10].
        Real corr = 0.0;
        for (int i = 10; i <= n + 10; ++i) corr += std::abs(z_best[(size_t)i] - z0[(size_t)i]);
        corr /= (Real)(n + 1);

        return corr;
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

            const Real coshdelta = std::cosh((Real)j * Yloc);
            const Real sinhdelta = std::sinh((Real)j * Yloc);

            const Real C = coshdelta + sinhdelta * Tanh[(size_t)j];
            const Real S = sinhdelta + coshdelta * Tanh[(size_t)j];

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
    Real mean_square_bed_orbital_velocity(int nph = 720) {
        if (!(d > 0.0) || !(T_target > 0.0)) return 0.0;
        const int Nph = std::max(36, nph);

        Real ub2 = 0.0;
        for (int i = 0; i < Nph; ++i) {
            const Real ph = (2.0 * Phys::PI) * (Real)i / (Real)Nph;
            Real u_abs, w_abs, ax;
            get_kinematics(0.0, ph, u_abs, w_abs, ax); // bed-origin coordinate z = 0
            const Real u_orb = u_abs - Uc;
            ub2 += u_orb * u_orb;
        }
        return ub2 / (Real)Nph;
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

        const Real c_dimless     = z[4] / kd_sqrt;
        const Real ce_dimless    = z[5] / kd_sqrt;
        const Real cs_dimless    = z[6] / kd_sqrt;
        const Real ubar_dimless  = z[7] / kd_sqrt;

        const Real kd_32 = kd * kd_sqrt;

        const Real Q_dimless = ubar_dimless - z[8] / kd_32;
        const Real R_dimless = 1.0 + z[9] / kd;

        const Real pulse = z[8] + kd * z[5];
        const Real ke = 0.5 * (z[4] * pulse - z[5] * Q_dimless * kd_32);

        Real pe = 0.0;
        for (int i = 1; i <= n; ++i) pe += 0.25 * (Y[(size_t)i] * Y[(size_t)i]);

        const Real ub2_alg = 2.0 * z[9] - z[4] * z[4];
        const Real q_term = z[7] * kd - z[8];

        const Real sxx = 4.0 * ke - 3.0 * pe + ub2_alg * kd + 2.0 * z[5] * q_term;
        const Real f = z[4] * (3.0 * ke - 2.0 * pe) + 0.5 * ub2_alg * (pulse + z[4] * kd) + z[4] * z[5] * q_term;

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
        EulerianCurrent = Uc;
        StokesCurrent = MassTransport;
        MeanFluidSpeed = ubar_dimless * std::sqrt(g * d);

        VolumeFluxQ = Q_dimless * std::sqrt(g * std::pow(d, 3.0));
        WaveVolumeFlux_q = MeanFluidSpeed * d - VolumeFluxQ;

        Bernoulli_r = (R_dimless - 1.0) * g * d;

        KineticEnergy   = Phys::RHO * g * (d * d) * KE_depth; // [J/m^2]
        PotentialEnergy = Phys::RHO * g * (d * d) * PE_depth; // [J/m^2]

        MeanSquareBedVelocity = mean_square_bed_orbital_velocity(720);

        Cg = (std::abs(EnergyDensity) > 1e-12) ? (Power / EnergyDensity) : 0.0;
    }

    // ----------------------------------------------------------------------
    // Core continuation, Newton solution of F(z)=0, and dimensional post-processing.
    // ----------------------------------------------------------------------
    void solve_internal() {
        const Real dhe = Height / (Real)nstep;
        const Real dho = MaxH   / (Real)nstep;

        // Continue in wave-height parameter from near-linear wave to target H.
        for (int ns = 1; ns <= nstep; ++ns) {
            height = (Real)ns * dhe;
            Hoverd = (Real)ns * dho;

            // Initial or extrapolated state vector for the current continuation step.
            if (ns == 1) {
                init_linear();
            } else {
                for (int i = 1; i <= num; ++i) z[(size_t)i] = 2.0 * SOL(i, 2) - SOL(i, 1);

                // Fallback if extrapolation yields an invalid Fenton state vector.
                if (!(z[1] > 0.0)) {
                    for (int i = 1; i <= num; ++i) z[(size_t)i] = SOL(i, 2);
                }
                if (!(z[1] > 0.0)) throw std::runtime_error("Invalid extrapolated start state for continuation step.");
                for (int i = 1; i <= num; ++i) {
                    if (!std::isfinite(z[(size_t)i])) {
                        for (int k2 = 1; k2 <= num; ++k2) z[(size_t)k2] = SOL(k2, 2);
                        break;
                    }
                }
                for (int i = 1; i <= num; ++i) {
                    if (!std::isfinite(z[(size_t)i])) throw std::runtime_error("Invalid extrapolated start state for continuation step.");
                }
            }

            bool step_converged = false;

            for (int it = 1; it <= number; ++it) {
                Real err = 0.0;

                try {
                    err = newton_step(it);
                } catch (const std::exception&) {
                    // Retry once from the last converged continuation state if the Newton step fails.
                    if (ns > 1 && it == 1) {
                        for (int i = 1; i <= num; ++i) z[(size_t)i] = SOL(i, 2);
                        err = newton_step(it);
                    } else {
                        throw;
                    }
                }

                if (!std::isfinite(err)) throw std::runtime_error("Non-finite Newton correction.");

                // Update continuation storage before testing convergence, preserving the next extrapolation state.
                if (ns == 1) {
                    for (int i = 1; i <= num; ++i) SOL(i, 2) = z[(size_t)i];
                } else {
                    for (int i = 1; i <= num; ++i) SOL(i, 1) = SOL(i, 2);
                    for (int i = 1; i <= num; ++i) SOL(i, 2) = z[(size_t)i];
                }

                // Guard the dense Newton solve against diverging or non-finite residuals.
                if (!(z[1] > 0.0)) throw std::runtime_error("Divergence: invalid state vector encountered.");
                for (int i = 1; i <= num; ++i) {
                    if (!std::isfinite(z[(size_t)i])) throw std::runtime_error("Divergence: non-finite state vector encountered.");
                }

                const Real criter = (ns == nstep) ? criter_final : crit;
                if (it > 1 && err < criter * std::abs(z[1])) {
                    step_converged = true;
                    break;
                }
            }

            if (!step_converged) {
                throw std::runtime_error("Newton did not converge within iteration budget.");
            }

            // Reconstruct free-surface harmonics and Fourier coefficients for this step.
            compute_Y_and_B();
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

        breaking_limit_miche = 0.142 * L * std::tanh(k * d);
        breaking_index = (breaking_limit_miche > 0.0) ? (H_target / breaking_limit_miche) : 0.0;
        is_breaking = (breaking_limit_miche > 0.0) && (H_target > breaking_limit_miche);

        calc_integral_props_cpp();

        // Kinematic diagnostics derived from the converged Fourier stream-function field.
        {
            Real w, ax;
            get_kinematics(0.0, 0.0, u_bed, w, ax);
        }

        // Quadratic bed-shear estimate retained as an engineering diagnostic, not part of F(z).
        const Real cf_est = 0.005;
        tau_bed = 0.5 * Phys::RHO * cf_est * (u_bed * u_bed);
        ExcursionBed = std::abs(u_bed) * T_target / (2.0 * Phys::PI);

        // Crest and trough horizontal velocities at the free surface.
        Real w, ax, u_trough;
        get_kinematics(d + eta_crest, 0.0, u_surf, w, ax);
        get_kinematics(d + eta_trough, Phys::PI, u_trough, w, ax);
        asymmetry = (std::abs(u_trough) > 0.0) ? std::abs(u_surf / u_trough) : 0.0;

        // Scan phase over the free surface for maximum vertical velocity and horizontal acceleration.
        acc_max = 0.0;
        w_max   = 0.0;
        for (int i = 0; i < 40; ++i) {
            const Real X = (Real)i * Phys::PI / (Real)39; // 40-point half-wave scan, X in [0,pi]
            const Real kEta = surface_keta(X);
            const Real z_surf = d * (1.0 + kEta / kd);

            Real u, ww, a;
            get_kinematics(z_surf, X, u, ww, a);
            acc_max = std::max(acc_max, std::abs(a));
            w_max   = std::max(w_max,   std::abs(ww));
        }
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
static std::string generate_output(double H_in, double T_in, double d_in, double Uc_in) {
    using namespace ReportFmt;

    const bool has_current = (Uc_in != 0.0);
    FentonStreamFunction solver0(H_in, T_in, d_in, 0.0);
    FentonStreamFunction solverC(H_in, T_in, d_in, Uc_in);

    if (has_current) {
#ifdef _OPENMP
#pragma omp parallel sections num_threads(2) if(openmp_thread_count() >= 2)
#endif
        {
#ifdef _OPENMP
#pragma omp section
#endif
            { solver0.solve(); }
#ifdef _OPENMP
#pragma omp section
#endif
            { solverC.solve(); }
        }
    } else {
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

    if (!issues.empty()) {
        box_title(out, "Numerical failure / non-convergence.");
        for (const auto& msg : issues) box_text(out, msg);
        hline(out, '-');
        return out.str();
    }

    // ------------------------------- report header ------------------------------
    box_title(out, "NONLINEAR WAVE HYDRODYNAMICS SOLVER (FENTON)");
    box_text(out, std::string("Wave height (H)             : ") + py_str_float(H_in) + " m");
    box_text(out, std::string("Wave period (τ)             : ") + py_str_float(T_in) + " s");
    box_text(out, std::string("Water depth (d)             : ") + py_str_float(d_in) + " m");
    box_text(out, std::string("Eulerian current ū₁         : ") + py_str_float(Uc_in) + " m/s (positive with wave propagation)");
    hline(out, '-');
    box_text(out, "Status: Full nonlinear system solved successfully.");
    hline(out, '-');
    out << "\n";

    // ---------------------------- hydrodynamic summary ---------------------------
    const Real g = solver0.g;
    const Real d = solver0.d;
    const Real sqrt_gd = std::sqrt(g * d);
    const Real sqrt_g_over_d = std::sqrt(g / d);

    auto wc_num = [&](double v) -> Cell { return has_current ? Cell::numv(v) : Cell::strv("-"); };
    auto wc_str = [&](const std::string& s) -> Cell { return has_current ? Cell::strv(s) : Cell::strv("-"); };
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
    rows.push_back(Row{ false, "", { Cell::strv("Wave height (H)"), Cell::numv(solver0.H_target), wc_num(solverC.H_target), Cell::strv("m") } });
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
    rows.push_back(Row{ false, "", { Cell::strv("Group velocity (C𝗀 = F/E)ㅤ"), Cell::numv(solver0.Cg), wc_num(solverC.Cg), Cell::strv("m/s") } });

    rows.push_back(Row{ true, "KINEMATICS (EXTREMES / BED ORBITAL MOTION)", {} });
    rows.push_back(Row{ false, "", { Cell::strv("Max surface horiz. vel |u|"), Cell::numv(solver0.u_surf), wc_num(solverC.u_surf), Cell::strv("m/s") } });
    rows.push_back(Row{ false, "", { Cell::strv("Max bed horiz. vel |u|"), Cell::numv(solver0.u_bed), wc_num(solverC.u_bed), Cell::strv("m/s") } });
    rows.push_back(Row{ false, "", { Cell::strv("Max horiz. accel |aₓ|"), Cell::numv(solver0.acc_max), wc_num(solverC.acc_max), Cell::strv("m/s²") } });
    rows.push_back(Row{ false, "", { Cell::strv("Velocity asymmetry |uc|/|ut|"), Cell::numv(solver0.asymmetry), wc_num(solverC.asymmetry), Cell::strv("-") } });
    rows.push_back(Row{ false, "", { Cell::strv("Mean square bed orbital vel ub²"), Cell::numv(solver0.MeanSquareBedVelocity), wc_num(solverC.MeanSquareBedVelocity), Cell::strv("m²/s²") } });
    rows.push_back(Row{ false, "", { Cell::strv("Bed orbital RMS velocity ub,rms"), Cell::numv(std::sqrt(std::max(0.0, solver0.MeanSquareBedVelocity))), wc_num(std::sqrt(std::max(0.0, solverC.MeanSquareBedVelocity))), Cell::strv("m/s") } });

    rows.push_back(Row{ true, "NONLINEARITY / BREAKING DIAGNOSTICS", {} });
    const std::string warn0 = solver0.is_breaking ? "BREAKING" : "STABLE";
    const std::string warnC = has_current ? (solverC.is_breaking ? "BREAKING" : "STABLE") : "-";
    rows.push_back(Row{ false, "", { Cell::strv("Miche breaking limit (Hmax)"), Cell::numv(solver0.breaking_limit_miche), wc_num(solverC.breaking_limit_miche), Cell::strv("m") } });
    rows.push_back(Row{ false, "", { Cell::strv("Saturation (H/Hmax)"), Cell::numv(solver0.breaking_index), wc_num(solverC.breaking_index), Cell::strv("-") } });
    rows.push_back(Row{ false, "", { Cell::strv("Breaking status"), Cell::strv(warn0), wc_str(warnC), Cell::strv("-") } });
    rows.push_back(Row{ false, "", { Cell::strv("Ursell number (U)"), Cell::numv(solver0.ursell), wc_num(solverC.ursell), Cell::strv("-") } });
    rows.push_back(Row{ false, "", { Cell::strv("Regime (by d/L)"), Cell::strv(solver0.regime), wc_str(solverC.regime), Cell::strv("-") } });

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
        rr.push_back(Row{ false, "", { Cell::numv(17), Cell::strv("Mean square of bed velocity"),     Cell::strv(f5(slv.MeanSquareBedVelocity)),Cell::strv("ub²/gd"),                          Cell::strv(f5(slv.MeanSquareBedVelocity / (g_ * d_))) } });
        rr.push_back(Row{ false, "", { Cell::numv(18), Cell::strv("Radiation stress"),                Cell::strv(f5(kN(slv.Sxx))),           Cell::strv("S_xx/ρgd²"),                       Cell::strv(f5(slv.Sxx_depth)) } });
        rr.push_back(Row{ false, "", { Cell::numv(19), Cell::strv("Wave power"),                      Cell::strv(f5(kW(slv.Power))),         Cell::strv("F/(ρg³ᐟ²d⁵ᐟ²)ㅤㅤ"),               Cell::strv(f5(slv.F_depth)) } });

        print_table(out, title, h, cw, al, rr);
    };

    print_solution_flat(solver0, "SOLUTION.RES (NO CURRENT)");
    if (has_current) print_solution_flat(solverC, "SOLUTION.RES (WITH CURRENT)");

    // --------------------------------- glossary --------------------------------
    const std::vector<std::string> gh = { "TERM / SYMBOL", "MEANING", "UNITS / NONDIM" };
    const std::vector<int> gcw = { 14, 64, 19 };
    const std::vector<std::string> gal = { "left", "left", "left" };

    std::vector<Row> terms;
    auto add_term = [&](const std::string& a, const std::string& b, const std::string& c) {
        terms.push_back(Row{ false, "", { Cell::strv(a), Cell::strv(b), Cell::strv(c) } });
    };

    add_term("d", "Still-water depth (bed to mean water level). Reference length scale.", "m ; d/d=1");
    add_term("H", "Wave height (crest-to-trough).", "m ; H/d");
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
    add_term("Cg", "Group velocity, defined here as Cg = F/E.", "m/s");
    add_term("Sₓₓ", "Radiation stress component in wave direction.", "kN/m ; Sₓₓ/(ρgd²)");
    add_term("ub²", "Mean square *orbital* bed velocity: ub² = <(ub(t) − ū₁)²>. Non-negative by definition; computed by phase averaging.", "m²/s² ; /gd");
    add_term("ub,rms", "Root-mean-square orbital bed velocity: ub,rms = √(ub²).", "m/s");
    add_term("usurf,max", "Maximum horizontal velocity at free surface (scanned over phase).", "m/s");
    add_term("ubed,max", "Maximum horizontal velocity at seabed (scanned over phase).", "m/s");
    add_term("a_x,max", "Maximum horizontal acceleration magnitude.", "m/s²");
    add_term("Asymmetry", "Velocity asymmetry indicator |uc|/|ut|.", "-");
    add_term("Hmax", "Miche breaking limit used as stability diagnostic.", "m");
    add_term("Ursell", "Ursell number, a shallow-water nonlinearity measure.", "-");
    add_term("Regime", "Depth regime based on d/L (deep/intermediate/shallow).", "-");

    print_table(out, "PARAMETER DEFINITIONS & GLOSSARY", gh, gcw, gal, terms);

    return out.str();
}
// Convert LF to CRLF for Win32 multiline EDIT control display.
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

        make_label(L"Wave Height (m):", y); g_hEditH = make_edit(L"3.0", IDC_EDIT_H, y); y += 35;
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

            if (!(H > 0.0) || !(T > 0.0) || !(d > 0.0)) {
                MessageBoxW(hwnd, L"Please enter positive values for H, T, and d.",
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
#ifdef _OPENMP
    omp_set_dynamic(0);
    omp_set_max_active_levels(1);
#endif
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
