/* ==============================================================================
 *  ENGINEERING TECHNICAL REFERENCE & THEORETICAL FORMULATION
 * ==============================================================================
 *  PROGRAM:      Nonlinear Wave Hydrodynamics Solver (Fenton's Stream Function)
 *  METHOD:       Fourier Approximation Method for Steady Water Waves (N=50)
 *  REFERENCE:    Fenton, J.D. (1999). "Numerical methods for nonlinear waves."
 *                In P.L.-F. Liu (Ed.), Advances in Coastal and Ocean Engineering
 *                (Vol. 5, pp. 241–324). World Scientific: Singapore.
 * ==============================================================================
 *
 *  1. INTRODUCTION & SCOPE
 *  -----------------------------------------------------------------------------
 *  This software calculates the hydrodynamics of steady, periodic surface gravity
 *  waves using high-order Stream Function theory. Unlike Linear (Airy) Theory,
 *  which assumes infinitesimal amplitudes, this method retains full nonlinearity
 *  in the boundary conditions.
 *
 *  Implementation Specifics (C++ Port):
 *  - Solver:    Dense nonlinear least-squares solver aligned with SciPy
 *               scipy.optimize.least_squares behavior:
 *               • TRF-style trust-region steps for robust continuation seeding.
 *               • MINPACK-style Levenberg–Marquardt for final rapid convergence.
 *  - Stability: Uses a Homotopy (continuation) method, stepping wave height
 *               incrementally from near-linear to target height to guarantee convergence.
 *               (Implemented as 4 steps: linspace(0.01, H_target, 4), matching Python.)
 *  - Regime:    Applicable to stable waves in shallow, intermediate, and deep
 *               water regimes up to the Miche breaking limit.
 *
 *  2. GOVERNING FIELD EQUATIONS
 *  -----------------------------------------------------------------------------
 *  The fluid is modeled as inviscid, incompressible, and irrotational.
 *  The flow is solved in a frame of reference moving with the wave celerity (c),
 *  rendering the flow steady.
 *
 *  A. Field Equation (Laplace):
 *     ∇²ψ = ∂²ψ/∂x² + ∂²ψ/∂z² = 0
 *     Where ψ(x,z) is the stream function. Velocities are defined as:
 *     u =  ∂ψ/∂z   (Horizontal)
 *     w = -∂ψ/∂x   (Vertical)
 *
 *  B. Bottom Boundary Condition (BBC) at z=0:
 *     The seabed is impermeable (a streamline).
 *     ψ(x, 0) = -Q
 *     Where Q is the volume flux per unit width in the moving frame.
 *
 *  3. FREE SURFACE BOUNDARY CONDITIONS
 *  -----------------------------------------------------------------------------
 *  The solution is constrained by two nonlinear conditions at the unknown
 *  free surface elevation z = η(x):
 *
 *  A. Kinematic Boundary Condition (KBC):
 *     The free surface is a streamline (constant ψ).
 *     ψ(x, η) = 0
 *
 *  B. Dynamic Boundary Condition (DBC - Bernoulli):
 *     Pressure is constant (atmospheric) along the surface.
 *     1/2 * [ (∂ψ/∂x)² + (∂ψ/∂z)² ] + gη = R
 *     Where R is the Bernoulli constant (Total Energy Head).
 *
 *  4. NUMERICAL SOLUTION (FOURIER ANSATZ)
 *  -----------------------------------------------------------------------------
 *  The stream function is approximated by a truncated Fourier series of order N
 *  (N=50) that analytically satisfies the Field Equation and Bottom BC:
 *
 *    ψ(x,z) = -(ū + c) z + Σ_{j=1..N} B_j * [sinh(jkz)/cosh(jkd)] * cos(jkx)
 *
 *  Deep Water Numerical Stability:
 *  To prevent floating-point overflow when kd >> 1, the code replaces the
 *  hyperbolic ratio with asymptotic exponentials when arguments > 25.0:
 *
 *    sinh(jkz)/cosh(jkd) ≈ exp(jk(z-d))
 *
 *  Optimization Vector (State Space):
 *  The solver minimizes residuals for the vector
 *
 *    X = [k, η_0...η_N, B_1...B_N, Q, R]
 *
 *  IMPORTANT (Overdetermined Residual System):
 *  The residual vector dimension is NOT equal to the number of unknowns.
 *  For N=50:
 *    - Unknowns: n = 1 + (N+1) + N + 2 = 2N + 4 = 104
 *    - Residuals: m = 3 + (N+1) + (N+1) = 2(N+1) + 3 = 105
 *  The implementation MUST allocate m=105 and never assume m==n, otherwise
 *  out-of-bounds writes will occur and results will become optimizer/flags dependent.
 *
 *  5. DERIVED PHYSICAL PARAMETERS & OUTPUT DEFINITIONS
 *  -----------------------------------------------------------------------------
 *  Upon convergence, the software calculates the following engineering parameters
 *  derived from the solved Fourier coefficients (B_j).
 *
 *  A. FUNDAMENTAL WAVE GEOMETRY & PHASE
 *  ------------------------------------
 *  1. Wavelength (L):
 *     Horizontal distance between crests. Solved via dispersion relation.
 *     L = c·T = 2π / k
 *
 *  2. Celerity (c):
 *     Phase velocity. c = L / T.
 *
 *  B. KINEMATICS (VELOCITIES & ACCELERATIONS)
 *  ------------------------------------------
 *  1. Horizontal Velocity (u):
 *     u(x,z) = c - ū + Σ_{j=1..N} jkB_j * [cosh(jkz)/cosh(jkd)] * cos(jkx)
 *
 *  2. Vertical Velocity (w):
 *     w(x,z) = Σ_{j=1..N} jkB_j * [sinh(jkz)/cosh(jkd)] * sin(jkx)
 *
 *  3. Max Acceleration (a_x):
 *     Total derivative (Convective acceleration).
 *     a_x = Du/Dt = u * ∂u/∂x + w * ∂u/∂z
 *
 *     NOTE (Python-parity detail):
 *     The vertical perturbation term uses +sin(j·phase) (not -sin). A sign error
 *     here distorts the convective term w·∂u/∂z and breaks Max Accel parity.
 *
 *  4. Velocity Asymmetry:
 *     Asymmetry = |u_crest| / |u_trough|
 *
 *  C. DYNAMICS (INTEGRAL PROPERTIES)
 *  ---------------------------------
 *  Computed using exact integral invariants (Fenton Eqs 14-16).
 *
 *  1. Impulse (I):
 *     Total wave momentum (kg·m/s).
 *     I = ρ(c d - Q)
 *
 *  2. Energy Density (E):
 *     Mean Energy (J/m²).
 *     PE = 1/2 ρ g mean(η²)
 *     KE = 1/2 (cI - Qρ U_c)
 *     E  = PE + KE
 *
 *  3. Power / Energy Flux (P):
 *     Rate of energy transfer (W/m).
 *     P = c(3KE - 2PE) + 1/2 mean(u_b²)(I + ρ c d) + 1/2 ρ Q U_c²
 *
 *     Note on mean(u_b²) (Mean Square Bed Velocity):
 *     To avoid deep-water integration errors, this is computed algebraically:
 *     mean(u_b²) = 2(R - g d) - c²
 *
 *  4. Radiation Stress (Sxx):
 *     Excess momentum flux (N/m).
 *     Sxx = 4KE - 3PE + ρ mean(u_b²) d + 2ρ I U_c
 *
 *  5. Mean Stokes Drift (U_drift):
 *     U_drift = I / (ρ d)
 *
 *  D. STABILITY & REGIME CLASSIFICATION
 *  ------------------------------------
 *  1. Ursell Number (U_r):
 *     U_r = H L² / d³ (Values > 26 indicate significant nonlinearity).
 *
 *  2. Miche Limit (H_max):
 *     Theoretical max height before breaking.
 *     H_max = 0.142 L tanh(kd)
 *
 *  3. Saturation (Breaking Index):
 *     Saturation = H / H_max
 *     - If > 1.0: Wave is BREAKING.
 *     - If < 1.0: Wave is STABLE.
 *
 *  4. Regime:
 *     - Shallow:      d/L < 0.05
 *     - Intermediate: 0.05 < d/L < 0.5
 *     - Deep:         d/L > 0.5
 *
 * ==============================================================================
 *  6. SOFTWARE USAGE & COMPILATION GUIDE (C++)
 * ==============================================================================
 *
 *  A. PREREQUISITES
 *  ----------------
 *  - Windows (Win32 API GUI build). MinGW-w64 (g++) or MSVC with C++17 support.
 *  - Basic familiarity with the command line (Terminal/CMD).
 *
 *  B. BUILDING FROM SOURCE
 *  ------------------------------------
 *  1. Ensure a C++17-capable toolchain is installed and on PATH (g++ recommended).
 *
 *  2. Compile (release-like build, GUI subsystem):
 *     > g++ fenton_gui.cpp -o fenton_gui.exe -O3 -std=gnu++17 -march=native ^
 *       -lgdi32 -luser32 -lkernel32 -lcomctl32 -static -mwindows -pthread
 *
 *     Build note (Python-parity / "no current" stability):
 *     This solver is path-sensitive: small floating-point differences can change trust-region
 *     acceptance and LM step damping, pushing the "no current" case to a different local minimum.
 *     On this system, Python-parity (including the "no current" scenario) is achieved ONLY when
 *     compiling with -march=native (enables CPU-specific instruction selection such as FMA and
 *     vectorization patterns). Removing -march=native has been observed to produce incorrect
 *     "no current" results even when all equations and tolerances are unchanged.
 *     For bitwise-stable results, keep the same CPU family, compiler version, and flags.
 *
 *  3. Run:
 *     > fenton_gui.exe
 *
 *  C. STANDALONE EXECUTABLE (.EXE)
 *  ------------------------------
 *  The build command above produces a standalone fenton_gui.exe (no Python required).
 *  The program writes results to the GUI output panel and also produces an "output.txt"
 *  file (UTF-8) to disk.
 *
 * ==============================================================================
 *  BIBLIOGRAPHY
 * ==============================================================================
 *
 *  1.  Fenton, J.D. (1999). "Numerical methods for nonlinear waves."
 *      In P.L.-F. Liu (Ed.), Advances in Coastal and Ocean Engineering (Vol. 5,
 *      pp. 241–324). World Scientific: Singapore.
 *      [Primary Source: Comprehensive review of fully-nonlinear methods including
 *      Fourier approximation, Boundary Integral Equation (BIE) methods, and
 *      Local Polynomial Approximation].
 *      URL: https://johndfenton.com/Papers/Fenton99Liu-Numerical-methods-for-nonlinear-waves.pdf
 *
 *  2.  Fenton, J.D. (1988). "The numerical solution of steady water wave problems."
 *      Computers & Geosciences, 14(3), 357–368.
 *      [The core algorithm for high-accuracy Stream Function Theory].
 *      URL: https://doi.org/10.1016/0098-3004(88)90066-0
 *
 *  3.  Fenton, J.D. (1985). "A fifth-order Stokes theory for steady waves."
 *      Journal of Waterway, Port, Coastal, and Ocean Engineering, 111(2), 216–234.
 *      [Standard analytical theory for deep/intermediate water pile design].
 *      URL: https://doi.org/10.1061/(ASCE)0733-950X(1985)111:2(216)
 *
 *  4.  Fenton, J.D. (1978). "Wave forces on vertical bodies of revolution."
 *      Journal of Fluid Mechanics, 85(2), 241–255.
 *      [Foundational diffraction theory for large diameter piles].
 *      URL: https://johndfenton.com/Papers/Fenton78-Waves-on-bodies-of-revolution.pdf
 *
 *  5.  Fenton, J.D. (1990). "Nonlinear wave theories." In B. Le Méhauté &
 *      D.M. Hanes (Eds.), The Sea: Ocean Engineering Science (Vol. 9, Part A).
 *      John Wiley & Sons.
 *      [Comprehensive guide for selecting wave theories: Stokes vs Cnoidal vs Stream].
 *      URL: https://www.johndfenton.com/Papers/Fenton90b-Nonlinear-wave-theories.pdf
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

// Continue with the rest of C/C++ includes
#include <algorithm>
#include <atomic>
#include <cmath>
#include <condition_variable>
#include <cstdint>
#include <functional>
#include <future>
#include <iomanip>
#include <limits>
#include <memory>
#include <mutex>
#include <queue>
#include <sstream>
#include <stdexcept>
#include <string>
#include <thread>
#include <utility>
#include <vector>
#include <fstream>

// ==============================================================================
//  THREAD POOL (Persistent, low overhead)
// ==============================================================================

class ThreadPool {
public:
    explicit ThreadPool(size_t threads = 0) : m_stop(false) {
        if (threads == 0) threads = std::thread::hardware_concurrency();
        if (threads == 0) threads = 2;

        m_workers.reserve(threads);
        for (size_t i = 0; i < threads; ++i) {
            m_workers.emplace_back([this]() {
                for (;;) {
                    std::function<void()> task;
                    {
                        std::unique_lock<std::mutex> lock(m_mutex);
                        m_cv.wait(lock, [this]() { return m_stop || !m_tasks.empty(); });
                        if (m_stop && m_tasks.empty()) return;
                        task = std::move(m_tasks.front());
                        m_tasks.pop();
                    }
                    task();
                }
            });
        }
    }

    template <class F, class... Args>
    auto enqueue(F&& f, Args&&... args)
        -> std::future<typename std::result_of<F(Args...)>::type>
    {
        using return_type = typename std::result_of<F(Args...)>::type;

        auto task = std::make_shared<std::packaged_task<return_type()>>(
            std::bind(std::forward<F>(f), std::forward<Args>(args)...)
        );

        std::future<return_type> res = task->get_future();
        {
            std::unique_lock<std::mutex> lock(m_mutex);
            if (m_stop) throw std::runtime_error("enqueue on stopped ThreadPool");
            m_tasks.emplace([task]() { (*task)(); });
        }
        m_cv.notify_one();
        return res;
    }

    ~ThreadPool() {
        {
            std::unique_lock<std::mutex> lock(m_mutex);
            m_stop = true;
        }
        m_cv.notify_all();
        for (std::thread& w : m_workers) {
            if (w.joinable()) w.join();
        }
    }

private:
    std::vector<std::thread> m_workers;
    std::queue<std::function<void()>> m_tasks;
    std::mutex m_mutex;
    std::condition_variable m_cv;
    bool m_stop;
};

// Lazy global pool (used for Jacobian construction)
static std::unique_ptr<ThreadPool> g_pool;

// ==============================================================================
//  PHYSICAL CONSTANTS (must match fenton_gui.py)
// ==============================================================================

namespace Phys {
    using Real = double;
    constexpr Real PI    = 3.141592653589793238462643383279502884;
    constexpr Real RHO   = 1025.0;
    constexpr Real G_STD = 9.80665;
}
using Phys::Real;

// ==============================================================================
//  SMALL LINEAR ALGEBRA (dense, for n ~ 100)
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
// Implemented dependency-free to support SciPy-parity TRF(tr_solver='exact')
// steps. For our dense n~100 systems, Jacobi is sufficiently robust.
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
//   SciPy's TRF('exact') path performs an SVD of J_h directly (LAPACK), which is
//   numerically more stable than forming J^T J (squares the condition number).
//   The earlier C++ parity implementation used eigen(J^T J) for convenience.
//   That is usually fine, but if you want to be extreme about accuracy
//   (especially near-breaking regimes, strong nonlinearity, or ill-conditioned
//   Jacobians), one-sided Jacobi SVD is a better match to SciPy's intent.
//
// Algorithm:
//   We orthogonalize the columns of A = J in-place by applying Jacobi rotations
//   to column pairs (p, q) until off-diagonal correlations are negligible.
//   Accumulated rotations are stored in V (right singular vectors). At the end,
//   the singular values are the column norms of the orthogonalized matrix.
//
// Notes:
//   - We do NOT explicitly construct U.
//   - TRF only needs (s, V) and uf = U^T f. We compute uf via
//       uf = diag(1/s) * V^T * (J^T f)
//     which remains valid for this decomposition.
//   - Output singular values are sorted descending, with V columns permuted.
// --------------------------------------------------------------------------------------
static bool svd_jacobi_onesided(int m, int n,
                                const std::vector<Real>& J,
                                std::vector<Real>& s,
                                std::vector<Real>& V)
{
    if (m < n || m <= 0 || n <= 0) return false;

    // Working copy of J (A in the algorithm).
    std::vector<Real> A = J; // row-major (m x n)

    // V starts as identity.
    V.assign((size_t)n * (size_t)n, 0.0);
    for (int i = 0; i < n; ++i) V[(size_t)i * (size_t)n + (size_t)i] = 1.0;

    auto col_dot = [&](int p, int q) -> Real {
        Real sum = 0.0;
        // Deterministic accumulation order; use compensated summation to reduce
        // round-off when m is large or columns are nearly parallel.
        Real c = 0.0;
        for (int i = 0; i < m; ++i) {
            const Real prod = A[(size_t)i * (size_t)n + (size_t)p] * A[(size_t)i * (size_t)n + (size_t)q];
            const Real y = prod - c;
            const Real t = sum + y;
            c = (t - sum) - y;
            sum = t;
        }
        return sum;
    };

    auto col_norm2 = [&](int p) -> Real {
        Real sum = 0.0;
        Real c = 0.0;
        for (int i = 0; i < m; ++i) {
            const Real v = A[(size_t)i * (size_t)n + (size_t)p];
            const Real prod = v * v;
            const Real y = prod - c;
            const Real t = sum + y;
            c = (t - sum) - y;
            sum = t;
        }
        return sum;
    };

    // Tolerance consistent with Jacobi sweeps.
    const Real eps = std::sqrt(std::numeric_limits<Real>::epsilon());

    // Sweeps: for n~104, 20-30 sweeps is more than enough for near-orthogonality.
    const int max_sweeps = 30;
    for (int sweep = 0; sweep < max_sweeps; ++sweep) {
        Real max_corr = 0.0;

        for (int p = 0; p < n; ++p) {
            const Real app = col_norm2(p);
            if (app <= 0.0) continue;
            for (int q = p + 1; q < n; ++q) {
                const Real aqq = col_norm2(q);
                if (aqq <= 0.0) continue;

                const Real apq = col_dot(p, q);
                const Real denom = std::sqrt(app * aqq);
                if (denom <= 0.0) continue;

                const Real corr = std::abs(apq) / denom;
                if (corr > max_corr) max_corr = corr;

                // If columns are nearly orthogonal, skip.
                if (corr < 10.0 * eps) continue;

                // Compute Jacobi rotation for columns p, q.
                const Real tau = (aqq - app) / (2.0 * apq);
                const Real t = std::copysign((Real)1.0, tau) /
                               (std::abs(tau) + std::sqrt((Real)1.0 + tau * tau));
                const Real c_rot = 1.0 / std::sqrt(1.0 + t * t);
                const Real s_rot = t * c_rot;

                // Rotate columns of A.
                for (int i = 0; i < m; ++i) {
                    const size_t off = (size_t)i * (size_t)n;
                    const Real aip = A[off + (size_t)p];
                    const Real aiq = A[off + (size_t)q];
                    A[off + (size_t)p] = c_rot * aip - s_rot * aiq;
                    A[off + (size_t)q] = s_rot * aip + c_rot * aiq;
                }

                // Accumulate rotation into V.
                for (int i = 0; i < n; ++i) {
                    const size_t off = (size_t)i * (size_t)n;
                    const Real vip = V[off + (size_t)p];
                    const Real viq = V[off + (size_t)q];
                    V[off + (size_t)p] = c_rot * vip - s_rot * viq;
                    V[off + (size_t)q] = s_rot * vip + c_rot * viq;
                }
            }
        }

        // Stop early if correlations are negligible.
        if (max_corr < 10.0 * eps) break;
    }

    // Singular values are norms of orthogonalized columns.
    s.assign((size_t)n, 0.0);
    for (int j = 0; j < n; ++j) {
        const Real n2 = col_norm2(j);
        s[(size_t)j] = std::sqrt(std::max((Real)0.0, n2));
    }

    // Sort descending singular values (like LAPACK/SciPy).
    std::vector<int> idx((size_t)n);
    for (int i = 0; i < n; ++i) idx[(size_t)i] = i;
    std::sort(idx.begin(), idx.end(), [&](int a, int b) { return s[(size_t)a] > s[(size_t)b]; });

    std::vector<Real> V_sorted((size_t)n * (size_t)n, 0.0);
    std::vector<Real> s_sorted((size_t)n, 0.0);
    for (int col = 0; col < n; ++col) {
        const int src = idx[(size_t)col];
        s_sorted[(size_t)col] = s[(size_t)src];
        for (int r = 0; r < n; ++r) {
            V_sorted[(size_t)r * (size_t)n + (size_t)col] = V[(size_t)r * (size_t)n + (size_t)src];
        }
    }
    V.swap(V_sorted);
    s.swap(s_sorted);

    return true;
}

} // namespace LinAlg

// ==============================================================================
//  FENTON STREAM FUNCTION SOLVER (Matches Python formulation)
// ==============================================================================


class FentonStreamFunction {
public:
    // ----------------------------- inputs -----------------------------------
    Real H_target;   // [m]
    Real T_target;   // [s]
    Real d;          // [m]
    Real Uc;         // [m/s] (Eulerian / lab-frame)
    Real g;          // [m/s^2]
    int  N;          // Fourier order (N=50)

    // ------------------------- solver controls ------------------------------
    int  nstep;          // continuation steps in wave height
    int  number;         // max Newton iterations per step
    Real crit;           // intermediate-step convergence factor
    Real criter_final;   // final-step convergence factor

    // Current criterion (1=Eulerian, 2=Stokes). GUI uses Eulerian current.
    int Current_criterion;

    // Derived nondimensional inputs (C++ Read_data equivalents)
    Real MaxH;     // H/d
    Real T_nd;     // T * sqrt(g/d)
    Real Height;   // (H/d) / (T_nd^2) = H/(g T^2)
    Real Current;  // Uc / sqrt(g d)

    // ----------------------------- outputs ----------------------------------
    Real k;        // [rad/m]
    Real L;        // [m]
    Real c;        // [m/s]

    std::vector<Real> eta_nodes; // size (N+1), absolute z from bed [m]
    std::vector<Real> Bj;        // size (N), B_1..B_N (depth scaling)

    Real eta_crest;    // [m] relative to SWL
    Real eta_trough;   // [m] relative to SWL

    Real steepness;
    Real rel_depth;
    Real ursell;
    std::string regime;

    Real breaking_limit_miche;
    Real breaking_index;
    bool is_breaking;

    // Integral properties (dimensional)
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

    // Additional invariants (depth-scaled; used in Solution-Flat reporting)
    Real E_depth;
    Real KE_depth;
    Real PE_depth;

    // Mean Stokes / mass-transport current (dimensional)
    Real MassTransport;

    // Bed orbital statistics / kinematics
    Real MeanSquareBedVelocity;  // ub^2 [m^2/s^2]
    Real u_bed;
    Real u_surf;
    Real acc_max;
    Real w_max;
    Real asymmetry;

    // Convenience (not printed in the report header but computed like Python)
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
        // Inputs screening (keep parity philosophy: fail-fast on invalid physical inputs)
        if (d > 0.0) {
            MaxH   = H_target / d;
            T_nd   = T_target * std::sqrt(g / d);
            Height = (T_nd > 0.0) ? (MaxH / (T_nd * T_nd)) : 0.0;
            Current = Uc / std::sqrt(g * d);
        }

        // Robustness for large ambient currents: increase continuation/iteration budget.
        if (std::abs(Current) >= 1.0) {
            nstep  = std::max(nstep, 8);
            number = std::max(number, 80);
        }

        // Internal arrays (C++-parity with Python's z[1..num] layout)
        n = N;
        num = 2 * n + 10;

        z.assign((size_t)num + 1, 0.0);
        rhs1.assign((size_t)num + 1, 0.0);
        rhs2.assign((size_t)num + 1, 0.0);

        coeff.assign((size_t)n + 1, 0.0);
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

    // GUI-facing kinematics: (u_abs, w_abs, a_x) at z from bed [m] and phase X=kx [rad].
    void get_kinematics(Real z_bed, Real phase, Real& u_abs, Real& w_abs, Real& ax) const {
        const Real kd = z[1];
        if (!(kd > 0.0) || !(d > 0.0)) {
            u_abs = w_abs = ax = 0.0;
            return;
        }

        const Real k_phys = kd / d;
        const Real X = phase;
        const Real Yloc = k_phys * (z_bed - d); // wave scaling: k(z-d)

        Real u_nd = 0.0, v_nd = 0.0, dudt_nd = 0.0;
        point(X, Yloc, u_nd, v_nd, dudt_nd);

        const Real scale_v = std::sqrt(g * d);
        u_abs = u_nd * scale_v;
        w_abs = v_nd * scale_v;
        ax    = dudt_nd * g;
    }

    // Accessor used for reporting (dimensionless wavenumber kd)
    Real kd_dimless() const { return (z.size() > 1) ? z[1] : 0.0; }

private:
    int n = 0;
    int num = 0;

    // 1-based vectors (index 0 unused)
    std::vector<Real> z, rhs1, rhs2, coeff, Tanh, B, Y;

    // trig tables
    std::vector<Real> cosa, sina;
    std::vector<Real> cos_nm, sin_nm; // (n+1) x n, row-major: [m*n + (j-1)]

    // continuation storage sol[i][1..2] flattened: sol[(i*3)+k]
    std::vector<Real> sol;

    // step variables
    Real height = 0.0; // stepped dimensionless height
    Real Hoverd = 0.0; // stepped H/d

private:
    inline Real& SOL(int i, int k) { return sol[(size_t)i * 3 + (size_t)k]; }
    inline Real  SOL(int i, int k) const { return sol[(size_t)i * 3 + (size_t)k]; }

    void init_trig_tables() {
        // cosa[k] = cos(k*pi/n), k=0..2n
        for (int i = 0; i <= 2 * n; ++i) {
            const Real ang = (Real)i * Phys::PI / (Real)n;
            cosa[(size_t)i] = std::cos(ang);
            sina[(size_t)i] = std::sin(ang);
        }

        // cos_nm[m,j] = cos((m*j)%2n * pi/n), for m=0..n, j=1..n
        for (int m = 0; m <= n; ++m) {
            for (int j = 1; j <= n; ++j) {
                const int nm = (m * j) % (2 * n);
                cos_nm[(size_t)m * (size_t)n + (size_t)(j - 1)] = cosa[(size_t)nm];
                sin_nm[(size_t)m * (size_t)n + (size_t)(j - 1)] = sina[(size_t)nm];
            }
        }
    }

    // ----------------------------------------------------------------------
    // Port of Python _init_linear()
    // ----------------------------------------------------------------------
    void init_linear() {
        const Real pi = Phys::PI;

        const Real sigma = (Hoverd > 0.0) ? (2.0 * pi * std::sqrt(height / Hoverd)) : 0.0;

        if (sigma > 0.0) {
            const Real t = std::tanh(std::pow(sigma, 1.5));
            // Fenton & McKee (1990) approximation (as used in the Python port)
            z[1] = (sigma * sigma) / std::pow(t, (Real)(2.0 / 3.0));
        } else {
            z[1] = 2.0 * pi * std::max(height, (Real)1e-12) / std::max(Hoverd, (Real)1e-12);
        }

        z[2] = z[1] * Hoverd;
        z[4] = std::sqrt(std::tanh(z[1]));
        z[3] = 2.0 * pi / z[4];

        // Current initialisation (finite depth)
        if (Current_criterion == 1) {
            z[5] = Current * std::sqrt(z[2]); // ce
            z[6] = 0.0;                       // cs
        } else {
            z[6] = Current * std::sqrt(z[2]);
            z[5] = 0.0;
        }

        z[7] = z[4];       // ubar (dimensionless)
        z[8] = 0.0;        // q-term
        z[9] = 0.5 * z[7] * z[7];

        z[10] = 0.5 * z[2];
        for (int i = 1; i <= n; ++i) {
            z[n + i + 10] = 0.0;                 // coeff
            z[i + 10] = 0.5 * z[2] * cosa[(size_t)i];
        }
        z[n + 11] = 0.5 * z[2] / z[7];

        // store sol[] for extrapolation
        for (int i = 1; i < 10; ++i) SOL(i, 1) = z[(size_t)i];
        for (int i = 10; i <= num; ++i) SOL(i, 1) = 0.0;
    }

    // ----------------------------------------------------------------------
    // Port of Python _eqns(): fills rhs_out[1..num] and returns sum(rhs^2)
    // ----------------------------------------------------------------------
    Real eqns(std::vector<Real>& rhs_out) {
        const Real pi = Phys::PI;
        rhs_out.assign((size_t)num + 1, 0.0);

        // Eqn 1
        rhs_out[1] = z[2] - z[1] * Hoverd;
        // Eqn 2 (Period case)
        rhs_out[2] = z[2] - height * z[3] * z[3];
        // Eqn 3
        rhs_out[3] = z[4] * z[3] - 2.0 * pi;
        // Eqn 4
        rhs_out[4] = z[5] + z[7] - z[4];
        // Eqn 5
        rhs_out[5] = z[1] * (z[6] + z[7] - z[4]) - z[8];

        // coeff and tanh tables
        for (int i = 1; i <= n; ++i) {
            coeff[(size_t)i] = z[(size_t)(n + i + 10)];
            Tanh[(size_t)i]  = std::tanh((Real)i * z[1]);
        }

        // Eqn 6 (finite depth; correction uses sqrt(z[1]))
        rhs_out[6] = z[(size_t)(Current_criterion + 4)] - Current * std::sqrt(z[1]);

        // Eqn 7 (mean free surface level; scaling constant irrelevant)
        rhs_out[7] = z[10] + z[(size_t)(n + 10)];
        for (int i = 1; i < n; ++i) rhs_out[7] += 2.0 * z[(size_t)(10 + i)];

        // Eqn 8 (wave height definition)
        rhs_out[8] = z[10] - z[(size_t)(n + 10)] - z[2];

        // Eqns 9.. and n+10.. : free surface BCs at nodes m=0..n
        for (int m = 0; m <= n; ++m) {
            const Real zsurf = z[(size_t)(10 + m)]; // k(eta-d)

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

                const Real tanh_jkd = Tanh[(size_t)j];

                const Real S = sinhkd + coshkd * tanh_jkd;
                const Real C = coshkd + sinhkd * tanh_jkd;

                const Real c_nm = cosrow[(size_t)(j - 1)];
                const Real s_nm = sinrow[(size_t)(j - 1)];

                const Real cj = coeff[(size_t)j];
                const Real j_cj = (Real)j * cj;

                psi += cj * S * c_nm;
                u   += j_cj * C * c_nm;
                v   += j_cj * S * s_nm;
            }

            rhs_out[(size_t)(m + 9)] = psi - z[8] - z[7] * z[(size_t)(m + 10)];
            rhs_out[(size_t)(n + m + 10)] = 0.5 * (((-z[7] + u) * (-z[7] + u)) + v * v)
                                           + z[(size_t)(m + 10)] - z[9];
        }

        // sum of squares
        Real ss = 0.0;
        for (int i = 1; i <= num; ++i) ss += rhs_out[(size_t)i] * rhs_out[(size_t)i];
        return ss;
    }

    // ----------------------------------------------------------------------
    // SVD solve with Press-style truncation: wmin = wmax*1e-12 (Python parity)
    // Solves A x = b for square A (num x num), returning x.
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
    // Port of Python _newton()
    // ----------------------------------------------------------------------
    Real newton_step(int /*iter_count*/) {
        const Real ss0 = eqns(rhs1);
        if (!std::isfinite(ss0)) {
            throw std::runtime_error("Non-finite residual norm at start of Newton step.");
        }

        const std::vector<Real> z0 = z;

        std::vector<Real> A((size_t)num * (size_t)num, 0.0);
        std::vector<Real> b((size_t)num, 0.0);

        // finite-difference Jacobian (column-wise)
        for (int i = 1; i <= num; ++i) {
            Real h = (Real)0.01 * z0[(size_t)i];
            if (std::abs(z0[(size_t)i]) < (Real)1e-4) h = (Real)1e-5;
            if (std::abs(h) > 1.0) h = std::copysign((Real)1.0, h);

            z[(size_t)i] = z0[(size_t)i] + h;
            eqns(rhs2);
            z[(size_t)i] = z0[(size_t)i];

            b[(size_t)(i - 1)] = -rhs1[(size_t)i];
            const Real inv_h = 1.0 / h;
            for (int r = 1; r <= num; ++r) {
                A[(size_t)(r - 1) * (size_t)num + (size_t)(i - 1)] =
                    (rhs2[(size_t)r] - rhs1[(size_t)r]) * inv_h;
            }
        }

        std::vector<Real> dx = svd_solve(A, b, num);
        for (Real v : dx) {
            if (!std::isfinite(v)) throw std::runtime_error("Non-finite Newton correction vector (dx).");
        }

        // Backtracking: prefer alpha=1, reduce if it worsens residuals or violates kd>0
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

        // mean absolute correction on the free surface nodes (indices 10..n+10 inclusive)
        Real corr = 0.0;
        for (int i = 10; i <= n + 10; ++i) corr += std::abs(z_best[(size_t)i] - z0[(size_t)i]);
        corr /= (Real)(n + 1);

        return corr;
    }

    // ----------------------------------------------------------------------
    // Port of Python _compute_Y_and_B()
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

        // Refresh Tanh[] for post-processing / kinematics.
        for (int i = 1; i <= n; ++i) Tanh[(size_t)i] = std::tanh((Real)i * z[1]);
    }

    // ----------------------------------------------------------------------
    // Port of Python _surface_keta()
    // ----------------------------------------------------------------------
    Real surface_keta(Real X) const {
        Real kEta = 0.0;
        for (int j = 1; j < n; ++j) kEta += Y[(size_t)j] * std::cos((Real)j * X);
        kEta += 0.5 * Y[(size_t)n] * std::cos((Real)n * X);
        return kEta;
    }

    // ----------------------------------------------------------------------
    // Port of Python _point(): returns (u_nd, v_nd, dudt_nd) in depth scaling.
    // ----------------------------------------------------------------------
    void point(Real X, Real Yloc, Real& u_out, Real& v_out, Real& dudt_out) const {
        const Real kd = z[1];
        const Real kd_sqrt = std::sqrt(kd);

        const Real c_nd  = z[4] / kd_sqrt; // c/√(gd)
        const Real ce_nd = z[5] / kd_sqrt; // ce/√(gd)

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
    // Mean square bed orbital velocity ub^2 = < (u_b - Uc)^2 >
    // ----------------------------------------------------------------------
    Real mean_square_bed_orbital_velocity(int nph = 720) {
        if (!(d > 0.0) || !(T_target > 0.0)) return 0.0;
        const int Nph = std::max(36, nph);

        Real ub2 = 0.0;
        for (int i = 0; i < Nph; ++i) {
            const Real ph = (2.0 * Phys::PI) * (Real)i / (Real)Nph;
            Real u_abs, w_abs, ax;
            get_kinematics(0.0, ph, u_abs, w_abs, ax); // bed: z_bed=0
            const Real u_orb = u_abs - Uc;
            ub2 += u_orb * u_orb;
        }
        return ub2 / (Real)Nph;
    }

    // ----------------------------------------------------------------------
    // Depth-scaled momentum flux S/(rho g d^2) in moving frame (phase-invariant).
    // ----------------------------------------------------------------------
    Real momentum_flux_S_depth(Real phase = 0.0, int npts = 1200) {
        const Real kd = z[1];
        if (!(kd > 0.0)) return 0.0;

        const Real c_nd = z[4] / std::sqrt(kd);
        const Real R_nd = 1.0 + z[9] / kd;

        const Real X = phase;
        const Real kEta = surface_keta(X);
        const Real eta_over_d = 1.0 + kEta / kd; // y at free surface

        const int Np = std::max(50, npts);
        const Real dy = eta_over_d / (Real)(Np - 1);

        Real integ = 0.0;
        for (int i = 0; i < Np; ++i) {
            const Real y = dy * (Real)i;
            const Real Yloc = kd * (y - 1.0); // Y = kd*(y-1)

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
    // Integral properties (Python _calc_integral_props_cpp parity)
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

        // dimensionalise
        EnergyDensity = Phys::RHO * g * (d * d) * E_depth;            // [J/m^2]
        Sxx           = Phys::RHO * g * (d * d) * Sxx_depth;          // [N/m]
        Power         = Phys::RHO * std::pow(g, 1.5) * std::pow(d, 2.5) * F_depth; // [W/m]

        MomentumFluxDepth = momentum_flux_S_depth(0.0, 1200);
        MomentumFlux      = Phys::RHO * g * (d * d) * MomentumFluxDepth;           // [N/m]

        Impulse = Phys::RHO * std::sqrt(g * std::pow(d, 3.0)) * I_depth;           // [kg/(m·s)]

        BernoulliR = R_dimless * g * d;
        MassTransport = cs_dimless * std::sqrt(g * d);

        // Convenience values for reporting
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
    // Core continuation + Newton loop + dimensional post-process (Python parity)
    // ----------------------------------------------------------------------
    void solve_internal() {
        const Real dhe = Height / (Real)nstep;
        const Real dho = MaxH   / (Real)nstep;

        // continuation in height
        for (int ns = 1; ns <= nstep; ++ns) {
            height = (Real)ns * dhe;
            Hoverd = (Real)ns * dho;

            // initial/extrapolated guess
            if (ns == 1) {
                init_linear();
            } else {
                for (int i = 1; i <= num; ++i) z[(size_t)i] = 2.0 * SOL(i, 2) - SOL(i, 1);

                // fallback if extrapolation yields invalid state
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
                    // Retry once from last converged state if first iteration failed.
                    if (ns > 1 && it == 1) {
                        for (int i = 1; i <= num; ++i) z[(size_t)i] = SOL(i, 2);
                        err = newton_step(it);
                    } else {
                        throw;
                    }
                }

                if (!std::isfinite(err)) throw std::runtime_error("Non-finite Newton correction.");

                // IMPORTANT: update continuation storage BEFORE convergence break (Python parity)
                if (ns == 1) {
                    for (int i = 1; i <= num; ++i) SOL(i, 2) = z[(size_t)i];
                } else {
                    for (int i = 1; i <= num; ++i) SOL(i, 1) = SOL(i, 2);
                    for (int i = 1; i <= num; ++i) SOL(i, 2) = z[(size_t)i];
                }

                // protect linear algebra calls on diverging states
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

            // update Y and B for this step (C++ does this each step)
            compute_Y_and_B();
        }

        // ------------------------- dimensional post-process --------------------
        const Real kd = z[1];
        if (!(kd > 0.0) || !std::isfinite(kd)) throw std::runtime_error("Invalid wavenumber (kd).");

        const Real k_phys = kd / d;
        const Real L_phys = 2.0 * Phys::PI / k_phys;

        const Real c_dimless = z[4] / std::sqrt(kd); // c / sqrt(g d)
        const Real c_phys = c_dimless * std::sqrt(g * d);

        if (!(L_phys > 0.0) || !std::isfinite(L_phys)) throw std::runtime_error("Invalid wavelength.");
        if (!std::isfinite(c_phys)) throw std::runtime_error("Invalid celerity.");

        // surface nodes correspond to m*pi/n (half wave: crest->trough)
        for (int m = 0; m <= n; ++m) {
            const Real kEta = z[(size_t)(10 + m)]; // k(eta-d) at node
            eta_nodes[(size_t)m] = d * (1.0 + kEta / kd);
        }

        k = k_phys;
        L = L_phys;
        c = c_phys;

        // store Bj as 0-based array
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

        // Kinematics summary
        {
            Real w, ax;
            get_kinematics(0.0, 0.0, u_bed, w, ax);
        }

        // Quadratic bed shear estimate kept as-is (engineering heuristic)
        const Real cf_est = 0.005;
        tau_bed = 0.5 * Phys::RHO * cf_est * (u_bed * u_bed);
        ExcursionBed = std::abs(u_bed) * T_target / (2.0 * Phys::PI);

        // crest/trough surface velocities
        Real w, ax, u_trough;
        get_kinematics(d + eta_crest, 0.0, u_surf, w, ax);
        get_kinematics(d + eta_trough, Phys::PI, u_trough, w, ax);
        asymmetry = (std::abs(u_trough) > 0.0) ? std::abs(u_surf / u_trough) : 0.0;

        // scan phases for max vertical velocity and horizontal acceleration on the surface
        acc_max = 0.0;
        w_max   = 0.0;
        for (int i = 0; i < 40; ++i) {
            const Real X = (Real)i * Phys::PI / (Real)39; // linspace(0,pi,40)
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
//  OUTPUT FORMATTING (Exact match to fenton_gui.py / output.txt)
// ==============================================================================

namespace ReportFmt {

static constexpr int W = 107; // report width (characters) including borders

// Count UTF-8 codepoints (good enough for monospace alignment in the report file)
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
    // C++ may print "inf"/"nan" in uppercase depending on locale; normalize:
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
        // fmt_cell will truncate it later to match Python behavior.
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

    // Header
    std::string line = "|";
    for (size_t i = 0; i < headers.size(); ++i) {
        Cell c = Cell::strv(headers[i]);
        line += " " + fmt_cell(c, col_w[i], aligns[i]) + " |";
    }
    out << line << "\n";
    table_sep(out, col_w);

    // Body
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
// Report generator (runs the solver for no-current and with-current cases).
// ----------------------------------------------------------------------------
static std::string generate_output(double H_in, double T_in, double d_in, double Uc_in) {
    using namespace ReportFmt;

    // Case A: No current
    FentonStreamFunction solver0(H_in, T_in, d_in, 0.0);
    solver0.solve();

    // Case B: With ambient current (Eulerian)
    FentonStreamFunction solverC(H_in, T_in, d_in, Uc_in);
    solverC.solve();

    const bool has_current = (Uc_in != 0.0);

    // Numerical sanity checks (mirror Python _solver_status behaviour)
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
    rows.push_back(Row{ false, "", { Cell::strv("Celerity / phase speed (c)"), Cell::numv(solver0.c), wc_num(solverC.c), Cell::strv("m/s") } });
    rows.push_back(Row{ false, "", { Cell::strv("c/√(gd)"), Cell::numv(solver0.c / sqrt_gd), wc_num(solverC.c / sqrt_gd), Cell::strv("-") } });
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
    rows.push_back(Row{ false, "", { Cell::strv(u8"Group velocity (C𝗀 = F/E)ㅤ"), Cell::numv(solver0.Cg), wc_num(solverC.Cg), Cell::strv("m/s") } });

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

    // ------------------------ SOLUTION-FLAT tables (exact set) -------------------
    auto print_solution_flat = [&](const FentonStreamFunction& slv, const std::string& title) {
        const std::vector<std::string> h = { "#", "PARAMETER", "value", "adim param", "adim value" };
        const std::vector<int> cw = { 2, 37, 13, 25, 14 };
        const std::vector<std::string> al = { "right", "left", "right", "left", "right" };

        const Real g_ = slv.g;
        const Real d_ = slv.d;
        const Real H_ = slv.H_target;
        const Real T_ = slv.T_target;
        const Real L_ = slv.L;
        const Real c_ = slv.c;

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
        rr.push_back(Row{ false, "", { Cell::numv(19), Cell::strv("Wave power"),                      Cell::strv(f5(kW(slv.Power))),         Cell::strv(u8"F/(ρg³ᐟ²d⁵ᐟ²)ㅤㅤ"),               Cell::strv(f5(slv.F_depth)) } });

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
    add_term("F", "Wave power (energy flux).", u8"kW/m ; F/(ρg³ᐟ²d⁵ᐟ²)ㅤ");
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
// Convert \n to \r\n for Win32 multiline EDIT control display.
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
//  GUI (Win32): Robust UI-thread updates (no direct control updates from worker)
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

// The formatted report uses a fixed line width (output.txt reference).
// We size the output area so the full line width is visible without horizontal clipping.
static constexpr int OUTPUT_COLS = 107;

static bool parse_double_w(const std::wstring& in, double& out_val) {
    std::wstring s = in;
    // accept decimal comma
    std::replace(s.begin(), s.end(), L',', L'.');

    // trim
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

    // Heavy compute
    std::string txt = generate_output(p->H, p->T, p->d, p->Uc);

    // Write output.txt (Python parity)
    try {
        std::ofstream f("output.txt", std::ios::binary);
        f.write(txt.data(), (std::streamsize)txt.size());
    } catch (...) {
        // ignore
    }

    // Display uses Windows newlines
    std::string display_txt = to_windows_newlines(txt);

    // Convert to UTF-16
    const int wlen = MultiByteToWideChar(CP_UTF8, 0, display_txt.c_str(), -1, nullptr, 0);
    std::wstring* wres = new std::wstring((size_t)wlen, L'\0');
    MultiByteToWideChar(CP_UTF8, 0, display_txt.c_str(), -1, &(*wres)[0], wlen);

    // Post back to UI thread (safe)
    if (!g_closing.load() && g_hMain && IsWindow(g_hMain)) {
        PostMessageW(g_hMain, WM_APP_RESULT, 0, (LPARAM)wres);
    } else {
        delete wres;
    }
    return 0;
}

static void resize_gui_to_fit_output(HWND hwnd) {
    if (!g_hOutput || !IsWindow(g_hOutput) || !g_hMonoFont) return;

    // Measure monospace glyph width in pixels for the active output font.
    int char_w = 10;
    {
        HDC hdc = GetDC(g_hOutput);
        if (hdc) {
            HFONT old = (HFONT)SelectObject(hdc, g_hMonoFont);
            SIZE sz{};
            if (GetTextExtentPoint32W(hdc, L"0", 1, &sz) && sz.cx > 0) {
                char_w = (int)sz.cx;
            }
            SelectObject(hdc, old);
            ReleaseDC(g_hOutput, hdc);
        }
    }

    // Remove internal edit margins so sizing is deterministic.
    SendMessageW(g_hOutput, EM_SETMARGINS, EC_LEFTMARGIN | EC_RIGHTMARGIN, MAKELPARAM(0, 0));

    // Account for the vertical scrollbar inside the EDIT control.
    const int vscroll_w = GetSystemMetrics(SM_CXVSCROLL);
    const int extra_px = 8; // cushion for borders/rounding

    const int desired_output_w = OUTPUT_COLS * char_w + vscroll_w + extra_px;

    // Layout constants used in WM_CREATE.
    const int output_x = 300;
    const int top_margin = 20;
    const int right_margin = 20;
    const int bottom_margin = 20;

    // Keep the current output height.
    RECT out_rc{};
    GetWindowRect(g_hOutput, &out_rc);
    const int out_h = (int)(out_rc.bottom - out_rc.top);

    // Resize output control.
    SetWindowPos(g_hOutput, nullptr,
        output_x, top_margin,
        desired_output_w, out_h,
        SWP_NOZORDER);

    // Resize the main window so the *client* width contains the output exactly.
    const int desired_client_w = output_x + desired_output_w + right_margin;
    const int desired_client_h = top_margin + out_h + bottom_margin;

    RECT rc{ 0, 0, desired_client_w, desired_client_h };
    const DWORD style = (DWORD)GetWindowLongPtrW(hwnd, GWL_STYLE);
    const DWORD exstyle = (DWORD)GetWindowLongPtrW(hwnd, GWL_EXSTYLE);
    AdjustWindowRectEx(&rc, style, FALSE, exstyle);
    const int win_w = rc.right - rc.left;
    const int win_h = rc.bottom - rc.top;

    SetWindowPos(hwnd, nullptr, 0, 0, win_w, win_h, SWP_NOMOVE | SWP_NOZORDER);
}

static LRESULT CALLBACK WndProc(HWND hwnd, UINT msg, WPARAM wParam, LPARAM lParam) {
    switch (msg) {
    case WM_CREATE: {
        g_hMain = hwnd;

        g_hUIFont = CreateFontW(19, 0, 0, 0, FW_NORMAL, FALSE, FALSE, FALSE,
            DEFAULT_CHARSET, OUT_DEFAULT_PRECIS, CLIP_DEFAULT_PRECIS,
            DEFAULT_QUALITY, DEFAULT_PITCH | FF_SWISS, L"Segoe UI");

        g_hMonoFont = CreateFontW(20, 0, 0, 0, FW_NORMAL, FALSE, FALSE, FALSE,
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

        g_hBtnCalc = CreateWindowW(L"BUTTON", L"CALCULATE HYDRODYNAMICS",
            WS_CHILD | WS_VISIBLE | BS_DEFPUSHBUTTON,
            20, y, 260, 40, hwnd, (HMENU)IDC_BTN_CALC, nullptr, nullptr);
        SendMessageW(g_hBtnCalc, WM_SETFONT, (WPARAM)g_hUIFont, TRUE);

        g_hOutput = CreateWindowW(L"EDIT", L"",
            WS_CHILD | WS_VISIBLE | WS_BORDER | ES_MULTILINE | WS_VSCROLL | ES_READONLY | WS_HSCROLL,
            300, 20, 800, 520, hwnd, (HMENU)IDC_OUTPUT, nullptr, nullptr);
        SendMessageW(g_hOutput, WM_SETFONT, (WPARAM)g_hMonoFont, TRUE);

        // Make the output area wide enough for the fixed-width report.
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
            SetWindowTextW(g_hBtnCalc, L"CALCULATE HYDRODYNAMICS");
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

    HWND hwnd = CreateWindowExW(
        0, L"FentonClassStable", L"Fenton Wave Solver",
        WS_OVERLAPPED | WS_CAPTION | WS_SYSMENU | WS_MINIMIZEBOX,
        CW_USEDEFAULT, CW_USEDEFAULT, 1140, 620,
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
