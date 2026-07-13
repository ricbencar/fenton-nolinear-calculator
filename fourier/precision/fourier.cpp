/*
================================================================================
fourier.cpp — production spectral solver for steady nonlinear gravity waves
================================================================================

PROGRAM PURPOSE
---------------
This program computes one steady, two-dimensional, periodic, finite-amplitude
gravity wave over a horizontal bed.  The free surface, wave celerity, mean-flow
quantities, stream-function coefficients, pressure field, velocity field, and
integral quantities are solved as one coupled nonlinear boundary-value problem.
An optional collinear current is prescribed using either the Eulerian-mean or
mass-transport definition selected in data.dat.

The executable is intentionally self-contained.  When data.dat is available,
it reads one annotated case from that file.  When data.dat is absent, it requests
every required numerical input interactively from the console; no physical input
is defaulted or inferred.  The program then solves the wave with a fixed 100-mode
spectral representation and writes solution.res, surface.res, and flowfield.res.
No command-line numerical parameters are required.

PHYSICAL MODEL
--------------
The calculation assumes:

  * homogeneous, incompressible, inviscid fluid;
  * irrotational two-dimensional motion;
  * a periodic wave of permanent form;
  * a horizontal impermeable bed;
  * constant atmospheric pressure on the free surface;
  * no breaking, viscosity, turbulence closure, vorticity, or surface tension.

The governing equations are written in a coordinate frame moving with the wave.
In that frame the free surface is steady and the fluid domain is periodic in the
phase coordinate X = k(x-c t), where k is the wavenumber and c is celerity.

DIMENSIONLESS VARIABLES
-----------------------
The nonlinear unknowns use wavenumber-based scaling:

  kd               dimensionless mean depth
  kH               dimensionless crest-to-trough height
  T*sqrt(gk)       dimensionless apparent period
  c*sqrt(k/g)      dimensionless celerity
  u_E*sqrt(k/g)    Eulerian mean current
  u_M*sqrt(k/g)    mass-transport mean current
  U*sqrt(k/g)      mean fluid speed in the moving frame
  q*sqrt(k^3/g)    wave-volume-flux variable
  r*k/g            Bernoulli constant

The input and engineering output are depth-scaled.  Conversion between the two
systems is performed only after the nonlinear state has converged.

SPECTRAL DISCRETIZATION
-----------------------
The half-wave interval is collocated at

  X_m = m*pi/N,  m = 0,...,N,

with N fixed at 100.  Wave symmetry makes the crest-to-trough half wave
sufficient.  The free surface is represented by N+1 nodal ordinates and the
stream function by N Fourier coefficients.

The finite-depth vertical basis is evaluated in the stable exponential form

  S_j(eta,kd) =
      [exp(j*eta) - exp(-j*eta - 2*j*kd)] / [1 + exp(-2*j*kd)]

  C_j(eta,kd) =
      [exp(j*eta) + exp(-j*eta - 2*j*kd)] / [1 + exp(-2*j*kd)]

rather than through direct hyperbolic functions.  The scaled representation
avoids overflow and cancellation over the admissible fluid domain.

STATE VECTOR AND RESIDUAL SYSTEM
--------------------------------
For N Fourier modes, the active state contains 2*N+10 entries:

   1  kd
   2  kH
   3  T*sqrt(gk)
   4  c*sqrt(k/g)
   5  Eulerian current
   6  mass-transport current
   7  mean wave-frame speed
   8  wave-volume-flux variable
   9  Bernoulli constant
  10 ... 10+N       free-surface nodal ordinates
  11+N ... 10+2*N   stream-function Fourier coefficients

The first eight residuals impose global wave, current, mean-level, and height
constraints.  Each of the N+1 collocation nodes then contributes one kinematic
and one dynamic free-surface equation.  The resulting system is square:

  8 + (N+1) + (N+1) = 2*N+10 equations.

NONLINEAR SOLUTION METHOD
-------------------------
The target wave is reached by adaptive continuation in wave height.  Each new
continuation point is predicted from the exact tangent of the converged branch
and corrected by a trust-region Newton iteration with an exact analytical
Jacobian.  Powell dogleg globalization combines the Newton and scaled Cauchy
directions.  Failed trials are rejected and the continuation step is reduced;
accepted inexpensive solves permit controlled step growth.

Positive state variables use logarithmic retractions, so depth, height, period,
and celerity remain positive without clipping.  Every candidate state is checked
for finite values, bed clearance, exponential safety, and residual decrease.
Unconverged states are never written to result files.

LINEAR ALGEBRA
--------------
Each Newton correction solves an equilibrated dense system.  The primary path is
partial-pivoted LU with long-double triangular solves and iterative refinement.
If the Jacobian is rank deficient or poorly conditioned, the solver falls back
to column-pivoted Householder QR and then to adaptively regularized least
squares.  The hierarchy preserves speed for normal cases while retaining a
usable descent direction near difficult states.

INPUT FILE
----------
The program first searches for data.dat in the current working directory and
then beside the executable.  When the file is found, it must contain these seven
non-empty lines:

  1. case title
  2. H/d
  3. Period or Wavelength
  4. T*sqrt(g/d) or L/d, respectively
  5. current criterion: 1 = Eulerian, 2 = mass transport
  6. dimensionless current magnitude u/sqrt(gd)
  7. FINISH

Only the leading token of lines 2 through 6 is interpreted.  Descriptive text
may follow values, for example:

  Test wave
  0.6        H/d
  Period     Measure of length: "Wavelength" or "Period"
  12.6042742750227
  1          Current criterion
  0.142808698122903
  FINISH

If data.dat is not found in either location, the program enters interactive
mode and requests, in order:

  1. H/d
  2. Period or Wavelength
  3. T*sqrt(g/d) or L/d, respectively
  4. current criterion: 1 = Eulerian, 2 = mass transport
  5. dimensionless current magnitude u/sqrt(gd)

Every field must be entered explicitly.  Invalid entries are rejected and
requested again; end-of-input aborts the run.  No default numerical values are
used.

The Fourier order is always 100.  The continuation step count is selected and
adapted internally from the requested H/d and observed nonlinear convergence.

OUTPUT FILES
------------
After successful convergence the program writes:

  solution.res   wave summary, integral quantities, invariants, and coefficients
  surface.res    free-surface samples and pressure consistency
  flowfield.res  velocity, acceleration, pressure, and Bernoulli profiles

USAGE
-----
Compile the source and run without arguments:

  ./fourier

On Windows the executable name is normally fourier.exe:

  fourier.exe

When data.dat is present in the working directory or beside the executable, it
is used.  Otherwise the program requests all required numerical inputs from the
console.  Pressing end-of-input before all fields are supplied aborts the run.

COMPILATION
-----------
Windows / MSYS2 / MinGW-w64, optimized build with OpenMP:

  g++ -std=c++20 fourier.cpp -O3 -march=native -flto=auto -fopenmp \
      -static-libgcc -static-libstdc++ -o fourier.exe

Linux, optimized build with OpenMP:

  g++ -std=c++20 fourier.cpp -O3 -march=native -flto=auto -fopenmp -o fourier

Portable build without OpenMP or architecture-specific optimization:

  g++ -std=c++20 fourier.cpp -O3 -o fourier

Recommended diagnostic build during development:

  g++ -std=c++20 fourier.cpp -O1 -g -Wall -Wextra -Wpedantic \
      -Wconversion -Wshadow -fsanitize=address,undefined -o fourier_debug

PROGRAM STRUCTURE AND ROUTINE MAP
---------------------------------
The translation unit is organized in the same order as the numerical workflow:

  1. OneBasedVector and DenseMatrix
     Small ownership types that centralize numerical indexing and contiguous
     storage.  They contain no wave physics.

  2. Input and configuration
     trim_line(), parse_leading_value(), ascii_iequals(), locate_data_file(),
     read_case_file(), read_case_interactively(), and configure_wave_case()
     implement the complete file and console input contracts and commit
     configuration only after validation succeeds.

  3. Spectral representation
     build_spectral_tables(), cosine_at_node(), sine_at_node(), and
     evaluate_vertical_basis() prepare the modal basis and its exact depth
     derivatives.

  4. Nonlinear equations
     evaluate_residuals() assembles all global and free-surface equations.
     assemble_jacobian() constructs their exact analytical Jacobian.
     state_is_admissible(), residual_is_finite(), and jacobian_is_finite()
     prevent invalid states from entering the nonlinear globalization logic.

  5. Scaled linear algebra
     build_scaled_linearization() equilibrates the system.  factorize_lu() and
     solve_lu() provide the fast path; solve_rank_revealing_qr() and
     solve_regularized_least_squares() provide robust fallbacks.  Iterative
     refinement checks the correction in the original scaled equations.

  6. Trust-region nonlinear correction
     build_cauchy_direction(), build_dogleg_step(), apply_retraction(), and
     solve_nonlinear_state() construct, constrain, test, accept, or reject each
     candidate state according to actual residual reduction.

  7. Adaptive continuation
     assemble_continuation_derivative(), compute_continuation_tangent(),
     predict_from_tangent(), predict_from_secant_history(), and
     continue_to_target() trace the physical solution branch from small wave
     height to the requested H/d with exact rollback after failed trials.

  8. Post-processing and output
     update_output_coefficients(), evaluate_flow_point(),
     compute_derived_quantities(), and the three dedicated result writers
     transform the converged state into engineering quantities and files.

  9. Application control
     allocate_workspace(), solve_wave(), run(), and main() enforce the sequence
     parse -> allocate -> solve -> verify -> write.  Result files are created
     only after a converged state exists.

================================================================================
*/

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <cstdio>
#include <exception>
#include <limits>
#include <numbers>
#include <vector>

#include <filesystem>
#include <fstream>
#include <iostream>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>

#if defined(_WIN32)

#ifndef NOMINMAX
#define NOMINMAX
#endif
#ifndef WIN32_LEAN_AND_MEAN
#define WIN32_LEAN_AND_MEAN
#endif
#include <windows.h>
#endif

#if defined(__linux__)

#include <limits.h>
#include <unistd.h>
#endif

namespace fs = std::filesystem;

namespace wave {

/*
Close a C output stream when its owning FilePtr leaves scope.  The deleter is
noexcept because it is used during stack unwinding; explicit end-of-run close
checks are performed separately when an I/O error must be reported.
*/
struct FileCloser {
  void operator()(std::FILE *f) const noexcept {
    if (f)
      std::fclose(f);
  }
};

using FilePtr = std::unique_ptr<std::FILE, FileCloser>;

/*
================================================================================
ONE-BASED NUMERICAL STORAGE
================================================================================

The nonlinear system is described with equation and unknown numbers beginning at
one. Rather than exposing repeated signed-to-unsigned conversions at every
access, OneBasedVector centralizes that convention in a small container adapter.

Element zero is intentionally allocated. It is not part of the mathematical
system, but retaining it makes the source equations read exactly like their
numbered definitions and eliminates error-prone index shifts throughout the
Jacobian, continuation, and output code.

The class exposes only the operations used by this translation unit. Storage
remains contiguous, ownership remains automatic, and the indexing convention is
confined to one auditable location.
================================================================================
*/
template <typename Value> class OneBasedVector {
public:
  using Storage = std::vector<Value>;
  using iterator = typename Storage::iterator;
  using const_iterator = typename Storage::const_iterator;

  // Construct an empty adapter. Storage is allocated later by assign().
  OneBasedVector() = default;

  // Allocate count contiguous elements, including the intentionally unused
  // index zero, and initialize every element to the supplied value.
  explicit OneBasedVector(std::size_t count,
                          const Value &initial_value = Value{})
      : values_(count, initial_value) {}

  // Replace the complete storage block. This is used when the fixed nonlinear
  // dimension is installed and guarantees value initialization of work arrays.
  void assign(std::size_t count, const Value &value) {
    values_.assign(count, value);
  }
  // Remove logical elements while retaining standard vector lifetime rules.
  void clear() noexcept { values_.clear(); }

  // Reserve storage before rebuilding integer pivot arrays during allocation.
  void reserve(std::size_t count) { values_.reserve(count); }

  // Report physical storage size. The active numerical indices are 1..size-1.
  [[nodiscard]] std::size_t size() const noexcept { return values_.size(); }
  // Expose the allocator limit for checked matrix extent calculations.
  [[nodiscard]] std::size_t max_size() const noexcept {
    return values_.max_size();
  }

  // Raw contiguous access is restricted to DenseMatrix row construction.
  Value *data() noexcept { return values_.data(); }
  const Value *data() const noexcept { return values_.data(); }

  // Access a numbered component. Assertions catch programming errors in debug
  // builds without adding bounds-check overhead to optimized numerical kernels.
  Value &operator[](int index) noexcept {
    assert(index >= 0);
    assert(static_cast<std::size_t>(index) < values_.size());
    return values_[static_cast<std::size_t>(index)];
  }

  // Const overload with the same one-based indexing contract.
  const Value &operator[](int index) const noexcept {
    assert(index >= 0);
    assert(static_cast<std::size_t>(index) < values_.size());
    return values_[static_cast<std::size_t>(index)];
  }

  // Iterators permit standard algorithms and range-for loops over the complete
  // allocated block, including element zero when bulk initialization is wanted.
  iterator begin() noexcept { return values_.begin(); }
  const_iterator begin() const noexcept { return values_.begin(); }
  iterator end() noexcept { return values_.end(); }
  const_iterator end() const noexcept { return values_.end(); }

private:
  Storage values_;
};

using Vector = OneBasedVector<double>;
using IndexVector = OneBasedVector<int>;

/*
================================================================================
DENSE SQUARE MATRIX
================================================================================

All nonlinear systems in this program are square. DenseMatrix owns one
contiguous row-major block and presents one-based row pointers to the numerical
kernels. Contiguous storage improves cache locality and avoids fragmented
per-row allocations.

The maximum matrix dimension is fixed by the allocated Fourier order. Individual
solves use the leading kSystemSize rows and columns.
================================================================================
*/
class DenseMatrix {
public:
  // Default construction creates a zero-dimensional matrix for later
  // allocation.
  DenseMatrix() = default;

  // Convenience construction for a known square dimension.
  explicit DenseMatrix(int dimension) { resize(dimension); }

  // Allocate a (dimension+1)^2 contiguous block so rows and columns may be
  // addressed from one. Existing contents are discarded and zero-initialized.
  void resize(int dimension) {
    if (dimension < 0) {
      throw std::invalid_argument("DenseMatrix dimension cannot be negative");
    }

    dimension_ = dimension;
    const std::size_t extent = static_cast<std::size_t>(dimension + 1);

    // Guard the row-major extent calculation explicitly. The solver's practical
    // dimensions are far below this limit, but the check prevents silent
    // wraparound if an invalid future configuration reaches this low-level
    // type.
    if (extent != 0U && extent > values_.max_size() / extent) {
      throw std::length_error(
          "DenseMatrix allocation exceeds addressable size");
    }

    values_.assign(extent * extent, 0.0);
  }

  // Return the beginning of one row. The returned pointer supports the natural
  // matrix[row][column] notation used throughout the dense kernels.
  double *operator[](int row) noexcept {
    assert(row >= 0 && row <= dimension_);
    return values_.data() + static_cast<std::size_t>(row) *
                                static_cast<std::size_t>(dimension_ + 1);
  }

  // Const row access for factorization back-solves and matrix-vector products.
  const double *operator[](int row) const noexcept {
    assert(row >= 0 && row <= dimension_);
    return values_.data() + static_cast<std::size_t>(row) *
                                static_cast<std::size_t>(dimension_ + 1);
  }

private:
  int dimension_ = 0;
  Vector values_;
};

class SteadyWaveApplication final {
private:
  /*
  ============================================================================
  NUMERICAL CONSTANTS AND ENUMERATIONS
  ============================================================================

  These constants define solver policy rather than wave physics. Keeping them in
  one block makes convergence behavior, continuation limits, and parallelization
  thresholds auditable without searching through algorithm bodies.
  ============================================================================
  */
  static constexpr double kPi = std::numbers::pi_v<double>;
  static constexpr int kConvergenceMaxIterations = 80;
  static constexpr double kConvergenceCriterion = 1.0e-10;
  static constexpr double kFinalResidualTolerance = 1.0e-12;
  static constexpr double kMinimumContinuationStep = 1.0e-7;
  static constexpr int kMaximumContinuationTrials = 512;
  static constexpr int kFourierOrder = 100;
  static constexpr double kSpectralDefectTolerance = 1.0e-9;
  static constexpr int kOpenMpMinimumOrder = 48;
  static constexpr int kSurfacePointsDefault = 50;
  static constexpr int kProfileCountDefault = 18;
  static constexpr int kVerticalPointCountDefault = 10;

  enum class LinearSolverKind { LU, QR, Regularized };
  enum class CurrentDefinition : int { Eulerian = 1, MassTransport = 2 };
  enum class WaveSpecification { Period, Wavelength };

  /*
  ============================================================================
  STATE AND RESIDUAL NUMBERING
  ============================================================================

  The first nine state entries are global wave quantities. They are followed by
  order+1 free-surface ordinates and order stream-function coefficients.

  The first eight residuals enforce global constraints. They are followed by
  order+1 kinematic conditions and order+1 dynamic conditions. The resulting
  square system therefore contains 2*order+10 equations and unknowns.
  ============================================================================
  */
  struct StateIndex {
    static constexpr int depth = 1;
    static constexpr int wave_height = 2;
    static constexpr int period = 3;
    static constexpr int celerity = 4;
    static constexpr int eulerian_current = 5;
    static constexpr int mass_transport_current = 6;
    static constexpr int wave_frame_speed = 7;
    static constexpr int flux = 8;
    static constexpr int bernoulli = 9;
    static constexpr int surface_begin = 10;
  };

  struct ResidualIndex {
    static constexpr int depth = 1;
    static constexpr int wave_specification = 2;
    static constexpr int dispersion = 3;
    static constexpr int eulerian_current = 4;
    static constexpr int mass_transport_current = 5;
    static constexpr int prescribed_current = 6;
    static constexpr int mean_surface = 7;
    static constexpr int wave_height = 8;
    static constexpr int kinematic_begin = 9;
  };

  /*
  ============================================================================
  OWNED PROGRAM STATE
  ============================================================================

  Data is grouped by responsibility so that each algorithm touches only the
  storage category it logically owns. This avoids a flat collection of unrelated
  members and makes lifetime, mutation, and rollback behavior explicit.
  ============================================================================
  */

  struct ProblemDefinition {
    std::string title;
    CurrentDefinition current_definition = CurrentDefinition::Eulerian;
    WaveSpecification wave_specification = WaveSpecification::Period;
    double current_value = 0.0;
    double target_relative_height = 0.0;
    double target_height_parameter = 0.0;
    double specified_wavelength = 0.0;
    double specified_period = 0.0;
  };

  // The nonlinear dimension is a compile-time invariant.  Keeping the order and
  // system size out of runtime storage prevents accidental zero-initialization,
  // stale configuration, or memory corruption from changing the discretization.
  static constexpr int kSystemSize =
      2 * kFourierOrder + StateIndex::surface_begin;

  struct Discretization {
    int continuation_steps = 0;
    int surface_point_count = kSurfacePointsDefault;
    int profile_count = kProfileCountDefault;
    int vertical_point_count = kVerticalPointCountDefault;
  };

  struct ContinuationState {
    double relative_height = 0.0;
    double height_parameter = 0.0;
  };

  struct SolverStatus {
    double residual_rms = std::numeric_limits<double>::infinity();
    double residual_max = std::numeric_limits<double>::infinity();
    double trust_radius = 1.0;
    double spectral_defect_max = std::numeric_limits<double>::infinity();
    LinearSolverKind linear_solver = LinearSolverKind::LU;
  };

  struct SpectralStorage {
    std::vector<double> cosine_table;
    std::vector<double> sine_table;
    Vector stream_coefficients;
    Vector surface_coefficients;
  };

  struct NonlinearStorage {
    Vector state;
    Vector residual;
    Vector trial_residual;
  };

  struct SolverWorkspace {
    DenseMatrix jacobian;
    DenseMatrix linear_matrix;
    Vector linear_rhs;
    Vector newton_direction;
    Vector cauchy_direction;
    Vector candidate_step;
    Vector base_state;
    Vector row_scale;
    Vector column_scale;
    Vector model_residual;
    IndexVector pivot;
    std::vector<long double> predictor_coefficients;
  };

  struct OutputHandles {
    // FilePtr objects in run() retain ownership for the entire write operation.
    // This lightweight bundle is passed explicitly, keeping output resources
    // out of persistent solver state and making write_outputs() referentially
    // clear.
    std::FILE *solution = nullptr;
    std::FILE *surface = nullptr;
    std::FILE *flowfield = nullptr;
  };

  /*
  Flush and close an owned output stream while preserving the close status.
  FilePtr normally closes automatically, but an explicit close at the end of a
  successful run is required to detect delayed write failures such as a full
  filesystem or a disconnected network volume.
  */
  [[nodiscard]] static bool close_output_file(FilePtr &file,
                                              const char *label) noexcept {
    std::FILE *stream = file.release();
    if (!stream)
      return true;

    bool success = true;
    if (std::fflush(stream) != 0 || std::ferror(stream) != 0) {
      std::fprintf(stderr, "Failed while writing %s.\n", label);
      success = false;
    }
    if (std::fclose(stream) != 0) {
      std::fprintf(stderr, "Failed while closing %s.\n", label);
      success = false;
    }
    return success;
  }

  ProblemDefinition problem_;
  Discretization grid_;
  ContinuationState continuation_;
  SolverStatus status_;
  SpectralStorage spectral_;
  NonlinearStorage nonlinear_;
  SolverWorkspace workspace_;

  /*
  ----------------------------------------------------------------------------
  INDEXING HELPERS
  ----------------------------------------------------------------------------

  These helpers are the only place where a Fourier mode or collocation node is
  converted into a state-vector or residual-vector position. Centralizing the
  layout prevents duplicated arithmetic and makes every later equation readable
  in terms of physical roles rather than raw offsets.

  For order N:
    - surface nodes occupy StateIndex::surface_begin ... +N;
    - stream coefficients occupy the following N positions;
    - kinematic residuals begin at ResidualIndex::kinematic_begin;
    - dynamic residuals follow the N+1 kinematic equations.
  ----------------------------------------------------------------------------
  */
  // Map collocation node m in [0,N] to its free-surface state component.
  int surface_state_index(int node) const noexcept {
    return StateIndex::surface_begin + node;
  }

  // Map Fourier mode j in [1,N] to its stream-function coefficient component.
  int coefficient_state_index(int mode) const noexcept {
    return StateIndex::surface_begin + kFourierOrder + mode;
  }

  // Map node m to the constant-stream-function residual row.
  int kinematic_residual_index(int node) const noexcept {
    return ResidualIndex::kinematic_begin + node;
  }

  // Map node m to the Bernoulli free-surface residual row.
  int dynamic_residual_index(int node) const noexcept {
    return StateIndex::surface_begin + kFourierOrder + node;
  }

  // Resolve the user-selected current convention once for residual/Jacobian
  // use.
  int selected_current_state_index() const noexcept {
    return problem_.current_definition == CurrentDefinition::Eulerian
               ? StateIndex::eulerian_current
               : StateIndex::mass_transport_current;
  }

  /*
  Return short human-readable labels used consistently in all output files.
  These functions contain presentation text only; no solver decision depends on
  the returned strings.
  */
  std::string method_description() const {
    return "# Solution by fixed " + std::to_string(kFourierOrder) +
           "-mode spectral collocation";
  }

  const char *current_name() const noexcept {
    return problem_.current_definition == CurrentDefinition::Eulerian
               ? "Eulerian"
               : "MassTransport";
  }

  /*
  ============================================================================
  SPECTRAL DISCRETIZATION AND NONLINEAR EQUATIONS
  ============================================================================
  */

  /*
  Precompute cos(m*j*pi/N) and sin(m*j*pi/N) for every collocation node m and
  Fourier mode j. The nonlinear residual and Jacobian evaluate these values
  O(N^2) times per iteration, so a contiguous lookup table avoids repeated
  transcendental calls in the dominant kernels.

  Long-double arithmetic is used for the angle before the result is stored as a
  double. This preserves accurate special values near crest and trough while
  keeping the hot table compact and cache-friendly.
  */
  void build_spectral_tables() {
    const std::size_t side = static_cast<std::size_t>(kFourierOrder + 1);
    spectral_.cosine_table.resize(side * side);
    spectral_.sine_table.resize(side * side);

    for (int node = 0; node <= kFourierOrder; ++node) {
      for (int mode = 0; mode <= kFourierOrder; ++mode) {
        const long double angle = std::numbers::pi_v<long double> *
                                  static_cast<long double>(node) *
                                  static_cast<long double>(mode) /
                                  static_cast<long double>(kFourierOrder);
        const std::size_t index = static_cast<std::size_t>(node) * side +
                                  static_cast<std::size_t>(mode);
        spectral_.cosine_table[index] = static_cast<double>(std::cos(angle));
        spectral_.sine_table[index] = static_cast<double>(std::sin(angle));
      }
    }
  }

  // Return cos(j*X_m) from the row-major table. Arguments are trusted internal
  // indices and are therefore not revalidated in this hot accessor.
  [[nodiscard]] inline double cosine_at_node(int node,
                                             int mode) const noexcept {
    return spectral_
        .cosine_table[static_cast<std::size_t>(node) *
                          static_cast<std::size_t>(kFourierOrder + 1) +
                      static_cast<std::size_t>(mode)];
  }

  // Return sin(j*X_m) using the same table layout as cosine_at_node().
  [[nodiscard]] inline double sine_at_node(int node, int mode) const noexcept {
    return spectral_.sine_table[static_cast<std::size_t>(node) *
                                    static_cast<std::size_t>(kFourierOrder + 1) +
                                static_cast<std::size_t>(mode)];
  }

  /*
  Construct a small-amplitude state at the current continuation parameter.

  The initialization has four responsibilities:
    1. obtain a positive first estimate of dimensionless depth;
    2. initialize period and celerity consistently;
    3. seed the selected current convention and global flow variables;
    4. create a cosine free surface and a first stream-function coefficient.

  This state is only a starting point for continuation. The final wave is
  defined exclusively by convergence of the complete nonlinear residual system.
  */
  void initialize_state() {
    /*
    --------------------------------------------------------------------------
    1. INITIAL DEPTH, PERIOD, AND CELERITY
    --------------------------------------------------------------------------

    Continuation begins at a small fraction of the requested H/d.  The first
    state must be close enough to the linear solution for the nonlinear
    corrector to identify the physical branch, but every quantity is
    subsequently adjusted by the complete residual system.

    Period input leaves kd unknown.  A smooth finite-depth dispersion estimate
    provides a positive starting value.  Wavelength input determines kd directly
    from L/d.  Neither relation is imposed as the final nonlinear solution; they
    are initialization formulas only.
    */
    if (problem_.wave_specification == WaveSpecification::Period) {
      const double sigma = 2.0 * kPi *
                           std::sqrt(continuation_.height_parameter /
                                     continuation_.relative_height);
      const double alpha = sigma * sigma;
      const double exp_argument = alpha * std::log(6.0 / 5.0);
      const double scale = std::exp(std::min(exp_argument, 20.0));
      const double tanh_argument = scale * std::sqrt(alpha);
      nonlinear_.state[StateIndex::depth] = alpha / std::tanh(tanh_argument);
    } else {
      nonlinear_.state[StateIndex::depth] = 2.0 * kPi *
                                            continuation_.height_parameter /
                                            continuation_.relative_height;
    }

    nonlinear_.state[StateIndex::wave_height] =
        nonlinear_.state[StateIndex::depth] * continuation_.relative_height;
    nonlinear_.state[StateIndex::celerity] =
        std::sqrt(std::tanh(nonlinear_.state[StateIndex::depth]));
    nonlinear_.state[StateIndex::period] =
        2.0 * kPi / nonlinear_.state[StateIndex::celerity];

    /*
    --------------------------------------------------------------------------
    2. CURRENT AND GLOBAL FLOW VARIABLES
    --------------------------------------------------------------------------

    The selected current definition receives the prescribed depth-scaled value,
    converted to native k-scaling through sqrt(kd).  The other current starts at
    zero and is determined by the coupled current and flux equations.  The
    moving-frame speed initially equals celerity, the wave-induced flux starts
    at zero, and the Bernoulli constant is initialized from the corresponding
    uniform-flow kinetic head.
    */
    nonlinear_.state[StateIndex::eulerian_current] = 0.0;
    nonlinear_.state[StateIndex::mass_transport_current] = 0.0;
    nonlinear_.state[selected_current_state_index()] =
        problem_.current_value * std::sqrt(nonlinear_.state[StateIndex::depth]);

    nonlinear_.state[StateIndex::wave_frame_speed] =
        nonlinear_.state[StateIndex::celerity];
    nonlinear_.state[StateIndex::flux] = 0.0;
    nonlinear_.state[StateIndex::bernoulli] =
        0.5 * nonlinear_.state[StateIndex::wave_frame_speed] *
        nonlinear_.state[StateIndex::wave_frame_speed];

    /*
    --------------------------------------------------------------------------
    3. LINEARIZED FREE-SURFACE AND STREAM-FUNCTION SEED
    --------------------------------------------------------------------------

    The first free-surface profile is a cosine with crest-to-trough height kH.
    All stream-function modes begin at zero except the fundamental coefficient,
    whose linear estimate supplies the required first-order motion.  The exact
    nonlinear solve immediately corrects these values.
    */
    build_spectral_tables();

    nonlinear_.state[surface_state_index(0)] =
        0.5 * nonlinear_.state[StateIndex::wave_height];
    for (int node = 1; node <= kFourierOrder; ++node) {
      nonlinear_.state[coefficient_state_index(node)] = 0.0;
      nonlinear_.state[surface_state_index(node)] =
          0.5 * nonlinear_.state[StateIndex::wave_height] *
          cosine_at_node(node, 1);
    }
    nonlinear_.state[coefficient_state_index(1)] =
        0.5 * nonlinear_.state[StateIndex::wave_height] /
        nonlinear_.state[StateIndex::wave_frame_speed];
  }

  struct VerticalBasisValues {
    double s = 0.0;
    double c = 0.0;
    double ds_dkd = 0.0;
    double dc_dkd = 0.0;
  };

  /*
  Evaluate one vertical mode without direct sinh/cosh calls.

  The expressions are rewritten with exponentials scaled by exp(-2*j*kd).
  This form avoids overflow and catastrophic cancellation over the admissible
  fluid domain while preserving exact derivatives with respect to kd.

  The routine returns false instead of propagating an overflow, invalid depth,
  or non-finite argument. Callers then reject the trial state through normal
  trust-region logic.
  */
  [[nodiscard]] static bool
  evaluate_vertical_basis(int mode, double eta, double kd,
                          VerticalBasisValues &basis) {
    /*
    Let a = j*eta and b = j*kd.  Direct evaluation with sinh(a), cosh(a), and
    tanh(b) can overflow even when the final basis is finite.  Multiplying the
    bed-reflected branch by exp(-2*b) produces an algebraically equivalent form
    in which every exponential has a controlled sign over the physical domain.
    */
    const double mode_value = static_cast<double>(mode);
    const double a = mode_value * eta;
    const double b = mode_value * kd;
    const double exponential_limit =
        0.98 * std::log(std::numeric_limits<double>::max());

    if (!std::isfinite(a) || !std::isfinite(b) || !(b > 0.0) ||
        a > exponential_limit) {
      return false;
    }

    const double depth_decay =
        (2.0 * b > exponential_limit) ? 0.0 : std::exp(-2.0 * b);
    const double reflected_exponent = -a - 2.0 * b;
    if (reflected_exponent > exponential_limit) {
      return false;
    }

    const double positive_branch = std::exp(a);
    const double reflected_branch = (reflected_exponent < -exponential_limit)
                                        ? 0.0
                                        : std::exp(reflected_exponent);
    const double denominator = 1.0 + depth_decay;

    basis.s = (positive_branch - reflected_branch) / denominator;
    basis.c = (positive_branch + reflected_branch) / denominator;

    /*
    Exact derivatives with respect to kd are required by the analytical
    Jacobian.  They are evaluated from the same scaled branches so the residual
    and Jacobian remain numerically consistent.
    */
    const double positive_times_decay = positive_branch * depth_decay;
    const double derivative_factor =
        2.0 * mode_value / (denominator * denominator);
    basis.ds_dkd =
        derivative_factor * (positive_times_decay + reflected_branch);
    basis.dc_dkd =
        derivative_factor * (positive_times_decay - reflected_branch);

    return std::isfinite(basis.s) && std::isfinite(basis.c) &&
           std::isfinite(basis.ds_dkd) && std::isfinite(basis.dc_dkd);
  }

  /*
  Compute sum(residual_i^2) using a scaled sum-of-squares algorithm.

  A naive sum can overflow even when the final norm is representable, or lose
  small components when magnitudes differ greatly. The scale/ssq recurrence is
  algebraically equivalent to the Euclidean norm calculation used in robust BLAS
  implementations.
  */
  [[nodiscard]] double stable_sum_squares(const Vector &residual) const {
    double scale = 0.0;
    double ssq = 1.0;
    for (int i = 1; i <= kSystemSize; ++i) {
      const double value = std::fabs(residual[i]);
      if (value == 0.0)
        continue;
      if (scale < value) {
        const double ratio = scale / value;
        ssq = 1.0 + ssq * ratio * ratio;
        scale = value;
      } else {
        const double ratio = value / scale;
        ssq += ratio * ratio;
      }
    }
    if (scale == 0.0)
      return 0.0;
    const long double result = static_cast<long double>(scale) *
                               static_cast<long double>(scale) *
                               static_cast<long double>(ssq);
    return result > static_cast<long double>(std::numeric_limits<double>::max())
               ? std::numeric_limits<double>::infinity()
               : static_cast<double>(result);
  }

  /*
  Assemble the complete square nonlinear system.

  Global equations enforce:
    - depth/height consistency;
    - period- or wavelength-based wave specification;
    - celerity-period consistency;
    - both mean-current definitions;
    - the selected prescribed current;
    - zero mean free-surface elevation;
    - crest-to-trough wave height.

  At every surface node the method then evaluates:
    - the constant-stream-function condition;
    - the free-surface Bernoulli condition.

  Long-double accumulators reduce cancellation in high-order modal sums. The
  returned value is the robust sum of squared residuals used as the merit
  function by the nonlinear solver.
  */
  [[nodiscard]] double evaluate_residuals(Vector &rhs) {
    /*
    Global equation R1: kH = kd*(H/d).  The continuation controller changes
    H/d while the nonlinear solver determines the corresponding kd and kH.
    */
    rhs[ResidualIndex::depth] =
        nonlinear_.state[StateIndex::wave_height] -
        nonlinear_.state[StateIndex::depth] * continuation_.relative_height;

    /*
    Global equation R2: impose the active continuation height in the selected
    input parameterization.
    */
    rhs[ResidualIndex::wave_specification] =
        problem_.wave_specification == WaveSpecification::Period
            ? nonlinear_.state[StateIndex::wave_height] -
                  continuation_.height_parameter *
                      nonlinear_.state[StateIndex::period] *
                      nonlinear_.state[StateIndex::period]
            : nonlinear_.state[StateIndex::wave_height] -
                  2.0 * kPi * continuation_.height_parameter;

    // Global equation R3: c*T = 2*pi in native scaling.
    rhs[ResidualIndex::dispersion] = nonlinear_.state[StateIndex::celerity] *
                                         nonlinear_.state[StateIndex::period] -
                                     2.0 * kPi;

    // Global equation R4: Eulerian-current definition.
    rhs[ResidualIndex::eulerian_current] =
        nonlinear_.state[StateIndex::eulerian_current] +
        nonlinear_.state[StateIndex::wave_frame_speed] -
        nonlinear_.state[StateIndex::celerity];

    // Global equation R5: mass-transport-current and flux definition.
    rhs[ResidualIndex::mass_transport_current] =
        nonlinear_.state[StateIndex::depth] *
            (nonlinear_.state[StateIndex::mass_transport_current] +
             nonlinear_.state[StateIndex::wave_frame_speed] -
             nonlinear_.state[StateIndex::celerity]) -
        nonlinear_.state[StateIndex::flux];

    /*
    Global equation R6: impose the selected depth-scaled current.  Conversion
    to native k-scaling introduces sqrt(kd).
    */
    rhs[ResidualIndex::prescribed_current] =
        nonlinear_.state[selected_current_state_index()] -
        problem_.current_value * std::sqrt(nonlinear_.state[StateIndex::depth]);

    /*
    Global equation R7: trapezoidal integral of the free surface over the
    symmetric half-wave. The common factor 1/N is omitted because the target is
    zero and constant scaling does not change the root.
    */
    rhs[ResidualIndex::mean_surface] =
        nonlinear_.state[surface_state_index(0)] +
        nonlinear_.state[surface_state_index(kFourierOrder)];
    for (int i = 1; i < kFourierOrder; ++i)
      rhs[ResidualIndex::mean_surface] =
          rhs[ResidualIndex::mean_surface] +
          nonlinear_.state[surface_state_index(i)] +
          nonlinear_.state[surface_state_index(i)];

    // Global equation R8: crest elevation minus trough elevation equals kH.
    rhs[ResidualIndex::wave_height] =
        nonlinear_.state[surface_state_index(0)] -
        nonlinear_.state[surface_state_index(kFourierOrder)] -
        nonlinear_.state[StateIndex::wave_height];

    /*
    Surface equations. Each collocation row is independent, so the outer loop is
    safe to parallelize. Modal sums remain serial within a row to preserve a
    deterministic long-double accumulation order.
    */
#if defined(_OPENMP)
#pragma omp parallel for schedule(static) if (kFourierOrder >=                   \
                                                  kOpenMpMinimumOrder)
#endif
    for (int m = 0; m <= kFourierOrder; ++m) {
      long double psi_acc = 0.0L;
      long double u_acc = 0.0L;
      long double v_acc = 0.0L;
      bool basis_ok = true;
      for (int j = 1; j <= kFourierOrder; ++j) {
        const double cosine = cosine_at_node(m, j);
        const double sine = sine_at_node(m, j);
        VerticalBasisValues basis;
        if (!evaluate_vertical_basis(
                j, nonlinear_.state[surface_state_index(m)],
                nonlinear_.state[StateIndex::depth], basis)) {
          basis_ok = false;
          break;
        }
        psi_acc += static_cast<long double>(
                       nonlinear_.state[coefficient_state_index(j)]) *
                   static_cast<long double>(basis.s) *
                   static_cast<long double>(cosine);
        u_acc += static_cast<long double>(j) *
                 static_cast<long double>(
                     nonlinear_.state[coefficient_state_index(j)]) *
                 static_cast<long double>(basis.c) *
                 static_cast<long double>(cosine);
        v_acc += static_cast<long double>(j) *
                 static_cast<long double>(
                     nonlinear_.state[coefficient_state_index(j)]) *
                 static_cast<long double>(basis.s) *
                 static_cast<long double>(sine);
      }
      if (!basis_ok) {
        rhs[kinematic_residual_index(m)] =
            std::numeric_limits<double>::infinity();
        rhs[dynamic_residual_index(m)] =
            std::numeric_limits<double>::infinity();
        continue;
      }

      // Constant-stream-function condition at surface node m.
      rhs[kinematic_residual_index(m)] = static_cast<double>(
          psi_acc -
          static_cast<long double>(nonlinear_.state[StateIndex::flux]) -
          static_cast<long double>(
              nonlinear_.state[StateIndex::wave_frame_speed]) *
              static_cast<long double>(
                  nonlinear_.state[surface_state_index(m)]));

      // Bernoulli condition at the same surface node.
      const long double wave_frame_u_acc =
          -static_cast<long double>(
              nonlinear_.state[StateIndex::wave_frame_speed]) +
          u_acc;
      rhs[dynamic_residual_index(m)] = static_cast<double>(
          0.5L * (wave_frame_u_acc * wave_frame_u_acc + v_acc * v_acc) +
          static_cast<long double>(nonlinear_.state[surface_state_index(m)]) -
          static_cast<long double>(nonlinear_.state[StateIndex::bernoulli]));
    }

    return stable_sum_squares(rhs);
  }

  /*
  Assemble the exact dense Jacobian of every residual with respect to every
  active unknown.

  The global rows are filled directly from their closed-form equations. Surface
  rows reuse the same vertical basis as evaluate_residuals and analytically
  differentiate with respect to depth, local surface elevation, global flow
  variables, and every Fourier coefficient.

  No finite differences are used. This removes differencing noise, permits
  quadratic local convergence, and makes the linear model reliable enough for
  trust-region acceptance tests.
  */
  void assemble_jacobian(DenseMatrix &jac) {
#if defined(_OPENMP)
#pragma omp parallel for schedule(static) if (kFourierOrder >=                   \
                                                  kOpenMpMinimumOrder)
#endif
    for (int i = 1; i <= kSystemSize; ++i) {
      for (int j = 1; j <= kSystemSize; ++j)
        jac[i][j] = 0.0;
    }

    // Derivatives of global equations R1 and R2.
    jac[ResidualIndex::depth][StateIndex::depth] =
        -continuation_.relative_height;
    jac[ResidualIndex::depth][StateIndex::wave_height] = 1.0;

    if (problem_.wave_specification != WaveSpecification::Period) {
      jac[ResidualIndex::wave_specification][StateIndex::wave_height] = 1.0;
    } else {
      jac[ResidualIndex::wave_specification][StateIndex::wave_height] = 1.0;
      jac[ResidualIndex::wave_specification][StateIndex::period] =
          -2.0 * continuation_.height_parameter *
          nonlinear_.state[StateIndex::period];
    }

    // Derivatives of global equations R3 through R6.
    jac[ResidualIndex::dispersion][StateIndex::period] =
        nonlinear_.state[StateIndex::celerity];
    jac[ResidualIndex::dispersion][StateIndex::celerity] =
        nonlinear_.state[StateIndex::period];

    jac[ResidualIndex::eulerian_current][StateIndex::celerity] = -1.0;
    jac[ResidualIndex::eulerian_current][StateIndex::eulerian_current] = 1.0;
    jac[ResidualIndex::eulerian_current][StateIndex::wave_frame_speed] = 1.0;

    jac[ResidualIndex::mass_transport_current][StateIndex::depth] =
        nonlinear_.state[StateIndex::mass_transport_current] +
        nonlinear_.state[StateIndex::wave_frame_speed] -
        nonlinear_.state[StateIndex::celerity];
    jac[ResidualIndex::mass_transport_current][StateIndex::celerity] =
        -nonlinear_.state[StateIndex::depth];
    jac[ResidualIndex::mass_transport_current]
       [StateIndex::mass_transport_current] =
           nonlinear_.state[StateIndex::depth];
    jac[ResidualIndex::mass_transport_current][StateIndex::wave_frame_speed] =
        nonlinear_.state[StateIndex::depth];
    jac[ResidualIndex::mass_transport_current][StateIndex::flux] = -1.0;

    jac[ResidualIndex::prescribed_current][selected_current_state_index()] =
        1.0;
    if (nonlinear_.state[StateIndex::depth] > 0.0) {
      jac[ResidualIndex::prescribed_current][StateIndex::depth] =
          -0.5 * problem_.current_value /
          std::sqrt(nonlinear_.state[StateIndex::depth]);
    }

    // Derivatives of mean-level and crest-to-trough constraints R7 and R8.
    jac[ResidualIndex::mean_surface][surface_state_index(0)] = 1.0;
    jac[ResidualIndex::mean_surface][surface_state_index(kFourierOrder)] = 1.0;
    for (int i = 1; i < kFourierOrder; ++i)
      jac[ResidualIndex::mean_surface][surface_state_index(i)] = 2.0;

    jac[ResidualIndex::wave_height][StateIndex::wave_height] = -1.0;
    jac[ResidualIndex::wave_height][surface_state_index(0)] = 1.0;
    jac[ResidualIndex::wave_height][surface_state_index(kFourierOrder)] = -1.0;

#if defined(_OPENMP)
#pragma omp parallel for schedule(static) if (kFourierOrder >=                   \
                                                  kOpenMpMinimumOrder)
#endif
    /*
    Each surface node contributes two independent Jacobian rows. Parallelizing
    over nodes is race-free because no two iterations write the same rows.
    */
    for (int m = 0; m <= kFourierOrder; ++m) {
      const int eta_index = surface_state_index(m);
      const int row_psi = kinematic_residual_index(m);
      const int row_dyn = dynamic_residual_index(m);
      const double eta = nonlinear_.state[eta_index];

      long double u_acc = 0.0L;
      long double v_acc = 0.0L;
      long double dpsi_dkd_acc = 0.0L;
      long double du_dkd_acc = 0.0L;
      long double dv_dkd_acc = 0.0L;
      long double du_deta_acc = 0.0L;
      long double dv_deta_acc = 0.0L;

      for (int j = 1; j <= kFourierOrder; ++j) {
        const int coeff_index = coefficient_state_index(j);
        const double bj = nonlinear_.state[coeff_index];
        const double jj = static_cast<double>(j);
        const double jj2 = jj * jj;
        const double cjx = cosine_at_node(m, j);
        const double sjx = sine_at_node(m, j);

        VerticalBasisValues basis;
        if (!evaluate_vertical_basis(
                j, eta, nonlinear_.state[StateIndex::depth], basis)) {
          jac[row_psi][eta_index] = std::numeric_limits<double>::quiet_NaN();
          jac[row_dyn][eta_index] = std::numeric_limits<double>::quiet_NaN();
          continue;
        }

        const double u_basis = jj * basis.c * cjx;
        const double v_basis = jj * basis.s * sjx;

        u_acc +=
            static_cast<long double>(bj) * static_cast<long double>(u_basis);
        v_acc +=
            static_cast<long double>(bj) * static_cast<long double>(v_basis);
        du_deta_acc += static_cast<long double>(jj2 * bj * cjx) *
                       static_cast<long double>(basis.s);
        dv_deta_acc += static_cast<long double>(jj2 * bj * sjx) *
                       static_cast<long double>(basis.c);

        dpsi_dkd_acc += static_cast<long double>(bj * cjx) *
                        static_cast<long double>(basis.ds_dkd);
        du_dkd_acc += static_cast<long double>(jj * bj * cjx) *
                      static_cast<long double>(basis.dc_dkd);
        dv_dkd_acc += static_cast<long double>(jj * bj * sjx) *
                      static_cast<long double>(basis.ds_dkd);

        /*
        Temporarily store S_j and d(u)/dB_j in the coefficient columns. Once the
        total velocity is known, these entries are converted in place to the
        final kinematic and dynamic derivatives. This removes a second
        exponential basis evaluation for every node/mode pair.
        */
        jac[row_psi][coeff_index] = basis.s;
        jac[row_dyn][coeff_index] = u_basis;
      }

      const double u_ = static_cast<double>(u_acc);
      const double dpsi_dkd = static_cast<double>(dpsi_dkd_acc);
      const double wave_frame_u =
          -nonlinear_.state[StateIndex::wave_frame_speed] + u_;

      jac[row_psi][StateIndex::depth] = dpsi_dkd;
      jac[row_psi][StateIndex::wave_frame_speed] = -eta;
      jac[row_psi][StateIndex::flux] = -1.0;
      jac[row_psi][eta_index] =
          u_ - nonlinear_.state[StateIndex::wave_frame_speed];

      jac[row_dyn][StateIndex::depth] = static_cast<double>(
          static_cast<long double>(wave_frame_u) * du_dkd_acc +
          v_acc * dv_dkd_acc);
      jac[row_dyn][StateIndex::wave_frame_speed] = -wave_frame_u;
      jac[row_dyn][StateIndex::bernoulli] = -1.0;
      jac[row_dyn][eta_index] = static_cast<double>(
          static_cast<long double>(wave_frame_u) * du_deta_acc +
          v_acc * dv_deta_acc + 1.0L);

      /*
      Finalize coefficient derivatives from the basis values cached during the
      first modal pass. No transcendental functions are evaluated here.
      */
      for (int j = 1; j <= kFourierOrder; ++j) {
        const int coeff_index = coefficient_state_index(j);
        const double jj = static_cast<double>(j);
        const double basis_s = jac[row_psi][coeff_index];
        const double du_db = jac[row_dyn][coeff_index];
        const double dv_db = jj * basis_s * sine_at_node(m, j);

        jac[row_psi][coeff_index] = basis_s * cosine_at_node(m, j);
        jac[row_dyn][coeff_index] =
            static_cast<double>(static_cast<long double>(wave_frame_u) *
                                    static_cast<long double>(du_db) +
                                v_acc * static_cast<long double>(dv_db));
      }
    }
  }

  /*
  ============================================================================
  STATE VALIDATION, NORMS, AND COORDINATE SCALING
  ============================================================================
  */

  // Return true only when every active residual is finite.
  [[nodiscard]] bool residual_is_finite(const Vector &residual) const {
    for (int i = 1; i <= kSystemSize; ++i) {
      if (!std::isfinite(residual[i]))
        return false;
    }
    return true;
  }

  /*
  Update RMS and maximum residual diagnostics from one evaluated residual
  vector. Centralizing this calculation keeps all failure and acceptance paths
  consistent.
  */
  void update_residual_diagnostics(const Vector &residual, double sum_squares) {
    if (!std::isfinite(sum_squares) || !residual_is_finite(residual)) {
      status_.residual_rms = std::numeric_limits<double>::infinity();
      status_.residual_max = std::numeric_limits<double>::infinity();
      return;
    }

    status_.residual_rms = std::sqrt(std::max(0.0, sum_squares) /
                                     static_cast<double>(kSystemSize));
    status_.residual_max = 0.0;
    for (int i = 1; i <= kSystemSize; ++i) {
      status_.residual_max =
          std::max(status_.residual_max, std::fabs(residual[i]));
    }
  }

  // Require both aggregate and worst-equation convergence. The maximum bound
  // prevents a single poorly satisfied boundary equation from being hidden by
  // an otherwise small RMS residual.
  [[nodiscard]] bool residual_converged(double tolerance) const noexcept {
    return status_.residual_rms <= tolerance &&
           status_.residual_max <= 10.0 * tolerance;
  }

  /*
  Validate a trial state before expensive residual acceptance tests.

  Positive global quantities must remain finite. In finite depth every surface
  node must stay above the bed by a scale-aware clearance. The largest modal
  exponent is also bounded so that all subsequent basis evaluations remain
  inside the safe floating-point range.
  */
  [[nodiscard]] bool state_is_admissible() const {
    for (int i = 1; i <= kSystemSize; ++i) {
      if (!std::isfinite(nonlinear_.state[i]))
        return false;
    }

    if (!(nonlinear_.state[StateIndex::depth] >
          64.0 * std::numeric_limits<double>::min())) {
      return false;
    }

    const double bed_clearance =
        64.0 * std::numeric_limits<double>::epsilon() *
        std::max(1.0, nonlinear_.state[StateIndex::depth]);
    for (int node = 0; node <= kFourierOrder; ++node) {
      if (!(nonlinear_.state[surface_state_index(node)] >
            -nonlinear_.state[StateIndex::depth] + bed_clearance)) {
        return false;
      }
    }

    if (!(nonlinear_.state[StateIndex::wave_height] > 0.0 &&
          nonlinear_.state[StateIndex::period] > 0.0 &&
          nonlinear_.state[StateIndex::celerity] > 0.0))
      return false;

    const double kSafeExpArgument =
        0.95 * std::log(std::numeric_limits<double>::max());
    for (int m = 0; m <= kFourierOrder; ++m) {
      if (static_cast<double>(kFourierOrder) *
              nonlinear_.state[surface_state_index(m)] >=
          kSafeExpArgument)
        return false;
    }

    return true;
  }

  // Reject a Jacobian immediately if any active entry is NaN or infinite.
  [[nodiscard]] bool jacobian_is_finite(const DenseMatrix &jac) const {
    for (int i = 1; i <= kSystemSize; ++i) {
      for (int j = 1; j <= kSystemSize; ++j) {
        if (!std::isfinite(jac[i][j]))
          return false;
      }
    }
    return true;
  }

  /*
  Stable Euclidean norm of the active one-based vector. The implementation uses
  the same scale/ssq recurrence as stable_sum_squares and therefore remains safe
  for highly disparate component magnitudes.
  */
  [[nodiscard]] double euclidean_norm(const Vector &vector) const {
    double scale = 0.0;
    double ssq = 1.0;
    for (int i = 1; i <= kSystemSize; ++i) {
      const double value = std::fabs(vector[i]);
      if (value == 0.0)
        continue;
      if (scale < value) {
        const double ratio = scale / value;
        ssq = 1.0 + ssq * ratio * ratio;
        scale = value;
      } else {
        const double ratio = value / scale;
        ssq += ratio * ratio;
      }
    }
    return scale == 0.0 ? 0.0 : scale * std::sqrt(ssq);
  }

  // Root-mean-square norm used for dimension-independent step diagnostics.
  [[nodiscard]] double rms_norm(const Vector &vector) const {
    return euclidean_norm(vector) /
           std::sqrt(static_cast<double>(kSystemSize));
  }

  // Long-double accumulated dot product for trust-region model calculations.
  [[nodiscard]] double dot_product(const Vector &a, const Vector &b) const {
    long double sum_ = 0.0L;
    for (int i = 1; i <= kSystemSize; ++i) {
      sum_ += static_cast<long double>(a[i]) * static_cast<long double>(b[i]);
    }
    return static_cast<double>(sum_);
  }

  // Identify strictly positive global variables updated in logarithmic rather
  // than additive coordinates by apply_retraction().
  [[nodiscard]] static bool uses_log_coordinate(int index) noexcept {
    return index >= StateIndex::depth && index <= StateIndex::celerity;
  }

  /*
  Map a scaled optimization step back to physical state variables.

  Depth, height, period, and celerity are updated multiplicatively in
  logarithmic coordinates, which preserves strict positivity for every finite
  step. Remaining variables use additive scaled coordinates. This retraction
  eliminates a common source of invalid Newton trials without clipping or
  changing the local derivative at zero step.
  */
  [[nodiscard]] bool apply_retraction(const Vector &base_state,
                                      const Vector &scaled_step) {
    for (int i = 1; i <= kSystemSize; ++i) {
      if (uses_log_coordinate(i)) {
        const double base = base_state[i];
        if (!(base > 0.0) || !std::isfinite(base))
          return false;
        const double exponent =
            (workspace_.column_scale[i] / base) * scaled_step[i];
        if (!std::isfinite(exponent) || std::fabs(exponent) > 60.0)
          return false;
        nonlinear_.state[i] = base * std::exp(exponent);
      } else {
        nonlinear_.state[i] =
            std::fma(workspace_.column_scale[i], scaled_step[i], base_state[i]);
      }
      if (!std::isfinite(nonlinear_.state[i]))
        return false;
    }
    return true;
  }

  // Advance one positive state component along a continuation tangent using a
  // multiplicative predictor. Invalid or excessive exponents safely retain the
  // previous value, allowing the nonlinear corrector to recover.
  [[nodiscard]] static double log_predict(double previous, double tangent,
                                          double delta) noexcept {
    if (!(previous > 0.0) || !std::isfinite(previous))
      return previous;
    const double exponent = delta * tangent / previous;
    if (!std::isfinite(exponent) || std::fabs(exponent) > 40.0)
      return previous;
    return previous * std::exp(exponent);
  }

  /*
  Select a characteristic scale for one state component.

  Global variables use magnitude-aware scales. Surface ordinates use wave-height
  scale. Modal coefficients additionally decay with mode number. These scales
  equilibrate heterogeneous physical quantities before linear solution and
  define the dimensionless trust-region norm.
  */
  [[nodiscard]] double state_scale(int index, const Vector &state) const {
    const double magnitude = std::fabs(state[index]);

    if (index == StateIndex::depth)
      return std::max(5.0e-2, magnitude);
    if (index == StateIndex::wave_height)
      return std::max(1.0e-3, magnitude);
    if (index == StateIndex::period || index == StateIndex::celerity)
      return std::max(1.0e-1, magnitude);
    if (index >= StateIndex::eulerian_current && index <= StateIndex::bernoulli)
      return std::max(1.0e-3, magnitude);

    const double wave_height_scale = std::fabs(state[StateIndex::wave_height]);
    if (index >= StateIndex::surface_begin &&
        index <= surface_state_index(kFourierOrder)) {
      return std::max({1.0e-4, wave_height_scale, magnitude});
    }

    const int mode = index - surface_state_index(kFourierOrder);
    const double spectral_scale =
        wave_height_scale / (2.0 * static_cast<double>(std::max(1, mode)));
    return std::max({1.0e-7, spectral_scale, magnitude});
  }

  /*
  Equilibrate the Jacobian and residual.

  Column scales express physically meaningful step sizes. Row scales normalize
  equation magnitudes. A second column adjustment reduces remaining RMS column
  imbalance. The resulting system has the same exact Newton solution as the
  original square system but substantially better numerical conditioning.
  */
  void build_scaled_linearization(const DenseMatrix &jac,
                                  const Vector &residual) {
    for (int j = 1; j <= kSystemSize; ++j) {
      workspace_.column_scale[j] = state_scale(j, workspace_.base_state);
    }

    for (int pass = 0; pass < 2; ++pass) {
      for (int i = 1; i <= kSystemSize; ++i) {
        double row_norm = 0.0;
        for (int j = 1; j <= kSystemSize; ++j) {
          row_norm = std::max(
              row_norm, std::fabs(jac[i][j] * workspace_.column_scale[j]));
        }
        workspace_.row_scale[i] = 1.0 / std::max(row_norm, 1.0e-14);
      }

      if (pass == 0) {
        for (int j = 1; j <= kSystemSize; ++j) {
          long double column_sum = 0.0L;
          for (int i = 1; i <= kSystemSize; ++i) {
            const long double value =
                static_cast<long double>(workspace_.row_scale[i]) *
                static_cast<long double>(jac[i][j]) *
                static_cast<long double>(workspace_.column_scale[j]);
            column_sum += value * value;
          }
          const double column_rms = std::sqrt(static_cast<double>(
              column_sum / static_cast<long double>(kSystemSize)));
          if (column_rms > 0.0 && std::isfinite(column_rms)) {
            const double adjustment = std::clamp(1.0 / column_rms, 0.25, 4.0);
            workspace_.column_scale[j] *= adjustment;
          }
        }
      }
    }

    for (int i = 1; i <= kSystemSize; ++i) {
      workspace_.linear_rhs[i] = -workspace_.row_scale[i] * residual[i];
      for (int j = 1; j <= kSystemSize; ++j) {
        workspace_.linear_matrix[i][j] =
            workspace_.row_scale[i] * jac[i][j] * workspace_.column_scale[j];
      }
    }
  }

  // Multiply by the fully row- and column-scaled Jacobian. This operation is
  // used for Cauchy directions and scaled linear residual verification; long-
  // double accumulation reduces cancellation without changing stored precision.
  void scaled_matrix_vector_product(const Vector &x, Vector &result) const {
    for (int i = 1; i <= kSystemSize; ++i) {
      long double sum_ = 0.0L;
      for (int j = 1; j <= kSystemSize; ++j) {
        sum_ += static_cast<long double>(workspace_.row_scale[i]) *
                static_cast<long double>(workspace_.jacobian[i][j]) *
                static_cast<long double>(workspace_.column_scale[j]) *
                static_cast<long double>(x[j]);
      }
      result[i] = static_cast<double>(sum_);
    }
  }

  /*
  Multiply a scaled-coordinate step by the physical Jacobian. Row scaling is not
  applied because the result is used in the physical residual-space merit model.
  */
  void physical_matrix_vector_product(const Vector &x, Vector &result) const {
    for (int i = 1; i <= kSystemSize; ++i) {
      long double sum_ = 0.0L;
      for (int j = 1; j <= kSystemSize; ++j) {
        sum_ += static_cast<long double>(workspace_.jacobian[i][j]) *
                static_cast<long double>(workspace_.column_scale[j]) *
                static_cast<long double>(x[j]);
      }
      result[i] = static_cast<double>(sum_);
    }
  }

  /*
  ============================================================================
  DENSE LINEAR ALGEBRA
  ============================================================================

  The primary Newton system is dense and moderate in size. A cache-friendly
  pivoted LU factorization is therefore faster than sparse or iterative methods
  for the supported Fourier orders. Rank-revealing QR and regularized normal
  equations provide progressively more defensive fallbacks.
  ============================================================================
  */

  /*
  In-place LU factorization with partial row pivoting.

  A scale-dependent pivot floor detects singular or numerically unusable
  factorizations before division. Row swaps are recorded for all subsequent
  right-hand sides and iterative-refinement corrections.
  */
  [[nodiscard]] bool factorize_lu(DenseMatrix &matrix, IndexVector &pivot) {
    pivot.assign(static_cast<std::size_t>(kSystemSize + 1), 0);
    double matrix_max = 0.0;
    for (int i = 1; i <= kSystemSize; ++i) {
      for (int j = 1; j <= kSystemSize; ++j) {
        matrix_max = std::max(matrix_max, std::fabs(matrix[i][j]));
      }
    }
    if (!(matrix_max > 0.0) || !std::isfinite(matrix_max))
      return false;

    const double pivot_floor = 64.0 * std::numeric_limits<double>::epsilon() *
                               static_cast<double>(kSystemSize) *
                               matrix_max;
    for (int k = 1; k <= kSystemSize; ++k) {
      int pivot_row = k;
      double pivot_abs = std::fabs(matrix[k][k]);
      for (int i = k + 1; i <= kSystemSize; ++i) {
        const double candidate = std::fabs(matrix[i][k]);
        if (candidate > pivot_abs) {
          pivot_abs = candidate;
          pivot_row = i;
        }
      }

      if (!(pivot_abs > pivot_floor) || !std::isfinite(pivot_abs))
        return false;

      pivot[k] = pivot_row;
      if (pivot_row != k) {
        for (int j = 1; j <= kSystemSize; ++j) {
          std::swap(matrix[k][j], matrix[pivot_row][j]);
        }
      }

      const double inverse_pivot = 1.0 / matrix[k][k];
      for (int i = k + 1; i <= kSystemSize; ++i) {
        matrix[i][k] *= inverse_pivot;
        const double multiplier = matrix[i][k];
        for (int j = k + 1; j <= kSystemSize; ++j) {
          matrix[i][j] = std::fma(-multiplier, matrix[k][j], matrix[i][j]);
        }
      }
    }

    return true;
  }

  /*
  Solve an already factorized LU system.

  Forward and backward substitutions accumulate in long double. This improves
  the correction quality at negligible cost relative to factorization and
  supports the iterative-refinement loop used for difficult Jacobians.
  */
  [[nodiscard]] bool solve_lu(const DenseMatrix &lu, const IndexVector &pivot,
                              const Vector &right_hand_side, Vector &solution) {
    for (int i = 1; i <= kSystemSize; ++i)
      solution[i] = right_hand_side[i];

    for (int k = 1; k <= kSystemSize; ++k) {
      if (pivot[k] != k)
        std::swap(solution[k], solution[pivot[k]]);
    }

    for (int i = 2; i <= kSystemSize; ++i) {
      long double value = static_cast<long double>(solution[i]);
      for (int j = 1; j < i; ++j) {
        value -= static_cast<long double>(lu[i][j]) *
                 static_cast<long double>(solution[j]);
      }
      solution[i] = static_cast<double>(value);
    }

    for (int i = kSystemSize; i >= 1; --i) {
      long double value = static_cast<long double>(solution[i]);
      for (int j = i + 1; j <= kSystemSize; ++j) {
        value -= static_cast<long double>(lu[i][j]) *
                 static_cast<long double>(solution[j]);
      }
      const double diagonal = lu[i][i];
      if (diagonal == 0.0 || !std::isfinite(diagonal))
        return false;
      solution[i] =
          static_cast<double>(value / static_cast<long double>(diagonal));
      if (!std::isfinite(solution[i]))
        return false;
    }

    return true;
  }

  // Validate every active component of a linear-system solution or search
  // direction before it participates in trust-region model arithmetic.
  [[nodiscard]] bool vector_is_finite(const Vector &vector) const {
    for (int i = 1; i <= kSystemSize; ++i) {
      if (!std::isfinite(vector[i]))
        return false;
    }
    return true;
  }

  /*
  Reconstruct the scaled Jacobian after LU has overwritten the working matrix.
  Fallback solvers require the original equilibrated coefficients.
  */
  void rebuild_scaled_matrix() {
    for (int i = 1; i <= kSystemSize; ++i) {
      for (int j = 1; j <= kSystemSize; ++j) {
        workspace_.linear_matrix[i][j] = workspace_.row_scale[i] *
                                         workspace_.jacobian[i][j] *
                                         workspace_.column_scale[j];
      }
    }
  }

  /*
  Column-pivoted Householder QR fallback.

  Column pivoting orders directions by remaining norm and exposes numerical
  rank. Reflectors and transformed right-hand sides are accumulated in long
  double. Rank-deficient trailing directions are assigned zero, producing a
  stable minimum-norm correction in the resolved subspace without forming J^T J.
  */
  [[nodiscard]] bool solve_rank_revealing_qr(DenseMatrix &matrix,
                                             const Vector &right_hand_side,
                                             Vector &solution) {
    const int dimension = kSystemSize;
    const std::size_t stride = static_cast<std::size_t>(dimension);
    std::vector<long double> qr(stride * stride, 0.0L);
    std::vector<long double> transformed_rhs(stride, 0.0L);
    std::vector<long double> householder(stride, 0.0L);
    std::vector<long double> column_norm_squared(stride, 0.0L);
    std::vector<long double> permuted_solution(stride, 0.0L);
    std::vector<int> permutation(stride, 0);

    auto at = [&](int row, int column) -> long double & {
      return qr[static_cast<std::size_t>(row) * stride +
                static_cast<std::size_t>(column)];
    };

    long double maximum_column_norm = 0.0L;
    for (int column = 0; column < dimension; ++column) {
      permutation[static_cast<std::size_t>(column)] = column;
      long double norm_squared = 0.0L;
      for (int row = 0; row < dimension; ++row) {
        const long double value =
            static_cast<long double>(matrix[row + 1][column + 1]);
        at(row, column) = value;
        norm_squared += value * value;
      }
      column_norm_squared[static_cast<std::size_t>(column)] = norm_squared;
      maximum_column_norm =
          std::max(maximum_column_norm, std::sqrt(norm_squared));
    }
    for (int row = 0; row < dimension; ++row) {
      transformed_rhs[static_cast<std::size_t>(row)] =
          static_cast<long double>(right_hand_side[row + 1]);
    }

    if (!(maximum_column_norm > 0.0L) || !std::isfinite(maximum_column_norm)) {
      return false;
    }

    const long double rank_threshold =
        64.0L *
        static_cast<long double>(std::numeric_limits<double>::epsilon()) *
        static_cast<long double>(dimension) * maximum_column_norm;

    int numerical_rank = 0;
    for (int k = 0; k < dimension; ++k) {
      int pivot_column = k;
      long double pivot_norm_squared =
          column_norm_squared[static_cast<std::size_t>(k)];
      for (int column = k + 1; column < dimension; ++column) {
        const long double candidate =
            column_norm_squared[static_cast<std::size_t>(column)];
        if (candidate > pivot_norm_squared) {
          pivot_norm_squared = candidate;
          pivot_column = column;
        }
      }

      if (pivot_column != k) {
        for (int row = 0; row < dimension; ++row) {
          std::swap(at(row, k), at(row, pivot_column));
        }
        std::swap(column_norm_squared[static_cast<std::size_t>(k)],
                  column_norm_squared[static_cast<std::size_t>(pivot_column)]);
        std::swap(permutation[static_cast<std::size_t>(k)],
                  permutation[static_cast<std::size_t>(pivot_column)]);
      }

      long double column_norm = 0.0L;
      for (int row = k; row < dimension; ++row) {
        column_norm = std::hypot(column_norm, at(row, k));
      }
      if (!(column_norm > rank_threshold) || !std::isfinite(column_norm))
        break;

      const long double diagonal = at(k, k);
      const long double reflector_diagonal =
          -std::copysign(column_norm, diagonal == 0.0L ? 1.0L : diagonal);

      long double reflector_norm = 0.0L;
      householder[static_cast<std::size_t>(k)] = diagonal - reflector_diagonal;
      reflector_norm =
          std::hypot(reflector_norm, householder[static_cast<std::size_t>(k)]);
      for (int row = k + 1; row < dimension; ++row) {
        householder[static_cast<std::size_t>(row)] = at(row, k);
        reflector_norm = std::hypot(reflector_norm,
                                    householder[static_cast<std::size_t>(row)]);
      }
      if (!(reflector_norm > 0.0L) || !std::isfinite(reflector_norm))
        break;

      for (int row = k; row < dimension; ++row) {
        householder[static_cast<std::size_t>(row)] /= reflector_norm;
      }

      for (int column = k; column < dimension; ++column) {
        long double projection = 0.0L;
        for (int row = k; row < dimension; ++row) {
          projection +=
              householder[static_cast<std::size_t>(row)] * at(row, column);
        }
        projection *= 2.0L;
        for (int row = k; row < dimension; ++row) {
          at(row, column) -=
              projection * householder[static_cast<std::size_t>(row)];
        }
      }

      long double rhs_projection = 0.0L;
      for (int row = k; row < dimension; ++row) {
        rhs_projection += householder[static_cast<std::size_t>(row)] *
                          transformed_rhs[static_cast<std::size_t>(row)];
      }
      rhs_projection *= 2.0L;
      for (int row = k; row < dimension; ++row) {
        transformed_rhs[static_cast<std::size_t>(row)] -=
            rhs_projection * householder[static_cast<std::size_t>(row)];
      }

      at(k, k) = reflector_diagonal;
      for (int row = k + 1; row < dimension; ++row)
        at(row, k) = 0.0L;
      ++numerical_rank;

      for (int column = k + 1; column < dimension; ++column) {
        long double norm_squared = 0.0L;
        for (int row = k + 1; row < dimension; ++row) {
          const long double value = at(row, column);
          norm_squared += value * value;
        }
        column_norm_squared[static_cast<std::size_t>(column)] = norm_squared;
      }
    }

    if (numerical_rank <= 0)
      return false;

    for (int row = numerical_rank - 1; row >= 0; --row) {
      long double value = transformed_rhs[static_cast<std::size_t>(row)];
      for (int column = row + 1; column < numerical_rank; ++column) {
        value -= at(row, column) *
                 permuted_solution[static_cast<std::size_t>(column)];
      }
      const long double diagonal = at(row, row);
      if (!(std::fabs(diagonal) > rank_threshold) || !std::isfinite(diagonal)) {
        return false;
      }
      permuted_solution[static_cast<std::size_t>(row)] = value / diagonal;
    }

    for (int i = 1; i <= dimension; ++i)
      solution[i] = 0.0;
    for (int column = 0; column < dimension; ++column) {
      const int original_column = permutation[static_cast<std::size_t>(column)];
      solution[original_column + 1] = static_cast<double>(
          permuted_solution[static_cast<std::size_t>(column)]);
    }

    return vector_is_finite(solution);
  }

  /*
  Final linear fallback using adaptive Tikhonov regularization.

  The routine forms J^T J + lambda I in long double and solves it by Cholesky
  factorization. Lambda is increased until the matrix is positive definite and
  the resulting direction is a descent direction. This path is intentionally
  used only when both LU and rank-revealing QR are inadequate.
  */
  [[nodiscard]] bool solve_regularized_least_squares(
      DenseMatrix &matrix, const Vector &right_hand_side, Vector &solution) {
    const int dimension = kSystemSize;
    const std::size_t stride = static_cast<std::size_t>(dimension);
    std::vector<long double> normal(stride * stride, 0.0L);
    std::vector<long double> gradient(stride, 0.0L);
    std::vector<long double> factor(stride * stride, 0.0L);
    std::vector<long double> work(stride, 0.0L);

    auto at = [&](std::vector<long double> &a, int row,
                  int column) -> long double & {
      return a[static_cast<std::size_t>(row) * stride +
               static_cast<std::size_t>(column)];
    };

    long double maximum_diagonal = 0.0L;
    for (int j = 0; j < dimension; ++j) {
      long double rhs_sum = 0.0L;
      for (int i = 0; i < dimension; ++i) {
        rhs_sum += static_cast<long double>(matrix[i + 1][j + 1]) *
                   static_cast<long double>(right_hand_side[i + 1]);
      }
      gradient[static_cast<std::size_t>(j)] = rhs_sum;

      for (int k = j; k < dimension; ++k) {
        long double sum_ = 0.0L;
        for (int i = 0; i < dimension; ++i) {
          sum_ += static_cast<long double>(matrix[i + 1][j + 1]) *
                  static_cast<long double>(matrix[i + 1][k + 1]);
        }
        at(normal, j, k) = sum_;
        at(normal, k, j) = sum_;
      }
      maximum_diagonal = std::max(maximum_diagonal, at(normal, j, j));
    }

    if (!(maximum_diagonal > 0.0L) || !std::isfinite(maximum_diagonal)) {
      return false;
    }

    long double lambda =
        64.0L *
        static_cast<long double>(std::numeric_limits<double>::epsilon()) *
        static_cast<long double>(dimension) * maximum_diagonal;
    lambda = std::max(lambda, 1.0e-28L * maximum_diagonal);

    for (int attempt = 0; attempt < 12; ++attempt) {
      factor = normal;
      for (int i = 0; i < dimension; ++i)
        at(factor, i, i) += lambda;

      bool positive_definite = true;
      for (int j = 0; j < dimension && positive_definite; ++j) {
        for (int k = 0; k <= j; ++k) {
          long double value = at(factor, j, k);
          for (int q = 0; q < k; ++q) {
            value -= at(factor, j, q) * at(factor, k, q);
          }
          if (j == k) {
            if (!(value > 0.0L) || !std::isfinite(value)) {
              positive_definite = false;
              break;
            }
            at(factor, j, k) = std::sqrt(value);
          } else {
            at(factor, j, k) = value / at(factor, k, k);
          }
        }
        for (int k = j + 1; k < dimension; ++k)
          at(factor, j, k) = 0.0L;
      }

      if (!positive_definite) {
        lambda *= 100.0L;
        continue;
      }

      for (int i = 0; i < dimension; ++i) {
        long double value = gradient[static_cast<std::size_t>(i)];
        for (int j = 0; j < i; ++j)
          value -= at(factor, i, j) * work[static_cast<std::size_t>(j)];
        work[static_cast<std::size_t>(i)] = value / at(factor, i, i);
      }
      for (int i = dimension - 1; i >= 0; --i) {
        long double value = work[static_cast<std::size_t>(i)];
        for (int j = i + 1; j < dimension; ++j) {
          value -= at(factor, j, i) * static_cast<long double>(solution[j + 1]);
        }
        solution[i + 1] = static_cast<double>(value / at(factor, i, i));
      }

      if (!vector_is_finite(solution)) {
        lambda *= 100.0L;
        continue;
      }

      long double directional_derivative = 0.0L;
      for (int i = 0; i < dimension; ++i) {
        directional_derivative -= gradient[static_cast<std::size_t>(i)] *
                                  static_cast<long double>(solution[i + 1]);
      }
      if (directional_derivative < 0.0L)
        return true;

      lambda *= 100.0L;
    }

    return false;
  }

  /*
  Solve the equilibrated Newton system through a tiered strategy.

  1. Pivoted LU for speed.
  2. Long-double iterative refinement against the unmodified Jacobian.
  3. Column-pivoted QR if LU is singular or insufficiently accurate.
  4. Regularized least squares as a last-resort descent direction.

  Every accepted direction is checked by its relative linear residual.
  */
  [[nodiscard]] bool solve_scaled_newton_system() {
    status_.linear_solver = LinearSolverKind::LU;

    bool lu_success =
        factorize_lu(workspace_.linear_matrix, workspace_.pivot) &&
        solve_lu(workspace_.linear_matrix, workspace_.pivot,
                 workspace_.linear_rhs, workspace_.newton_direction);

    if (lu_success) {
      for (int refinement = 0; refinement < 4; ++refinement) {
        scaled_matrix_vector_product(workspace_.newton_direction,
                                     workspace_.model_residual);
        for (int i = 1; i <= kSystemSize; ++i) {
          workspace_.model_residual[i] =
              workspace_.linear_rhs[i] - workspace_.model_residual[i];
        }

        const double residual_rms = rms_norm(workspace_.model_residual);
        const double rhs_rms = rms_norm(workspace_.linear_rhs);
        if (residual_rms <= 32.0 * std::numeric_limits<double>::epsilon() *
                                static_cast<double>(kSystemSize) *
                                (1.0 + rhs_rms)) {
          break;
        }

        if (!solve_lu(workspace_.linear_matrix, workspace_.pivot,
                      workspace_.model_residual, workspace_.cauchy_direction)) {
          lu_success = false;
          break;
        }
        for (int i = 1; i <= kSystemSize; ++i) {
          workspace_.newton_direction[i] += workspace_.cauchy_direction[i];
        }
      }

      if (lu_success) {
        scaled_matrix_vector_product(workspace_.newton_direction,
                                     workspace_.model_residual);
        for (int i = 1; i <= kSystemSize; ++i) {
          workspace_.model_residual[i] =
              workspace_.linear_rhs[i] - workspace_.model_residual[i];
        }
        const double relative_residual =
            rms_norm(workspace_.model_residual) /
            (1.0 + rms_norm(workspace_.linear_rhs));
        lu_success = vector_is_finite(workspace_.newton_direction) &&
                     relative_residual <= 1.0e-9;
      }
    }

    if (lu_success)
      return true;

    rebuild_scaled_matrix();
    if (solve_rank_revealing_qr(workspace_.linear_matrix, workspace_.linear_rhs,
                                workspace_.newton_direction)) {
      scaled_matrix_vector_product(workspace_.newton_direction,
                                   workspace_.model_residual);
      for (int i = 1; i <= kSystemSize; ++i) {
        workspace_.model_residual[i] =
            workspace_.linear_rhs[i] - workspace_.model_residual[i];
      }
      const double relative_residual = rms_norm(workspace_.model_residual) /
                                       (1.0 + rms_norm(workspace_.linear_rhs));
      if (relative_residual <= 1.0e-8) {
        status_.linear_solver = LinearSolverKind::QR;
        return true;
      }
    }

    status_.linear_solver = LinearSolverKind::Regularized;
    rebuild_scaled_matrix();
    return solve_regularized_least_squares(workspace_.linear_matrix,
                                           workspace_.linear_rhs,
                                           workspace_.newton_direction);
  }

  /*
  ============================================================================
  NONLINEAR GLOBALIZATION AND CONTINUATION
  ============================================================================
  */

  /*
  Assemble the explicit derivative of the residual with respect to normalized
  continuation fraction. Solving J dz/df = -dF/df gives the local
  solution-branch tangent used by the predictor.
  */
  void assemble_continuation_derivative(Vector &derivative) {
    for (int i = 1; i <= kSystemSize; ++i)
      derivative[i] = 0.0;

    derivative[ResidualIndex::depth] =
        -nonlinear_.state[StateIndex::depth] * problem_.target_relative_height;
    derivative[2] = problem_.wave_specification == WaveSpecification::Period
                        ? -problem_.target_height_parameter *
                              nonlinear_.state[StateIndex::period] *
                              nonlinear_.state[StateIndex::period]
                        : -2.0 * kPi * problem_.target_height_parameter;
  }

  /*
  Compute and validate the tangent to the nonlinear solution branch.

  The same exact Jacobian, equilibration, and robust linear solver used by
  Newton are reused here. The final tangent is transformed back from scaled
  coordinates and rejected when its normalized RMS indicates a near-singular or
  unreliable branch direction.
  */
  [[nodiscard]] bool compute_continuation_tangent(Vector &tangent) {
    tangent.assign(static_cast<std::size_t>(kSystemSize + 1), 0.0);
    for (int i = 1; i <= kSystemSize; ++i)
      workspace_.base_state[i] = nonlinear_.state[i];

    assemble_jacobian(workspace_.jacobian);
    if (!jacobian_is_finite(workspace_.jacobian))
      return false;

    assemble_continuation_derivative(nonlinear_.trial_residual);
    build_scaled_linearization(workspace_.jacobian, nonlinear_.trial_residual);
    if (!solve_scaled_newton_system())
      return false;

    long double scaled_norm_squared = 0.0L;
    for (int i = 1; i <= kSystemSize; ++i) {
      tangent[i] = workspace_.column_scale[i] * workspace_.newton_direction[i];
      const long double scale = static_cast<long double>(
          std::max(1.0e-12, state_scale(i, nonlinear_.state)));
      const long double value = static_cast<long double>(tangent[i]) / scale;
      scaled_norm_squared += value * value;
    }

    const double scaled_rms = static_cast<double>(std::sqrt(
        scaled_norm_squared / static_cast<long double>(kSystemSize)));
    return std::isfinite(scaled_rms) && scaled_rms < 1.0e6;
  }

  /*
  Build the steepest-descent minimizer of the linearized least-squares model.

  This direction guarantees descent when the Newton direction lies outside the
  trust region or is unreliable. It forms the first segment of the Powell
  dogleg path.
  */
  [[nodiscard]] bool build_cauchy_direction() {
    for (int j = 1; j <= kSystemSize; ++j) {
      long double gradient = 0.0L;
      for (int i = 1; i <= kSystemSize; ++i) {
        gradient += static_cast<long double>(workspace_.jacobian[i][j]) *
                    static_cast<long double>(workspace_.column_scale[j]) *
                    static_cast<long double>(nonlinear_.residual[i]);
      }
      workspace_.cauchy_direction[j] = static_cast<double>(gradient);
    }

    const double gradient_squared =
        dot_product(workspace_.cauchy_direction, workspace_.cauchy_direction);
    if (!(gradient_squared > 0.0) || !std::isfinite(gradient_squared))
      return false;

    physical_matrix_vector_product(workspace_.cauchy_direction,
                                   workspace_.model_residual);
    const double ag_squared =
        dot_product(workspace_.model_residual, workspace_.model_residual);
    if (!(ag_squared > 0.0) || !std::isfinite(ag_squared))
      return false;

    const double alpha = gradient_squared / ag_squared;
    for (int j = 1; j <= kSystemSize; ++j) {
      workspace_.cauchy_direction[j] = -alpha * workspace_.cauchy_direction[j];
    }
    return vector_is_finite(workspace_.cauchy_direction);
  }

  /*
  Select a step on the Powell dogleg path.

  The choice is:
    - full Newton when it fits inside the radius;
    - a radius-scaled Cauchy step when steepest descent reaches the boundary;
    - the intersection of the Cauchy-to-Newton segment with the boundary
      otherwise.
  */
  [[nodiscard]] bool build_dogleg_step(double trust_radius) {
    const double target_norm =
        trust_radius * std::sqrt(static_cast<double>(kSystemSize));
    const double newton_norm = euclidean_norm(workspace_.newton_direction);
    const double cauchy_norm = euclidean_norm(workspace_.cauchy_direction);

    if (!std::isfinite(newton_norm) || !std::isfinite(cauchy_norm) ||
        !(target_norm > 0.0))
      return false;

    if (newton_norm <= target_norm) {
      for (int i = 1; i <= kSystemSize; ++i)
        workspace_.candidate_step[i] = workspace_.newton_direction[i];
      return true;
    }

    if (!(cauchy_norm > 0.0))
      return false;
    if (cauchy_norm >= target_norm) {
      const double factor = target_norm / cauchy_norm;
      for (int i = 1; i <= kSystemSize; ++i) {
        workspace_.candidate_step[i] = factor * workspace_.cauchy_direction[i];
      }
      return true;
    }

    long double direction_squared = 0.0L;
    long double cauchy_dot_direction = 0.0L;
    for (int i = 1; i <= kSystemSize; ++i) {
      const long double direction =
          static_cast<long double>(workspace_.newton_direction[i]) -
          static_cast<long double>(workspace_.cauchy_direction[i]);
      direction_squared += direction * direction;
      cauchy_dot_direction +=
          static_cast<long double>(workspace_.cauchy_direction[i]) * direction;
    }

    if (!(direction_squared > 0.0L))
      return false;
    const long double cauchy_squared = static_cast<long double>(cauchy_norm) *
                                       static_cast<long double>(cauchy_norm);
    const long double target_squared = static_cast<long double>(target_norm) *
                                       static_cast<long double>(target_norm);
    const long double discriminant =
        cauchy_dot_direction * cauchy_dot_direction +
        direction_squared * (target_squared - cauchy_squared);
    if (!(discriminant >= 0.0L))
      return false;

    const double tau = static_cast<double>(
        (-cauchy_dot_direction + std::sqrt(discriminant)) / direction_squared);
    if (!(tau >= 0.0 && tau <= 1.0) || !std::isfinite(tau))
      return false;

    for (int i = 1; i <= kSystemSize; ++i) {
      workspace_.candidate_step[i] = workspace_.cauchy_direction[i] +
                                     tau * (workspace_.newton_direction[i] -
                                            workspace_.cauchy_direction[i]);
    }
    return vector_is_finite(workspace_.candidate_step);
  }

  /*
  ============================================================================
  TRUST-REGION STEP
  ============================================================================

  One call performs one globally safeguarded nonlinear iteration.

  1. Evaluate the residual at the current state.
  2. Assemble and equilibrate the exact Jacobian.
  3. Compute a Newton direction and a Cauchy descent direction.
  4. Select a Powell-dogleg step inside the current trust region.
  5. Retract positive variables through logarithmic coordinates.
  6. Accept the trial only when the actual merit reduction agrees sufficiently
     with the local linear model.
  7. Expand or contract the trust radius according to model quality.

  The returned object contains the acceptance flag and the RMS norm of the
  accepted step in scaled coordinates. Residual diagnostics remain in status_
  because they are displayed and reused by the nonlinear convergence controller.
  ============================================================================
  */
  struct TrustRegionStepResult {
    bool accepted = false;
    double scaled_step_rms = std::numeric_limits<double>::infinity();
  };

  TrustRegionStepResult trust_region_iteration(int iteration) {
    constexpr double kMinimumTrustRadius = 1.0e-12;
    constexpr double kMaximumTrustRadius = 100.0;
    constexpr double kAcceptanceRatio = 1.0e-4;
    constexpr int kMaximumTrustTrials = 18;

    if (iteration <= 1)
      status_.trust_radius = 1.0;

    const double base_sum_squares = evaluate_residuals(nonlinear_.residual);
    if (!std::isfinite(base_sum_squares) ||
        !residual_is_finite(nonlinear_.residual)) {
      update_residual_diagnostics(nonlinear_.residual, base_sum_squares);
      return {};
    }

    // Preserve an exact rollback point before any trial state is formed.
    for (int i = 1; i <= kSystemSize; ++i)
      workspace_.base_state[i] = nonlinear_.state[i];

    assemble_jacobian(workspace_.jacobian);
    if (!jacobian_is_finite(workspace_.jacobian)) {
      update_residual_diagnostics(nonlinear_.residual, base_sum_squares);
      return {};
    }

    build_scaled_linearization(workspace_.jacobian, nonlinear_.residual);
    if (!solve_scaled_newton_system() || !build_cauchy_direction()) {
      update_residual_diagnostics(nonlinear_.residual, base_sum_squares);
      return {};
    }

    const double base_merit = 0.5 * base_sum_squares;
    double accepted_sum_squares = base_sum_squares;
    bool accepted = false;

    auto restore_base_state = [&]() {
      for (int i = 1; i <= kSystemSize; ++i)
        nonlinear_.state[i] = workspace_.base_state[i];
    };

    /*
    Near a solution, an unrestricted Newton step provides quadratic convergence.
    Attempt it first only when both the residual and scaled direction are
    already moderate. The trial still passes the same admissibility and merit
    checks used by the general trust-region path.
    */
    const double base_residual_rms =
        std::sqrt(base_sum_squares / static_cast<double>(kSystemSize));
    const double full_newton_rms = rms_norm(workspace_.newton_direction);
    if (base_residual_rms < 1.0e-7 && std::isfinite(full_newton_rms) &&
        full_newton_rms < 2.0) {
      for (int i = 1; i <= kSystemSize; ++i)
        workspace_.candidate_step[i] = workspace_.newton_direction[i];

      if (apply_retraction(workspace_.base_state, workspace_.candidate_step) &&
          state_is_admissible()) {
        const double trial_sum_squares =
            evaluate_residuals(nonlinear_.trial_residual);
        if (std::isfinite(trial_sum_squares) &&
            residual_is_finite(nonlinear_.trial_residual) &&
            trial_sum_squares < 0.8 * base_sum_squares) {
          accepted_sum_squares = trial_sum_squares;
          status_.trust_radius =
              std::clamp(std::max(status_.trust_radius, 1.25 * full_newton_rms),
                         kMinimumTrustRadius, kMaximumTrustRadius);
          update_residual_diagnostics(nonlinear_.trial_residual,
                                      accepted_sum_squares);
          return {true, full_newton_rms};
        }
      }
      restore_base_state();
    }

    /*
    General dogleg globalization. Rejected trials never modify the accepted
    state: every attempt starts from base_state and failure contracts the
    radius.
    */
    for (int trial = 0; trial < kMaximumTrustTrials; ++trial) {
      restore_base_state();

      if (!build_dogleg_step(status_.trust_radius)) {
        status_.trust_radius =
            std::max(kMinimumTrustRadius, 0.25 * status_.trust_radius);
        continue;
      }

      const double step_norm = rms_norm(workspace_.candidate_step);
      if (!apply_retraction(workspace_.base_state, workspace_.candidate_step) ||
          !state_is_admissible()) {
        status_.trust_radius =
            std::max(kMinimumTrustRadius, 0.25 * status_.trust_radius);
        continue;
      }

      const double trial_sum_squares =
          evaluate_residuals(nonlinear_.trial_residual);
      if (!std::isfinite(trial_sum_squares) ||
          !residual_is_finite(nonlinear_.trial_residual)) {
        status_.trust_radius =
            std::max(kMinimumTrustRadius, 0.25 * status_.trust_radius);
        continue;
      }

      /*
      The model residual is F + J*step in physical residual coordinates.
      Comparing its predicted reduction with the actual nonlinear reduction is
      the standard trust-region quality measure.
      */
      physical_matrix_vector_product(workspace_.candidate_step,
                                     workspace_.model_residual);
      for (int i = 1; i <= kSystemSize; ++i)
        workspace_.model_residual[i] += nonlinear_.residual[i];

      const double model_sum_squares =
          stable_sum_squares(workspace_.model_residual);
      const double actual_reduction = base_merit - 0.5 * trial_sum_squares;
      const double predicted_reduction = base_merit - 0.5 * model_sum_squares;
      const double ratio = predicted_reduction > 0.0
                               ? actual_reduction / predicted_reduction
                               : -std::numeric_limits<double>::infinity();

      if (actual_reduction > 0.0 && ratio > kAcceptanceRatio) {
        accepted = true;
        accepted_sum_squares = trial_sum_squares;

        if (ratio < 0.25) {
          status_.trust_radius =
              std::max(kMinimumTrustRadius, 0.25 * status_.trust_radius);
        } else if (ratio > 0.75 && step_norm >= 0.8 * status_.trust_radius) {
          status_.trust_radius =
              std::min(kMaximumTrustRadius, 2.0 * status_.trust_radius);
        }
        break;
      }

      status_.trust_radius =
          std::max(kMinimumTrustRadius, 0.25 * status_.trust_radius);
    }

    if (!accepted) {
      restore_base_state();
      update_residual_diagnostics(nonlinear_.residual, base_sum_squares);
      return {};
    }

    update_residual_diagnostics(nonlinear_.trial_residual,
                                accepted_sum_squares);
    return {true, rms_norm(workspace_.candidate_step)};
  }

  /*
  Result returned by one fixed-parameter nonlinear correction stage. Keeping
  convergence and iteration count together prevents partially initialized status
  flags from escaping the solver.
  */
  struct NonlinearSolveResult {
    bool converged = false;
    int iterations = 0;
  };

  [[nodiscard]] const char *linear_solver_name() const noexcept {
    switch (status_.linear_solver) {
    case LinearSolverKind::LU:
      return "LU";
    case LinearSolverKind::QR:
      return "QR";
    case LinearSolverKind::Regularized:
      return "REG";
    }
    return "?";
  }

  /*
  Iterate trust-region corrections until both RMS and maximum residual tests
  pass.

  Stagnation detection terminates unproductive solves early so continuation can
  reduce its step and retry from a verified state. A solve never reports success
  based solely on a small state correction.
  */
  [[nodiscard]] NonlinearSolveResult
  solve_nonlinear_system(double tolerance, const char *stage_name) {
    NonlinearSolveResult result;

    /*
    A predictor can occasionally satisfy the target equations immediately.
    Recognize that state before requesting a Newton/Cauchy direction, because
    the zero-gradient least-squares model has no meaningful descent step.
    */
    const double initial_sum_squares = evaluate_residuals(nonlinear_.residual);
    update_residual_diagnostics(nonlinear_.residual, initial_sum_squares);
    if (residual_converged(tolerance)) {
      result.converged = true;
      return result;
    }

    int stagnant_iterations = 0;
    double previous_residual_rms = std::numeric_limits<double>::infinity();

    for (int iter = 1; iter <= kConvergenceMaxIterations; ++iter) {
      result.iterations = iter;
      std::fprintf(stdout, "\n%s iteration %3d:", stage_name, iter);

      const TrustRegionStepResult step = trust_region_iteration(iter);
      std::fprintf(
          stdout,
          " Scaled step RMS: %8.1e  Residual RMS: %8.1e  Residual max: %8.1e"
          "  Trust radius: %8.1e  Linear: %s",
          step.scaled_step_rms, status_.residual_rms, status_.residual_max,
          status_.trust_radius, linear_solver_name());

      if (!step.accepted || !std::isfinite(step.scaled_step_rms) ||
          !std::isfinite(status_.residual_rms) ||
          !std::isfinite(status_.residual_max)) {
        return result;
      }

      if (residual_converged(tolerance)) {
        result.converged = true;
        return result;
      }

      if (status_.residual_rms >= 0.9999 * previous_residual_rms) {
        ++stagnant_iterations;
      } else {
        stagnant_iterations = 0;
      }
      previous_residual_rms = status_.residual_rms;

      if (stagnant_iterations >= 5 ||
          (status_.trust_radius < 1.0e-12 &&
           status_.residual_rms > 100.0 * tolerance)) {
        return result;
      }
    }

    return result;
  }

  /*
  Transform nodal surface elevations to cosine coefficients and copy the solved
  stream coefficients into dedicated post-processing storage.

  Long-double summation is used for the discrete cosine transform. Separating
  post-processing coefficients from nonlinear workspace prevents output routines
  from depending on transient trial states.
  */
  void update_output_coefficients() {
    spectral_.surface_coefficients[0] = 0.0;
    for (int mode = 1; mode <= kFourierOrder; ++mode) {
      spectral_.stream_coefficients[mode] =
          nonlinear_.state[coefficient_state_index(mode)];
      long double surface_sum =
          0.5L *
          (static_cast<long double>(nonlinear_.state[surface_state_index(0)]) +
           static_cast<long double>(
               nonlinear_.state[surface_state_index(kFourierOrder)]) *
               ((mode & 1) ? -1.0L : 1.0L));
      for (int node = 1; node <= kFourierOrder - 1; ++node) {
        surface_sum += static_cast<long double>(
                           nonlinear_.state[surface_state_index(node)]) *
                       static_cast<long double>(cosine_at_node(node, mode));
      }
      spectral_.surface_coefficients[mode] = static_cast<double>(
          2.0L * surface_sum / static_cast<long double>(kFourierOrder));
    }
  }

  /*
  Apply a smooth high-mode filter to continuation predictors only.

  Extrapolation can amplify small unresolved tail coefficients before Newton has
  a chance to restore the equations. The exponential filter damps only the upper
  spectral band of the predicted state; converged states and governing residuals
  are never filtered.
  */
  void filter_predictor_spectrum() {
    if (kFourierOrder < 16)
      return;

    const int cutoff =
        std::max(4, static_cast<int>(0.72 * static_cast<double>(kFourierOrder)));
    const double denominator =
        static_cast<double>(std::max(1, kFourierOrder - cutoff));
    auto &cosine_coefficients = workspace_.predictor_coefficients;
    std::fill(cosine_coefficients.begin(),
              cosine_coefficients.begin() + kFourierOrder + 1, 0.0L);

    for (int mode = 0; mode <= kFourierOrder; ++mode) {
      long double sum_ = 0.5L * static_cast<long double>(
                                    nonlinear_.state[surface_state_index(0)]);
      sum_ += 0.5L *
              static_cast<long double>(
                  nonlinear_.state[surface_state_index(kFourierOrder)]) *
              ((mode & 1) ? -1.0L : 1.0L);
      for (int node = 1; node < kFourierOrder; ++node) {
        sum_ += static_cast<long double>(
                    nonlinear_.state[surface_state_index(node)]) *
                static_cast<long double>(cosine_at_node(node, mode));
      }
      cosine_coefficients[static_cast<std::size_t>(mode)] =
          2.0L * sum_ / static_cast<long double>(kFourierOrder);
    }

    auto filter = [&](int mode) -> long double {
      if (mode <= cutoff)
        return 1.0L;
      const long double xi = static_cast<long double>(mode - cutoff) /
                             static_cast<long double>(denominator);
      const long double xi2 = xi * xi;
      const long double xi4 = xi2 * xi2;
      const long double xi8 = xi4 * xi4;
      return std::exp(-36.0L * xi8 * xi4);
    };

    for (int mode = cutoff + 1; mode <= kFourierOrder; ++mode) {
      const long double sigma = filter(mode);
      cosine_coefficients[static_cast<std::size_t>(mode)] *= sigma;
      nonlinear_.state[coefficient_state_index(mode)] *=
          static_cast<double>(sigma);
    }

    for (int node = 0; node <= kFourierOrder; ++node) {
      long double value = 0.5L * cosine_coefficients[0];
      for (int mode = 1; mode < kFourierOrder; ++mode) {
        value += cosine_coefficients[static_cast<std::size_t>(mode)] *
                 static_cast<long double>(cosine_at_node(node, mode));
      }
      value += 0.5L *
               cosine_coefficients[static_cast<std::size_t>(kFourierOrder)] *
               ((node & 1) ? -1.0L : 1.0L);
      nonlinear_.state[surface_state_index(node)] = static_cast<double>(value);
    }
  }

  /*
  Continuation statistics are returned as one value object. accepted_steps
  counts verified branch advances; trials includes rejected attempts and
  therefore measures the actual globalization work.
  */
  struct ContinuationResult {
    bool converged = false;
    int accepted_steps = 0;
    int trials = 0;
  };

  /*
  Predict from the exact continuation tangent at the most recent converged
  state. Positive variables are extrapolated in logarithmic coordinates; all
  others use fused multiply-add in physical coordinates.
  */
  void predict_from_tangent(const Vector &previous_state, const Vector &tangent,
                            double delta) {
    for (int i = 1; i <= kSystemSize; ++i) {
      nonlinear_.state[i] =
          uses_log_coordinate(i)
              ? log_predict(previous_state[i], tangent[i], delta)
              : std::fma(delta, tangent[i], previous_state[i]);
    }
  }

  /*
  Predict from two converged branch states when a reliable tangent is
  unavailable or produces an inadmissible trial. The extrapolation ratio is
  bounded to avoid excessive branch jumps. Positive variables use a log-secant
  extrapolation.
  */
  void predict_from_secant_history(const Vector &current_state,
                                   const Vector &older_state,
                                   double current_fraction,
                                   double older_fraction,
                                   double trial_fraction) {
    const double fraction_interval = current_fraction - older_fraction;
    const double raw_ratio =
        fraction_interval > 0.0
            ? (trial_fraction - current_fraction) / fraction_interval
            : 0.0;
    const double ratio = std::clamp(raw_ratio, 0.0, 1.5);

    for (int i = 1; i <= kSystemSize; ++i) {
      const double current = current_state[i];
      const double older = older_state[i];
      if (uses_log_coordinate(i) && current > 0.0 && older > 0.0) {
        nonlinear_.state[i] =
            current * std::exp(ratio * (std::log(current) - std::log(older)));
      } else {
        nonlinear_.state[i] = std::fma(ratio, current - older, current);
      }
    }
  }

  /*
  Convert nonlinear iteration count into the next continuation-step multiplier.
  Fast corrections permit controlled growth; expensive corrections reduce the
  next increment before an actual failure occurs.
  */
  [[nodiscard]] static double continuation_growth(int iterations) noexcept {
    if (iterations <= 4)
      return 1.8;
    if (iterations <= 7)
      return 1.4;
    if (iterations <= 12)
      return 1.15;
    if (iterations >= 25)
      return 0.65;
    if (iterations >= 18)
      return 0.8;
    return 1.0;
  }

  /*
  Advance from the linear limit to the requested wave height.

  The controller combines:
    - tangent prediction when a reliable branch derivative is available;
    - secant history prediction as a fallback;
    - exact rollback after every failed nonlinear correction;
    - adaptive step growth from observed Newton iteration counts;
    - conservative maximum growth to avoid branch jumping.

  All history states are converged states. Rejected iterates are never inserted
  into the predictor.
  */
  [[nodiscard]] ContinuationResult continue_to_target(int initial_step_count,
                                                      const char *stage_name,
                                                      bool report_trials) {
    ContinuationResult result;
    const double initial_step =
        1.0 / static_cast<double>(std::max(1, initial_step_count));
    double step = initial_step;
    double fraction = 0.0;
    double previous_previous_fraction = 0.0;

    Vector previous_state(nonlinear_.state.size(), 0.0);
    Vector previous_previous_state(nonlinear_.state.size(), 0.0);
    Vector tangent(nonlinear_.state.size(), 0.0);
    bool have_previous_state = false;
    bool have_previous_previous_state = false;
    bool have_tangent = false;

    while (fraction < 1.0 - 16.0 * std::numeric_limits<double>::epsilon()) {
      if (++result.trials > kMaximumContinuationTrials)
        return result;

      step = std::min(step, 1.0 - fraction);
      const double trial_fraction = fraction + step;
      continuation_.height_parameter =
          trial_fraction * problem_.target_height_parameter;
      continuation_.relative_height =
          trial_fraction * problem_.target_relative_height;

      if (report_trials) {
        std::fprintf(stdout,
                     "\n\nContinuation trial %d: target fraction %.10f, "
                     "increment %.3e\n",
                     result.trials, trial_fraction, step);
      }

      if (!have_previous_state) {
        initialize_state();
      } else if (have_tangent) {
        predict_from_tangent(previous_state, tangent,
                             trial_fraction - fraction);
        if (!state_is_admissible() && have_previous_previous_state) {
          predict_from_secant_history(previous_state, previous_previous_state,
                                      fraction, previous_previous_fraction,
                                      trial_fraction);
        }
      } else if (have_previous_previous_state) {
        predict_from_secant_history(previous_state, previous_previous_state,
                                    fraction, previous_previous_fraction,
                                    trial_fraction);
      } else {
        std::copy(previous_state.begin(), previous_state.end(),
                  nonlinear_.state.begin());
      }

      if (have_previous_state)
        filter_predictor_spectrum();

      const bool final_target =
          trial_fraction >= 1.0 - 16.0 * std::numeric_limits<double>::epsilon();
      const double tolerance =
          final_target ? kFinalResidualTolerance : kConvergenceCriterion;
      const NonlinearSolveResult nonlinear =
          solve_nonlinear_system(tolerance, stage_name);

      if (!nonlinear.converged) {
        if (have_previous_state) {
          std::copy(previous_state.begin(), previous_state.end(),
                    nonlinear_.state.begin());
        }
        step *= 0.5;
        if (step < kMinimumContinuationStep)
          return result;
        if (report_trials) {
          std::fprintf(stderr,
                       "\nContinuation rollback at fraction %.10f; retrying "
                       "with increment %.3e.\n",
                       trial_fraction, step);
        }
        continue;
      }

      if (have_previous_state) {
        previous_previous_state = previous_state;
        previous_previous_fraction = fraction;
        have_previous_previous_state = true;
      }
      previous_state = nonlinear_.state;
      have_previous_state = true;
      fraction = trial_fraction;
      ++result.accepted_steps;
      have_tangent = compute_continuation_tangent(tangent);

      const double growth = continuation_growth(nonlinear.iterations);

      if (fraction < 1.0 - 16.0 * std::numeric_limits<double>::epsilon()) {
        const double remaining = 1.0 - fraction;
        const double maximum_step = std::min(1.5 * initial_step, remaining);
        step = std::min(std::max(step * growth, kMinimumContinuationStep),
                        maximum_step);
      }
    }

    continuation_.height_parameter = problem_.target_height_parameter;
    continuation_.relative_height = problem_.target_relative_height;
    update_output_coefficients();
    result.converged = true;
    return result;
  }

  /*
  ============================================================================
  INPUT, CASE CONFIGURATION, AND TEXT SUMMARIES
  ============================================================================
  */

  // Write the normalized interpretation of the input before numerical work
  // starts.
  void write_input_summary(std::FILE *file) const {
    std::fprintf(file, "# %s", problem_.title.c_str());
    std::fprintf(file, "\n\n# Input data");
    std::fprintf(file, "\n\n# Height/Depth:%6.8f",
                 problem_.target_relative_height);
    if (problem_.wave_specification != WaveSpecification::Period) {
      std::fprintf(file, "\n# Length/Depth:%7.8f",
                   problem_.specified_wavelength);
    } else {
      std::fprintf(file, "\n# Dimensionless Period T*sqrt(g/d):%8.8f",
                   problem_.specified_period);
    }
    std::fprintf(file, "\n# Current criterion: %s,  Dimensionless value:%6.8lf",
                 current_name(), problem_.current_value);
    const std::string method = method_description();
    std::fprintf(file, "\n%s\n", method.c_str());
  }

  /*
  Remove surrounding horizontal whitespace and a possible Windows carriage
  return. Empty or whitespace-only lines become an empty string.
  */
  [[nodiscard]] static std::string trim_line(std::string line) {
    if (!line.empty() && line.back() == '\r')
      line.pop_back();
    const auto first = line.find_first_not_of(" \t");
    if (first == std::string::npos)
      return {};
    const auto last = line.find_last_not_of(" \t");
    return line.substr(first, last - first + 1);
  }

  /*
  Parse the first token of a data-file line. Trailing comments are intentionally
  accepted, allowing human-readable annotations after each value.
  */
  template <typename Value>
  [[nodiscard]] static bool parse_leading_value(const std::string &line,
                                                Value &value) {
    std::istringstream stream(line);
    stream >> value;
    return !stream.fail();
  }

  /*
  ASCII-only case-insensitive comparison for command keywords. The accepted
  vocabulary is deliberately restricted to English protocol tokens, so locale-
  dependent case conversion would add ambiguity without benefit.
  */
  [[nodiscard]] static bool ascii_iequals(std::string_view left,
                                          std::string_view right) noexcept {
    if (left.size() != right.size())
      return false;
    for (std::size_t i = 0; i < left.size(); ++i) {
      const unsigned char a = static_cast<unsigned char>(left[i]);
      const unsigned char b = static_cast<unsigned char>(right[i]);
      const unsigned char lower_a =
          (a >= 'A' && a <= 'Z') ? static_cast<unsigned char>(a + ('a' - 'A'))
                                 : a;
      const unsigned char lower_b =
          (b >= 'A' && b <= 'Z') ? static_cast<unsigned char>(b + ('a' - 'A'))
                                 : b;
      if (lower_a != lower_b)
        return false;
    }
    return true;
  }

  [[nodiscard]] static const char *program_name(int argc,
                                                char **argv) noexcept {
    return (argc > 0 && argv && argv[0]) ? argv[0] : "fourier";
  }

  /*
  Print the complete runtime contract.  The executable intentionally accepts no
  numerical command-line parameters: one data.dat file defines one calculation.
  A help switch is retained because it is useful in installed and scripted
  environments and does not create a second configuration path.
  */
  static void print_usage(std::FILE *file, const char *executable) {
    std::fprintf(
        file,
        "Usage:\n"
        "  %s\n"
        "      Use data.dat when available; otherwise request every required\n"
        "      numerical value interactively. Solve one fixed 100-mode wave\n"
        "      case and write solution.res, surface.res, and flowfield.res.\n\n"
        "  %s --help\n"
        "      Display this message.\n\n"
        "data.dat contains seven non-empty lines:\n"
        "  title\n"
        "  H/d\n"
        "  Period or Wavelength\n"
        "  T*sqrt(g/d) or L/d\n"
        "  current criterion: 1 or 2\n"
        "  dimensionless current magnitude\n"
        "  FINISH\n\n"
        "If data.dat is absent, the same five numerical/keyword fields are\n"
        "requested from the console. No default input values are used.\n",
        executable, executable);
  }

  /*
  Determine the directory containing the executable without assuming that the
  process was launched from that directory.  Native operating-system facilities
  are preferred.  argv[0] is used only as a portable fallback.
  */
  [[nodiscard]] static fs::path
  resolve_executable_directory(const fs::path &fallback_directory, int argc,
                               char **argv) {
    fs::path directory = fallback_directory;
#if defined(_WIN32)
    std::wstring buffer(32768, L'\0');
    const DWORD length = ::GetModuleFileNameW(
        nullptr, buffer.data(), static_cast<DWORD>(buffer.size()));
    if (length > 0 && static_cast<std::size_t>(length) < buffer.size()) {
      buffer.resize(length);
      directory = fs::path(buffer).parent_path();
    }
#elif defined(__linux__)
    char buffer[PATH_MAX] = {0};
    const ssize_t length =
        ::readlink("/proc/self/exe", buffer, sizeof(buffer) - 1);
    if (length > 0) {
      buffer[length] = '\0';
      directory = fs::path(buffer).parent_path();
    }
#endif
    if (directory.empty() && argc > 0 && argv && argv[0]) {
      directory = fs::path(argv[0]).parent_path();
    }
    return directory;
  }

  /*
  Search for data.dat in deterministic priority order:

    1. current working directory;
    2. directory containing the executable.

  Returning an empty path is not a fatal condition.  It explicitly selects the
  interactive input path and prevents accidental use of an unrelated file from
  another search location.
  */
  [[nodiscard]] static fs::path
  locate_data_file(const fs::path &current_directory,
                   const fs::path &executable_directory) {
    const fs::path local = current_directory / "data.dat";
    if (fs::exists(local)) {
      return local;
    }

    const fs::path beside_executable = executable_directory / "data.dat";
    if (fs::exists(beside_executable)) {
      return beside_executable;
    }

    return {};
  }

  /*
  Estimate the first continuation increment from H/d.

  This value is only a seed.  The continuation controller subsequently enlarges
  or reduces the increment according to actual Newton iteration counts and rolls
  back exactly after failed trials.  The quadratic dependence gives gentle
  waves a short path while retaining conservative branch tracking for steep
  waves.
  */
  [[nodiscard]] static int
  choose_initial_continuation_step_count(double relative_height) noexcept {
    const double normalized = std::clamp(relative_height / 0.70, 0.0, 1.5);
    const double estimate = 4.0 + 12.0 * normalized * normalized;
    return std::clamp(static_cast<int>(std::ceil(estimate)), 5, 32);
  }

  /*
  Validate and normalize one finite-depth wave specification.

  Period and wavelength are alternative closures for the same nonlinear system.
  Both are converted to a continuation height parameter so the continuation and
  Newton routines do not contain duplicated case logic.
  */
  [[nodiscard]] bool configure_wave_case(double relative_height,
                                         WaveSpecification wave_specification,
                                         double specified_value,
                                         std::string &error) {
    if (!(relative_height > 0.0) || !std::isfinite(relative_height)) {
      error = "H/d must be a finite positive number";
      return false;
    }
    if (!(specified_value > 0.0) || !std::isfinite(specified_value)) {
      error = "specified period or wavelength must be finite and positive";
      return false;
    }

    ProblemDefinition candidate = problem_;
    candidate.wave_specification = wave_specification;
    if (wave_specification == WaveSpecification::Period) {
      candidate.specified_period = specified_value;
      candidate.specified_wavelength = 0.0;
      candidate.target_height_parameter =
          relative_height / (specified_value * specified_value);
    } else {
      candidate.specified_wavelength = specified_value;
      candidate.specified_period = 0.0;
      candidate.target_height_parameter = relative_height / specified_value;
    }

    candidate.target_relative_height = relative_height;
    if (!(candidate.target_height_parameter > 0.0) ||
        !std::isfinite(candidate.target_height_parameter)) {
      error = "normalized continuation height is not finite and positive";
      return false;
    }

    problem_ = std::move(candidate);
    return true;
  }

  /*
  Read and validate one complete seven-line case from data.dat.

  Blank lines are ignored.  Numeric and keyword lines may contain descriptive
  text after the leading token.  All values are parsed into local temporaries
  and validated before the application state is committed, so a malformed file
  cannot leave a partially configured solver.
  */
  [[nodiscard]] bool read_case_file(const fs::path &path, std::string &error) {
    std::ifstream input(path);
    if (!input) {
      error = "cannot open input file";
      return false;
    }

    std::vector<std::string> lines;
    std::string line;
    while (std::getline(input, line)) {
      line = trim_line(std::move(line));
      if (!line.empty()) {
        lines.push_back(std::move(line));
      }
    }

    if (lines.size() != 7U) {
      error = "data.dat must contain exactly seven non-empty lines";
      return false;
    }
    if (ascii_iequals(lines.front(), "FINISH")) {
      error = "FINISH appears before the wave case";
      return false;
    }
    if (!ascii_iequals(lines.back(), "FINISH")) {
      error = "the final non-empty line must be FINISH";
      return false;
    }

    double relative_height = 0.0;
    double specified_value = 0.0;
    double current_value = 0.0;
    int current_definition = 0;
    std::string specification_token;
    if (!parse_leading_value(lines[1], relative_height) ||
        !parse_leading_value(lines[2], specification_token) ||
        !parse_leading_value(lines[3], specified_value) ||
        !parse_leading_value(lines[4], current_definition) ||
        !parse_leading_value(lines[5], current_value)) {
      error = "one or more numeric fields are invalid";
      return false;
    }
    if (!(current_definition == 1 || current_definition == 2)) {
      error = "current criterion must be 1 or 2";
      return false;
    }
    if (!std::isfinite(current_value)) {
      error = "current magnitude must be finite";
      return false;
    }
    WaveSpecification wave_specification;
    if (ascii_iequals(specification_token, "Period")) {
      wave_specification = WaveSpecification::Period;
    } else if (ascii_iequals(specification_token, "Wavelength")) {
      wave_specification = WaveSpecification::Wavelength;
    } else {
      error = "measure of length must begin with Period or Wavelength";
      return false;
    }
    if (!configure_wave_case(relative_height, wave_specification,
                             specified_value, error)) {
      return false;
    }

    problem_.title = lines[0];
    problem_.current_definition =
        static_cast<CurrentDefinition>(current_definition);
    problem_.current_value = current_value;

    grid_.continuation_steps =
        choose_initial_continuation_step_count(relative_height);
    return true;
  }


  /*
  Read one non-empty line from standard input after displaying and flushing a
  prompt.  Empty entries are not interpreted as defaults: they are returned to
  the caller as invalid input and must be entered again.  Returning false means
  that standard input reached end-of-file or suffered an unrecoverable read
  failure before the requested field was supplied.
  */
  [[nodiscard]] static bool prompt_line(std::string_view prompt,
                                        std::string &line) {
    std::fputs(prompt.data(), stdout);
    std::fflush(stdout);
    if (!std::getline(std::cin, line)) {
      return false;
    }
    line = trim_line(std::move(line));
    return true;
  }

  /*
  Request every required case value directly from the console.

  No numerical defaults are embedded in this path.  Each prompt remains active
  until the user supplies a syntactically and physically valid value.  Parsing
  accepts only the leading token, matching annotated data.dat lines and allowing
  a user to paste either a bare value or a value followed by explanatory text.

  The complete candidate case is committed only after all five fields have been
  entered and configure_wave_case() has validated their coupled normalization.
  This transaction-like behavior prevents a failed or interrupted interactive
  session from leaving a partially configured solver.
  */
  [[nodiscard]] bool read_case_interactively(std::string &error) {
    double relative_height = 0.0;
    double specified_value = 0.0;
    double current_value = 0.0;
    int current_definition = 0;
    WaveSpecification wave_specification = WaveSpecification::Period;
    std::string line;

    for (;;) {
      if (!prompt_line("H/d: ", line)) {
        error = "end-of-input while reading H/d";
        return false;
      }
      if (parse_leading_value(line, relative_height) &&
          relative_height > 0.0 && std::isfinite(relative_height)) {
        break;
      }
      std::fprintf(stderr, "Invalid H/d. Enter a finite positive number.\n");
    }

    for (;;) {
      if (!prompt_line(
              "Measure of length (enter Period or Wavelength): ", line)) {
        error = "end-of-input while reading the measure of length";
        return false;
      }

      std::string token;
      if (parse_leading_value(line, token) && ascii_iequals(token, "Period")) {
        wave_specification = WaveSpecification::Period;
        break;
      }
      if (parse_leading_value(line, token) &&
          ascii_iequals(token, "Wavelength")) {
        wave_specification = WaveSpecification::Wavelength;
        break;
      }
      std::fprintf(stderr,
                   "Invalid measure of length. Enter Period or Wavelength.\n");
    }

    const char *specified_value_prompt =
        wave_specification == WaveSpecification::Period
            ? "Dimensionless period T*sqrt(g/d): "
            : "Dimensionless wavelength L/d: ";
    for (;;) {
      if (!prompt_line(specified_value_prompt, line)) {
        error = "end-of-input while reading the specified period or wavelength";
        return false;
      }
      if (parse_leading_value(line, specified_value) &&
          specified_value > 0.0 && std::isfinite(specified_value)) {
        break;
      }
      std::fprintf(stderr,
                   "Invalid value. Enter a finite positive period or "
                   "wavelength.\n");
    }

    for (;;) {
      if (!prompt_line("Current criterion (1 for Euler, or 2 for Stokes): ", line)) {
        error = "end-of-input while reading the current criterion";
        return false;
      }
      if (parse_leading_value(line, current_definition) &&
          (current_definition == 1 || current_definition == 2)) {
        break;
      }
      std::fprintf(stderr, "Invalid current criterion. Enter 1 or 2.\n");
    }

    for (;;) {
      if (!prompt_line(
              "Current magnitude, dimensionless ubar/sqrt(gd): ", line)) {
        error = "end-of-input while reading the current magnitude";
        return false;
      }
      if (parse_leading_value(line, current_value) &&
          std::isfinite(current_value)) {
        break;
      }
      std::fprintf(stderr,
                   "Invalid current magnitude. Enter a finite number.\n");
    }

    if (!configure_wave_case(relative_height, wave_specification,
                             specified_value, error)) {
      return false;
    }

    problem_.title = "Interactive wave case";
    problem_.current_definition =
        static_cast<CurrentDefinition>(current_definition);
    problem_.current_value = current_value;
    grid_.continuation_steps =
        choose_initial_continuation_step_count(relative_height);
    return true;
  }

  /*
  Write the common converged-wave summary used by the console and result files.
  All values are derived from the final nonlinear state; this routine performs
  no mutation and cannot influence convergence.
  */
  void write_finite_summary(std::FILE *file) const {
    const double wavelength = 2.0 * kPi / nonlinear_.state[StateIndex::depth];
    const double maximum_relative_height =
        (0.0077829 * wavelength * wavelength * wavelength +
         0.0095721 * wavelength * wavelength + 0.141063 * wavelength) /
        (0.0093407 * wavelength * wavelength * wavelength +
         0.0317567 * wavelength * wavelength + 0.078834 * wavelength + 1.0);
    std::fprintf(file, "# %s", problem_.title.c_str());
    {
      const std::string method = method_description();
      std::fprintf(file, "\n%s\n", method.c_str());
    }
    std::fprintf(file,
                 "\n# Height/Depth:%6.5f, %3.0lf%% of the maximum of H/d "
                 "=%6.5f for this length:",
                 nonlinear_.state[StateIndex::wave_height] /
                     nonlinear_.state[StateIndex::depth],
                 nonlinear_.state[StateIndex::wave_height] /
                     nonlinear_.state[StateIndex::depth] /
                     maximum_relative_height * 100.0,
                 maximum_relative_height);
    std::fprintf(file, "\n# Length/Depth:%7.5f",
                 2 * kPi / nonlinear_.state[StateIndex::depth]);
    std::fprintf(file, "\n# Dimensionless Period T*sqrt(g/d):%7.5f",
                 nonlinear_.state[StateIndex::period] /
                     std::sqrt(nonlinear_.state[StateIndex::depth]));
    std::fprintf(file,
                 "\n# Current criterion: %s,  Dimensionless value:%6.5lf\n",
                 current_name(), problem_.current_value);
  }

  /*
  Write one scalar row with both native wavenumber-based and depth-based values.
  The two-column contract is fixed because every supported calculation has a
  finite mean depth.
  */
  void write_result_row(std::FILE *file, const char *description,
                        double native_value, double depth_value) const {
    std::fprintf(file, "\n%s % .16e % .16e", description, native_value,
                 depth_value);
  }

  /*
  ============================================================================
  POST-PROCESSING AND OUTPUT
  ============================================================================
  */

  /*
  Evaluate the smooth cosine interpolant of the converged nodal free surface at
  an arbitrary phase.  Post-processing uses this representation instead of
  nearest-node lookup so surface and flow profiles remain smooth between
  collocation points.
  */
  [[nodiscard]] double surface_elevation(double phase) const {
    long double elevation = 0.0L;
    for (int mode = 1; mode < kFourierOrder; ++mode) {
      elevation +=
          static_cast<long double>(spectral_.surface_coefficients[mode]) *
          std::cos(static_cast<long double>(mode) *
                   static_cast<long double>(phase));
    }
    elevation +=
        0.5L *
        static_cast<long double>(spectral_.surface_coefficients[kFourierOrder]) *
        std::cos(static_cast<long double>(kFourierOrder) *
                 static_cast<long double>(phase));
    return static_cast<double>(elevation);
  }

  /*
  Evaluate boundary conditions between collocation nodes.

  Collocation residuals can be extremely small even when the fixed spectrum is
  under-resolved. Half-shifted verification points therefore provide an
  independent internal defect estimate. The check never changes the prescribed
  100-mode order and is intentionally omitted from result files.
  */
  void update_spectral_defect() {
    const int sample_count = std::max(2 * kFourierOrder, 32);
    long double maximum = 0.0L;

    for (int sample = 0; sample < sample_count; ++sample) {
      const long double x = (static_cast<long double>(sample) + 0.5L) *
                            static_cast<long double>(kPi) /
                            static_cast<long double>(sample_count);
      const double eta = surface_elevation(static_cast<double>(x));
      long double psi_acc = 0.0L;
      long double u_acc = 0.0L;
      long double v_acc = 0.0L;
      bool valid = true;

      for (int mode = 1; mode <= kFourierOrder; ++mode) {
        VerticalBasisValues basis;
        if (!evaluate_vertical_basis(
                mode, eta, nonlinear_.state[StateIndex::depth], basis)) {
          valid = false;
          break;
        }
        const long double phase = static_cast<long double>(mode) * x;
        const long double cosine = std::cos(phase);
        const long double sine = std::sin(phase);
        const long double coefficient =
            static_cast<long double>(spectral_.stream_coefficients[mode]);
        psi_acc += coefficient * static_cast<long double>(basis.s) * cosine;
        u_acc += static_cast<long double>(mode) * coefficient *
                 static_cast<long double>(basis.c) * cosine;
        v_acc += static_cast<long double>(mode) * coefficient *
                 static_cast<long double>(basis.s) * sine;
      }

      if (!valid) {
        status_.spectral_defect_max = std::numeric_limits<double>::infinity();
        return;
      }

      const long double eta_ld = static_cast<long double>(eta);
      const long double kinematic =
          psi_acc -
          static_cast<long double>(nonlinear_.state[StateIndex::flux]) -
          static_cast<long double>(
              nonlinear_.state[StateIndex::wave_frame_speed]) *
              eta_ld;
      const long double relative_u =
          -static_cast<long double>(
              nonlinear_.state[StateIndex::wave_frame_speed]) +
          u_acc;
      const long double dynamic =
          0.5L * (relative_u * relative_u + v_acc * v_acc) + eta_ld -
          static_cast<long double>(nonlinear_.state[StateIndex::bernoulli]);

      maximum = std::max(maximum, std::fabs(kinematic));
      maximum = std::max(maximum, std::fabs(dynamic));
    }

    status_.spectral_defect_max = static_cast<double>(maximum);
  }

  /*
  Complete local kinematic and dynamic state at one phase/elevation pair.

  The structure is returned by value so profile generation cannot observe stale
  mutable intermediates. valid is false whenever basis evaluation or a derived
  quantity becomes non-finite.
  */
  struct FlowPoint {
    bool valid = true;
    double horizontal_velocity = 0.0;
    double vertical_velocity = 0.0;
    double horizontal_gradient = 0.0;
    double vertical_gradient = 0.0;
    double potential_time_derivative = 0.0;
    double horizontal_time_derivative = 0.0;
    double vertical_time_derivative = 0.0;
    double pressure = 0.0;
    double bernoulli_residual = 0.0;
  };

  /*
  Evaluate potential, stream function, velocity gradients, accelerations,
  pressure, and a Bernoulli consistency residual from the converged spectrum.

  Modal sums are accumulated in long double and then converted to depth-based
  nondimensionalization. Time derivatives follow from
  the steady travelling-wave relation.
  */
  [[nodiscard]] FlowPoint
  evaluate_flow_point(double phase, double vertical_coordinate) const {
    FlowPoint point;
    long double horizontal_velocity = 0.0L;
    long double vertical_velocity = 0.0L;
    long double horizontal_gradient = 0.0L;
    long double vertical_gradient = 0.0L;

    /*
    Evaluate the native k-scaled spectral field.  The basis satisfies the bed
    condition analytically, while the trigonometric factors provide horizontal
    periodicity.  Derivative sums are accumulated alongside the field values so
    no mode is evaluated more than once at a point.
    */
    for (int mode = 1; mode <= kFourierOrder; ++mode) {
      VerticalBasisValues basis;
      if (!evaluate_vertical_basis(mode, vertical_coordinate,
                                   nonlinear_.state[StateIndex::depth],
                                   basis)) {
        point.valid = false;
        return point;
      }

      const long double mode_value = static_cast<long double>(mode);
      const long double coefficient =
          static_cast<long double>(spectral_.stream_coefficients[mode]);
      const long double angle = mode_value * static_cast<long double>(phase);
      const long double cosine = std::cos(angle);
      const long double sine = std::sin(angle);

      horizontal_velocity +=
          mode_value * coefficient * static_cast<long double>(basis.c) * cosine;
      vertical_velocity +=
          mode_value * coefficient * static_cast<long double>(basis.s) * sine;
      horizontal_gradient -= mode_value * mode_value * coefficient *
                             static_cast<long double>(basis.c) * sine;
      vertical_gradient += mode_value * mode_value * coefficient *
                           static_cast<long double>(basis.s) * cosine;
    }

    /*
    Convert the native solution to depth-based engineering variables.  The
    physical output coordinate uses y/d=0 at the bed and y/d=1 at mean water
    level, while vertical_coordinate is k(y-d).
    */
    const double depth = nonlinear_.state[StateIndex::depth];
    if (!(depth > 0.0)) {
      point.valid = false;
      return point;
    }

    const double sqrt_depth = std::sqrt(depth);
    const double celerity = nonlinear_.state[StateIndex::celerity] / sqrt_depth;
    const double eulerian_current =
        nonlinear_.state[StateIndex::eulerian_current] / sqrt_depth;
    const double bernoulli_constant =
        1.0 + nonlinear_.state[StateIndex::bernoulli] / depth;
    const double elevation = 1.0 + vertical_coordinate / depth;

    point.horizontal_velocity =
        eulerian_current +
        static_cast<double>(horizontal_velocity) / sqrt_depth;
    point.vertical_velocity =
        static_cast<double>(vertical_velocity) / sqrt_depth;
    point.horizontal_gradient =
        static_cast<double>(horizontal_gradient) * sqrt_depth;
    point.vertical_gradient =
        static_cast<double>(vertical_gradient) * sqrt_depth;

    /*
    A permanent travelling wave satisfies d/dt = -c*d/dx.  This identity gives
    local time derivatives without numerical differencing.  Material
    accelerations then follow from the Eulerian time derivative plus convective
    terms.  Irrotationality and incompressibility provide the cross derivatives.
    */
    point.potential_time_derivative = -celerity * point.horizontal_velocity;
    point.horizontal_time_derivative = -celerity * point.horizontal_gradient;
    point.vertical_time_derivative = -celerity * point.vertical_gradient;
    /*
    Pressure is recovered from Bernoulli's equation in the moving frame.  The
    final residual recomputes the same constant in the stationary frame and is
    used as a local post-processing consistency check.
    */
    point.pressure = bernoulli_constant - elevation -
                     0.5 * ((point.horizontal_velocity - celerity) *
                                (point.horizontal_velocity - celerity) +
                            point.vertical_velocity * point.vertical_velocity);
    point.bernoulli_residual =
        point.potential_time_derivative + point.pressure + elevation +
        0.5 * (point.horizontal_velocity * point.horizontal_velocity +
               point.vertical_velocity * point.vertical_velocity) -
        (bernoulli_constant - 0.5 * celerity * celerity);

    point.valid = std::isfinite(point.pressure) &&
                  std::isfinite(point.bernoulli_residual);
    return point;
  }

  /*
  Compute phase-averaged squared orbital velocity at the bed.

  The steady Eulerian component is removed before squaring, so the result is the
  non-negative oscillatory contribution rather than total current speed. Uniform
  phase sampling and long-double accumulation make the diagnostic independent of
  the output profile grid.
  */
  [[nodiscard]] double
  mean_square_bed_orbital_velocity(int phase_samples = 720) const {
    if (!(nonlinear_.state[StateIndex::depth] > 0.0)) {
      return std::numeric_limits<double>::quiet_NaN();
    }
    const int sample_count = std::max(36, phase_samples);
    const double depth = nonlinear_.state[StateIndex::depth];
    const double sqrt_depth = std::sqrt(depth);
    long double sum_squares = 0.0L;

    for (int sample = 0; sample < sample_count; ++sample) {
      const double phase = 2.0 * kPi * static_cast<double>(sample) /
                           static_cast<double>(sample_count);
      const FlowPoint point = evaluate_flow_point(phase, -depth);
      if (!point.valid)
        return std::numeric_limits<double>::quiet_NaN();
      const double absolute_velocity = point.horizontal_velocity * sqrt_depth;
      const double orbital_velocity =
          absolute_velocity - nonlinear_.state[StateIndex::eulerian_current];
      sum_squares += static_cast<long double>(orbital_velocity) *
                     static_cast<long double>(orbital_velocity);
    }
    return static_cast<double>(sum_squares /
                               static_cast<long double>(sample_count));
  }

  /*
  Scalar quantities derived once from the converged state and reused by all
  output sections. Centralizing these formulas avoids repeated conversions and
  ensures console, tables, and profile headers share identical values.
  */
  struct DerivedQuantities {
    double depth = 1.0;
    double wavelength = 0.0;
    double wave_height = 0.0;
    double period = 0.0;
    double celerity = 0.0;
    double eulerian_current = 0.0;
    double mass_transport_current = 0.0;
    double wave_frame_speed = 0.0;
    double volume_flux = 0.0;
    double bernoulli_constant = 0.0;
    double momentum_flux = 0.0;
    double impulse = 0.0;
    double kinetic_energy = 0.0;
    double potential_energy = 0.0;
    double mean_square_bed_velocity = 0.0;
    double radiation_stress = 0.0;
    double wave_power = 0.0;
  };

  /*
  Convert the native spectral solution into depth-based engineering quantities
  and integral invariants. The routine is pure with respect to solver state.
  */
  [[nodiscard]] DerivedQuantities compute_derived_quantities() const {
    DerivedQuantities result;
    result.depth = nonlinear_.state[StateIndex::depth];
    const double sqrt_depth = std::sqrt(result.depth);
    const double depth_power_3_2 = result.depth * sqrt_depth;
    result.wavelength = 2.0 * kPi / result.depth;
    result.wave_height =
        nonlinear_.state[StateIndex::wave_height] / result.depth;
    result.period = nonlinear_.state[StateIndex::period] / sqrt_depth;
    result.celerity = nonlinear_.state[StateIndex::celerity] / sqrt_depth;
    result.eulerian_current =
        nonlinear_.state[StateIndex::eulerian_current] / sqrt_depth;
    result.mass_transport_current =
        nonlinear_.state[StateIndex::mass_transport_current] / sqrt_depth;
    result.wave_frame_speed =
        nonlinear_.state[StateIndex::wave_frame_speed] / sqrt_depth;
    result.volume_flux = result.wave_frame_speed -
                         nonlinear_.state[StateIndex::flux] / depth_power_3_2;
    result.bernoulli_constant =
        1.0 + nonlinear_.state[StateIndex::bernoulli] / result.depth;
    result.impulse =
        nonlinear_.state[StateIndex::flux] +
        result.depth * nonlinear_.state[StateIndex::eulerian_current];
    result.kinetic_energy =
        0.5 * (nonlinear_.state[StateIndex::celerity] * result.impulse -
               nonlinear_.state[StateIndex::eulerian_current] *
                   result.volume_flux * depth_power_3_2);
    for (int mode = 1; mode <= kFourierOrder; ++mode) {
      result.potential_energy += 0.25 * spectral_.surface_coefficients[mode] *
                                 spectral_.surface_coefficients[mode];
    }
    const double bed_velocity_identity =
        2.0 * nonlinear_.state[StateIndex::bernoulli] -
        nonlinear_.state[StateIndex::celerity] *
            nonlinear_.state[StateIndex::celerity];
    result.radiation_stress =
        4.0 * result.kinetic_energy - 3.0 * result.potential_energy +
        bed_velocity_identity * result.depth +
        2.0 * nonlinear_.state[StateIndex::eulerian_current] *
            (nonlinear_.state[StateIndex::wave_frame_speed] * result.depth -
             nonlinear_.state[StateIndex::flux]);
    result.wave_power =
        nonlinear_.state[StateIndex::celerity] *
            (3.0 * result.kinetic_energy - 2.0 * result.potential_energy) +
        0.5 * bed_velocity_identity *
            (result.impulse +
             nonlinear_.state[StateIndex::celerity] * result.depth) +
        nonlinear_.state[StateIndex::celerity] *
            nonlinear_.state[StateIndex::eulerian_current] *
            (nonlinear_.state[StateIndex::wave_frame_speed] * result.depth -
             nonlinear_.state[StateIndex::flux]);
    result.mean_square_bed_velocity = mean_square_bed_orbital_velocity();
    result.momentum_flux =
        result.radiation_stress -
        2.0 * nonlinear_.state[StateIndex::celerity] * result.impulse +
        (nonlinear_.state[StateIndex::celerity] *
             nonlinear_.state[StateIndex::celerity] +
         0.5 * result.depth) *
            result.depth;
    return result;
  }

  /*
  Write the scalar invariants and converged spectral coefficients. Keeping this
  file-specific formatter separate from surface and field sampling makes each
  output contract independently reviewable and prevents unrelated formatting
  changes from becoming entangled.
  */
  void write_solution_file(std::FILE *file,
                           const DerivedQuantities &quantities) const {
    write_finite_summary(file);

    /*
    The Stokes-Ursell number provides a compact nonlinearity/dispersion measure
    for the converged solution.  It is reported only as a diagnostic; it does
    not enter the nonlinear solve.
    */
    std::fprintf(file, "\n# Stokes-Ursell number %7.3f",
                 0.5 * nonlinear_.state[StateIndex::wave_height] /
                     std::pow(nonlinear_.state[StateIndex::depth], 3));

    std::fprintf(file, "\n\n# Integral quantities");
    std::fprintf(file,
                 "\n# (1) Quantity, (2) symbol, solution non-dimensionalised "
                 "by (3) g & wavenumber, and (4) g & mean depth\n");

    write_result_row(file, "# Water depth                        (d)",
                     nonlinear_.state[StateIndex::depth], 1.0);
    write_result_row(file, "# Wave length                   (lambda)",
                     2.0 * kPi, quantities.wavelength);
    write_result_row(file, "# Wave height                        (H)",
                     nonlinear_.state[StateIndex::wave_height],
                     quantities.wave_height);
    write_result_row(file, "# Wave period                      (tau)",
                     nonlinear_.state[StateIndex::period], quantities.period);
    write_result_row(file, "# Wave speed                         (c)",
                     nonlinear_.state[StateIndex::celerity],
                     quantities.celerity);
    write_result_row(file, "# Eulerian current                  (u1)",
                     nonlinear_.state[StateIndex::eulerian_current],
                     quantities.eulerian_current);
    write_result_row(file, "# Stokes current                    (u2)",
                     nonlinear_.state[StateIndex::mass_transport_current],
                     quantities.mass_transport_current);
    write_result_row(file, "# Mean fluid speed in frame of wave  (U)",
                     nonlinear_.state[StateIndex::wave_frame_speed],
                     quantities.wave_frame_speed);
    write_result_row(file, "# Volume flux due to waves           (q)",
                     nonlinear_.state[StateIndex::flux],
                     nonlinear_.state[StateIndex::flux] /
                         std::pow(quantities.depth, 1.5));
    write_result_row(file, "# Bernoulli constant                 (r)",
                     nonlinear_.state[StateIndex::bernoulli],
                     nonlinear_.state[StateIndex::bernoulli] /
                         quantities.depth);

    write_result_row(file, "# Volume flux                        (Q)",
                     quantities.volume_flux * std::pow(quantities.depth, 1.5),
                     quantities.volume_flux);
    write_result_row(file, "# Bernoulli constant                 (R)",
                     quantities.bernoulli_constant * quantities.depth,
                     quantities.bernoulli_constant);
    write_result_row(file, "# Momentum flux                      (S)",
                     quantities.momentum_flux,
                     quantities.momentum_flux /
                         (quantities.depth * quantities.depth));
    write_result_row(file, "# Impulse                            (I)",
                     quantities.impulse,
                     quantities.impulse / std::pow(quantities.depth, 1.5));
    write_result_row(file, "# Kinetic energy                     (T)",
                     quantities.kinetic_energy,
                     quantities.kinetic_energy /
                         (quantities.depth * quantities.depth));
    write_result_row(file, "# Potential energy                   (V)",
                     quantities.potential_energy,
                     quantities.potential_energy /
                         (quantities.depth * quantities.depth));
    write_result_row(file, "# Mean square of bed velocity      (ub2)",
                     quantities.mean_square_bed_velocity,
                     quantities.mean_square_bed_velocity / quantities.depth);
    write_result_row(file, "# Radiation stress                 (Sxx)",
                     quantities.radiation_stress,
                     quantities.radiation_stress /
                         (quantities.depth * quantities.depth));
    write_result_row(file, "# Wave power                         (F)",
                     quantities.wave_power,
                     quantities.wave_power / std::pow(quantities.depth, 2.5));

    /*
    Both coefficient sets are emitted with all 100 modes.  B_j describes the
    stream-function representation and E_j the cosine representation of the
    free surface used for continuous post-processing.
    */
    std::fprintf(file, "\n\n# Dimensionless coefficients in Fourier series");
    std::fprintf(file, "\n# Potential/Streamfn\tSurface elevations");
    std::fprintf(file, "\n# j, B[j], & E[j], j=1..N\n");
    for (int mode = 1; mode <= kFourierOrder; ++mode) {
      std::fprintf(file, "\n%4d  % .16e  % .16e", mode,
                   spectral_.stream_coefficients[mode],
                   spectral_.surface_coefficients[mode]);
    }

    std::fprintf(file,
                 "\n\n# Numerical verification"
                 "\n# Final Fourier order: %d"
                 "\n# Collocation residual RMS: %.16e"
                 "\n# Collocation residual max: %.16e\n",
                 kFourierOrder, status_.residual_rms, status_.residual_max);
  }

  /*
  Sample the free surface on a crest-clustered quadratic phase grid. Every row
  is evaluated directly from the final spectral coefficients; no interpolation
  of previously sampled output is used.
  */
  void write_surface_file(std::FILE *file,
                          const DerivedQuantities &quantities) const {
    std::fprintf(file, "# %s\n", problem_.title.c_str());
    const std::string method = method_description();
    std::fprintf(file, "%s\n", method.c_str());
    std::fprintf(file, "\n# Surface of wave - trough-crest-trough,");
    std::fprintf(file, " note quadratic point spacing clustered around crest");
    std::fprintf(file, "\n# Non-dimensionalised with respect to depth");
    std::fprintf(file, "\n# X/d, eta/d, & check of surface pressure\n");
    std::fprintf(file, "\n0.\t0.\t0. # Dummy point to scale plot\n");

    for (int sample = -grid_.surface_point_count / 2;
         sample <= grid_.surface_point_count / 2; ++sample) {
      const double sample_value = static_cast<double>(sample);
      const double point_count = static_cast<double>(grid_.surface_point_count);
      const double phase = 4.0 * kPi * sample_value * std::fabs(sample_value) /
                           (point_count * point_count);
      const double eta = surface_elevation(phase);
      const FlowPoint point = evaluate_flow_point(phase, eta);
      std::fprintf(file, "\n%8.4lf\t%7.4f\t%7.0e", phase / quantities.depth,
                   1.0 + eta / quantities.depth, point.pressure);
    }
    std::fprintf(file, "\n\n");
  }

  /*
  Write vertical kinematic profiles at uniformly spaced half-wave phases. Local
  velocity, acceleration, pressure, and Bernoulli diagnostics are evaluated from
  the same immutable converged state used by the solution summary.
  */
  void write_flowfield_file(std::FILE *file,
                            const DerivedQuantities &quantities) const {
    std::fprintf(file, "# %s\n", problem_.title.c_str());
    const std::string method = method_description();
    std::fprintf(file, "%s\n", method.c_str());
    std::fprintf(
        file, "\n# Velocity and acceleration profiles and Bernoulli checks\n");
    std::fprintf(
        file,
        "\n# All quantities are dimensionless with respect to g and/or d\n");
    std::fprintf(file, "\n#****************************************************"
                       "***************************");
    std::fprintf(file,
                 "\n# y        u       v    dphi/dt   du/dt   dv/dt  du/dx   "
                 "du/dy Bernoulli check  ");
    std::fprintf(file, "\n# -     -------------   -------  ------   -----  "
                       "------------- ---------------  ");
    std::fprintf(file, "\n# d        sqrt(gd)       gd        g       g       "
                       "sqrt(g/d)        gd         ");
    std::fprintf(file, "\n#****************************************************"
                       "***************************");
    std::fprintf(file, "\n# Increasing X/d and Phase below moves from crest to "
                       "trough over half a wave.");
    std::fprintf(
        file,
        "\n# At a fixed physical x, X=x-c*t decreases as time increases.\n");

    for (int profile = 0; profile <= grid_.profile_count; ++profile) {
      const double phase = kPi * static_cast<double>(profile) /
                           static_cast<double>(grid_.profile_count);
      const double eta = surface_elevation(phase);
      std::fprintf(file, "\n\n# X/d = %8.7f, Phase = %6.1f\n",
                   phase / quantities.depth, phase * 180.0 / kPi);

      for (int point_index = 0; point_index <= grid_.vertical_point_count;
           ++point_index) {
        const double vertical_output =
            static_cast<double>(point_index) * (1.0 + eta / quantities.depth) /
            static_cast<double>(grid_.vertical_point_count);
        const FlowPoint point = evaluate_flow_point(
            phase, quantities.depth * (vertical_output - 1.0));
        std::fprintf(file,
                     "\n% 10.6f % 10.6f % 10.6f % 10.6f % 10.6f % 10.6f % "
                     "10.6f % 10.6f % 10.6f",
                     vertical_output, point.horizontal_velocity,
                     point.vertical_velocity, point.potential_time_derivative,
                     point.horizontal_time_derivative,
                     point.vertical_time_derivative, point.horizontal_gradient,
                     point.vertical_gradient, point.bernoulli_residual);
      }
    }
    std::fprintf(file, "\n\n");
  }

  /*
  Generate solution.res, surface.res, and flowfield.res from the verified state.

  Numerical solving is complete before this routine begins. Output generation is
  therefore read-only and cannot alter convergence state or spectral resolution.
  */
  void write_outputs(const OutputHandles &output) const {
    const DerivedQuantities quantities = compute_derived_quantities();

    std::fprintf(stdout, "\n\n# Solution summary:\n\n");
    write_finite_summary(stdout);
    write_solution_file(output.solution, quantities);
    write_surface_file(output.surface, quantities);
    write_flowfield_file(output.flowfield, quantities);
  }

  /*
  ============================================================================
  MEMORY MANAGEMENT AND APPLICATION CONTROL
  ============================================================================
  */

  /*
  Allocate every persistent numerical array exactly once for the fixed 100-mode
  system.  No dynamic allocation occurs inside continuation, Newton iteration,
  residual evaluation, Jacobian assembly, or output sampling.

  One extra element is retained because the mathematical state and residual
  numbering begin at one.  Dense matrices use the corresponding one-based row
  and column convention through DenseMatrix.
  */
  void allocate_workspace() {

    const std::size_t vector_size =
        static_cast<std::size_t>(kSystemSize + 1);
    const std::size_t spectral_size = static_cast<std::size_t>(kFourierOrder + 1);
    const std::size_t table_side = spectral_size;

    spectral_.surface_coefficients.assign(spectral_size, 0.0);
    spectral_.stream_coefficients.assign(spectral_size, 0.0);
    spectral_.cosine_table.reserve(table_side * table_side);
    spectral_.sine_table.reserve(table_side * table_side);

    nonlinear_.state.assign(vector_size, 0.0);
    nonlinear_.residual.assign(vector_size, 0.0);
    nonlinear_.trial_residual.assign(vector_size, 0.0);

    workspace_.jacobian.resize(kSystemSize);
    workspace_.linear_matrix.resize(kSystemSize);
    workspace_.linear_rhs.assign(vector_size, 0.0);
    workspace_.newton_direction.assign(vector_size, 0.0);
    workspace_.cauchy_direction.assign(vector_size, 0.0);
    workspace_.candidate_step.assign(vector_size, 0.0);
    workspace_.base_state.assign(vector_size, 0.0);
    workspace_.row_scale.assign(vector_size, 0.0);
    workspace_.column_scale.assign(vector_size, 0.0);
    workspace_.model_residual.assign(vector_size, 0.0);
    workspace_.pivot.clear();
    workspace_.pivot.reserve(vector_size);
    workspace_.predictor_coefficients.assign(spectral_size, 0.0L);
  }

  /*
  Execute adaptive continuation to the exact target H/d.  The fixed spectral
  order never changes.  After convergence, the off-grid defect is evaluated as
  an internal resolution check; an elevated defect produces a warning but does
  not alter the user-specified production discretization.
  */
  [[nodiscard]] bool solve_wave() {
    const ContinuationResult continuation = continue_to_target(
        grid_.continuation_steps, "spectral continuation", true);
    if (!continuation.converged) {
      std::fprintf(stderr,
                   "\nFatal: adaptive continuation failed before reaching "
                   "the target wave height. No unconverged result files were "
                   "written.\n");
      return false;
    }

    std::fprintf(stdout,
                 "\nAdaptive continuation completed in %d accepted steps and "
                 "%d total trials.\n",
                 continuation.accepted_steps, continuation.trials);

    update_spectral_defect();
    if (!std::isfinite(status_.spectral_defect_max) ||
        status_.spectral_defect_max > kSpectralDefectTolerance) {
      std::fprintf(stderr,
                   "Warning: the fixed 100-mode solution has an elevated "
                   "internal between-node boundary defect.\n");
    }
    return true;
  }

public:
  /*
  Application entry point.

  Responsibilities are deliberately limited to:
    - locating and parsing input;
    - allocating solver storage;
    - invoking the numerical pipeline;
    - opening output files with RAII ownership;
    - dispatching read-only result generation;
    - returning a meaningful process status.
  */
  int run(int argc, char **argv) {
    /*
    The only accepted option is --help. Numerical command-line arguments are
    intentionally rejected: cases come from data.dat when it exists and from
    explicit interactive prompts otherwise.
    */
    if (argc == 2 && argv && argv[1] &&
        (std::string_view(argv[1]) == "-h" ||
         std::string_view(argv[1]) == "--help")) {
      print_usage(stdout, program_name(argc, argv));
      return 0;
    }
    if (argc != 1) {
      std::fprintf(
          stderr,
          "This executable accepts no numerical command-line arguments.\n\n");
      print_usage(stderr, program_name(argc, argv));
      return 2;
    }

    /*
    Locate and parse input before allocating the large dense work matrices.  The
    output directory is the directory containing the selected data.dat, which
    keeps all files for one calculation together.
    */
    const fs::path current_directory = fs::current_path();
    const fs::path executable_directory =
        resolve_executable_directory(current_directory, argc, argv);
    const fs::path data_path =
        locate_data_file(current_directory, executable_directory);
    fs::path io_directory = current_directory;
    std::string input_error;

    if (!data_path.empty()) {
      io_directory = data_path.parent_path();
      std::fprintf(stdout, "Using input file: %s\n",
                   fs::absolute(data_path).string().c_str());
      if (!read_case_file(data_path, input_error)) {
        std::fprintf(stderr, "Invalid data.dat: %s.\n", input_error.c_str());
        return 2;
      }
    } else {
      std::fprintf(stdout,
                   "data.dat was not found in the working directory or beside "
                   "the executable.\n"
                   "Enter every required input value manually; no defaults "
                   "will be used.\n\n");
      if (!read_case_interactively(input_error)) {
        std::fprintf(stderr, "Interactive input failed: %s.\n",
                     input_error.c_str());
        return 2;
      }
    }

    write_input_summary(stdout);
    std::fprintf(stdout,
                 "Fixed Fourier order: %d\n"
                 "Adaptive continuation seed: %d increments\n",
                 kFourierOrder, grid_.continuation_steps);

    /*
    Allocate, solve, and verify before creating any result file.  This ordering
    guarantees that a failed nonlinear calculation cannot leave apparently
    valid but unconverged output.
    */
    allocate_workspace();
    if (!solve_wave()) {
      return 3;
    }

    const fs::path solution_path = io_directory / "solution.res";
    const fs::path surface_path = io_directory / "surface.res";
    const fs::path flowfield_path = io_directory / "flowfield.res";

    FilePtr solution_file{std::fopen(solution_path.string().c_str(), "w")};
    FilePtr surface_file{std::fopen(surface_path.string().c_str(), "w")};
    FilePtr flowfield_file{std::fopen(flowfield_path.string().c_str(), "w")};
    if (!solution_file || !surface_file || !flowfield_file) {
      std::perror("output file");
      return 4;
    }

    const OutputHandles output{solution_file.get(), surface_file.get(),
                               flowfield_file.get()};
    write_outputs(output);

    /*
    Explicitly flush and close each stream so delayed write failures are visible
    to the process exit status.  RAII still protects every earlier return path.
    */
    const bool solution_closed =
        close_output_file(solution_file, "solution.res");
    const bool surface_closed = close_output_file(surface_file, "surface.res");
    const bool flowfield_closed =
        close_output_file(flowfield_file, "flowfield.res");
    if (!(solution_closed && surface_closed && flowfield_closed)) {
      return 4;
    }

    std::fprintf(stdout,
                 "\nFinished. Output files:\n"
                 "  %s\n"
                 "  %s\n"
                 "  %s\n",
                 fs::absolute(solution_path).string().c_str(),
                 fs::absolute(surface_path).string().c_str(),
                 fs::absolute(flowfield_path).string().c_str());
    return 0;
  }
};

} // namespace wave

/*
Process entry point.  All resources and numerical state are owned by the
application object; main only converts uncaught exceptions into a nonzero exit
status and a concise diagnostic suitable for scripts and batch execution.
*/
int main(int argc, char **argv) {
  try {
    return wave::SteadyWaveApplication{}.run(argc, argv);
  } catch (const std::exception &error) {
    std::fprintf(stderr, "Fatal error: %s\n", error.what());
    return 1;
  } catch (...) {
    std::fprintf(stderr, "Fatal error: unknown exception.\n");
    return 1;
  }
}
