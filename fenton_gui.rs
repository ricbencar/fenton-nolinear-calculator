// ============================================================================
// FENTON NONLINEAR WAVE SOLVER - COMPLETE ENGINEERING AND BUILD REFERENCE
// ============================================================================
// PROGRAM      : Fenton Nonlinear Wave Solver, native Windows graphical program
// SOURCE       : fenton_gui.rs (single source file, Rust 2021 edition)
// FORMULATION  : Rienecker-Fenton steady nonlinear stream-function/Fourier
//                collocation method for periodic gravity waves in finite depth.
// IMPLEMENTED  : Nonlinear wavelength, wave/current interaction, hydrodynamics,
//                kinematics, integral quantities, accuracy checks and Win32 GUI.
// BUILD HOST   : Windows 10 x64 (1809+) or Windows 11 x64, MSYS2 UCRT64 shell.
// TARGET       : x86_64-pc-windows-gnu (GNU/MinGW-w64, not Microsoft MSVC).
// DEPENDENCIES : Rust standard library, GNU linker, native Windows system DLLs.
// EXTERNAL RUST CRATES: None. Cargo and Visual Studio are not required.
//
// IMPORTANT PLATFORM DISTINCTION
// ------------------------------
// The instructions in this header are for MODERN Windows using the actively
// supported MSYS2 UCRT64 environment. Current MSYS2 installers/packages do not
// support Windows 7 as a build host. An older Windows 7 toolchain is a DIFFERENT
// configuration and is not covered by the commands below. As of 2026, MSYS2
// recommends UCRT64 rather than the legacy/deprecated MINGW64 environment.
// The Windows runtime DLLs imported by the resulting executable should always
// be inspected before distributing it to a different machine.
//
// ============================================================================
//  1. INSTALLATION AND COMPILATION - MSYS2 UCRT64 + PACMAN + RUSTC
// ============================================================================
//
// STEP 1 - INSTALL MSYS2
// ---------------------
// Download the current x86_64 installer from https://www.msys2.org/ .
// Install MSYS2 to C:\msys64 (the recommended short, ASCII-only path).
// Open the "MSYS2 UCRT64" shortcut from the Windows Start menu. Its prompt
// identifies UCRT64 and its toolchain executables reside in /ucrt64/bin.
// DO NOT use the plain "MSYS2 MSYS" shell for this native Windows GNU build.
// DO NOT mix /mingw64/bin and /ucrt64/bin tools or libraries in one build.
//
// STEP 2 - UPDATE MSYS2 (in the MSYS2 UCRT64 BASH shell)
// ------------------------------------------------------
//     pacman -Syu
//
// If pacman requests termination/restart of the terminal, close it, reopen
// "MSYS2 UCRT64", and run 'pacman -Syu' again to finish updating.
// Finish updates BEFORE installing compiler packages. An update can change the
// base runtime and pacman database, so avoid a partially updated MSYS2 system.
//
// STEP 3 - INSTALL GNU COMPILER, LINKER, RUST AND CARGO
// ----------------------------------------------------
// Still inside the *UCRT64* Bash shell:
//
//     pacman -S --needed base-devel mingw-w64-ucrt-x86_64-toolchain mingw-w64-ucrt-x86_64-rust
//
// Accept the package selection and installation confirmations. This command
// installs GCC, binutils, MinGW-w64 headers/runtime, rustc and Cargo. Both
// rustc and GCC must come from the SAME UCRT64 prefix. A systemwide rustup
// installation is not needed and should not shadow the MSYS2 pacman compiler.
//
// STEP 4 - CHECK THE ACTIVE ENVIRONMENT (MSYS2 UCRT64 BASH)
// --------------------------------------------------------
//     echo "$MSYSTEM"                       # expected: UCRT64
//     command -v rustc                     # expected: /ucrt64/bin/rustc
//     command -v cargo                     # expected: /ucrt64/bin/cargo
//     command -v gcc                       # expected: /ucrt64/bin/gcc
//     rustc --version
//     cargo --version
//     gcc --version
//     rustc -vV                            # host: x86_64-pc-windows-gnu
//     rustc --print target-libdir --target x86_64-pc-windows-gnu
//
// If rustc comes from ~/.cargo/bin or a different prefix, fix PATH before
// compiling. Do not combine a separately installed MSVCRT-based Rust stdlib
// with a UCRT64 linker, or mix conflicting versions of Rust stdlib/linker.
//
// STEP 5 - ENTER THE SOURCE DIRECTORY (MSYS2 UCRT64 BASH)
// -------------------------------------------------------
// Example for a project stored in the same Windows path as this application:
//
//     cd /c/fenton-nolinear-calculator-main/rust
//     ls -l fenton_gui.rs
//
// In MSYS2, the Windows path C:\folder\file becomes /c/folder/file. Quotes
// are required if any directory name contains spaces. All source code is in
// fenton_gui.rs: no Cargo.toml and no external Rust dependencies are required.
//
// STEP 6 - COMPILE THE NATIVE 64-BIT WINDOWS GUI (CPU-OPTIMIZED)
// ----------------------------------------------
// For the fastest candidate on THIS build computer, select target-cpu=native.
// It enables instruction selection and scheduling for the local x86-64 CPU.
// Explicitly disable hardware FMA in the application code to avoid
// changing the separate multiplication/addition rounding assumptions used
// by the outward-rounded interval arithmetic. Keep IEEE binary64 (f64),
// -C panic=unwind and every existing numerical acceptance threshold. These
// optimizations do NOT imply bit-for-bit identical math-library results on
// all CPUs or that native will be faster in every case; run the provided
// benchmark/regression checks. This EXE is not portable to older CPUs.
//
// Paste the following COMPLETE command into MSYS2 UCRT64 Bash on ONE LINE:
//
//     rustc --edition=2021 --target x86_64-pc-windows-gnu -C linker=gcc -C opt-level=3 -C codegen-units=1 -C lto=fat -C target-cpu=native -C target-feature=-fma -C prefer-dynamic=no -C panic=unwind -C strip=symbols "fenton_gui.rs" -o "fenton_gui.exe"
//
// Equivalently, use the following Bash continuation syntax (BACKSLASH, not a
// PowerShell backtick; no trailing whitespace after each backslash):
//
//     rustc --edition=2021 --target x86_64-pc-windows-gnu \
//         -C linker=gcc \
//         -C opt-level=3 \
//         -C codegen-units=1 \
//         -C lto=fat \
//         -C target-cpu=native \
//         -C target-feature=-fma \
//         -C prefer-dynamic=no \
//         -C panic=unwind \
//         -C strip=symbols \
//         "fenton_gui.rs" -o "fenton_gui.exe"
//
// PORTABLE ALTERNATIVE (baseline x86-64, no CPU-specific ISA requirements):
//     rustc --edition=2021 --target x86_64-pc-windows-gnu -C linker=gcc -C opt-level=3 -C codegen-units=1 -C lto=fat -C target-cpu=x86-64 -C prefer-dynamic=no -C panic=unwind -C strip=symbols "fenton_gui.rs" -o "fenton_gui_portable.exe"
//
// For convenience, if the companion build_msys2_ucrt64.sh is available:
//
//     bash build_msys2_ucrt64.sh gui         # local-CPU optimized GUI
//     bash build_msys2_ucrt64.sh portable    # baseline x86-64 GUI
//     bash build_msys2_ucrt64.sh all         # GUI, CLI and unit tests
//     bash build_msys2_ucrt64.sh bench       # comparative numeric/timing run
//
// 'linker=gcc' explicitly selects the installed UCRT64 GNU linker driver.
// The old Windows 7 command forced '-C link-self-contained=yes'. That option
// is deliberately OMITTED here: pacman supplies a complete UCRT64 GNU linker
// and native libraries. rustc can select the appropriate native link objects
// using the active environment. Avoid blindly mixing the rustc-supplied
// self-contained MSVCRT objects with UCRT64 compiler/linker packages.
//
// STEP 7 - RUN THE PROGRAM AND CHECK NATIVE DLL IMPORTS
// -----------------------------------------------------
//     test -f fenton_gui.exe && echo 'GUI build OK'
//     ./fenton_gui.exe
//     objdump -p fenton_gui.exe | grep -i 'DLL Name'
//
// GUI execution uses the Windows subsystem, without an extra console window.
// Enter wave height H [m], period T [s], depth d [m], current Uc [m/s], and
// press CALCULATE. The report appears in the interface and is also written to
// output.txt in the CURRENT WORKING DIRECTORY. The disk report is UTF-8;
// the Win32 EDIT control receives UTF-16 text with Windows CRLF line breaks.
//
// A native Windows program normally imports Windows-owned DLLs such as
// KERNEL32.dll, USER32.dll, GDI32.dll, ucrtbase.dll, etc. These are NOT bundled.
// 'prefer-dynamic=no' prefers static Rust libraries but is not a promise that
// every GNU support library is statically linked. Inspect the import list for
// non-system DLLs (e.g. libwinpthread-1.dll or libgcc_s_seh-1.dll). If present,
// package compatible DLLs or investigate static native linking and RETEST.
// Do not assume that a successful compilation proves stand-alone deployment.
//
// OPTIONAL - ATTEMPT STATIC LINKAGE OF GNU SUPPORT LIBRARIES
// ----------------------------------------------------------
// The following native-linker flags may reduce dependencies on MinGW DLLs,
// when corresponding static libraries are available in your UCRT64 toolchain:
//
//     rustc --edition=2021 --target x86_64-pc-windows-gnu -C linker=gcc -C opt-level=3 -C codegen-units=1 -C lto=fat -C target-cpu=x86-64 -C prefer-dynamic=no -C panic=unwind -C strip=symbols -C link-arg=-static -C link-arg=-static-libgcc "fenton_gui.rs" -o "fenton_gui.exe"
//
// This optional command is not guaranteed to work with every GCC/Rust package
// combination. It cannot remove essential Windows operating-system imports.
// After using it, inspect 'objdump -p' and run the application again.
//
// OPTIONAL - TEST EXECUTABLE (COMPILED FROM THE SAME FILE)
// -------------------------------------------------------
//     rustc --edition=2021 --test --target x86_64-pc-windows-gnu -C linker=gcc -C panic=unwind "fenton_gui.rs" -o "fenton_tests.exe"
//     ./fenton_tests.exe
//     ./fenton_tests.exe --ignored           # slower numerical regression tests
//
// The test build is a console-mode program. The normal GUI build uses the
// 'windows' subsystem attribute; no WinMain wrapper or third-party GUI crate
// is necessary. These checks do not substitute for accuracy validation.
//
// OPTIONAL - COMPILE FROM AN EXISTING WINDOWS POWERSHELL SESSION
// --------------------------------------------------------------
// The recommended build terminal is MSYS2 UCRT64. To run its Bash environment
// from PowerShell, while preserving the current working directory:
//
//     $env:MSYSTEM = 'UCRT64'
//     $env:CHERE_INVOKING = 'yes'
//     & 'C:\msys64\usr\bin\bash.exe' -lc 'bash build_msys2_ucrt64.sh'
//
// Alternatively, start the MSYS2 UCRT64 terminal and use the direct 'rustc'
// command above. A companion build_msys2_ucrt64.bat can launch the same build
// from PowerShell or cmd.exe when MSYS2 lives at C:\msys64.
//
// MEANING OF THE RUSTC OPTIONS
// ----------------------------
// --edition=2021             Rust 2021 syntax, as used by this source file.
// --target x86_64-pc-windows-gnu
//                            Windows x64 binary using the GNU ABI.
// -C linker=gcc              Use the currently selected UCRT64 GCC driver.
// -C opt-level=3             Optimized release machine code.
// -C codegen-units=1         Single LLVM code-generation unit, potentially
//                            improving optimization at the cost of compile time.
// -C lto=fat                 Whole-program link-time optimization.
// -C target-cpu=native        Generate code tuned to the exact build computer.
//                            Do NOT redistribute to an incompatible processor.
// -C target-feature=-fma
//                            Avoid hardware fused multiply-add so interval
//                            rounding remains based on separate operations.
// -C target-cpu=x86-64        Portable baseline alternative (see above).
// -C prefer-dynamic=no       Prefer static linkage of Rust dependencies.
// -C panic=unwind            Retain catch_unwind used to handle numerical
//                            domain/interval evaluation errors gracefully.
// -C strip=symbols           Remove debug symbols from the release executable.
// -C link-arg=-static        OPTIONAL GCC argument for static native support.
// -o fenton_gui.exe          Name the resulting graphical Windows executable.
//
// For numerical reproducibility, do not enable fast-math or relaxed floating-
// point reassociation. On CPUs with FMA, -fma keeps multiplication and
// addition separate. This avoids changing the interval-rounding model.
// IEEE-754 binary64 interval tests depend on arithmetic
// assumptions; verify reproducibility across compiler versions and CPUs.
//
// DIAGNOSTICS AND TROUBLESHOOTING
// -------------------------------
//  - 'rustc: command not found': verify UCRT64, run 'command -v rustc', and
//    install mingw-w64-ucrt-x86_64-rust through pacman.
//  - 'linker gcc not found': install mingw-w64-ucrt-x86_64-toolchain and check
//    'command -v gcc' reports /ucrt64/bin/gcc.
//  - 'can't find crate for std': inspect rustc -vV, rustc --print target-libdir,
//    and ensure the GNU standard library matches the installed Rust toolchain.
//  - '.cargo/bin/rustc' is first on PATH: remove the conflicting PATH entry
//    from the UCRT64 session or deliberately select a single toolchain family.
//  - 'undefined reference' / runtime mismatch: check linker resolution and
//    prevent mixing MSVCRT-based MINGW64 with UCRT64 libraries.
//  - 'pacman -Syu' requests a restart: close and reopen UCRT64 before
//    running package installation; repeat full update if necessary.
//  - 'MSVCP...dll' or 'vcruntime...dll' dependency: check whether the binary
//    was accidentally built with MSVC or depends on third-party libraries.
//  - 'libwinpthread-1.dll missing': inspect native imports, then include a
//    compatible redistributable or rebuild with appropriate static linking.
//  - Preserved certificate metadata has narrowly scoped #[allow(dead_code)]
//    attributes: no warnings are suppressed for the numerical solver itself.
//  - Compiler 'error:' or a nonzero status means compilation failed.
//  - Compile with panic=unwind, not panic=abort: the solver uses catch_unwind.
//  - Use target-cpu=x86-64 for portability; target-cpu=native is a
//    machine-specific optimization, not a universal speed guarantee.
//  - Current UCRT64 packages target modern Windows; do not assume that the
//    resulting binary or MSYS2 itself can run on Windows 7 or 8.1.
//
// ALTERNATIVE ENVIRONMENT - LEGACY MINGW64 (NOT THE RECOMMENDED CHOICE)
// --------------------------------------------------------------------
// Existing MINGW64 installations historically use:
//
//     pacman -S --needed base-devel mingw-w64-x86_64-toolchain mingw-w64-x86_64-rust
//
// MINGW64 uses /mingw64/bin and MSVCRT rather than /ucrt64/bin and UCRT.
// Its packages are being phased out. Use either MINGW64 consistently OR UCRT64
// consistently, never a mixed runtime/linker environment. For new Windows
// 10/11 installations, follow the UCRT64 instructions above.
//
// Official references: https://www.msys2.org/docs/environments/
//                      https://www.msys2.org/docs/windows_support/
//                      https://packages.msys2.org/packages/mingw-w64-ucrt-x86_64-rust
//                      https://doc.rust-lang.org/rustc/codegen-options/
//
// ============================================================================
//  2. PHYSICAL SCOPE, COORDINATES AND INPUT CONVENTIONS
// ============================================================================
//
// This program treats a two-dimensional, inviscid, irrotational, steady
// progressive regular gravity wave propagating over a horizontal, impermeable
// bottom. It is not an irregular-sea, breaking-wave, turbulent or 3-D model.
// The primary engineering relationship computed is:
//
//     L = F(H, T, d, Uc),            k = 2*pi/L,          c = L/T
//
// with H [m] the crest-to-trough height of ONE regular wave (not Hs),
// T [s] the fixed-frame wave period, d [m] the still-water depth, and
// Uc [m/s] the prescribed collinear Eulerian current. Positive Uc follows
// the direction of propagation; negative Uc opposes the wave.
//
// Introduce the travelling-wave coordinate X and vertical coordinate Y:
//
//     X = k (x - c t),       Y = k y,       y = 0 at still-water level,
//     y = -d at the bed,    Y_bed = -kd,  zeta(X) = k eta(X).
//
// The numerical surface is symmetric: X_m = m*pi/N for m=0,...,N. One
// complete wavelength corresponds to a phase increase of 2*pi.
// The solver uses g = 9.80665 m/s^2 and rho = 1025.0 kg/m^3.
// The nondimensional stream function is Psi = psi sqrt(k^3/g), with the
// corresponding velocity scale sqrt(g/k). The Fenton state is 1-indexed:
//
//   z[1]    = kd                         nondimensional water depth
//   z[2]    = kH                         nondimensional wave height
//   z[3]    = T sqrt(gk)                 nondimensional wave period
//   z[4]    = c sqrt(k/g)                nondimensional phase velocity
//   z[5]    = ubar_1 sqrt(k/g)           Eulerian current variable
//   z[6]    = ubar_2 sqrt(k/g)           Stokes/mass-transport current variable
//   z[7]    = Ubar sqrt(k/g)             mean wave-frame velocity variable
//   z[8]    = q sqrt(k^3/g)             volume-flux-related variable
//   z[9]    = r k/g                      reduced Bernoulli constant variable
//   z[10+m] = zeta_m = k eta(X_m),       m=0,...,N
//   z[N+10+j] = B_j,                    j=1,...,N
//
// There are 2*N+10 active unknowns and as many nonlinear equations. z[0] is
// deliberately unused to keep the implementation close to Fenton's notation.
// The index ranges MUST NOT be shifted when modifying equations or Jacobians.
//
// ============================================================================
//  3. FOURIER STREAM FUNCTION AND FINITE-DEPTH BASIS
// ============================================================================
//
// The steady stream-function representation is
//
//     Psi(X,Y) = -Ubar* (Y+kd)
//              + sum_{j=1..N} B_j S_j(Y) cos(jX),
//     Ubar*    = Ubar sqrt(k/g),
//
//     S_j(Y) = sinh[j(Y+kd)] / cosh(jkd)
//            = [exp(jY)-exp(-jY-2jkd)]/[1+exp(-2jkd)],
//
//     C_j(Y) = cosh[j(Y+kd)] / cosh(jkd)
//            = [exp(jY)+exp(-jY-2jkd)]/[1+exp(-2jkd)].
//
// At Y=-kd the S_j terms vanish and the impermeable-bed condition is
// satisfied analytically. The exponential ratios avoid overflow and loss of
// significance that can occur in naive evaluation of large hyperbolic terms.
// The free surface zeta(X) is reconstructed from cosine coefficients after
// each accepted continuation step; its highest harmonic |E_N| is a useful
// (but not sufficient) diagnostic of spectral resolution.
//
// ============================================================================
//  4. FULL NONLINEAR RESIDUAL SYSTEM
// ============================================================================
//
// Let h_s = H/(g*T^2) and h_r = H/d. For the active continuation fraction,
// the code evaluates stepped versions of h_s and h_r. Global residuals are:
//
//    F1 = z2 - z1*h_r
//    F2 = z2 - h_s*z3^2
//    F3 = z4*z3 - 2*pi
//    F4 = z5 + z7 - z4
//    F5 = z1*(z6+z7-z4) - z8
//    F6 = z[current_criterion+4] - (Uc/sqrt(gd))*sqrt(z1)
//    F7 = z10 + z[N+10] + 2*sum_{m=1..N-1} z[10+m]
//    F8 = z10 - z[N+10] - z2
//
// F7 imposes the half-wave mean-level constraint and F8 the wave height.
// For every collocation node m = 0,...,N calculate Fourier sums psi_m,
// u_m, v_m from S_j, C_j and the harmonic trigonometric table, then impose:
//
//    F[9+m]    = psi_m - z8 - z7*z[10+m]
//    F[N+10+m] = 0.5*((u_m-z7)^2+v_m^2) + z[10+m] - z9.
//
// These are the nonlinear surface streamline and Bernoulli conditions,
// coupled to the dispersion/period/current/flux conditions above. The
// explicit analytic Jacobian J=dF/dz is assembled with identical row/column
// indexing. This is a finite-amplitude nonlinear solution, not an empirical
// correction applied to the linear Airy dispersion relation.
//
// ============================================================================
//  5. CURRENTS, MASS FLUX, BERNOULLI AND DIMENSIONAL SCALES
// ============================================================================
//
// Current criterion 1 specifies Eulerian current (GUI default); criterion 2
// permits the alternate Stokes-current closure in the solver object. The
// Fenton variables are distinct:
//
//    ubar_1 = c - Ubar,            ubar_2 = c - Q/d,
//    q      = Ubar*d - Q,          r = R - g*d.
//
// The code maintains z[5] and z[6] rather than conflating Eulerian and Stokes
// current. The report includes velocity, volume flux, Bernoulli constants,
// energy, stress and wave power, using the following dimensional scales:
//
//    velocity:  sqrt(gd)          energy/area: rho*g*d^2
//    flux:      sqrt(g*d^3)       force/length: rho*g*d^2
//    impulse:   rho*sqrt(g*d^3)  power/length: rho*g^(3/2)*d^(5/2)
//
// Numerical energy, radiation-stress, impulse and power formulas are derived
// from the accepted Fourier state, and several are checked against independent
// phase/depth integrations. Reported engineering units are m, s, m/s, m/s^2,
// kN/m, kJ/m^2 and kW/m as indicated in the report tables.
//
// ============================================================================
//  6. NUMERICAL SOLUTION METHOD
// ============================================================================
//
// (a) Initialize kd using a Fenton–McKee-like approximation to seed the
//     nonlinear solve; no linear approximation is used as the final answer.
// (b) Increase H in adaptive continuation increments, initially four steps.
//     A failed step is retried with a smaller height increment.
// (c) Solve J*delta_z=-F at each Newton step with a column-equilibrated,
//     one-sided Jacobi singular-value decomposition (SVD). The relative
//     singular-value cutoff is 1e-15. Equilibration reduces column scaling
//     disparities; SVD is more robust than naive Gaussian elimination in
//     ill-conditioned wave problems.
// (d) Backtrack alpha=1,1/2,1/4,... until the trial decreases sum(F_i^2),
//     has positive kd and contains only finite state variables.
// (e) Require max|F_i| <= 1e-8 at intermediate continuation stages and
//     <= 1e-10 at final wave height, plus a sufficiently small surface
//     correction. Do not lower these limits merely to force convergence.
// (f) Estimate Jacobian conditioning by SVD. Above 1e11, reduce the Fourier
//     order to a minimum of 12 and repeat, retaining report diagnostics.
// (g) Convert the accepted nondimensional state to dimensional outputs and
//     re-evaluate all equations to ensure normalization has not changed the
//     solution beyond its residual tolerance.
//
// All primary scalar arithmetic uses IEEE-754 binary64 (Rust f64). Bitwise
// results can depend on the compiler, target math library, CPU and flags;
// neither LTO nor -O3 should be treated as a substitute for regression tests.
//
// ============================================================================
//  7. VELOCITY, ACCELERATION AND ENGINEERING DIAGNOSTICS
// ============================================================================
//
// At arbitrary (X,Y), the Fourier sums yield horizontal and vertical velocity
// in the fixed frame. Velocity derivatives give the material horizontal
// acceleration a_x = partial_t(u) + u*partial_x(u) + w*partial_y(u).
// The travelling-wave relation gives partial_t(u)=-c*partial_x(u).
// Kinematics are converted from dimensionless to [m/s] and [m/s^2].
//
// The implementation searches over phases 0..pi, refines local candidates
// with a golden-section method and reports surface/bed peak absolute velocity
// and horizontal material acceleration. The mean-square ORBITAL bed velocity
// is evaluated via orthogonality of the Fourier harmonics, separately from
// the Eulerian-current contribution to mean-square absolute bed velocity.
//
// The Miche criterion, Fenton (1990) highest-steady-wave fit, 0.98*Hm
// indicator, Ursell number, surface-tail coefficient and crest stagnation
// margin are diagnostics. They are not stability or physical wave-breaking
// proofs and should not be interpreted as measured wave-breaking thresholds.
//
// ============================================================================
//  8. VALIDATION, LIMITS OF INTERVAL CHECKING AND NUMERICAL TESTS
// ============================================================================
//
// The acceptance pipeline includes:
//    (i)  a final nonlinear collocation residual <= 1e-10;
//   (ii)  interpolation of the free-surface equations at >=2048 phases,
//         with a maximum surface residual <= 1e-6;
//  (iii)  spectral-vs-integral energy/momentum/power identities within
//         1e-7 + 1e-8*abs(reference);
//   (iv)  two independent phase/depth Gauss–Legendre quadratures of different
//         resolution, using 64/96 vertical nodes;
//    (v)  binary64 interval enclosures for all finite-dimensional residuals
//         and nine selected output formulas.
//
// The interval module expands endpoints by one adjacent binary64 number and
// evaluates conservative bounds with explicit Taylor remainders for exp and
// trigonometric collocation. These checks bound a finite algebraic evaluation
// under the assumptions of the implemented arithmetic. They are NOT a proof
// that a nearby exact solution of the water-wave boundary-value problem exists;
// neither full physical stability nor all machine/external-library behavior
// is certified. A failed accuracy check rejects the affected wave output.
//
// Unit tests cover input validation, a challenging SVD matrix, initialization
// and fixed-width report rendering. Ignored, more expensive reference tests
// exercise nonlinear solutions and selected wavelength/power values.
//
// ============================================================================
//  9. PROGRAM STRUCTURE AND WINDOWS EXECUTION
// ============================================================================
//
//   linalg      Dense equilibrated one-sided Jacobi SVD Newton solve.
//   verified    Binary64 interval operations and finite-system certificates.
//   solver      Full nonlinear wave, continuation, kinematics and checks.
//   report      Fixed-width engineering tables and glossary.
//   windows_gui Native Win32 GUI (only when building for Windows GUI mode).
//   main        GUI or command-line entry point under compile-time cfg flags.
//
// Win32 interoperation is declared directly using extern "system" imports from
// user32.dll, gdi32.dll and kernel32.dll, with no third-party GUI framework.
// The UI sends heavy wave computations to a worker thread, preventing the
// message loop from blocking. On completion the worker posts a message to the
// UI thread, which installs the report into a read-only multiline EDIT control.
// A thread panicking during report generation is handled by catch_unwind;
// writing output.txt can still fail for filesystem/permission reasons.
//
// ============================================================================
// 10. REFERENCES AND ENGINEERING USE
// ============================================================================
//
// Rienecker, M. M. & Fenton, J. D. (1981), A Fourier approximation method for
// steady water waves, Journal of Fluid Mechanics, 104, 119–137.
// Fenton, J. D. (1990), Nonlinear wave theories, The Sea, Vol. 9, Part A.
// Fenton, J. D. & McKee, W. D. (1990), On calculating the lengths of water waves.
//
// Reported H denotes regular-wave height, not significant wave height Hs.
// Verify applicability of the steady periodic inviscid model and compare with
// independent solutions when designing critical coastal or marine structures.
// ============================================================================
//
#![cfg_attr(all(windows, not(any(feature = "headless", test))), windows_subsystem = "windows")]

// ===========================================================================
// MODULE: linalg
// ===========================================================================
// ============================================================================
// MODULE linalg — SVD STABILIZATION FOR THE NONLINEAR NEWTON STEP
// ============================================================================
// The Jacobian is a dense square matrix of size 2*N+10. The solver uses a
// column-equilibrated one-sided Jacobi SVD rather than normal equations,
// avoiding unnecessary squaring of the condition number. The rotation sequence
// is explicitly ordered: changing it can alter floating-point convergence.
mod linalg {
//! Column-equilibrated one-sided Jacobi SVD for the nonlinear Newton system.
//!
//! Given the row-major Jacobian A and a right-hand side b, compute a truncated
//! pseudo-inverse solution to A*x=b. Column equilibration acts before the Jacobi
//! sweeps; column magnitudes and rotations are accumulated in column-major
//! workspaces W and V so that each column is contiguous in memory. The function
//! additionally returns an estimated condition ratio sigma_max/sigma_min.

/// Solve A*x = b by column-equilibrated one-sided Jacobi SVD.
///
/// `matrix` is n-by-n and stored in **row-major** order; the returned x is
/// n-dimensional. Internally W and V store **columns contiguously**.
///
/// Steps: validate shapes -> scale columns -> orthogonalize column pairs ->
/// compute singular values -> apply a relative singular-value cutoff -> form
/// the least-squares correction -> undo scaling. Condition estimates are
/// informative but do not establish a rigorous bound on the Jacobian inverse.
/// A non-converged Jacobi sweep or a non-finite datum returns `Err`.
pub fn svd_solve(matrix: &[f64], rhs: &[f64], n: usize) -> Result<(Vec<f64>, f64), String> {
    if matrix.len() != n * n || rhs.len() != n {
        return Err("Invalid matrix dimensions".to_string());
    }
    // W = A * diag(1/scale): each column is equilibrated independently.
    // V accumulates the right singular-vector rotations, starting as identity.
    // The vectors are initialized once and reused throughout Jacobi sweeps.
    let mut w = vec![0.0; n * n];
    let mut v = vec![0.0; n * n];
    let mut scale = vec![1.0; n];
    for j in 0..n {
        let mut maximum = 0.0_f64;
        for i in 0..n {
            maximum = maximum.max(matrix[i * n + j].abs());
        }
        if maximum > 0.0 { scale[j] = maximum; }
    }
    for j in 0..n {
        v[j * n + j] = 1.0;
        for i in 0..n {
            let value = matrix[i * n + j];
            if !value.is_finite() { return Err("Non-finite Jacobian".into()); }
            w[j * n + i] = value / scale[j];
        }
        if !rhs[j].is_finite() { return Err("Non-finite linear right-hand side".into()); }
    }
    // Small dot products relative to the column norms indicate numerical
    // orthogonality. Use a tolerance tied to binary64 machine epsilon, not a
    // user-specified Newton residual. At most 100 full sweeps are attempted.
    let tolerance = 32.0 * f64::EPSILON;
    let mut orthogonal = false;
    for _ in 0..100 {
        orthogonal = true;
        for p in 0..n {
            for q in (p + 1)..n {
                // Contiguous column slices permit bounds-check elimination and
                // efficient machine-code vectorization without changing the
                // summation or rotation order of the Jacobi algorithm.
                // Rust split_at_mut gives nonoverlapping writable slices for
                // the two matrix columns. This avoids copying a column and
                // preserves the exact order of floating-point operations.
                let (w_left, w_right) = w.split_at_mut(q * n);
                let wp = &mut w_left[p * n..(p + 1) * n];
                let wq = &mut w_right[..n];
                let (mut pp, mut qq, mut pq) = (0.0, 0.0, 0.0);
                for i in 0..n {
                    let a = wp[i];
                    let b = wq[i];
                    pp += a * a;
                    qq += b * b;
                    pq += a * b;
                }
                if pp <= 0.0 || qq <= 0.0 || pq.abs() <= tolerance * pp.sqrt() * qq.sqrt() {
                    continue;
                }
                orthogonal = false;
                // Stable tangent formula for the Jacobi plane rotation:
                // tau = (||q||^2-||p||^2)/(2 p.q), t = sign(tau)/
                // (|tau|+sqrt(1+tau^2)), c=1/sqrt(1+t^2), s=t*c.
                // Apply the same orthogonal rotation to W and V.
                let tau = (qq - pp) / (2.0 * pq);
                let t = 1.0_f64.copysign(tau) / (tau.abs() + 1.0_f64.hypot(tau));
                let c = 1.0 / 1.0_f64.hypot(t);
                let sn = t * c;
                let (v_left, v_right) = v.split_at_mut(q * n);
                let vp = &mut v_left[p * n..(p + 1) * n];
                let vq = &mut v_right[..n];
                for i in 0..n {
                    let a = wp[i];
                    let b = wq[i];
                    wp[i] = c * a - sn * b;
                    wq[i] = sn * a + c * b;
                    let va = vp[i];
                    let vb = vq[i];
                    vp[i] = c * va - sn * vb;
                    vq[i] = sn * va + c * vb;
                }
            }
        }
        if orthogonal { break; }
    }
    // An incomplete SVD is not a valid Newton linear correction. Stop rather
    // than silently using a partially orthogonalized Jacobian factorization.
    if !orthogonal { return Err("Jacobi SVD did not converge".into()); }
    // Each W column norm is an approximate singular value. hypot accumulates
    // Euclidean norms without the avoidable overflow/underflow of naive sum
    // of squares. The condition ratio is a post-rotation diagnostic.
    let mut sigma: Vec<f64> = vec![0.0; n];
    let mut x = vec![0.0; n];
    let mut largest = 0.0_f64;
    for j in 0..n {
        for i in 0..n { sigma[j] = sigma[j].hypot(w[j * n + i]); }
        largest = largest.max(sigma[j]);
    }
    if !largest.is_finite() { return Err("Non-finite singular value".into()); }
    let smallest = sigma.iter().copied().fold(f64::INFINITY, f64::min);
    let condition = if smallest > largest * 1e-300 { largest / smallest } else { f64::INFINITY };
    // Reconstruct x = D^{-1} V Sigma^+ U^T b. Singular values less than
    // 1e-15 of the maximum are omitted to avoid unstable Newton corrections.
    // The omitted modes constitute a truncated pseudo-inverse, not a proof
    // that the original Jacobian is full rank.
    for j in 0..n {
        let sj = sigma[j];
        if sj <= largest * 1e-15 { continue; }
        let mut projection = 0.0;
        for i in 0..n { projection += (w[j * n + i] / sj) * rhs[i]; }
        for i in 0..n { x[i] += v[j * n + i] * (projection / sj); }
    }
    // Undo the original column scaling to recover the correction in the
    // physical ordering of the 1-based Fenton state (converted to 0-based dx).
    for j in 0..n { x[j] /= scale[j]; }
    Ok((x, condition))
}
}

// ===========================================================================
// MODULE: verified
// ===========================================================================
// ============================================================================
// MODULE verified — ENCLOSURES FOR THE FINITE-DIMENSIONAL PROBLEM
// ============================================================================
// The purpose is to bound floating-point evaluations of residuals and nine
// derived formulas after the nonlinear solver has converged. This module is
// separate from the SVD and Newton solver: it does not find a new root.
//
// All interval operations expect finite endpoints. Invalid domains, division
// by an interval containing zero, overflow or a failed square-root bound
// cause a panic. The caller contains these panics using catch_unwind.
mod verified {
//! Binary64 interval enclosures for the finite residuals and nine derived outputs.
//! This is a *finite-system check*, not a computer-assisted proof of root existence.
//!
//! Endpoint expansion uses adjacent IEEE-754 binary64 numbers. Exponentials
//! are enclosed by a reduced-argument Taylor sum plus explicit remainder,
//! and collocation trigonometric values by finite Taylor bounds. These checks
//! have scope limited to the implemented finite algebraic expressions.
use std::collections::BTreeMap;
use std::ops::{Add, Div, Mul, Neg, Sub};

/// Adjacent IEEE-754 binary64 number in the direction of increasing values.
/// Move to the adjacent representable `f64` toward positive infinity.
/// For ±0, the next strictly positive subnormal number is returned; NaN and
/// positive infinity are propagated. This is one-ULP outward expansion.
fn next_up(x: f64) -> f64 {
    if x.is_nan() || x == f64::INFINITY { return x; }
    if x == 0.0 { return f64::from_bits(1); }
    let b = x.to_bits();
    f64::from_bits(if x > 0.0 { b + 1 } else { b - 1 })
}
/// Move to the adjacent representable `f64` toward negative infinity.
/// The sign-aware integer-bit adjustment works on finite positive and negative
/// IEEE-754 numbers; both signed zeros map to the smallest negative subnormal.
fn next_down(x: f64) -> f64 {
    if x.is_nan() || x == f64::NEG_INFINITY { return x; }
    if x == 0.0 { return -f64::from_bits(1); }
    let b = x.to_bits();
    f64::from_bits(if x > 0.0 { b - 1 } else { b + 1 })
}

/// Closed binary64 interval [lo, hi] with finite, ordered endpoints.
/// `Copy` permits small values to be passed by value without aliasing.
#[derive(Debug, Copy, Clone, Default)]
pub struct Interval { pub lo: f64, pub hi: f64 }
impl Interval {
    // Zero-radius interval describing an exact representable input value.
    // It is not the mathematical real-number uncertainty of a decimal literal.
    pub fn point(x: f64) -> Self { Self::new(x, x) }
    // All constructors enforce finite, ascending endpoints. NaN/Inf arithmetic
    // is deliberately treated as an invalid certificate rather than accepted.
    pub fn new(lo: f64, hi: f64) -> Self {
        assert!(lo.is_finite() && hi.is_finite() && lo <= hi, "Invalid interval endpoints");
        Self { lo, hi }
    }
    pub fn down(x: f64) -> f64 { next_down(x) }
    // Largest absolute magnitude over the closed interval; useful for a
    // conservative infinity norm of an equation residual.
    pub fn magnitude(self) -> f64 { self.lo.abs().max(self.hi.abs()) }
    // Square is monotone only on nonnegative or nonpositive input intervals.
    // If zero lies inside, the mathematical lower bound is exactly zero.
    // Enlarge floating-point lower/upper products by adjacent f64 numbers.
    pub fn square(self) -> Self {
        let lower = if self.lo <= 0.0 && self.hi >= 0.0 { 0.0 } else { self.lo.abs().min(self.hi.abs()) };
        let upper = self.magnitude();
        Self::new(0.0_f64.max(next_down(lower * lower)), next_up(upper * upper))
    }
    // Integer exponentiation by repeated squaring, O(log n) multiplications.
    // Used to compute Taylor terms and the high-order remainder estimates.
    pub fn power(self, mut n: u32) -> Self {
        let mut value = Self::point(1.0);
        let mut base = self;
        while n != 0 {
            if n & 1 != 0 { value = value * base; }
            n >>= 1;
            if n != 0 { base = base.square(); }
        }
        value
    }
    // Enclose exp(x) for one binary64 argument. Reduce |x| to <=0.125 by
    // repeated halving, sum exp(r)=sum(r^j/j!) through j=32, then bound the
    // positive Taylor remainder by the next power with an exp(r)<2 factor.
    // Use reciprocity for x<0 and square back once per reduction step.
    // Domain and overflow errors propagate rather than reporting false checks.
    fn exponential_point(x: f64) -> Self {
        if x == 0.0 { return Self::point(1.0); }
        let mut reduced = x.abs();
        let mut halvings = 0;
        while reduced > 0.125 { reduced *= 0.5; halvings += 1; }
        // The reduced argument is itself enclosed one ULP in both directions.
        // The lower bound is clamped at zero because reduction uses abs(x).
        let r = Self::new(0.0_f64.max(next_down(reduced)), next_up(reduced));
        let mut sum = Self::point(1.0);
        let mut term = Self::point(1.0);
        for j in 1..=32 {
            term = term * r / Self::point(j as f64);
            sum = sum + term;
        }
        // First omitted Taylor term is r^33/33!. For 0 <= r <= 0.125,
        // multiplying by 2 is a conservative bound on the remaining tail.
        let mut remainder = r.power(33) * Self::point(2.0);
        for j in 2..=33 { remainder = remainder / Self::point(j as f64); }
        sum = Self::new(sum.lo, next_up(sum.hi + remainder.hi));
        if x < 0.0 { sum = Self::point(1.0) / sum; }
        for _ in 0..halvings { sum = sum.square(); }
        sum
    }
    // exp is increasing, so bounding the two endpoints suffices to enclose
    // the result on the whole input interval.
    pub fn exp(self) -> Self {
        Self::new(Self::exponential_point(self.lo).lo, Self::exponential_point(self.hi).hi)
    }
    // sqrt is defined only for nonnegative intervals. Starting from the
    // hardware sqrt, adjust endpoints outward and verify by squaring. The
    // finite adjustment loop fails explicitly if the bracket cannot be shown.
    pub fn sqrt(self) -> Self {
        assert!(self.lo >= 0.0, "Negative square-root interval");
        if self.hi == 0.0 { return Self::point(0.0); }
        let mut lower = if self.lo == 0.0 { 0.0 } else { 0.0_f64.max(next_down(self.lo.sqrt())) };
        let mut upper = next_up(self.hi.sqrt());
        for _ in 0..16 {
            if lower <= 0.0 || next_up(lower * lower) <= self.lo { break; }
            lower = 0.0_f64.max(next_down(lower));
        }
        for _ in 0..16 {
            if next_down(upper * upper) >= self.hi { break; }
            upper = next_up(upper);
        }
        assert!((lower == 0.0 || next_up(lower * lower) <= self.lo) && next_down(upper * upper) >= self.hi,
            "Square-root enclosure failed");
        Self::new(lower, upper)
    }
}
// Elementary interval operators. Addition/subtraction enlarge the exact
// endpoint sums/differences. Multiplication encloses all four endpoint
// products. Division requires the denominator not to contain zero and
// multiplies by an outward-expanded reciprocal interval. Unary negation
// reverses endpoints without additional rounding.
impl Add for Interval { type Output=Self; fn add(self,b:Self)->Self { Self::new(next_down(self.lo+b.lo),next_up(self.hi+b.hi)) } }
impl Neg for Interval { type Output=Self; fn neg(self)->Self { Self::new(-self.hi,-self.lo) } }
impl Sub for Interval { type Output=Self; fn sub(self,b:Self)->Self { self + (-b) } }
impl Mul for Interval {
    type Output=Self;
        // Four cross-products are necessary because the extrema may occur
        // at any corner of the input rectangle [a.lo,a.hi] x [b.lo,b.hi].
    fn mul(self,b:Self)->Self {
        let p=[self.lo*b.lo,self.lo*b.hi,self.hi*b.lo,self.hi*b.hi];
        assert!(p.iter().all(|v|v.is_finite()),"Interval product overflow");
        let low=p.iter().copied().fold(f64::INFINITY,f64::min);
        let high=p.iter().copied().fold(f64::NEG_INFINITY,f64::max);
        Self::new(next_down(low),next_up(high))
    }
}
impl Div for Interval {
    type Output=Self;
        // A denominator spanning zero makes division unbounded or undefined:
        // such a certificate must fail instead of returning a finite number.
    fn div(self,b:Self)->Self {
        assert!(b.lo > 0.0 || b.hi < 0.0,"Interval divisor contains zero");
        self * Self::new(next_down(1.0/b.hi),next_up(1.0/b.lo))
    }
}
/// Interval enclosing mathematical pi around Rust's f64 PI constant.
/// Subsequent collocation trigonometric bounds are based on this enclosure.
pub fn pi() -> Interval { Interval::new(next_down(std::f64::consts::PI),next_up(std::f64::consts::PI)) }

/// Sin/cos enclosures at m*pi/N; Taylor enclosures match the source formulation.
/// Compute conservative (sin, cos) intervals at the 2N phase indices.
///
/// Symmetry reduces phase indices after N to the negative interval; angles
/// are obtained from m*pi/N. Taylor series through 32 sine/cosine terms are
/// supplemented by the first neglected-term bounds, with interval arithmetic
/// on every operation. Lookup indices are used in the residual certificate.
pub fn collocation(n: usize) -> Vec<(Interval,Interval)> {
    let i = Interval::point;
    let mut table=Vec::with_capacity(2*n);
    for index in 0..2*n {
        // Fourier phases modulo 2*pi: indices N+1..2N-1 correspond to the
        // negative phase partners, keeping |angle| <= pi for Taylor accuracy.
        let q=if index<=n { index as i32 } else {index as i32 - 2*n as i32};
        let x=pi()*i(q as f64)/i(n as f64);
        let square=x.square();
        let(mut sine,mut cosine,mut st,mut ct)=(i(0.0),i(0.0),x,i(1.0));
        for m in 0..32 {
            sine=sine+st; cosine=cosine+ct;
            st = -st*square/i(((2*m+2)*(2*m+3)) as f64);
            ct = -ct*square/i(((2*m+1)*(2*m+2)) as f64);
        }
        // The first omitted sine term has magnitude at most |x|^65/65!,
        // and the first omitted cosine term at most |x|^64/64! as evaluated
        // by this series layout. Symmetric error intervals expand both sums.
        let mut sin_error=i(x.magnitude()).power(65);
        let mut cos_error=i(x.magnitude()).power(64);
        for j in 2..=65 { sin_error=sin_error/i(j as f64); }
        for j in 2..=64 { cos_error=cos_error/i(j as f64); }
        table.push((sine+Interval::new(-sin_error.hi,sin_error.hi),
                    cosine+Interval::new(-cos_error.hi,cos_error.hi)));
    }
    table
}

#[derive(Clone,Debug)]
/// One checked engineering output: an interval result, a reported scalar,
/// its conservative discrepancy and its permitted error budget.
pub struct OutputBound {
    // Retained for forensic inspection of the full outward-rounded output bound.
    #[allow(dead_code)]
    pub enclosure: Interval,
    // Original engineering scalar retained alongside the conservative bound.
    #[allow(dead_code)]
    pub reported: f64,
    pub error_upper: f64,
    pub budget_lower: f64,
}
#[derive(Clone,Debug)]
/// Complete finite-system checking record for one nonlinear wave solution.
/// `root_existence_proved` is deliberately false: this is a residual and
/// output-formula enclosure, not a fixed-point or interval-Newton existence
/// theorem. The report exposes only selected diagnostic information.
pub struct Certificate {
    // Human-readable provenance of the IEEE binary64 enclosure method.
    #[allow(dead_code)]
    pub arithmetic: String,
    // Explicit scope of the certified finite equations and output formulas.
    #[allow(dead_code)]
    pub scope: String,
    // Deliberately false: finite residual checks do not prove a nearby root.
    #[allow(dead_code)]
    pub root_existence_proved: bool,
    pub passed: bool,
    pub order: usize,
    pub equation_count: usize,
    pub residual_intervals: Vec<Interval>,
    pub state: Vec<f64>,
    pub surface_coefficients: Vec<f64>,
    pub residual_upper: f64,
    pub residual_tolerance: f64,
    pub outputs: BTreeMap<String,OutputBound>,
}
// Reset certificates before each independent wave solve or verification.
// This avoids reusing bounds from an earlier wave/current scenario.
impl Default for Certificate {
    fn default() -> Self {
        Self { arithmetic:"Outward-rounded IEEE binary64 intervals".into(),
            scope:"Finite-state residual and nine selected finite-formula outputs".into(),
            root_existence_proved:false,passed:false,order:0,equation_count:0,
            residual_intervals:Vec::new(),state:Vec::new(),surface_coefficients:Vec::new(),
            residual_upper:0.0,residual_tolerance:1e-10,outputs:BTreeMap::new() }
    }
}
}

// ===========================================================================
// MODULE: solver
// ===========================================================================
// ============================================================================
// MODULE solver — RIENECKER–FENTON FOURIER COLLOCATION AND ENGINEERING OUTPUT
// ============================================================================
// All physical inputs and report quantities live in FentonStreamFunction.
// The numerical state is intentionally 1-indexed and uses 2N+10 equations.
// The solver solves at zero current and (if requested) at a prescribed current;
// each case has independent continuation, verification and diagnostics.
mod solver {
//! Steady periodic finite-amplitude waves: Fenton stream-function Fourier collocation.
//! The 1-based z vector, equations, analytic Jacobian, Newton continuation, integral
//! diagnostics and finite-system verification follow the Fourier formulation.
//!
//! The solution path is: construct nondimensional inputs -> precompute Fourier
//! collocation factors -> initialize a small-height wave -> height continuation
//! and damped Newton/SVD -> independent output verification -> report.
//! No physical quantity is silently accepted when its numerical checks fail.
use crate::linalg;
use crate::verified::{self, Certificate, Interval as I, OutputBound};
use std::collections::BTreeMap;
use std::f64::consts::PI;

// Standard gravitational acceleration [m/s^2] and nominal seawater density
// [kg/m^3]. Constants are held fixed throughout the model and post-processing.
pub const G_STD: f64=9.80665;
pub const RHO: f64=1025.0;
// Numerical acceptance controls (not merely formatting tolerances).
// Intermediate continuation may be less strict; final wave-height residual
// and output verification remain stringent. The 0.98 limiting-wave factor
// is a diagnostic comparison, not a proof of convergence or stability.
const FINAL_RESIDUAL_TOLERANCE: f64=1e-10;
const INTERMEDIATE_RESIDUAL_TOLERANCE: f64=1e-8;
const FOURIER_APPLICABILITY_FRACTION: f64=0.98;

type ResultT<T> = Result<T,String>;

#[derive(Default,Clone)]
/// One complete steady nonlinear wave problem, including internal state,
/// physical outputs, continuation controls and verification certificate.
///
/// Physical dimensional quantities are stored in SI units. Internal z[] is
/// the dimensionless Fenton state with element zero intentionally unused:
/// z[1..=2*N+10]. Solver methods use analytic derivatives and a separate
/// interval module for post-solution finite-system checks.
pub struct FentonStreamFunction {
    // ------------------------- INPUTS AND PHYSICAL CONSTANTS -------------------
    // h_target [m]: crest-to-trough wave height; t_target [s]: observed period;
    // d [m]: depth from bed to still water level; uc [m/s]: Eulerian current.
    pub h_target: f64, pub t_target: f64, pub d: f64, pub uc: f64,
    // Gravity [m/s^2] and active Fourier order N (may be smaller than
    // requested_order if the Jacobian conditioning requires order reduction).
    pub g: f64, pub n: usize,
    // ----------------------- NEWTON / CONTINUATION CONTROL ---------------------
    // nstep: starting number of height increments; number: Newton iteration
    // ceiling per height step; crit and criter_final: mean surface-correction
    // tolerances for intermediate and final stages respectively.
    pub nstep: usize, pub number: usize, pub crit: f64, pub criter_final: f64,
    // Fenton closure selector: 1 = Eulerian-current constraint (GUI default);
    // 2 = Stokes/mass-transport current constraint.
    pub current_criterion: usize,
    // ------------------------- NONDIMENSIONAL INPUTS ---------------------------
    // max_h = H/d; t_nd = T*sqrt(g/d); height = H/(g*T^2);
    // current = Uc/sqrt(g*d). 'height' is not the dimensional wave height.
    pub max_h: f64, pub t_nd: f64, pub height: f64, pub current: f64,
    // ----------------------- DISPERSION / SURFACE GEOMETRY --------------------
    // k [rad/m] = 2*pi/L, l [m] is wavelength, c [m/s] = L/T.
    pub k: f64, pub l: f64, pub c: f64,
    // eta_nodes holds bed-origin free-surface heights [m] at half-wave
    // collocation nodes; bj is the user-facing zero-based copy of B_1..B_N.
    pub eta_nodes: Vec<f64>, pub bj: Vec<f64>,
    // Crest/trough elevation relative to still-water level [m]; wave
    // steepness is H/L. The crest and trough correspond to X=0 and X=pi.
    pub eta_crest: f64, pub eta_trough: f64, pub steepness: f64,
    // Relative depth d/L; Ursell nonlinearity measure H*L^2/d^3; and a
    // qualitative shallow/intermediate/deep label based on d/L.
    pub rel_depth: f64, pub ursell: f64, pub regime: String,
    // --------------------- PHYSICAL APPLICABILITY SCREENS ----------------------
    // Miche limiting-height approximation Hmax = 0.142 L tanh(kd), and H/Hmax.
    // These are screening indicators, not an embedded breaking-wave model.
    pub breaking_limit_miche: f64, pub breaking_index: f64, pub is_breaking: bool,
    // Highest steady wave Hm estimated from Fenton (1990) Eq. 32, and H/Hm.
    pub highest_wave_limit_fenton: f64, pub highest_wave_ratio: f64,
    // 0.98*Hm diagnostic and magnitude of the last free-surface Fourier
    // coefficient. A small last coefficient is a useful, nonexclusive
    // indicator of adequate harmonic resolution.
    pub fourier_accuracy_limit: f64, pub surface_tail_coefficient: f64,
    // ||F||_infinity in the dimensionless nonlinear equations; and c-u_crest
    // [m/s], a physical near-stagnation / near-limiting-wave indicator.
    pub max_residual: f64, pub crest_stagnation_margin: f64,
    // Flags triggered by the analytical height-ratio and H/d screens.
    pub near_limiting_wave: bool, pub practical_breaking_risk: bool,
    // ---------------------------- MEAN-FLOW OUTPUTS ----------------------------
    // Eulerian current, Stokes/mass-transport current and mean fluid speed
    // in m/s. They are distinct quantities in the Fenton current closure.
    pub eulerian_current: f64, pub stokes_current: f64, pub mean_fluid_speed: f64,
    // q = Ubar*d-Q and Q, both dimensional area flows [m^2/s] per unit width.
    pub wave_volume_flux_q: f64, pub volume_flux_q: f64,
    // Full R and reduced r=R-g*d Bernoulli quantities [m^2/s^2].
    pub bernoulli_r_cap: f64, pub bernoulli_r: f64,
    // ----------------------- HYDRODYNAMIC INTEGRAL OUTPUTS ---------------------
    // Kinetic, potential and total energy density [J/m^2] (report in kJ/m^2).
    pub kinetic_energy: f64, pub potential_energy: f64, pub energy_density: f64,
    // Momentum flux [N/m] and its dimensionless depth-scaled counterpart.
    pub momentum_flux: f64, pub momentum_flux_depth: f64,
    // Radiation stress Sxx [N/m] and impulse [kg/(m*s)], with their
    // nondimensional depth-scaled representations.
    pub sxx: f64, pub sxx_depth: f64, pub impulse: f64, pub i_depth: f64,
    // Wave power [W/m], dimensionless wave-power function and energy-flux
    // velocity cg=power/energy_density [m/s] where E is sufficiently nonzero.
    pub power: f64, pub f_depth: f64, pub cg: f64,
    // Dimensionless total/kinetic/potential energy divided by rho*g*d^2.
    pub e_depth: f64, pub ke_depth: f64, pub pe_depth: f64,
    // Stokes mass transport [m/s], copied into stokes_current for reporting.
    pub mass_transport: f64,
    // ---------------------------- KINEMATICS OUTPUTS ---------------------------
    // Bed mean-square fixed-frame velocity and mean-square orbital velocity
    // [m^2/s^2]. Their distinction matters when mean current is nonzero.
    pub mean_square_bed_absolute_velocity: f64, pub mean_square_bed_velocity: f64,
    // Bed/surface maximum absolute horizontal velocities, plus the
    // crest-phase signed bed and free-surface horizontal velocities [m/s].
    pub u_bed: f64, pub u_surf: f64, pub u_bed_crest: f64, pub u_surf_crest: f64,
    // Phase [rad] of each searched maximum; acc_max is the maximum magnitude
    // of horizontal MATERIAL acceleration at the free surface [m/s^2].
    pub extrema_phases: BTreeMap<String,f64>, pub acc_max: f64,
    // Bed max |a_x| [m/s^2], max surface |w| [m/s] and velocity asymmetry.
    pub acc_bed_max: f64, pub w_max: f64, pub asymmetry: f64,
    // Optional estimated bed shear [Pa] from fixed cf=0.005 and orbital
    // excursion length [m]. Values are secondary engineering approximations.
    pub tau_bed: f64, pub excursion_bed: f64,
    // ------------------------- NUMERICAL STATE / STATUS ------------------------
    // Requested order and actual (equilibrated) Jacobian condition estimate.
    pub requested_order: usize, pub jacobian_condition: f64,
    // Accepted target-height fraction and explanatory numerical diagnostics.
    pub continuation_fraction: f64, pub numerical_note: String,
    // `converged` indicates solver acceptance, not just a small single Newton
    // update. Accuracy flags are set only after independent checks are passed.
    pub converged: bool, pub accuracy_verified: bool,
    // Residual/output interval-check flag and saved finite-system certificate.
    pub finite_system_verified: bool, pub finite_system_certificate: Certificate,
    // Named ratios of measured numerical discrepancies to tolerance budgets.
    pub accuracy_metrics: BTreeMap<String,f64>,
    // Machine-readable acceptance status and human-readable failure details.
    pub solution_status: String, pub last_error: String,
    // Private 1-based state vector, with index 0 unused.
    // ---------------------------- INTERNAL WORKSPACE ---------------------------
    // num = 2*N+10; private `z` and associated arrays use 1-based indexing.
    num: usize,
    // z: nonlinear state; tanh_values: harmonic tanh(jkd); b: B_j Fourier
    // amplitudes; y: cosine coefficients for zeta(X), all working arrays.
    z: Vec<f64>, tanh_values: Vec<f64>, b: Vec<f64>, y: Vec<f64>,
    // Cached cos/sin(j X_m). cos_nm and sin_nm flatten an (N+1) by N
    // collocation table; each value is at m*N+(j-1) in memory.
    cosa: Vec<f64>, sina: Vec<f64>, cos_nm: Vec<f64>, sin_nm: Vec<f64>,
    // Precomputed [1+exp(-2 j kd)] denominators for physical postprocessing.
    depth_denominator: Vec<f64>,
    // Current continuation-stage H/(gT^2) and H/d, distinct from targets.
    height_step: f64, hoverd: f64,
}
impl FentonStreamFunction {
    /// Construct one wave with SI inputs (H [m], T [s], d [m], current [m/s]).
    ///
    /// Allocates 2N+10 unknowns with an additional unused slot at index zero;
    /// prepares Fourier lookup tables and nondimensional input quantities.
    /// The initial Fourier order may be reduced later if the Jacobian is ill
    /// conditioned. Large nondimensional currents receive larger continuation
    /// and Newton iteration budgets, without loosening acceptance tolerances.
    pub fn new(h: f64,t: f64,depth: f64,current: f64,fourier_order: usize) -> Self {
        let n=fourier_order.max(1);
        let num=2*n+10;
        // Rust struct update syntax initializes only selected fields here;
        // unspecified numeric diagnostics use Default, ensuring fresh state.
        let mut s=Self {
            h_target:h,t_target:t,d:depth,uc:current,
            g:G_STD,n,nstep:4,number:40,crit:1e-8,criter_final:1e-10,
            current_criterion:1, requested_order:n,
            eta_nodes:vec![depth;n+1],bj:vec![0.0;n],
            num,z:vec![0.0;num+1],tanh_values:vec![0.0;n+1],
            b:vec![0.0;n+1],y:vec![0.0;num+1],
            cosa:vec![0.0;2*n+1],sina:vec![0.0;2*n+1],
            cos_nm:vec![0.0;(n+1)*n],sin_nm:vec![0.0;(n+1)*n],
            depth_denominator:vec![1.0;n+1],
            solution_status:"not_solved".into(),
            ..Self::default()
        };
        if depth>0.0 {
            s.max_h=h/depth;
            s.t_nd=t*(G_STD/depth).sqrt();
            s.height=if s.t_nd>0.0{s.max_h/(s.t_nd*s.t_nd)}else{0.0};
            s.current=current/(G_STD*depth).sqrt();
        }
        if s.current.abs()>=1.0 {s.nstep=s.nstep.max(8);s.number=s.number.max(80);}
        s.init_trig_tables();
        s
    }
    /// Access the state variable kd without exposing mutable state storage.
    /// Returns zero when the state has not been allocated or solved.
    pub fn kd_dimless(&self)->f64{self.z.get(1).copied().unwrap_or(0.0)}

    /// Solve a nonlinear wave and independently verify accepted outputs.
    ///
    /// Restarts accuracy/certificate state, attempts a solve at current N,
    /// estimates the final Jacobian condition, and may retry at lower N.
    /// When the collocation solution converges, strict output verification is
    /// performed under catch_unwind so interval domain errors are recorded.
    /// The solver never declares success solely because Newton has stopped.
    pub fn solve(&mut self) {
        self.accuracy_verified=false;
        self.finite_system_verified=false;
        self.finite_system_certificate=Certificate::default();
        self.accuracy_metrics.clear();
        self.solution_status="not_solved".into();
        // Keep the requested Fourier order for reporting. Rebuilding the
        // solver at smaller N discards old state and recomputes all stages.
        let requested=self.requested_order;
        let mut order=self.n;
        loop {
            self.solve_at_order();
            if self.converged {
                match linalg::svd_solve(&self.jacobian(),&vec![0.0;self.num],self.num) {
                    Ok((_,condition))=>self.jacobian_condition=condition,
                    Err(e)=>{self.converged=false;self.last_error=e;self.jacobian_condition=f64::INFINITY;}
                }
            }
            // Condition-based order backoff: high Fourier orders can make
            // the dense collocation Jacobian too ill-conditioned for f64.
            // Keep the same residual limits at every attempted order.
            if self.jacobian_condition<=1e11 || order<=12 {break;}
            order=(order*8/10).max(12);
            let mut next=Self::new(self.h_target,self.t_target,self.d,self.uc,order);
            next.nstep=self.nstep;next.number=self.number;
            next.crit=self.crit;next.criter_final=self.criter_final;
            next.current_criterion=self.current_criterion;
            next.requested_order=requested;
            *self=next;
        }
        if self.converged && self.jacobian_condition>1e11 {
            self.converged=false;
            self.last_error="The Fourier system is too ill-conditioned for a reliable binary64 solution.".into();
        }
        if order!=requested {
            self.numerical_note=format!("Fourier order used: {order} (requested {requested}); limited by numerical conditioning.");
        }
        if self.converged {
            // Accuracy checks include quadrature, independent identities
            // and outward interval arithmetic. catch_unwind needs the
            // release compiler flag -C panic=unwind (not panic=abort).
            // Catch interval domain/overflow panics rather than terminating the GUI worker.
            let result=std::panic::catch_unwind(std::panic::AssertUnwindSafe(||self.verify_accuracy()));
            match result {
                Ok(Ok(()))=>{},
                Ok(Err(e))=>{self.converged=false;self.accuracy_verified=false;
                    self.solution_status="accuracy_failure".into();self.last_error=format!("Accuracy verification failed: {e}");},
                Err(_)=>{self.converged=false;self.accuracy_verified=false;
                    self.solution_status="accuracy_failure".into();
                    self.last_error="Accuracy verification failed: interval enclosure domain or overflow error".into();}
            }
        } else {
            self.solution_status="nonconvergence".into();
        }
    }
    /// Validate the physical inputs and solve exactly one active Fourier order.
    /// Failure is captured in `last_error` rather than exposing an invalid
    /// dimensional report. Both selected current-closure options are supported.
    fn solve_at_order(&mut self) {
        self.converged=false;
        self.last_error.clear();
        if !(self.h_target.is_finite()&&self.t_target.is_finite()&&self.d.is_finite()&&self.uc.is_finite()
            &&self.h_target>0.0&&self.t_target>0.0&&self.d>0.0) {
            self.last_error="Invalid inputs: H, T, and d must be finite and > 0; current must be finite.".into();
            return;
        }
        self.max_residual=f64::INFINITY;
        if self.current_criterion!=1 && self.current_criterion!=2 {
            self.last_error="Current criterion must be 1 (Eulerian) or 2 (Stokes).".into();return;
        }
        match self.solve_internal() {
            Ok(())=>self.converged=true,
            Err(e)=>{self.converged=false;self.last_error=format!("Solver error: {e}");}
        }
    }
    /// Precompute sin/cos for the half-wave collocation grid.
    /// X_m = m*pi/N, m=0..N, j=1..N. Because m*j is an integer, modulo 2N
    /// lookup avoids redundant trigonometric calls within every Newton step.
    fn init_trig_tables(&mut self) {
        for i in 0..=2*self.n {
            let a=i as f64*PI/self.n as f64;
            self.cosa[i]=a.cos();self.sina[i]=a.sin();
        }
        for m in 0..=self.n {
            for j in 1..=self.n {
                let nm=m*j%(2*self.n);
                self.cos_nm[m*self.n+j-1]=self.cosa[nm];
                self.sin_nm[m*self.n+j-1]=self.sina[nm];
            }
        }
    }
    /// Seed the first continuation step with a small-amplitude approximation.
    ///
    /// sigma = 2*pi*sqrt([H/(gT^2)]/[H/d]) and the Fenton–McKee-like kd
    /// estimate provide only a starting point for the fully nonlinear solve.
    /// Initialize phase speed, current variables, Bernoulli variable, free-
    /// surface nodes and Fourier coefficients consistently with this guess.
    fn init_linear(&mut self) {
        let sigma=if self.hoverd>0.0 {2.0*PI*(self.height_step/self.hoverd).sqrt()}else{0.0};
        if sigma>0.0 {
            let t=sigma.powf(1.5).tanh();
            self.z[1]=sigma*sigma/t.powf(2.0/3.0);
        } else {
            self.z[1]=2.0*PI*self.height_step.max(1e-12)/self.hoverd.max(1e-12);
        }
        // z2=kH; z3=T*sqrt(gk); z4=c*sqrt(k/g).
        // The nonlinear solver subsequently enforces z3*z4=2*pi exactly.
        self.z[2]=self.z[1]*self.hoverd;
        self.z[4]=self.z[1].tanh().sqrt();
        self.z[3]=2.0*PI/self.z[4];
        if self.current_criterion==1 {
            self.z[5]=self.current*self.z[1].sqrt();self.z[6]=0.0;
        }else{self.z[6]=self.current*self.z[1].sqrt();self.z[5]=0.0;}
        self.z[7]=self.z[4];self.z[8]=0.0;self.z[9]=0.5*self.z[7]*self.z[7];
        // Seed crest/trough Fourier-surface nodes from a sinusoidal wave.
        // The first B_j is set to an Airy-like harmonic estimate; all higher
        // harmonics initially vanish. Newton computes the nonlinear harmonics.
        self.z[10]=0.5*self.z[2];
        for i in 1..=self.n {
            self.z[self.n+i+10]=0.0;
            self.z[i+10]=0.5*self.z[2]*self.cosa[i];
        }
        self.z[self.n+11]=0.5*self.z[2]/self.z[7];
    }
    /// Convenience wrapper for the nonlinear equation assembly.
    /// Returns the 1-based residual vector and its squared Euclidean norm.
    /// A reusable-buffer implementation below avoids repeated allocations
    /// during Newton line searches.
    fn eqns_for_state(&self,s:&[f64])->ResultT<(Vec<f64>,f64)> {
        let mut r=vec![0.0;self.num+1];
        let ss=self.eqns_for_state_into(s,&mut r)?;
        Ok((r,ss))
    }
    /// Assemble the exact finite nonlinear residual equations F(z).
    ///
    /// Rows 1..8 enforce depth, height, period, current and mean-level
    /// relations. Rows 9..N+9 give N+1 streamline equations; rows N+10..
    /// 2N+10 give the corresponding Bernoulli equations. The returned scalar
    /// is sum(F_i^2), used as the line-search merit function. Exponentials
    /// with unsafe arguments are rejected to detect Newton divergence.
    fn eqns_for_state_into(&self,s:&[f64],r:&mut [f64])->ResultT<f64> {
        let n=self.n;
        debug_assert_eq!(r.len(),self.num+1);
        r.fill(0.0);
        // The denominator depends on kd, not on the collocation phase.
        // S_j and C_j share the phase-independent finite-depth denominator
        // 1 + exp(-2*j*kd). Compute it once for this candidate state.
        let mut den=vec![0.0;n+1];
        for j in 1..=n {den[j]=1.0+(-2.0*j as f64*s[1]).exp();}
        // Global constraints F1..F8 (see source-level theoretical reference):
        // relative depth, period, wave speed, current closure, flux, zero-mean
        // surface elevation and prescribed crest-to-trough wave height.
        r[1]=s[2]-s[1]*self.hoverd;
        r[2]=s[2]-self.height_step*s[3]*s[3];
        r[3]=s[4]*s[3]-2.0*PI;
        r[4]=s[5]+s[7]-s[4];
        r[5]=s[1]*(s[6]+s[7]-s[4])-s[8];
        r[6]=s[self.current_criterion+4]-self.current*s[1].sqrt();
        r[7]=s[10]+s[n+10];
        for i in 1..n {r[7]+=2.0*s[10+i];}
        r[8]=s[10]-s[n+10]-s[2];
        // Enforce the two nonlinear free-surface conditions at each of the
        // N+1 nodes from crest (m=0) to trough (m=N).
        for m in 0..=n {
            let zsurf=s[10+m];
            let(mut psi,mut u,mut v)=(0.0,0.0,0.0);
            for j in 1..=n {
                let x=j as f64*zsurf;
                // Refuse an unphysically large exponential argument in a
                // trial iterate. A line search can then try a smaller step.
                if x>60.0 || x< -60.0 {return Err("Divergence: exp(j*zsurf) out of safe range.".into());}
                let ep=x.exp();
                let em=(-x-2.0*j as f64*s[1]).exp();
                let sh=(ep-em)/den[j];
                let ch=(ep+em)/den[j];
                let cos=self.cos_nm[m*n+j-1];let sin=self.sin_nm[m*n+j-1];
                let cj=s[n+j+10];
                psi+=cj*sh*cos;
                u+=j as f64*cj*ch*cos;
                v+=j as f64*cj*sh*sin;
            }
            // Kinematic/streamline condition and dynamic Bernoulli condition.
            // Both are expressed in the nondimensional moving-wave frame.
            r[m+9]=psi-s[8]-s[7]*s[m+10];
            r[n+m+10]=0.5*((-s[7]+u)*(-s[7]+u)+v*v)+s[m+10]-s[9];
        }
        let ss=(1..=self.num).map(|i|r[i]*r[i]).sum();
        Ok(ss)
    }
    /// Infinity norm ||F||_inf = max_i |F_i| over active equations only.
    /// The final acceptance threshold is independent of wave height/order.
    fn residual_inf_norm(&self,r:&[f64])->f64 {
        (1..=self.num).fold(0.0_f64,|v,i|v.max(r[i].abs()))
    }
    /// Analytic dense Jacobian J_ij = partial F_i / partial z_j.
    ///
    /// The row-major n-by-n matrix uses a 1-based row/column accessor `ix`.
    /// Scalar constraints are differentiated directly. For each free-surface
    /// phase, S_j and C_j plus their derivatives wrt kd and zeta give the
    /// Jacobian columns for wavenumber, surface shape, Bernoulli and Fourier
    /// coefficients. Analytic derivatives avoid expensive finite-difference
    /// evaluations and their step-size sensitivity.
    fn jacobian(&self)->Vec<f64> {
        let n=self.n;let num=self.num;let z=&self.z;
        let mut a=vec![0.0;num*num];
        // Translate documented 1-based equation/state indices to Rust's
        // 0-based contiguous matrix layout; do not swap row/column order.
        let ix=|r:usize,c:usize|(r-1)*num+c-1;
        a[ix(1,1)]=-self.hoverd;a[ix(1,2)]=1.0;
        a[ix(2,2)]=1.0;a[ix(2,3)]=-2.0*self.height_step*z[3];
        a[ix(3,3)]=z[4];a[ix(3,4)]=z[3];
        a[ix(4,5)]=1.0;a[ix(4,7)]=1.0;a[ix(4,4)]=-1.0;
        a[ix(5,1)]=z[6]+z[7]-z[4];
        a[ix(5,6)]=z[1];a[ix(5,7)]=z[1];a[ix(5,4)]=-z[1];a[ix(5,8)]=-1.0;
        a[ix(6,self.current_criterion+4)]=1.0;
        a[ix(6,1)]=-0.5*self.current/z[1].sqrt();
        for m in 0..=n {a[ix(7,10+m)]=if m==0 || m==n {1.0}else{2.0};}
        a[ix(8,10)]=1.0;a[ix(8,n+10)]=-1.0;a[ix(8,2)]=-1.0;
        // Precompute exp(-2jkd) for all j and reuse in every collocation row.
        // This reduces repeated transcendental evaluations in the Jacobian.
        let mut depth_q=vec![0.0;n+1];
        let mut depth_den=vec![0.0;n+1];
        for j in 1..=n {
            let q=(-2.0*j as f64*z[1]).exp();
            depth_q[j]=q;depth_den[j]=1.0+q;
        }
        let mut sh=vec![0.0;n];
        let mut ch=vec![0.0;n];
        for m in 0..=n {
            let sr=9+m;let dr=n+10+m;let surface=z[10+m];
            let(mut u,mut v,mut duk,mut dvk,mut duy,mut dvy)=(0.0,0.0,0.0,0.0,0.0,0.0);
            for j in 1..=n {
                let jf=j as f64;
                let idx=m*n+j-1;
                let co=self.cos_nm[idx];let si=self.sin_nm[idx];let bj=z[n+10+j];
                let ep=(jf*surface).exp();let q=depth_q[j];
                let em=(-jf*surface-2.0*jf*z[1]).exp();let den=depth_den[j];
                let sj=(ep-em)/den;let cj=(ep+em)/den;
                // Differentiate the stable hyperbolic basis wrt kd.
                // ds = dS_j/d(kd), dc = dC_j/d(kd), with zeta fixed.
                let ds=2.0*jf*(em+q*ep)/(den*den);
                let dc=2.0*jf*(q*ep-em)/(den*den);
                sh[j-1]=sj;ch[j-1]=cj;
                u+=jf*bj*cj*co;v+=jf*bj*sj*si;
                a[ix(sr,1)]+=bj*ds*co;
                duk+=jf*bj*dc*co;dvk+=jf*bj*ds*si;
                duy+=jf*jf*bj*sj*co;dvy+=jf*jf*bj*cj*si;
                a[ix(sr,n+10+j)]=sj*co;
            }
            // In the moving frame the horizontal velocity is u-z7; the
            // Bernoulli residual contains 0.5[(u-z7)^2+v^2] + zeta-z9.
            // Apply the chain rule to all Jacobian contributions.
            let uw=u-z[7];
            a[ix(sr,7)]=-surface;a[ix(sr,8)]=-1.0;a[ix(sr,10+m)]=uw;
            a[ix(dr,1)]=uw*duk+v*dvk;
            a[ix(dr,7)]=-uw;a[ix(dr,9)]=-1.0;
            a[ix(dr,10+m)]=1.0+uw*duy+v*dvy;
            for j in 1..=n {
                let idx=m*n+j-1;let jf=j as f64;
                a[ix(dr,n+10+j)]=jf*(uw*ch[j-1]*self.cos_nm[idx]+v*sh[j-1]*self.sin_nm[idx]);
            }
        }
        a
    }
    /// Perform one damped Newton/SVD correction to the Fenton state.
    ///
    /// Compute F, J and delta=-J^+F, then backtrack alpha from 1 by halves.
    /// Only trials with kd>0, finite entries and non-increasing sum(F_i^2)
    /// are accepted. Return mean absolute correction of the N+1 free-surface
    /// ordinates, which is combined with the independent residual criterion.
    fn newton_step(&mut self)->ResultT<f64> {
        let mut r0=vec![0.0;self.num+1];
        let ss0=self.eqns_for_state_into(&self.z,&mut r0)?;
        if !ss0.is_finite(){return Err("Non-finite initial residual.".into());}
        let z0=self.z.clone();
        let a=self.jacobian();
        let rhs=(1..=self.num).map(|i|-r0[i]).collect::<Vec<_>>();
        let (dx,condition)=linalg::svd_solve(&a,&rhs,self.num)?;
        self.jacobian_condition=condition;
        if dx.iter().any(|v|!v.is_finite()) {return Err("Non-finite Newton step.".into());}
        let mut best=z0.clone();
        let mut trial=z0.clone();
        let mut workspace=vec![0.0;self.num+1];
        // Armijo-like monotone residual safeguard without an additional
        // slope parameter: reject each trial if its squared residual grows.
        // The smallest attempted factor is approximately 1e-4.
        let mut alpha=1.0;
        while alpha>=1e-4 {
            trial.copy_from_slice(&z0);
            for i in 1..=self.num {trial[i]+=alpha*dx[i-1];}
            if trial[1]>0.0 && trial.iter().all(|v|v.is_finite()) {
                if let Ok(ss)=self.eqns_for_state_into(&trial,&mut workspace) {
                    if ss.is_finite() && ss<=ss0 {best.copy_from_slice(&trial);break;}
                }
            }
            alpha*=0.5;
        }
        let correction=(10..=self.n+10).map(|i|(best[i]-z0[i]).abs()).sum::<f64>()/(self.n+1) as f64;
        self.z=best;
        Ok(correction)
    }
    /// Reconstruct cosine coefficients of zeta(X) from collocation values.
    ///
    /// The discrete symmetric cosine transformation handles crest/trough
    /// endpoint weights separately and stores Y_j in `y[j]`; original stream-
    /// function coefficients are copied to `b[j]`. Hyperbolic denominator
    /// factors are refreshed for later kinematic and integral evaluations.
    fn compute_y_and_b(&mut self) {
        self.y.fill(0.0);
        for j in 1..=self.n {
            self.b[j]=self.z[j+self.n+10];
            let mut s=0.5*(self.z[10]+self.z[self.n+10]*(if j%2==0{1.0}else{-1.0}));
            for m in 1..self.n {
                let idx=m*j%(2*self.n);
                s+=self.z[10+m]*self.cosa[idx];
            }
            self.y[j]=2.0*s/self.n as f64;
        }
        for i in 1..=self.n {
            self.tanh_values[i]=(i as f64*self.z[1]).tanh();
            self.depth_denominator[i]=1.0+(-2.0*i as f64*self.z[1]).exp();
        }
    }
    /// Reconstruct k*eta(X) from the accepted finite cosine series.
    /// The final mode N receives the half weight dictated by the chosen
    /// symmetric finite Fourier convention.
    fn surface_keta(&self,x:f64)->f64 {
        let mut eta=0.0;
        for j in 1..self.n {eta+=self.y[j]*(j as f64*x).cos();}
        eta+=0.5*self.y[self.n]*(self.n as f64*x).cos();
        eta
    }
    /// Returns depth-scaled fixed-frame (u,w) and depth-scaled material a_x/g.
    /// Return depth-scaled fixed-frame (u,w,a_x) at dimensionless (X,Y).
    ///
    /// S_j and C_j provide the velocity field and its spatial derivatives.
    /// Compute local time derivative from travelling-wave phase (-c*u_x),
    /// then add advective terms (u*u_x+w*u_y) for MATERIAL acceleration.
    /// The caller later scales velocity by sqrt(gd) and acceleration by g.
    fn point(&self,x:f64,yloc:f64)->(f64,f64,f64){
        let kd=self.z[1];let kd_sqrt=kd.sqrt();
        let c_nd=self.z[4]/kd_sqrt;let ce_nd=self.z[5]/kd_sqrt;
        let(mut u,mut v,mut ux,mut vx)=(0.0,0.0,0.0,0.0);
        for j in 1..=self.n {
            let jf=j as f64;
            let cos=(jf*x).cos();let sin=(jf*x).sin();
            let a=jf*yloc;let b=jf*kd;
            let ep=a.exp();let em=(-a-2.0*b).exp();let den=self.depth_denominator[j];
            let ch=(ep+em)/den;let sh=(ep-em)/den;
            let bj=self.b[j];
            u+=jf*bj*ch*cos;
            v+=jf*bj*sh*sin;
            ux+=-jf*jf*bj*ch*sin;
            vx+= jf*jf*bj*sh*cos;
        }
        // Convert harmonic velocity to depth-scaled units and add the
        // Eulerian-current component. Derivatives scale differently because
        // dimensionless phase X contains the wave number k.
        u=u/kd_sqrt+ce_nd;v/=kd_sqrt;ux*=kd_sqrt;vx*=kd_sqrt;
        let ut=-c_nd*ux;let uy=vx;
        (u,v,ut+u*ux+v*uy)
    }
    /// Dimensional (u,w,a_x) with coordinate z measured upwards from seabed.
    /// Physical kinematics at elevation z_bed [m] above the seabed and
    /// phase X [rad]. Returns (u,w,a_x) with units (m/s,m/s,m/s^2).
    /// `z_bed=0` selects the bed; the free surface is at d+eta(X).
    pub fn get_kinematics(&self,z_bed:f64,phase:f64)->(f64,f64,f64) {
        if self.z[1]<=0.0 || self.d<=0.0 {return (0.0,0.0,0.0);}
        let kd=self.z[1];let k_phys=kd/self.d;
        let yloc=k_phys*(z_bed-self.d);
        let (u,v,ax)=self.point(phase,yloc);
        let scale_v=(self.g*self.d).sqrt();
        (u*scale_v,v*scale_v,ax*self.g)
    }
    /// Mean square *orbital* bed velocity [m^2/s^2] from Fourier orthogonality.
    ///
    /// The constant Eulerian current is subtracted before phase averaging:
    /// <(u_b-Uc)^2> = (g*d/(2kd)) sum_j (j B_j sech(jkd))^2.
    /// This avoids an expensive phase-quadrature loop and removes the mean-
    /// current contribution from the orbital RMS velocity.
    fn mean_square_bed_orbital_velocity(&self)->f64{
        if self.d<=0.0 || self.t_target<=0.0 || self.z[1]<=0.0 {return 0.0;}
        let mut sum=0.0;
        for j in 1..=self.n {
            let e=(-(j as f64)*self.z[1]).exp();
            let sech=2.0*e/(1.0+e*e);
            let amplitude=j as f64*self.b[j]*sech;
            sum+=amplitude*amplitude;
        }
        0.5*self.g*self.d/self.z[1]*sum
    }
    /// Depth-normalized momentum flux S/(rho*g*d^2) at a selected phase.
    /// The integrand combines moving-frame pressure and relative horizontal
    /// momentum flux; the vertical integration uses a trapezoidal rule with
    /// at least 50 points (the normal solver requests 1200).
    fn momentum_flux_s_depth(&self,phase:f64,npts:usize)->f64{
        let kd=self.z[1];if kd<=0.0{return 0.0;}
        let c_nd=self.z[4]/kd.sqrt();
        let r_nd=1.0+self.z[9]/kd;
        let keta=self.surface_keta(phase);
        let eta_over_d=1.0+keta/kd;
        let np=npts.max(50);
        let dy=eta_over_d/(np-1) as f64;
        let mut integ=0.0;
        for i in 0..np {
            let y=dy*i as f64;let yloc=kd*(y-1.0);
            let (u,v,_ax)=self.point(phase,yloc);
            let urel=u-c_nd;let p=r_nd-y-0.5*(urel*urel+v*v);
            let f=p+urel*urel;
            integ+= (if i==0 || i==np-1 {0.5}else{1.0})*f;
        }
        integ*dy
    }
    /// Derive wave-integrated engineering quantities from the Fourier state.
    ///
    /// Convert depth-scaled kinetic/potential energy, impulse, radiation
    /// stress and wave power to SI units. Reconstruct the current variables,
    /// Bernoulli offsets and orbital bed velocity. These analytic spectral
    /// expressions are checked later against independent defining integrals.
    fn calc_integral_props(&mut self) {
        let kd=self.z[1];
        if kd<=0.0 {self.power=0.0;self.energy_density=0.0;self.sxx=0.0;
            self.impulse=0.0;self.cg=0.0;self.mass_transport=0.0;self.bernoulli_r_cap=0.0;return;}
        let kd_sqrt=kd.sqrt();
        let ce_dimless=self.z[5]/kd_sqrt;
        let cs_dimless=self.z[6]/kd_sqrt;
        let ubar_dimless=self.z[7]/kd_sqrt;
        let kd_32=kd*kd_sqrt;
        // Flux variables depend on (kd)^(3/2). The original Fenton notation
        // distinguishes wave transport and total volume flux; the report
        // reproduces both rather than treating them as interchangeable.
        let q_dimless=ubar_dimless-self.z[8]/kd_32;
        let r_dimless=1.0+self.z[9]/kd;
        let pulse=self.z[8]+kd*self.z[5];
        let ke=0.5*(self.z[4]*pulse-self.z[5]*q_dimless*kd_32);
        let mut pe=0.0;
        for i in 1..self.n {pe+=0.25*self.y[i]*self.y[i];}
        pe+=self.y[self.n]*self.y[self.n]/16.0;
        let beta=2.0*self.z[9]-self.z[4]*self.z[4]+2.0*self.z[4]*self.z[5];
        let sxx=4.0*ke-3.0*pe+beta*kd-2.0*self.z[5]*pulse;
        let f=self.z[4]*(3.0*ke-2.0*pe)+0.5*beta*(pulse+self.z[4]*kd)-2.0*self.z[4]*self.z[5]*pulse;
        let kd2=kd*kd;let kd25=kd2*kd_sqrt;
        self.e_depth=(ke+pe)/kd2;self.ke_depth=ke/kd2;self.pe_depth=pe/kd2;
        self.sxx_depth=sxx/kd2;self.f_depth=f/kd25;self.i_depth=pulse/kd_32;
        // Dimensional scales: rho*g*d^2 for E,Sxx; rho*g^(3/2)*d^(5/2)
        // for energy flux per crest width; rho*sqrt(g*d^3) for impulse.
        self.energy_density=RHO*self.g*self.d*self.d*self.e_depth;
        self.sxx=RHO*self.g*self.d*self.d*self.sxx_depth;
        self.power=RHO*self.g.powf(1.5)*self.d.powf(2.5)*self.f_depth;
        self.momentum_flux_depth=self.momentum_flux_s_depth(0.0,1200);
        self.momentum_flux=RHO*self.g*self.d*self.d*self.momentum_flux_depth;
        self.impulse=RHO*(self.g*self.d.powi(3)).sqrt()*self.i_depth;
        self.bernoulli_r_cap=r_dimless*self.g*self.d;
        self.mass_transport=cs_dimless*(self.g*self.d).sqrt();
        self.eulerian_current=ce_dimless*(self.g*self.d).sqrt();
        self.stokes_current=self.mass_transport;
        self.mean_fluid_speed=ubar_dimless*(self.g*self.d).sqrt();
        self.volume_flux_q=q_dimless*(self.g*self.d.powi(3)).sqrt();
        self.wave_volume_flux_q=self.mean_fluid_speed*self.d-self.volume_flux_q;
        self.bernoulli_r=(r_dimless-1.0)*self.g*self.d;
        self.kinetic_energy=RHO*self.g*self.d*self.d*self.ke_depth;
        self.potential_energy=RHO*self.g*self.d*self.d*self.pe_depth;
        // Absolute bed mean-square speed adds Uc^2 to the Fourier orbital
        // contribution because the phase-average orbital velocity is zero.
        self.mean_square_bed_velocity=self.mean_square_bed_orbital_velocity();
        self.mean_square_bed_absolute_velocity=self.mean_square_bed_velocity+self.eulerian_current*self.eulerian_current;
        self.cg=if self.energy_density.abs()>1e-12 {self.power/self.energy_density}else{0.0};
    }
    /// Evaluate |quantity| at bed or instantaneous free surface and phase X.
    /// The named quantity may be horizontal/vertical velocity, horizontal
    /// material acceleration, or orbital bed speed with Eulerian mean removed.
    fn phase_quantity(&self,phase:f64,quantity:&str)->ResultT<f64>{
        let bed=quantity=="bed"||quantity=="bed_orbital"||quantity=="bed_acceleration";
        let top=if bed{0.0}else{self.d*(1.0+self.surface_keta(phase)/self.z[1])};
        let(u,w,acc)=self.get_kinematics(top,phase);
        let value=match quantity {
            "bed_orbital"=>u-self.eulerian_current,
            "vertical"=>w,
            "acceleration"|"bed_acceleration"=>acc,
            _=>u,
        };
        if !value.is_finite() {return Err("Non-finite phase evaluation".into());}
        Ok(value.abs())
    }
    /// Search for the absolute maximum over the symmetric half-period 0..pi.
    ///
    /// Sample at >=257 phases, identify candidate discrete local maxima,
    /// refine each candidate by golden-section optimization, and also retain
    /// sampled endpoints. Return (peak absolute value, phase in radians).
    /// This is a numerical extremum search, not a mathematical global proof.
    fn phase_maximum(&self,quantity:&str)->ResultT<(f64,f64)>{
        let count=(8*self.n+1).max(257);
        let mut phases=vec![0.0;count];let mut values=vec![0.0;count];
        let(mut best,mut phase)=(-1.0,0.0);
        for i in 0..count {
            phases[i]=PI*i as f64/(count-1) as f64;
            values[i]=self.phase_quantity(phases[i],quantity)?;
            if values[i]>best {best=values[i];phase=phases[i];}
        }
        // Golden-section ratio yields nested brackets without derivatives.
        // Continue until the phase bracket is shorter than 1e-12 rad or
        // 64 refinement iterations are exhausted.
        let golden=(5.0_f64.sqrt()-1.0)/2.0;
        for i in 1..count-1 {
            if !(values[i]>=values[i-1]&&values[i]>=values[i+1]
                &&(values[i]>values[i-1]||values[i]>values[i+1])) {continue;}
            let(mut left,mut right)=(phases[i-1],phases[i+1]);
            let mut x1=right-golden*(right-left);
            let mut x2=left+golden*(right-left);
            let mut f1=self.phase_quantity(x1,quantity)?;
            let mut f2=self.phase_quantity(x2,quantity)?;
            for _ in 0..64 {
                if right-left<1e-12 {break;}
                if f1<f2 {left=x1;x1=x2;f1=f2;
                    x2=left+golden*(right-left);f2=self.phase_quantity(x2,quantity)?;}
                else {right=x2;x2=x1;f2=f1;
                    x1=right-golden*(right-left);f1=self.phase_quantity(x1,quantity)?;}
                if f1>best {best=f1;phase=x1;}
                if f2>best {best=f2;phase=x2;}
            }
        }
        Ok((best,phase))
    }
    /// Populate peak bed/surface speeds and accelerations, plus the phase of
    /// each extremum. Separate bed and surface evaluations prevent accidental
    /// comparison of kinematics at different vertical elevations.
    fn compute_extrema(&mut self)->ResultT<()> {
        self.extrema_phases.clear();
        for key in ["bed","surface","vertical","acceleration","bed_acceleration"] {
            let(value,phase)=self.phase_maximum(key)?;
            match key {
                "bed"=>self.u_bed=value,"surface"=>self.u_surf=value,
                "vertical"=>self.w_max=value,"acceleration"=>self.acc_max=value,
                "bed_acceleration"=>self.acc_bed_max=value,_=>{}
            }
            self.extrema_phases.insert(key.into(),phase);
        }
        Ok(())
    }
    /// Wavelength/depth-dependent highest steady wave ratio Hm/d.
    /// Implements the rational fit attributed to Fenton (1990), Eq. (32).
    /// For x=L/d>1, evaluate in reciprocal powers to reduce overflow risk.
    /// It is an applicability diagnostic, not a breaking simulation.
    fn fenton_highest_hoverd(x:f64)->f64 {
        if !x.is_finite() || x<=0.0 {return f64::NAN;}
        if x<=1.0 {
            let numerator=x*(0.141063+x*(0.0095721+0.0077829*x));
            let denominator=1.0+x*(0.0788340+x*(0.0317567+0.0093407*x));
            numerator/denominator
        }else{
            let inv=1.0/x;let inv2=inv*inv;let inv3=inv2*inv;
            let numerator=0.0077829+0.0095721*inv+0.141063*inv2;
            let denominator=0.0093407+0.0317567*inv+0.0788340*inv2+inv3;
            numerator/denominator
        }
    }
    /// Adaptive wave-height continuation and final dimensional postprocessing.
    ///
    /// Fraction starts at zero, advances toward one, and only accepted states
    /// are retained. For a new stage, extrapolate from the previous two states
    /// when available; otherwise restore the last accepted state. Failure
    /// halves the height increment down to a fixed floor and attempt budget.
    /// After the final stage, recompute the residual, rebuild dimensional
    /// kinematics/fluxes and verify normalized state consistency.
    fn solve_internal(&mut self)->ResultT<()> {
        let(mut fraction,mut increment,mut previous_fraction)=(0.0,1.0/self.nstep.max(1) as f64,0.0);
        let mut accepted:Vec<f64>=Vec::new();
        let mut previous:Vec<f64>=Vec::new();
        let mut attempts=0;
        // The continuation parameter multiplies H/(gT^2) and H/d together,
        // so the physical wave grows while water depth and period remain fixed.
        while fraction<1.0 {
            attempts+=1;
            let target=1.0_f64.min(fraction+increment);
            self.height_step=target*self.height;self.hoverd=target*self.max_h;
            if accepted.is_empty() {self.init_linear();}
            else {
                self.z=accepted.clone();
                // Secant extrapolation in the continuation parameter uses
                // two previously accepted states to improve initial guesses.
                if !previous.is_empty() {
                    let ratio=(target-fraction)/(fraction-previous_fraction);
                    for i in 1..=self.num {self.z[i]+=ratio*(accepted[i]-previous[i]);}
                }
                if self.z[1]<=0.0 || self.z.iter().any(|v|!v.is_finite()) {self.z=accepted.clone();}
            }
            let mut step_converged=false;
            // At the exact target wave height, enforce final 1e-10 residual
            // and surface correction requirements. Earlier steps use 1e-8.
            let residual_limit=if target==1.0{FINAL_RESIDUAL_TOLERANCE}else{INTERMEDIATE_RESIDUAL_TOLERANCE};
            let correction_limit=if target==1.0{self.criter_final}else{self.crit};
            for it in 1..=self.number {
                let err=match self.newton_step(){Ok(v)=>v,Err(_)=>break};
                let (rhs,ss)=match self.eqns_for_state(&self.z){Ok(v)=>v,Err(_)=>break};
                if !ss.is_finite() {break;}
                self.max_residual=self.residual_inf_norm(&rhs);
                if it>1 && err<correction_limit*self.z[1].abs() && self.max_residual<=residual_limit {
                    step_converged=true;break;
                }
                if err==0.0 && self.max_residual>residual_limit {break;}
            }
            // A failure in Newton, SVD, finite residual or iteration budget
            // triggers smaller continuation steps; no unverified state is used.
            if !step_converged {
                increment*=0.5;
                if increment<1.0/4096.0 || attempts>=80 {
                    return Err(format!("No converged solution at H={} m, T={} s, Uc={} m/s. Height continuation reached H={} m; max|F|={:.3e}. The requested wave may be near or beyond the limiting wave height; reduce H or increase T.",
                       self.h_target,self.t_target,self.uc,fraction*self.h_target,self.max_residual));
                }
                continue;
            }
            previous_fraction=fraction;previous=accepted;
            fraction=target;accepted=self.z.clone();self.continuation_fraction=fraction;
            increment=(increment*1.5).min(1.0/self.nstep.max(1) as f64);
            self.compute_y_and_b();
        }
        // First final full-equation residual check, before any dimensional
        // consistency normalization or engineering output computation.
        let (r,ss)=self.eqns_for_state(&self.z)?;
        if !ss.is_finite(){return Err("Final nonlinear residual norm is non-finite".into());}
        self.max_residual=self.residual_inf_norm(&r);
        if self.max_residual>FINAL_RESIDUAL_TOLERANCE {
            return Err(format!("False convergence rejected: maximum residual = {:.3e}, tolerance = {:.3e}.",
                 self.max_residual,FINAL_RESIDUAL_TOLERANCE));
        }
        // Highest retained cosine coefficient |E_N| is a spectral-resolution
        // indicator; independent physical/accuracy checks remain necessary.
        self.surface_tail_coefficient=self.y[self.n].abs();
        let kd=self.z[1];
        if kd<=0.0 || !kd.is_finite() {return Err("Invalid wavenumber (kd)".into());}
        // Nonlinear wavenumber kd -> physical k, L and phase speed c.
        // These values are fixed by the converged nonlinear period relation.
        let k_phys=kd/self.d;let l_phys=2.0*PI/k_phys;let c_phys=l_phys/self.t_target;
        if l_phys<=0.0 || !l_phys.is_finite(){return Err("Invalid wavelength".into());}
        if c_phys<=0.0 || !c_phys.is_finite(){return Err("Invalid celerity".into());}
        for m in 0..=self.n {self.eta_nodes[m]=self.d*(1.0+self.z[10+m]/kd);}
        self.k=k_phys;self.l=l_phys;self.c=c_phys;
        // Rebuild the period/phase-current state from dimensional values as a
        // consistency normalization. This must preserve F(z) within tolerance.
        let z3_period=self.t_target*(self.g*k_phys).sqrt();
        let z4_period=c_phys*(k_phys/self.g).sqrt();
        if z3_period.is_finite()&&z4_period.is_finite()&&z3_period>0.0&&z4_period>0.0 {
            self.z[3]=z3_period;self.z[4]=z4_period;
            if self.current_criterion==1 {
                self.z[5]=self.uc*(k_phys/self.g).sqrt();
                self.z[7]=self.z[4]-self.z[5];
                self.z[6]=self.z[8]/self.z[1]-self.z[7]+self.z[4];
            }else{
                self.z[6]=self.uc*(k_phys/self.g).sqrt();
                self.z[7]=self.z[8]/self.z[1]-self.z[6]+self.z[4];
                self.z[5]=self.z[4]-self.z[7];
            }
        }
        let (r,ss)=self.eqns_for_state(&self.z)?;
        if !ss.is_finite(){return Err("Non-finite residual after dimensional consistency normalization".into());}
        self.max_residual=self.residual_inf_norm(&r);
        if self.max_residual>FINAL_RESIDUAL_TOLERANCE {
            // Refuse an apparently converged wave if numerical normalization
            // disturbed the nonlinear equations beyond the original tolerance.
            return Err("Dimensional consistency normalization invalidated the nonlinear solution.".into());
        }
        for j in 1..=self.n {self.bj[j-1]=self.b[j];}
        self.eta_crest=self.eta_nodes[0]-self.d;
        self.eta_trough=self.eta_nodes[self.n]-self.d;
        // Compact regular-wave shape diagnostics: H/L, d/L and Ursell number.
        self.steepness=self.h_target/self.l;
        self.rel_depth=self.d/self.l;
        self.ursell=self.h_target*self.l*self.l/self.d.powi(3);
        self.regime=if self.rel_depth<0.05{"Shallow"}else if self.rel_depth<0.5{"Intermediate"}else{"Deep"}.into();
        let highest_hoverd=Self::fenton_highest_hoverd(self.l/self.d);
        if highest_hoverd<=0.0 || !highest_hoverd.is_finite(){return Err("Could not evaluate the theoretical highest-wave limit".into());}
        self.highest_wave_limit_fenton=highest_hoverd*self.d;
        self.fourier_accuracy_limit=FOURIER_APPLICABILITY_FRACTION*self.highest_wave_limit_fenton;
        self.highest_wave_ratio=self.h_target/self.highest_wave_limit_fenton;
        self.near_limiting_wave=self.highest_wave_ratio>=0.90;
        self.practical_breaking_risk=self.max_h>0.55;
        // Miche's regular-wave limiting approximation, independent of the
        // Fenton highest-steady-wave rational fit above. Neither supersedes
        // a detailed engineering assessment of near-breaking wave behavior.
        self.breaking_limit_miche=0.142*self.l*(self.k*self.d).tanh();
        self.breaking_index=if self.breaking_limit_miche>0.0{self.h_target/self.breaking_limit_miche}else{0.0};
        self.is_breaking=self.breaking_limit_miche>0.0&&self.h_target>self.breaking_limit_miche;
        self.calc_integral_props();
        self.u_bed_crest=self.get_kinematics(0.0,0.0).0;
        self.u_surf_crest=self.get_kinematics(self.d+self.eta_crest,0.0).0;
        let u_trough=self.get_kinematics(self.d+self.eta_trough,PI).0;
        self.asymmetry=if u_trough.abs()>0.0 {(self.u_surf_crest/u_trough).abs()}else{0.0};
        self.crest_stagnation_margin=self.c-self.u_surf_crest;
        self.compute_extrema()?;
        // Secondary bed-shear indicator with fixed empirical cf=0.005:
        // tau_bed = 0.5*rho*cf*(peak bed speed)^2. Not a calibrated drag law.
        self.tau_bed=0.5*RHO*0.005*self.u_bed*self.u_bed;
        self.excursion_bed=self.phase_maximum("bed_orbital")?.0*self.t_target/(2.0*PI);
        for value in [self.l,self.k,self.c,self.power,self.sxx,self.energy_density,self.cg,self.impulse,
            self.potential_energy,self.kinetic_energy,self.mean_square_bed_velocity,
            self.momentum_flux,self.u_bed,self.u_surf,self.acc_max,self.acc_bed_max,self.w_max] {
            if !value.is_finite(){return Err("Non-finite dimensional output".into());}
        }
        let (r,ss)=self.eqns_for_state(&self.z)?;
        if !ss.is_finite(){return Err("Non-finite final residual".into());}
        self.max_residual=self.residual_inf_norm(&r);
        if self.max_residual>FINAL_RESIDUAL_TOLERANCE{return Err("Final residual exceeds tolerance".into());}
        Ok(())
    }
    /// Generate symmetric Gauss–Legendre nodes and weights on [-1,1].
    ///
    /// Roots of P_n are refined by Newton iteration using the polynomial
    /// three-term recurrence. Reflection symmetry fills the paired nodes;
    /// weights are 2/[(1-x_i^2) (P_n'(x_i))^2]. Failed root refinement
    /// rejects the numerical verification run.
    fn gauss_legendre(count: usize)->ResultT<(Vec<f64>,Vec<f64>)>{
        let mut nodes=vec![0.0;count];let mut weights=vec![0.0;count];
        for i in 0..(count+1)/2 {
            let mut root=(PI*(i as f64+0.75)/(count as f64+0.5)).cos();
            let mut derivative=0.0;
            let mut ready=false;
            for _ in 0..100 {
                let(mut p0,mut p1)=(1.0,root);
                for j in 2..=count {
                    let next=((2*j-1) as f64*root*p1-(j-1) as f64*p0)/j as f64;
                    p0=p1;p1=next;
                }
                derivative=count as f64*(root*p1-p0)/(root*root-1.0);
                let delta=p1/derivative;
                if delta.abs()<=4.0*f64::EPSILON {ready=true;break;}
                root-=delta;
            }
            if !ready {return Err("Gauss-Legendre node iteration did not converge".into());}
            let weight=2.0/((1.0-root*root)*derivative*derivative);
            nodes[i]=-root;nodes[count-1-i]=root;
            weights[i]=weight;weights[count-1-i]=weight;
        }
        Ok((nodes,weights))
    }
    /// Independently integrate the defining hydrodynamic quantities.
    ///
    /// Outer integration is a uniform 0..2*pi phase average. For each phase,
    /// map vertical Gauss–Legendre nodes from [-1,1] to [-d,eta(X)] and sum
    /// kinetic energy, impulse, radiation stress and energy transport using
    /// reconstructed dimensional velocities and Bernoulli pressure. Use two
    /// levels of phase/vertical resolution to detect inadequate quadrature.
    fn defining_integrals(&self,phase_count:usize,vertical_count:usize)->ResultT<BTreeMap<String,f64>>{
        let (nodes,weights)=Self::gauss_legendre(vertical_count)?;
        let(mut impulse,mut kinetic,mut potential,mut stress,mut power)=(0.0,0.0,0.0,0.0,0.0);
        let mut cos=vec![0.0;self.n+1];let mut sin=vec![0.0;self.n+1];
        // Depth ratios and spectral amplitudes are fixed throughout the
        // two-dimensional integration; only the local vertical coordinate
        // and phase change between quadrature nodes.
        let mut amplitudes=vec![0.0;self.n+1];
        // Cache the phase-invariant harmonic amplitudes to avoid repeatedly
        // calculating j*B_j*sqrt(g/k) in the nested integration loops.
        let velocity_scale=(self.g/self.k).sqrt();
        for j in 1..=self.n {
            amplitudes[j]=j as f64*self.b[j]*velocity_scale;
        }
        for i in 0..phase_count {
            let phase=2.0*PI*i as f64/phase_count as f64;
            let eta=self.surface_keta(phase)/self.k;
            let thickness=(self.d+eta)/2.0;
            if thickness<=0.0 || !thickness.is_finite(){return Err("The interpolated surface intersects the bed or is non-finite".into());}
            potential+=0.5*RHO*self.g*eta*eta/phase_count as f64;
            for j in 1..=self.n {
                cos[j]=(j as f64*phase).cos();sin[j]=(j as f64*phase).sin();
            }
            for m in 0..vertical_count {
                let y=-self.d+thickness*(nodes[m]+1.0);
                let(mut u,mut v)=(self.eulerian_current,0.0);
                for j in 1..=self.n {
                    let jf=j as f64;
                    let a=jf*self.k*y;let b=jf*self.z[1];
                    let positive=a.exp();let reflected=(-a-2.0*b).exp();
                    let amplitude=amplitudes[j];
                    u+=(positive+reflected)/self.depth_denominator[j]*amplitude*cos[j];
                    v+=(positive-reflected)/self.depth_denominator[j]*amplitude*sin[j];
                }
                // Apply reduced Bernoulli pressure in the travelling frame.
                // Integrands for kinetic energy, momentum flux and power use
                // the same physical velocity field but different invariants.
                let pressure=RHO*(self.bernoulli_r-self.g*y-0.5*((u-self.c)*(u-self.c)+v*v));
                let weight=weights[m]*thickness/phase_count as f64;
                impulse+=weight*RHO*u;
                kinetic+=weight*0.5*RHO*(u*u+v*v);
                stress+=weight*(pressure+RHO*u*u);
                power+=weight*(pressure+RHO*self.g*y+0.5*RHO*(u*u+v*v))*u;
            }
        }
        Ok(BTreeMap::from([
            ("Impulse".into(),impulse),("KineticEnergy".into(),kinetic),
            ("PotentialEnergy".into(),potential),("Sxx".into(),stress-0.5*RHO*self.g*self.d*self.d),
            ("Power".into(),power)]))
    }
    /// Enforce independent numerical/physical consistency requirements.
    ///
    /// Evaluate the boundary conditions on a dense phase grid, compare two
    /// spectral integral identities, compare coarse/fine defining integrals,
    /// and call the interval finite-system checker. Every discrepancy uses
    /// an absolute-plus-relative tolerance budget. Failed checks are listed
    /// in `last_error`, the solution status is set to resolution_failure, and
    /// the case cannot be reported as accepted.
    fn verify_accuracy(&mut self)->ResultT<()> {
        self.accuracy_verified=false;
        self.finite_system_verified=false;
        self.finite_system_certificate=Certificate::default();
        self.accuracy_metrics.clear();
        let mut failures:Vec<String>=Vec::new();
        // Continuous Fourier surface conditions are tested at many phases,
        // rather than only at the N+1 nodes used to solve the equations.
        let phase_count=2048_usize.max(32*self.n);
        let mut boundary=0.0_f64;
        for i in 0..phase_count {
            let phase=2.0*PI*i as f64/phase_count as f64;
            let eta=self.surface_keta(phase);
            let(mut psi,mut u,mut v)=(0.0,0.0,0.0);
            for j in 1..=self.n {
                let jf=j as f64;let a=jf*eta;let b=jf*self.z[1];
                let positive=a.exp();let reflected=(-a-2.0*b).exp();
                let sh=(positive-reflected)/self.depth_denominator[j];
                let ch=(positive+reflected)/self.depth_denominator[j];
                let co=(jf*phase).cos();let si=(jf*phase).sin();
                psi+=sh*co*self.b[j];
                u+=ch*co*jf*self.b[j];
                v+=sh*si*jf*self.b[j];
            }
            let streamline=psi-self.z[8]-self.z[7]*eta;
            let bernoulli=0.5*((u-self.z[7]).powi(2)+v*v)+eta-self.z[9];
            if !streamline.is_finite() || !bernoulli.is_finite(){return Err("Non-finite interpolated-surface residual".into());}
            boundary=boundary.max(streamline.abs().max(bernoulli.abs()));
        }
        self.accuracy_metrics.insert("surface_residual".into(),boundary);
        if boundary>1e-6 {failures.push("interpolated-surface residual exceeds 1e-6".into());}
        // Independent spectral/bed identities for wave power and radiation
        // stress; use the already computed dimensional energy and impulse.
        let beta=self.mean_square_bed_absolute_velocity;
        let spectral=BTreeMap::from([
            ("Power".to_string(),self.c*(3.0*self.kinetic_energy-2.0*self.potential_energy)
                +0.5*beta*(self.impulse+RHO*self.c*self.d)-2.0*self.c*self.eulerian_current*self.impulse),
            ("Sxx".to_string(),4.0*self.kinetic_energy-3.0*self.potential_energy+RHO*beta*self.d
                -2.0*self.eulerian_current*self.impulse)]);
        let values=BTreeMap::from([
            ("Power".to_string(),self.power),("Sxx".to_string(),self.sxx),
            ("Impulse".to_string(),self.impulse),("KineticEnergy".to_string(),self.kinetic_energy),
            ("PotentialEnergy".to_string(),self.potential_energy)]);
        for (key,reference) in &spectral {
            let error=(values[key]-*reference).abs();
            let budget=1e-7+1e-8*reference.abs();
            self.accuracy_metrics.insert(format!("{key}_spectral_ratio"),error/budget);
            if !error.is_finite() || error>budget {
                failures.push(format!("{key} spectral-bed comparison exceeds its budget"));
            }
        }
        // Independently integrate at two resolutions and compare both their
        // mutual agreement and their agreement with spectral expressions.
        let coarse=self.defining_integrals(256.max(8*self.n),64)?;
        let fine=self.defining_integrals(512.max(16*self.n),96)?;
        for (key,reference) in &fine {
            let budget=1e-7+1e-8*reference.abs();
            let difference=(*reference-coarse[key]).abs();
            let error=(values[key]-*reference).abs();
            self.accuracy_metrics.insert(format!("{key}_quadrature_ratio"),difference/budget);
            self.accuracy_metrics.insert(format!("{key}_integral_ratio"),error/budget);
            if !difference.is_finite() || !error.is_finite() || error.max(difference)>budget {
                failures.push(format!("{key} defining-integral comparison exceeds its budget"));
            }
        }
        // Output formulas and all finite nonlinear residuals are enclosed
        // separately by interval arithmetic; failure rejects this solution.
        if !self.verify_finite_system()? {
            failures.push("finite-system interval verification exceeds its budget".into());
        }
        self.accuracy_metrics.insert("finite_residual_ratio".into(),
            self.finite_system_certificate.residual_upper/self.finite_system_certificate.residual_tolerance);
        let mut output_ratio=0.0_f64;
        for bound in self.finite_system_certificate.outputs.values(){
            output_ratio=output_ratio.max(bound.error_upper/bound.budget_lower);
        }
        self.accuracy_metrics.insert("finite_output_ratio".into(),output_ratio);
        if !failures.is_empty() {
            self.converged=false;
            self.solution_status="resolution_failure".into();
            self.last_error=format!("The collocation system converged, but the available Fourier resolution does not satisfy the output accuracy budgets: {}.",failures.join("; "));
            return Ok(());
        }
        self.accuracy_verified=true;
        self.solution_status="verified_numerically".into();
        Ok(())
    }
    /// Evaluate the complete accepted state using interval-valued expressions.
    ///
    /// Constructs all 2N+10 residual intervals and nine selected derived
    /// output intervals (L,c,impulse,kinetic/potential/total energy,Sxx,power,
    /// mean-square orbital bed velocity). Compare a conservative upper bound
    /// on residual magnitude with 1e-10 and each reported output with its
    /// own absolute-plus-relative error budget. The resulting certificate
    /// explicitly records that exact root existence was NOT proven.
    fn verify_finite_system(&mut self)->ResultT<bool>{
        self.finite_system_verified=false;
        self.finite_system_certificate=Certificate::default();
        // Check binary64 radix, mantissa precision and gradual underflow.
        // `black_box` avoids compile-time folding of the subnormal probe.
        let smallest=std::hint::black_box(f64::from_bits(1));
        if f64::RADIX!=2 || f64::MANTISSA_DIGITS!=53 ||
            std::hint::black_box(smallest*2.0)!=f64::from_bits(2) {
            return Err("Interval verification requires IEEE binary64 with gradual underflow".into());
        }
        let i=I::point;
        let state=self.z.iter().map(|&v|i(v)).collect::<Vec<_>>();
        let h=i(self.h_target);let t=i(self.t_target);
        let depth=i(self.d);let current=i(self.uc);
        let gravity=i(self.g);let density=i(RHO);
        let n=self.n;
        // Interval form of the very same Fenton residual system used in
        // Newton, but stored zero-based here: f[0] corresponds to F1.
        let mut f=vec![i(0.0);2*n+10];
        f[0]=state[2]-state[1]*h/depth;
        f[1]=state[2]-h/(gravity*t.square())*state[3].square();
        f[2]=state[4]*state[3]-i(2.0)*verified::pi();
        f[3]=state[5]+state[7]-state[4];
        f[4]=state[1]*(state[6]+state[7]-state[4])-state[8];
        f[5]=state[self.current_criterion+4]-current/(gravity*depth).sqrt()*state[1].sqrt();
        f[6]=state[10]+state[n+10];
        for j in 11..n+10 {f[6]=f[6]+i(2.0)*state[j];}
        f[7]=state[10]-state[n+10]-state[2];
        // Certified trigonometric lookup for each harmonic at X=m*pi/N.
        // The normalized hyperbolic basis is evaluated through interval exp.
        let table=verified::collocation(n);
        let mut den=Vec::with_capacity(n);
        for j in 1..=n {den.push(i(1.0)+(-i((2*j) as f64)*state[1]).exp());}
        for m in 0..=n {
            let ordinate=state[m+10];
            let(mut psi,mut u,mut v)=(i(0.0),i(0.0),i(0.0));
            for j in 1..=n {
                let jf=i(j as f64);
                let ep=(jf*ordinate).exp();
                let em=(-jf*ordinate-i((2*j) as f64)*state[1]).exp();
                let sh=(ep-em)/den[j-1];let ch=(ep+em)/den[j-1];
                let trig=table[m*j%(2*n)];
                let coefficient=state[n+10+j];
                psi=psi+coefficient*sh*trig.1;
                u=u+jf*coefficient*ch*trig.1;
                v=v+jf*coefficient*sh*trig.0;
            }
            f[8+m]=psi-state[8]-state[7]*ordinate;
            f[n+9+m]=((u-state[7]).square()+v.square())/i(2.0)+ordinate-state[9];
        }
        // Form the nine selected output formulas independently from the
        // accepted state and surface Fourier coefficients, with correct SI
        // scaling for velocity, energy, impulse, stress and wave power.
        let wave_number=state[1]/depth;
        let speed_scale=(gravity/wave_number).sqrt();
        let energy_scale=density*gravity/wave_number.square();
        let pulse=state[8]+state[1]*state[5];
        let flux_scaled=state[7]*state[1]-state[8];
        let kinetic=(state[4]*pulse-state[5]*flux_scaled)/i(2.0);
        let mut potential=i(0.0);
        for j in 1..n {potential=potential+i(self.y[j]).square()/i(4.0);}
        potential=potential+i(self.y[n]).square()/i(16.0);
        let beta=i(2.0)*state[9]-state[4].square()+i(2.0)*state[4]*state[5];
        let stress=i(4.0)*kinetic-i(3.0)*potential+beta*state[1]-i(2.0)*state[5]*pulse;
        let power=state[4]*(i(3.0)*kinetic-i(2.0)*potential)
            +beta*(pulse+state[4]*state[1])/i(2.0)-i(2.0)*state[4]*state[5]*pulse;
        let mut bed_sum=i(0.0);
        for j in 1..=n {
            let inverse_cosh=i(2.0)*(-i(j as f64)*state[1]).exp()/den[j-1];
            bed_sum=bed_sum+(i(j as f64)*state[n+10+j]*inverse_cosh).square()/i(2.0);
        }
        let outputs=BTreeMap::from([
            ("L".to_string(),i(2.0)*verified::pi()/wave_number),
            ("c".to_string(),state[4]*speed_scale),
            ("Impulse".to_string(),density*speed_scale/wave_number*pulse),
            ("KineticEnergy".to_string(),energy_scale*kinetic),
            ("PotentialEnergy".to_string(),energy_scale*potential),
            ("EnergyDensity".to_string(),energy_scale*(kinetic+potential)),
            ("Sxx".to_string(),energy_scale*stress),
            ("Power".to_string(),energy_scale*speed_scale*power),
            ("MeanSquareBedVelocity".to_string(),gravity/wave_number*bed_sum)
        ]);
        let reported=BTreeMap::from([
            ("L".to_string(),self.l),("c".to_string(),self.c),
            ("Impulse".to_string(),self.impulse),
            ("KineticEnergy".to_string(),self.kinetic_energy),
            ("PotentialEnergy".to_string(),self.potential_energy),
            ("EnergyDensity".to_string(),self.energy_density),
            ("Sxx".to_string(),self.sxx),("Power".to_string(),self.power),
            ("MeanSquareBedVelocity".to_string(),self.mean_square_bed_velocity)
        ]);
        // Save equation bounds and all output comparisons. The conservative
        // residual limit is rounded toward minus infinity so that accepting
        // the interval comparison does not loosen the nominal tolerance.
        let cert=&mut self.finite_system_certificate;
        cert.state=self.z.clone();
        cert.surface_coefficients=self.y[..=n].to_vec();
        cert.order=n;cert.equation_count=f.len();cert.residual_intervals=f;
        cert.residual_tolerance=I::down(1e-10);
        for residual in &cert.residual_intervals {cert.residual_upper=cert.residual_upper.max(residual.magnitude());}
        let mut valid=cert.residual_upper<=cert.residual_tolerance;
        for (key,interval) in outputs {
            let min_abs=if interval.lo<=0.0 && interval.hi>=0.0 {0.0}
                else{interval.lo.abs().min(interval.hi.abs())};
            let check=OutputBound{
                enclosure:interval,reported:reported[&key],
                error_upper:(i(reported[&key])-interval).magnitude(),
                budget_lower:I::down(I::down(1e-7)+I::down(I::down(1e-8)*min_abs))
            };
            // Each accepted output must fit the budget; a passed global
            // residual alone is insufficient for acceptance of report values.
            valid=valid && check.error_upper<=check.budget_lower;
            cert.outputs.insert(key,check);
        }
        cert.passed=valid;self.finite_system_verified=valid;
        Ok(valid)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    // Reject a negative physical wave height without allowing later
    // calculations to consume invalid nondimensional input values.
    #[test]
    fn input_validation() {
        let mut s=FentonStreamFunction::new(-1.0,9.0,5.0,0.0,12);
        s.solve();
        assert!(!s.converged);
        assert!(s.last_error.contains("Invalid inputs"));
    }
    // Regression for a nearly rank-sensitive matrix with equal-norm columns.
    // Confirms finite SVD correction and a small linear-system residual.
    #[test]
    fn jacobi_equal_norm_columns() {
        // Equal column norms must still admit a nonzero Jacobi rotation.
        let matrix=[1.0,1.0,0.0,0.0,1.0,1.0,1.0,0.0,0.0];
        let rhs=[1.0,1.0,1.0];
        let (solution,condition)=crate::linalg::svd_solve(&matrix,&rhs,3).unwrap();
        assert!(solution.iter().all(|x|x.is_finite()));
        assert!(condition.is_finite() && condition>=1.0);
        for row in 0..3 {
            let actual=(0..3).map(|col|matrix[row*3+col]*solution[col]).sum::<f64>();
            assert!((actual-rhs[row]).abs()<1e-12,"Jacobi SVD residual is excessive");
        }
    }
    // Exercise the approximate initialization independently of full Newton
    // continuation, validating state dimensions and finite residuals.
    #[test]
    fn initial_values_are_finite() {
        let mut s=FentonStreamFunction::new(0.5,8.0,12.0,0.0,12);
        s.height_step=s.height/4.0;s.hoverd=s.max_h/4.0;
        s.init_linear();
        assert!(s.z[1].is_finite() && s.z[1]>0.0);
        assert!(s.z[3].is_finite());
        let (r,ss)=s.eqns_for_state(&s.z).unwrap();
        assert_eq!(r.len(),2*12+11);
        assert!(ss.is_finite());
    }
    // Expensive nonlinear regression data: wavelength and wave power are
    // compared to rounded five-decimal reference cases, including currents.
    // Kept ignored by default because it exercises full verification work.
    #[test]
    #[ignore = "Computational reference test; run with --ignored"]
    fn reference_wavelengths() {
        // Independent double-precision reference solutions, rounded to five decimals.
        for (h,t,d,uc,wavelength,power_kw) in [
            (0.5,8.0,12.0,0.0,75.88257,2.31615),
            (1.0,9.0,15.0,0.5,101.38371,13.08299),
            (1.0,5.0,10.0,-0.5,32.59377,3.54351),
            (3.0,9.0,5.0,1.0,78.82722,73.50878),
        ] {
            let mut s=FentonStreamFunction::new(h,t,d,uc,50);
            s.solve();
            assert!(s.converged,"H={h}, T={t}, d={d}, Uc={uc}: {}",s.last_error);
            assert!(s.accuracy_verified && s.finite_system_verified);
            assert!((s.l-wavelength).abs()<0.00006,
                "Reference wavelength mismatch H={h}, T={t}, d={d}, Uc={uc}: {} vs {wavelength}",s.l);
            assert!((s.power/1000.0-power_kw).abs()<0.00010,
                "Reference power mismatch H={h}, T={t}, d={d}, Uc={uc}: {} vs {power_kw}",s.power/1000.0);
        }
    }
    // Additional full nonlinear integration check of solution convergence,
    // final ||F||_inf and k*L=2*pi. Enable explicitly when required.
    #[test]
    #[ignore = "Full nonlinear integration test; run explicitly with cargo test -- --ignored"]
    fn nonlinear_small_wave() {
        let mut s=FentonStreamFunction::new(0.5,8.0,12.0,0.0,24);
        s.solve();
        assert!(s.converged,"{}",s.last_error);
        assert!(s.max_residual<=1e-10);
        assert!((s.k*s.l-2.0*PI).abs()<1e-10);
    }
}
}

// ===========================================================================
// MODULE: report
// ===========================================================================
// ============================================================================
// MODULE report — ENGINEERING TABLES, UNITS, STATUS AND SOLUTION.RES FORMAT
// ============================================================================
// One shared report generator handles a zero-current comparison and (if Uc is
// nonzero) an additional current-bearing wave. Only numerically accepted cases
// are reported. The fixed-column format is designed for a monospaced Win32
// multiline EDIT control and a UTF-8 plain-text output.txt file.
mod report {
//! Fixed-width 107-character engineering report, hydrodynamics tables and glossary.
//! The 2026-11-02 Python `output.txt` is the layout reference: labels, row
//! ordering, numerical precision, spacing, wrapping and glossary wording.
//! Numerical values continue to be computed from the Rust solver, not copied
//! from the single supplied reference wave case.
//!
//! Reported quantities have SI engineering units. Missing/failed wave cases
//! are not filled with plausible numerical guesses: their columns are blanked
//! and explanatory diagnostic messages are retained. The output uses Unicode
//! scientific labels but calculates widths by Unicode scalar values.
use crate::solver::{FentonStreamFunction as Solver,G_STD};
use std::fmt::Write as _;
use std::thread;

// Exact total report line width including the two outer border characters.
// A monospaced font is essential; some non-ASCII glyphs may occupy different
// terminal display widths even when their scalar-value count is one.
const WIDTH:usize=107;
#[derive(Clone)]
// A table cell is absent, numeric, or textual. Explicit Cell::None prevents
// reporting unverified values for a failed current or zero-current case.
enum Cell { None, Num(f64), Text(String) }
impl Cell {
    fn n(v:f64)->Self{Self::Num(v)}
    fn s(v:impl Into<String>)->Self{Self::Text(v.into())}
}
#[derive(Clone)]
// Rows can be section dividers (centered text) or regular four/five-column
// data. The section divider does not participate in numeric alignment.
struct Row {section:Option<String>,cells:Vec<Cell>}
impl Row {
    fn section(s:&str)->Self {Self{section:Some(s.into()),cells:Vec::new()}}
    fn cells(cells:Vec<Cell>)->Self{Self{section:None,cells}}
}
// Unicode-aware padding functions avoid cutting through UTF-8 code units.
// The reference report widths are measured in characters rather than bytes.
fn char_len(s:&str)->usize{s.chars().count()}
fn trunc(s:&str,w:usize)->String{s.chars().take(w).collect()}
fn pad_left(s:&str,w:usize)->String{let t=trunc(s,w);format!("{}{}"," ".repeat(w.saturating_sub(char_len(&t))),t)}
fn pad_right(s:&str,w:usize)->String{let t=trunc(s,w);format!("{}{}",t," ".repeat(w.saturating_sub(char_len(&t))))}
fn pad_center(s:&str,w:usize,fill:char)->String {
    let t=trunc(s,w);let left=w.saturating_sub(char_len(&t))/2;
    let right=w.saturating_sub(char_len(&t)+left);
    format!("{}{}{}",fill.to_string().repeat(left),t,fill.to_string().repeat(right))
}
/// Render an input with a decimal marker even for integral-valued f64.
/// This produces a stable, human-readable representation of GUI inputs.
fn py_str_float(v:f64)->String {
    if !v.is_finite(){return "nan".into();}
    let mut s=v.to_string();
    if !s.contains('.')&&!s.contains('e')&&!s.contains('E'){s.push_str(".0");}
    s
}
/// Fit finite numbers into a fixed-width table cell.
/// Try 5 down to 0 fixed decimals; if these do not fit, use scientific
/// notation with decreasing precision. This affects presentation only:
/// underlying solver values remain f64 and are not rounded before checks.
fn fmt_float(v:f64,w:usize)->String{
    if !v.is_finite(){return pad_left("nan",w);}
    for dec in [5,4,3,2,1,0] {
        let s=format!("{v:.dec$}");
        if char_len(&s)<=w{return pad_left(&s,w);}
    }
    for dec in [6,5,4,3] {
        let s=format!("{v:.dec$e}");
        if char_len(&s)<=w{return pad_left(&s,w);}
    }
    pad_left(&trunc(&format!("{v:.2e}"),w),w)
}
/// Render one table cell using requested width/alignment.
/// Text cells strip line breaks and surrounding whitespace; numeric cells
/// are right-aligned; missing cells display a dash in the relevant column.
fn fmt_cell(c:&Cell,w:usize,align:&str)->String {
    match c {
        Cell::Num(v)=>fmt_float(*v,w),
        Cell::Text(s)=>{
            let s=s.replace('\n'," ");let s=s.trim();
            match align {"left"=>pad_right(s,w),"center"=>pad_center(s,w,' '),_=>pad_left(s,w)}
        },
        Cell::None=>pad_left("-",w)
    }
}
// Reusable table-border and title helpers. All separator/title lines have
// total length WIDTH, including the '+' and '|' boundary characters.
fn hline(out:&mut String,ch:char){let _=writeln!(out,"+{}+",ch.to_string().repeat(WIDTH-2));}
fn box_title(out:&mut String,title:&str){hline(out,'-');let _=writeln!(out,"|{}|",pad_center(title,WIDTH-2,' '));hline(out,'-');}
fn box_text(out:&mut String,text:&str){
    let mut s=text.replace('\n'," ").trim_end().to_string();
    if char_len(&s)>WIDTH-2 {s=format!("{}...",trunc(&s,WIDTH-5));}
    let _=writeln!(out,"|{}|",pad_right(&s,WIDTH-2));
}
fn table_sep(out:&mut String,col_w:&[usize]){
    out.push('|');
    for (i,w) in col_w.iter().enumerate(){
        out.push_str(&"-".repeat(w+2));
        out.push(if i+1<col_w.len(){'+'}else{'|'});
    }
    out.push('\n');
}
/// Wrap cell text on whitespace, preserving whole words where practical.
/// Strings longer than the column width are ultimately truncated in fmt_cell
/// to prevent a malformed report table.
fn wrap_text(text:&str,width:usize)->Vec<String>{
    let s=text.replace('\n'," ").trim_end().to_string();
    if char_len(&s)<=width{return vec![s];}
    let mut lines=Vec::new();let mut line=String::new();
    for word in s.split_whitespace(){
        if line.is_empty(){line=word.into();}
        else if char_len(&line)+1+char_len(word)<=width {line.push(' ');line.push_str(word);}
        else {lines.push(line);line=word.into();}
        if !line.is_empty()&&char_len(&line)>width {lines.push(line);line=String::new();}
    }
    if !line.is_empty(){lines.push(line);}
    if lines.is_empty(){lines.push(String::new());}
    lines
}
/// Write one complete formatted report table, including title, headings,
/// separator, section rows, wrapped cell rows and final border.
/// `headers`, `widths`, and `aligns` are parallel arrays; content producers
/// must provide matching column counts. Section rows span the table width.
fn print_table(out:&mut String,title:&str,headers:&[&str],widths:&[usize],aligns:&[&str],rows:&[Row]){
    box_title(out,title);
    out.push('|');
    for i in 0..headers.len() {let _=write!(out," {} |",fmt_cell(&Cell::s(headers[i]),widths[i],aligns[i]));}
    out.push('\n');table_sep(out,widths);
    for row in rows {
        if let Some(title)=&row.section {
            let _=writeln!(out,"|{}|",pad_center(&format!(" {title} "),WIDTH-2,'-'));
            continue;
        }
        // Calculate the maximum number of physical lines required by this
        // logical row; later lines leave numeric cells empty to avoid repeats.
        let mut wrapped=Vec::with_capacity(row.cells.len());let mut nlines=1;
        for (i,c) in row.cells.iter().enumerate(){
            if let Cell::Text(s)=c {
                let lines=wrap_text(s,widths[i]);nlines=nlines.max(lines.len());wrapped.push(lines);
            }else{wrapped.push(vec![String::new()]);}
        }
        for li in 0..nlines {
            out.push('|');
            for (ci,c) in row.cells.iter().enumerate(){
                let cell=match c {
                    Cell::Text(_)=>Cell::s(wrapped[ci].get(li).map(|v|v.as_str()).unwrap_or("")),
                    _ if li>0=>Cell::s(""),
                    _=>c.clone(),
                };
                let _=write!(out," {} |",fmt_cell(&cell,widths[ci],aligns[ci]));
            }
            out.push('\n');
        }
    }
    hline(out,'-');out.push('\n');
}
// Report diagnostics use scientific notation for small residual values;
// values are formatted only after acceptance calculations are complete.
/// Translate solver acceptance flags and finite-value checks into a diagnostic.
/// A nonconverged or nonfinite wave is never silently displayed as valid.
fn issue(s:&Solver,label:&str)->Option<String>{
    if !s.converged{return Some(format!("[{label}] {}",if s.last_error.is_empty(){"Did not converge."}else{&s.last_error}));}
    if [s.l,s.k,s.c].iter().any(|v|!v.is_finite()){return Some(format!("[{label}] Non-finite result."));}
    None
}

/// Solve the input regular wave and assemble the complete engineering report.
///
/// A zero-current case is always evaluated. If Uc != 0, the current-bearing
/// case runs concurrently in a second scoped thread. Each solver object has
/// a separate Jacobian, Fourier state, residual and certificate, so there is
/// no shared mutable numerical state between the two computations.
/// The returned String is later saved as UTF-8 and displayed in the GUI.
pub fn generate_output(h:f64,t:f64,d:f64,uc:f64)->String {
    let has_current=uc!=0.0;
    let mut solver0=Solver::new(h,t,d,0.0,50);
    let mut solver_c=Solver::new(h,t,d,uc,50);
    // The independent no-current and with-current Fenton problems can be
    // calculated in parallel. Scoped threads borrow local solver objects
    // without forcing global storage or shared-mutable solver state.
    if has_current {
        thread::scope(|scope|{
            let j=scope.spawn(||solver0.solve());
            solver_c.solve();
            let _=j.join();
        });
    }else{solver0.solve();}
    let mut issues=Vec::new();
    if let Some(e)=issue(&solver0,"No current"){issues.push(e);}
    if has_current {if let Some(e)=issue(&solver_c,"With current"){issues.push(e);}}
    let valid0=issue(&solver0,"No current").is_none();
    let validc=has_current&&issue(&solver_c,"With current").is_none();
    let mut out=String::new();
    // If neither scenario meets all solver/verification requirements, emit
    // only an explanatory error header; do not manufacture result tables.
    if !valid0 && !validc {
        box_title(&mut out,"INVALID, OUT-OF-DOMAIN, OR NON-CONVERGED WAVE");
        for msg in &issues {let _=writeln!(out,"{msg}");}
        hline(&mut out,'-');return out;
    }
    // Show numerical results only for successfully verified cases.
    // Fallback references simplify shared report-scaling calculations, but
    // failed columns are still explicitly masked before the table is emitted.
    let s0=if valid0{&solver0}else{&solver_c};
    let sc=if validc{&solver_c}else{&solver0};
    box_title(&mut out,"NONLINEAR WAVE HYDRODYNAMICS SOLVER (FENTON)");
    box_text(&mut out,&format!("Wave height (H)             : {} m",py_str_float(h)));
    box_text(&mut out,&format!("Wave period (τ)             : {} s",py_str_float(t)));
    box_text(&mut out,&format!("Water depth (d)             : {} m",py_str_float(d)));
    box_text(&mut out,&format!("Eulerian current ū₁         : {} m/s (positive with wave propagation)",py_str_float(uc)));
    hline(&mut out,'-');
    box_text(&mut out,if issues.is_empty(){"Status: Finite-system and output accuracy checks passed."}
                      else{"Status: Results shown for accepted cases only."});
    for (valid,s,label) in [(valid0,s0,"No current"),(validc,sc,"With current")] {
        if !valid{continue;}
        box_text(&mut out,&format!("{label}: N={}; max|F|={:.3e}",s.n,s.max_residual));
        if s.finite_system_verified {
            box_text(&mut out,&format!("Finite-system interval residual bound: {:.3e}; nine output bounds checked.",
                s.finite_system_certificate.residual_upper));
        }
        if !s.numerical_note.is_empty(){box_text(&mut out,&s.numerical_note);}
    }
    hline(&mut out,'-');
    for e in &issues {let _=writeln!(out,"{e}");}
    out.push('\n');
    // Physical conversion scales sqrt(gd) [m/s] and sqrt(g/d) [1/s] used by
    // the depth-normalized dispersion and mean-flow table.
    let sqrt_gd=(G_STD*d).sqrt();
    let sqrt_g_over_d=(G_STD/d).sqrt();
    let wc_num=|v:f64|if has_current{Cell::n(v)}else{Cell::s("-")};
    let wc_str=|v:&str|if has_current{Cell::s(v)}else{Cell::s("-")};
    let c0=s0.l/s0.t_target;let cc=sc.l/sc.t_target;
    let mut rows=Vec::new();
    // Local row-construction macros guarantee consistent parameter labels,
    // units and paired zero/current columns throughout this large table.
    macro_rules! pair_num {
        ($label:expr,$a:expr,$b:expr,$unit:expr) => {
            rows.push(Row::cells(vec![Cell::s($label),Cell::n($a),wc_num($b),Cell::s($unit)]));
        }
    }
    macro_rules! pair_str {
        ($label:expr,$a:expr,$b:expr,$unit:expr) => {
            rows.push(Row::cells(vec![Cell::s($label),Cell::s($a),wc_str($b),Cell::s($unit)]));
        }
    }
    // Engineering summary is deliberately partitioned into inputs, geometry,
    // mean flows, flux/energy, bed/free-surface kinematics and breaking screens.
    rows.push(Row::section("INPUTS & REFERENCE SCALES"));
    pair_num!("Water depth (d)",s0.d,sc.d,"m");
    pair_num!("Wave height (H)",s0.h_target,sc.h_target,"m");
    pair_num!("Wave period (τ)",s0.t_target,sc.t_target,"s");
    pair_num!("H/d",s0.h_target/s0.d,sc.h_target/sc.d,"-");
    pair_num!("τ√(g/d)",s0.t_target*sqrt_g_over_d,sc.t_target*sqrt_g_over_d,"-");
    rows.push(Row::section("DISPERSION & PHASE (GEOMETRY)"));
    pair_num!("Wavelength (L)",s0.l,sc.l,"m");
    pair_num!("Wave number (k)",s0.k,sc.k,"rad/m");
    pair_num!("kd",s0.kd_dimless(),sc.kd_dimless(),"-");
    pair_num!("Angular frequency (ω)",2.0*std::f64::consts::PI/s0.t_target,2.0*std::f64::consts::PI/sc.t_target,"rad/s");
    pair_num!("Celerity / phase speed (c)",c0,cc,"m/s");
    pair_num!("c/√(gd)",c0/sqrt_gd,cc/sqrt_gd,"-");
    pair_num!("Crest elevation (ηc)",s0.eta_crest,sc.eta_crest,"m");
    pair_num!("Trough elevation (ηt)",s0.eta_trough,sc.eta_trough,"m");
    rows.push(Row::section("MEAN FLOWS (FENTON SOLUTION-FLAT)"));
    pair_num!("Eulerian current (ū₁)",s0.eulerian_current,sc.eulerian_current,"m/s");
    pair_num!("Stokes current (ū₂)",s0.stokes_current,sc.stokes_current,"m/s");
    pair_num!("Mean fluid speed (Ū)",s0.mean_fluid_speed,sc.mean_fluid_speed,"m/s");
    rows.push(Row::section("FLUXES & BERNOULLI CONSTANTS"));
    pair_num!("Wave volume flux (q)",s0.wave_volume_flux_q,sc.wave_volume_flux_q,"m²/s");
    pair_num!("Volume flux (Q)",s0.volume_flux_q,sc.volume_flux_q,"m²/s");
    pair_num!("Bernoulli constant (R)",s0.bernoulli_r_cap,sc.bernoulli_r_cap,"m²/s²");
    pair_num!("Reduced Bernoulli (r = R−g d)",s0.bernoulli_r,sc.bernoulli_r,"m²/s²");
    rows.push(Row::section("INTEGRAL QUANTITIES (PER UNIT CREST WIDTH)"));
    pair_num!("Kinetic energy (T)",s0.kinetic_energy/1000.0,sc.kinetic_energy/1000.0,"kJ/m²");
    pair_num!("Potential energy (V)",s0.potential_energy/1000.0,sc.potential_energy/1000.0,"kJ/m²");
    pair_num!("Total energy (E = T+V)",s0.energy_density/1000.0,sc.energy_density/1000.0,"kJ/m²");
    pair_num!("Momentum flux (S)",s0.momentum_flux/1000.0,sc.momentum_flux/1000.0,"kN/m");
    pair_num!("Radiation stress (Sₓₓ)",s0.sxx/1000.0,sc.sxx/1000.0,"kN/m");
    pair_num!("Impulse (I)",s0.impulse/1000.0,sc.impulse/1000.0,"10³ kg/(m·s)");
    pair_num!("Wave power (F)",s0.power/1000.0,sc.power/1000.0,"kW/m");
    pair_num!("Energy-flux velocity (F/E)ㅤ",s0.cg,sc.cg,"m/s");
    rows.push(Row::section("KINEMATICS (EXTREMES / BED ORBITAL MOTION)"));
    pair_num!("Surface max |u|",s0.u_surf,sc.u_surf,"m/s");
    pair_num!("Bed max |u|",s0.u_bed,sc.u_bed,"m/s");
    pair_num!("Surface max |aₓ|",s0.acc_max,sc.acc_max,"m/s²");
    pair_num!("Bed max |aₓ|",s0.acc_bed_max,sc.acc_bed_max,"m/s²");
    pair_num!("Velocity asymmetry |uc|/|ut|",s0.asymmetry,sc.asymmetry,"-");
    pair_num!("Mean square bed orbital vel ub²",s0.mean_square_bed_velocity,sc.mean_square_bed_velocity,"m²/s²");
    pair_num!("Bed orbital RMS velocity ub,rms",s0.mean_square_bed_velocity.max(0.0).sqrt(),sc.mean_square_bed_velocity.max(0.0).sqrt(),"m/s");
    rows.push(Row::section("NONLINEARITY / BREAKING DIAGNOSTICS"));
    pair_num!("Miche breaking limit (Hmax)",s0.breaking_limit_miche,sc.breaking_limit_miche,"m");
    pair_num!("Saturation (H/Hmax)",s0.breaking_index,sc.breaking_index,"-");
    pair_str!("Breaking status",if s0.is_breaking{"BREAKING"}else{"BELOW MICHE"},
              if sc.is_breaking{"BREAKING"}else{"BELOW MICHE"},"-");
    // The reference report has only the Miche three-row breaking screen,
    // followed by the Ursell number and the depth-regime classification.
    // Fenton limiting-wave information and spectral/residual diagnostics remain
    // available to the solver, but are deliberately not printed in this table.
    pair_num!("Ursell number (U)",s0.ursell,sc.ursell,"-");
    pair_str!("Regime (by d/L)",s0.regime.as_str(),sc.regime.as_str(),"-");
    // Suppress all values from unsuccessful current scenarios, independently
    // of fallback references used for report construction.
    for row in &mut rows {
        if row.cells.len()==4 {
            if !valid0 {row.cells[1]=Cell::None;}
            if !validc {row.cells[2]=Cell::None;}
        }
    }
    print_table(&mut out,"CALCULATED HYDRODYNAMIC PARAMETERS",
        &["PARAMETER","NO CURRENT","WITH CURRENT","UNIT"],&[42,16,16,20],
        &["left","right","right","left"],&rows);
    if valid0 {print_solution_flat(&mut out,&solver0,"SOLUTION.RES (NO CURRENT)");}
    if validc {print_solution_flat(&mut out,&solver_c,"SOLUTION.RES (WITH CURRENT)");}
    print_glossary(&mut out);
    out
}
/// Print exactly nineteen SOLUTION.RES entries in the reference order.
/// Entry 17 wraps its 38-character name over two physical rows, with the
/// numeric and other cells deliberately blank on the continuation row.
/// All scalar values use five decimals and remain f64 until formatting.
fn print_solution_flat(out:&mut String,s:&Solver,title:&str) {
    let (g,d,h,t,l)=(s.g,s.d,s.h_target,s.t_target,s.l);
    let c=l/t;let sqrt_gd=(g*d).sqrt();let sqrt_gd3=(g*d.powi(3)).sqrt();
    let f5=|v:f64|format!("{v:.5}");
    let mut rows=Vec::new();
    // Keep numbering, descriptions and five-decimal text aligned across all
    // 19 reference-formulation rows without duplicating table formatting code.
    macro_rules! flat {
        ($number:expr,$name:expr,$value:expr,$dimensionless:expr,$value_nd:expr) => {
            rows.push(Row::cells(vec![Cell::n($number as f64),Cell::s($name),
                Cell::s(f5($value)),Cell::s($dimensionless),Cell::s(f5($value_nd))]));
        }
    }
    flat!(1,"Water depth",d,"d/d = 1",1.0);
    flat!(2,"Wave length",l,"λ/d",l/d);
    flat!(3,"Wave height",h,"H/d",h/d);
    flat!(4,"Wave period",t,"τ√(g/d)",t*(g/d).sqrt());
    flat!(5,"Wave speed",c,"c/√(gd)",c/sqrt_gd);
    flat!(6,"Eulerian current",s.eulerian_current,"ū₁/√(gd)",s.eulerian_current/sqrt_gd);
    flat!(7,"Stokes current",s.stokes_current,"ū₂/√(gd)",s.stokes_current/sqrt_gd);
    flat!(8,"Mean fluid speed",s.mean_fluid_speed,"Ū/√(gd)",s.mean_fluid_speed/sqrt_gd);
    flat!(9,"Wave volume flux, q = Ū d − Q",s.wave_volume_flux_q,"q/√(gd³)",s.wave_volume_flux_q/sqrt_gd3);
    flat!(10,"Bernoulli constant, r = R − gd",s.bernoulli_r,"r/gd",s.bernoulli_r/(g*d));
    flat!(11,"Volume flux",s.volume_flux_q,"Q/√(gd³)",s.volume_flux_q/sqrt_gd3);
    flat!(12,"Bernoulli constant",s.bernoulli_r_cap,"R/gd",s.bernoulli_r_cap/(g*d));
    flat!(13,"Momentum flux",s.momentum_flux/1000.0,"S/ρgd²",s.momentum_flux_depth);
    flat!(14,"Impulse",s.impulse/1000.0,"I/(ρ√(gd³))",s.i_depth);
    flat!(15,"Kinetic energy",s.kinetic_energy/1000.0,"T/ρgd²",s.ke_depth);
    flat!(16,"Potential energy",s.potential_energy/1000.0,"V/ρgd²",s.pe_depth);
    flat!(17,"Mean square of bed velocity",s.mean_square_bed_absolute_velocity,
        "ub²/gd",s.mean_square_bed_absolute_velocity/(g*d));
    flat!(18,"Radiation stress",s.sxx/1000.0,"S_xx/ρgd²",s.sxx_depth);
    flat!(19,"Wave power",s.power/1000.0,"F/(ρg³ᐟ²d⁵ᐟ²)ㅤㅤ",s.f_depth);
    print_table(out,title,&["#","PARAMETER","value","adim param","adim value"],
        &[2,37,13,25,14],&["right","left","right","left","right"],&rows);
}
/// Append the glossary with the exact sequence, wording, scientific characters
/// and units of the supplied Python output reference. Some definitions wrap
/// to a second physical line to preserve the fixed 107-character layout.
fn print_glossary(out:&mut String) {
    // Labels, definitions, units and their order reproduce the supplied Python report.
    let terms:[(&str,&str,&str);33]=[
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
        ("F", "Wave power (energy flux).", "kW/m ; F/(ρg³ᐟ²d⁵ᐟ²)ㅤ"),
        ("Cg", "Energy-flux / energy-density ratio F/E.", "m/s"),
        ("Sₓₓ", "Radiation stress component in wave direction.", "kN/m ; Sₓₓ/(ρgd²)"),
        ("ub²", "Mean square *orbital* bed velocity: ub² = <(ub(t) − ū₁)²>. Non-negative by definition; computed by phase averaging.", "m²/s² ; /gd"),
        ("ub,rms", "Root-mean-square orbital bed velocity: ub,rms = √(ub²).", "m/s"),
        ("usurf,max", "Maximum absolute surface horizontal velocity over phase.", "m/s"),
        ("ubed,max", "Maximum absolute bed horizontal velocity over phase.", "m/s"),
        ("a_x,surf,max", "Maximum absolute surface material horizontal acceleration over phase.", "m/s²"),
        ("a_x,bed,max", "Maximum absolute bed material horizontal acceleration over phase.", "m/s²"),
        ("Asymmetry", "Velocity asymmetry indicator |uc|/|ut|.", "-"),
        ("Hmax", "Miche breaking limit used as stability diagnostic.", "m"),
        ("Ursell", "Ursell number, a shallow-water nonlinearity measure.", "-"),
        ("Regime", "Depth regime based on d/L (deep/intermediate/shallow).", "-"),
    ];
    let rows=terms.iter().map(|(sym,meaning,units)|Row::cells(vec![Cell::s(*sym),Cell::s(*meaning),Cell::s(*units)]))
        .collect::<Vec<_>>();
    print_table(out,"PARAMETER DEFINITIONS & GLOSSARY",&["TERM / SYMBOL","MEANING","UNITS / NONDIM"],
        &[14,64,19],&["left","left","left"],&rows);
}

#[cfg(test)]
mod tests {
    use super::*;
    // Formatting regression: every border/header/body line must occupy the
    // expected fixed-width 107-character record, including Unicode labels.
    #[test]
    fn reference_glossary_layout_and_terms() {
        let mut out=String::new();
        print_glossary(&mut out);
        assert!(out.contains("| H              | Wave height (crest-to-trough)."));
        assert!(out.contains("| Hmax           | Miche breaking limit used as stability diagnostic."));
        assert!(!out.contains("| Hm             |"));
        assert!(out.lines().filter(|s|s.starts_with('|')||s.starts_with('+'))
            .all(|s|s.chars().count()==WIDTH));
    }
    #[test]
    fn reference_solution_row_seventeen_wraps() {
        let mut out=String::new();
        let rows=vec![Row::cells(vec![Cell::n(17.0),
            Cell::s("Mean square of bed velocity"),
            Cell::s("1.80068"),Cell::s("ub²/gd"),Cell::s("0.03672")])];
        print_table(&mut out,"TEST",&["#","PARAMETER","value","adim param","adim value"],
            &[2,37,13,25,14],&["right","left","right","left","right"],&rows);
        assert!(out.contains("| 17 | Mean square of bed velocity"));
        assert!(out.lines().filter(|s|s.starts_with('|')||s.starts_with('+'))
            .all(|s|s.chars().count()==WIDTH));
    }
    #[test]
    fn borders_are_107_characters() {
        let mut out=String::new();
        print_table(&mut out,"TEST",&["a","b"],&[40,60],&["left","right"],
            &[Row::cells(vec![Cell::s("Value"),Cell::n(3.14)])]);
        assert!(out.lines().filter(|s|s.starts_with('|')||s.starts_with('+')).all(|s|s.chars().count()==WIDTH));
    }
}
}

// ===========================================================================
// MODULE: windows_gui
// ===========================================================================
// ============================================================================
// MODULE windows_gui — NATIVE WIN32 USER INTERFACE AND WORKER THREAD
// ============================================================================
// This module is compiled only on Windows when the optional 'headless' cfg
// is absent. It calls the Windows APIs directly without a Rust GUI framework.
// The window controls use UTF-16 text, and the report textbox is read-only
// with horizontal and vertical scroll bars. Input and report GUI state live
// on the main message-loop thread; heavy computations run on a worker thread.
#[cfg(all(windows, not(feature = "headless")))]
mod windows_gui {
//! Dependency-free native Win32 GUI (UTF-16 controls, editable inputs, read-only report).
//! The numerical code is pure Rust; Windows functions are invoked via system DLLs.
//!
//! Win32 handles are represented as pointer-sized integers in the direct ABI
//! declarations. The unsafe boundary is confined to GUI utility functions,
//! callback dispatch and system FFI; mathematical functions remain safe Rust.
//! GUI messages must be handled on the owning UI thread.
#![cfg(windows)]
use crate::report::generate_output;
use std::ffi::c_void;
use std::ptr::{null,null_mut};
use std::sync::atomic::{AtomicBool,AtomicIsize,Ordering};

// ---------------------------- WIN32 MESSAGE CODES ---------------------------
// WM_CREATE: construct controls; WM_COMMAND: handle Calculate; WM_DESTROY:
// exit message loop; WM_APP_RESULT: private message conveying completed report.
// WM_SETFONT selects the UI/report fonts; EM_SETMARGINS adjusts edit padding.
const WM_CREATE:u32=0x0001;
const WM_DESTROY:u32=0x0002;
const WM_COMMAND:u32=0x0111;
const WM_SETFONT:u32=0x0030;
const WM_APP_RESULT:u32=0x8000+1;
const EM_SETMARGINS:u32=0x00d3;
// --------------------------- WINDOW STYLE FLAGS -----------------------------
// These are ordinary Win32 style bitmasks. WS_CHILD|WS_VISIBLE creates child
// controls; ES_READONLY and ES_MULTILINE make the report a text viewer.
// WS_CAPTION|WS_SYSMENU|WS_MINIMIZEBOX prevents arbitrary frame resizing.
const EC_LEFTMARGIN:usize=1;
const EC_RIGHTMARGIN:usize=2;
const WS_CHILD:u32=0x4000_0000;
const WS_VISIBLE:u32=0x1000_0000;
const WS_BORDER:u32=0x0080_0000;
const WS_VSCROLL:u32=0x0020_0000;
const WS_HSCROLL:u32=0x0010_0000;
const WS_CAPTION:u32=0x00c0_0000;
const WS_SYSMENU:u32=0x0008_0000;
const WS_MINIMIZEBOX:u32=0x0002_0000;
const ES_MULTILINE:u32=0x0004;
const ES_READONLY:u32=0x0800;
const BS_DEFPUSHBUTTON:u32=0x0001;
const SWP_NOMOVE:u32=0x0002;
const SWP_NOZORDER:u32=0x0004;
const GWL_STYLE:i32=-16;
const GWL_EXSTYLE:i32=-20;
// Screen metrics, system color, child-control identifiers and font traits.
// Control IDs are used by WM_COMMAND to identify the Calculate button.
const SM_CXSCREEN:i32=0;
const SM_CYSCREEN:i32=1;
const COLOR_WINDOW:isize=5;
const IDC_EDIT_H:isize=101;
const IDC_EDIT_T:isize=102;
const IDC_EDIT_D:isize=103;
const IDC_EDIT_UC:isize=104;
const IDC_BTN_CALC:isize=110;
const IDC_OUTPUT:isize=111;
const MB_ICONERROR:u32=0x10;
const FW_NORMAL:i32=400;
const DEFAULT_CHARSET:u32=1;
const FF_SWISS:u32=32;
const FF_MODERN:u32=48;
const FIXED_PITCH:u32=1;

// ----------------------------- WIN32 ABI TYPES ------------------------------
// #[repr(C)] guarantees C-compatible field layout for Win32 RECT, POINT,
// MSG and WNDCLASSEXW. Using the wrong field order/width would corrupt API
// communication. These definitions target the 64-bit Windows ABI used here.
#[repr(C)]
#[derive(Copy,Clone,Default)]
struct Rect{left:i32,top:i32,right:i32,bottom:i32}
#[repr(C)]
#[derive(Copy,Clone,Default)]
struct Point{x:i32,y:i32}
#[repr(C)]
#[derive(Default)]
struct Msg{hwnd:isize,message:u32,wparam:usize,lparam:isize,time:u32,pt:Point,lprivate:u32}
#[repr(C)]
struct WndClassEx {
    cb_size:u32,style:u32,proc:Option<unsafe extern "system" fn(isize,u32,usize,isize)->isize>,
    cls_extra:i32,wnd_extra:i32,instance:isize,icon:isize,cursor:isize,background:isize,
    menu_name:*const u16,class_name:*const u16,icon_sm:isize,
}
// ----------------------------- WIN32 SYSTEM DLLS ----------------------------
// user32: window classes, edit controls, message handling, geometry, dialog.
// gdi32: fonts and GDI resource cleanup.
// kernel32: module-handle access.
// These are operating-system DLL imports; they are not third-party Rust crates.
// `extern "system"` selects Windows' native calling convention.
#[link(name="user32")]
extern "system" {
    fn RegisterClassExW(class:*const WndClassEx)->u16;
    fn CreateWindowExW(ex_style:u32,class:*const u16,name:*const u16,style:u32,
        x:i32,y:i32,w:i32,h:i32,parent:isize,menu:isize,instance:isize,param:*mut c_void)->isize;
    fn DefWindowProcW(hwnd:isize,msg:u32,wparam:usize,lparam:isize)->isize;
    fn LoadIconW(instance:isize,id:*const u16)->isize;
    fn LoadCursorW(instance:isize,id:*const u16)->isize;
    fn GetSystemMetrics(index:i32)->i32;
    fn GetWindowLongPtrW(hwnd:isize,index:i32)->isize;
    fn AdjustWindowRectEx(rect:*mut Rect,style:u32,has_menu:i32,exstyle:u32)->i32;
    fn SetWindowPos(hwnd:isize,insert_after:isize,x:i32,y:i32,w:i32,h:i32,flags:u32)->i32;
    fn SendMessageW(hwnd:isize,msg:u32,wparam:usize,lparam:isize)->isize;
    fn SetWindowTextW(hwnd:isize,text:*const u16)->i32;
    fn GetWindowTextW(hwnd:isize,text:*mut u16,cap:i32)->i32;
    fn MessageBoxW(parent:isize,text:*const u16,caption:*const u16,flags:u32)->i32;
    fn EnableWindow(hwnd:isize,enabled:i32)->i32;
    fn IsWindow(hwnd:isize)->i32;
    fn PostMessageW(hwnd:isize,msg:u32,wparam:usize,lparam:isize)->i32;
    fn PostQuitMessage(exit_code:i32);
    fn GetMessageW(msg:*mut Msg,hwnd:isize,min:u32,max:u32)->i32;
    fn TranslateMessage(msg:*const Msg)->i32;
    fn DispatchMessageW(msg:*const Msg)->isize;
    fn ShowWindow(hwnd:isize,cmd:i32)->i32;
    fn UpdateWindow(hwnd:isize)->i32;
}
#[link(name="gdi32")]
extern "system" {
    fn CreateFontW(height:i32,width:i32,escapement:i32,orientation:i32,weight:i32,
        italic:u32,underline:u32,strike_out:u32,char_set:u32,output_prec:u32,clip_prec:u32,
        quality:u32,pitch_family:u32,face:*const u16)->isize;
    fn DeleteObject(obj:isize)->i32;
}
#[link(name="kernel32")]
extern "system" {fn GetModuleHandleW(module:*const u16)->isize;}

// --------------------------- GLOBAL GUI HANDLES -----------------------------
// AtomicIsize stores HWND and HFONT handle values without a mutable static.
// Relaxed ordering suffices for ordinary handle lookup after UI creation;
// CLOSING uses Acquire/Release to signal a worker that the GUI is exiting.
// The Win32 window itself is owned by the message-loop thread.
static MAIN:AtomicIsize=AtomicIsize::new(0);
static INPUT_H:AtomicIsize=AtomicIsize::new(0);
static INPUT_T:AtomicIsize=AtomicIsize::new(0);
static INPUT_D:AtomicIsize=AtomicIsize::new(0);
static INPUT_UC:AtomicIsize=AtomicIsize::new(0);
static OUTPUT:AtomicIsize=AtomicIsize::new(0);
static BUTTON:AtomicIsize=AtomicIsize::new(0);
static UI_FONT:AtomicIsize=AtomicIsize::new(0);
static MONO_FONT:AtomicIsize=AtomicIsize::new(0);
static CLOSING:AtomicBool=AtomicBool::new(false);
/// Read a pointer-sized Win32 handle from its atomic storage.
fn get(a:&AtomicIsize)->isize{a.load(Ordering::Relaxed)}
/// Convert UTF-8 Rust text to a NUL-terminated UTF-16 Win32 string.
/// The returned allocation must remain alive for the duration of the FFI call.
fn wide(s:&str)->Vec<u16>{s.encode_utf16().chain(std::iter::once(0)).collect()}
/// Parse user numbers accepting either decimal comma or decimal point.
/// Trimming whitespace and rejecting NaN/Infinity keeps UI input validation
/// consistent with the numerical solver's finite-value requirements.
fn parse_input(text:&str)->Option<f64>{
    let t=text.trim().replace(',', ".");
    let value=t.parse::<f64>().ok()?;
    if value.is_finite(){Some(value)}else{None}
}
/// Read numeric text from a Win32 EDIT handle and parse it as `f64`.
/// Uses a fixed UTF-16 buffer of 128 code units; an empty edit returns None.
unsafe fn read_double(hwnd:isize)->Option<f64>{
    let mut buf=[0u16;128];
    let len=GetWindowTextW(hwnd,buf.as_mut_ptr(),buf.len() as i32);
    if len<=0{return None;}
    parse_input(&String::from_utf16_lossy(&buf[..len as usize]))
}
/// Display a modal Win32 error dialog with a UTF-16 message.
/// All invalid physical or numeric inputs are rejected before worker creation.
unsafe fn show_error(hwnd:isize,text:&str){
    MessageBoxW(hwnd,wide(text).as_ptr(),wide("Input Error").as_ptr(),MB_ICONERROR);
}
/// Create one child control, set its location/size and optional shared font.
/// `class` is a standard Win32 control class (STATIC, EDIT or BUTTON).
/// The HWND is returned as a pointer-sized value and retained by the UI.
unsafe fn control(hwnd:isize,class:&str,title:&str,style:u32,id:isize,x:i32,y:i32,w:i32,h:i32,font:isize)->isize{
    let window=CreateWindowExW(0,wide(class).as_ptr(),wide(title).as_ptr(),style,
         x,y,w,h,hwnd,id,GetModuleHandleW(null()),null_mut());
    if window!=0 && font!=0{SendMessageW(window,WM_SETFONT,font as usize,1);}
    window
}
/// Dimension the main frame and report area from the available screen size.
///
/// Desired client width/height are limited to 1600x950 with screen margins;
/// AdjustWindowRectEx accounts for non-client caption/border dimensions.
/// The text report control is laid out from x=300, giving input controls a
/// fixed panel width. Horizontal scroll is available for 107-column lines.
unsafe fn resize_gui(hwnd:isize) {
    let width=(GetSystemMetrics(SM_CXSCREEN)-80).min(1600).max(600);
    let height=(GetSystemMetrics(SM_CYSCREEN)-120).min(950).max(400);
    let mut rc=Rect{left:0,top:0,right:width,bottom:height};
    let style=GetWindowLongPtrW(hwnd,GWL_STYLE) as u32;
    let exstyle=GetWindowLongPtrW(hwnd,GWL_EXSTYLE) as u32;
    AdjustWindowRectEx(&mut rc,style,0,exstyle);
    SetWindowPos(hwnd,0,0,0,rc.right-rc.left,rc.bottom-rc.top,SWP_NOMOVE|SWP_NOZORDER);
    let output=get(&OUTPUT);
    if output!=0 && IsWindow(output)!=0 {
        SendMessageW(output,EM_SETMARGINS,EC_LEFTMARGIN|EC_RIGHTMARGIN,0);
        SetWindowPos(output,0,300,20,(width-320).max(100),(height-40).max(100),SWP_NOZORDER);
    }
}
/// Spawn a detached numerical worker for the selected wave parameters.
///
/// The worker computes the full report and writes output.txt as UTF-8. It
/// passes ownership of a heap-allocated String to the UI thread by posting
/// WM_APP_RESULT. On posting failure, the worker reclaims the allocation.
/// Numerical panics are caught so a failed worker supplies an error message.
/// No Win32 EDIT control is modified directly from this worker thread.
fn run_calculation(hwnd:isize,h:f64,t:f64,d:f64,uc:f64){
    std::thread::spawn(move ||{
        let result=std::panic::catch_unwind(||generate_output(h,t,d,uc));
        // Panic recovery protects the GUI from interval domain failures or
        // other unexpected panics during numerical result construction.
        let report=result.unwrap_or_else(|_|"Calculation failed: unexpected numerical panic.\n".to_string());
        let _=std::fs::write("output.txt",report.as_bytes());
        // A posted message carries a raw pointer to a Box<String>. The GUI
        // message handler reconstructs the Box exactly once and frees it.
        // If the main window closes, avoid posting when the close flag is set.
        if !CLOSING.load(Ordering::Acquire){
            let ptr=Box::into_raw(Box::new(report));
            unsafe {
                if PostMessageW(hwnd,WM_APP_RESULT,0,ptr as isize)==0 {
                    drop(Box::from_raw(ptr));
                }
            }
        }
    });
}
/// Central Win32 window procedure executed by the UI message loop.
///
/// WM_CREATE creates inputs, button and monospaced read-only output control.
/// WM_COMMAND validates inputs and launches the solver worker.
/// WM_APP_RESULT receives an owned report and re-enables the button.
/// WM_DESTROY releases fonts and requests the message loop to terminate.
/// Unhandled messages are forwarded to DefWindowProcW.
unsafe extern "system" fn wnd_proc(hwnd:isize,msg:u32,wparam:usize,lparam:isize)->isize {
    match msg {
        // -------------------------- WINDOW INITIALIZATION ---------------------
        // Register the HWND, construct Segoe UI/Consolas fonts and create four
        // initial-value EDIT controls for regular H, period, depth and current.
        WM_CREATE=>{
            MAIN.store(hwnd,Ordering::Relaxed);
            let ui=CreateFontW(19,0,0,0,FW_NORMAL,0,0,0,DEFAULT_CHARSET,0,0,0,
                FF_SWISS,wide("Segoe UI").as_ptr());
            let mono=CreateFontW(22,0,0,0,FW_NORMAL,0,0,0,DEFAULT_CHARSET,0,0,0,
                FIXED_PITCH|FF_MODERN,wide("Consolas").as_ptr());
            UI_FONT.store(ui,Ordering::Relaxed);MONO_FONT.store(mono,Ordering::Relaxed);
            let mut y=20;
            // The four input rows use consistent label positions and numeric
            // edit field widths; factory helper applies shared GUI font.
            for (label,value,id,store) in [
                ("Regular Wave H (m):","3.0",IDC_EDIT_H,&INPUT_H),
                ("Wave Period (s):","9.0",IDC_EDIT_T,&INPUT_T),
                ("Water Depth (m):","5.0",IDC_EDIT_D,&INPUT_D),
                ("Current (m/s):","1.0",IDC_EDIT_UC,&INPUT_UC)
            ] {
                control(hwnd,"STATIC",label,WS_CHILD|WS_VISIBLE,0,20,y,165,25,ui);
                let input=control(hwnd,"EDIT",value,WS_CHILD|WS_VISIBLE|WS_BORDER,id,190,y,90,25,ui);
                store.store(input,Ordering::Relaxed);
                y+=35;
            }
            y+=10;
            let button=control(hwnd,"BUTTON","CALCULATE",WS_CHILD|WS_VISIBLE|BS_DEFPUSHBUTTON,
                IDC_BTN_CALC,20,y,260,40,ui);
            BUTTON.store(button,Ordering::Relaxed);
            // Multiline read-only report: vertical/horizontal scroll bars,
            // Consolas fixed-pitch font and large text for engineering tables.
            let output=control(hwnd,"EDIT","",WS_CHILD|WS_VISIBLE|WS_BORDER|ES_MULTILINE|
                WS_VSCROLL|ES_READONLY|WS_HSCROLL,IDC_OUTPUT,300,20,800,520,mono);
            OUTPUT.store(output,Ordering::Relaxed);
            resize_gui(hwnd);
            0
        },
        // ----------------------------- CALCULATE ------------------------------
        // LOWORD(wParam) identifies the button ID. All four inputs must be
        // finite, and H,T,d strictly positive; current can be signed or zero.
        WM_COMMAND if (wparam&0xffff) as isize==IDC_BTN_CALC=>{
            let inputs=(read_double(get(&INPUT_H)),read_double(get(&INPUT_T)),
                        read_double(get(&INPUT_D)),read_double(get(&INPUT_UC)));
            match inputs {
                (Some(h),Some(t),Some(d),Some(uc))=>{
                    if h<=0.0 || t<=0.0 || d<=0.0 {
                        show_error(hwnd,"H, T, and d must be finite and positive; current must be finite.");
                    }else{
                        // Disable repeated calculation clicks until the worker
                        // sends WM_APP_RESULT. The UI thread remains responsive.
                        EnableWindow(get(&BUTTON),0);
                        SetWindowTextW(get(&BUTTON),wide("SOLVING...").as_ptr());
                        run_calculation(hwnd,h,t,d,uc);
                    }
                },
                _=>show_error(hwnd,"Invalid numeric input.\n\nAccepted formats: 3.5 or 3,5 (decimal comma supported)."),
            }
            0
        },
        // ---------------------- WORKER COMPLETION MESSAGE ---------------------
        // WM_APP_RESULT lparam owns a Box<String> posted by run_calculation.
        // Convert LF to CRLF for Windows EDIT control and release Box on scope
        // exit; re-enable the calculate button regardless of solution status.
        WM_APP_RESULT=>{
            if lparam!=0 {
                let report=Box::from_raw(lparam as *mut String);
                let display=report.replace('\n',"\r\n");
                let output=get(&OUTPUT);
                if output!=0 && IsWindow(output)!=0 {
                    SetWindowTextW(output,wide(&display).as_ptr());
                }
            }
            let button=get(&BUTTON);
            if button!=0 && IsWindow(button)!=0 {
                EnableWindow(button,1);
                SetWindowTextW(button,wide("CALCULATE").as_ptr());
            }
            0
        },
        // ---------------------------- GUI SHUTDOWN ----------------------------
        // Mark the main window as closing, free GDI font resources and request
        // GetMessageW to terminate via PostQuitMessage(0).
        WM_DESTROY=>{
            CLOSING.store(true,Ordering::Release);
            for f in [get(&UI_FONT),get(&MONO_FONT)]{
                if f!=0 {DeleteObject(f);}
            }
            PostQuitMessage(0);
            0
        },
        _=>DefWindowProcW(hwnd,msg,wparam,lparam),
    }
}
/// Register the application window class and enter the native message loop.
///
/// The class callback points to wnd_proc; Win32 provides the message pump,
/// dispatch and keyboard translation. The outer window has a caption, system
/// menu and minimize control. Its size is determined from current display
/// dimensions and the report uses a high-visibility monospaced text control.
pub fn run(){unsafe{
    let instance=GetModuleHandleW(null());
    let class_name=wide("FentonClassRust");
    // The class name identifies the application's own Win32 window class;
    // the standard icon/cursor are loaded from the operating system.
    let wc=WndClassEx{
        cb_size:std::mem::size_of::<WndClassEx>() as u32,style:0,proc:Some(wnd_proc),
        cls_extra:0,wnd_extra:0,instance,
        icon:LoadIconW(0,32512usize as *const u16),
        cursor:LoadCursorW(0,32512usize as *const u16),
        background:COLOR_WINDOW+1,menu_name:null(),class_name:class_name.as_ptr(),icon_sm:0,
    };
    if RegisterClassExW(&wc)==0{return;}
    let style=WS_CAPTION|WS_SYSMENU|WS_MINIMIZEBOX;
    let screen_width=GetSystemMetrics(SM_CXSCREEN);
    let screen_height=GetSystemMetrics(SM_CYSCREEN);
    let w=(screen_width-80).min(1600).max(600);
    let h=(screen_height-120).min(950).max(400);
    let mut rect=Rect{left:0,top:0,right:w,bottom:h};
    AdjustWindowRectEx(&mut rect,style,0,0);
    let hwnd=CreateWindowExW(0,class_name.as_ptr(),wide("Fenton Wave Solver").as_ptr(),style,
        i32::MIN,i32::MIN,rect.right-rect.left,rect.bottom-rect.top,0,0,instance,null_mut());
    if hwnd==0 {return;}
    ShowWindow(hwnd,1);
    UpdateWindow(hwnd);
    let mut msg=Msg::default();
    // Standard message pump: GetMessageW blocks until a message arrives.
    // Zero/negative status ends dispatch and returns control to main().
    loop {
        let status=GetMessageW(&mut msg,0,0,0);
        if status<=0 {break;}
        TranslateMessage(&msg);
        DispatchMessageW(&msg);
    }
}}
}

// ===========================================================================
// ENTRY POINT
// ===========================================================================

// ============================================================================
// PROGRAM ENTRY POINT — GUI WHEN TARGET WINDOWS, CONSOLE WHEN HEADLESS
// ============================================================================
// Normal Windows GNU compilation (command at top of file) enters run() and
// uses the Windows subsystem with no extra console window. Non-Windows builds
// and explicitly headless builds parse H, T, d, Uc from command-line arguments
// and print the same engineering report while writing output.txt.
// Do not pass CLI arguments to a GUI release executable and expect them to be
// read: the GUI reads its four native input fields instead.
fn main() {
    #[cfg(all(windows,not(feature="headless")))]
    windows_gui::run();

    // Cross-platform command-line build; on Windows use --cfg feature="headless".
    // Optional command-line execution path: exactly four finite numbers
    // define one wave; without arguments, a test wave is used. An invalid
    // value exits with an error status. The GUI build does not enter this arm.
    #[cfg(any(not(windows),feature="headless"))]
    {
        let args=std::env::args().skip(1).collect::<Vec<_>>();
        let vals=if args.len()==4 {
            let mut parsed=Vec::new();
            for s in &args {
                match s.replace(',', ".").parse::<f64>() {
                    Ok(v) if v.is_finite()=>parsed.push(v),
                    _=>{eprintln!("Invalid numeric argument: {s}");std::process::exit(2);}
                }
            }
            parsed
        }else if args.is_empty(){vec![3.0,9.0,5.0,1.0]}
        else{eprintln!("Usage: fenton_gui [H T d Uc]\nExample: fenton_gui 3.0 9.0 5.0 1.0");std::process::exit(2);};
        if vals[0]<=0.0||vals[1]<=0.0||vals[2]<=0.0 {
            eprintln!("H, T, d must be positive.");std::process::exit(2);
        }
        let output=report::generate_output(vals[0],vals[1],vals[2],vals[3]);
        let _=std::fs::write("output.txt",&output);
        print!("{output}");
    }
}
