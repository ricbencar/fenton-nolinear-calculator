! -----------------------------------------------------------------------------
! fenton88.f90
!
! Production Fortran 2008 single-source implementation of the steady water-wave
! formulation described in: Fenton, J. D. (1988). The numerical solution of
! steady water wave problems. Computers & Geosciences, 14(3), 357-368.
!
! This source is adapted from the method and Appendix program in Fenton (1988).
! It is not a verbatim reproduction of the printed code. The numerical model,
! notation and output terminology follow the paper, while the implementation is
! organised as one program with internal procedures, uses allocatable arrays,
! explicit file I/O and an internal linear-system solver.
!
! Physical problem
! ----------------
! Periodic waves of wavelength lambda and crest-to-trough height H propagate
! steadily over water of mean depth d. A coordinate system is fixed under a
! wave crest and moves with the wave, so the flow is steady in that frame. The
! method applies Fenton's wavenumber non-dimensionalisation using k = 2*pi/lambda,
! together with gravitational acceleration g. For finite depth, kd is solved as
! part of the nonlinear system. For deep water, the paper's dummy convention
! kd = -1 is retained so that the same variable ordering can be used.
!
! Numerical formulation
! ---------------------
! The unknown vector z has 2N+10 components, following Fenton's Table 1:
!
!   z(1)             kd
!   z(2)             kH
!   z(3)             tau*(g*k)**0.5
!   z(4)             c*(k/g)**0.5
!   z(5)             c_E*(k/g)**0.5
!   z(6)             c_S*(k/g)**0.5
!   z(7)             ubar*(k/g)**0.5
!   z(8)             q*(k**3/g)**0.5
!   z(9)             r*k/g, where r = R - g*d
!   z(10:N+10)       k*eta_m at N+1 surface collocation points
!   z(N+11:2N+10)    stream-function Fourier coefficients B_j
!
! The residual vector f(i), stored in rhs(i), contains the same 2N+10
! equations and uses Fenton's variable ordering:
!
!   f(1)              depth/height residual:
!                     finite: z(2) - (H/d)*z(1)
!                     deep:   z(1) + 1, with dummy kd = -1
!                     variables: z(1)=kd, z(2)=kH
!
!   f(2)              wavelength or period residual:
!                     wavelength: z(2) - 2*pi*(H/lambda)
!                     period:     z(2) - [H/(g*tau**2)]*z(3)**2
!                     variables: z(2)=kH, z(3)=tau*(g*k)**0.5
!
!   f(3)              wave-speed definition residual:
!                     z(4)*z(3) - 2*pi
!                     variables: z(3)=tau*(g*k)**0.5, z(4)=c*(k/g)**0.5
!
!   f(4)              Eulerian-current frame residual:
!                     z(5) + z(7) - z(4)
!                     variables: z(4)=c*, z(5)=c_E*, z(7)=ubar*
!
!   f(5)              Stokes-current / volume-flux residual:
!                     finite: z(6) + z(7) - z(4) - z(8)/z(1)
!                     deep:   z(6) + z(7) - z(4)
!                     variables: z(1)=kd, z(4)=c*, z(6)=c_S*,
!                                z(7)=ubar*, z(8)=q*
!
!   f(6)              prescribed-current residual:
!                     z(5) - input*sqrt(z(2)) for an Eulerian input, or
!                     z(6) - input*sqrt(z(2)) for a Stokes input
!                     variables: z(2)=kH, z(5)=c_E*, z(6)=c_S*
!
!   f(7)              zero-mean free-surface residual:
!                     z(10) + z(N+10) + 2*sum[z(11:N+9)]
!                     variables: z(10:N+10)=k*eta_m
!
!   f(8)              wave-height residual:
!                     z(10) - z(N+10) - z(2)
!                     variables: z(2)=kH, z(10)=k*eta_0, z(N+10)=k*eta_N
!
!   f(9+m)            kinematic free-surface residual at collocation point m,
!                     m=0,...,N: psi_m - z(8) - z(7)*z(10+m)
!                     variables: z(7)=ubar*, z(8)=q*, z(10+m)=k*eta_m,
!                                z(N+11:2N+10)=B_j
!
!   f(N+10+m)         dynamic free-surface Bernoulli residual at point m,
!                     m=0,...,N:
!                     0.5*((-z(7)+u_m)**2 + v_m**2) + z(10+m) - z(9)
!                     variables: z(7)=ubar*, z(9)=r*k/g, z(10+m)=k*eta_m,
!                                z(N+11:2N+10)=B_j
!
! Here c*, c_E*, c_S*, ubar* and q* denote the dimensionless variables defined
! in the z-vector table above. u_m and v_m are the dimensionless velocity-series
! contributions at the free-surface collocation point m.
!
! The stream function satisfies incompressibility, irrotationality, periodicity
! and the bottom condition analytically. The nonlinear free-surface conditions
! are imposed at equally spaced points from crest to trough. Newton iteration is
! used, and the Jacobian is obtained by numerical differentiation as in Fenton's
! simplified coding strategy.
!
! Wave-height continuation
! ------------------------
! For high or long waves, a purely linear first estimate can be inadequate. The
! target wave height is therefore reached in M steps. The first step starts from
! linear wave theory; later steps are extrapolated from the two preceding
! converged solutions. This follows the procedure described in the paper.
!
! Input file: input.txt
! ---------------------
! Four free-format lines are read:
!
!   1) deep | finite,        H/d
!   2) wavelength | period,  H/lambda or H/(g*tau**2)
!   3) euler | stokes,       c_E/(gH)**0.5 or c_S/(gH)**0.5
!   4) N, M                 Fourier coefficients and height-continuation steps
!
! Example from Fenton's Table 2(a):
!
!   deep 0.
!   wavelength 0.09762055
!   euler 0.
!   20 5
!
! Example finite-depth period-specified case:
!
!   finite 0.60000000000000
!   period 0.00377672671473
!   euler 0.184365236507856
!   50 5
!
! Output file: output.txt
! -----------------------
! The output contains the iteration history, final solution variables, surface
! elevations, stream-function Fourier coefficients, integral quantities and, for
! finite depth, volume flux Q, Bernoulli constant R and momentum flux S. Small
! high-order Fourier coefficients are printed in scientific notation so that
! truncation behaviour is visible instead of being rounded to zero.
!
! Mean squared bed velocity
! -------------------------
! The reported mean squared bed velocity is evaluated as a phase average of the
! orbital horizontal velocity at the bed:
!
!   ub2 = < (U_b - c_E)**2 >,
!
! where U_b is the physical-frame horizontal velocity at y = -d, c_E is the
! Eulerian current component, and < > denotes averaging over one complete wave
! cycle. In the code this is kept in Fenton's k-based scaling, so the reported
! k-based value is ub2*k/g and the d-based value is ub2/(g*d). This numerical
! phase-average definition is non-negative by construction. The original
! Cokelet/Fenton algebraic expression 2*r*k/g - c**2*k/g is retained internally
! only where it is part of the published radiation-stress and wave-power formulae.
!
! Compilation
! -----------
!   gfortran -std=f2008 -O3 fenton88.f90 -o fenton88.exe
!
! Static standalone example:
!   gfortran -std=f2008 -O3 -static -static-libgcc -static-libgfortran \
!            -static-libquadmath fenton88.f90 -o fenton88.exe
! -----------------------------------------------------------------------------

! -----------------------------------------------------------------------------
! Basic types, constants and keyword handling.
!
! The program stores the input choices in wave_problem and the nonlinear
! solution arrays in wave_workspace. The names are longer than in the Appendix
! program, but the physical quantities are the same as Fenton's z-vector.
! -----------------------------------------------------------------------------
program steady
  use, intrinsic :: iso_fortran_env, only : real64, error_unit
  implicit none


  ! Floating-point kind used throughout the solver. Fenton's printed program
  ! used double precision; real64 gives the same intent explicitly.
  integer, parameter :: dp = real64
  ! Constants controlling the Newton iteration. The tolerance follows the
  ! Appendix program: 1.0E-3 during continuation and 1.0E-5 on the final step.
  real(dp), parameter :: pi_value = acos(-1.0_dp)
  real(dp), parameter :: newton_tolerance = 1.0e-9_dp
  integer, parameter :: max_newton_iterations = 100

  ! User-specified problem definition read from input.txt.
  type :: wave_problem
     ! 'deep' selects the infinite-depth equations; 'finite' solves for kd.
     character(len=16) :: depth_kind = ''
     ! 'wavelength' prescribes H/lambda; 'period' prescribes H/(g*tau**2).
     character(len=16) :: specified_kind = ''
     ! 'euler' prescribes c_E; 'stokes' prescribes c_S.
     character(len=16) :: current_kind = ''
     ! Target H/d. Used directly for finite depth; ignored for deep water.
     real(dp) :: target_height_over_depth = 0.0_dp
     ! Target H/lambda or H/(g*tau**2), according to specified_kind.
     real(dp) :: target_height_parameter = 0.0_dp
     ! Prescribed current magnitude, scaled with (gH)**0.5.
     real(dp) :: current_value = 0.0_dp
     ! N is both the number of Fourier coefficients and the number of
     ! crest-to-trough intervals used for the free-surface collocation grid.
     integer :: n = 0
     ! M in the paper: number of continuation steps in wave height.
     integer :: nstep = 0
     ! Total number of unknowns and residual equations: 2N+10.
     integer :: num = 0
  end type wave_problem

  ! Working arrays used by the nonlinear solver and output routines.
  type :: wave_workspace
     ! Solution vector z(1:2N+10), using Fenton's variable ordering.
     real(dp), allocatable :: z(:)
     ! Discrete Fourier transform of the final free-surface elevations.
     real(dp), allocatable :: y(:)
     ! Trigonometric tables for m*j*pi/N at collocation points.
     real(dp), allocatable :: cos_table(:)
     real(dp), allocatable :: sin_table(:)
     ! Last two converged continuation solutions, used for height extrapolation.
     real(dp), allocatable :: previous_solution(:,:)
  end type wave_workspace


  call run_solver('input.txt', 'output.txt')

contains


  ! Convert input keywords to lower case and left-justify them. This makes the
  ! four-line input format tolerant of capitalisation without altering the
  ! numerical interpretation of Fenton's data choices.
  pure subroutine normalize_keyword(word)
    character(len=*), intent(inout) :: word
    integer :: i, c

    word = adjustl(word)
    do i = 1, len(word)
       c = iachar(word(i:i))
       if (c >= iachar('A') .and. c <= iachar('Z')) then
          word(i:i) = achar(c + iachar('a') - iachar('A'))
       end if
    end do
  end subroutine normalize_keyword


  ! Solve A*x=b in place. On return, b contains x.
  !
  ! info = 0   successful solution
  ! info > 0   singular or numerically singular pivot detected
  ! info = -1  inconsistent matrix/vector dimensions
  subroutine solve_linear_system(a, b, info)
    real(dp), intent(inout) :: a(:,:)
    real(dp), intent(inout) :: b(:)
    integer, intent(out) :: info

    integer :: i, j, k, n, pivot_row
    real(dp) :: pivot_abs, candidate_abs, factor, tmp
    real(dp) :: row_tmp(size(b))

    n = size(b)
    info = 0

    if (size(a,1) /= n .or. size(a,2) /= n) then
       info = -1
       return
    end if

    ! Forward elimination with partial pivoting.
    do k = 1, n - 1
       pivot_row = k
       pivot_abs = abs(a(k,k))
       do i = k + 1, n
          candidate_abs = abs(a(i,k))
          if (candidate_abs > pivot_abs) then
             pivot_abs = candidate_abs
             pivot_row = i
          end if
       end do

       if (pivot_abs <= tiny(1.0_dp)) then
          info = k
          return
       end if

       if (pivot_row /= k) then
          row_tmp(:) = a(k,:)
          a(k,:) = a(pivot_row,:)
          a(pivot_row,:) = row_tmp(:)
          tmp = b(k)
          b(k) = b(pivot_row)
          b(pivot_row) = tmp
       end if

       do i = k + 1, n
          factor = a(i,k)/a(k,k)
          a(i,k) = 0.0_dp
          do j = k + 1, n
             a(i,j) = a(i,j) - factor*a(k,j)
          end do
          b(i) = b(i) - factor*b(k)
       end do
    end do

    if (abs(a(n,n)) <= tiny(1.0_dp)) then
       info = n
       return
    end if

    ! Back substitution.
    do i = n, 1, -1
       tmp = b(i)
       do j = i + 1, n
          tmp = tmp - a(i,j)*b(j)
       end do
       b(i) = tmp/a(i,i)
    end do

  end subroutine solve_linear_system


  ! Main driver used by program steady. It deliberately keeps the user-facing
  ! workflow simple: read input.txt, solve the steady-wave problem, write
  ! output.txt.
  subroutine run_solver(input_file, output_file)
    character(len=*), intent(in) :: input_file, output_file

    type(wave_problem) :: problem
    type(wave_workspace) :: work
    integer :: input_unit, output_unit_local, ios

    ! Read the four-line Fenton-style input file.
    open(newunit=input_unit, file=input_file, status='old', action='read', &
         iostat=ios)
    if (ios /= 0) then
       write(error_unit,'(a)') 'ERROR: cannot open input.txt in the current directory.'
       stop 1
    end if

    call read_problem(input_unit, problem)
    close(input_unit)

    call validate_problem(problem)
    call allocate_workspace(problem, work)

    ! Replace output.txt on every run, matching the usual batch-solver workflow.
    open(newunit=output_unit_local, file=output_file, status='replace', &
         action='write', iostat=ios)
    if (ios /= 0) then
       write(error_unit,'(a)') 'ERROR: cannot create output.txt in the current directory.'
       stop 1
    end if

    call write_header(output_unit_local)
    call solve_problem(problem, work, output_unit_local)
    close(output_unit_local)
  end subroutine run_solver

  ! Parse the four free-format input lines. The numerical values are not scaled
  ! here; wave-height continuation applies the scaling later.
  subroutine read_problem(unit_in, problem)
    integer, intent(in) :: unit_in
    type(wave_problem), intent(out) :: problem
    integer :: ios

    read(unit_in, *, iostat=ios) problem%depth_kind, &
         problem%target_height_over_depth
    call assert_read_ok(ios, 'input line 1')

    read(unit_in, *, iostat=ios) problem%specified_kind, &
         problem%target_height_parameter
    call assert_read_ok(ios, 'input line 2')

    read(unit_in, *, iostat=ios) problem%current_kind, problem%current_value
    call assert_read_ok(ios, 'input line 3')

    read(unit_in, *, iostat=ios) problem%n, problem%nstep
    call assert_read_ok(ios, 'input line 4')

    call normalize_keyword(problem%depth_kind)
    call normalize_keyword(problem%specified_kind)
    call normalize_keyword(problem%current_kind)

    problem%num = 2*problem%n + 10
  end subroutine read_problem

  ! Abort early on malformed input, before any allocation or numerical work.
  subroutine assert_read_ok(ios, label)
    integer, intent(in) :: ios
    character(len=*), intent(in) :: label

    if (ios /= 0) then
       write(error_unit,'(a,a)') 'ERROR: invalid or incomplete ', trim(label)
       stop 1
    end if
  end subroutine assert_read_ok

  ! Validate only the structural input choices. Physical feasibility and
  ! convergence limits are determined by the nonlinear wave solution itself.
  subroutine validate_problem(problem)
    type(wave_problem), intent(in) :: problem

    select case (trim(problem%depth_kind))
    case ('deep', 'finite')
    case default
       write(error_unit,'(a)') "ERROR: depth keyword must be 'deep' or 'finite'."
       stop 1
    end select

    select case (trim(problem%specified_kind))
    case ('wavelength', 'period')
    case default
       write(error_unit,'(a)') "ERROR: case keyword must be 'wavelength' or 'period'."
       stop 1
    end select

    select case (trim(problem%current_kind))
    case ('euler', 'stokes')
    case default
       write(error_unit,'(a)') "ERROR: current keyword must be 'euler' or 'stokes'."
       stop 1
    end select

    if (problem%n < 1) then
       write(error_unit,'(a)') 'ERROR: N must be at least 1.'
       stop 1
    end if

    if (problem%nstep < 1) then
       write(error_unit,'(a)') 'ERROR: number of height steps must be at least 1.'
       stop 1
    end if

    if (trim(problem%depth_kind) == 'finite' .and. &
        problem%target_height_over_depth <= 0.0_dp) then
       write(error_unit,'(a)') 'ERROR: finite-depth H/d must be positive.'
       stop 1
    end if

    if (problem%target_height_parameter <= 0.0_dp) then
       write(error_unit,'(a)') 'ERROR: wave-height parameter must be positive.'
       stop 1
    end if
  end subroutine validate_problem

  ! Allocate arrays after N is known. The original printed program used fixed
  ! array dimensions; allocatable arrays avoid arbitrary hard limits on N.
  subroutine allocate_workspace(problem, work)
    type(wave_problem), intent(in) :: problem
    type(wave_workspace), intent(out) :: work
    integer :: m

    allocate(work%z(problem%num))
    allocate(work%y(problem%n))
    allocate(work%cos_table(0:2*problem%n))
    allocate(work%sin_table(0:2*problem%n))
    allocate(work%previous_solution(problem%num,2))

    work%z = 0.0_dp
    work%y = 0.0_dp
    work%previous_solution = 0.0_dp

    ! Store cos(m*pi/N) and sin(m*pi/N). Products m*j are reduced modulo 2N.
    do m = 0, 2*problem%n
       work%cos_table(m) = cos(real(m,dp)*pi_value/real(problem%n,dp))
       work%sin_table(m) = sin(real(m,dp)*pi_value/real(problem%n,dp))
    end do
  end subroutine allocate_workspace

  ! Bibliographic header written at the top of output.txt.
  subroutine write_header(unit_out)
    integer, intent(in) :: unit_out

    write(unit_out,'(a)') 'FENTON88 STEADY WATER WAVE SOLVER'
    write(unit_out,'(a)') 'Fenton, J. D. (1988). The numerical solution of steady water wave problems.'
    write(unit_out,'(a)') 'Computers & Geosciences, 14(3), 357-368.'
    write(unit_out,'(a)') '--------------------------------------------------------------------------'
    write(unit_out,'()')
  end subroutine write_header

  ! Advance from zero/small wave height to the requested wave height using the
  ! continuation strategy described by Fenton.
  subroutine solve_problem(problem, work, unit_out)
    type(wave_problem), intent(in) :: problem
    type(wave_workspace), intent(inout) :: work
    integer, intent(in) :: unit_out

    integer :: ns, i
    real(dp) :: step_height_parameter, step_height_over_depth
    real(dp) :: height_parameter, height_over_depth

    ! The prescribed height parameter and H/d are scaled in equal increments.
    step_height_parameter = problem%target_height_parameter/real(problem%nstep,dp)
    step_height_over_depth = problem%target_height_over_depth/real(problem%nstep,dp)

    call write_input_summary(problem, unit_out)

    do ns = 1, problem%nstep
       height_parameter = real(ns,dp)*step_height_parameter
       height_over_depth = real(ns,dp)*step_height_over_depth
       write(unit_out,23) ns, problem%nstep

       if (ns == 1) then
          ! First height level: linear wave theory supplies z^(0).
          call initialize_linear_solution(problem, work, height_parameter, &
                                          height_over_depth, unit_out)
       else
          ! Later height levels: linear extrapolation from the two previous
          ! converged nonlinear solutions.
          do i = 1, problem%num
             work%z(i) = 2.0_dp*work%previous_solution(i,2) - &
                         work%previous_solution(i,1)
          end do
       end if

       call newton_solve(problem, work, height_parameter, height_over_depth, &
                         ns, unit_out)

       if (ns == 1) then
          work%previous_solution(:,2) = work%z(:)
       else
          work%previous_solution(:,1) = work%previous_solution(:,2)
          work%previous_solution(:,2) = work%z(:)
       end if
    end do

    call output_solution(problem, work, unit_out)

23  format(/,'Height step ',i2,' of ',i2)
  end subroutine solve_problem

  ! Echo the input choices using Fenton's terminology.
  subroutine write_input_summary(problem, unit_out)
    type(wave_problem), intent(in) :: problem
    integer, intent(in) :: unit_out

    write(unit_out,'(a)') 'Input data'
    write(unit_out,'(a,a,", Height/Depth ",f10.6)') 'Depth: ', &
         trim(problem%depth_kind), problem%target_height_over_depth

    if (trim(problem%specified_kind) == 'period') then
       write(unit_out,'(a,f16.12)') 'Period parameter H/(g*T**2): ', &
            problem%target_height_parameter
    else
       write(unit_out,'(a,f16.12)') 'Wave steepness H/L:          ', &
            problem%target_height_parameter
    end if

    write(unit_out,'(a,a,a,f16.12)') 'Current criterion ', &
         trim(problem%current_kind), ', Magnitude ', problem%current_value
    write(unit_out,'(a,i0)') 'Fourier components N:        ', problem%n
    write(unit_out,'(a,i0)') 'Wave-height steps:           ', problem%nstep
    write(unit_out,'()')
  end subroutine write_input_summary

  ! Build the initial estimate from linear wave theory.
  !
  ! For finite-depth period input, kd is estimated with Eckart's approximation
  ! followed by one Newton refinement, as described in the paper. For wavelength
  ! input, kd = 2*pi*d/lambda follows directly. Deep water uses kd = -1 and the
  ! deep-water linear relation.
  subroutine initialize_linear_solution(problem, work, height_parameter, &
                                        height_over_depth, unit_out)
    type(wave_problem), intent(in) :: problem
    type(wave_workspace), intent(inout) :: work
    real(dp), intent(in) :: height_parameter, height_over_depth
    integer, intent(in) :: unit_out

    integer :: i
    real(dp) :: alpha, beta, tanh_beta

    work%z = 0.0_dp

    ! z(1)=kd, z(2)=kH and z(4)=c*(k/g)**0.5 are the primary
    ! linear-theory quantities required to initialise the remaining variables.
    if (trim(problem%depth_kind) == 'finite') then
       if (trim(problem%specified_kind) == 'period') then
          ! alpha = 4*pi**2*d/(g*tau**2). With height_parameter=H/(g*tau**2)
          ! and height_over_depth=H/d, alpha is obtained as shown below.
          alpha = 4.0_dp*pi_value*pi_value*height_parameter/height_over_depth
          beta = alpha/sqrt(tanh(alpha))
          tanh_beta = tanh(beta)
          work%z(1) = beta + (alpha - beta*tanh_beta) / &
               (tanh_beta + beta*(1.0_dp - tanh_beta*tanh_beta))
       else
          work%z(1) = 2.0_dp*pi_value*height_parameter/height_over_depth
       end if
       work%z(2) = work%z(1)*height_over_depth
       work%z(4) = sqrt(tanh(work%z(1)))
    else
       work%z(1) = -1.0_dp
       work%z(4) = 1.0_dp
       if (trim(problem%specified_kind) == 'period') then
          ! Deep-water linear period relation gives kH directly from
          ! H/(g*tau**2): kH = 4*pi**2*H/(g*tau**2).
          work%z(2) = 4.0_dp*pi_value*pi_value*height_parameter
       else
          work%z(2) = 2.0_dp*pi_value*height_parameter
       end if
    end if

    ! Linear wave period in the physical frame: tau*(g*k)**0.5 = 2*pi/c*.
    work%z(3) = 2.0_dp*pi_value/work%z(4)

    ! Only one current criterion is prescribed by the input. The other is set to
    ! zero as an initial estimate and is subsequently solved by the full system.
    if (trim(problem%current_kind) == 'euler') then
       work%z(5) = problem%current_value*sqrt(work%z(2))
       work%z(6) = 0.0_dp
    else
       work%z(6) = problem%current_value*sqrt(work%z(2))
       work%z(5) = 0.0_dp
    end if

    ! Linear estimates for mean relative speed, wave-induced volume flux,
    ! Bernoulli constant, surface ordinates and stream-function coefficients.
    work%z(7) = work%z(4)
    work%z(8) = 0.0_dp
    work%z(9) = 0.5_dp*work%z(7)**2
    work%z(10) = 0.5_dp*work%z(2)

    do i = 1, problem%n
       work%z(10+i) = 0.5_dp*work%z(2)*work%cos_table(i)
       work%z(problem%n+10+i) = 0.0_dp
    end do

    work%z(problem%n+11) = 0.5_dp*work%z(2)/work%z(7)

    write(unit_out,2) (work%z(i), i = 1, problem%num)
2   format('Initial linear solution',20(/,6(1x,e13.6)),/)

    work%previous_solution = 0.0_dp
    work%previous_solution(1:9,1) = work%z(1:9)
  end subroutine initialize_linear_solution

  ! Newton solution of f(z)=0 for one continuation height.
  !
  ! The Jacobian is approximated column by column by perturbing each unknown.
  ! This is the central simplification of Fenton's 1988 Appendix code relative
  ! to approaches that require an analytically coded Jacobian.
  subroutine newton_solve(problem, work, height_parameter, height_over_depth, &
                          height_step, unit_out)
    type(wave_problem), intent(in) :: problem
    type(wave_workspace), intent(inout) :: work
    real(dp), intent(in) :: height_parameter, height_over_depth
    integer, intent(in) :: height_step, unit_out

    integer :: i, j, iter, info
    real(dp) :: h, correction_sum, tolerance
    real(dp), allocatable :: rhs1(:), rhs2(:), jacobian(:,:), correction(:)

    allocate(rhs1(problem%num), rhs2(problem%num))
    allocate(jacobian(problem%num, problem%num), correction(problem%num))

    do iter = 1, max_newton_iterations
       write(unit_out,24) iter

       ! Residual vector f(z) at the current Newton iterate.
       call evaluate_equations(problem, work, height_parameter, &
                               height_over_depth, rhs1)

       ! Numerical differentiation: delta_i=z_i/100 unless z_i is very small,
       ! in which case the absolute increment 1.0E-5 is used.
       do i = 1, problem%num
          h = 0.01_dp*work%z(i)
          if (abs(work%z(i)) < 1.0e-4_dp) h = 1.0e-5_dp

          work%z(i) = work%z(i) + h
          call evaluate_equations(problem, work, height_parameter, &
                                  height_over_depth, rhs2)
          work%z(i) = work%z(i) - h

          correction(i) = -rhs1(i)
          do j = 1, problem%num
             jacobian(j,i) = (rhs2(j) - rhs1(j))/h
          end do
       end do

       ! Solve J*delta_z = -f. The correction vector overwrites correction(:).
       call solve_linear_system(jacobian, correction, info)
       if (info /= 0) then
          write(unit_out,'(a,i0)') 'Matrix singular. Pivot/status code: ', info
          stop
       end if

       ! Fenton's convergence measure is the sum of absolute corrections.
       correction_sum = sum(abs(correction))
       work%z(:) = work%z(:) + correction(:)

       write(unit_out,25) (work%z(i), i = 1, problem%num)

       tolerance = newton_tolerance
       if (height_step == problem%nstep) tolerance = 0.01_dp*newton_tolerance
       if (correction_sum < tolerance) then
          deallocate(rhs1, rhs2, jacobian, correction)
          return
       end if
    end do

    write(unit_out,26) max_newton_iterations
    stop

24  format(/,'Iteration ',i3)
25  format('Solution vector',20(/,6(1x,e13.6)),/)
26  format('Did not converge sufficiently after ',i3,' iterations.')
  end subroutine newton_solve

  ! Evaluate the 2N+10 nonlinear residual equations.
  !
  ! Residuals 1-8 are geometric/kinematic constraints and current criteria.
  ! Residuals 9 to N+9 impose that the free surface is a streamline. Residuals
  ! N+10 to 2N+10 impose constant pressure on the free surface through
  ! Bernoulli's equation.
  subroutine evaluate_equations(problem, work, height_parameter, &
                                height_over_depth, rhs)
    type(wave_problem), intent(in) :: problem
    type(wave_workspace), intent(in) :: work
    real(dp), intent(in) :: height_parameter, height_over_depth
    real(dp), intent(out) :: rhs(:)

    integer :: i, current_index, j, m, nm
    real(dp) :: b_j, e_term, psi, s_ratio, c_ratio
    real(dp) :: u_velocity, v_velocity
    real(dp) :: tanh_jkd(problem%n)

    rhs = 0.0_dp

    ! f1: relation between kH and kd through the specified H/d. For deep water
    ! the dummy equation kd+1=0 retains Fenton's variable ordering.
    if (trim(problem%depth_kind) == 'finite') then
       rhs(1) = work%z(2) - work%z(1)*height_over_depth
    else
       rhs(1) = work%z(1) + 1.0_dp
    end if

    ! f2: either prescribed steepness H/lambda or prescribed period parameter.
    if (trim(problem%specified_kind) == 'wavelength') then
       rhs(2) = work%z(2) - 2.0_dp*pi_value*height_parameter
    else
       rhs(2) = work%z(2) - height_parameter*work%z(3)**2
    end if

    ! f3: definition c=lambda/tau in dimensionless form.
    rhs(3) = work%z(4)*work%z(3) - 2.0_dp*pi_value

    ! f4 and f5: transform between the steady wave frame and the physical
    ! frame for Eulerian mean velocity and mean mass-transport velocity.
    rhs(4) = work%z(5) + work%z(7) - work%z(4)
    rhs(5) = work%z(6) + work%z(7) - work%z(4)

    ! Finite depth includes the wave-induced volume flux q/kd term.  The
    ! hyperbolic ratios appearing in the stream-function series are evaluated
    ! through tanh(jkd), avoiding the large cosh(jkd) denominator used in the
    ! printed Appendix program. This is algebraically equivalent to
    !
    !   sinh[j(kd+k*eta)]/cosh(jkd)
    !   cosh[j(kd+k*eta)]/cosh(jkd)
    !
    ! but is numerically safer for moderate-to-large kd or high Fourier order.
    if (trim(problem%depth_kind) == 'finite') then
       rhs(5) = rhs(5) - work%z(8)/work%z(1)
       do j = 1, problem%n
          tanh_jkd(j) = tanh(real(j,dp)*work%z(1))
       end do
    else
       tanh_jkd = 0.0_dp
    end if

    ! f6: enforce the selected input current criterion, c_E or c_S.
    current_index = 6
    if (trim(problem%current_kind) == 'euler') current_index = 5
    rhs(6) = work%z(current_index) - problem%current_value*sqrt(work%z(2))

    ! f7: zero mean free-surface elevation using the Fourier/trapezoidal sum.
    rhs(7) = work%z(10) + work%z(problem%n+10)
    do i = 1, problem%n - 1
       rhs(7) = rhs(7) + 2.0_dp*work%z(10+i)
    end do

    ! f8: crest-to-trough elevation difference equals kH.
    rhs(8) = work%z(10) - work%z(problem%n+10) - work%z(2)

    ! Surface collocation from crest m=0 to trough m=N. At each point the
    ! stream-function value gives the kinematic residual and Bernoulli's equation
    ! gives the dynamic residual.
    do m = 0, problem%n
       psi = 0.0_dp
       u_velocity = 0.0_dp
       v_velocity = 0.0_dp

       if (trim(problem%depth_kind) == 'finite') then
          do j = 1, problem%n
             nm = mod(m*j, 2*problem%n)
             b_j = work%z(problem%n+10+j)
             call finite_depth_shape_functions(real(j,dp), work%z(10+m), &
                                               tanh_jkd(j), s_ratio, c_ratio)

             psi = psi + b_j*s_ratio*work%cos_table(nm)
             u_velocity = u_velocity + real(j,dp)*b_j*c_ratio*work%cos_table(nm)
             v_velocity = v_velocity + real(j,dp)*b_j*s_ratio*work%sin_table(nm)
          end do
       else
          ! Deep-water terms reduce to exp(j*k*eta).
          do j = 1, problem%n
             nm = mod(m*j, 2*problem%n)
             b_j = work%z(problem%n+10+j)
             e_term = exp(real(j,dp)*work%z(10+m))
             psi = psi + b_j*e_term*work%cos_table(nm)
             u_velocity = u_velocity + real(j,dp)*b_j*e_term*work%cos_table(nm)
             v_velocity = v_velocity + real(j,dp)*b_j*e_term*work%sin_table(nm)
          end do
       end if

       ! Kinematic condition: free surface is a streamline.
       rhs(m+9) = psi - work%z(8) - work%z(7)*work%z(m+10)

       ! Dynamic condition: pressure is constant on the free surface.
       rhs(problem%n+m+10) = &
            0.5_dp*((-work%z(7) + u_velocity)**2 + v_velocity**2) + &
            work%z(m+10) - work%z(9)
    end do

  end subroutine evaluate_equations

  ! Stable finite-depth hyperbolic ratios used in the stream-function series:
  !
  !   S_j = sinh[j(kd+ky)]/cosh(jkd)
  !   C_j = cosh[j(kd+ky)]/cosh(jkd)
  !
  ! with ky = k*y measured from the still-water level.  The identities below
  ! avoid explicitly forming sinh(jkd) or cosh(jkd), which can overflow for
  ! large jkd even though the final ratios are well conditioned.
  pure subroutine finite_depth_shape_functions(j_value, ky, tanh_jkd, &
                                               s_ratio, c_ratio)
    real(dp), intent(in) :: j_value, ky, tanh_jkd
    real(dp), intent(out) :: s_ratio, c_ratio

    real(dp) :: theta, sinh_theta, cosh_theta

    theta = j_value*ky
    sinh_theta = sinh(theta)
    cosh_theta = cosh(theta)

    s_ratio = sinh_theta + cosh_theta*tanh_jkd
    c_ratio = cosh_theta + sinh_theta*tanh_jkd
  end subroutine finite_depth_shape_functions

  ! Final report. Quantities are written using Fenton's terminology and
  ! dimensionless groups. Finite-depth output includes both k-based and d-based
  ! values because these are the forms most convenient for engineering use.
  subroutine output_solution(problem, work, unit_out)
    type(wave_problem), intent(in) :: problem
    type(wave_workspace), intent(inout) :: work
    integer, intent(in) :: unit_out

    integer :: i, j, m, iend
    real(dp) :: sqkd, kd32, kd52, pulse, kinetic_energy
    real(dp) :: potential_energy, bed_velocity_squared, fenton_bed_velocity_squared
    real(dp) :: radiation_stress, wave_power
    real(dp) :: volume_flux, bernoulli_constant, momentum_flux

    sqkd = 0.0_dp
    kd32 = 0.0_dp
    kd52 = 0.0_dp

    ! Compute Y_j, the Fourier coefficients of the surface elevation. These are
    ! used by point_values for local velocity/pressure evaluation.
    call compute_surface_fourier_coefficients(problem, work)

    write(unit_out,'()')
    if (trim(problem%depth_kind) == 'finite') then
       sqkd = sqrt(work%z(1))
       kd32 = work%z(1)*sqkd
       kd52 = work%z(1)*work%z(1)*sqkd

       write(unit_out,'(a)') &
            'Solution, non-dimensionalized by wavenumber and by water depth'
       write(unit_out,'(a)') &
            'Variable                               k-based value       d-based value'
       write(unit_out,3) 'Water depth', work%z(1), 1.0_dp
       write(unit_out,3) 'Wave height', work%z(2), work%z(2)/work%z(1)
       write(unit_out,3) 'Wave period', work%z(3), work%z(3)/sqkd
       write(unit_out,3) 'Wavelength', 2.0_dp*pi_value, 2.0_dp*pi_value/work%z(1)
       write(unit_out,3) 'Wave speed', work%z(4), work%z(4)/sqkd
       write(unit_out,3) 'Mean Eulerian fluid velocity', work%z(5), work%z(5)/sqkd
       write(unit_out,3) 'Mean mass transport velocity', work%z(6), work%z(6)/sqkd
       write(unit_out,3) 'Mean fluid speed relative to wave', work%z(7), &
            work%z(7)/sqkd
       write(unit_out,3) 'Volume flux due to waves', work%z(8), work%z(8)/kd32
       write(unit_out,3) 'Bernoulli constant', work%z(9), work%z(9)/work%z(1)
    else
       write(unit_out,'(a)') 'Solution, non-dimensionalized by wavenumber'
       write(unit_out,'(a)') &
            'Variable                               k-based value       d-based value'
       write(unit_out,42) 'Wave height', work%z(2), 'not applicable'
       write(unit_out,42) 'Wave period', work%z(3), 'not applicable'
       write(unit_out,42) 'Wavelength', 2.0_dp*pi_value, 'not applicable'
       write(unit_out,42) 'Wave speed', work%z(4), 'not applicable'
       write(unit_out,42) 'Mean Eulerian fluid velocity', work%z(5), 'not applicable'
       write(unit_out,42) 'Mean mass transport velocity', work%z(6), &
            'not applicable'
       write(unit_out,42) 'Mean fluid speed relative to wave', work%z(7), &
            'not applicable'
       write(unit_out,42) 'Volume flux due to waves', work%z(8), 'not applicable'
       write(unit_out,42) 'Bernoulli constant', work%z(9), 'not applicable'
    end if

    write(unit_out,'()')
    write(unit_out,'(a)') 'Surface elevations - crest to trough'
    do i = 10, problem%n + 10, 10
       iend = min(i + 9, problem%n + 10)
       write(unit_out,'(10(1x,f9.4))') (work%z(j), j = i, iend)
    end do

    write(unit_out,'()')
    write(unit_out,'(a)') 'Fourier coefficients'
    do i = 1, problem%n, 4
       iend = min(i + 3, problem%n)
       write(unit_out,'(4(i4,1x,es16.8,3x))') &
            (j, work%z(problem%n+10+j), j = i, iend)
    end do

    ! Integral quantities in the notation of the paper. The variable name
    ! 'pulse' is the dimensionless impulse I/(rho*(g/k**3)**0.5).
    pulse = work%z(8) + work%z(1)*work%z(5)
    kinetic_energy = 0.5_dp*(work%z(4)*pulse + &
         work%z(5)*(work%z(8) - work%z(7)*work%z(1)))

    potential_energy = 0.5_dp*(work%z(10)**2 + work%z(problem%n+10)**2)
    do m = 1, problem%n - 1
       potential_energy = potential_energy + work%z(10+m)**2
    end do
    potential_energy = potential_energy/(2.0_dp*real(problem%n,dp))

    ! Mean squared bed orbital velocity used for the printed output row.
    ! This follows the GUI convention: ub2 = <(U_b - c_E)**2>, evaluated
    ! by uniform phase averaging over one complete wave cycle. The result is
    ! therefore non-negative and represents the oscillatory bed velocity after
    ! removing the prescribed Eulerian current component.
    bed_velocity_squared = mean_square_bed_orbital_velocity(problem, work)

    ! Algebraic bed-velocity term appearing in the Cokelet/Fenton integral
    ! formulae for radiation stress and wave power. It is kept separate from
    ! the reported phase-averaged orbital ub2 so no other integral formula is
    ! altered by the output-definition correction above.
    fenton_bed_velocity_squared = 2.0_dp*work%z(9) - work%z(4)**2

    radiation_stress = 4.0_dp*kinetic_energy - 3.0_dp*potential_energy + &
         fenton_bed_velocity_squared*work%z(1) + &
         2.0_dp*work%z(5)*(work%z(7)*work%z(1) - work%z(8))
    wave_power = work%z(4)*(3.0_dp*kinetic_energy - 2.0_dp*potential_energy) + &
         0.5_dp*fenton_bed_velocity_squared*(pulse + work%z(4)*work%z(1)) + &
         work%z(4)*work%z(5)*(work%z(7)*work%z(1) - work%z(8))

    write(unit_out,'()')
    if (trim(problem%depth_kind) == 'finite') then
       write(unit_out,'(a,/,a)') 'Integral quantities', &
            'Quantity                              k-based value       d-based value'
       write(unit_out,31) 'Impulse (I)', pulse, pulse/kd32
       write(unit_out,31) 'Kinetic energy (T)', kinetic_energy, &
            kinetic_energy/(work%z(1)*work%z(1))
       write(unit_out,31) 'Potential energy (V)', potential_energy, &
            potential_energy/(work%z(1)*work%z(1))
       write(unit_out,31) 'Mean squared bed velocity (ub2)', bed_velocity_squared, &
            bed_velocity_squared/work%z(1)
       write(unit_out,31) 'Radiation stress (Sxx)', radiation_stress, &
            radiation_stress/(work%z(1)*work%z(1))
       write(unit_out,31) 'Wave power (F)', wave_power, wave_power/kd52
    else
       write(unit_out,9) pulse, kinetic_energy, potential_energy, &
            bed_velocity_squared, radiation_stress, wave_power
    end if

    if (trim(problem%depth_kind) == 'finite') then
       ! Fenton reports the finite-depth invariants Q, R and S separately.
       volume_flux = work%z(7)*work%z(1) - work%z(8)
       bernoulli_constant = work%z(9) + work%z(1)
       momentum_flux = radiation_stress - 2.0_dp*work%z(4)*pulse + &
            (work%z(4)**2 + 0.5_dp*work%z(1))*work%z(1)

       write(unit_out,'()')
       write(unit_out,'(a,/,a)') 'Invariants for finite depth', &
            'Quantity                              k-based value       d-based value'
       write(unit_out,31) 'Volume flux (Q)', volume_flux, volume_flux/kd32
       write(unit_out,31) 'Bernoulli constant (R)', bernoulli_constant, &
            bernoulli_constant/work%z(1)
       write(unit_out,31) 'Momentum flux (S)', momentum_flux, &
            momentum_flux/(work%z(1)*work%z(1))
    end if

    call write_final_residuals(problem, work, unit_out)

3   format(a35,1x,f16.6,1x,f16.6)
31  format(a35,1x,f16.6,1x,f16.6)
42  format(a35,1x,f16.6,1x,a16)
9   format('Integral quantities',/, &
           'Impulse (I)                    ',f16.6,/, &
           'Kinetic energy (T)             ',f16.6,/, &
           'Potential energy (V)           ',f16.6,/, &
           'Mean squared bed velocity      ',f16.6,/, &
           'Radiation stress (Sxx)         ',f16.6,/, &
           'Wave power (F)                 ',f16.6)
  end subroutine output_solution

  ! Print the final residual vector f(z) after all engineering output sections.
  !
  ! This is a complete diagnostic dump of the nonlinear equations actually solved
  ! by Newton's method.  The values should be close to zero at convergence. The
  ! labels intentionally follow the same f(i) numbering used by evaluate_equations
  ! so that the output can be checked directly against the code.
  subroutine write_final_residuals(problem, work, unit_out)
    type(wave_problem), intent(in) :: problem
    type(wave_workspace), intent(in) :: work
    integer, intent(in) :: unit_out

    integer :: i
    real(dp), allocatable :: residual(:)
    character(len=10) :: equation_text
    character(len=3) :: point_text
    character(len=18) :: condition_text
    character(len=14) :: variable_text

    allocate(residual(problem%num))
    call evaluate_equations(problem, work, problem%target_height_parameter, &
                            problem%target_height_over_depth, residual)

    write(unit_out,'()')
    write(unit_out,'(a)') 'Final residual vector f(z)'
    write(unit_out,'(a)') 'Complete residuals of the 2N+10 equations solved by Newton iteration.'
    write(unit_out,'(a,1x,es12.5)') 'Maximum absolute residual:', maxval(abs(residual))
    write(unit_out,'(a,1x,es12.5)') 'Root-mean-square residual:', &
         sqrt(sum(residual*residual)/real(problem%num,dp))
    write(unit_out,'()')
    write(unit_out,'(a)') 'Variable notation:'
    write(unit_out,'(a)') '  z1=kd, z2=kH, z3=tau*(gk)^0.5, z4=c*, z5=cE*'
    write(unit_out,'(a)') '  z6=cS*, z7=ubar*, z8=q*, z9=rk/g'
    write(unit_out,'(a)') '  eta_m means z(10+m)=k*eta_m; B_j means z(N+10+j).'
    write(unit_out,'()')
    write(unit_out,'(a)') &
         '  i  residual    m  equation           variables         value'
    write(unit_out,'(a)') repeat('-', 70)

    do i = 1, problem%num
       call residual_row_fields(problem, i, equation_text, point_text, &
                                condition_text, variable_text)
       write(unit_out,'(i3,2x,a10,2x,a3,2x,a18,2x,a14,2x,es12.4)') &
            i, equation_text, point_text, condition_text, variable_text, residual(i)
    end do

    deallocate(residual)
  end subroutine write_final_residuals

  ! Build fixed-width fields for one residual output row.  The output is
  ! intentionally split into short columns so that the residual table remains
  ! readable in plain text output.txt, even for N=50 or larger.
  subroutine residual_row_fields(problem, index_value, equation_text, point_text, &
                                 condition_text, variable_text)
    type(wave_problem), intent(in) :: problem
    integer, intent(in) :: index_value
    character(len=*), intent(out) :: equation_text
    character(len=*), intent(out) :: point_text
    character(len=*), intent(out) :: condition_text
    character(len=*), intent(out) :: variable_text

    integer :: m

    equation_text = ''
    point_text = '--'
    condition_text = ''
    variable_text = ''

    select case (index_value)
    case (1)
       equation_text = 'f1'
       if (trim(problem%depth_kind) == 'finite') then
          condition_text = 'Depth/height'
          variable_text = 'z1,z2'
       else
          condition_text = 'Deep-water dummy'
          variable_text = 'z1'
       end if

    case (2)
       equation_text = 'f2'
       if (trim(problem%specified_kind) == 'wavelength') then
          condition_text = 'Wavelength'
          variable_text = 'z2,H/lambda'
       else
          condition_text = 'Period'
          variable_text = 'z2,z3'
       end if

    case (3)
       equation_text = 'f3'
       condition_text = 'Wave speed'
       variable_text = 'z3,z4'

    case (4)
       equation_text = 'f4'
       condition_text = 'Eulerian frame'
       variable_text = 'z4,z5,z7'

    case (5)
       equation_text = 'f5'
       if (trim(problem%depth_kind) == 'finite') then
          condition_text = 'Stokes/flux'
          variable_text = 'z1,z4,z6,z7,z8'
       else
          condition_text = 'Stokes current'
          variable_text = 'z4,z6,z7'
       end if

    case (6)
       equation_text = 'f6'
       if (trim(problem%current_kind) == 'euler') then
          condition_text = 'Input Eulerian'
          variable_text = 'z2,z5'
       else
          condition_text = 'Input Stokes'
          variable_text = 'z2,z6'
       end if

    case (7)
       equation_text = 'f7'
       condition_text = 'Zero-mean eta'
       variable_text = 'z10:zN+10'

    case (8)
       equation_text = 'f8'
       condition_text = 'Wave height'
       variable_text = 'z2,z10,zN+10'

    case default
       if (index_value >= 9 .and. index_value <= problem%n + 9) then
          m = index_value - 9
          equation_text = 'f(9+m)'
          write(point_text,'(i0)') m
          condition_text = 'Kinematic FS'
          variable_text = 'eta_m,B,z7,z8'
       else if (index_value >= problem%n + 10 .and. index_value <= problem%num) then
          m = index_value - (problem%n + 10)
          equation_text = 'f(N+10+m)'
          write(point_text,'(i0)') m
          condition_text = 'Dynamic Bernoulli'
          variable_text = 'eta_m,B,z7,z9'
       else
          equation_text = 'unknown'
          point_text = '--'
          condition_text = 'Unassigned index'
          variable_text = '--'
       end if
    end select
  end subroutine residual_row_fields

  ! Discrete Fourier transform of k*eta_m values from crest to trough. The two
  ! endpoints receive half weight, consistent with the interpolation formula in
  ! the paper.
  subroutine compute_surface_fourier_coefficients(problem, work)
    type(wave_problem), intent(in) :: problem
    type(wave_workspace), intent(inout) :: work

    integer :: j, m
    real(dp) :: sum_value

    do j = 1, problem%n
       sum_value = 0.5_dp*(work%z(10) + &
            work%z(problem%n+10)*(-1.0_dp)**j)
       do m = 1, problem%n - 1
          sum_value = sum_value + work%z(10+m)* &
               work%cos_table(mod(m*j, 2*problem%n))
       end do
       work%y(j) = 2.0_dp*sum_value/real(problem%n,dp)
    end do
  end subroutine compute_surface_fourier_coefficients

  ! Mean squared bed orbital velocity in the output convention used by the
  ! associated GUI program. The calculation samples the physical-frame horizontal
  ! velocity at the bed over one complete phase cycle and removes the Eulerian
  ! current component before squaring.
  !
  ! Returned value:
  !   k-based ub2 = <(U_b - c_E)**2>*k/g.
  !
  ! For deep water there is no finite bed in the computational domain; the
  ! corresponding limiting orbital velocity at the bed is therefore zero.
  function mean_square_bed_orbital_velocity(problem, work) result(ub2)
    type(wave_problem), intent(in) :: problem
    type(wave_workspace), intent(in) :: work
    real(dp) :: ub2

    integer, parameter :: phase_count = 720
    integer :: iphase
    real(dp) :: kx, orbital_velocity

    if (trim(problem%depth_kind) == 'deep') then
       ub2 = 0.0_dp
       return
    end if

    ub2 = 0.0_dp

    do iphase = 0, phase_count - 1
       kx = 2.0_dp*pi_value*real(iphase,dp)/real(phase_count,dp)

       ! The helper returns the physical-frame horizontal velocity scaled as
       ! U*(k/g)**0.5 at y=-d.  z(5) is c_E*(k/g)**0.5.  Their difference is
       ! the orbital component at the bed, in the same k-based scaling.
       orbital_velocity = bed_horizontal_velocity(problem, work, kx) - work%z(5)
       ub2 = ub2 + orbital_velocity*orbital_velocity
    end do

    ub2 = ub2/real(phase_count,dp)
  end function mean_square_bed_orbital_velocity

  ! Horizontal velocity at the finite-depth bed in the physical frame.
  !
  ! Returned value:
  !   U_b*(k/g)**0.5, evaluated at y=-d.
  !
  ! The mean-squared-bed-velocity output uses this focused helper instead of the
  ! general point routine, avoiding unnecessary pressure, elevation and vertical
  ! velocity calculations during the phase average.
  function bed_horizontal_velocity(problem, work, kx) result(u_bed)
    type(wave_problem), intent(in) :: problem
    type(wave_workspace), intent(in) :: work
    real(dp), intent(in) :: kx
    real(dp) :: u_bed

    integer :: j
    real(dp) :: b_j, c_ratio, s_ratio, tanh_jkd

    u_bed = work%z(4) - work%z(7)

    if (trim(problem%depth_kind) == 'deep') return

    ! At the bed, ky=-kd. The stable ratio C_j is exactly sech(jkd), so the
    ! expression remains well conditioned even when jkd is large.
    do j = 1, problem%n
       b_j = work%z(problem%n+10+j)
       tanh_jkd = tanh(real(j,dp)*work%z(1))
       call finite_depth_shape_functions(real(j,dp), -work%z(1), tanh_jkd, &
                                         s_ratio, c_ratio)
       u_bed = u_bed + real(j,dp)*b_j*c_ratio*cos(real(j,dp)*kx)
    end do
  end function bed_horizontal_velocity

  ! Local velocities, pressure and surface elevation at a point.
  !
  ! Inputs kx and ky are dimensionless coordinates k*(X-c*t) and kY. Returned
  ! velocities are in the physical frame through which the wave travels:
  ! U*(k/g)**0.5 and V*(k/g)**0.5. Pressure is p/(rho*g/k).
  !
  ! The main report uses the same velocity representation indirectly for the
  ! phase-averaged bed-velocity row, and this general routine remains available
  ! for engineering calculations at arbitrary points in the flow.
  subroutine point_values(problem, work, kx, ky, u, v, press, elevn)
    type(wave_problem), intent(in) :: problem
    type(wave_workspace), intent(in) :: work
    real(dp), intent(in) :: kx, ky
    real(dp), intent(out) :: u, v, press, elevn

    integer :: j
    real(dp) :: b_j, e_term, s_ratio, c_ratio, tanh_jkd

    elevn = 0.5_dp*work%y(problem%n)*cos(real(problem%n,dp)*kx)
    do j = 1, problem%n - 1
       elevn = elevn + work%y(j)*cos(real(j,dp)*kx)
    end do

    u = work%z(4) - work%z(7)
    v = 0.0_dp

    if (trim(problem%depth_kind) == 'finite') then
       do j = 1, problem%n
          b_j = work%z(problem%n+10+j)
          tanh_jkd = tanh(real(j,dp)*work%z(1))
          call finite_depth_shape_functions(real(j,dp), ky, tanh_jkd, &
                                            s_ratio, c_ratio)
          u = u + real(j,dp)*b_j*c_ratio*cos(real(j,dp)*kx)
          v = v + real(j,dp)*b_j*s_ratio*sin(real(j,dp)*kx)
       end do
    else
       do j = 1, problem%n
          b_j = work%z(problem%n+10+j)
          e_term = exp(real(j,dp)*ky)
          u = u + real(j,dp)*b_j*e_term*cos(real(j,dp)*kx)
          v = v + real(j,dp)*b_j*e_term*sin(real(j,dp)*kx)
       end do
    end if

    press = work%z(9) - ky - 0.5_dp*((u - work%z(4))**2 + v*v)
  end subroutine point_values


end program steady
