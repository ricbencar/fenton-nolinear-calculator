! -----------------------------------------------------------------------------
! fenton88_steady.f90
!
! Clean transcription of the FORTRAN program "STEADY" printed in the Appendix of:
! J. D. Fenton (1988), "The numerical solution of steady water wave problems",
! Computers & Geosciences, 14(3), 357-368.
!
! Notes on this transcription:
!   - The algorithm, variable numbering, equations, input format and output format
!     follow the printed program.
!   - The printed program calls LINPACK routines DGEFA and DGESL. To make this file
!     standalone, compatible replacement subroutines with the same interfaces are
!     included at the end.
!
! Default files:
!   - input : input.txt
!   - output: output.txt
!   - maximum/default tested Fourier components: N = 50
!   - output.txt is spaced with blank lines between the main result sections.
!
! Example input from Fenton 1988 Table 2(a), to be saved as input.txt:
!   deep 0.
!   wavelength 0.09762055
!   euler 0.
!   10 1
!
! Additional finite-depth, period-specified example with Eulerian current:
!   finite 0.600000000000000
!   period 0.003776726714733
!   euler 0.184365236507856
!   50 5
!
! Run, without shell redirection:
!   fenton88_steady.exe
!
! Static standalone executable compilation, for example:
!   gfortran -std=legacy -O3 -static -static-libgcc -static-libgfortran -static-libquadmath fenton88_steady.f90 -o fenton88_steady.exe
! -----------------------------------------------------------------------------

program steady
  implicit double precision(a-h,k-l,o-z)
  character*10 depth, case, currnt
  integer, parameter :: maxn = 50
  integer, parameter :: maxvar = 2*maxn + 10

  common /one/ n, num, pi, hoverd, height, value, depth, case, currnt
  common /two/ z(maxvar), cosa(0:maxvar), sina(0:maxvar), coeff(maxvar), sol(maxvar,2), y(maxvar)

  dimension rhs1(maxvar), rhs2(maxvar), a(maxvar,maxvar), b(maxvar), ipvt(maxvar)

  ! DEFAULT FILES.
  ! The executable is intended to be run directly, without shell redirection.
  ! It reads input.txt from the working directory and writes output.txt there.
  open(unit=5, file='input.txt', status='old', action='read', iostat=ios)
  if (ios .ne. 0) then
     write(*,'(a)') 'ERROR: cannot open input.txt in the current working directory.'
     stop 1
  endif

  open(unit=6, file='output.txt', status='replace', action='write', iostat=ios)
  if (ios .ne. 0) then
     write(*,'(a)') 'ERROR: cannot create output.txt in the current working directory.'
     stop 1
  endif

  ! INPUT DATA
  ! "depth"  is either 'deep' or 'finite'.
  ! "hoverd" is wave height / depth.
  read(5,*) depth, hoverd

  ! "case" is either 'period' or 'wavelength'.
  ! "height" is height/length if "case" is 'wavelength'.
  ! "height" is height/(g*T**2) if "case" is 'period'.
  read(5,*) case, height

  ! "currnt" is either 'euler' or 'stokes'.
  ! "value" is the magnitude of the mean Eulerian or Stokes velocity,
  ! non-dimensionalized with respect to wave height.
  read(5,*) currnt, value

  ! "n" is the number of terms in the Fourier series and the number
  ! of intervals in half a wavelength.
  ! "nstep" is the number of steps in wave height.
  read(5,*) n, nstep

  ! "number" is the number of iterations for each wave-height step.
  number = 9

  ! "crit" is the convergence criterion. If the sum of magnitudes of
  ! corrections is smaller than crit, the iteration stops.
  crit = 1.d-3

  call write_header
  call write_blank

  write(6,20) depth, hoverd
  if (case .eq. 'period') then
     write(6,'(A,F16.12)') 'Period parameter H/(g*T**2): ', height
  else
     write(6,'(A,F16.12)') 'Wave steepness H/L:          ', height
  endif
  write(6,22) currnt, value
  write(6,'(A,I0)') 'Fourier components N:        ', n
  write(6,'(A,I0)') 'Wave-height steps:           ', nstep
  call write_blank

  num = 2*n + 10
  if (n .gt. maxn) then
     write(6,'(A,I0,A,I0)') 'N too large: ', n, '. Maximum N is ', maxn
     stop
  endif

  pi = 4.d0*atan(1.d0)
  dhe = height/nstep
  dho = hoverd/nstep

  ! COMMENCE STEPPING THROUGH STEPS IN WAVE HEIGHT.
  do 1 ns = 1, nstep
     write(6,23) ns, nstep
     height = ns*dhe
     hoverd = ns*dho

     ! CALCULATE INITIAL LINEAR SOLUTION.
     if (ns .le. 1) then
        call init
     else
        ! OR, EXTRAPOLATE FOR NEXT WAVE HEIGHT, IF NECESSARY.
        do 3 i = 1, num
3          z(i) = 2.d0*sol(i,2) - sol(i,1)
     endif

     ! COMMENCE ITERATIVE SOLUTION.
     do 4 iter = 1, number
        write(6,24) iter

        ! CALCULATE RIGHT SIDES OF EQUATIONS AND DIFFERENTIATE NUMERICALLY
        ! TO OBTAIN JACOBIAN MATRIX.
        call eqns(rhs1)
        do 5 i = 1, num
           h = 0.01d0*z(i)
           if (abs(z(i)) .lt. 1.d-4) h = 1.d-5
           z(i) = z(i) + h
           call eqns(rhs2)
           z(i) = z(i) - h
           b(i) = -rhs1(i)
           do 6 j = 1, num
6             a(j,i) = (rhs2(j) - rhs1(j))/h
5       continue

        ! SOLVE MATRIX EQUATION AND CORRECT VARIABLES, USING LINPACK-STYLE ROUTINES.
        ! The matrix equation [a(i,j)][correction vector] = [b(i)] is solved.
        call dgefa(a, maxvar, num, ipvt, info)
        if (info .ne. 0) then
           write(6,27)
           stop
        endif
        call dgesl(a, maxvar, num, ipvt, b, 0)

        ! The b(i) are now the corrections to each variable.
        sum = 0.d0
        do 7 i = 1, num
           sum = sum + abs(b(i))
7          z(i) = z(i) + b(i)

        write(6,25) (z(i), i = 1, num)
        criter = crit
        if (ns .eq. nstep) criter = 0.01d0*crit
        if (sum .lt. criter) goto 8
4    continue

     write(6,26) number
     stop

8    if (ns .eq. 1) then
        do 9 i = 1, num
9          sol(i,2) = z(i)
     else
        do 10 i = 1, num
           sol(i,1) = sol(i,2)
10         sol(i,2) = z(i)
     endif
1 continue

  ! OUTPUT OF RESULTS.
  call output

20 format('Input data',/,'Depth: ',a6,', Height/Depth ',f10.6)
22 format('Current criterion ',a6,', Magnitude ',f16.12)
23 format(/,'Height step ',i2,' of ',i2)
24 format(/,'Iteration ',i3)
25 format('Solution vector',20(/,6e13.6),/)
26 format('Did not converge sufficiently after ',i3,' iterations.')
27 format('Matrix singular')

  stop
end program steady


! -----------------------------------------------------------------------------
! SUBROUTINE TO WRITE PROGRAM IDENTIFICATION HEADER.
! -----------------------------------------------------------------------------
subroutine write_header
  implicit none

  write(6,'(a)') 'FENTON88 STEADY WATER WAVE SOLVER'
  write(6,'(a)') 'Fenton, J. D. (1988). The numerical solution of steady water wave problems.'
  write(6,'(a)') 'Computers & Geosciences, 14(3), 357-368.'
  write(6,'(a)') '--------------------------------------------------------------------------'

  return
end subroutine write_header

! -----------------------------------------------------------------------------
! SUBROUTINE TO WRITE ONE BLANK LINE IN OUTPUT.TXT.
! -----------------------------------------------------------------------------
subroutine write_blank
  implicit none

  write(6,'()')

  return
end subroutine write_blank

! -----------------------------------------------------------------------------
! SUBROUTINE TO CALCULATE INITIAL SOLUTION FROM LINEAR WAVE THEORY.
! -----------------------------------------------------------------------------
subroutine init
  implicit double precision(a-h,k-l,o-z)
  character*10 depth, case, currnt

  integer, parameter :: maxn = 50
  integer, parameter :: maxvar = 2*maxn + 10

  common /one/ n, num, pi, hoverd, height, value, depth, case, currnt
  common /two/ z(maxvar), cosa(0:maxvar), sina(0:maxvar), coeff(maxvar), sol(maxvar,2), y(maxvar)

  if (depth .eq. 'finite') then
     if (case .eq. 'period') then
        a = 4.d0*pi*pi*height/hoverd
        b = a/sqrt(tanh(a))
        t = tanh(b)
        z(1) = b + (a - b*t)/(t + b*(1.d0 - t*t))
     else
        z(1) = 2.d0*pi*height/hoverd
     endif
     z(2) = z(1)*hoverd
     z(4) = sqrt(tanh(z(1)))
  else
     z(1) = -1.d0
     z(4) = 1.d0
     if (case .eq. 'period') then
        z(2) = 4.d0*pi*pi*height
     else
        z(2) = 2.d0*pi*height
     endif
  endif

  z(3) = 2.d0*pi/z(4)

  if (currnt .eq. 'euler') then
     z(5) = value*sqrt(z(2))
     z(6) = 0.d0
  else
     z(6) = value*sqrt(z(2))
     z(5) = 0.d0
  endif

  z(7) = z(4)
  z(8) = 0.d0
  z(9) = 0.5d0*z(7)**2

  cosa(0) = 1.d0
  sina(0) = 0.d0
  z(10) = 0.5d0*z(2)

  do 1 i = 1, n
     cosa(i)   = cos(i*pi/n)
     cosa(i+n) = cos((i+n)*pi/n)
     sina(i)   = sin(i*pi/n)
     sina(i+n) = sin((i+n)*pi/n)
     z(n+i+10) = 0.d0
1    z(i+10)   = 0.5d0*z(2)*cosa(i)

  z(n+11) = 0.5d0*z(2)/z(7)

  write(6,2) (z(i), i = 1, num)
2 format('Initial linear solution',20(/,6e13.6),/)

  do 3 i = 1, 9
3    sol(i,1) = z(i)
  sol(i,2) = 0.d0
  do 4 i = 10, num
4    sol(i,1) = 0.d0

  return
end subroutine init

! -----------------------------------------------------------------------------
! SUBROUTINE FOR EVALUATION OF EQUATIONS.
! -----------------------------------------------------------------------------
subroutine eqns(rhs)
  implicit double precision(a-h,k-l,o-z)
  character*10 depth, case, currnt

  integer, parameter :: maxn = 50
  integer, parameter :: maxvar = 2*maxn + 10

  common /one/ n, num, pi, hoverd, height, value, depth, case, currnt
  common /two/ z(maxvar), cosa(0:maxvar), sina(0:maxvar), coeff(maxvar), sol(maxvar,2), y(maxvar)

  dimension rhs(maxvar)

  if (depth .eq. 'finite') then
     rhs(1) = z(2) - z(1)*hoverd
  else
     rhs(1) = z(1) + 1.d0
  endif

  if (case .eq. 'wavelength') then
     rhs(2) = z(2) - 2.d0*pi*height
  else
     rhs(2) = z(2) - height*z(3)**2
  endif

  rhs(3) = z(4)*z(3) - pi - pi
  rhs(4) = z(5) + z(7) - z(4)
  rhs(5) = z(6) + z(7) - z(4)

  if (depth .eq. 'finite') then
     rhs(5) = rhs(5) - z(8)/z(1)
     do 2 i = 1, n
2       coeff(i) = z(n+i+10)/cosh(i*z(1))
  endif

  it = 6
  if (currnt .eq. 'euler') it = 5
  rhs(6) = z(it) - value*sqrt(z(2))

  rhs(7) = z(10) + z(n+10)
  do 1 i = 1, n-1
1    rhs(7) = rhs(7) + z(10+i) + z(10+i)

  rhs(8) = z(10) - z(n+10) - z(2)

  do 3 m = 0, n
     psi = 0.d0
     u = 0.d0
     v = 0.d0

     if (depth .eq. 'finite') then
        do 4 j = 1, n
           nm = mod(m*j, n+n)
           e = exp(j*(z(1) + z(10+m)))
           s = 0.5d0*(e - 1.d0/e)
           c = 0.5d0*(e + 1.d0/e)
           psi = psi + coeff(j)*s*cosa(nm)
           u = u + j*coeff(j)*c*cosa(nm)
4          v = v + j*coeff(j)*s*sina(nm)
     else
        do 5 j = 1, n
           nm = mod(m*j, n+n)
           e = exp(j*z(10+m))
           psi = psi + z(n+j+10)*e*cosa(nm)
           u = u + j*z(n+j+10)*e*cosa(nm)
5          v = v + j*z(n+j+10)*e*sina(nm)
     endif

     rhs(m+9) = psi - z(8) - z(7)*z(m+10)
     rhs(n+m+10) = 0.5d0*((-z(7) + u)**2 + v**2) + z(m+10) - z(9)
3 continue

  return
end subroutine eqns

! -----------------------------------------------------------------------------
! SUBROUTINE FOR OUTPUT OF RESULTS.
! -----------------------------------------------------------------------------
subroutine output
  implicit double precision(a-h,k-l,o-z)
  character*10 depth, case, currnt

  integer, parameter :: maxn = 50
  integer, parameter :: maxvar = 2*maxn + 10

  common /one/ n, num, pi, hoverd, height, value, depth, case, currnt
  common /two/ z(maxvar), cosa(0:maxvar), sina(0:maxvar), coeff(maxvar), sol(maxvar,2), y(maxvar)

  ! Calculate Fourier coefficients of surface elevation.
  do 10 j = 1, n
     sum = 0.5d0*(z(10) + z(n+10)*(-1.d0)**j)
     do 11 m = 1, n-1
11      sum = sum + z(10+m)*cosa(mod(m*j, n+n))
10   y(j) = 2.d0*sum/n

  call write_blank
  write(6,1)
1 format('Solution, non-dimensionalized by wavenumber')

  if (depth .eq. 'finite') write(6,2) z(1)
2 format('Water depth                         ',f10.6)

  write(6,3) (z(i), i = 2, 9)
3 format('Wave height                         ',f10.6,/, &
         'Wave period                         ',f10.6,/, &
         'Wave speed                          ',f10.6,/, &
         'Mean Eulerian fluid velocity        ',f10.6,/, &
         'Mean mass transport velocity        ',f10.6,/, &
         'Mean fluid speed relative to wave   ',f10.6,/, &
         'Volume flux due to waves            ',f10.6,/, &
         'Bernoulli constant                  ',f10.6)

  call write_blank
  write(6,'(A)') 'Surface elevations - crest to trough'
  do 20 i = 10, n+10, 10
     iend = min(i+9, n+10)
     write(6,'(10(1x,f9.4))') (z(j), j = i, iend)
20 continue

  call write_blank
  write(6,'(A)') 'Fourier coefficients'
  do 21 i = 1, n, 5
     iend = min(i+4, n)
     write(6,'(5(i4,1x,f10.6,3x))') (j, z(j+n+10), j = i, iend)
21 continue

  call write_blank

  pulse = z(8) + z(1)*z(5)
  ke = 0.5d0*(z(4)*pulse + z(5)*(z(8) - z(7)*z(1)))

  pe = 0.5d0*(z(10)**2 + z(n+10)**2)
  do 7 i = 1, n-1
7    pe = pe + z(10+i)**2
  pe = pe/(2.d0*n)

  ! Historical Fenton (1988) bed-velocity term, retained exactly.
  ! Under some current/reference-frame cases this value may be negative.
  ub2 = 2.d0*z(9) - z(4)**2
  sxx = 4.d0*ke - 3.d0*pe + ub2*z(1) + 2.d0*z(5)*(z(7)*z(1) - z(8))
  f = z(4)*(3.d0*ke - 2.d0*pe) + 0.5d0*ub2*(pulse + z(4)*z(1)) &
      + z(4)*z(5)*(z(7)*z(1) - z(8))

  write(6,8) pulse, ke, pe, ub2, sxx, f
8 format('Integral quantities',/, &
         'Impulse (I)                    ',e13.6,/, &
         'Kinetic energy (T)             ',e13.6,/, &
         'Potential energy (V)           ',e13.6,/, &
         'Mean squared bed velocity      ',e13.6,/, &
         'Radiation stress (Sxx)         ',e13.6,/, &
         'Wave power (F)                 ',e13.6)

  if (depth .eq. 'finite') then
     q = z(7)*z(1) - z(8)
     r = z(9) + z(1)
     s = sxx - 2.d0*z(4)*pulse + (z(4)**2 + 0.5d0*z(1))*z(1)
     call write_blank
     write(6,9) q, r, s
9    format('Invariants for finite depth',/, &
            'Volume flux (Q)               ',e13.6,/, &
            'Bernoulli constant (R)         ',e13.6,/, &
            'Momentum flux (S)             ',e13.6)
  endif

  return
end subroutine output

! -----------------------------------------------------------------------------
! SUBROUTINE FOR CALCULATION OF FREE SURFACE ELEVATION elevn AT kx,
! AND VELOCITY (u,v) AND PRESSURE press AT POINT (kx,ky).
! -----------------------------------------------------------------------------
subroutine point(kx, ky, u, v, press, elevn)
  implicit double precision(a-h,k-l,o-z)
  character*10 depth, case, currnt

  integer, parameter :: maxn = 50
  integer, parameter :: maxvar = 2*maxn + 10

  common /one/ n, num, pi, hoverd, height, value, depth, case, currnt
  common /two/ z(maxvar), cosa(0:maxvar), sina(0:maxvar), coeff(maxvar), sol(maxvar,2), y(maxvar)

  elevn = 0.5d0*y(n)*cos(n*kx)
  do 1 j = 1, n-1
1    elevn = elevn + y(j)*cos(j*kx)

  u = z(4) - z(7)
  v = 0.d0

  if (depth .eq. 'finite') then
     do 4 j = 1, n
        e = exp(j*(z(1) + ky))
        s = 0.5d0*(e - 1.d0/e)
        c = 0.5d0*(e + 1.d0/e)
        b = z(n+j+10)/cosh(j*z(1))
        u = u + j*b*c*cos(j*kx)
4       v = v + j*b*s*sin(j*kx)
  else
     do 5 j = 1, n
        e = exp(j*ky)
        u = u + j*z(n+j+10)*e*cos(j*kx)
5       v = v + j*z(n+j+10)*e*sin(j*kx)
  endif

  press = z(9) - ky - 0.5d0*((u - z(4))**2 + v*v)
  return
end subroutine point

! -----------------------------------------------------------------------------
! Minimal LINPACK-compatible LU factorisation routine.
! Interface matches DGEFA(A,LDA,N,IPVT,INFO).
! -----------------------------------------------------------------------------
subroutine dgefa(a, lda, n, ipvt, info)
  implicit double precision(a-h,o-z)
  integer lda, n, ipvt(n), info
  double precision a(lda,n), maxa, t

  info = 0
  if (n .eq. 0) return

  do 10 k = 1, n-1
     kp1 = k + 1
     l = k
     maxa = abs(a(k,k))
     do 20 i = kp1, n
        if (abs(a(i,k)) .gt. maxa) then
           maxa = abs(a(i,k))
           l = i
        endif
20   continue
     ipvt(k) = l

     if (a(l,k) .eq. 0.d0) then
        info = k
     else
        if (l .ne. k) then
           t = a(l,k)
           a(l,k) = a(k,k)
           a(k,k) = t
        endif

        do 30 i = kp1, n
30         a(i,k) = a(i,k)/a(k,k)

        do 50 j = kp1, n
           if (l .ne. k) then
              t = a(l,j)
              a(l,j) = a(k,j)
              a(k,j) = t
           endif
           do 40 i = kp1, n
40            a(i,j) = a(i,j) - a(i,k)*a(k,j)
50      continue
     endif
10 continue

  ipvt(n) = n
  if (a(n,n) .eq. 0.d0) info = n
  return
end subroutine dgefa

! -----------------------------------------------------------------------------
! Minimal LINPACK-compatible solve routine.
! Interface matches DGESL(A,LDA,N,IPVT,B,JOB). Only JOB=0 is needed here.
! -----------------------------------------------------------------------------
subroutine dgesl(a, lda, n, ipvt, b, job)
  implicit double precision(a-h,o-z)
  integer lda, n, ipvt(n), job
  double precision a(lda,n), b(n), t

  if (job .ne. 0) then
     write(6,*) 'DGESL replacement only supports JOB=0.'
     stop
  endif

  ! Solve L*y = b.
  do 10 k = 1, n-1
     l = ipvt(k)
     if (l .ne. k) then
        t = b(l)
        b(l) = b(k)
        b(k) = t
     endif
     do 20 i = k+1, n
20      b(i) = b(i) - a(i,k)*b(k)
10 continue

  ! Solve U*x = y.
  do 30 kb = 1, n
     k = n + 1 - kb
     b(k) = b(k)/a(k,k)
     t = -b(k)
     do 40 i = 1, k-1
40      b(i) = b(i) + t*a(i,k)
30 continue

  return
end subroutine dgesl
