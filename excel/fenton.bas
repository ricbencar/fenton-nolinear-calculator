Attribute VB_Name = "fenton"
' -----------------------------------------------------------------------
' WAVELENGTH CALCULATOR (Fenton Nonlinear Fourier Solver)
' =====================================================================
' PURPOSE:
'   Computes nonlinear wavelength for steady periodic gravity waves in
'   finite depth with optional Eulerian current, following the Fenton
'   Fourier / stream-function solution strategy.
'
' PUBLIC ENTRY POINT:
'   L_fenton(H, T, d, [u], [UseHomotopy])
'
' INPUTS:
'   H : Wave Height [m]
'   T : Wave Period [s]
'   d : Water Depth [m]
'   u : Current Velocity [m/s] (+ following, - opposing)
'
' RETURNS:
'   L_fenton : Nonlinear wavelength [m] or VBA error text when the solve
'              does not converge or the input state is invalid.
'
' NUMERICAL METHOD:
'   Uses continuation, Newton iteration, dense linear algebra, and the
'   same solver architecture adopted in the reference Python/C++ Fenton
'   implementation used elsewhere in the project.
'
' ROLE IN WORKBOOK:
'   Highest-fidelity physical wavelength model available in the Excel
'   workbook and the main nonlinear benchmark for the surrogate modules.
' -----------------------------------------------------------------------

Option Explicit

' ------------------------------ constants -------------------------------------
Private Const PI As Double = 3.14159265358979
Private Const G_STD As Double = 9.80665
Private Const N_FOURIER As Long = 20
Private Const CURRENT_CRITERION_EULERIAN As Long = 1

' ==============================================================================
'  Public API
' ==============================================================================
Public Function L_fenton(ByVal H As Double, ByVal T As Double, ByVal d As Double, Optional ByVal U As Double = 0#, Optional ByVal UseHomotopy As Boolean = True) As Variant
    ' UseHomotopy is kept for backward compatibility;
    ' the solver always uses the function.py continuation strategy.
    On Error GoTo Fail

    Dim errMsg As String
    Dim L As Double

    L = SolveFentonWavelength(H, T, d, U, errMsg)
    If Len(errMsg) <> 0 Then
        L_fenton = "Err: " & errMsg
    Else
        L_fenton = L
    End If
    Exit Function

Fail:
    L_fenton = "Err: " & err.Description
End Function

' ==============================================================================
'  Core solver (function.py parity)
' ==============================================================================
Private Function SolveFentonWavelength(ByVal H As Double, ByVal T As Double, ByVal d As Double, ByVal U As Double, ByRef errMsg As String) As Double
    errMsg = vbNullString

    If (H <= 0#) Or (T <= 0#) Or (d <= 0#) Then
        errMsg = "Invalid inputs: H, T, and d must be > 0."
        SolveFentonWavelength = 0#
        Exit Function
    End If

    Dim n As Long
    n = N_FOURIER
    Dim num As Long
    num = 2 * n + 10

    ' Unknown vector (1-based; index 0 unused).
    Dim z() As Double
    ReDim z(0 To num)
    Dim rhs1() As Double
    ReDim rhs1(0 To num)
    Dim rhs2() As Double
    ReDim rhs2(0 To num)
    Dim coeff() As Double
    ReDim coeff(0 To n)
    Dim tanhT() As Double
    ReDim tanhT(0 To n)
    Dim b() As Double
    ReDim b(0 To n)
    Dim y() As Double
    ReDim y(0 To num)

    ' Continuation storage sol[i][1..2].
    Dim sol() As Double
    ReDim sol(0 To num, 0 To 2)

    ' Trig tables as in function.py init().
    Dim cosa() As Double
    ReDim cosa(0 To 2 * n)
    Dim sina() As Double
    ReDim sina(0 To 2 * n)
    Dim cosnm() As Double
    ReDim cosnm(0 To n, 1 To n)
    Dim sinnm() As Double
    ReDim sinnm(0 To n, 1 To n)

    Call PrecomputeTrigTables(n, cosa, sina, cosnm, sinnm)

    ' Non-dimensional input groups (match function.py).
    Dim MaxH As Double
    MaxH = H / d
    Dim Tnd As Double
    Tnd = T * Sqr(G_STD / d)
    Dim Height As Double
    If Tnd > 0# Then
        Height = MaxH / (Tnd * Tnd)     ' == H / (g*T^2)
    Else
        Height = 0#
    End If
    Dim Current As Double
    Current = U / Sqr(G_STD * d)

    Dim CurrentCriterion As Long
    CurrentCriterion = CURRENT_CRITERION_EULERIAN

    ' Solver control (match function.py defaults).
    Dim nstep As Long
    nstep = 4
    Dim number As Long
    number = 40
    Dim crit As Double
    crit = 0.00000001
    Dim criterFinal As Double
    criterFinal = 0.0000000001

    ' Large-current robustness budget (equations unchanged).
    If Abs(Current) >= 1# Then
        If nstep < 8 Then nstep = 8
        If number < 80 Then number = 80
    End If

    Dim dhe As Double
    dhe = Height / CDbl(nstep)
    Dim dho As Double
    dho = MaxH / CDbl(nstep)

    Dim ns As Long, it As Long, i As Long
    Dim heightStep As Double, Hoverd As Double
    Dim err As Double, criter As Double
    Dim stepConverged As Boolean

    For ns = 1 To nstep
        heightStep = CDbl(ns) * dhe
        Hoverd = CDbl(ns) * dho

        If ns = 1 Then
            Call InitLinearState(z, sol, n, num, heightStep, Hoverd, Current, CurrentCriterion, cosa)
        Else
            ' Extrapolation: z = 2*sol(:,2) - sol(:,1)
            For i = 1 To num
                z(i) = 2# * sol(i, 2) - sol(i, 1)
            Next i

            ' Fallback to last converged state if extrapolation is invalid.
            If (Not IsFiniteVec(z, 1, num)) Or (z(1) <= 0#) Then
                For i = 1 To num
                    z(i) = sol(i, 2)
                Next i
            End If
            If (Not IsFiniteVec(z, 1, num)) Or (z(1) <= 0#) Then
                errMsg = "Invalid extrapolated start state for continuation step."
                SolveFentonWavelength = 0#
                Exit Function
            End If
        End If

        stepConverged = False
        For it = 1 To number
            err = NewtonUpdateWithRetry(z, rhs1, rhs2, coeff, tanhT, _
                                        cosnm, sinnm, sol, n, num, ns, it, _
                                        Hoverd, heightStep, Current, _
                                        CurrentCriterion)

            ' IMPORTANT: update continuation storage BEFORE convergence break.
            If ns = 1 Then
                For i = 1 To num
                    sol(i, 2) = z(i)
                Next i
            Else
                For i = 1 To num
                    sol(i, 1) = sol(i, 2)
                    sol(i, 2) = z(i)
                Next i
            End If

            If (Not IsFiniteVec(z, 1, num)) Or (z(1) <= 0#) Then
                errMsg = "Divergence: non-finite/invalid state vector encountered."
                SolveFentonWavelength = 0#
                Exit Function
            End If

            If ns = nstep Then
                criter = criterFinal
            Else
                criter = crit
            End If

            If (it > 1) And (err < criter * Abs(z(1))) Then
                stepConverged = True
                Exit For
            End If
        Next it

        If Not stepConverged Then
            errMsg = "Newton did not converge within " & CStr(number) & " iterations at continuation step " & CStr(ns) & "/" & CStr(nstep) & "."
            SolveFentonWavelength = 0#
            Exit Function
        End If

        ' Update Y and B (as in function.py / C++ reference).
        Call ComputeYandB(z, b, y, cosa, n)
    Next ns

    Dim kd As Double
    kd = z(1)
    If (kd <= 0#) Or (Not IsFinite(kd)) Then
        errMsg = "Invalid wavenumber (kd)."
        SolveFentonWavelength = 0#
        Exit Function
    End If

    SolveFentonWavelength = 2# * PI * d / kd
End Function

' ==============================================================================
'  Newton wrapper with the same first-iteration retry logic as function.py
' ==============================================================================
Private Function NewtonUpdateWithRetry(ByRef z() As Double, _
                                       ByRef rhs1() As Double, _
                                       ByRef rhs2() As Double, _
                                       ByRef coeff() As Double, _
                                       ByRef tanhT() As Double, _
                                       ByRef cosnm() As Double, _
                                       ByRef sinnm() As Double, _
                                       ByRef sol() As Double, _
                                       ByVal n As Long, _
                                       ByVal num As Long, _
                                       ByVal ns As Long, _
                                       ByVal it As Long, _
                                       ByVal Hoverd As Double, _
                                       ByVal heightStep As Double, _
                                       ByVal Current As Double, _
                                       ByVal CurrentCriterion As Long) As Double
    On Error GoTo RetryNewton

    NewtonUpdateWithRetry = NewtonUpdate(z, rhs1, rhs2, coeff, tanhT, _
                                         cosnm, sinnm, n, num, Hoverd, _
                                         heightStep, Current, _
                                         CurrentCriterion)
    Exit Function

RetryNewton:
    If (ns > 1) And (it = 1) Then
        Dim i As Long
        For i = 1 To num
            z(i) = sol(i, 2)
        Next i

        err.Clear
        NewtonUpdateWithRetry = NewtonUpdate(z, rhs1, rhs2, coeff, tanhT, _
                                             cosnm, sinnm, n, num, Hoverd, _
                                             heightStep, Current, _
                                             CurrentCriterion)
        Exit Function
    End If

    err.Raise err.number, err.Source, err.Description
End Function

' ==============================================================================
'  Initial state (function.py _init_linear)
' ==============================================================================
Private Sub InitLinearState(ByRef z() As Double, ByRef sol() As Double, ByVal n As Long, ByVal num As Long, ByVal heightStep As Double, ByVal Hoverd As Double, ByVal Current As Double, ByVal CurrentCriterion As Long, ByRef cosa() As Double)
    Dim sigma As Double
    If Hoverd > 0# Then
        sigma = 2# * PI * Sqr(heightStep / Hoverd)
    Else
        sigma = 0#
    End If

    If sigma > 0# Then
        z(1) = (sigma * sigma) / (Tanh_S(sigma ^ 1.5) ^ (2# / 3#))
    Else
        z(1) = 2# * PI * MaxD(heightStep, 0.000000000001) / MaxD(Hoverd, 0.000000000001)
    End If

    z(2) = z(1) * Hoverd
    z(4) = Sqr(Tanh_S(z(1)))
    z(3) = 2# * PI / z(4)

    If CurrentCriterion = 1 Then
        z(5) = Current * Sqr(z(2))
        z(6) = 0#
    Else
        z(6) = Current * Sqr(z(2))
        z(5) = 0#
    End If

    z(7) = z(4)
    z(8) = 0#
    z(9) = 0.5 * z(7) * z(7)

    z(10) = 0.5 * z(2)

    Dim i As Long
    For i = 1 To n
        z(n + i + 10) = 0#
        z(10 + i) = 0.5 * z(2) * cosa(i)
    Next i

    z(n + 11) = 0.5 * z(2) / z(7)

    ' Extrapolation seed (function.py): sol[1..9,1]=z; sol[10..,1]=0.
    For i = 1 To 9
        sol(i, 1) = z(i)
    Next i
    For i = 10 To num
        sol(i, 1) = 0#
    Next i
End Sub

' ==============================================================================
'  Residual equations Eqns() (finite depth, period case)
' ==============================================================================
Private Function EqnsSS(ByRef z() As Double, ByRef rhs() As Double, ByRef coeff() As Double, ByRef tanhT() As Double, ByRef cosnm() As Double, ByRef sinnm() As Double, ByVal n As Long, ByVal num As Long, ByVal Hoverd As Double, ByVal heightStep As Double, ByVal Current As Double, ByVal CurrentCriterion As Long, ByRef ok As Boolean) As Double
    ok = True

    Dim kd As Double
    kd = z(1)
    Dim i As Long

    ' Eqn 1..5
    rhs(1) = z(2) - z(1) * Hoverd
    rhs(2) = z(2) - heightStep * z(3) * z(3)
    rhs(3) = z(4) * z(3) - 2# * PI
    rhs(4) = z(5) + z(7) - z(4)
    rhs(5) = z(1) * (z(6) + z(7) - z(4)) - z(8)

    ' coeff and tanh tables
    For i = 1 To n
        coeff(i) = z(n + i + 10)
        tanhT(i) = Tanh_S(CDbl(i) * kd)
    Next i

    ' Eqn 6 (finite depth; uses sqrt(kd))
    rhs(6) = z(CurrentCriterion + 4) - Current * Sqr(kd)

    ' Eqn 7 (mean level)
    Dim rhs7 As Double
    rhs7 = z(10) + z(n + 10)
    For i = 1 To n - 1
        rhs7 = rhs7 + 2# * z(10 + i)
    Next i
    rhs(7) = rhs7

    ' Eqn 8 (wave height)
    rhs(8) = z(10) - z(n + 10) - z(2)

    ' Eqns 9..
    Dim m As Long, jj As Long
    Dim zsurf As Double, x As Double, e As Double, invE As Double
    Dim sinhkd As Double, coshkd As Double, s As Double, c As Double
    Dim psi As Double, U As Double, v As Double, jcj As Double
    Dim ccos As Double, ssin As Double, tj As Double

    For m = 0 To n
        zsurf = z(10 + m)

        psi = 0#
        U = 0#
        v = 0#

        For jj = 1 To n
            x = CDbl(jj) * zsurf
            If (x > 60#) Or (x < -60#) Then
                ok = False
                EqnsSS = 1E+308
                Exit Function
            End If

            e = Exp(x)
            invE = 1# / e
            sinhkd = 0.5 * (e - invE)
            coshkd = 0.5 * (e + invE)

            tj = tanhT(jj)
            s = sinhkd + coshkd * tj
            c = coshkd + sinhkd * tj

            ccos = cosnm(m, jj)
            ssin = sinnm(m, jj)

            psi = psi + coeff(jj) * s * ccos
            jcj = CDbl(jj) * coeff(jj)
            U = U + jcj * c * ccos
            v = v + jcj * s * ssin
        Next jj

        rhs(m + 9) = psi - z(8) - z(7) * z(10 + m)
        rhs(n + m + 10) = 0.5 * ((-z(7) + U) ^ 2 + v ^ 2) + z(10 + m) - z(9)
    Next m

    ' Sum of squares
    Dim ss As Double
    ss = 0#
    For i = 1 To num
        ss = ss + rhs(i) * rhs(i)
    Next i
    EqnsSS = ss
End Function

' ==============================================================================
'  Newton update (finite-difference Jacobian + SVD solve + backtracking)
' ==============================================================================
Private Function NewtonUpdate(ByRef z() As Double, ByRef rhs1() As Double, ByRef rhs2() As Double, ByRef coeff() As Double, ByRef tanhT() As Double, ByRef cosnm() As Double, ByRef sinnm() As Double, ByVal n As Long, ByVal num As Long, ByVal Hoverd As Double, ByVal heightStep As Double, ByVal Current As Double, ByVal CurrentCriterion As Long) As Double
    Dim ok As Boolean
    Dim ss0 As Double
    ss0 = EqnsSS(z, rhs1, coeff, tanhT, cosnm, sinnm, n, num, Hoverd, heightStep, Current, CurrentCriterion, ok)
    If Not ok Then err.Raise vbObjectError + 7001, , "Non-finite residual."

    Dim z0() As Double
    ReDim z0(0 To num)
    Dim i As Long
    For i = 1 To num
        z0(i) = z(i)
    Next i

    Dim a() As Double
    ReDim a(1 To num, 1 To num)
    Dim b() As Double
    ReDim b(1 To num)

    Dim H As Double, r As Long
    For i = 1 To num
        H = 0.01 * z0(i)
        If Abs(z0(i)) < 0.0001 Then H = 0.00001
        If Abs(H) > 1# Then H = SignD(1#, H)

        z(i) = z0(i) + H
        Call EqnsSS(z, rhs2, coeff, tanhT, cosnm, sinnm, n, num, Hoverd, heightStep, Current, CurrentCriterion, ok)
        z(i) = z0(i)
        If Not ok Then err.Raise vbObjectError + 7002, , "Divergence in Jacobian FD."

        b(i) = -rhs1(i)
        For r = 1 To num
            a(r, i) = (rhs2(r) - rhs1(r)) / H
        Next r
    Next i

    Dim dx() As Double
    dx = SvdSolveSquare(a, b, num)
    If Not IsFiniteVec(dx, 1, num) Then
        err.Raise vbObjectError + 7003, , "Non-finite Newton correction vector."
    End If

    Dim alpha As Double
    alpha = 1#
    Dim ssBest As Double
    ssBest = ss0
    Dim zBest() As Double
    ReDim zBest(0 To num)
    Dim zTry() As Double
    ReDim zTry(0 To num)

    For i = 1 To num
        zBest(i) = z0(i)
    Next i

    Dim ss1 As Double
    Do While alpha >= 0.0001
        For i = 1 To num
            zTry(i) = z0(i) + alpha * dx(i)
        Next i

        If (zTry(1) <= 0#) Or (Not IsFiniteVec(zTry, 1, num)) Then
            alpha = alpha * 0.5
            GoTo ContinueLineSearch
        End If

        For i = 1 To num
            z(i) = zTry(i)
        Next i

        ss1 = EqnsSS(z, rhs2, coeff, tanhT, cosnm, sinnm, n, num, Hoverd, heightStep, Current, CurrentCriterion, ok)

        If ok And IsFinite(ss1) And (ss1 <= ssBest) Then
            ssBest = ss1
            For i = 1 To num
                zBest(i) = zTry(i)
            Next i

            If ss1 <= ss0 Then Exit Do
        End If

        alpha = alpha * 0.5
ContinueLineSearch:
    Loop

    ' Commit best found
    For i = 1 To num
        z(i) = zBest(i)
    Next i

    ' err = mean(abs(zBest[10..n+10] - z0[10..n+10]))
    Dim sumAbs As Double
    sumAbs = 0#
    Dim cnt As Long
    cnt = 0
    For i = 10 To (n + 10)
        sumAbs = sumAbs + Abs(zBest(i) - z0(i))
        cnt = cnt + 1
    Next i
    If cnt > 0 Then
        NewtonUpdate = sumAbs / CDbl(cnt)
    Else
        NewtonUpdate = 0#
    End If
End Function

' ==============================================================================
'  Post-convergence Fourier transform block (Compute Y and B)
' ==============================================================================
Private Sub ComputeYandB(ByRef z() As Double, ByRef b() As Double, ByRef y() As Double, ByRef cosa() As Double, ByVal n As Long)
    Dim j As Long, m As Long
    Dim twoN As Long
    twoN = 2 * n

    For j = 0 To UBound(y)
        y(j) = 0#
    Next j

    For j = 1 To n
        b(j) = z(j + n + 10)

        Dim sign As Double
        If (j Mod 2) = 1 Then
            sign = -1#
        Else
            sign = 1#
        End If

        Dim s As Double
        s = 0.5 * (z(10) + z(n + 10) * sign)
        For m = 1 To n - 1
            s = s + z(10 + m) * cosa((m * j) Mod twoN)
        Next m
        y(j) = 2# * s / CDbl(n)
    Next j
End Sub

' ==============================================================================
'  Trig precomputation (match function.py init())
' ==============================================================================
Private Sub PrecomputeTrigTables(ByVal n As Long, ByRef cosa() As Double, ByRef sina() As Double, ByRef cosnm() As Double, ByRef sinnm() As Double)
    Dim k As Long
    For k = 0 To 2 * n
        cosa(k) = Cos(CDbl(k) * PI / CDbl(n))
        sina(k) = Sin(CDbl(k) * PI / CDbl(n))
    Next k

    Dim m As Long, j As Long, idx As Long
    For m = 0 To n
        For j = 1 To n
            idx = (m * j) Mod (2 * n)
            cosnm(m, j) = cosa(idx)
            sinnm(m, j) = sina(idx)
        Next j
    Next m
End Sub

' ==============================================================================
'  Linear solve: SVD pseudo-inverse (Press-style truncation)
' ==============================================================================
Private Function SvdSolveSquare(ByRef a() As Double, ByRef b() As Double, ByVal n As Long) As Double()
    Dim U() As Double
    ReDim U(0 To n, 0 To n)
    Dim v() As Double
    ReDim v(0 To n, 0 To n)
    Dim w() As Double
    ReDim w(0 To n)
    Dim x() As Double
    ReDim x(0 To n)

    Dim i As Long, j As Long
    For i = 1 To n
        For j = 1 To n
            U(i, j) = a(i, j)
        Next j
    Next i

    Call SVDCMP(U, n, n, w, v)

    Dim wMax As Double
    wMax = 0#
    For i = 1 To n
        If w(i) > wMax Then wMax = w(i)
    Next i
    Dim wMin As Double
    wMin = wMax * 0.000000000001
    For i = 1 To n
        If w(i) <= wMin Then w(i) = 0#
    Next i

    Call SVBKSB(U, w, v, n, n, b, x)
    SvdSolveSquare = x
End Function

' ==============================================================================
'  SVD routines (Numerical Recipes style, adapted for VBA Double arrays)
' ==============================================================================
Private Sub SVBKSB(ByRef U() As Double, ByRef w() As Double, ByRef v() As Double, ByVal m As Long, ByVal n As Long, ByRef b() As Double, ByRef x() As Double)
    Dim jj As Long, j As Long, i As Long
    Dim s As Double
    Dim tmp() As Double
    ReDim tmp(0 To n)

    For j = 1 To n
        s = 0#
        If w(j) <> 0# Then
            For i = 1 To m
                s = s + U(i, j) * b(i)
            Next i
            s = s / w(j)
        End If
        tmp(j) = s
    Next j

    For j = 1 To n
        s = 0#
        For jj = 1 To n
            s = s + v(j, jj) * tmp(jj)
        Next jj
        x(j) = s
    Next j
End Sub

Private Sub SVDCMP(ByRef a() As Double, ByVal m As Long, ByVal n As Long, ByRef w() As Double, ByRef v() As Double)
    Dim flag As Long, i As Long, its As Long, j As Long, jj As Long
    Dim k As Long, L As Long, nm As Long
    Dim anorm As Double
    Dim c As Double
    Dim f As Double
    Dim G As Double
    Dim H As Double
    Dim s As Double
    Dim sc As Double
    Dim x As Double
    Dim y As Double
    Dim z As Double
    Dim rv1() As Double
    ReDim rv1(1 To n)

    G = 0#
    sc = 0#
    anorm = 0#

    For i = 1 To n
        L = i + 1
        rv1(i) = sc * G
        G = 0#
        s = 0#
        sc = 0#

        If i <= m Then
            For k = i To m
                sc = sc + Abs(a(k, i))
            Next k
            If sc <> 0# Then
                For k = i To m
                    a(k, i) = a(k, i) / sc
                    s = s + a(k, i) * a(k, i)
                Next k
                f = a(i, i)
                G = -SignD(Sqr(s), f)
                H = f * G - s
                a(i, i) = f - G
                If i <> n Then
                    For j = L To n
                        s = 0#
                        For k = i To m
                            s = s + a(k, i) * a(k, j)
                        Next k
                        f = s / H
                        For k = i To m
                            a(k, j) = a(k, j) + f * a(k, i)
                        Next k
                    Next j
                End If
                For k = i To m
                    a(k, i) = a(k, i) * sc
                Next k
            End If
        End If

        w(i) = sc * G
        G = 0#
        s = 0#
        sc = 0#

        If (i <= m) And (i <> n) Then
            For k = L To n
                sc = sc + Abs(a(i, k))
            Next k
            If sc <> 0# Then
                For k = L To n
                    a(i, k) = a(i, k) / sc
                    s = s + a(i, k) * a(i, k)
                Next k
                f = a(i, L)
                G = -SignD(Sqr(s), f)
                H = f * G - s
                a(i, L) = f - G
                For k = L To n
                    rv1(k) = a(i, k) / H
                Next k
                If i <> m Then
                    For j = L To m
                        s = 0#
                        For k = L To n
                            s = s + a(j, k) * a(i, k)
                        Next k
                        For k = L To n
                            a(j, k) = a(j, k) + s * rv1(k)
                        Next k
                    Next j
                End If
                For k = L To n
                    a(i, k) = a(i, k) * sc
                Next k
            End If
        End If

        anorm = MaxD(anorm, Abs(w(i)) + Abs(rv1(i)))
    Next i

    ' Accumulation of right-hand transformations.
    For i = n To 1 Step -1
        If i < n Then
            If G <> 0# Then
                For j = L To n
                    v(j, i) = (a(i, j) / a(i, L)) / G
                Next j
                For j = L To n
                    s = 0#
                    For k = L To n
                        s = s + a(i, k) * v(k, j)
                    Next k
                    For k = L To n
                        v(k, j) = v(k, j) + s * v(k, i)
                    Next k
                Next j
            End If
            For j = L To n
                v(i, j) = 0#
                v(j, i) = 0#
            Next j
        End If
        v(i, i) = 1#
        G = rv1(i)
        L = i
    Next i

    ' Accumulation of left-hand transformations.
    Dim mn As Long
    mn = IIf(m < n, m, n)
    For i = mn To 1 Step -1
        L = i + 1
        G = w(i)
        If i < n Then
            For j = L To n
                a(i, j) = 0#
            Next j
        End If
        If G <> 0# Then
            G = 1# / G
            If i <> n Then
                For j = L To n
                    s = 0#
                    For k = L To m
                        s = s + a(k, i) * a(k, j)
                    Next k
                    f = (s / a(i, i)) * G
                    For k = i To m
                        a(k, j) = a(k, j) + f * a(k, i)
                    Next k
                Next j
            End If
            For j = i To m
                a(j, i) = a(j, i) * G
            Next j
        Else
            For j = i To m
                a(j, i) = 0#
            Next j
        End If
        a(i, i) = a(i, i) + 1#
    Next i

    ' Diagonalisation of the bidiagonal form.
    For k = n To 1 Step -1
        For its = 1 To 30
            flag = 1
            For L = k To 1 Step -1
                nm = L - 1
                If (Abs(rv1(L)) + anorm) = anorm Then
                    flag = 0
                    Exit For
                End If
                If (Abs(w(nm)) + anorm) = anorm Then Exit For
            Next L

            If flag <> 0 Then
                c = 0#
                s = 1#
                For i = L To k
                    f = s * rv1(i)
                    rv1(i) = c * rv1(i)
                    If (Abs(f) + anorm) = anorm Then Exit For
                    G = w(i)
                    H = Pythag(f, G)
                    w(i) = H
                    H = 1# / H
                    c = G * H
                    s = -f * H
                    For j = 1 To m
                        y = a(j, nm)
                        z = a(j, i)
                        a(j, nm) = y * c + z * s
                        a(j, i) = z * c - y * s
                    Next j
                Next i
            End If

            z = w(k)
            If L = k Then
                If z < 0# Then
                    w(k) = -z
                    For j = 1 To n
                        v(j, k) = -v(j, k)
                    Next j
                End If
                Exit For
            End If

            If its = 30 Then err.Raise vbObjectError + 7101, , "SVD no converge."

            x = w(L)
            nm = k - 1
            y = w(nm)
            G = rv1(nm)
            H = rv1(k)
            f = ((y - z) * (y + z) + (G - H) * (G + H)) / (2# * H * y)
            G = Pythag(f, 1#)
            f = ((x - z) * (x + z) + H * (y / (f + SignD(G, f)) - H)) / x

            c = 1#
            s = 1#
            For j = L To nm
                i = j + 1
                G = rv1(i)
                y = w(i)
                H = s * G
                G = c * G
                z = Pythag(f, H)
                rv1(j) = z
                c = f / z
                s = H / z
                f = x * c + G * s
                G = G * c - x * s
                H = y * s
                y = y * c

                For jj = 1 To n
                    x = v(jj, j)
                    z = v(jj, i)
                    v(jj, j) = x * c + z * s
                    v(jj, i) = z * c - x * s
                Next jj

                z = Pythag(f, H)
                w(j) = z
                If z <> 0# Then
                    z = 1# / z
                    c = f * z
                    s = H * z
                End If
                f = c * G + s * y
                x = c * y - s * G

                For jj = 1 To m
                    y = a(jj, j)
                    z = a(jj, i)
                    a(jj, j) = y * c + z * s
                    a(jj, i) = z * c - y * s
                Next jj
            Next j
            rv1(L) = 0#
            rv1(k) = f
            w(k) = x
        Next its
    Next k
End Sub

' ==============================================================================
'  Numerics helpers
' ==============================================================================
Private Function Pythag(ByVal a As Double, ByVal b As Double) As Double
    Dim absa As Double
    absa = Abs(a)
    Dim absb As Double
    absb = Abs(b)
    If absa > absb Then
        Pythag = absa * Sqr(1# + (absb / absa) ^ 2)
    ElseIf absb = 0# Then
        Pythag = 0#
    Else
        Pythag = absb * Sqr(1# + (absa / absb) ^ 2)
    End If
End Function

Private Function SignD(ByVal a As Double, ByVal b As Double) As Double
    If b >= 0# Then
        SignD = Abs(a)
    Else
        SignD = -Abs(a)
    End If
End Function

Private Function MaxD(ByVal a As Double, ByVal b As Double) As Double
    If a >= b Then
        MaxD = a
    Else
        MaxD = b
    End If
End Function

Private Function IsFinite(ByVal x As Double) As Boolean
    ' True for normal finite doubles; False for NaN/Inf/overflowed values.
    IsFinite = (x = x) And (Abs(x) < 1E+308)
End Function

Private Function IsFiniteVec(ByRef v() As Double, ByVal lo As Long, ByVal hi As Long) As Boolean
    Dim i As Long
    For i = lo To hi
        If Not IsFinite(v(i)) Then
            IsFiniteVec = False
            Exit Function
        End If
    Next i
    IsFiniteVec = True
End Function

Private Function Tanh_S(ByVal x As Double) As Double
    ' Stable tanh approximation for VBA (matches typical libm behaviour).
    If x > 20# Then
        Tanh_S = 1#
    ElseIf x < -20# Then
        Tanh_S = -1#
    Else
        Dim e2 As Double
        e2 = Exp(2# * x)
        Tanh_S = (e2 - 1#) / (e2 + 1#)
    End If
End Function

