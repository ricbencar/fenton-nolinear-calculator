Attribute VB_Name = "fenton"
Option Explicit

' Finite-depth steady nonlinear waves with prescribed Eulerian current.
' L_fenton(H, T, d, [U]) returns wavelength in metres.
' Inputs use metres, seconds and metres per second; g = 9.80665 m/s^2.
' A failed, inadmissible or under-resolved state returns an "Err: " message.
' Acceptance checks the finite residual and selected outputs; it is not a Lean program proof.
' The Fourier order is limited by the column-equilibrated Jacobian condition.

Private Type FentonInterval
    lo As Double
    hi As Double
End Type

Private Const RHO As Double = 1025#
Private Const PI As Double = 3.1415926535897931
Private Const G_STD As Double = 9.80665
Private Const N_FOURIER As Long = 50
Private Const MIN_FOURIER As Long = 12
Private Const FINAL_TOL As Double = 1E-10
Private Const CONDITION_LIMIT As Double = 1E+11
Private Const FAILURE_VALUE As Double = 1E+300
Private Const MAX_ITERATIONS As Long = 80
Private Const MAX_TRIALS As Long = 160
Private Const MIN_INCREMENT As Double = 0.000244140625

Public Function L_fenton(ByVal H As Double, ByVal T As Double, ByVal d As Double, _
                         Optional ByVal U As Double = 0#) As Variant
    On Error GoTo Failed
    Dim message As String, wavelength As Double
    Dim diagnostics() As Double, values() As Double, bounds() As Double
    wavelength = SolveFentonWavelength(H, T, d, U, message, diagnostics, values, bounds)
    If Len(message) = 0 Then
        L_fenton = wavelength
    Else
        L_fenton = "Err: " & message
    End If
    Exit Function
Failed:
    L_fenton = "Err: " & Err.Description
End Function

Private Function SolveFentonWavelength(ByVal H As Double, ByVal T As Double, _
                                      ByVal d As Double, ByVal U As Double, _
                                      ByRef message As String, ByRef diagnostics() As Double, _
                                      ByRef values() As Double, ByRef bounds() As Double) As Double
    message = vbNullString
    ReDim diagnostics(1 To 11): ReDim values(1 To 9)
    ReDim bounds(1 To 1, 1 To 4)
    If Not IsFinite(H) Or Not IsFinite(T) Or Not IsFinite(d) Or Not IsFinite(U) Then
        message = "Inputs must be finite."
        Exit Function
    End If
    If H <= 0# Or T <= 0# Or d <= 0# Then
        message = "H, T and d must be positive."
        Exit Function
    End If

    Dim relativeHeight As Double, heightParameter As Double, current As Double
    relativeHeight = H / d
    heightParameter = H / G_STD / T / T
    current = U / Sqr(G_STD * d)
    If Not IsFinite(relativeHeight) Or Not IsFinite(heightParameter) Then
        message = "Invalid dimensionless input groups."
        Exit Function
    End If
    If relativeHeight <= 0# Or heightParameter <= 0# Then
        message = "Dimensionless height is outside the numerical range."
        Exit Function
    End If

    Dim n As Long, num As Long, solved As Boolean
    Dim condition As Double, residual As Double
    Dim z() As Double, cosnm() As Double, sinnm() As Double
    n = N_FOURIER
    Do
        num = 2 * n + 10
        ReDim z(0 To num)
        ReDim cosnm(0 To n, 1 To n)
        ReDim sinnm(0 To n, 1 To n)
        PrecomputeTrigTables n, cosnm, sinnm
        solved = ContinueToTarget(z, cosnm, sinnm, n, num, relativeHeight, _
                                  heightParameter, current, residual)
        condition = StateCondition(z, cosnm, sinnm, n, num, relativeHeight, _
                                   heightParameter, current)
        If condition <= CONDITION_LIMIT Or n <= MIN_FOURIER Then Exit Do
        n = CLng(Int(0.8 * n))
        If n < MIN_FOURIER Then n = MIN_FOURIER
    Loop
    diagnostics(1) = n: diagnostics(2) = condition: diagnostics(3) = residual
    If Not solved Then
        message = "No converged solution at the requested height; N=" & CStr(n) & _
                  ", max|F|=" & Format$(residual, "0.000E+00") & "."
        Exit Function
    End If
    If condition > CONDITION_LIMIT Then
        message = "The Fourier system exceeds the numerical conditioning limit."
        Exit Function
    End If
    Dim kPhysical As Double, cPhysical As Double
    kPhysical = z(1) / d: cPhysical = 2# * PI / kPhysical / T
    z(3) = T * Sqr(G_STD * kPhysical)
    z(4) = cPhysical * Sqr(kPhysical / G_STD)
    z(5) = U * Sqr(kPhysical / G_STD)
    z(7) = z(4) - z(5)
    z(6) = z(8) / z(1) - z(7) + z(4)
    Dim rhs() As Double, merit As Double, ok As Boolean
    ReDim rhs(0 To num)
    merit = EqnsSS(z, rhs, cosnm, sinnm, n, num, relativeHeight, heightParameter, current, ok)
    If Not ok Then
        message = "Non-finite final residual."
        Exit Function
    End If
    diagnostics(3) = MaxResidual(rhs, num)
    If diagnostics(3) > FINAL_TOL Then
        message = "The final residual exceeds tolerance."
        Exit Function
    End If
    If Not VerifyAccuracy(z, n, H, T, d, U, values, diagnostics, bounds, message) Then Exit Function
    SolveFentonWavelength = values(1)
End Function

Private Function ContinueToTarget(ByRef z() As Double, ByRef cosnm() As Double, _
                                  ByRef sinnm() As Double, ByVal n As Long, ByVal num As Long, _
                                  ByVal relativeHeight As Double, ByVal heightParameter As Double, _
                                  ByVal current As Double, ByRef residual As Double) As Boolean
    On Error GoTo Failed
    Dim accepted() As Double, previous() As Double, rhs() As Double
    ReDim accepted(0 To num)
    ReDim previous(0 To num)
    ReDim rhs(0 To num)
    Dim fraction As Double, previousFraction As Double, increment As Double
    Dim target As Double, ratio As Double, trialHeight As Double, trialParameter As Double
    Dim trial As Long, iteration As Long, i As Long
    Dim havePrevious As Boolean, converged As Boolean, ok As Boolean, merit As Double
    increment = 0.25
    If Abs(current) >= 1# Then increment = 0.125
    residual = FAILURE_VALUE
    For trial = 1 To MAX_TRIALS
        target = fraction + increment
        If target > 1# Then target = 1#
        trialHeight = target * relativeHeight
        trialParameter = target * heightParameter
        If fraction = 0# Then
            InitLinearState z, n, num, trialParameter, trialHeight, current, cosnm
        Else
            ratio = 0#
            If havePrevious Then ratio = (target - fraction) / (fraction - previousFraction)
            For i = 1 To num
                z(i) = accepted(i) + ratio * (accepted(i) - previous(i))
            Next i
            If Not Admissible(z, n, num) Then CopyState accepted, z, num
        End If
        converged = False
        For iteration = 0 To MAX_ITERATIONS
            merit = EqnsSS(z, rhs, cosnm, sinnm, n, num, trialHeight, trialParameter, current, ok)
            residual = FAILURE_VALUE
            If ok Then residual = MaxResidual(rhs, num)
            If residual <= FINAL_TOL Then
                converged = True
                Exit For
            End If
            If Not ok Or iteration = MAX_ITERATIONS Then Exit For
            If Not NewtonUpdate(z, cosnm, sinnm, n, num, trialHeight, trialParameter, current) Then Exit For
        Next iteration
        If converged Then
            CopyState accepted, previous, num
            previousFraction = fraction
            havePrevious = (fraction > 0#)
            CopyState z, accepted, num
            fraction = target
            If fraction = 1# Then
                ContinueToTarget = True
                Exit Function
            End If
            If iteration <= 8 Then increment = increment * 1.5
            If increment > 0.25 Then increment = 0.25
        Else
            If fraction > 0# Then CopyState accepted, z, num
            increment = increment * 0.5
            If increment < MIN_INCREMENT Then Exit Function
        End If
    Next trial
    Exit Function
Failed:
    If fraction > 0# Then CopyState accepted, z, num
    residual = FAILURE_VALUE
End Function

Private Sub CopyState(ByRef source() As Double, ByRef destination() As Double, ByVal num As Long)
    Dim i As Long
    For i = 1 To num
        destination(i) = source(i)
    Next i
End Sub

Private Sub PrecomputeTrigTables(ByVal n As Long, ByRef cosnm() As Double, ByRef sinnm() As Double)
    Dim m As Long, j As Long, phase As Double
    For m = 0 To n
        For j = 1 To n
            phase = CDbl((m * j) Mod (2 * n)) * PI / n
            cosnm(m, j) = Cos(phase)
            sinnm(m, j) = Sin(phase)
        Next j
    Next m
End Sub

Private Sub InitLinearState(ByRef z() As Double, ByVal n As Long, ByVal num As Long, _
                            ByVal heightParameter As Double, ByVal relativeHeight As Double, _
                            ByVal current As Double, ByRef cosnm() As Double)
    Dim sigma As Double, i As Long
    For i = 1 To num
        z(i) = 0#
    Next i
    sigma = 2# * PI * Sqr(heightParameter / relativeHeight)
    z(1) = sigma * sigma / (Tanh_S(sigma ^ 1.5) ^ (2# / 3#))
    z(2) = z(1) * relativeHeight
    z(7) = Sqr(Tanh_S(z(1)))
    z(5) = current * Sqr(z(1))
    z(6) = z(5)
    z(4) = z(7) + z(5)
    If z(4) <= 0# Then z(4) = z(7)
    z(3) = 2# * PI / z(4)
    z(9) = 0.5 * z(7) * z(7)
    For i = 0 To n
        z(10 + i) = 0.5 * z(2) * cosnm(i, 1)
    Next i
    z(n + 11) = 0.5 * z(2) / z(7)
End Sub

Private Function ExpSafe(ByVal exponent As Double) As Double
    If exponent < -745# Then
        ExpSafe = 0#
    Else
        ExpSafe = Exp(exponent)
    End If
End Function

Private Sub VerticalBasis(ByVal mode As Long, ByVal eta As Double, ByVal kd As Double, _
                          ByRef stream As Double, ByRef velocity As Double)
    Dim forward As Double, reflected As Double, denominator As Double
    forward = ExpSafe(mode * eta)
    reflected = ExpSafe(-mode * (2# * kd + eta))
    denominator = 1# + ExpSafe(-2# * mode * kd)
    stream = (forward - reflected) / denominator
    velocity = (forward + reflected) / denominator
End Sub

Private Function Admissible(ByRef z() As Double, ByVal n As Long, ByVal num As Long) As Boolean
    If Not IsFiniteVec(z, 1, num) Then Exit Function
    If z(1) <= 0# Or z(2) <= 0# Or z(3) <= 0# Or z(4) <= 0# Or z(7) <= 0# Then Exit Function
    Dim m As Long
    For m = 0 To n
        If z(1) + z(10 + m) <= 0# Then Exit Function
        If n * Abs(z(10 + m)) > 600# Then Exit Function
    Next m
    Admissible = True
End Function

Private Function EqnsSS(ByRef z() As Double, ByRef rhs() As Double, _
                        ByRef cosnm() As Double, ByRef sinnm() As Double, _
                        ByVal n As Long, ByVal num As Long, ByVal relativeHeight As Double, _
                        ByVal heightParameter As Double, ByVal current As Double, ByRef ok As Boolean) As Double
    On Error GoTo Failed
    ok = False
    EqnsSS = FAILURE_VALUE
    If Not Admissible(z, n, num) Then Exit Function
    rhs(1) = z(2) - z(1) * relativeHeight
    rhs(2) = z(2) - heightParameter * z(3) * z(3)
    rhs(3) = z(4) * z(3) - 2# * PI
    rhs(4) = z(5) + z(7) - z(4)
    rhs(5) = z(1) * (z(6) + z(7) - z(4)) - z(8)
    rhs(6) = z(5) - current * Sqr(z(1))
    rhs(7) = z(10) + z(n + 10)
    Dim i As Long, m As Long, j As Long
    For i = 1 To n - 1
        rhs(7) = rhs(7) + 2# * z(10 + i)
    Next i
    rhs(8) = z(10) - z(n + 10) - z(2)
    Dim psi As Double, uSum As Double, vSum As Double, stream As Double, velocity As Double
    Dim coefficient As Double, eta As Double, total As Double
    For m = 0 To n
        psi = 0#: uSum = 0#: vSum = 0#
        eta = z(10 + m)
        For j = 1 To n
            VerticalBasis j, eta, z(1), stream, velocity
            coefficient = z(n + 10 + j)
            psi = psi + coefficient * stream * cosnm(m, j)
            uSum = uSum + j * coefficient * velocity * cosnm(m, j)
            vSum = vSum + j * coefficient * stream * sinnm(m, j)
        Next j
        rhs(m + 9) = psi - z(8) - z(7) * eta
        rhs(n + m + 10) = 0.5 * ((uSum - z(7)) ^ 2 + vSum ^ 2) + eta - z(9)
    Next m
    For i = 1 To num
        total = total + rhs(i) * rhs(i)
    Next i
    If Not IsFinite(total) Then Exit Function
    EqnsSS = total
    ok = True
    Exit Function
Failed:
    ok = False
    EqnsSS = FAILURE_VALUE
End Function

Private Function MaxResidual(ByRef rhs() As Double, ByVal num As Long) As Double
    Dim i As Long, maximum As Double
    For i = 1 To num
        If Not IsFinite(rhs(i)) Then
            MaxResidual = FAILURE_VALUE
            Exit Function
        End If
        If Abs(rhs(i)) > maximum Then maximum = Abs(rhs(i))
    Next i
    MaxResidual = maximum
End Function

' Analytic derivatives of the global constraints and surface equations.
Private Sub AssembleJacobian(ByRef z() As Double, ByRef jac() As Double, _
                             ByRef cosnm() As Double, ByRef sinnm() As Double, _
                             ByVal n As Long, ByVal num As Long, ByVal relativeHeight As Double, _
                             ByVal heightParameter As Double, ByVal current As Double)
    ReDim jac(1 To num, 1 To num)
    jac(1, 1) = -relativeHeight: jac(1, 2) = 1#
    jac(2, 2) = 1#: jac(2, 3) = -2# * heightParameter * z(3)
    jac(3, 3) = z(4): jac(3, 4) = z(3)
    jac(4, 5) = 1#: jac(4, 7) = 1#: jac(4, 4) = -1#
    jac(5, 1) = z(6) + z(7) - z(4)
    jac(5, 6) = z(1): jac(5, 7) = z(1): jac(5, 4) = -z(1): jac(5, 8) = -1#
    jac(6, 5) = 1#: jac(6, 1) = -0.5 * current / Sqr(z(1))
    jac(7, 10) = 1#: jac(7, n + 10) = 1#
    Dim m As Long, j As Long, row As Long
    For m = 1 To n - 1
        jac(7, 10 + m) = 2#
    Next m
    jac(8, 10) = 1#: jac(8, n + 10) = -1#: jac(8, 2) = -1#
    Dim eta As Double, stream As Double, velocity As Double, coefficient As Double
    Dim uSum As Double, vSum As Double, uySum As Double, vySum As Double
    Dim pd As Double, ud As Double, vd As Double, e As Double, forward As Double
    Dim reflected As Double, factor As Double, sd As Double, cd As Double, relativeU As Double
    For m = 0 To n
        eta = z(10 + m)
        uSum = 0#: vSum = 0#: uySum = 0#: vySum = 0#
        pd = 0#: ud = 0#: vd = 0#
        For j = 1 To n
            VerticalBasis j, eta, z(1), stream, velocity
            coefficient = z(n + 10 + j)
            uSum = uSum + j * coefficient * velocity * cosnm(m, j)
            vSum = vSum + j * coefficient * stream * sinnm(m, j)
            uySum = uySum + CDbl(j) * j * coefficient * stream * cosnm(m, j)
            vySum = vySum + CDbl(j) * j * coefficient * velocity * sinnm(m, j)
            e = ExpSafe(-2# * j * z(1))
            forward = ExpSafe(j * eta)
            reflected = ExpSafe(-j * (2# * z(1) + eta))
            factor = 2# * j / ((1# + e) * (1# + e))
            sd = factor * (reflected + e * forward)
            cd = factor * (e * forward - reflected)
            pd = pd + coefficient * sd * cosnm(m, j)
            ud = ud + j * coefficient * cd * cosnm(m, j)
            vd = vd + j * coefficient * sd * sinnm(m, j)
            jac(m + 9, n + 10 + j) = stream * cosnm(m, j)
        Next j
        relativeU = uSum - z(7)
        jac(m + 9, 1) = pd
        jac(m + 9, 7) = -eta: jac(m + 9, 8) = -1#
        jac(m + 9, 10 + m) = relativeU
        row = n + m + 10
        jac(row, 1) = relativeU * ud + vSum * vd
        jac(row, 7) = -relativeU: jac(row, 9) = -1#
        jac(row, 10 + m) = 1# + relativeU * uySum + vSum * vySum
        For j = 1 To n
            VerticalBasis j, eta, z(1), stream, velocity
            jac(row, n + 10 + j) = j * (relativeU * velocity * cosnm(m, j) + vSum * stream * sinnm(m, j))
        Next j
    Next m
End Sub

Private Function NewtonUpdate(ByRef z() As Double, ByRef cosnm() As Double, ByRef sinnm() As Double, _
                              ByVal n As Long, ByVal num As Long, ByVal relativeHeight As Double, _
                              ByVal heightParameter As Double, ByVal current As Double) As Boolean
    On Error GoTo Failed
    Dim z0() As Double, rhs() As Double, b() As Double, jac() As Double, dx() As Double
    ReDim z0(0 To num): ReDim rhs(0 To num): ReDim b(1 To num)
    CopyState z, z0, num
    Dim baseReady As Boolean
    baseReady = True
    Dim ok As Boolean, ss0 As Double, ss1 As Double, condition As Double, alpha As Double, i As Long
    ss0 = EqnsSS(z, rhs, cosnm, sinnm, n, num, relativeHeight, heightParameter, current, ok)
    If Not ok Then Exit Function
    For i = 1 To num
        b(i) = -rhs(i)
    Next i
    AssembleJacobian z, jac, cosnm, sinnm, n, num, relativeHeight, heightParameter, current
    dx = SvdSolveSquare(jac, b, num, condition)
    If Not IsFiniteVec(dx, 1, num) Then Exit Function
    alpha = 1#
    Do While alpha >= 0.00000095367431640625
        For i = 1 To num
            z(i) = z0(i) + alpha * dx(i)
        Next i
        ss1 = EqnsSS(z, rhs, cosnm, sinnm, n, num, relativeHeight, heightParameter, current, ok)
        If ok Then
            If ss1 < ss0 * (1# - 0.0001 * alpha) Or MaxResidual(rhs, num) <= FINAL_TOL Then
                NewtonUpdate = True
                Exit Function
            End If
        End If
        alpha = alpha * 0.5
    Loop
    CopyState z0, z, num
    Exit Function
Failed:
    If baseReady Then CopyState z0, z, num
End Function

Private Function StateCondition(ByRef z() As Double, ByRef cosnm() As Double, ByRef sinnm() As Double, _
                                ByVal n As Long, ByVal num As Long, ByVal relativeHeight As Double, _
                                ByVal heightParameter As Double, ByVal current As Double) As Double
    On Error GoTo Failed
    StateCondition = FAILURE_VALUE
    If Not Admissible(z, n, num) Then Exit Function
    Dim jac() As Double, b() As Double, dx() As Double, condition As Double
    ReDim b(1 To num)
    ' Use the height constraints associated with the retained state.
    relativeHeight = z(2) / z(1)
    heightParameter = z(2) / z(3) / z(3)
    AssembleJacobian z, jac, cosnm, sinnm, n, num, relativeHeight, heightParameter, current
    dx = SvdSolveSquare(jac, b, num, condition)
    StateCondition = condition
    Exit Function
Failed:
    StateCondition = FAILURE_VALUE
End Function

Private Function SvdSolveSquare(ByRef a() As Double, ByRef b() As Double, ByVal n As Long, _
                                ByRef condition As Double) As Double()
    Dim uMat() As Double, vMat() As Double, w() As Double, x() As Double, scales() As Double
    ReDim uMat(0 To n, 0 To n): ReDim vMat(0 To n, 0 To n)
    ReDim w(0 To n): ReDim x(0 To n): ReDim scales(1 To n)
    Dim i As Long, j As Long, maximum As Double, wMax As Double, wMin As Double, smallest As Double
    For j = 1 To n
        maximum = 0#
        For i = 1 To n
            If Not IsFinite(a(i, j)) Then Err.Raise vbObjectError + 7002, , "Non-finite Jacobian."
            If Abs(a(i, j)) > maximum Then maximum = Abs(a(i, j))
        Next i
        If maximum = 0# Then maximum = 1#
        scales(j) = maximum
        For i = 1 To n
            uMat(i, j) = a(i, j) / maximum
        Next i
    Next j
    SVDCMP uMat, n, n, w, vMat
    smallest = FAILURE_VALUE
    For i = 1 To n
        If Not IsFinite(w(i)) Then Err.Raise vbObjectError + 7003, , "Non-finite singular value."
        If w(i) > wMax Then wMax = w(i)
        If w(i) < smallest Then smallest = w(i)
    Next i
    condition = FAILURE_VALUE
    If smallest > wMax * 1E-300 And smallest > 0# Then condition = wMax / smallest
    wMin = wMax * 1E-15
    For i = 1 To n
        If w(i) <= wMin Then w(i) = 0#
    Next i
    SVBKSB uMat, w, vMat, n, n, b, x
    For i = 1 To n
        x(i) = x(i) / scales(i)
    Next i
    SvdSolveSquare = x
End Function

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
        For its = 1 To 100
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

            If its = 100 Then err.Raise vbObjectError + 7101, , "SVD did not converge."

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
    ' Saturated branches avoid overflow in the exponential.
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


Private Sub SurfaceCoefficients(ByRef z() As Double, ByVal n As Long, ByRef coeff() As Double)
    ReDim coeff(0 To n)
    Dim j As Long, m As Long, total As Double, phase As Double
    For j = 1 To n
        total = z(10) + z(n + 10) * Cos(j * PI)
        For m = 1 To n - 1
            phase = CDbl((m * j) Mod (2 * n)) * PI / n
            total = total + 2# * z(10 + m) * Cos(phase)
        Next m
        coeff(j) = total / n
    Next j
End Sub

Private Function SurfaceOrdinate(ByRef coeff() As Double, ByVal n As Long, ByVal phase As Double) As Double
    Dim j As Long, total As Double
    For j = 1 To n - 1
        total = total + coeff(j) * Cos(j * phase)
    Next j
    SurfaceOrdinate = total + 0.5 * coeff(n) * Cos(n * phase)
End Function

Private Sub SelectedOutputs(ByRef z() As Double, ByRef coeff() As Double, ByVal n As Long, _
                            ByVal d As Double, ByRef values() As Double, ByRef bedAbsolute As Double)
    ReDim values(1 To 9)
    Dim k As Double, speed As Double, energyScale As Double, pulse As Double, flux As Double
    Dim kinetic As Double, potential As Double, beta As Double, bed As Double, inverseCosh As Double
    Dim j As Long
    k = z(1) / d: speed = Sqr(G_STD / k): energyScale = RHO * G_STD / (k * k)
    pulse = z(8) + z(1) * z(5): flux = z(7) * z(1) - z(8)
    kinetic = 0.5 * (z(4) * pulse - z(5) * flux)
    For j = 1 To n - 1
        potential = potential + 0.25 * coeff(j) * coeff(j)
    Next j
    potential = potential + coeff(n) * coeff(n) / 16#
    beta = 2# * z(9) - z(4) * z(4) + 2# * z(4) * z(5)
    For j = 1 To n
        inverseCosh = 2# * ExpSafe(-j * z(1)) / (1# + ExpSafe(-2# * j * z(1)))
        bed = bed + 0.5 * (j * z(n + 10 + j) * inverseCosh) ^ 2
    Next j
    values(1) = 2# * PI / k
    values(2) = z(4) * speed
    values(3) = RHO * speed / k * pulse
    values(4) = energyScale * kinetic
    values(5) = energyScale * potential
    values(6) = energyScale * (kinetic + potential)
    values(7) = energyScale * (4# * kinetic - 3# * potential + beta * z(1) - 2# * z(5) * pulse)
    values(8) = energyScale * speed * (z(4) * (3# * kinetic - 2# * potential) _
                + 0.5 * beta * (pulse + z(4) * z(1)) - 2# * z(4) * z(5) * pulse)
    values(9) = G_STD / k * bed
    bedAbsolute = values(9) + z(5) * z(5) * G_STD / k
    If Not IsFiniteVec(values, 1, 9) Then Err.Raise vbObjectError + 7210, , "Non-finite dimensional output."
End Sub

Private Function SurfaceResidual(ByRef z() As Double, ByRef coeff() As Double, ByVal n As Long) As Double
    Dim count As Long, m As Long, j As Long, phase As Double, eta As Double
    Dim stream As Double, velocity As Double, psi As Double, uSum As Double, vSum As Double
    Dim co As Double, si As Double, amplitude As Double, maximum As Double
    count = CLng(MaxD(2048#, 32# * n))
    For m = 0 To count - 1
        phase = m * 2# * PI / count
        eta = SurfaceOrdinate(coeff, n, phase)
        psi = 0#: uSum = 0#: vSum = 0#
        For j = 1 To n
            VerticalBasis j, eta, z(1), stream, velocity
            co = Cos(j * phase): si = Sin(j * phase): amplitude = z(n + 10 + j)
            psi = psi + amplitude * stream * co
            uSum = uSum + j * amplitude * velocity * co
            vSum = vSum + j * amplitude * stream * si
        Next j
        maximum = MaxD(maximum, Abs(psi - z(8) - z(7) * eta))
        maximum = MaxD(maximum, Abs(0.5 * ((uSum - z(7)) ^ 2 + vSum ^ 2) + eta - z(9)))
    Next m
    If Not IsFinite(maximum) Then Err.Raise vbObjectError + 7211, , "Non-finite surface residual."
    SurfaceResidual = maximum
End Function

Private Sub LegendreValue(ByVal n As Long, ByVal x As Double, ByRef value As Double, ByRef derivative As Double)
    Dim j As Long, p0 As Double, p1 As Double, p2 As Double
    p0 = 1#: p1 = x
    For j = 2 To n
        p2 = ((2# * j - 1#) * x * p1 - (j - 1#) * p0) / j
        p0 = p1: p1 = p2
    Next j
    value = p1
    derivative = n * (x * p1 - p0) / (x * x - 1#)
End Sub

Private Sub GaussLegendre(ByVal n As Long, ByRef nodes() As Double, ByRef weights() As Double)
    ReDim nodes(0 To n - 1): ReDim weights(0 To n - 1)
    Dim i As Long, iteration As Long, x As Double, nextX As Double, p As Double, derivative As Double
    Dim converged As Boolean, weight As Double
    For i = 0 To CLng(Int((n - 1) / 2#))
        x = Cos(PI * (i + 0.75) / (n + 0.5))
        converged = False
        For iteration = 1 To 100
            LegendreValue n, x, p, derivative
            nextX = x - p / derivative
            If Abs(nextX - x) <= 8.8817841970012523E-16 Then
                x = nextX: converged = True
                Exit For
            End If
            x = nextX
        Next iteration
        If Not converged Then Err.Raise vbObjectError + 7212, , "Gauss-Legendre nodes did not converge."
        LegendreValue n, x, p, derivative
        weight = 2# / ((1# - x * x) * derivative * derivative)
        nodes(i) = -x: nodes(n - 1 - i) = x
        weights(i) = weight: weights(n - 1 - i) = weight
    Next i
End Sub

Private Sub DefiningIntegrals(ByRef z() As Double, ByRef coeff() As Double, ByVal n As Long, _
                             ByVal d As Double, ByVal phaseCount As Long, ByVal verticalCount As Long, _
                             ByRef totals() As Double)
    ReDim totals(1 To 5)
    Dim nodes() As Double, weights() As Double, co() As Double, si() As Double, amplitude() As Double
    ReDim co(1 To n): ReDim si(1 To n): ReDim amplitude(1 To n)
    GaussLegendre verticalCount, nodes, weights
    Dim k As Double, speed As Double, c As Double, current As Double, bernoulli As Double
    Dim m As Long, j As Long, v As Long, phase As Double, eta As Double, thickness As Double
    Dim y As Double, weight As Double, stream As Double, velocity As Double
    Dim uSum As Double, vSum As Double, pressure As Double, kinetic As Double
    k = z(1) / d: speed = Sqr(G_STD / k)
    c = z(4) * speed: current = z(5) * speed: bernoulli = z(9) * G_STD / k
    For j = 1 To n
        amplitude(j) = j * z(n + 10 + j) * speed
    Next j
    For m = 0 To phaseCount - 1
        phase = 2# * PI * m / phaseCount
        eta = SurfaceOrdinate(coeff, n, phase) / k
        thickness = 0.5 * (d + eta)
        If thickness <= 0# Or Not IsFinite(thickness) Then Err.Raise vbObjectError + 7213, , "The interpolated surface intersects the bed."
        For j = 1 To n
            co(j) = Cos(j * phase): si(j) = Sin(j * phase)
        Next j
        totals(3) = totals(3) + 0.5 * RHO * G_STD * eta * eta / phaseCount
        For v = 0 To verticalCount - 1
            y = -d + thickness * (nodes(v) + 1#)
            weight = weights(v) * thickness / phaseCount
            uSum = current: vSum = 0#
            For j = 1 To n
                VerticalBasis j, k * y, z(1), stream, velocity
                uSum = uSum + amplitude(j) * velocity * co(j)
                vSum = vSum + amplitude(j) * stream * si(j)
            Next j
            pressure = RHO * (bernoulli - G_STD * y - 0.5 * ((uSum - c) ^ 2 + vSum ^ 2))
            kinetic = 0.5 * RHO * (uSum * uSum + vSum * vSum)
            totals(1) = totals(1) + weight * RHO * uSum
            totals(2) = totals(2) + weight * kinetic
            totals(4) = totals(4) + weight * (pressure + RHO * uSum * uSum)
            totals(5) = totals(5) + weight * (pressure + RHO * G_STD * y + kinetic) * uSum
        Next v
    Next m
    totals(4) = totals(4) - 0.5 * RHO * G_STD * d * d
    If Not IsFiniteVec(totals, 1, 5) Then Err.Raise vbObjectError + 7214, , "Non-finite defining integral."
End Sub

Private Function OutputBudget(ByVal reference As Double) As Double
    OutputBudget = 1E-7 + 1E-8 * Abs(reference)
End Function

Private Function VerifyAccuracy(ByRef z() As Double, ByVal n As Long, ByVal H As Double, _
                                ByVal T As Double, ByVal d As Double, ByVal current As Double, _
                                ByRef values() As Double, ByRef diagnostics() As Double, _
                                ByRef bounds() As Double, ByRef message As String) As Boolean
    Dim coeff() As Double, coarse() As Double, fine() As Double, indices() As Double
    Dim bedAbsolute As Double, c As Double, eulerian As Double, expected As Double, ratio As Double
    Dim j As Long, index As Long, failures As String
    For j = 4 To 11
        diagnostics(j) = 0#
    Next j
    SurfaceCoefficients z, n, coeff
    SelectedOutputs z, coeff, n, d, values, bedAbsolute
    diagnostics(4) = SurfaceResidual(z, coeff, n)
    If diagnostics(4) > 1E-6 Then failures = failures & "interpolated-surface residual exceeds 1e-6; "
    c = values(2): eulerian = z(5) * Sqr(G_STD * d / z(1))
    expected = c * (3# * values(4) - 2# * values(5)) + 0.5 * bedAbsolute * (values(3) + RHO * c * d) _
               - 2# * c * eulerian * values(3)
    diagnostics(5) = Abs(values(8) - expected) / OutputBudget(expected)
    expected = 4# * values(4) - 3# * values(5) + RHO * bedAbsolute * d - 2# * eulerian * values(3)
    diagnostics(5) = MaxD(diagnostics(5), Abs(values(7) - expected) / OutputBudget(expected))
    If diagnostics(5) > 1# Then failures = failures & "spectral-bed comparison exceeds its budget; "
    DefiningIntegrals z, coeff, n, d, CLng(MaxD(256#, 8# * n)), 64, coarse
    DefiningIntegrals z, coeff, n, d, CLng(MaxD(512#, 16# * n)), 96, fine
    ReDim indices(1 To 5)
    indices(1) = 3: indices(2) = 4: indices(3) = 5: indices(4) = 7: indices(5) = 8
    For j = 1 To 5
        index = CLng(indices(j))
        ratio = Abs(fine(j) - coarse(j)) / OutputBudget(fine(j))
        diagnostics(6) = MaxD(diagnostics(6), ratio)
        ratio = Abs(values(index) - fine(j)) / OutputBudget(fine(j))
        diagnostics(7) = MaxD(diagnostics(7), ratio)
    Next j
    If diagnostics(6) > 1# Or diagnostics(7) > 1# Then failures = failures & "defining-integral comparison exceeds its budget; "
    If Not VerifyFiniteSystem(z, coeff, n, H, T, d, current, values, diagnostics, bounds) Then
        failures = failures & "finite-system interval verification exceeds its budget; "
    End If
    If Len(failures) > 0 Then
        message = "The collocation system converged, but its Fourier resolution or output accuracy does not meet the budgets: " & failures
        Exit Function
    End If
    diagnostics(11) = 1#
    VerifyAccuracy = True
End Function

' Finite interval endpoints are enlarged after every arithmetic operation.
' The margin exceeds binary64 rounding, including subnormal results.
Private Function RoundDown(ByVal x As Double) As Double
    If Not IsFinite(x) Then Err.Raise vbObjectError + 7200, , "Non-finite interval endpoint."
    RoundDown = (x - Abs(x) * 8.8817841970012523E-16) - 2.2250738585072014E-308
End Function

Private Function RoundUp(ByVal x As Double) As Double
    RoundUp = -RoundDown(-x)
End Function

Private Function IBounds(ByVal lo As Double, ByVal hi As Double) As FentonInterval
    If Not IsFinite(lo) Or Not IsFinite(hi) Or lo > hi Then Err.Raise vbObjectError + 7201, , "Invalid interval bounds."
    Dim result As FentonInterval
    result.lo = lo
    result.hi = hi
    IBounds = result
End Function

Private Function IP(ByVal x As Double) As FentonInterval
    IP = IBounds(x, x)
End Function

Private Function IAdd(ByRef a As FentonInterval, ByRef b As FentonInterval) As FentonInterval
    IAdd = IBounds(RoundDown(a.lo + b.lo), RoundUp(a.hi + b.hi))
End Function

Private Function INeg(ByRef a As FentonInterval) As FentonInterval
    INeg = IBounds(-a.hi, -a.lo)
End Function

Private Function ISub(ByRef a As FentonInterval, ByRef b As FentonInterval) As FentonInterval
    ISub = IAdd(a, INeg(b))
End Function

Private Function IMul(ByRef a As FentonInterval, ByRef b As FentonInterval) As FentonInterval
    Dim p1 As Double, p2 As Double, p3 As Double, p4 As Double
    p1 = a.lo * b.lo: p2 = a.lo * b.hi
    p3 = a.hi * b.lo: p4 = a.hi * b.hi
    IMul = IBounds(RoundDown(MinD(MinD(p1, p2), MinD(p3, p4))), _
                    RoundUp(MaxD(MaxD(p1, p2), MaxD(p3, p4))))
End Function

Private Function IDiv(ByRef a As FentonInterval, ByRef b As FentonInterval) As FentonInterval
    If b.lo <= 0# And b.hi >= 0# Then Err.Raise vbObjectError + 7202, , "Interval divisor contains zero."
    IDiv = IMul(a, IBounds(RoundDown(1# / b.hi), RoundUp(1# / b.lo)))
End Function

Private Function IMagnitude(ByRef a As FentonInterval) As Double
    IMagnitude = MaxD(Abs(a.lo), Abs(a.hi))
End Function

Private Function ISquare(ByRef a As FentonInterval) As FentonInterval
    Dim lower As Double, upper As Double
    If a.lo > 0# Or a.hi < 0# Then lower = MinD(Abs(a.lo), Abs(a.hi))
    upper = IMagnitude(a)
    ISquare = IBounds(MaxD(0#, RoundDown(lower * lower)), RoundUp(upper * upper))
End Function

Private Function IPower(ByRef a As FentonInterval, ByVal exponent As Long) As FentonInterval
    Dim result As FentonInterval, factor As FentonInterval
    result = IP(1#): factor = a
    Do While exponent > 0
        If exponent Mod 2 = 1 Then result = IMul(result, factor)
        exponent = CLng(Int(exponent / 2#))
        If exponent > 0 Then factor = ISquare(factor)
    Loop
    IPower = result
End Function

Private Function IExpPoint(ByVal x As Double) As FentonInterval
    If x = 0# Then
        IExpPoint = IP(1#)
        Exit Function
    End If
    Dim reduced As Double, halvings As Long, j As Long
    Dim r As FentonInterval, total As FentonInterval, term As FentonInterval, remainder As FentonInterval
    reduced = Abs(x)
    Do While reduced > 0.125
        reduced = reduced * 0.5
        halvings = halvings + 1
    Loop
    r = IBounds(MaxD(0#, RoundDown(reduced)), RoundUp(reduced))
    total = IP(1#): term = IP(1#)
    For j = 1 To 32
        term = IDiv(IMul(term, r), IP(CDbl(j)))
        total = IAdd(total, term)
    Next j
    ' For 0 <= r <= 1/4, exp(r) < 2 bounds the Taylor remainder.
    remainder = IMul(IPower(r, 33), IP(2#))
    For j = 2 To 33
        remainder = IDiv(remainder, IP(CDbl(j)))
    Next j
    total = IBounds(total.lo, RoundUp(total.hi + remainder.hi))
    If x < 0# Then total = IDiv(IP(1#), total)
    For j = 1 To halvings
        total = ISquare(total)
    Next j
    IExpPoint = total
End Function

Private Function IExp(ByRef a As FentonInterval) As FentonInterval
    Dim lower As FentonInterval, upper As FentonInterval
    lower = IExpPoint(a.lo): upper = IExpPoint(a.hi)
    IExp = IBounds(lower.lo, upper.hi)
End Function

Private Function ISqrt(ByRef a As FentonInterval) As FentonInterval
    If a.lo < 0# Then Err.Raise vbObjectError + 7203, , "Negative square-root interval."
    If a.hi = 0# Then
        ISqrt = IP(0#)
        Exit Function
    End If
    Dim lower As Double, upper As Double, j As Long
    lower = MaxD(0#, RoundDown(Sqr(a.lo)))
    upper = RoundUp(Sqr(a.hi))
    For j = 1 To 16
        If lower = 0# Or RoundUp(lower * lower) <= a.lo Then Exit For
        lower = MaxD(0#, RoundDown(lower))
    Next j
    For j = 1 To 16
        If RoundDown(upper * upper) >= a.hi Then Exit For
        upper = RoundUp(upper)
    Next j
    If lower > 0# And RoundUp(lower * lower) > a.lo Then Err.Raise vbObjectError + 7204, , "Square-root enclosure failed."
    If RoundDown(upper * upper) < a.hi Then Err.Raise vbObjectError + 7204, , "Square-root enclosure failed."
    ISqrt = IBounds(lower, upper)
End Function

Private Function IPi() As FentonInterval
    IPi = IBounds(RoundDown(PI), RoundUp(PI))
End Function

Private Sub IntervalTrig(ByVal n As Long, ByRef sine() As FentonInterval, ByRef cosine() As FentonInterval)
    ReDim sine(0 To 2 * n - 1): ReDim cosine(0 To 2 * n - 1)
    Dim index As Long, q As Long, m As Long, j As Long
    Dim x As FentonInterval, xx As FentonInterval, st As FentonInterval, ct As FentonInterval
    Dim se As FentonInterval, ce As FentonInterval
    For index = 0 To 2 * n - 1
        q = index
        If index > n Then q = index - 2 * n
        x = IDiv(IMul(IPi(), IP(CDbl(q))), IP(CDbl(n)))
        xx = ISquare(x): st = x: ct = IP(1#)
        sine(index) = IP(0#): cosine(index) = IP(0#)
        For m = 0 To 31
            sine(index) = IAdd(sine(index), st)
            cosine(index) = IAdd(cosine(index), ct)
            st = IDiv(IMul(INeg(st), xx), IP(CDbl((2 * m + 2) * (2 * m + 3))))
            ct = IDiv(IMul(INeg(ct), xx), IP(CDbl((2 * m + 1) * (2 * m + 2))))
        Next m
        se = IPower(IP(IMagnitude(x)), 65)
        ce = IPower(IP(IMagnitude(x)), 64)
        For j = 2 To 65
            se = IDiv(se, IP(CDbl(j)))
        Next j
        For j = 2 To 64
            ce = IDiv(ce, IP(CDbl(j)))
        Next j
        sine(index) = IAdd(sine(index), IBounds(-se.hi, se.hi))
        cosine(index) = IAdd(cosine(index), IBounds(-ce.hi, ce.hi))
    Next index
End Sub

Private Function MinD(ByVal a As Double, ByVal b As Double) As Double
    If a <= b Then
        MinD = a
    Else
        MinD = b
    End If
End Function

Private Function VerifyFiniteSystem(ByRef state() As Double, ByRef coeff() As Double, ByVal n As Long, _
                                    ByVal H As Double, ByVal T As Double, ByVal d As Double, ByVal current As Double, _
                                    ByRef values() As Double, ByRef diagnostics() As Double, ByRef bounds() As Double) As Boolean
    Dim z() As FentonInterval, y() As FentonInterval, f() As FentonInterval
    Dim sine() As FentonInterval, cosine() As FentonInterval, denominator() As FentonInterval
    Dim enclosed() As FentonInterval
    ReDim z(0 To 2 * n + 10): ReDim y(0 To n): ReDim f(1 To 2 * n + 10)
    ReDim denominator(1 To n): ReDim enclosed(1 To 9)
    ReDim bounds(1 To 2 * n + 19, 1 To 4)
    Dim ih As FentonInterval, it As FentonInterval, depth As FentonInterval, uc As FentonInterval
    Dim gravity As FentonInterval, density As FentonInterval, piBound As FentonInterval
    Dim middle As FentonInterval, ordinate As FentonInterval, psi As FentonInterval
    Dim us As FentonInterval, vs As FentonInterval, ep As FentonInterval, em As FentonInterval
    Dim stream As FentonInterval, velocity As FentonInterval, coefficient As FentonInterval
    Dim sn As FentonInterval, cs As FentonInterval, k As FentonInterval, speed As FentonInterval
    Dim energyScale As FentonInterval, pulse As FentonInterval, flux As FentonInterval
    Dim kinetic As FentonInterval, potential As FentonInterval, beta As FentonInterval
    Dim stress As FentonInterval, power As FentonInterval, bed As FentonInterval, inverseCosh As FentonInterval
    Dim j As Long, m As Long, index As Long, row As Long
    Dim maximum As Double, errorBound As Double, budget As Double, minimum As Double, passed As Boolean
    diagnostics(8) = 0#: diagnostics(9) = 0#: diagnostics(10) = 0#
    ih = IP(H): it = IP(T): depth = IP(d): uc = IP(current)
    gravity = IP(G_STD): density = IP(RHO): piBound = IPi()
    For j = 0 To 2 * n + 10
        z(j) = IP(state(j))
    Next j
    For j = 0 To n
        y(j) = IP(coeff(j))
    Next j
    f(1) = ISub(z(2), IDiv(IMul(z(1), ih), depth))
    f(2) = ISub(z(2), IMul(IDiv(ih, IMul(gravity, ISquare(it))), ISquare(z(3))))
    f(3) = ISub(IMul(z(4), z(3)), IMul(IP(2#), piBound))
    f(4) = ISub(IAdd(z(5), z(7)), z(4))
    f(5) = ISub(IMul(z(1), ISub(IAdd(z(6), z(7)), z(4))), z(8))
    f(6) = ISub(z(5), IMul(IDiv(uc, ISqrt(IMul(gravity, depth))), ISqrt(z(1))))
    middle = IP(0#)
    For j = 11 To n + 9
        middle = IAdd(middle, z(j))
    Next j
    f(7) = IAdd(IAdd(z(10), z(n + 10)), IMul(IP(2#), middle))
    f(8) = ISub(ISub(z(10), z(n + 10)), z(2))
    IntervalTrig n, sine, cosine
    For j = 1 To n
        denominator(j) = IAdd(IP(1#), IExp(IMul(IMul(INeg(IP(2#)), IP(CDbl(j))), z(1))))
    Next j
    For m = 0 To n
        ordinate = z(m + 10): psi = IP(0#): us = IP(0#): vs = IP(0#)
        For j = 1 To n
            ep = IExp(IMul(IP(CDbl(j)), ordinate))
            em = IExp(ISub(IMul(INeg(IP(CDbl(j))), ordinate), IMul(IMul(IP(2#), IP(CDbl(j))), z(1))))
            stream = IDiv(ISub(ep, em), denominator(j))
            velocity = IDiv(IAdd(ep, em), denominator(j))
            index = (m * j) Mod (2 * n)
            sn = sine(index): cs = cosine(index): coefficient = z(n + 10 + j)
            psi = IAdd(psi, IMul(IMul(coefficient, stream), cs))
            us = IAdd(us, IMul(IMul(IMul(IP(CDbl(j)), coefficient), velocity), cs))
            vs = IAdd(vs, IMul(IMul(IMul(IP(CDbl(j)), coefficient), stream), sn))
        Next j
        f(9 + m) = ISub(ISub(psi, z(8)), IMul(z(7), ordinate))
        f(n + 10 + m) = ISub(IAdd(IDiv(IAdd(ISquare(ISub(us, z(7))), ISquare(vs)), IP(2#)), ordinate), z(9))
    Next m
    k = IDiv(z(1), depth)
    speed = ISqrt(IDiv(gravity, k))
    energyScale = IDiv(IMul(density, gravity), ISquare(k))
    pulse = IAdd(z(8), IMul(z(1), z(5)))
    flux = ISub(IMul(z(7), z(1)), z(8))
    kinetic = IDiv(ISub(IMul(z(4), pulse), IMul(z(5), flux)), IP(2#))
    potential = IP(0#)
    For j = 1 To n - 1
        potential = IAdd(potential, IDiv(ISquare(y(j)), IP(4#)))
    Next j
    potential = IAdd(potential, IDiv(ISquare(y(n)), IP(16#)))
    beta = IAdd(ISub(IMul(IP(2#), z(9)), ISquare(z(4))), IMul(IMul(IP(2#), z(4)), z(5)))
    stress = ISub(IAdd(ISub(IMul(IP(4#), kinetic), IMul(IP(3#), potential)), IMul(beta, z(1))), IMul(IMul(IP(2#), z(5)), pulse))
    power = ISub(IAdd(IMul(z(4), ISub(IMul(IP(3#), kinetic), IMul(IP(2#), potential))), IDiv(IMul(beta, IAdd(pulse, IMul(z(4), z(1)))), IP(2#))), IMul(IMul(IMul(IP(2#), z(4)), z(5)), pulse))
    bed = IP(0#)
    For j = 1 To n
        inverseCosh = IDiv(IMul(IP(2#), IExp(IMul(INeg(IP(CDbl(j))), z(1)))), denominator(j))
        bed = IAdd(bed, IDiv(ISquare(IMul(IMul(IP(CDbl(j)), z(n + 10 + j)), inverseCosh)), IP(2#)))
    Next j
    enclosed(1) = IDiv(IMul(IP(2#), piBound), k)
    enclosed(2) = IMul(z(4), speed)
    enclosed(3) = IMul(IDiv(IMul(density, speed), k), pulse)
    enclosed(4) = IMul(energyScale, kinetic)
    enclosed(5) = IMul(energyScale, potential)
    enclosed(6) = IMul(energyScale, IAdd(kinetic, potential))
    enclosed(7) = IMul(energyScale, stress)
    enclosed(8) = IMul(IMul(energyScale, speed), power)
    enclosed(9) = IMul(IDiv(gravity, k), bed)
    For j = 1 To 2 * n + 10
        bounds(j, 1) = f(j).lo: bounds(j, 2) = f(j).hi
        maximum = MaxD(maximum, IMagnitude(f(j)))
    Next j
    diagnostics(8) = maximum
    passed = (maximum <= RoundDown(FINAL_TOL))
    For j = 1 To 9
        errorBound = IMagnitude(ISub(IP(values(j)), enclosed(j)))
        minimum = 0#
        If enclosed(j).lo > 0# Or enclosed(j).hi < 0# Then minimum = MinD(Abs(enclosed(j).lo), Abs(enclosed(j).hi))
        budget = RoundDown(RoundDown(1E-7) + RoundDown(RoundDown(1E-8) * minimum))
        diagnostics(9) = MaxD(diagnostics(9), RoundUp(errorBound / budget))
        If errorBound > budget Then passed = False
        row = 2 * n + 10 + j
        bounds(row, 1) = enclosed(j).lo: bounds(row, 2) = enclosed(j).hi
        bounds(row, 3) = errorBound: bounds(row, 4) = budget
    Next j
    If passed Then diagnostics(10) = 1#
    VerifyFiniteSystem = passed
End Function
