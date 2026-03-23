Attribute VB_Name = "effective"
' -----------------------------------------------------------------------
' WAVELENGTH CALCULATOR (Effective Period / Effective Depth Model)
' =====================================================================
' PURPOSE:
'   Computes wavelength by first transforming the physical inputs into an
'   effective period and an effective depth, then solving the linear
'   dispersion relation with those modified variables.
'
' PUBLIC ENTRY POINT:
'   L_effective(H, T, d, U)
'
' INPUTS:
'   H : Wave Height [m]
'   T : Wave Period [s]
'   d : Water Depth [m]
'   U : Current Velocity [m/s] (+ following, - opposing)
'
' RETURNS:
'   L_effective : Predicted wavelength [m]
'
' MODEL STRUCTURE:
'   T_eff(H,T,U,d) = ((U^2) - H / 1819 / 708) * (-93 / 943) + U * (500 / 739) + T
'   d_eff(H,T,U,d) = Sqr(H / Sqr(d / H)) + U * (986 / 873) * tanh(H) + d
'   L              = Linear_Dispersion_Relation(T_eff, d_eff)
'
' ROLE IN WORKBOOK:
'   Compact explicit alternative model that preserves a physical linear
'   dispersion solve while embedding empirically evolved corrections.
' -----------------------------------------------------------------------

Public Function L_effective(H As Double, T As Double, d As Double, U As Double) As Double

    Const G As Double = 9.80665
    Const PI As Double = 3.14159265359

    If d <= 0# Then
        L_effective = 0#
        Exit Function
    End If

    If H <= 0# Then H = 0.000001

    Dim T_eff As Double
    Dim d_eff As Double

    ' Exact evolved formula, with internal variable mapping H, T, U, d
    T_eff = (((U ^ 2) - H / 1819# / 708#) * (-93# / 943#)) + U * (500# / 739#) + T
    d_eff = Sqr(H / Sqr(d / H)) + U * (986# / 873#) * FnTanh(H) + d

    If Abs(T_eff) < 0.0001 Then
        L_effective = 0#
        Exit Function
    End If

    Dim T_abs As Double
    Dim d_abs As Double

    T_abs = Abs(T_eff)
    d_abs = Abs(d_eff)

    If T_abs < 0.1 Then T_abs = 0.1
    If d_abs < 0.01 Then d_abs = 0.01

    Dim omega As Double
    Dim k As Double
    Dim kd As Double
    Dim th As Double
    Dim f As Double
    Dim df As Double
    Dim diff As Double
    Dim j As Integer

    omega = (2# * PI) / T_abs
    k = (omega * omega) / G
    If k < 0.000001 Then k = 0.000001

    For j = 1 To 5
        If k > 500# Then k = 500#

        kd = k * d_abs
        th = FnTanh(kd)

        f = G * k * th - (omega * omega)
        df = G * th + G * kd * (1# - th * th)

        If Abs(df) < 0.000000001 Then Exit For

        diff = f / df
        k = k - diff

        If Abs(diff) < 0.00001 Then Exit For
    Next j

    If k <= 0# Then
        L_effective = 0#
    Else
        L_effective = (2# * PI) / k
    End If

End Function

Private Function FnTanh(x As Double) As Double
    If x > 20# Then
        FnTanh = 1#
    ElseIf x < -20# Then
        FnTanh = -1#
    Else
        Dim ex2 As Double
        ex2 = Exp(2# * x)
        FnTanh = (ex2 - 1#) / (ex2 + 1#)
    End If
End Function

