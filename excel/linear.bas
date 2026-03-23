Attribute VB_Name = "linear"
' -----------------------------------------------------------------------
' WAVELENGTH CALCULATOR (Linear Dispersion with Current)
' =====================================================================
' PURPOSE:
'   Computes the linear Airy-theory wavelength for a regular wave in
'   finite depth, including Doppler interaction with a steady uniform
'   ambient current.
'
' PUBLIC ENTRY POINT:
'   L_linear(T, d, [u])
'
' INPUTS:
'   T : Wave Period [s]
'   d : Water Depth [m]
'   u : Current Velocity [m/s] (+ following, - opposing)
'
' RETURNS:
'   L_linear : Linear predicted wavelength [m]
'
' NUMERICAL METHOD:
'   Solves the Doppler-shifted linear dispersion relation by Newton-
'   Raphson iteration on wavenumber, then converts k to wavelength.
'
' ROLE IN WORKBOOK:
'   Baseline physical solver and reference initializer for the higher-
'   order nonlinear and surrogate wavelength modules.
' -----------------------------------------------------------------------

Option Explicit

Public Function L_linear(T As Double, d As Double, Optional U As Double = 0#) As Variant
    On Error GoTo ErrHandler

    ' --- 1. Physics Baseline ---
    Const G As Double = 9.80665
    Const PI As Double = 3.14159265359
    Const Tolerance As Double = 0.000000000001 ' 1E-12 Precision
    Const MaxIter As Integer = 100
    
    Dim omega As Double, c_shallow As Double
    omega = 2 * PI / T
    c_shallow = Sqr(G * d)
    
    ' --- 2. Initial Guess Logic (User Snippet) ---
    Dim k As Double, L0 As Double
    If U < -0.5 * c_shallow Then
        ' Guard for strong opposing currents
        k = omega / (c_shallow * 0.6)
    Else
        ' Standard deep water guess
        L0 = (G * T ^ 2) / (2 * PI)
        k = 2 * PI / L0
    End If
    
    ' --- 3. Newton-Raphson Iteration ---
    Dim i As Integer
    Dim sig As Double, th As Double, f As Double, dfdk As Double
    Dim k_new As Double
    
    For i = 1 To MaxIter
        ' Intrinsic frequency (Doppler shifted)
        sig = omega - k * U
        
        ' Hyperbolic Tangent with stability guard
        th = Tanh_S(k * d)
        
        ' Function: f(k) = g*k*tanh(kd) - (omega - kU)^2
        f = G * k * th - sig ^ 2
        
        ' Derivative: df/dk = g*tanh(kd) + g*k*d*sech^2(kd) + 2*(omega - kU)*U
        ' Note: sech^2(x) = 1 - tanh^2(x)
        dfdk = G * th + G * k * d * (1 - th ^ 2) + 2 * sig * U
        
        ' Safety check for derivative
        If Abs(dfdk) < 0.000000000000001 Then Exit For
        
        ' Update k
        k_new = k - f / dfdk
        
        ' Check Convergence
        If Abs(k_new - k) < Tolerance Then
            k = k_new
            GoTo Success
        End If
        
        k = k_new
    Next i

Success:
    If k <= 0 Then
        L_linear = "Err: No Conv"
    Else
        L_linear = 2 * PI / k
    End If
    Exit Function

ErrHandler:
    L_linear = "Err: " & err.Description
End Function

' --- Helper Math Function for Stability ---
Private Function Tanh_S(x As Double) As Double
    If x > 20 Then
        Tanh_S = 1
    ElseIf x < -20 Then
        Tanh_S = -1
    Else
        Dim ex As Double, emx As Double
        ex = Exp(x)
        emx = Exp(-x)
        Tanh_S = (ex - emx) / (ex + emx)
    End If
End Function

