Attribute VB_Name = "genetic"
' -----------------------------------------------------------------------
' WAVELENGTH CALCULATOR (Genetic Symbolic-Regression Surrogate)
' =====================================================================
' PURPOSE:
'   Computes wavelength from an explicit gene-expression-programming
'   formula built as a correction to a linear physics-based baseline.
'
' PUBLIC ENTRY POINT:
'   L_genetic(H, T, d, u)
'
' INPUTS:
'   H : Wave Height [m]
'   T : Wave Period [s]
'   d : Water Depth [m]
'   u : Current Velocity [m/s] (+ following, - opposing)
'
' RETURNS:
'   L_genetic : Predicted wavelength [m]
'
' MODEL STRUCTURE:
'   1. Solve the linear baseline wavelength.
'   2. Construct evolved dimensionless features.
'   3. Apply the Hall-of-Fame symbolic-regression multiplier.
'
' ROLE IN WORKBOOK:
'   Transparent compact surrogate that preserves an explicit engineering
'   formula while still capturing nonlinear and current-driven effects.
' -----------------------------------------------------------------------

Public Function L_genetic(H As Double, T As Double, d As Double, U As Double) As Double
    
    ' 1. CONSTANTS
    ' ---------------------------------------------------------
    Const G As Double = 9.80665
    Const PI As Double = 3.14159265359
    
    ' Safety checks to prevent errors
    If T <= 0 Or d <= 0 Then
        L_genetic = 0
        Exit Function
    End If
    If H <= 0 Then H = 0.000001 ' Avoid Log(0) errors
    
    ' 2. CALCULATE LINEAR BASELINE (L_linear)
    ' ---------------------------------------------------------
    ' Solves the implicit equation: L = (gT^2 / 2pi) * tanh(2pi * d / L)
    ' using fixed-point iteration.
    
    Dim L_linear As Double
    Dim L0 As Double
    Dim k As Double
    Dim i As Integer
    
    ' Deep water guess (L0)
    L0 = (G * T ^ 2) / (2 * PI)
    L_linear = L0
    
    ' Iterative Solver (matches C++ Newton-Raphson logic)
    ' Converges L_linear for the dispersion relation without current
    For i = 1 To 20
        Dim k_w As Double
        k_w = (2 * PI) / L_linear
        Dim next_L As Double
        next_L = L0 * FnTanh(k_w * d)
        
        ' Check convergence
        If Abs(next_L - L_linear) < 1 Then
            L_linear = next_L
            Exit For
        End If
        L_linear = next_L
    Next i
    
    ' 3. CALCULATE DIMENSIONLESS PARAMETERS (x0 - x8)
    ' ---------------------------------------------------------
    ' These variables map to the features defined in the C++ code.
    
    ' x2 = ln(H/d) : Relative Height (Logarithmic)
    Dim x2 As Double
    x2 = Log(H / d)
    
    ' x5 = (U * T) / L_linear : Doppler Proxy
    Dim x5 As Double
    x5 = (U * T) / L_linear
    
    ' x6 = U / C0 : Velocity Ratio
    ' C0 (Deep water phase speed) = L0 / T = gT / 2pi
    Dim C0 As Double
    C0 = (G * T) / (2 * PI)
    Dim x6 As Double
    x6 = U / C0
    
    ' x8 = ln(T * sqrt(g/d)) : Dimensionless Period (Logarithmic)
    Dim x8 As Double
    x8 = Log(T * Sqr(G / d))
    
    ' 4. APPLY GEP CORRECTION FORMULA
    ' ---------------------------------------------------------
    ' Raw Formula: -5/67 / x2 * tanh(x8 * x5 + 2117/961) + exp(x5 * 345/509) + x6 / (x6 + 435/526)
    
    Dim Term1 As Double
    Dim Term2 As Double
    Dim Term3 As Double
    
    ' Term 1: (-5/67 / x2) * tanh(...)
    ' Note: VBA performs division/mult left-to-right.
    Dim inner_tanh As Double
    inner_tanh = (x8 * x5) + (2117# / 961#)
    Term1 = ((-5# / 67#) / x2) * FnTanh(inner_tanh)
    
    ' Term 2: exp(x5 * 345/509)
    Term2 = Exp(x5 * (345# / 509#))
    
    ' Term 3: x6 / (x6 + 435/526)
    Term3 = x6 / (x6 + (435# / 526#))
    
    Dim Multiplier As Double
    Multiplier = Term1 + Term2 + Term3
    
    ' 5. FINAL CALCULATION
    ' ---------------------------------------------------------
    L_genetic = L_linear * Multiplier

End Function

' Helper Function: Tanh (Hyperbolic Tangent) for VBA
' Required because older Excel VBA versions do not have Tanh built-in.
Private Function FnTanh(x As Double) As Double
    ' Safe Tanh implementation for VBA
    If x > 20 Then
        FnTanh = 1
    ElseIf x < -20 Then
        FnTanh = -1
    Else
        Dim ex2 As Double
        ex2 = Exp(2 * x)
        FnTanh = (ex2 - 1) / (ex2 + 1)
    End If
End Function

