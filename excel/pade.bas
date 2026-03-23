Attribute VB_Name = "pade"
' -----------------------------------------------------------------------
' WAVELENGTH CALCULATOR (Padé Regime-Based Surrogate)
' =====================================================================
' PURPOSE:
'   Computes wavelength using a compact regime-dependent Padé / rational
'   approximation calibrated against the reference nonlinear dataset.
'
' PUBLIC ENTRY POINTS:
'   L_pade(H, T, d, Uc)
'   L(H, T, d, Uc)     ' Alias to L_pade
'
' INPUTS:
'   H  : Wave Height [m]
'   T  : Wave Period [s]
'   d  : Water Depth [m]
'   Uc : Current Velocity [m/s] (+ following, - opposing)
'
' RETURNS:
'   L_pade : Predicted wavelength [m]
'
' MODEL STRUCTURE:
'   1. Solve a linear Doppler baseline wavelength.
'   2. Classify relative-depth regime.
'   3. Apply the corresponding Padé correction formula.
'
' ROLE IN WORKBOOK:
'   Fast explicit surrogate with low runtime cost and clear engineering
'   behaviour across shallow, intermediate, and deep-water regimes.
' -----------------------------------------------------------------------

Private Const G As Double = 9.80665
Private Const PI As Double = 3.14159265358979
Private Const REL_SHALLOW As Double = 0.05
Private Const REL_DEEP As Double = 0.5
Private Const EPS_K As Double = 0.0001
Private Const EPS_LOG As Double = 0.0000001
Private Const EPS_DEN As Double = 0.001
Private Const LOG_CORR_MIN As Double = -2.5
Private Const LOG_CORR_MAX As Double = 2.5

Private Const REG_SHALLOW As Long = 0
Private Const REG_INTER As Long = 1
Private Const REG_DEEP As Long = 2

Public Function L_pade(ByVal H As Double, ByVal T As Double, _
                       ByVal d As Double, ByVal Uc As Double) As Double
    Dim L_lin As Double
    Dim L_base As Double
    Dim rel_depth As Double
    Dim reg As Long

    If T <= 0# Or d <= 0# Then
        L_pade = 0#
        Exit Function
    End If

    L_lin = SolveLinearDoppler(T, d, Uc)
    L_base = MaxVal(L_lin, 0.1)
    rel_depth = d / MaxVal(L_lin, EPS_K)

    If rel_depth < REL_SHALLOW Then
        reg = REG_SHALLOW
    ElseIf rel_depth < REL_DEEP Then
        reg = REG_INTER
    Else
        reg = REG_DEEP
    End If

    L_pade = PredictRegime(H, T, d, Uc, L_base, reg)
End Function

Public Function L(ByVal H As Double, ByVal T As Double, _
                  ByVal d As Double, ByVal Uc As Double) As Double
    L = L_pade(H, T, d, Uc)
End Function

Private Function SolveLinearDoppler(ByVal T As Double, ByVal d As Double, _
                                    ByVal Uc As Double) As Double
    Dim L0 As Double
    Dim omega As Double
    Dim c_shallow As Double
    Dim k As Double
    Dim kd As Double
    Dim th As Double
    Dim sigma As Double
    Dim f As Double
    Dim sech2 As Double
    Dim d_sigma As Double
    Dim df As Double
    Dim k_new As Double
    Dim diff_val As Double
    Dim it As Long

    omega = 2# * PI / T
    c_shallow = Sqr(G * d)

    If Uc < -0.5 * c_shallow Then
        k = omega / (c_shallow * 0.6)
    Else
        L0 = (G * T * T) / (2# * PI)
        k = 2# * PI / L0
    End If

    For it = 1 To 40
        If k < EPS_K Then k = EPS_K
        If k > 200# Then k = 200#

        kd = k * d
        th = FnTanh(kd)
        sigma = Sqr(G * k * th)
        f = sigma + k * Uc - omega
        sech2 = 1# - th * th

        If sigma > 0.000000001 Then
            d_sigma = (G * th + G * kd * sech2) / (2# * sigma)
        Else
            d_sigma = 0#
        End If

        df = d_sigma + Uc
        If Abs(df) < 0.000000001 Then Exit For

        k_new = k - f / df
        diff_val = Abs(k_new - k)
        k = 0.8 * k + 0.2 * k_new

        If diff_val < 0.0000001 Then Exit For
    Next it

    If k <= EPS_K Then k = EPS_K
    SolveLinearDoppler = 2# * PI / k
End Function

Private Sub Features(ByVal H As Double, ByVal T As Double, ByVal d As Double, _
                     ByVal Uc As Double, ByVal L_base As Double, _
                     ByRef x0_raw As Double, ByRef x1_raw As Double, _
                     ByRef x2_raw As Double, ByRef x3_raw As Double)
    Dim L_feat As Double
    Dim d_safe As Double

    L_feat = MaxVal(L_base, 0.1)
    x0_raw = Log(MaxVal(H / L_feat, EPS_LOG))
    x1_raw = Log(MaxVal(d / L_feat, EPS_LOG))
    x2_raw = (Uc * T) / L_feat
    d_safe = MaxVal(d, 0.1)
    x3_raw = Log(MaxVal((H * L_feat * L_feat) / _
                        (d_safe * d_safe * d_safe), EPS_LOG))
End Sub

Private Sub ComputePowers(ByVal max_degree As Long, ByVal x0 As Double, _
                          ByVal x1 As Double, ByVal x2 As Double, _
                          ByVal x3 As Double, ByRef p0() As Double, _
                          ByRef p1() As Double, ByRef p2() As Double, _
                          ByRef p3() As Double)
    Dim ip As Long

    p0(0) = 1#: p1(0) = 1#: p2(0) = 1#: p3(0) = 1#
    For ip = 1 To max_degree
        p0(ip) = p0(ip - 1) * x0
        p1(ip) = p1(ip - 1) * x1
        p2(ip) = p2(ip - 1) * x2
        p3(ip) = p3(ip - 1) * x3
    Next ip
End Sub

Private Function PredictRegime(ByVal H As Double, ByVal T As Double, _
                               ByVal d As Double, ByVal Uc As Double, _
                               ByVal L_base As Double, ByVal reg As Long) As Double
    Dim x0_raw As Double, x1_raw As Double
    Dim x2_raw As Double, x3_raw As Double
    Dim x0 As Double, x1 As Double, x2 As Double, x3 As Double
    Dim m0 As Double, m1 As Double, m2 As Double, m3 As Double
    Dim s0 As Double, s1 As Double, s2 As Double, s3 As Double
    Dim p0(0 To 5) As Double
    Dim p1(0 To 5) As Double
    Dim p2(0 To 5) As Double
    Dim p3(0 To 5) As Double
    Dim poly As Double
    Dim P As Double, Q As Double
    Dim log_corr As Double

    Call Features(H, T, d, Uc, L_base, x0_raw, x1_raw, x2_raw, x3_raw)

    Select Case reg
    Case REG_SHALLOW
        m0 = -3.84276294002735
        m1 = -3.22269430421023
        m2 = 0.172137743622236
        m3 = 5.82531997260334
        s0 = 0.391277290256598
        s1 = 0.188117617935001
        s2 = 0.142509166671896
        s3 = 0.623346387573241

        x0 = (x0_raw - m0) / s0
        x1 = (x1_raw - m1) / s1
        x2 = (x2_raw - m2) / s2
        x3 = (x3_raw - m3) / s3

        Call ComputePowers(3, x0, x1, x2, x3, p0, p1, p2, p3)

        poly = 0#
        poly = poly + 0.122658891590006
        poly = poly + -0.223162432469857 * p0(1)
        poly = poly + 0.122154375864433 * p1(1)
        poly = poly + -0.032557341198453 * p2(1)
        poly = poly + 2.57536434260184E-02 * p3(1)
        poly = poly + 0.257331371712484 * p0(2)
        poly = poly + -0.308771404759481 * p0(1) * p1(1)
        poly = poly + 0.138200968009481 * p0(1) * p2(1)
        poly = poly + -0.220137192146353 * p0(1) * p3(1)
        poly = poly + 6.82539989430941E-02 * p1(2)
        poly = poly + -0.139487186293543 * p1(1) * p2(1)
        poly = poly + -1.73676555918994E-03 * p1(1) * p3(1)
        poly = poly + -8.65028847534429E-03 * p2(2)
        poly = poly + -0.130879558323783 * p2(1) * p3(1)
        poly = poly + -6.62448649411293E-02 * p3(2)
        poly = poly + 0.186956878005965 * p0(3)
        poly = poly + -0.12925232447195 * p0(2) * p1(1)
        poly = poly + 0.216879015407666 * p0(2) * p2(1)
        poly = poly + -0.312039395521106 * p0(2) * p3(1)
        poly = poly + -0.154345478551834 * p0(1) * p1(2)
        poly = poly + -0.373278850829006 * p0(1) * p1(1) * p2(1)
        poly = poly + -0.124310452436107 * p0(1) * p1(1) * p3(1)
        poly = poly + 7.20011253746945E-02 * p0(1) * p2(2)
        poly = poly + -0.287588778092046 * p0(1) * p2(1) * p3(1)
        poly = poly + -7.60134194128335E-03 * p0(1) * p3(2)
        poly = poly + -6.47783993781433E-02 * p1(3)
        poly = poly + 7.98497398531808E-02 * p1(2) * p2(1)
        poly = poly + -0.11129694402868 * p1(2) * p3(1)
        poly = poly + -0.104170470140354 * p1(1) * p2(2)
        poly = poly + 9.87890014999881E-03 * p1(1) * p2(1) * p3(1)
        poly = poly + 0.10984275632681 * p1(1) * p3(2)
        poly = poly + -1.76417179219626E-03 * p2(3)
        poly = poly + -0.107352889557145 * p2(2) * p3(1)
        poly = poly + -8.35623109010058E-02 * p2(1) * p3(2)
        poly = poly + 0.167229374489742 * p3(3)
        P = poly

        poly = 1#
        poly = poly + -1.45145515254481 * p0(1)
        poly = poly + 0.126886479878433 * p1(1)
        poly = poly + -8.79634413123627E-02 * p2(1)
        poly = poly + -1.02596276940039 * p3(1)
        poly = poly + 0.962884351101783 * p0(2)
        poly = poly + -7.31762766321935E-02 * p0(1) * p1(1)
        poly = poly + 0.139802357215793 * p0(1) * p2(1)
        poly = poly + 0.670657645121159 * p0(1) * p3(1)
        poly = poly + -2.98881255268956E-02 * p1(2)
        poly = poly + 1.36052608122462E-02 * p1(1) * p2(1)
        poly = poly + -1.88735614204462E-02 * p1(1) * p3(1)
        poly = poly + -6.90961025695448E-02 * p2(2)
        poly = poly + 7.54369012639411E-02 * p2(1) * p3(1)
        poly = poly + 0.438062142925391 * p3(2)
        poly = poly + -0.135196555657141 * p0(3)
        poly = poly + 1.53048375533114E-02 * p0(2) * p1(1)
        poly = poly + 4.55031931454175E-02 * p0(2) * p2(1)
        poly = poly + -9.87198639308659E-02 * p0(2) * p3(1)
        poly = poly + -8.70034410285395E-04 * p0(1) * p1(2)
        poly = poly + -1.50650346005766E-02 * p0(1) * p1(1) * p2(1)
        poly = poly + 1.03946085289015E-02 * p0(1) * p1(1) * p3(1)
        poly = poly + 4.03177660803925E-02 * p0(1) * p2(2)
        poly = poly + 4.22018351380404E-02 * p0(1) * p2(1) * p3(1)
        poly = poly + -7.13777583062031E-02 * p0(1) * p3(2)
        poly = poly + 2.79733351513631E-03 * p1(3)
        poly = poly + -8.98159249807842E-03 * p1(2) * p2(1)
        poly = poly + -3.07871762783684E-03 * p1(2) * p3(1)
        poly = poly + -1.88057339551477E-02 * p1(1) * p2(2)
        poly = poly + -1.32481527820806E-03 * p1(1) * p2(1) * p3(1)
        poly = poly + 9.31208956672424E-03 * p1(1) * p3(2)
        poly = poly + -1.08533474918853E-02 * p2(3)
        poly = poly + 4.23335989032786E-02 * p2(2) * p3(1)
        poly = poly + 2.76897135965507E-02 * p2(1) * p3(2)
        poly = poly + -0.053234928586283 * p3(3)
        Q = poly

    Case REG_INTER
        m0 = -3.53130688174989
        m1 = -1.67429460487877
        m2 = 1.50335096000887E-03
        m3 = 1.49157693288641
        s0 = 0.711777902074339
        s1 = 0.48853216529815
        s2 = 0.134203739319999
        s3 = 1.50697715170332

        x0 = (x0_raw - m0) / s0
        x1 = (x1_raw - m1) / s1
        x2 = (x2_raw - m2) / s2
        x3 = (x3_raw - m3) / s3

        Call ComputePowers(3, x0, x1, x2, x3, p0, p1, p2, p3)

        poly = 0#
        poly = poly + 1.07601956662038E-02
        poly = poly + -9.12423459356319E-02 * p0(1)
        poly = poly + 0.183722838767907 * p1(1)
        poly = poly + -3.43747435979402E-04 * p2(1)
        poly = poly + 0.189606763941413 * p3(1)
        poly = poly + -0.232308855215674 * p0(2)
        poly = poly + 0.133775070931928 * p0(1) * p1(1)
        poly = poly + 0.069301346894819 * p0(1) * p2(1)
        poly = poly + 0.326359969729259 * p0(1) * p3(1)
        poly = poly + 0.695263570031413 * p1(2)
        poly = poly + -0.148967216618729 * p1(1) * p2(1)
        poly = poly + 1.0279626010453 * p1(1) * p3(1)
        poly = poly + 5.50583989197696E-04 * p2(2)
        poly = poly + -0.151670762076735 * p2(1) * p3(1)
        poly = poly + 0.327297909931705 * p3(2)
        poly = poly + -0.138610707375903 * p0(3)
        poly = poly + 0.13154177010249 * p0(2) * p1(1)
        poly = poly + -0.141827502253434 * p0(2) * p2(1)
        poly = poly + 0.360721953982323 * p0(2) * p3(1)
        poly = poly + 0.231290304289201 * p0(1) * p1(2)
        poly = poly + 0.372892379904205 * p0(1) * p1(1) * p2(1)
        poly = poly + 0.250000264163254 * p0(1) * p1(1) * p3(1)
        poly = poly + 0.190120974338761 * p0(1) * p2(2)
        poly = poly + 0.193519010104878 * p0(1) * p2(1) * p3(1)
        poly = poly + -1.95276475309372E-02 * p0(1) * p3(2)
        poly = poly + 0.185080974822618 * p1(3)
        poly = poly + -0.159912570724368 * p1(2) * p2(1)
        poly = poly + 0.090822668623141 * p1(2) * p3(1)
        poly = poly + -0.391486546148403 * p1(1) * p2(2)
        poly = poly + 6.15543088340233E-02 * p1(1) * p2(1) * p3(1)
        poly = poly + -0.337838393859852 * p1(1) * p3(2)
        poly = poly + -1.30149043702889E-04 * p2(3)
        poly = poly + -0.402625505319557 * p2(2) * p3(1)
        poly = poly + 0.231764538322574 * p2(1) * p3(2)
        poly = poly + -0.242437725567292 * p3(3)
        P = poly

        poly = 1#
        poly = poly + -1.28414928855244 * p0(1)
        poly = poly + 7.84924386226821E-02 * p1(1)
        poly = poly + 5.91047329224304E-02 * p2(1)
        poly = poly + -0.41981950654197 * p3(1)
        poly = poly + 0.766408761399911 * p0(2)
        poly = poly + -0.662675151056155 * p0(1) * p1(1)
        poly = poly + -0.13406101427693 * p0(1) * p2(1)
        poly = poly + 0.285949522416194 * p0(1) * p3(1)
        poly = poly + 0.385395629988835 * p1(2)
        poly = poly + 8.51644632323136E-02 * p1(1) * p2(1)
        poly = poly + -0.265573227686085 * p1(1) * p3(1)
        poly = poly + 0.01726608092209 * p2(2)
        poly = poly + 8.85871425726676E-02 * p2(1) * p3(1)
        poly = poly + -0.565194869680856 * p3(2)
        poly = poly + -0.207092632594496 * p0(3)
        poly = poly + 0.695047685172946 * p0(2) * p1(1)
        poly = poly + 7.02492430451155E-02 * p0(2) * p2(1)
        poly = poly + -8.67232798226327E-02 * p0(2) * p3(1)
        poly = poly + -1.43336324220046 * p0(1) * p1(2)
        poly = poly + 0.124478545225257 * p0(1) * p1(1) * p2(1)
        poly = poly + 5.89129000085894E-03 * p0(1) * p1(1) * p3(1)
        poly = poly + -0.125829410744164 * p0(1) * p2(2)
        poly = poly + 3.48423448156924E-03 * p0(1) * p2(1) * p3(1)
        poly = poly + 8.49475502503417E-02 * p0(1) * p3(2)
        poly = poly + 1.12971147774955 * p1(3)
        poly = poly + -0.442404371971004 * p1(2) * p2(1)
        poly = poly + 0.362391153998164 * p1(2) * p3(1)
        poly = poly + 0.250900869616184 * p1(1) * p2(2)
        poly = poly + -0.639661961537757 * p1(1) * p2(1) * p3(1)
        poly = poly + 1.54567546687247E-02 * p1(1) * p3(2)
        poly = poly + -2.52384985675274E-04 * p2(3)
        poly = poly + 0.25298188495864 * p2(2) * p3(1)
        poly = poly + -0.192457999314718 * p2(1) * p3(2)
        poly = poly + 0.851433854411603 * p3(3)
        Q = poly

    Case REG_DEEP
        m0 = -2.74456246031991
        m1 = -0.174860770385788
        m2 = 8.5264201624294E-04
        m3 = -2.21998014916254
        s0 = 0.570366357511055
        s1 = 0.440791669311128
        s2 = 0.180567751120582
        s3 = 1.27954215870031

        x0 = (x0_raw - m0) / s0
        x1 = (x1_raw - m1) / s1
        x2 = (x2_raw - m2) / s2
        x3 = (x3_raw - m3) / s3

        Call ComputePowers(5, x0, x1, x2, x3, p0, p1, p2, p3)

        poly = 0#
        poly = poly + 0.034426196370823
        poly = poly + -1.30728970965385 * p0(1)
        poly = poly + 2.811483877117 * p1(1)
        poly = poly + -7.07647377280387E-03 * p2(1)
        poly = poly + 2.72606501754838 * p3(1)
        poly = poly + 1.04324143340433 * p0(2)
        poly = poly + -2.01883222802041 * p0(1) * p1(1)
        poly = poly + 0.587645474165238 * p0(1) * p2(1)
        poly = poly + -2.33852137910641 * p0(1) * p3(1)
        poly = poly + -0.509765982740888 * p1(2)
        poly = poly + -1.29573198773404 * p1(1) * p2(1)
        poly = poly + -0.127889296716919 * p1(1) * p3(1)
        poly = poly + -1.45090030224523E-03 * p2(2)
        poly = poly + -1.2559431956521 * p2(1) * p3(1)
        poly = poly + 0.358394020366298 * p3(2)
        poly = poly + 0.217064780533775 * p0(3)
        poly = poly + -0.656137330252744 * p0(2) * p1(1)
        poly = poly + 1.26620131293536 * p0(2) * p2(1)
        poly = poly + -0.431757821517899 * p0(2) * p3(1)
        poly = poly + 0.420560524047101 * p0(1) * p1(2)
        poly = poly + -3.0121545955335 * p0(1) * p1(1) * p2(1)
        poly = poly + 0.185818527558822 * p0(1) * p1(1) * p3(1)
        poly = poly + 0.151638742655726 * p0(1) * p2(2)
        poly = poly + -3.13256065065272 * p0(1) * p2(1) * p3(1)
        poly = poly + -5.69755074883839E-02 * p0(1) * p3(2)
        poly = poly + -0.214325086407473 * p1(3)
        poly = poly + 8.94399592031439E-04 * p1(2) * p2(1)
        poly = poly + -0.199587007293985 * p1(2) * p3(1)
        poly = poly + -0.355639408096815 * p1(1) * p2(2)
        poly = poly + 0.511226847014798 * p1(1) * p2(1) * p3(1)
        poly = poly + -0.177011132857736 * p1(1) * p3(2)
        poly = poly + 6.71202810507883E-04 * p2(3)
        poly = poly + -0.345320358136189 * p2(2) * p3(1)
        poly = poly + 0.493149163810967 * p2(1) * p3(2)
        poly = poly + -0.178616342551935 * p3(3)
        poly = poly + -0.452083806632251 * p0(4)
        poly = poly + 0.929590433854566 * p0(3) * p1(1)
        poly = poly + -0.227463161677271 * p0(3) * p2(1)
        poly = poly + 1.37821278706062 * p0(3) * p3(1)
        poly = poly + 0.20110524805254 * p0(2) * p1(2)
        poly = poly + 0.751117817359393 * p0(2) * p1(1) * p2(1)
        poly = poly + -0.865280078507087 * p0(2) * p1(1) * p3(1)
        poly = poly + 0.548835233760403 * p0(2) * p2(2)
        poly = poly + 0.58658928165099 * p0(2) * p2(1) * p3(1)
        poly = poly + -0.877588541138445 * p0(2) * p3(2)
        poly = poly + -4.89004238119476E-02 * p0(1) * p1(3)
        poly = poly + -0.526233676992449 * p0(1) * p1(2) * p2(1)
        poly = poly + 0.187717078324943 * p0(1) * p1(2) * p3(1)
        poly = poly + -1.21952592110831 * p0(1) * p1(1) * p2(2)
        poly = poly + -0.592863832796373 * p0(1) * p1(1) * p2(1) * p3(1)
        poly = poly + 0.186759925386856 * p0(1) * p1(1) * p3(2)
        poly = poly + -2.79606389893333E-02 * p0(1) * p2(3)
        poly = poly + -1.22413191383729 * p0(1) * p2(2) * p3(1)
        poly = poly + -0.219524884189707 * p0(1) * p2(1) * p3(2)
        poly = poly + -1.78262510086378E-03 * p0(1) * p3(3)
        poly = poly + -1.11381040379119E-02 * p1(4)
        poly = poly + -4.48194420526048E-03 * p1(3) * p2(1)
        poly = poly + -2.88024431873925E-02 * p1(3) * p3(1)
        poly = poly + -6.68561867553627E-02 * p1(2) * p2(2)
        poly = poly + -0.206102204695054 * p1(2) * p2(1) * p3(1)
        poly = poly + -3.20082095194857E-02 * p1(2) * p3(2)
        poly = poly + 6.33739806546677E-02 * p1(1) * p2(3)
        poly = poly + -2.16843871585281E-02 * p1(1) * p2(2) * p3(1)
        poly = poly + -6.55504444159509E-02 * p1(1) * p2(1) * p3(2)
        poly = poly + -0.108839292687417 * p1(1) * p3(3)
        poly = poly + 3.5738115909278E-04 * p2(4)
        poly = poly + 6.18374982238066E-02 * p2(3) * p3(1)
        poly = poly + 0.041693635667715 * p2(2) * p3(2)
        poly = poly + 0.125596499838677 * p2(1) * p3(3)
        poly = poly + -9.17298728161615E-02 * p3(4)
        poly = poly + -6.68914077560428E-02 * p0(5)
        poly = poly + 7.05228313949367E-02 * p0(4) * p1(1)
        poly = poly + 0.240708004857002 * p0(4) * p2(1)
        poly = poly + 0.283532148921935 * p0(4) * p3(1)
        poly = poly + 0.244973577421466 * p0(3) * p1(2)
        poly = poly + -0.756182868996171 * p0(3) * p1(1) * p2(1)
        poly = poly + -5.18317098147288E-02 * p0(3) * p1(1) * p3(1)
        poly = poly + 0.614161630192344 * p0(3) * p2(2)
        poly = poly + -1.04399492339536 * p0(3) * p2(1) * p3(1)
        poly = poly + -0.220476664856711 * p0(3) * p3(2)
        poly = poly + -0.065423390992346 * p0(2) * p1(3)
        poly = poly + 0.478398700971647 * p0(2) * p1(2) * p2(1)
        poly = poly + -5.77362022954065E-02 * p0(2) * p1(2) * p3(1)
        poly = poly + -1.53347030888513 * p0(2) * p1(1) * p2(2)
        poly = poly + 1.97910814934406 * p0(2) * p1(1) * p2(1) * p3(1)
        poly = poly + -0.170041232476451 * p0(2) * p1(1) * p3(2)
        poly = poly + 0.091825984261077 * p0(2) * p2(3)
        poly = poly + -1.22336381862551 * p0(2) * p2(2) * p3(1)
        poly = poly + 0.980173204097067 * p0(2) * p2(1) * p3(2)
        poly = poly + -0.100757026702261 * p0(2) * p3(3)
        poly = poly + 4.61805560194696E-02 * p0(1) * p1(4)
        poly = poly + 0.154628561826729 * p0(1) * p1(3) * p2(1)
        poly = poly + 0.025980224621631 * p0(1) * p1(3) * p3(1)
        poly = poly + 0.197116417364951 * p0(1) * p1(2) * p2(2)
        poly = poly + -0.401806059350958 * p0(1) * p1(2) * p2(1) * p3(1)
        poly = poly + -3.20667756912502E-02 * p0(1) * p1(2) * p3(2)
        poly = poly + -0.208477222572714 * p0(1) * p1(1) * p2(3)
        poly = poly + -0.203951140797432 * p0(1) * p1(1) * p2(2) * p3(1)
        poly = poly + -0.034920586889727 * p0(1) * p1(1) * p2(1) * p3(2)
        poly = poly + -3.43241401273053E-02 * p0(1) * p1(1) * p3(3)
        poly = poly + -1.02698105873091E-02 * p0(1) * p2(4)
        poly = poly + -0.205956987347702 * p0(1) * p2(3) * p3(1)
        poly = poly + -0.385869899828892 * p0(1) * p2(2) * p3(2)
        poly = poly + 0.492053632828282 * p0(1) * p2(1) * p3(3)
        poly = poly + -6.94249231245914E-02 * p0(1) * p3(4)
        poly = poly + 1.29787107522251E-02 * p1(5)
        poly = poly + -4.16664398088651E-02 * p1(4) * p2(1)
        poly = poly + -1.53023079636677E-02 * p1(4) * p3(1)
        poly = poly + 3.56814664868936E-02 * p1(3) * p2(2)
        poly = poly + -5.62049556365616E-02 * p1(3) * p2(1) * p3(1)
        poly = poly + -4.12202061905592E-03 * p1(3) * p3(2)
        poly = poly + -1.15195869238608E-02 * p1(2) * p2(3)
        poly = poly + 5.38315719093637E-02 * p1(2) * p2(2) * p3(1)
        poly = poly + -3.47463864198555E-02 * p1(2) * p2(1) * p3(2)
        poly = poly + -2.38078878604471E-02 * p1(2) * p3(3)
        poly = poly + 2.28967189000334E-02 * p1(1) * p2(4)
        poly = poly + -1.38977232702885E-02 * p1(1) * p2(3) * p3(1)
        poly = poly + 1.29254019245706E-02 * p1(1) * p2(2) * p3(2)
        poly = poly + -8.27073179940645E-02 * p1(1) * p2(1) * p3(3)
        poly = poly + -3.24023387057993E-02 * p1(1) * p3(4)
        poly = poly + 9.02380413437505E-06 * p2(5)
        poly = poly + 2.22296261631702E-02 * p2(4) * p3(1)
        poly = poly + -2.55835935828475E-03 * p2(3) * p3(2)
        poly = poly + -5.54315251358381E-03 * p2(2) * p3(3)
        poly = poly + -6.18804608283115E-02 * p2(1) * p3(4)
        poly = poly + 1.16251392997031E-02 * p3(5)
        P = poly

        poly = 1#
        poly = poly + -3.20798197234749 * p0(1)
        poly = poly + -0.724320562130214 * p1(1)
        poly = poly + 4.77249031861665E-02 * p2(1)
        poly = poly + -0.681416844855583 * p3(1)
        poly = poly + 4.37614336217539 * p0(2)
        poly = poly + 0.976400831526362 * p0(1) * p1(1)
        poly = poly + -0.14532356502234 * p0(1) * p2(1)
        poly = poly + 0.941615629305515 * p0(1) * p3(1)
        poly = poly + 0.182852057386867 * p1(2)
        poly = poly + -1.03715038438628E-02 * p1(1) * p2(1)
        poly = poly + 0.24626558399216 * p1(1) * p3(1)
        poly = poly + -1.67270832230894E-02 * p2(2)
        poly = poly + -0.05406046899581 * p2(1) * p3(1)
        poly = poly + 0.165223485905615 * p3(2)
        poly = poly + -2.97774030252351 * p0(3)
        poly = poly + -0.636717022659393 * p0(2) * p1(1)
        poly = poly + 0.252488918702817 * p0(2) * p2(1)
        poly = poly + -0.669320818624025 * p0(2) * p3(1)
        poly = poly + -0.156992867909834 * p0(1) * p1(2)
        poly = poly + -1.62839621614379E-04 * p0(1) * p1(1) * p2(1)
        poly = poly + -0.121573582743352 * p0(1) * p1(1) * p3(1)
        poly = poly + -6.28207316470531E-02 * p0(1) * p2(2)
        poly = poly + 0.112717294102726 * p0(1) * p2(1) * p3(1)
        poly = poly + -0.172711946430928 * p0(1) * p3(2)
        poly = poly + -4.30961662873439E-02 * p1(3)
        poly = poly + -1.82950163918499E-02 * p1(2) * p2(1)
        poly = poly + -2.54420388626298E-02 * p1(2) * p3(1)
        poly = poly + -1.13768707185896E-02 * p1(1) * p2(2)
        poly = poly + 1.88348560427273E-02 * p1(1) * p2(1) * p3(1)
        poly = poly + -2.78987021507144E-02 * p1(1) * p3(2)
        poly = poly + 4.25441534306874E-02 * p2(3)
        poly = poly + -1.62451457362522E-02 * p2(2) * p3(1)
        poly = poly + 3.07792929188282E-02 * p2(1) * p3(2)
        poly = poly + -4.81551392843235E-02 * p3(3)
        poly = poly + 1.02901341629659 * p0(4)
        poly = poly + 0.193964119110369 * p0(3) * p1(1)
        poly = poly + -0.220885354911628 * p0(3) * p2(1)
        poly = poly + 0.258234028192201 * p0(3) * p3(1)
        poly = poly + 7.16725871555428E-02 * p0(2) * p1(2)
        poly = poly + 8.9965521406317E-03 * p0(2) * p1(1) * p2(1)
        poly = poly + 1.23892554030648E-02 * p0(2) * p1(1) * p3(1)
        poly = poly + 0.178805703479562 * p0(2) * p2(2)
        poly = poly + -0.107759157391984 * p0(2) * p2(1) * p3(1)
        poly = poly + 0.102305940324051 * p0(2) * p3(2)
        poly = poly + 3.08843447820079E-02 * p0(1) * p1(3)
        poly = poly + 1.98689955511547E-02 * p0(1) * p1(2) * p2(1)
        poly = poly + 3.04363223820002E-05 * p0(1) * p1(2) * p3(1)
        poly = poly + 3.94530518348764E-02 * p0(1) * p1(1) * p2(2)
        poly = poly + -1.65238280800464E-02 * p0(1) * p1(1) * p2(1) * p3(1)
        poly = poly + 5.49115386528278E-03 * p0(1) * p1(1) * p3(2)
        poly = poly + -0.101863397256816 * p0(1) * p2(3)
        poly = poly + 3.89303494524459E-02 * p0(1) * p2(2) * p3(1)
        poly = poly + -3.09575589183528E-02 * p0(1) * p2(1) * p3(2)
        poly = poly + 3.99287382804525E-02 * p0(1) * p3(3)
        poly = poly + 8.36216870529288E-03 * p1(4)
        poly = poly + 6.13794379641066E-03 * p1(3) * p2(1)
        poly = poly + 5.12485435001769E-03 * p1(3) * p3(1)
        poly = poly + 6.78140695727817E-03 * p1(2) * p2(2)
        poly = poly + 2.5133537540195E-03 * p1(2) * p2(1) * p3(1)
        poly = poly + -5.28283961849323E-03 * p1(2) * p3(2)
        poly = poly + -2.54596795930958E-02 * p1(1) * p2(3)
        poly = poly + 0.01057810070633 * p1(1) * p2(2) * p3(1)
        poly = poly + -9.96311361407819E-03 * p1(1) * p2(1) * p3(2)
        poly = poly + 7.90741559390111E-03 * p1(1) * p3(3)
        poly = poly + 0.011019866632944 * p2(4)
        poly = poly + -1.90944966314051E-02 * p2(3) * p3(1)
        poly = poly + 6.42131407996161E-03 * p2(2) * p3(2)
        poly = poly + -3.50295496003817E-03 * p2(1) * p3(3)
        poly = poly + 9.62644723524048E-03 * p3(4)
        poly = poly + -0.139612935865398 * p0(5)
        poly = poly + -1.61227519738391E-02 * p0(4) * p1(1)
        poly = poly + 6.57582738550188E-02 * p0(4) * p2(1)
        poly = poly + -4.55711499310408E-02 * p0(4) * p3(1)
        poly = poly + -1.16596229987947E-02 * p0(3) * p1(2)
        poly = poly + -6.57658204949983E-03 * p0(3) * p1(1) * p2(1)
        poly = poly + 4.86308156601893E-03 * p0(3) * p1(1) * p3(1)
        poly = poly + -0.10036373145077 * p0(3) * p2(2)
        poly = poly + 3.61090242124856E-02 * p0(3) * p2(1) * p3(1)
        poly = poly + -2.53395893066711E-02 * p0(3) * p3(2)
        poly = poly + -7.35753658704463E-03 * p0(2) * p1(3)
        poly = poly + -7.98573896074787E-03 * p0(2) * p1(2) * p2(1)
        poly = poly + 2.40645651190969E-03 * p0(2) * p1(2) * p3(1)
        poly = poly + -2.38632874835867E-02 * p0(2) * p1(1) * p2(2)
        poly = poly + 5.32150051019692E-03 * p0(2) * p1(1) * p2(1) * p3(1)
        poly = poly + -3.19251525955071E-04 * p0(2) * p1(1) * p3(2)
        poly = poly + 6.35477992389113E-02 * p0(2) * p2(3)
        poly = poly + -2.00758420025715E-02 * p0(2) * p2(2) * p3(1)
        poly = poly + 1.05962458910113E-02 * p0(2) * p2(1) * p3(2)
        poly = poly + -1.09653884466126E-02 * p0(2) * p3(3)
        poly = poly + -2.4897251686087E-03 * p0(1) * p1(4)
        poly = poly + -2.28404714959298E-03 * p0(1) * p1(3) * p2(1)
        poly = poly + -7.0661287854459E-04 * p0(1) * p1(3) * p3(1)
        poly = poly + -5.07384675059637E-03 * p0(1) * p1(2) * p2(2)
        poly = poly + -1.1991970570974E-03 * p0(1) * p1(2) * p2(1) * p3(1)
        poly = poly + 1.80296720890438E-03 * p0(1) * p1(2) * p3(2)
        poly = poly + 1.67509237597409E-02 * p0(1) * p1(1) * p2(3)
        poly = poly + -5.39355635192798E-03 * p0(1) * p1(1) * p2(2) * p3(1)
        poly = poly + 3.61145320924182E-03 * p0(1) * p1(1) * p2(1) * p3(2)
        poly = poly + -2.0056280554991E-03 * p0(1) * p1(1) * p3(3)
        poly = poly + -1.15009649201076E-02 * p0(1) * p2(4)
        poly = poly + 1.10152869027075E-02 * p0(1) * p2(3) * p3(1)
        poly = poly + -3.3748620549211E-03 * p0(1) * p2(2) * p3(2)
        poly = poly + 9.91023657467516E-04 * p0(1) * p2(1) * p3(3)
        poly = poly + -2.81515051342284E-03 * p0(1) * p3(4)
        poly = poly + -2.26489248901108E-04 * p1(5)
        poly = poly + 7.0944459020728E-05 * p1(4) * p2(1)
        poly = poly + -8.75745646046383E-04 * p1(4) * p3(1)
        poly = poly + -9.70882179776355E-04 * p1(3) * p2(2)
        poly = poly + -1.0914563456733E-03 * p1(3) * p2(1) * p3(1)
        poly = poly + 5.90084871048974E-04 * p1(3) * p3(2)
        poly = poly + 3.75068689549917E-03 * p1(2) * p2(3)
        poly = poly + -1.25832908070351E-03 * p1(2) * p2(2) * p3(1)
        poly = poly + 5.9343586031131E-04 * p1(2) * p2(1) * p3(2)
        poly = poly + 1.93849122566873E-04 * p1(2) * p3(3)
        poly = poly + -2.98043857440057E-03 * p1(1) * p2(4)
        poly = poly + 3.5906195176163E-03 * p1(1) * p2(3) * p3(1)
        poly = poly + -1.10377808289296E-03 * p1(1) * p2(2) * p3(2)
        poly = poly + 9.965284864964E-04 * p1(1) * p2(1) * p3(3)
        poly = poly + -1.09436506364567E-03 * p1(1) * p3(4)
        poly = poly + 7.41023220590129E-04 * p2(5)
        poly = poly + -2.04643992208591E-03 * p2(4) * p3(1)
        poly = poly + 1.19933870524047E-03 * p2(3) * p3(2)
        poly = poly + -3.63646615897484E-04 * p2(2) * p3(3)
        poly = poly + -5.88134651195864E-04 * p2(1) * p3(4)
        poly = poly + -1.2387611635104E-04 * p3(5)
        Q = poly
    End Select

    If Abs(Q) < EPS_DEN Then
        If Q < 0# Then
            Q = -EPS_DEN
        Else
            Q = EPS_DEN
        End If
    End If

    log_corr = P / Q
    If log_corr < LOG_CORR_MIN Then log_corr = LOG_CORR_MIN
    If log_corr > LOG_CORR_MAX Then log_corr = LOG_CORR_MAX

    PredictRegime = MaxVal(L_base, 0.1) * Exp(log_corr)
End Function

Private Function MaxVal(ByVal a As Double, ByVal b As Double) As Double
    If a >= b Then
        MaxVal = a
    Else
        MaxVal = b
    End If
End Function

Private Function FnTanh(ByVal x As Double) As Double
    If x > 20# Then
        FnTanh = 1#
    ElseIf x < -20# Then
        FnTanh = -1#
    Else
        FnTanh = (Exp(x) - Exp(-x)) / (Exp(x) + Exp(-x))
    End If
End Function

