Attribute VB_Name = "pade"
Option Explicit

Private Const G As Double = 9.80665#
Private Const PI As Double = 3.1415926535897931#
Private Const REL_SHALLOW As Double = 0.05#
Private Const REL_DEEP As Double = 0.5#
Private Const EPS_K As Double = 0.0001#
Private Const EPS_LOG As Double = 0.0000001#
Private Const EPS_DEN As Double = 0.001#
Private Const LOG_CORR_MIN As Double = -2.5#
Private Const LOG_CORR_MAX As Double = 2.5#

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
    L_base = MaxVal(L_lin, 0.1#)
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
    Dim C_shallow As Double
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
    C_shallow = Sqr(G * d)

    If Uc < -0.5# * C_shallow Then
        k = omega / (C_shallow * 0.6#)
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

        If sigma > 0.000000001# Then
            d_sigma = (G * th + G * kd * sech2) / (2# * sigma)
        Else
            d_sigma = 0#
        End If

        df = d_sigma + Uc
        If Abs(df) < 0.000000001# Then Exit For

        k_new = k - f / df
        diff_val = Abs(k_new - k)
        k = 0.8# * k + 0.2# * k_new

        If diff_val < 0.0000001# Then Exit For
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

    L_feat = MaxVal(L_base, 0.1#)
    x0_raw = Log(MaxVal(H / L_feat, EPS_LOG))
    x1_raw = Log(MaxVal(d / L_feat, EPS_LOG))
    x2_raw = (Uc * T) / L_feat
    d_safe = MaxVal(d, 0.1#)
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
        m0 = -3.8427629400273480E+00
        m1 = -3.2226943042102336E+00
        m2 = 1.7213774362223636E-01
        m3 = 5.8253199726033387E+00
        s0 = 3.9127729025659802E-01
        s1 = 1.8811761793500070E-01
        s2 = 1.4250916667189623E-01
        s3 = 6.2334638757324090E-01

        x0 = (x0_raw - m0) / s0
        x1 = (x1_raw - m1) / s1
        x2 = (x2_raw - m2) / s2
        x3 = (x3_raw - m3) / s3

        Call ComputePowers(3, x0, x1, x2, x3, p0, p1, p2, p3)

        poly = 0#
        poly = poly + 1.2265889159000636E-01
        poly = poly + -2.2316243246985726E-01 * p0(1)
        poly = poly + 1.2215437586443346E-01 * p1(1)
        poly = poly + -3.2557341198453042E-02 * p2(1)
        poly = poly + 2.5753643426018411E-02 * p3(1)
        poly = poly + 2.5733137171248410E-01 * p0(2)
        poly = poly + -3.0877140475948062E-01 * p0(1) * p1(1)
        poly = poly + 1.3820096800948112E-01 * p0(1) * p2(1)
        poly = poly + -2.2013719214635302E-01 * p0(1) * p3(1)
        poly = poly + 6.8253998943094141E-02 * p1(2)
        poly = poly + -1.3948718629354254E-01 * p1(1) * p2(1)
        poly = poly + -1.7367655591899365E-03 * p1(1) * p3(1)
        poly = poly + -8.6502884753442890E-03 * p2(2)
        poly = poly + -1.3087955832378312E-01 * p2(1) * p3(1)
        poly = poly + -6.6244864941129342E-02 * p3(2)
        poly = poly + 1.8695687800596483E-01 * p0(3)
        poly = poly + -1.2925232447194968E-01 * p0(2) * p1(1)
        poly = poly + 2.1687901540766594E-01 * p0(2) * p2(1)
        poly = poly + -3.1203939552110610E-01 * p0(2) * p3(1)
        poly = poly + -1.5434547855183373E-01 * p0(1) * p1(2)
        poly = poly + -3.7327885082900597E-01 * p0(1) * p1(1) * p2(1)
        poly = poly + -1.2431045243610674E-01 * p0(1) * p1(1) * p3(1)
        poly = poly + 7.2001125374694511E-02 * p0(1) * p2(2)
        poly = poly + -2.8758877809204642E-01 * p0(1) * p2(1) * p3(1)
        poly = poly + -7.6013419412833500E-03 * p0(1) * p3(2)
        poly = poly + -6.4778399378143323E-02 * p1(3)
        poly = poly + 7.9849739853180821E-02 * p1(2) * p2(1)
        poly = poly + -1.1129694402867968E-01 * p1(2) * p3(1)
        poly = poly + -1.0417047014035401E-01 * p1(1) * p2(2)
        poly = poly + 9.8789001499988115E-03 * p1(1) * p2(1) * p3(1)
        poly = poly + 1.0984275632681019E-01 * p1(1) * p3(2)
        poly = poly + -1.7641717921962629E-03 * p2(3)
        poly = poly + -1.0735288955714541E-01 * p2(2) * p3(1)
        poly = poly + -8.3562310901005799E-02 * p2(1) * p3(2)
        poly = poly + 1.6722937448974157E-01 * p3(3)
        P = poly

        poly = 1#
        poly = poly + -1.4514551525448056E+00 * p0(1)
        poly = poly + 1.2688647987843332E-01 * p1(1)
        poly = poly + -8.7963441312362736E-02 * p2(1)
        poly = poly + -1.0259627694003886E+00 * p3(1)
        poly = poly + 9.6288435110178328E-01 * p0(2)
        poly = poly + -7.3176276632193518E-02 * p0(1) * p1(1)
        poly = poly + 1.3980235721579284E-01 * p0(1) * p2(1)
        poly = poly + 6.7065764512115933E-01 * p0(1) * p3(1)
        poly = poly + -2.9888125526895609E-02 * p1(2)
        poly = poly + 1.3605260812246237E-02 * p1(1) * p2(1)
        poly = poly + -1.8873561420446178E-02 * p1(1) * p3(1)
        poly = poly + -6.9096102569544834E-02 * p2(2)
        poly = poly + 7.5436901263941092E-02 * p2(1) * p3(1)
        poly = poly + 4.3806214292539092E-01 * p3(2)
        poly = poly + -1.3519655565714139E-01 * p0(3)
        poly = poly + 1.5304837553311426E-02 * p0(2) * p1(1)
        poly = poly + 4.5503193145417530E-02 * p0(2) * p2(1)
        poly = poly + -9.8719863930865931E-02 * p0(2) * p3(1)
        poly = poly + -8.7003441028539463E-04 * p0(1) * p1(2)
        poly = poly + -1.5065034600576578E-02 * p0(1) * p1(1) * p2(1)
        poly = poly + 1.0394608528901499E-02 * p0(1) * p1(1) * p3(1)
        poly = poly + 4.0317766080392510E-02 * p0(1) * p2(2)
        poly = poly + 4.2201835138040378E-02 * p0(1) * p2(1) * p3(1)
        poly = poly + -7.1377758306203076E-02 * p0(1) * p3(2)
        poly = poly + 2.7973335151363072E-03 * p1(3)
        poly = poly + -8.9815924980784240E-03 * p1(2) * p2(1)
        poly = poly + -3.0787176278368376E-03 * p1(2) * p3(1)
        poly = poly + -1.8805733955147730E-02 * p1(1) * p2(2)
        poly = poly + -1.3248152782080562E-03 * p1(1) * p2(1) * p3(1)
        poly = poly + 9.3120895667242378E-03 * p1(1) * p3(2)
        poly = poly + -1.0853347491885299E-02 * p2(3)
        poly = poly + 4.2333598903278591E-02 * p2(2) * p3(1)
        poly = poly + 2.7689713596550716E-02 * p2(1) * p3(2)
        poly = poly + -5.3234928586282954E-02 * p3(3)
        Q = poly

    Case REG_INTER
        m0 = -3.5313068817498912E+00
        m1 = -1.6742946048787659E+00
        m2 = 1.5033509600088672E-03
        m3 = 1.4915769328864064E+00
        s0 = 7.1177790207433878E-01
        s1 = 4.8853216529815047E-01
        s2 = 1.3420373931999938E-01
        s3 = 1.5069771517033199E+00

        x0 = (x0_raw - m0) / s0
        x1 = (x1_raw - m1) / s1
        x2 = (x2_raw - m2) / s2
        x3 = (x3_raw - m3) / s3

        Call ComputePowers(3, x0, x1, x2, x3, p0, p1, p2, p3)

        poly = 0#
        poly = poly + 1.0760195666203811E-02
        poly = poly + -9.1242345935631852E-02 * p0(1)
        poly = poly + 1.8372283876790654E-01 * p1(1)
        poly = poly + -3.4374743597940160E-04 * p2(1)
        poly = poly + 1.8960676394141349E-01 * p3(1)
        poly = poly + -2.3230885521567424E-01 * p0(2)
        poly = poly + 1.3377507093192773E-01 * p0(1) * p1(1)
        poly = poly + 6.9301346894818955E-02 * p0(1) * p2(1)
        poly = poly + 3.2635996972925924E-01 * p0(1) * p3(1)
        poly = poly + 6.9526357003141348E-01 * p1(2)
        poly = poly + -1.4896721661872914E-01 * p1(1) * p2(1)
        poly = poly + 1.0279626010452982E+00 * p1(1) * p3(1)
        poly = poly + 5.5058398919769565E-04 * p2(2)
        poly = poly + -1.5167076207673469E-01 * p2(1) * p3(1)
        poly = poly + 3.2729790993170499E-01 * p3(2)
        poly = poly + -1.3861070737590328E-01 * p0(3)
        poly = poly + 1.3154177010248952E-01 * p0(2) * p1(1)
        poly = poly + -1.4182750225343368E-01 * p0(2) * p2(1)
        poly = poly + 3.6072195398232298E-01 * p0(2) * p3(1)
        poly = poly + 2.3129030428920050E-01 * p0(1) * p1(2)
        poly = poly + 3.7289237990420526E-01 * p0(1) * p1(1) * p2(1)
        poly = poly + 2.5000026416325422E-01 * p0(1) * p1(1) * p3(1)
        poly = poly + 1.9012097433876143E-01 * p0(1) * p2(2)
        poly = poly + 1.9351901010487779E-01 * p0(1) * p2(1) * p3(1)
        poly = poly + -1.9527647530937230E-02 * p0(1) * p3(2)
        poly = poly + 1.8508097482261843E-01 * p1(3)
        poly = poly + -1.5991257072436768E-01 * p1(2) * p2(1)
        poly = poly + 9.0822668623141034E-02 * p1(2) * p3(1)
        poly = poly + -3.9148654614840284E-01 * p1(1) * p2(2)
        poly = poly + 6.1554308834023304E-02 * p1(1) * p2(1) * p3(1)
        poly = poly + -3.3783839385985154E-01 * p1(1) * p3(2)
        poly = poly + -1.3014904370288875E-04 * p2(3)
        poly = poly + -4.0262550531955738E-01 * p2(2) * p3(1)
        poly = poly + 2.3176453832257379E-01 * p2(1) * p3(2)
        poly = poly + -2.4243772556729160E-01 * p3(3)
        P = poly

        poly = 1#
        poly = poly + -1.2841492885524408E+00 * p0(1)
        poly = poly + 7.8492438622682070E-02 * p1(1)
        poly = poly + 5.9104732922430375E-02 * p2(1)
        poly = poly + -4.1981950654197026E-01 * p3(1)
        poly = poly + 7.6640876139991054E-01 * p0(2)
        poly = poly + -6.6267515105615460E-01 * p0(1) * p1(1)
        poly = poly + -1.3406101427693026E-01 * p0(1) * p2(1)
        poly = poly + 2.8594952241619392E-01 * p0(1) * p3(1)
        poly = poly + 3.8539562998883459E-01 * p1(2)
        poly = poly + 8.5164463232313567E-02 * p1(1) * p2(1)
        poly = poly + -2.6557322768608504E-01 * p1(1) * p3(1)
        poly = poly + 1.7266080922089984E-02 * p2(2)
        poly = poly + 8.8587142572667560E-02 * p2(1) * p3(1)
        poly = poly + -5.6519486968085608E-01 * p3(2)
        poly = poly + -2.0709263259449648E-01 * p0(3)
        poly = poly + 6.9504768517294568E-01 * p0(2) * p1(1)
        poly = poly + 7.0249243045115470E-02 * p0(2) * p2(1)
        poly = poly + -8.6723279822632671E-02 * p0(2) * p3(1)
        poly = poly + -1.4333632422004592E+00 * p0(1) * p1(2)
        poly = poly + 1.2447854522525711E-01 * p0(1) * p1(1) * p2(1)
        poly = poly + 5.8912900008589366E-03 * p0(1) * p1(1) * p3(1)
        poly = poly + -1.2582941074416365E-01 * p0(1) * p2(2)
        poly = poly + 3.4842344815692408E-03 * p0(1) * p2(1) * p3(1)
        poly = poly + 8.4947550250341702E-02 * p0(1) * p3(2)
        poly = poly + 1.1297114777495512E+00 * p1(3)
        poly = poly + -4.4240437197100363E-01 * p1(2) * p2(1)
        poly = poly + 3.6239115399816446E-01 * p1(2) * p3(1)
        poly = poly + 2.5090086961618407E-01 * p1(1) * p2(2)
        poly = poly + -6.3966196153775734E-01 * p1(1) * p2(1) * p3(1)
        poly = poly + 1.5456754668724705E-02 * p1(1) * p3(2)
        poly = poly + -2.5238498567527429E-04 * p2(3)
        poly = poly + 2.5298188495864032E-01 * p2(2) * p3(1)
        poly = poly + -1.9245799931471816E-01 * p2(1) * p3(2)
        poly = poly + 8.5143385441160258E-01 * p3(3)
        Q = poly

    Case REG_DEEP
        m0 = -2.7445624603199126E+00
        m1 = -1.7486077038578807E-01
        m2 = 8.5264201624293985E-04
        m3 = -2.2199801491625442E+00
        s0 = 5.7036635751105458E-01
        s1 = 4.4079166931112801E-01
        s2 = 1.8056775112058210E-01
        s3 = 1.2795421587003140E+00

        x0 = (x0_raw - m0) / s0
        x1 = (x1_raw - m1) / s1
        x2 = (x2_raw - m2) / s2
        x3 = (x3_raw - m3) / s3

        Call ComputePowers(5, x0, x1, x2, x3, p0, p1, p2, p3)

        poly = 0#
        poly = poly + 3.4426196370822976E-02
        poly = poly + -1.3072897096538516E+00 * p0(1)
        poly = poly + 2.8114838771170039E+00 * p1(1)
        poly = poly + -7.0764737728038667E-03 * p2(1)
        poly = poly + 2.7260650175483843E+00 * p3(1)
        poly = poly + 1.0432414334043347E+00 * p0(2)
        poly = poly + -2.0188322280204134E+00 * p0(1) * p1(1)
        poly = poly + 5.8764547416523849E-01 * p0(1) * p2(1)
        poly = poly + -2.3385213791064121E+00 * p0(1) * p3(1)
        poly = poly + -5.0976598274088802E-01 * p1(2)
        poly = poly + -1.2957319877340383E+00 * p1(1) * p2(1)
        poly = poly + -1.2788929671691907E-01 * p1(1) * p3(1)
        poly = poly + -1.4509003022452308E-03 * p2(2)
        poly = poly + -1.2559431956521041E+00 * p2(1) * p3(1)
        poly = poly + 3.5839402036629792E-01 * p3(2)
        poly = poly + 2.1706478053377504E-01 * p0(3)
        poly = poly + -6.5613733025274368E-01 * p0(2) * p1(1)
        poly = poly + 1.2662013129353633E+00 * p0(2) * p2(1)
        poly = poly + -4.3175782151789943E-01 * p0(2) * p3(1)
        poly = poly + 4.2056052404710093E-01 * p0(1) * p1(2)
        poly = poly + -3.0121545955335036E+00 * p0(1) * p1(1) * p2(1)
        poly = poly + 1.8581852755882183E-01 * p0(1) * p1(1) * p3(1)
        poly = poly + 1.5163874265572602E-01 * p0(1) * p2(2)
        poly = poly + -3.1325606506527159E+00 * p0(1) * p2(1) * p3(1)
        poly = poly + -5.6975507488383946E-02 * p0(1) * p3(2)
        poly = poly + -2.1432508640747333E-01 * p1(3)
        poly = poly + 8.9439959203143867E-04 * p1(2) * p2(1)
        poly = poly + -1.9958700729398529E-01 * p1(2) * p3(1)
        poly = poly + -3.5563940809681532E-01 * p1(1) * p2(2)
        poly = poly + 5.1122684701479837E-01 * p1(1) * p2(1) * p3(1)
        poly = poly + -1.7701113285773568E-01 * p1(1) * p3(2)
        poly = poly + 6.7120281050788279E-04 * p2(3)
        poly = poly + -3.4532035813618939E-01 * p2(2) * p3(1)
        poly = poly + 4.9314916381096680E-01 * p2(1) * p3(2)
        poly = poly + -1.7861634255193468E-01 * p3(3)
        poly = poly + -4.5208380663225139E-01 * p0(4)
        poly = poly + 9.2959043385456552E-01 * p0(3) * p1(1)
        poly = poly + -2.2746316167727146E-01 * p0(3) * p2(1)
        poly = poly + 1.3782127870606189E+00 * p0(3) * p3(1)
        poly = poly + 2.0110524805254043E-01 * p0(2) * p1(2)
        poly = poly + 7.5111781735939276E-01 * p0(2) * p1(1) * p2(1)
        poly = poly + -8.6528007850708732E-01 * p0(2) * p1(1) * p3(1)
        poly = poly + 5.4883523376040344E-01 * p0(2) * p2(2)
        poly = poly + 5.8658928165099011E-01 * p0(2) * p2(1) * p3(1)
        poly = poly + -8.7758854113844509E-01 * p0(2) * p3(2)
        poly = poly + -4.8900423811947605E-02 * p0(1) * p1(3)
        poly = poly + -5.2623367699244927E-01 * p0(1) * p1(2) * p2(1)
        poly = poly + 1.8771707832494319E-01 * p0(1) * p1(2) * p3(1)
        poly = poly + -1.2195259211083069E+00 * p0(1) * p1(1) * p2(2)
        poly = poly + -5.9286383279637300E-01 * p0(1) * p1(1) * p2(1) * p3(1)
        poly = poly + 1.8675992538685626E-01 * p0(1) * p1(1) * p3(2)
        poly = poly + -2.7960638989333344E-02 * p0(1) * p2(3)
        poly = poly + -1.2241319138372908E+00 * p0(1) * p2(2) * p3(1)
        poly = poly + -2.1952488418970723E-01 * p0(1) * p2(1) * p3(2)
        poly = poly + -1.7826251008637806E-03 * p0(1) * p3(3)
        poly = poly + -1.1138104037911873E-02 * p1(4)
        poly = poly + -4.4819442052604847E-03 * p1(3) * p2(1)
        poly = poly + -2.8802443187392456E-02 * p1(3) * p3(1)
        poly = poly + -6.6856186755362715E-02 * p1(2) * p2(2)
        poly = poly + -2.0610220469505397E-01 * p1(2) * p2(1) * p3(1)
        poly = poly + -3.2008209519485700E-02 * p1(2) * p3(2)
        poly = poly + 6.3373980654667678E-02 * p1(1) * p2(3)
        poly = poly + -2.1684387158528124E-02 * p1(1) * p2(2) * p3(1)
        poly = poly + -6.5550444415950856E-02 * p1(1) * p2(1) * p3(2)
        poly = poly + -1.0883929268741739E-01 * p1(1) * p3(3)
        poly = poly + 3.5738115909278038E-04 * p2(4)
        poly = poly + 6.1837498223806610E-02 * p2(3) * p3(1)
        poly = poly + 4.1693635667714994E-02 * p2(2) * p3(2)
        poly = poly + 1.2559649983867738E-01 * p2(1) * p3(3)
        poly = poly + -9.1729872816161517E-02 * p3(4)
        poly = poly + -6.6891407756042812E-02 * p0(5)
        poly = poly + 7.0522831394936744E-02 * p0(4) * p1(1)
        poly = poly + 2.4070800485700181E-01 * p0(4) * p2(1)
        poly = poly + 2.8353214892193518E-01 * p0(4) * p3(1)
        poly = poly + 2.4497357742146608E-01 * p0(3) * p1(2)
        poly = poly + -7.5618286899617071E-01 * p0(3) * p1(1) * p2(1)
        poly = poly + -5.1831709814728809E-02 * p0(3) * p1(1) * p3(1)
        poly = poly + 6.1416163019234360E-01 * p0(3) * p2(2)
        poly = poly + -1.0439949233953567E+00 * p0(3) * p2(1) * p3(1)
        poly = poly + -2.2047666485671136E-01 * p0(3) * p3(2)
        poly = poly + -6.5423390992346037E-02 * p0(2) * p1(3)
        poly = poly + 4.7839870097164733E-01 * p0(2) * p1(2) * p2(1)
        poly = poly + -5.7736202295406526E-02 * p0(2) * p1(2) * p3(1)
        poly = poly + -1.5334703088851322E+00 * p0(2) * p1(1) * p2(2)
        poly = poly + 1.9791081493440581E+00 * p0(2) * p1(1) * p2(1) * p3(1)
        poly = poly + -1.7004123247645145E-01 * p0(2) * p1(1) * p3(2)
        poly = poly + 9.1825984261076993E-02 * p0(2) * p2(3)
        poly = poly + -1.2233638186255120E+00 * p0(2) * p2(2) * p3(1)
        poly = poly + 9.8017320409706699E-01 * p0(2) * p2(1) * p3(2)
        poly = poly + -1.0075702670226144E-01 * p0(2) * p3(3)
        poly = poly + 4.6180556019469639E-02 * p0(1) * p1(4)
        poly = poly + 1.5462856182672940E-01 * p0(1) * p1(3) * p2(1)
        poly = poly + 2.5980224621631038E-02 * p0(1) * p1(3) * p3(1)
        poly = poly + 1.9711641736495114E-01 * p0(1) * p1(2) * p2(2)
        poly = poly + -4.0180605935095781E-01 * p0(1) * p1(2) * p2(1) * p3(1)
        poly = poly + -3.2066775691250209E-02 * p0(1) * p1(2) * p3(2)
        poly = poly + -2.0847722257271398E-01 * p0(1) * p1(1) * p2(3)
        poly = poly + -2.0395114079743193E-01 * p0(1) * p1(1) * p2(2) * p3(1)
        poly = poly + -3.4920586889726969E-02 * p0(1) * p1(1) * p2(1) * p3(2)
        poly = poly + -3.4324140127305267E-02 * p0(1) * p1(1) * p3(3)
        poly = poly + -1.0269810587309119E-02 * p0(1) * p2(4)
        poly = poly + -2.0595698734770188E-01 * p0(1) * p2(3) * p3(1)
        poly = poly + -3.8586989982889225E-01 * p0(1) * p2(2) * p3(2)
        poly = poly + 4.9205363282828174E-01 * p0(1) * p2(1) * p3(3)
        poly = poly + -6.9424923124591420E-02 * p0(1) * p3(4)
        poly = poly + 1.2978710752225082E-02 * p1(5)
        poly = poly + -4.1666439808865072E-02 * p1(4) * p2(1)
        poly = poly + -1.5302307963667744E-02 * p1(4) * p3(1)
        poly = poly + 3.5681466486893562E-02 * p1(3) * p2(2)
        poly = poly + -5.6204955636561603E-02 * p1(3) * p2(1) * p3(1)
        poly = poly + -4.1220206190559243E-03 * p1(3) * p3(2)
        poly = poly + -1.1519586923860837E-02 * p1(2) * p2(3)
        poly = poly + 5.3831571909363683E-02 * p1(2) * p2(2) * p3(1)
        poly = poly + -3.4746386419855536E-02 * p1(2) * p2(1) * p3(2)
        poly = poly + -2.3807887860447099E-02 * p1(2) * p3(3)
        poly = poly + 2.2896718900033398E-02 * p1(1) * p2(4)
        poly = poly + -1.3897723270288528E-02 * p1(1) * p2(3) * p3(1)
        poly = poly + 1.2925401924570562E-02 * p1(1) * p2(2) * p3(2)
        poly = poly + -8.2707317994064541E-02 * p1(1) * p2(1) * p3(3)
        poly = poly + -3.2402338705799338E-02 * p1(1) * p3(4)
        poly = poly + 9.0238041343750525E-06 * p2(5)
        poly = poly + 2.2229626163170185E-02 * p2(4) * p3(1)
        poly = poly + -2.5583593582847470E-03 * p2(3) * p3(2)
        poly = poly + -5.5431525135838074E-03 * p2(2) * p3(3)
        poly = poly + -6.1880460828311543E-02 * p2(1) * p3(4)
        poly = poly + 1.1625139299703095E-02 * p3(5)
        P = poly

        poly = 1#
        poly = poly + -3.2079819723474947E+00 * p0(1)
        poly = poly + -7.2432056213021356E-01 * p1(1)
        poly = poly + 4.7724903186166540E-02 * p2(1)
        poly = poly + -6.8141684485558252E-01 * p3(1)
        poly = poly + 4.3761433621753891E+00 * p0(2)
        poly = poly + 9.7640083152636203E-01 * p0(1) * p1(1)
        poly = poly + -1.4532356502233956E-01 * p0(1) * p2(1)
        poly = poly + 9.4161562930551523E-01 * p0(1) * p3(1)
        poly = poly + 1.8285205738686719E-01 * p1(2)
        poly = poly + -1.0371503843862819E-02 * p1(1) * p2(1)
        poly = poly + 2.4626558399216045E-01 * p1(1) * p3(1)
        poly = poly + -1.6727083223089419E-02 * p2(2)
        poly = poly + -5.4060468995810046E-02 * p2(1) * p3(1)
        poly = poly + 1.6522348590561511E-01 * p3(2)
        poly = poly + -2.9777403025235110E+00 * p0(3)
        poly = poly + -6.3671702265939334E-01 * p0(2) * p1(1)
        poly = poly + 2.5248891870281737E-01 * p0(2) * p2(1)
        poly = poly + -6.6932081862402515E-01 * p0(2) * p3(1)
        poly = poly + -1.5699286790983411E-01 * p0(1) * p1(2)
        poly = poly + -1.6283962161437859E-04 * p0(1) * p1(1) * p2(1)
        poly = poly + -1.2157358274335248E-01 * p0(1) * p1(1) * p3(1)
        poly = poly + -6.2820731647053110E-02 * p0(1) * p2(2)
        poly = poly + 1.1271729410272635E-01 * p0(1) * p2(1) * p3(1)
        poly = poly + -1.7271194643092833E-01 * p0(1) * p3(2)
        poly = poly + -4.3096166287343946E-02 * p1(3)
        poly = poly + -1.8295016391849939E-02 * p1(2) * p2(1)
        poly = poly + -2.5442038862629787E-02 * p1(2) * p3(1)
        poly = poly + -1.1376870718589595E-02 * p1(1) * p2(2)
        poly = poly + 1.8834856042727284E-02 * p1(1) * p2(1) * p3(1)
        poly = poly + -2.7898702150714421E-02 * p1(1) * p3(2)
        poly = poly + 4.2544153430687383E-02 * p2(3)
        poly = poly + -1.6245145736252221E-02 * p2(2) * p3(1)
        poly = poly + 3.0779292918828204E-02 * p2(1) * p3(2)
        poly = poly + -4.8155139284323523E-02 * p3(3)
        poly = poly + 1.0290134162965947E+00 * p0(4)
        poly = poly + 1.9396411911036895E-01 * p0(3) * p1(1)
        poly = poly + -2.2088535491162770E-01 * p0(3) * p2(1)
        poly = poly + 2.5823402819220148E-01 * p0(3) * p3(1)
        poly = poly + 7.1672587155542811E-02 * p0(2) * p1(2)
        poly = poly + 8.9965521406316987E-03 * p0(2) * p1(1) * p2(1)
        poly = poly + 1.2389255403064844E-02 * p0(2) * p1(1) * p3(1)
        poly = poly + 1.7880570347956204E-01 * p0(2) * p2(2)
        poly = poly + -1.0775915739198381E-01 * p0(2) * p2(1) * p3(1)
        poly = poly + 1.0230594032405126E-01 * p0(2) * p3(2)
        poly = poly + 3.0884344782007900E-02 * p0(1) * p1(3)
        poly = poly + 1.9868995551154709E-02 * p0(1) * p1(2) * p2(1)
        poly = poly + 3.0436322382000216E-05 * p0(1) * p1(2) * p3(1)
        poly = poly + 3.9453051834876372E-02 * p0(1) * p1(1) * p2(2)
        poly = poly + -1.6523828080046418E-02 * p0(1) * p1(1) * p2(1) * p3(1)
        poly = poly + 5.4911538652827841E-03 * p0(1) * p1(1) * p3(2)
        poly = poly + -1.0186339725681588E-01 * p0(1) * p2(3)
        poly = poly + 3.8930349452445938E-02 * p0(1) * p2(2) * p3(1)
        poly = poly + -3.0957558918352800E-02 * p0(1) * p2(1) * p3(2)
        poly = poly + 3.9928738280452493E-02 * p0(1) * p3(3)
        poly = poly + 8.3621687052928811E-03 * p1(4)
        poly = poly + 6.1379437964106586E-03 * p1(3) * p2(1)
        poly = poly + 5.1248543500176866E-03 * p1(3) * p3(1)
        poly = poly + 6.7814069572781676E-03 * p1(2) * p2(2)
        poly = poly + 2.5133537540195042E-03 * p1(2) * p2(1) * p3(1)
        poly = poly + -5.2828396184932290E-03 * p1(2) * p3(2)
        poly = poly + -2.5459679593095832E-02 * p1(1) * p2(3)
        poly = poly + 1.0578100706329993E-02 * p1(1) * p2(2) * p3(1)
        poly = poly + -9.9631136140781850E-03 * p1(1) * p2(1) * p3(2)
        poly = poly + 7.9074155939011098E-03 * p1(1) * p3(3)
        poly = poly + 1.1019866632944009E-02 * p2(4)
        poly = poly + -1.9094496631405108E-02 * p2(3) * p3(1)
        poly = poly + 6.4213140799616068E-03 * p2(2) * p3(2)
        poly = poly + -3.5029549600381748E-03 * p2(1) * p3(3)
        poly = poly + 9.6264472352404797E-03 * p3(4)
        poly = poly + -1.3961293586539816E-01 * p0(5)
        poly = poly + -1.6122751973839070E-02 * p0(4) * p1(1)
        poly = poly + 6.5758273855018826E-02 * p0(4) * p2(1)
        poly = poly + -4.5571149931040847E-02 * p0(4) * p3(1)
        poly = poly + -1.1659622998794661E-02 * p0(3) * p1(2)
        poly = poly + -6.5765820494998262E-03 * p0(3) * p1(1) * p2(1)
        poly = poly + 4.8630815660189331E-03 * p0(3) * p1(1) * p3(1)
        poly = poly + -1.0036373145077002E-01 * p0(3) * p2(2)
        poly = poly + 3.6109024212485630E-02 * p0(3) * p2(1) * p3(1)
        poly = poly + -2.5339589306671070E-02 * p0(3) * p3(2)
        poly = poly + -7.3575365870446312E-03 * p0(2) * p1(3)
        poly = poly + -7.9857389607478658E-03 * p0(2) * p1(2) * p2(1)
        poly = poly + 2.4064565119096869E-03 * p0(2) * p1(2) * p3(1)
        poly = poly + -2.3863287483586657E-02 * p0(2) * p1(1) * p2(2)
        poly = poly + 5.3215005101969158E-03 * p0(2) * p1(1) * p2(1) * p3(1)
        poly = poly + -3.1925152595507143E-04 * p0(2) * p1(1) * p3(2)
        poly = poly + 6.3547799238911337E-02 * p0(2) * p2(3)
        poly = poly + -2.0075842002571474E-02 * p0(2) * p2(2) * p3(1)
        poly = poly + 1.0596245891011255E-02 * p0(2) * p2(1) * p3(2)
        poly = poly + -1.0965388446612573E-02 * p0(2) * p3(3)
        poly = poly + -2.4897251686087037E-03 * p0(1) * p1(4)
        poly = poly + -2.2840471495929843E-03 * p0(1) * p1(3) * p2(1)
        poly = poly + -7.0661287854459039E-04 * p0(1) * p1(3) * p3(1)
        poly = poly + -5.0738467505963688E-03 * p0(1) * p1(2) * p2(2)
        poly = poly + -1.1991970570974016E-03 * p0(1) * p1(2) * p2(1) * p3(1)
        poly = poly + 1.8029672089043777E-03 * p0(1) * p1(2) * p3(2)
        poly = poly + 1.6750923759740879E-02 * p0(1) * p1(1) * p2(3)
        poly = poly + -5.3935563519279840E-03 * p0(1) * p1(1) * p2(2) * p3(1)
        poly = poly + 3.6114532092418180E-03 * p0(1) * p1(1) * p2(1) * p3(2)
        poly = poly + -2.0056280554990955E-03 * p0(1) * p1(1) * p3(3)
        poly = poly + -1.1500964920107618E-02 * p0(1) * p2(4)
        poly = poly + 1.1015286902707481E-02 * p0(1) * p2(3) * p3(1)
        poly = poly + -3.3748620549211013E-03 * p0(1) * p2(2) * p3(2)
        poly = poly + 9.9102365746751592E-04 * p0(1) * p2(1) * p3(3)
        poly = poly + -2.8151505134228383E-03 * p0(1) * p3(4)
        poly = poly + -2.2648924890110786E-04 * p1(5)
        poly = poly + 7.0944459020728034E-05 * p1(4) * p2(1)
        poly = poly + -8.7574564604638342E-04 * p1(4) * p3(1)
        poly = poly + -9.7088217977635489E-04 * p1(3) * p2(2)
        poly = poly + -1.0914563456732977E-03 * p1(3) * p2(1) * p3(1)
        poly = poly + 5.9008487104897432E-04 * p1(3) * p3(2)
        poly = poly + 3.7506868954991675E-03 * p1(2) * p2(3)
        poly = poly + -1.2583290807035138E-03 * p1(2) * p2(2) * p3(1)
        poly = poly + 5.9343586031130994E-04 * p1(2) * p2(1) * p3(2)
        poly = poly + 1.9384912256687338E-04 * p1(2) * p3(3)
        poly = poly + -2.9804385744005655E-03 * p1(1) * p2(4)
        poly = poly + 3.5906195176162990E-03 * p1(1) * p2(3) * p3(1)
        poly = poly + -1.1037780828929606E-03 * p1(1) * p2(2) * p3(2)
        poly = poly + 9.9652848649640045E-04 * p1(1) * p2(1) * p3(3)
        poly = poly + -1.0943650636456694E-03 * p1(1) * p3(4)
        poly = poly + 7.4102322059012900E-04 * p2(5)
        poly = poly + -2.0464399220859134E-03 * p2(4) * p3(1)
        poly = poly + 1.1993387052404672E-03 * p2(3) * p3(2)
        poly = poly + -3.6364661589748429E-04 * p2(2) * p3(3)
        poly = poly + -5.8813465119586365E-04 * p2(1) * p3(4)
        poly = poly + -1.2387611635104035E-04 * p3(5)
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

    PredictRegime = MaxVal(L_base, 0.1#) * Exp(log_corr)
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
