# SOLVER AUDIT REPORT
# Case: Test wave
# Solution by fixed 50-mode spectral collocation
# This file records the complete nonlinear state, complete residual vector, continuation history, and independent final convergence checks.

# Solver configuration
Fourier order N                         50
Nonlinear unknowns                     110
Nonlinear equations                    110
Intermediate residual tolerance        1.0000000000000000e-10
Final residual RMS tolerance           9.9999999999999998e-13
Final residual max tolerance           9.9999999999999994e-12
Off-grid spectral defect tolerance     1.0000000000000001e-09
Target H/d                             5.9999999999999998e-01
Current definition                     Eulerian
Prescribed current / sqrt(gd)          1.4280869812290301e-01

# Adaptive continuation history
# trial       fraction       increment status        residual_rms       residual_max    next_increment tangent
      1  7.692307692308e-02  7.692307692308e-02 ACCEPT    5.222078424032e-13  1.219153332076e-12  1.153846153846e-01 YES
      2  1.923076923077e-01  1.153846153846e-01 ACCEPT    9.410304270164e-12  5.698222449446e-11  1.153846153846e-01 YES
      3  3.076923076923e-01  1.153846153846e-01 ACCEPT    8.065271988354e-12  4.883765462513e-11  1.153846153846e-01 YES
      4  4.230769230769e-01  1.153846153846e-01 ACCEPT    1.680845475431e-14  1.017860277173e-13  1.153846153846e-01 YES
      5  5.384615384615e-01  1.153846153846e-01 ACCEPT    2.678203813846e-16  1.422684652913e-15  1.153846153846e-01 YES
      6  6.538461538462e-01  1.153846153846e-01 ACCEPT    6.273994883589e-17  5.412337245048e-16  1.153846153846e-01 YES
      7  7.692307692308e-01  1.153846153846e-01 ACCEPT    1.847860371264e-16  1.443289932013e-15  1.153846153846e-01 YES
      8  8.846153846154e-01  1.153846153846e-01 ACCEPT    2.502004560702e-16  2.609024107869e-15  1.153846153846e-01 YES
      9  1.000000000000e+00  1.153846153846e-01 ACCEPT    1.378636209957e-16  1.054711873394e-15  1.153846153846e-01 YES
# Accepted continuation steps: 9
# Total continuation trials:   9

# Independent final convergence proof
Target continuation fraction reached       PASS
Final state admissible                  PASS
All final residuals finite              PASS
All final Jacobian entries finite       PASS
Final residual sum of squares           2.0907015793435818e-30
Final residual L2 norm                  1.4459258554101526e-15
Final residual RMS                      1.3786362099565123e-16
Required residual RMS <=                9.9999999999999998e-13  PASS
Final maximum absolute residual         1.0547118733938987e-15
Required residual max <=                9.9999999999999994e-12  PASS
Index of maximum absolute residual      F[7]
Equation with maximum residual          trapezoidal zero-mean free-surface constraint
Combined residual convergence test      PASS
Final scaled Newton correction solved   PASS
Final scaled Newton correction RMS      3.0176226864166552e-16
Final linear solver                     LU
Maximum between-node spectral defect    1.9822312275754803e-11
Spectral resolution check               PASS
Overall nonlinear convergence            PASS
# The nonlinear convergence result is based on an independent re-evaluation of all 110 equations at the final state. The spectral defect is a separate truncation-resolution check.

# Residual groups
Global constraints                 F[  1..  8] count=  8 RMS= 5.0658939237266969e-16 MAX= 1.0547118733938987e-15 at F[7]
Kinematic surface equations        F[  9.. 59] count= 51 RMS= 7.8658926012321446e-18 MAX= 2.4936649967166602e-17 at F[10]
Dynamic surface equations          F[ 60..110] count= 51 RMS= 2.6002863356801221e-17 MAX= 4.3340981845108040e-17 at F[86]

# Complete final nonlinear state vector z[1..110]
# index              value  variable and definition
z[  1] =  3.98541591251242111e-01  kd: mean water depth scaled by wavenumber
z[  2] =  2.39124954750745300e-01  kH: crest-to-trough wave height scaled by wavenumber
z[  3] =  7.95709732959653326e+00  T*sqrt(gk): apparent wave period in native scaling
z[  4] =  7.89632833044430771e-01  c*sqrt(k/g): wave celerity in native scaling
z[  5] =  9.01553461692551628e-02  u_E*sqrt(k/g): Eulerian-mean current
z[  6] =  1.07619351116591974e-01  u_M*sqrt(k/g): mass-transport mean current
z[  7] =  6.99477486875175525e-01  U*sqrt(k/g): mean fluid speed relative to the wave
z[  8] =  6.96013232133114567e-03  q*sqrt(k^3/g): wave-volume-flux variable
z[  9] =  2.47888342046438698e-01  r*k/g: Bernoulli constant
z[ 10] =  1.98381053610311636e-01  eta[0]: free-surface ordinate at X=0*pi/N
z[ 11] =  1.92785807180771374e-01  eta[1]: free-surface ordinate at X=1*pi/N
z[ 12] =  1.78140347346803118e-01  eta[2]: free-surface ordinate at X=2*pi/N
z[ 13] =  1.58579885393630549e-01  eta[3]: free-surface ordinate at X=3*pi/N
z[ 14] =  1.37373701960425421e-01  eta[4]: free-surface ordinate at X=4*pi/N
z[ 15] =  1.16409932920982673e-01  eta[5]: free-surface ordinate at X=5*pi/N
z[ 16] =  9.66456365930722372e-02  eta[6]: free-surface ordinate at X=6*pi/N
z[ 17] =  7.85198302308441204e-02  eta[7]: free-surface ordinate at X=7*pi/N
z[ 18] =  6.21931409082371670e-02  eta[8]: free-surface ordinate at X=8*pi/N
z[ 19] =  4.76750693504148157e-02  eta[9]: free-surface ordinate at X=9*pi/N
z[ 20] =  3.48922639147840191e-02  eta[10]: free-surface ordinate at X=10*pi/N
z[ 21] =  2.37267819562240641e-02  eta[11]: free-surface ordinate at X=11*pi/N
z[ 22] =  1.40386521146659126e-02  eta[12]: free-surface ordinate at X=12*pi/N
z[ 23] =  5.67976123401355611e-03  eta[13]: free-surface ordinate at X=13*pi/N
z[ 24] = -1.49739077953406764e-03  eta[14]: free-surface ordinate at X=14*pi/N
z[ 25] = -7.63418701439231608e-03  eta[15]: free-surface ordinate at X=15*pi/N
z[ 26] = -1.28625664375079472e-02  eta[16]: free-surface ordinate at X=16*pi/N
z[ 27] = -1.73032045743832623e-02  eta[17]: free-surface ordinate at X=17*pi/N
z[ 28] = -2.10647489873358185e-02  eta[18]: free-surface ordinate at X=18*pi/N
z[ 29] = -2.42437833732603476e-02  eta[19]: free-surface ordinate at X=19*pi/N
z[ 30] = -2.69252777338798675e-02  eta[20]: free-surface ordinate at X=20*pi/N
z[ 31] = -2.91833431052354042e-02  eta[21]: free-surface ordinate at X=21*pi/N
z[ 32] = -3.10821566742463980e-02  eta[22]: free-surface ordinate at X=22*pi/N
z[ 33] = -3.26769606130282650e-02  eta[23]: free-surface ordinate at X=23*pi/N
z[ 34] = -3.40150675094438373e-02  eta[24]: free-surface ordinate at X=24*pi/N
z[ 35] = -3.51368280943858413e-02  eta[25]: free-surface ordinate at X=25*pi/N
z[ 36] = -3.60765341196268643e-02  eta[26]: free-surface ordinate at X=26*pi/N
z[ 37] = -3.68632417301404947e-02  eta[27]: free-surface ordinate at X=27*pi/N
z[ 38] = -3.75215094461581952e-02  eta[28]: free-surface ordinate at X=28*pi/N
z[ 39] = -3.80720507700349278e-02  eta[29]: free-surface ordinate at X=29*pi/N
z[ 40] = -3.85323051844824044e-02  eta[30]: free-surface ordinate at X=30*pi/N
z[ 41] = -3.89169335026839508e-02  eta[31]: free-surface ordinate at X=31*pi/N
z[ 42] = -3.92382446318091868e-02  eta[32]: free-surface ordinate at X=32*pi/N
z[ 43] = -3.95065611729188138e-02  eta[33]: free-surface ordinate at X=33*pi/N
z[ 44] = -3.97305311631543923e-02  eta[34]: free-surface ordinate at X=34*pi/N
z[ 45] = -3.99173928577966491e-02  eta[35]: free-surface ordinate at X=35*pi/N
z[ 46] = -4.00731988810002979e-02  eta[36]: free-surface ordinate at X=36*pi/N
z[ 47] = -4.02030054332456269e-02  eta[37]: free-surface ordinate at X=37*pi/N
z[ 48] = -4.03110315884811499e-02  eta[38]: free-surface ordinate at X=38*pi/N
z[ 49] = -4.04007930803730944e-02  eta[39]: free-surface ordinate at X=39*pi/N
z[ 50] = -4.04752143856329516e-02  eta[40]: free-surface ordinate at X=40*pi/N
z[ 51] = -4.05367223733776871e-02  eta[41]: free-surface ordinate at X=41*pi/N
z[ 52] = -4.05873243064315131e-02  eta[42]: free-surface ordinate at X=42*pi/N
z[ 53] = -4.06286725526570705e-02  eta[43]: free-surface ordinate at X=43*pi/N
z[ 54] = -4.06621179884786230e-02  eta[44]: free-surface ordinate at X=44*pi/N
z[ 55] = -4.06887537479668146e-02  eta[45]: free-surface ordinate at X=45*pi/N
z[ 56] = -4.07094506837530415e-02  eta[46]: free-surface ordinate at X=46*pi/N
z[ 57] = -4.07248856549784002e-02  eta[47]: free-surface ordinate at X=47*pi/N
z[ 58] = -4.07355635368796890e-02  eta[48]: free-surface ordinate at X=48*pi/N
z[ 59] = -4.07418336511125578e-02  eta[49]: free-surface ordinate at X=49*pi/N
z[ 60] = -4.07439011404336290e-02  eta[50]: free-surface ordinate at X=50*pi/N
z[ 61] =  1.05067973627037167e-01  B[1]: stream-function Fourier coefficient
z[ 62] =  3.61089684232196034e-02  B[2]: stream-function Fourier coefficient
z[ 63] =  1.41109518426180448e-02  B[3]: stream-function Fourier coefficient
z[ 64] =  5.50838843096240587e-03  B[4]: stream-function Fourier coefficient
z[ 65] =  2.06115220910717224e-03  B[5]: stream-function Fourier coefficient
z[ 66] =  7.20957606320987276e-04  B[6]: stream-function Fourier coefficient
z[ 67] =  2.30648236511205116e-04  B[7]: stream-function Fourier coefficient
z[ 68] =  6.61227756266817660e-05  B[8]: stream-function Fourier coefficient
z[ 69] =  1.70224193680355176e-05  B[9]: stream-function Fourier coefficient
z[ 70] =  4.54672633026073837e-06  B[10]: stream-function Fourier coefficient
z[ 71] =  1.95218288227848842e-06  B[11]: stream-function Fourier coefficient
z[ 72] =  1.38191405574222841e-06  B[12]: stream-function Fourier coefficient
z[ 73] =  1.03309247504867816e-06  B[13]: stream-function Fourier coefficient
z[ 74] =  6.92804516938986030e-07  B[14]: stream-function Fourier coefficient
z[ 75] =  4.13759723861588701e-07  B[15]: stream-function Fourier coefficient
z[ 76] =  2.24193087131496623e-07  B[16]: stream-function Fourier coefficient
z[ 77] =  1.12406979118317915e-07  B[17]: stream-function Fourier coefficient
z[ 78] =  5.32216073445646407e-08  B[18]: stream-function Fourier coefficient
z[ 79] =  2.43906438491294681e-08  B[19]: stream-function Fourier coefficient
z[ 80] =  1.11693069588464767e-08  B[20]: stream-function Fourier coefficient
z[ 81] =  5.29549289421307377e-09  B[21]: stream-function Fourier coefficient
z[ 82] =  2.66625514776606305e-09  B[22]: stream-function Fourier coefficient
z[ 83] =  1.42672744807999562e-09  B[23]: stream-function Fourier coefficient
z[ 84] =  7.93903752686030888e-10  B[24]: stream-function Fourier coefficient
z[ 85] =  4.46831766826857594e-10  B[25]: stream-function Fourier coefficient
z[ 86] =  2.49185854643780601e-10  B[26]: stream-function Fourier coefficient
z[ 87] =  1.36280931264682798e-10  B[27]: stream-function Fourier coefficient
z[ 88] =  7.29564371286465321e-11  B[28]: stream-function Fourier coefficient
z[ 89] =  3.83796209934269338e-11  B[29]: stream-function Fourier coefficient
z[ 90] =  1.99885325311793773e-11  B[30]: stream-function Fourier coefficient
z[ 91] =  1.03987694435806905e-11  B[31]: stream-function Fourier coefficient
z[ 92] =  5.44813491510778719e-12  B[32]: stream-function Fourier coefficient
z[ 93] =  2.88973595416260614e-12  B[33]: stream-function Fourier coefficient
z[ 94] =  1.55348743817930694e-12  B[34]: stream-function Fourier coefficient
z[ 95] =  8.44250242138076212e-13  B[35]: stream-function Fourier coefficient
z[ 96] =  4.61661314691138055e-13  B[36]: stream-function Fourier coefficient
z[ 97] =  2.52805749531635616e-13  B[37]: stream-function Fourier coefficient
z[ 98] =  1.38145542219812019e-13  B[38]: stream-function Fourier coefficient
z[ 99] =  7.52072588925628973e-14  B[39]: stream-function Fourier coefficient
z[100] =  4.07949561863947681e-14  B[40]: stream-function Fourier coefficient
z[101] =  2.20831044161123785e-14  B[41]: stream-function Fourier coefficient
z[102] =  1.19555605022286690e-14  B[42]: stream-function Fourier coefficient
z[103] =  6.48981105166934509e-15  B[43]: stream-function Fourier coefficient
z[104] =  3.53965510031210346e-15  B[44]: stream-function Fourier coefficient
z[105] =  1.94444249060438622e-15  B[45]: stream-function Fourier coefficient
z[106] =  1.07799296662741714e-15  B[46]: stream-function Fourier coefficient
z[107] =  6.06371006882315942e-16  B[47]: stream-function Fourier coefficient
z[108] =  3.51294591518102954e-16  B[48]: stream-function Fourier coefficient
z[109] =  2.23817690067850842e-16  B[49]: stream-function Fourier coefficient
z[110] =  9.32859244414695949e-17  B[50]: stream-function Fourier coefficient

# Complete final residual vector F[1..110]
# index             residual       absolute_value  fraction_of_max_limit  equation
F[  1] =  4.21560917563122331e-17   4.21560917563122331e-17   4.215609176e-06  kH - kd*(H/d)
F[  2] =  1.06785044722761168e-16   1.06785044722761168e-16   1.067850447e-05  kH - height_parameter*T^2
F[  3] = -9.56010382435361781e-16   9.56010382435361781e-16   9.560103824e-05  c*T - 2*pi
F[  4] = -1.11022302462515654e-16   1.11022302462515654e-16   1.110223025e-05  u_E + U - c
F[  5] = -1.79444234892111576e-17   1.79444234892111576e-17   1.794442349e-06  kd*(u_M + U - c) - q
F[  6] =  9.53963883014665368e-18   9.53963883014665368e-18   9.539638830e-07  u_E - prescribed_current*sqrt(kd)
F[  7] = -1.05471187339389871e-15   1.05471187339389871e-15   1.054711873e-04  trapezoidal zero-mean free-surface constraint
F[  8] = -2.77555756156289135e-17   2.77555756156289135e-17   2.775557562e-06  eta[0] - eta[N] - kH
F[  9] =  2.01119502996061073e-17   2.01119502996061073e-17   2.011195030e-06  kinematic free-surface equation at node m=0
F[ 10] =  2.49366499671666020e-17   2.49366499671666020e-17   2.493664997e-06  kinematic free-surface equation at node m=1
F[ 11] =  1.52330405234213373e-17   1.52330405234213373e-17   1.523304052e-06  kinematic free-surface equation at node m=2
F[ 12] = -7.86046575051990715e-18   7.86046575051990715e-18   7.860465751e-07  kinematic free-surface equation at node m=3
F[ 13] =  2.67662411332358907e-18   2.67662411332358907e-18   2.676624113e-07  kinematic free-surface equation at node m=4
F[ 14] =  3.54398585131199262e-18   3.54398585131199262e-18   3.543985851e-07  kinematic free-surface equation at node m=5
F[ 15] =  1.93462325152882197e-17   1.93462325152882197e-17   1.934623252e-06  kinematic free-surface equation at node m=6
F[ 16] =  1.04354459101729802e-18   1.04354459101729802e-18   1.043544591e-07  kinematic free-surface equation at node m=7
F[ 17] = -9.21571846612678769e-19   9.21571846612678769e-19   9.215718466e-08  kinematic free-surface equation at node m=8
F[ 18] =  7.15573433840432926e-18   7.15573433840432926e-18   7.155734338e-07  kinematic free-surface equation at node m=9
F[ 19] =  1.67255125764834145e-17   1.67255125764834145e-17   1.672551258e-06  kinematic free-surface equation at node m=10
F[ 20] = -6.20705743747951288e-18   6.20705743747951288e-18   6.207057437e-07  kinematic free-surface equation at node m=11
F[ 21] =  5.91567810362403357e-18   5.91567810362403357e-18   5.915678104e-07  kinematic free-surface equation at node m=12
F[ 22] =  1.35173752887577520e-17   1.35173752887577520e-17   1.351737529e-06  kinematic free-surface equation at node m=13
F[ 23] =  3.50248123689653190e-18   3.50248123689653190e-18   3.502481237e-07  kinematic free-surface equation at node m=14
F[ 24] =  3.13952761899806421e-18   3.13952761899806421e-18   3.139527619e-07  kinematic free-surface equation at node m=15
F[ 25] =  4.43760561066527948e-18   4.43760561066527948e-18   4.437605611e-07  kinematic free-surface equation at node m=16
F[ 26] =  5.84113920426565514e-18   5.84113920426565514e-18   5.841139204e-07  kinematic free-surface equation at node m=17
F[ 27] = -1.99984478846740310e-18   1.99984478846740310e-18   1.999844788e-07  kinematic free-surface equation at node m=18
F[ 28] =  2.42759642683082477e-18   2.42759642683082477e-18   2.427596427e-07  kinematic free-surface equation at node m=19
F[ 29] = -1.43995601033231058e-19   1.43995601033231058e-19   1.439956010e-08  kinematic free-surface equation at node m=20
F[ 30] = -1.69406589450860068e-18   1.69406589450860068e-18   1.694065895e-07  kinematic free-surface equation at node m=21
F[ 31] =  2.49535906261116880e-18   2.49535906261116880e-18   2.495359063e-07  kinematic free-surface equation at node m=22
F[ 32] =  3.32545135092038313e-18   3.32545135092038313e-18   3.325451351e-07  kinematic free-surface equation at node m=23
F[ 33] = -3.04931861011548122e-18   3.04931861011548122e-18   3.049318610e-07  kinematic free-surface equation at node m=24
F[ 34] = -4.40457132572236176e-20   4.40457132572236176e-20   4.404571326e-09  kinematic free-surface equation at node m=25
F[ 35] =  5.53959547504312422e-19   5.53959547504312422e-19   5.539595475e-08  kinematic free-surface equation at node m=26
F[ 36] =  5.87671458805033575e-18   5.87671458805033575e-18   5.876714588e-07  kinematic free-surface equation at node m=27
F[ 37] =  3.28648783534668532e-19   3.28648783534668532e-19   3.286487835e-08  kinematic free-surface equation at node m=28
F[ 38] = -2.37169225231204095e-19   2.37169225231204095e-19   2.371692252e-08  kinematic free-surface equation at node m=29
F[ 39] = -2.58853268680914184e-18   2.58853268680914184e-18   2.588532687e-07  kinematic free-surface equation at node m=30
F[ 40] = -2.99002630380768020e-18   2.99002630380768020e-18   2.990026304e-07  kinematic free-surface equation at node m=31
F[ 41] = -3.18314981578166067e-18   3.18314981578166067e-18   3.183149816e-07  kinematic free-surface equation at node m=32
F[ 42] = -7.53181696698523861e-18   7.53181696698523861e-18   7.531816967e-07  kinematic free-surface equation at node m=33
F[ 43] = -4.05220561966457282e-18   4.05220561966457282e-18   4.052205620e-07  kinematic free-surface equation at node m=34
F[ 44] =  1.20448085099561508e-18   1.20448085099561508e-18   1.204480851e-07  kinematic free-surface equation at node m=35
F[ 45] = -8.30092288309214332e-20   8.30092288309214332e-20   8.300922883e-09  kinematic free-surface equation at node m=36
F[ 46] = -4.03187682893046961e-19   4.03187682893046961e-19   4.031876829e-08  kinematic free-surface equation at node m=37
F[ 47] = -2.64443686132792566e-18   2.64443686132792566e-18   2.644436861e-07  kinematic free-surface equation at node m=38
F[ 48] =  4.77895988840876251e-18   4.77895988840876251e-18   4.778959888e-07  kinematic free-surface equation at node m=39
F[ 49] =  3.41523684332933897e-18   3.41523684332933897e-18   3.415236843e-07  kinematic free-surface equation at node m=40
F[ 50] = -8.80067232197218052e-18   8.80067232197218052e-18   8.800672322e-07  kinematic free-surface equation at node m=41
F[ 51] =  3.55076211489002702e-18   3.55076211489002702e-18   3.550762115e-07  kinematic free-surface equation at node m=42
F[ 52] = -3.56939683972962163e-18   3.56939683972962163e-18   3.569396840e-07  kinematic free-surface equation at node m=43
F[ 53] = -1.86516654985396935e-18   1.86516654985396935e-18   1.865166550e-07  kinematic free-surface equation at node m=44
F[ 54] = -6.55942314353730183e-18   6.55942314353730183e-18   6.559423144e-07  kinematic free-surface equation at node m=45
F[ 55] =  1.42470941728173317e-18   1.42470941728173317e-18   1.424709417e-07  kinematic free-surface equation at node m=46
F[ 56] = -9.23265912507187370e-18   9.23265912507187370e-18   9.232659125e-07  kinematic free-surface equation at node m=47
F[ 57] = -7.54028729645778162e-18   7.54028729645778162e-18   7.540287296e-07  kinematic free-surface equation at node m=48
F[ 58] = -1.22802836692928463e-17   1.22802836692928463e-17   1.228028367e-06  kinematic free-surface equation at node m=49
F[ 59] = -1.07522362324460885e-17   1.07522362324460885e-17   1.075223623e-06  kinematic free-surface equation at node m=50
F[ 60] = -1.26038502551439890e-18   1.26038502551439890e-18   1.260385026e-07  dynamic Bernoulli free-surface equation at node m=0
F[ 61] =  2.03152382069471393e-17   2.03152382069471393e-17   2.031523821e-06  dynamic Bernoulli free-surface equation at node m=1
F[ 62] = -1.43792313125890026e-17   1.43792313125890026e-17   1.437923131e-06  dynamic Bernoulli free-surface equation at node m=2
F[ 63] = -1.00966327312712600e-17   1.00966327312712600e-17   1.009663273e-06  dynamic Bernoulli free-surface equation at node m=3
F[ 64] = -1.94207714146465982e-17   1.94207714146465982e-17   1.942077141e-06  dynamic Bernoulli free-surface equation at node m=4
F[ 65] = -1.42166009867161769e-17   1.42166009867161769e-17   1.421660099e-06  dynamic Bernoulli free-surface equation at node m=5
F[ 66] = -1.02999206386122921e-17   1.02999206386122921e-17   1.029992064e-06  dynamic Bernoulli free-surface equation at node m=6
F[ 67] = -1.62494800601264977e-17   1.62494800601264977e-17   1.624948006e-06  dynamic Bernoulli free-surface equation at node m=7
F[ 68] = -2.29850860566926940e-17   2.29850860566926940e-17   2.298508606e-06  dynamic Bernoulli free-surface equation at node m=8
F[ 69] = -3.49384150083453804e-17   3.49384150083453804e-17   3.493841501e-06  dynamic Bernoulli free-surface equation at node m=9
F[ 70] = -2.38117902132128911e-17   2.38117902132128911e-17   2.381179021e-06  dynamic Bernoulli free-surface equation at node m=10
F[ 71] = -2.13858878522765750e-17   2.13858878522765750e-17   2.138588785e-06  dynamic Bernoulli free-surface equation at node m=11
F[ 72] = -3.59413020178944720e-17   3.59413020178944720e-17   3.594130202e-06  dynamic Bernoulli free-surface equation at node m=12
F[ 73] = -2.99510850149120600e-17   2.99510850149120600e-17   2.995108501e-06  dynamic Bernoulli free-surface equation at node m=13
F[ 74] = -3.51281503885303437e-17   3.51281503885303437e-17   3.512815039e-06  dynamic Bernoulli free-surface equation at node m=14
F[ 75] = -3.17671236538252799e-17   3.17671236538252799e-17   3.176712365e-06  dynamic Bernoulli free-surface equation at node m=15
F[ 76] = -3.76218153852470039e-17   3.76218153852470039e-17   3.762181539e-06  dynamic Bernoulli free-surface equation at node m=16
F[ 77] = -1.94072188874905294e-17   1.94072188874905294e-17   1.940721889e-06  dynamic Bernoulli free-surface equation at node m=17
F[ 78] = -3.49113099540332428e-17   3.49113099540332428e-17   3.491130995e-06  dynamic Bernoulli free-surface equation at node m=18
F[ 79] = -3.97631146759058751e-17   3.97631146759058751e-17   3.976311468e-06  dynamic Bernoulli free-surface equation at node m=19
F[ 80] = -2.11148373091551989e-17   2.11148373091551989e-17   2.111483731e-06  dynamic Bernoulli free-surface equation at node m=20
F[ 81] = -2.20093041014557400e-17   2.20093041014557400e-17   2.200930410e-06  dynamic Bernoulli free-surface equation at node m=21
F[ 82] = -3.07371315899640507e-17   3.07371315899640507e-17   3.073713159e-06  dynamic Bernoulli free-surface equation at node m=22
F[ 83] = -2.82434665932473905e-17   2.82434665932473905e-17   2.824346659e-06  dynamic Bernoulli free-surface equation at node m=23
F[ 84] = -2.68068987147040971e-17   2.68068987147040971e-17   2.680689871e-06  dynamic Bernoulli free-surface equation at node m=24
F[ 85] = -1.46909394371785851e-17   1.46909394371785851e-17   1.469093944e-06  dynamic Bernoulli free-surface equation at node m=25
F[ 86] = -4.33409818451080397e-17   4.33409818451080397e-17   4.334098185e-06  dynamic Bernoulli free-surface equation at node m=26
F[ 87] = -3.59413020178944720e-17   3.59413020178944720e-17   3.594130202e-06  dynamic Bernoulli free-surface equation at node m=27
F[ 88] = -3.45589442479754538e-17   3.45589442479754538e-17   3.455894425e-06  dynamic Bernoulli free-surface equation at node m=28
F[ 89] = -1.76182853028894471e-17   1.76182853028894471e-17   1.761828530e-06  dynamic Bernoulli free-surface equation at node m=29
F[ 90] = -2.53974358904729414e-17   2.53974358904729414e-17   2.539743589e-06  dynamic Bernoulli free-surface equation at node m=30
F[ 91] = -2.37711326317446847e-17   2.37711326317446847e-17   2.377113263e-06  dynamic Bernoulli free-surface equation at node m=31
F[ 92] = -2.95174041459178582e-17   2.95174041459178582e-17   2.951740415e-06  dynamic Bernoulli free-surface equation at node m=32
F[ 93] = -2.86229373536173171e-17   2.86229373536173171e-17   2.862293735e-06  dynamic Bernoulli free-surface equation at node m=33
F[ 94] = -3.05473962097790874e-17   3.05473962097790874e-17   3.054739621e-06  dynamic Bernoulli free-surface equation at node m=34
F[ 95] = -1.08149166705429067e-17   1.08149166705429067e-17   1.081491667e-06  dynamic Bernoulli free-surface equation at node m=35
F[ 96] = -3.07371315899640507e-17   3.07371315899640507e-17   3.073713159e-06  dynamic Bernoulli free-surface equation at node m=36
F[ 97] = -1.98137947021725935e-17   1.98137947021725935e-17   1.981379470e-06  dynamic Bernoulli free-surface equation at node m=37
F[ 98] = -2.04914210599760338e-17   2.04914210599760338e-17   2.049142106e-06  dynamic Bernoulli free-surface equation at node m=38
F[ 99] = -1.28206906896410899e-17   1.28206906896410899e-17   1.282069069e-06  dynamic Bernoulli free-surface equation at node m=39
F[100] = -1.75369701399530342e-17   1.75369701399530342e-17   1.753697014e-06  dynamic Bernoulli free-surface equation at node m=40
F[101] = -1.96511643762997679e-17   1.96511643762997679e-17   1.965116438e-06  dynamic Bernoulli free-surface equation at node m=41
F[102] = -3.51823604971546189e-17   3.51823604971546189e-17   3.518236050e-06  dynamic Bernoulli free-surface equation at node m=42
F[103] = -2.49908600757908772e-17   2.49908600757908772e-17   2.499086008e-06  dynamic Bernoulli free-surface equation at node m=43
F[104] = -3.06287113727155003e-17   3.06287113727155003e-17   3.062871137e-06  dynamic Bernoulli free-surface equation at node m=44
F[105] = -1.24954300378954386e-17   1.24954300378954386e-17   1.249543004e-06  dynamic Bernoulli free-surface equation at node m=45
F[106] = -1.91090632900570156e-17   1.91090632900570156e-17   1.910906329e-06  dynamic Bernoulli free-surface equation at node m=46
F[107] = -1.33085816672595669e-17   1.33085816672595669e-17   1.330858167e-06  dynamic Bernoulli free-surface equation at node m=47
F[108] = -1.65882932390282178e-17   1.65882932390282178e-17   1.658829324e-06  dynamic Bernoulli free-surface equation at node m=48
F[109] = -4.27988807588652875e-17   4.27988807588652875e-17   4.279888076e-06  dynamic Bernoulli free-surface equation at node m=49
F[110] = -2.36627124144961343e-17   2.36627124144961343e-17   2.366271241e-06  dynamic Bernoulli free-surface equation at node m=50
