include("../src/Psychro.jl")
using Base.Test


# Testing the polynomial evaluation macro:
p = randn(6)
x = 0.33
@test @evalpoly(x, p[1], p[2]) ≈ Psychro.@polyeval(x, p, 2) atol=1e-13
@test @evalpoly(x, p[1], p[2], p[3]) ≈ Psychro.@polyeval(x, p, 3) atol=1e-13
@test @evalpoly(x, p[1], p[2], p[3], p[4]) ≈ Psychro.@polyeval(x, p, 4) atol=1e-13
@test @evalpoly(x, p[1], p[2], p[3], p[4], p[5]) ≈ Psychro.@polyeval(x, p, 5) atol=1e-13
@test @evalpoly(x, p[1], p[2], p[3], p[4], p[5], p[6]) ≈ Psychro.@polyeval(x, p, 6) atol=1e-13



#=
Testing the virial coefficients. Table 1 from reference [2].
=#
# Baa
@test -55.94 ≈ Psychro.Baa(173.15)*1e6 atol=0.01
@test -13.15 ≈ Psychro.Baa(273.15)*1e6 atol=0.01
@test 3.72 ≈ Psychro.Baa(373.15)*1e6 atol=0.01
@test 12.31 ≈ Psychro.Baa(473.15)*1e6 atol=0.01

# dBaa/dT
@test 0.7240 ≈ Psychro.dBaa(173.15)*1e6 atol=0.0001
@test 0.2460 ≈ Psychro.dBaa(273.15)*1e6 atol=0.0001
@test 0.1146 ≈ Psychro.dBaa(373.15)*1e6 atol=0.0001
@test 0.0640 ≈ Psychro.dBaa(473.15)*1e6 atol=0.0001

# Caaa
@test 2267.0 ≈ Psychro.Caaa(173.15)*1e12 atol=1.0
@test 1409.0 ≈ Psychro.Caaa(273.15)*1e12 atol=1.0
@test 1202.0 ≈ Psychro.Caaa(373.15)*1e12 atol=1.0
@test 1139.0 ≈ Psychro.Caaa(473.15)*1e12 atol=1.0

# dCaaa/dT
@test -18.00 ≈ Psychro.dCaaa(173.15)*1e12 atol=0.01
@test -3.65 ≈ Psychro.dCaaa(273.15)*1e12 atol=0.01
@test -1.06 ≈ Psychro.dCaaa(373.15)*1e12 atol=0.01
@test -0.34 ≈ Psychro.dCaaa(473.15)*1e12 atol=0.01

# Baw
@test -93.3 ≈ Psychro.Baw(173.15)*1e6 atol=0.1
@test -36.4 ≈ Psychro.Baw(273.15)*1e6 atol=0.1
@test -14.5 ≈ Psychro.Baw(373.15)*1e6 atol=0.1
@test -3.1 ≈ Psychro.Baw(473.15)*1e6 atol=0.1

# dBaw/dT
@test 1.01 ≈ Psychro.dBaw(173.15)*1e6 atol=0.01
@test 0.318 ≈ Psychro.dBaw(273.15)*1e6 atol=0.001
@test 0.151 ≈ Psychro.dBaw(373.15)*1e6 atol=0.001
@test 0.0869 ≈ Psychro.dBaw(473.15)*1e6 atol=0.0001

# Caaw
@test 1023 ≈ Psychro.Caaw(173.15)*1e12 atol=1.0
@test 861 ≈ Psychro.Caaw(273.15)*1e12 atol=1.0
@test 696 ≈ Psychro.Caaw(373.15)*1e12 atol=1.0
@test 627 ≈ Psychro.Caaw(473.15)*1e12 atol=1.0

# dCaaw/dT
@test  5.56 ≈ Psychro.dCaaw(173.15)*1e12 atol=0.01
@test -2.44 ≈ Psychro.dCaaw(273.15)*1e12 atol=0.01
@test -1.02 ≈ Psychro.dCaaw(373.15)*1e12 atol=0.01
@test -0.457 ≈ Psychro.dCaaw(473.15)*1e12 atol=0.001

# Caww
@test -2.0e7 ≈ Psychro.Caww(173.15)*1e12 atol=0.1e7
@test -2.2e5 ≈ Psychro.Caww(273.15)*1e12 atol=0.1e5
@test -3.0e4 ≈ Psychro.Caww(373.15)*1e12 atol=0.1e4
@test -8.4e3 ≈ Psychro.Caww(473.15)*1e12 atol=0.1e3

# dCaww/dT
@test 1.61e6 ≈ Psychro.dCaww(173.15)*1e12 atol=0.01e6
@test 6.05e3 ≈ Psychro.dCaww(273.15)*1e12 atol=0.01e3
@test    456 ≈ Psychro.dCaww(373.15)*1e12 atol=1.0
@test   86.9 ≈ Psychro.dCaww(473.15)*1e12 atol=0.1
# Note that on table 1, the last line is 8.69 instead of 86.9. Probably a typo in the paper.


# Specific volume of dry air
# Table 3
P = 0.1e6
@test 0.49510 ≈ Psychro.volumeair(173.15, P) atol=0.00003
@test 0.78364 ≈ Psychro.volumeair(273.15, P) atol=0.0001
@test 1.0713 ≈ Psychro.volumeair(373.15, P) atol=0.0001
@test 1.3586 ≈ Psychro.volumeair(473.15, P) atol=0.0001

P = 0.5e6
@test 0.09747 ≈ Psychro.volumeair(173.15, P) atol=0.00003
@test 0.15637 ≈ Psychro.volumeair(273.15, P) atol=0.00003
@test 0.21436 ≈ Psychro.volumeair(373.15, P) atol=0.00003
@test 0.27207 ≈ Psychro.volumeair(473.15, P) atol=0.00003

P = 1e6
@test 0.04775 ≈ Psychro.volumeair(173.15, P) atol=0.00003
@test 0.07797 ≈ Psychro.volumeair(273.15, P) atol=0.00003
@test 0.10725 ≈ Psychro.volumeair(373.15, P) atol=0.00003
@test 0.13625 ≈ Psychro.volumeair(473.15, P) atol=0.00003

P = 5e6
@test 0.00795 ≈ Psychro.volumeair(173.15, P) atol=0.00003
@test 0.01533 ≈ Psychro.volumeair(273.15, P) atol=0.00003
@test 0.02162 ≈ Psychro.volumeair(373.15, P) atol=0.00003
@test 0.02763 ≈ Psychro.volumeair(473.15, P) atol=0.00003

# Specific enthalpy of dry air
# Table 3
P = 0.1e6
@test -100.63 ≈ Psychro.enthalpyair(173.15, P)/1000 atol=0.01
@test    0.004≈ Psychro.enthalpyair(273.15, P)/1000 atol=0.01
@test  100.79 ≈ Psychro.enthalpyair(373.15, P)/1000 atol=0.01
@test  202.55 ≈ Psychro.enthalpyair(473.15, P)/1000 atol=0.01

P = 0.5e6
@test -103.17 ≈ Psychro.enthalpyair(173.15, P)/1000 atol=0.01
@test   -1.10 ≈ Psychro.enthalpyair(273.15, P)/1000 atol=0.01
@test  100.26 ≈ Psychro.enthalpyair(373.15, P)/1000 atol=0.01
@test  202.30 ≈ Psychro.enthalpyair(473.15, P)/1000 atol=0.01

P = 1e6
@test -106.42 ≈ Psychro.enthalpyair(173.15, P)/1000 atol=0.01
@test   -2.48 ≈ Psychro.enthalpyair(273.15, P)/1000 atol=0.01
@test   99.60 ≈ Psychro.enthalpyair(373.15, P)/1000 atol=0.01
@test  202.00 ≈ Psychro.enthalpyair(473.15, P)/1000 atol=0.01

P = 5e6
@test -135.55 ≈ Psychro.enthalpyair(173.15, P)/1000 atol=0.01
@test  -13.15 ≈ Psychro.enthalpyair(273.15, P)/1000 atol=0.01
@test   94.63 ≈ Psychro.enthalpyair(373.15, P)/1000 atol=0.01
@test  199.8 ≈ Psychro.enthalpyair(473.15, P)/1000 atol=0.1

# Specific entropy of dry air
P = 0.1e6
@test -0.4550 ≈ Psychro.entropyair(173.15, P)/1000 atol=0.0001
@test  0.0038 ≈ Psychro.entropyair(273.15, P)/1000 atol=0.0001
@test  0.3182 ≈ Psychro.entropyair(373.15, P)/1000 atol=0.0001
@test  0.5597 ≈ Psychro.entropyair(473.15, P)/1000 atol=0.0001


# Table 4, appendix of reference [2]
# Testing the volume of saturated moist air
volmoist(Tk, P) = Psychro.volumemoist(Tk, P, Psychro.efactor(Tk,P)*Psychro.Pws(Tk)/P)
P = 0.1e6
@test 0.49510 ≈ volmoist(173.15, P) atol=0.00005
@test 0.78846 ≈ volmoist(273.15, P) atol=0.00005
@test 3.4974 ≈ volmoist(363.15, P) atol=0.0002
P = 0.5e6
@test 0.097466 ≈ volmoist(173.15, P) atol=0.000006
@test 0.15657 ≈ volmoist(273.15, P) atol=0.00005
@test 0.24273 ≈ volmoist(363.15, P) atol=0.00003
@test 0.26904 ≈ volmoist(373.15, P) atol=0.00003
P=1e6
@test 0.047752 ≈ volmoist(173.15, P) atol=0.000006
@test 0.078022 ≈ volmoist(273.15, P) atol=0.000005
@test 0.11226 ≈ volmoist(363.15, P) atol=0.00001
@test 0.11942 ≈ volmoist(373.15, P) atol=0.00001
@test 0.23281 ≈ volmoist(423.15, P) atol=0.00002
@test 0.32711 ≈ volmoist(433.15, P) atol=0.00002
@test 0.61489 ≈ volmoist(443.15, P) atol=0.00005

P=5e6
@test 0.007951 ≈ volmoist(173.15, P) atol=0.000006
@test 0.015331 ≈ volmoist(273.15, P) atol=0.000005
@test 0.021310 ≈ volmoist(363.15, P) atol=0.00001
@test 0.022073 ≈ volmoist(373.15, P) atol=0.00001
@test 0.02727  ≈ volmoist(423.15, P) atol=0.00001
@test 0.02885  ≈ volmoist(433.15, P) atol=0.00001
@test 0.03077  ≈ volmoist(443.15, P) atol=0.00001
@test 0.03315  ≈ volmoist(453.15, P) atol=0.00001
@test 0.0402  ≈ volmoist(473.15, P) atol=0.0001

hmoist(Tk, P) = Psychro.enthalpymoist(Tk, P, Psychro.efactor(Tk,P)*Psychro.Pws(Tk)/P)/1000
P = 0.1e6
@test -100.627 ≈ hmoist(173.15, P) atol=0.005
@test    9.602 ≈ hmoist(273.15, P) atol=0.005
@test   4035.0 ≈ hmoist(363.15, P) atol=0.3 # In the paper the number is 4305!

P = 0.5e6
@test -103.168 ≈ hmoist(173.15, P) atol=0.005
@test    0.831 ≈ hmoist(273.15, P) atol=0.005
@test   365.55 ≈ hmoist(363.15, P) atol=0.03
@test   533.27 ≈ hmoist(373.15, P) atol=0.03
@test   35957 ≈ hmoist(423.15, P) atol=2.0

P = 1e6
@test -106.415 ≈ hmoist(173.15, P) atol=0.005
@test   -1.498 ≈ hmoist(273.15, P) atol=0.005
@test   217.77 ≈ hmoist(363.15, P) atol=0.03
@test   293.17 ≈ hmoist(373.15, P) atol=0.03
@test   1788.8 ≈ hmoist(423.15, P) atol=0.5
@test   3113.3 ≈ hmoist(433.15, P) atol=0.5
@test   7206.0 ≈ hmoist(443.15, P) atol=2.0

P = 5e6
@test -135.547 ≈ hmoist(173.15, P) atol=0.005
@test  -12.929 ≈ hmoist(273.15, P) atol=0.005
@test   109.89 ≈ hmoist(363.15, P) atol=0.03
@test   132.42 ≈ hmoist(373.15, P) atol=0.03
@test   347.3 ≈ hmoist(423.15, P) atol=0.5
@test   428.8 ≈ hmoist(433.15, P) atol=0.5
@test   534.1 ≈ hmoist(443.15, P) atol=0.5
@test   672.0 ≈ hmoist(453.15, P) atol=1.0
@test   1113.0 ≈ hmoist(473.15, P) atol=1.0

# Entropy of saturated moist air
# Table 4 in the appendix of reference [2].
smoist(Tk, P) = Psychro.entropymoist(Tk, P, Psychro.efactor(Tk,P)*Psychro.Pws(Tk)/P)/1000
P = 0.1e6
@test  -0.45499 ≈ smoist(173.15, P) atol=0.00003
@test   0.04070 ≈ smoist(273.15, P) atol=0.00003
@test  11.7287 ≈ smoist(363.15, P) atol=0.0005

P = 0.5e6
@test  -0.92718 ≈ smoist(173.15, P) atol=0.00005
@test  -0.45417 ≈ smoist(273.15, P) atol=0.00005
@test   0.64517 ≈ smoist(363.15, P) atol=0.0005
@test   1.1102 ≈ smoist(373.15, P) atol=0.0005
@test  90.01 ≈ smoist(423.15, P) atol=0.05

P = 1e6
@test  -1.13929 ≈ smoist(173.15, P) atol=0.00005 # In the table the signal is +
@test  -0.66104 ≈ smoist(273.15, P) atol=0.00005
@test   0.00805 ≈ smoist(363.15, P) atol=0.00005
@test   0.21676 ≈ smoist(373.15, P) atol=0.0005
@test   4.0478  ≈ smoist(423.15, P) atol=0.0005
@test   7.3103  ≈ smoist(433.15, P) atol=0.0005
@test   17.221  ≈ smoist(443.15, P) atol=0.005

P = 5e6
@test  -1.72380 ≈ smoist(173.15, P) atol=0.0001 # In the table the signal is +
@test  -1.15924 ≈ smoist(273.15, P) atol=0.00005
@test  -0.77519 ≈ smoist(363.15, P) atol=0.00005
@test  -0.71331 ≈ smoist(373.15, P) atol=0.00005
@test  -0.1629  ≈ smoist(423.15, P) atol=0.0002
@test   0.0364  ≈ smoist(433.15, P) atol=0.0002
@test   0.2901  ≈ smoist(443.15, P) atol=0.0002
@test   0.6185  ≈ smoist(453.15, P) atol=0.0002
@test   1.642  ≈ smoist(473.15, P) atol=0.002


# Appendix of ref. [1]. First table
# Specific volume of saturated ice:
@test 1.0768e-3 ≈ Psychro.volumeice(173.15) atol=0.0001e-3
@test 1.0829e-3 ≈ Psychro.volumeice(223.15) atol=0.0001e-3
@test 1.0909e-3 ≈ Psychro.volumeice(273.16) atol=0.0001e-3

# Specific volume of saturated liquid water
@test 1.00021e-3 ≈ Psychro.volumewater(273.16) atol=0.00001e-3
@test 1.01215e-3 ≈ Psychro.volumewater(323.15) atol=0.00001e-3
@test 1.04346e-3 ≈ Psychro.volumewater(373.15) atol=0.00001e-3 # In the paper there is a typo in the power of 10: 1.04346e10
@test 1.09050e-3 ≈ Psychro.volumewater(423.15) atol=0.00001e-3
@test 1.15653e-3 ≈ Psychro.volumewater(473.15) atol=0.00001e-3


# Saturation enthalpy of saturated ice
@test -507.215e3 ≈ Psychro.enthalpyice(173.15) atol=0.001e3
@test -429.413e3 ≈ Psychro.enthalpyice(223.15) atol=0.001e3
@test -333.429e3 ≈ Psychro.enthalpyice(273.15) atol=0.001e3  # The table presents -333.409. Typo???

# Saturation enthalpy of saturated liquid water
@test 0.0 ≈ Psychro.enthalpywater(273.16) atol=1.0
@test 209.330e3 ≈ Psychro.enthalpywater(323.15) atol=5
@test 419.158e3 ≈ Psychro.enthalpywater(373.15) atol=5
@test 632.210e3 ≈ Psychro.enthalpywater(423.15) atol=9
@test 852.329e3 ≈ Psychro.enthalpywater(473.15) atol=15


# Table 2 of reference [2]
# Enhancement factor
P = 0.1e6
@test 1.0105 ≈ Psychro.efactor(173.15, P) atol=0.0001
@test 1.0039 ≈ Psychro.efactor(273.15, P) atol=0.0001
@test 1.0039 ≈ Psychro.efactor(363.15, P) atol=0.0001
P = 0.5e6
@test 1.054 ≈ Psychro.efactor(173.15, P) atol=0.001
@test 1.0177 ≈ Psychro.efactor(273.15, P) atol=0.0006  # I think the table has a typo!!?
@test 1.0180 ≈ Psychro.efactor(363.15, P) atol=0.0001
@test 1.0188 ≈ Psychro.efactor(373.15, P) atol=0.0001
@test 1.0022 ≈ Psychro.efactor(423.15, P) atol=0.0001
P = 1.0e6
@test 1.113 ≈ Psychro.efactor(173.15, P) atol=0.001
@test 1.0353 ≈ Psychro.efactor(273.15, P) atol=0.0015  # I think the table has a typo!!?
@test 1.0284 ≈ Psychro.efactor(363.15, P) atol=0.0001
@test 1.0295 ≈ Psychro.efactor(373.15, P) atol=0.0001
@test 1.0288 ≈ Psychro.efactor(423.15, P) atol=0.0001
@test 1.0235 ≈ Psychro.efactor(433.15, P) atol=0.0001
@test 1.0142 ≈ Psychro.efactor(443.15, P) atol=0.0001


P = 5.0e6
@test 1.820 ≈ Psychro.efactor(173.15, P) atol=0.001
@test 1.191 ≈ Psychro.efactor(273.15, P) atol=0.008  # I think the table has a typo!!?
@test 1.1102 ≈ Psychro.efactor(363.15, P) atol=0.0001
@test 1.1082 ≈ Psychro.efactor(373.15, P) atol=0.0001
@test 1.111 ≈ Psychro.efactor(423.15, P) atol=0.001
@test 1.114 ≈ Psychro.efactor(433.15, P) atol=0.001
@test 1.116 ≈ Psychro.efactor(443.15, P) atol=0.001
@test 1.117 ≈ Psychro.efactor(453.15, P) atol=0.001
@test 1.116 ≈ Psychro.efactor(473.15, P) atol=0.001

# Saturation pressure of water vapor over ice
# Second table in the appendix of reference [1]
@test 1.40510e-3 ≈ Psychro.Pws_s(173.15) atol=0.00001e-3
@test 6.11153e2 ≈ Psychro.Pws_s(273.15) atol=0.00001e2

# Saturation pressure of water vapor over liquid water
# Second table in the appendix of reference [1]
@test 6.11213e2 ≈ Psychro.Pws_l(273.15) atol=0.00001e2
@test 1.01419e5 ≈ Psychro.Pws_l(373.15) atol=0.00001e5
@test 1.55507e6 ≈ Psychro.Pws_l(473.15) atol=0.00001e6



# Testing saturation temperature function:
@test Psychro.Tws(Psychro.Pws_s(173.15)) ≈ 173.15 atol=1e-5
@test Psychro.Tws(Psychro.Pws_s(223.15)) ≈ 223.15 atol=1e-5
@test Psychro.Tws(Psychro.Pws_s(273.15)) ≈ 273.15 atol=1e-5
@test Psychro.Tws(Psychro.Pws_s(273.16)) ≈ 273.16 atol=1e-5
@test Psychro.Tws(Psychro.Pws_l(273.16)) ≈ 273.16 atol=1e-5
@test Psychro.Tws(Psychro.Pws_l(323.15)) ≈ 323.15 atol=1e-5
@test Psychro.Tws(Psychro.Pws_l(373.15)) ≈ 373.15 atol=1e-5
@test Psychro.Tws(Psychro.Pws_l(473.15)) ≈ 473.15 atol=1e-5


# First virial coefficient of saturated vapor B'
# Second table in the appendix of reference [1]
@test -3.2939e-5 ≈ Psychro.Blin(173.15) atol=0.0001e-5
@test -8.3497e-7 ≈ Psychro.Blin(273.15) atol=0.0001e-7
@test -1.4658e-7 ≈ Psychro.Blin(373.15) atol=0.0001e-5
@test -5.0508e-8 ≈ Psychro.Blin(473.15) atol=0.0001e-5

# Second virial coefficient of saturated vapor C'
# Second table in the appendix of reference [1]
@test -4.6563e-9 ≈ Psychro.Clin(173.15) atol=0.0001e-9
@test -2.0928e-12 ≈ Psychro.Clin(273.15) atol=0.0001e-12
@test -5.7548e-14 ≈ Psychro.Clin(373.15) atol=0.0001e-14
@test -6.3933e-15 ≈ Psychro.Clin(473.15) atol=0.0001e-15

# Specific enthalpy of saturated water vapor
# First table of the appendix of ref. [1]
@test 2315.87 ≈ Psychro.enthalpyvapor(173.15)/1000 atol=0.01
@test 2408.41 ≈ Psychro.enthalpyvapor(223.15)/1000 atol=0.01
@test 2500.81 ≈ Psychro.enthalpyvapor(273.16)/1000 atol=0.01
@test 2591.29 ≈ Psychro.enthalpyvapor(323.15)/1000 atol=0.01
@test 2675.46 ≈ Psychro.enthalpyvapor(373.15)/1000 atol=0.01
@test 2746.15 ≈ Psychro.enthalpyvapor(423.15)/1000 atol=0.01
@test 2793.11 ≈ Psychro.enthalpyvapor(473.15)/1000 atol=0.01

# Specific volume of saturated water vapor
# First table of the appendix of ref. [1]
@test 5.6873e7 ≈ Psychro.volumevapor(173.15) atol=0.0001e7
@test 2.6146e4 ≈ Psychro.volumevapor(223.15) atol=0.0001e4
@test 2.0601e2 ≈ Psychro.volumevapor(273.16) atol=0.0001e2
@test 1.2030e1 ≈ Psychro.volumevapor(323.15) atol=0.0001e1
@test 1.6718e0 ≈ Psychro.volumevapor(373.15) atol=0.0001e0
@test 3.9253e-1 ≈ Psychro.volumevapor(423.15) atol=0.0001e-1
@test 1.2722e-1 ≈ Psychro.volumevapor(473.15) atol=0.0001e0


