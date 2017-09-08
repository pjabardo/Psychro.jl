# =================================================================================================
# =                                         Hyland83.jl                                           =
# =================================================================================================

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


# Specific entropy of saturated water vapor
# First table of the appendix of ref. [1]
@test 14.30387 ≈ Psychro.entropyvapor(173.15)/1000 atol=0.00005
@test 11.10965 ≈ Psychro.entropyvapor(223.15)/1000 atol=0.00005
@test  9.15510 ≈ Psychro.entropyvapor(273.16)/1000 atol=0.00005
@test  8.07477 ≈ Psychro.entropyvapor(323.15)/1000 atol=0.00005
@test  7.35365 ≈ Psychro.entropyvapor(373.15)/1000 atol=0.00005
@test  6.83731 ≈ Psychro.entropyvapor(423.15)/1000 atol=0.00005
@test  6.43218 ≈ Psychro.entropyvapor(473.15)/1000 atol=0.00005

