include("../src/Psychro.jl")
using Base.Test

# write your own tests here

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
