
using Unitful

# Testing th units part.
@testset "Units" begin
Tf = 70.0u"°F"
Tc = uconvert(u"°C", Tf)
Tk = uconvert(u"K", Tf)
Pa = 1.0u"atm"
Pl = uconvert(u"lbf/inch^2", Pa)
Pk = uconvert(u"kPa", Pa)
Pp = uconvert(u"Pa", Pa)

v1 = volume(DryAir, Tf, Pa, u"inch^3/lb")
v2 = volume(DryAir, Tk.val, Pp.val)
@test uconvert(u"m^3/kg", v1).val ≈ v2 rtol=1e-8

v1 = volumem(DryAir, Tf, Pa, u"inch^3/kmol")
v2 = volumem(DryAir, Tk.val, Pp.val)
@test uconvert(u"m^3/mol", v1).val ≈ v2 rtol=1e-8

v1 = density(DryAir, Tf, Pa, u"lb/inch^3")
v2 = density(DryAir, Tk.val, Pp.val)
@test uconvert(u"kg/m^3", v1).val ≈ v2 rtol=1e-8


v1 = enthalpy(DryAir, Tf, Pa, u"inch^2/hr^2")
v2 = enthalpy(DryAir, Tk.val, Pp.val)
@test uconvert(u"J/kg", v1).val ≈ v2 rtol=1e-8

v1 = enthalpym(DryAir, Tf, Pa, u"lb*inch^2/hr^2/kmol")
v2 = enthalpym(DryAir, Tk.val, Pp.val)
@test uconvert(u"J/mol", v1).val ≈ v2 rtol=1e-8

v1 = entropy(DryAir, Tf, Pa, u"inch^2/hr^2/Ra")
v2 = entropy(DryAir, Tk.val, Pp.val)
@test uconvert(u"J/kg/K", v1).val ≈ v2 rtol=1e-8

v1 = entropym(DryAir, Tf, Pa, u"lb*inch^2/hr^2/Ra/kmol")
v2 = entropym(DryAir, Tk.val, Pp.val)
@test uconvert(u"J/mol/K", v1).val ≈ v2 rtol=1e-8

v1 = compressfactor(DryAir, Tf, Pa)
v2 = compressfactor(DryAir, Tk.val, Pp.val)
@test v1 ≈ v2 rtol=1e-8


v1 = volume(Vapor, Tf, u"inch^3/lb")
v2 = volume(Vapor, Tk.val)
@test uconvert(u"m^3/kg", v1).val ≈ v2 rtol=1e-8


v1 = volumem(Vapor, Tf, u"inch^3/kmol")
v2 = volumem(Vapor, Tk.val)
@test uconvert(u"m^3/mol", v1).val ≈ v2 rtol=1e-8

v1 = density(Vapor, Tf, u"lb/inch^3")
v2 = density(Vapor, Tk.val)
@test uconvert(u"kg/m^3", v1).val ≈ v2 rtol=1e-8


v1 = enthalpy(Vapor, Tf, u"inch^2/hr^2")
v2 = enthalpy(Vapor, Tk.val)
@test uconvert(u"J/kg", v1).val ≈ v2 rtol=1e-8

v1 = enthalpym(Vapor, Tf, u"lb*inch^2/hr^2/kmol")
v2 = enthalpym(Vapor, Tk.val)
@test uconvert(u"J/mol", v1).val ≈ v2 rtol=1e-8

v1 = entropy(Vapor, Tf, u"inch^2/hr^2/Ra")
v2 = entropy(Vapor, Tk.val)
@test uconvert(u"J/kg/K", v1).val ≈ v2 rtol=1e-8

v1 = entropym(Vapor, Tf, u"lb*inch^2/hr^2/Ra/kmol")
v2 = entropym(Vapor
, Tk.val)
@test uconvert(u"J/mol/K", v1).val ≈ v2 rtol=1e-8

v1 = compressfactor(Vapor, Tf)
v2 = compressfactor(Vapor, Tk.val)
@test v1 ≈ v2 rtol=1e-8

Df = 60.0u"°F"
Dc = uconvert(u"°C", Df)
Dk = uconvert(u"K", Df)

r = Psychro.relhum(MoistAir, Tk.val, DewPoint, Dk.val, Pp.val)
B = Psychro.wetbulb(MoistAir, Tk.val, DewPoint, Dk.val, Pp.val)
w = Psychro.humrat(MoistAir, Tk.val, DewPoint, Dk.val, Pp.val)
q = Psychro.spechum(MoistAir, Tk.val, DewPoint, Dk.val, Pp.val)
x = Psychro.molarfrac(Tk.val, DewPoint, Dk.val, Pp.val)

Bf = uconvert(u"Ra", B*u"K")

x1 = Psychro.molarfrac(Tf, DewPoint, Df, Pl)
x2 = Psychro.molarfrac(Tf, WetBulb, Bf, Pl)
x3 = Psychro.molarfrac(Tf, HumRat, w, Pl)
x4 = Psychro.molarfrac(Tf, SpecHum, q, Pl)
x5 = Psychro.molarfrac(Tf, RelHum, r, Pl)
@test x1 ≈ x rtol=1e-8
@test x2 ≈ x rtol=1e-8
@test x3 ≈ x rtol=1e-8
@test x4 ≈ x rtol=1e-8
@test x5 ≈ x rtol=1e-8

x  = volume(MoistAir, Tk.val, DewPoint, Dk.val, Pp.val)
x1 = volume(MoistAir, Tf, DewPoint, Df, Pl, u"inch^3/lb")
x2 = volume(MoistAir, Tf, WetBulb, Bf, Pl, u"inch^3/lb")
x3 = volume(MoistAir, Tf, HumRat, w, Pl, u"inch^3/lb")
x4 = volume(MoistAir, Tf, SpecHum, q, Pl, u"inch^3/lb")
x5 = volume(MoistAir, Tf, RelHum, r, Pl, u"inch^3/lb")
@test uconvert(u"m^3/kg", x1).val ≈ x rtol=1e-8
@test uconvert(u"m^3/kg", x2).val ≈ x rtol=1e-8
@test uconvert(u"m^3/kg", x3).val ≈ x rtol=1e-8
@test uconvert(u"m^3/kg", x4).val ≈ x rtol=1e-8
@test uconvert(u"m^3/kg", x5).val ≈ x rtol=1e-8

x  = volumem(MoistAir, Tk.val, DewPoint, Dk.val, Pp.val)
x1 = volumem(MoistAir, Tf, DewPoint, Df, Pl, u"inch^3/kmol")
x2 = volumem(MoistAir, Tf, WetBulb, Bf, Pl, u"inch^3/kmol")
x3 = volumem(MoistAir, Tf, HumRat, w, Pl, u"inch^3/kmol")
x4 = volumem(MoistAir, Tf, SpecHum, q, Pl, u"inch^3/kmol")
x5 = volumem(MoistAir, Tf, RelHum, r, Pl, u"inch^3/kmol")
@test uconvert(u"m^3/mol", x1).val ≈ x rtol=1e-8
@test uconvert(u"m^3/mol", x2).val ≈ x rtol=1e-8
@test uconvert(u"m^3/mol", x3).val ≈ x rtol=1e-8
@test uconvert(u"m^3/mol", x4).val ≈ x rtol=1e-8
@test uconvert(u"m^3/mol", x5).val ≈ x rtol=1e-8


x  = density(MoistAir, Tk.val, DewPoint, Dk.val, Pp.val)
x1 = density(MoistAir, Tf, DewPoint, Df, Pl, u"lb/inch^3")
x2 = density(MoistAir, Tf, WetBulb, Bf, Pl, u"lb/inch^3")
x3 = density(MoistAir, Tf, HumRat, w, Pl, u"lb/inch^3")
x4 = density(MoistAir, Tf, SpecHum, q, Pl, u"lb/inch^3")
x5 = density(MoistAir, Tf, RelHum, r, Pl, u"lb/inch^3")
@test uconvert(u"kg/m^3", x1).val ≈ x rtol=1e-8
@test uconvert(u"kg/m^3", x2).val ≈ x rtol=1e-8
@test uconvert(u"kg/m^3", x3).val ≈ x rtol=1e-8
@test uconvert(u"kg/m^3", x4).val ≈ x rtol=1e-8
@test uconvert(u"kg/m^3", x5).val ≈ x rtol=1e-8


x  = enthalpy(MoistAir, Tk.val, DewPoint, Dk.val, Pp.val)
x1 = enthalpy(MoistAir, Tf, DewPoint, Df, Pl, u"inch^2/hr^2")
x2 = enthalpy(MoistAir, Tf, WetBulb, Bf, Pl, u"inch^2/hr^2")
x3 = enthalpy(MoistAir, Tf, HumRat, w, Pl, u"inch^2/hr^2")
x4 = enthalpy(MoistAir, Tf, SpecHum, q, Pl, u"inch^2/hr^2")
x5 = enthalpy(MoistAir, Tf, RelHum, r, Pl, u"inch^2/hr^2")
@test uconvert(u"J/kg", x1).val ≈ x rtol=1e-8
@test uconvert(u"J/kg", x2).val ≈ x rtol=1e-8
@test uconvert(u"J/kg", x3).val ≈ x rtol=1e-8
@test uconvert(u"J/kg", x4).val ≈ x rtol=1e-8
@test uconvert(u"J/kg", x5).val ≈ x rtol=1e-8


x  = enthalpym(MoistAir, Tk.val, DewPoint, Dk.val, Pp.val)
x1 = enthalpym(MoistAir, Tf, DewPoint, Df, Pl, u"inch^2*lb/hr^2/kmol")
x2 = enthalpym(MoistAir, Tf, WetBulb, Bf, Pl, u"inch^2*lb/hr^2/kmol")
x3 = enthalpym(MoistAir, Tf, HumRat, w, Pl, u"inch^2*lb/hr^2/kmol")
x4 = enthalpym(MoistAir, Tf, SpecHum, q, Pl, u"inch^2*lb/hr^2/kmol")
x5 = enthalpym(MoistAir, Tf, RelHum, r, Pl, u"inch^2*lb/hr^2/kmol")
@test uconvert(u"J/mol", x1).val ≈ x rtol=1e-8
@test uconvert(u"J/mol", x2).val ≈ x rtol=1e-8
@test uconvert(u"J/mol", x3).val ≈ x rtol=1e-8
@test uconvert(u"J/mol", x4).val ≈ x rtol=1e-8
@test uconvert(u"J/mol", x5).val ≈ x rtol=1e-8


x  = entropy(MoistAir, Tk.val, DewPoint, Dk.val, Pp.val)
x1 = entropy(MoistAir, Tf, DewPoint, Df, Pl, u"inch^2/hr^2/Ra")
x2 = entropy(MoistAir, Tf, WetBulb, Bf, Pl, u"inch^2/hr^2/Ra")
x3 = entropy(MoistAir, Tf, HumRat, w, Pl, u"inch^2/hr^2/Ra")
x4 = entropy(MoistAir, Tf, SpecHum, q, Pl, u"inch^2/hr^2/Ra")
x5 = entropy(MoistAir, Tf, RelHum, r, Pl, u"inch^2/hr^2/Ra")
@test uconvert(u"J/kg/K", x1).val ≈ x rtol=1e-8
@test uconvert(u"J/kg/K", x2).val ≈ x rtol=1e-8
@test uconvert(u"J/kg/K", x3).val ≈ x rtol=1e-8
@test uconvert(u"J/kg/K", x4).val ≈ x rtol=1e-8
@test uconvert(u"J/kg/K", x5).val ≈ x rtol=1e-8


x  = entropym(MoistAir, Tk.val, DewPoint, Dk.val, Pp.val)
x1 = entropym(MoistAir, Tf, DewPoint, Df, Pl, u"inch^2/hr^2/Ra*lb/kmol")
x2 = entropym(MoistAir, Tf, WetBulb, Bf, Pl, u"inch^2/hr^2/Ra*lb/kmol")
x3 = entropym(MoistAir, Tf, HumRat, w, Pl, u"inch^2/hr^2/Ra*lb/kmol")
x4 = entropym(MoistAir, Tf, SpecHum, q, Pl, u"inch^2/hr^2/Ra*lb/kmol")
x5 = entropym(MoistAir, Tf, RelHum, r, Pl, u"inch^2/hr^2/Ra*lb/kmol")
@test uconvert(u"J/mol/K", x1).val ≈ x rtol=1e-8
@test uconvert(u"J/mol/K", x2).val ≈ x rtol=1e-8
@test uconvert(u"J/mol/K", x3).val ≈ x rtol=1e-8
@test uconvert(u"J/mol/K", x4).val ≈ x rtol=1e-8
@test uconvert(u"J/mol/K", x5).val ≈ x rtol=1e-8

x  = compressfactor(MoistAir, Tk.val, DewPoint, Dk.val, Pp.val)
x1 = compressfactor(MoistAir, Tf, DewPoint, Df, Pl)
x2 = compressfactor(MoistAir, Tf, WetBulb, Bf, Pl)
x3 = compressfactor(MoistAir, Tf, HumRat, w, Pl)
x4 = compressfactor(MoistAir, Tf, SpecHum, q, Pl)
x5 = compressfactor(MoistAir, Tf, RelHum, r, Pl)
@test x1 ≈ x rtol=1e-8
@test x2 ≈ x rtol=1e-8
@test x3 ≈ x rtol=1e-8
@test x4 ≈ x rtol=1e-8
@test x5 ≈ x rtol=1e-8

#xv = Psychro.molarfrac(Tk.val, DewPoint, Dk.val, Pp.val)
#xv1 = Psychro.molarfrac(Tk, DewPoint, Df, Pl)


x  = dewpoint(MoistAir, Tk.val, DewPoint, Dk.val, Pp.val)
x1 = dewpoint(MoistAir, Tf, DewPoint, Df, Pl, u"Ra")
x2 = dewpoint(MoistAir, Tf, WetBulb, Bf, Pl, u"Ra")
x3 = dewpoint(MoistAir, Tf, HumRat, w, Pl, u"Ra")
x4 = dewpoint(MoistAir, Tf, SpecHum, q, Pl, u"Ra")
x5 = dewpoint(MoistAir, Tf, RelHum, r, Pl, u"Ra")
@test uconvert(u"K", x1).val ≈ x rtol=1e-8
@test uconvert(u"K", x2).val ≈ x rtol=1e-8
@test uconvert(u"K", x3).val ≈ x rtol=1e-8
@test uconvert(u"K", x4).val ≈ x rtol=1e-8
@test uconvert(u"K", x5).val ≈ x rtol=1e-8


x  = wetbulb(MoistAir, Tk.val, DewPoint, Dk.val, Pp.val)
x1 = wetbulb(MoistAir, Tf, DewPoint, Df, Pl, u"Ra")
x2 = wetbulb(MoistAir, Tf, WetBulb, Bf, Pl, u"Ra")
x3 = wetbulb(MoistAir, Tf, HumRat, w, Pl, u"Ra")
x4 = wetbulb(MoistAir, Tf, SpecHum, q, Pl, u"Ra")
x5 = wetbulb(MoistAir, Tf, RelHum, r, Pl, u"Ra")
@test uconvert(u"K", x1).val ≈ x rtol=1e-8
@test uconvert(u"K", x2).val ≈ x rtol=1e-8
@test uconvert(u"K", x3).val ≈ x rtol=1e-8
@test uconvert(u"K", x4).val ≈ x rtol=1e-8
@test uconvert(u"K", x5).val ≈ x rtol=1e-8


x  = humrat(MoistAir, Tk.val, DewPoint, Dk.val, Pp.val)
x1 = humrat(MoistAir, Tf, DewPoint, Df, Pl)
x2 = humrat(MoistAir, Tf, WetBulb, Bf, Pl)
x3 = humrat(MoistAir, Tf, HumRat, w, Pl)
x4 = humrat(MoistAir, Tf, SpecHum, q, Pl)
x5 = humrat(MoistAir, Tf, RelHum, r, Pl)
@test x1 ≈ x rtol=1e-8
@test x2 ≈ x rtol=1e-8
@test x3 ≈ x rtol=1e-8
@test x4 ≈ x rtol=1e-8
@test x5 ≈ x rtol=1e-8


x  = relhum(MoistAir, Tk.val, DewPoint, Dk.val, Pp.val)
x1 = relhum(MoistAir, Tf, DewPoint, Df, Pl)
x2 = relhum(MoistAir, Tf, WetBulb, Bf, Pl)
x3 = relhum(MoistAir, Tf, HumRat, w, Pl)
x4 = relhum(MoistAir, Tf, SpecHum, q, Pl)
x5 = relhum(MoistAir, Tf, RelHum, r, Pl)
@test x1 ≈ x rtol=1e-8
@test x2 ≈ x rtol=1e-8
@test x3 ≈ x rtol=1e-8
@test x4 ≈ x rtol=1e-8
@test x5 ≈ x rtol=1e-8


x  = spechum(MoistAir, Tk.val, DewPoint, Dk.val, Pp.val)
x1 = spechum(MoistAir, Tf, DewPoint, Df, Pl)
x2 = spechum(MoistAir, Tf, WetBulb, Bf, Pl)
x3 = spechum(MoistAir, Tf, HumRat, w, Pl)
x4 = spechum(MoistAir, Tf, SpecHum, q, Pl)
x5 = spechum(MoistAir, Tf, RelHum, r, Pl)
@test x1 ≈ x rtol=1e-8
@test x2 ≈ x rtol=1e-8
@test x3 ≈ x rtol=1e-8
@test x4 ≈ x rtol=1e-8
@test x5 ≈ x rtol=1e-8
end