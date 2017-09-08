#=
This file actually implements the user interface of the Psychro library.
=#


"""
    volume(DryAir, T, P[, out_unit=u"m^3/kg"])
    volume(Vapor, T[, out_unit=u"m^3/kg"])
    volume(MoistAir, T, HumidityType, hum, P[, out_unit="m^3/kg"])

Calculates the specific volume of a gas. In the case of moist air, parameters specifying the humidity shoud be added and the the property is based on the mass of dry air not mass of the gas (as is the case for `DryAir` or `Vapor`). 

If `Unitful` units are used, all parameters of the function (except dimensionless) should have units associated. When no units are provided, all parameters should use SI units (and so does the return type). For further information, checkout the doc for `Psychro` module.

Temperature should be in the range `173.15 K < T < 473.15` and pressure should be below 5 MPa.

# Examples
```julia-repl
julia> volume(DryAir, 293.15, 101325.0)
0.8302353995092402

julia> volume(DryAir, 20.0u"°C", 1.0u"atm")
0.8302353995092402 kg^-1 m^3

julia> volume(Vapor, 293.15)
57.77492045936827

julia> volume(Vapor, 20.0u"°C")
57.77492045936827 kg^-1 m^3

julia> volume(MoistAir, 293.15, RelHum, 0.7, 101325.0)
0.843889817602806

julia> volume(MoistAir, 20.0u"°C", WetBulb, 18.0u"°C", 1.0u"atm")
0.8464079202783964 kg^-1 m^3

julia> volume(MoistAir, 20.0u"°C", WetBulb, 18.0u"°C", 1.0u"atm", u"inch^3/lb")
23428.490579267214 in^3 lb^-1
```
"""
function volume(::Type{DryAir}, Tk, P)
    @assert 173.1 < Tk < 473.2 "Out of range error: Temperature should be between 173.15 K and 473.15 K!"
    @assert 0 <= P < 5e6 "Out of range error: Pressure should be below 5×10⁶ Pa"
    volumeair(Tk, P)
end


"""
    volumem(DryAir, T, P[, out_unit=u"m^3/kg"])
    volumem(Vapor, T[, out_unit=u"m^3/kg"])
    volumem(MoistAir, T, HumidityType, hum, P[, out_unit="m^3/kg"])

Calculates the molar volume of a gas. In the case of moist air, parameters specifying the humidity shoud be added. 

If `Unitful` units are used, all parameters of the function (except dimensionless) should have units associated. When no units are provided, all parameters should use SI units (and so does the return type). For further information, checkout the doc for `Psychro` module.

Temperature should be in the range `173.15 K < T < 473.15` and pressure should be below 5 MPa.

# Examples
```julia-repl
julia> volumem(DryAir, 293.15, 101325.0)
0.024046522993685877

julia> volumem(DryAir, 20.0u"°C", 1.0u"atm")
0.024046522993685877 m^3 mol^-1

julia> volumem(Vapor, 293.15)
1.040831369053248

julia> volumem(Vapor, 20.0u"°C")
1.040831369053248 m^3 mol^-1

julia> volumem(MoistAir, 293.15, RelHum, 0.7, 101325.0)
0.02404547233948706

julia> volumem(MoistAir, 20.0u"°C", WetBulb, 18.0u"°C", 1.0u"atm")
0.024045209859305458 m^3 mol^-1

julia> volumem(MoistAir, 20.0u"°C", WetBulb, 18.0u"°C", 1.0u"atm", u"inch^3/mol")
1467.32873315839 in^3 mol^-1
```
"""
function volumem(::Type{DryAir}, Tk, P)
    @assert 173.1 < Tk < 473.2 "Out of range error: Temperature should be between 173.15 K and 473.15 K!"
    @assert 0 <= P < 5e6 "Out of range error: Pressure should be below 5×10⁶ Pa"
    molarvolumeair(Tk, P)
end
    
   
"""
    density(DryAir, T, P[, out_unit=u"m^3/kg"])
    density(Vapor, T[, out_unit=u"m^3/kg"])
    density(MoistAir, T, HumidityType, hum, P[, out_unit="m^3/kg"])

Calculates the density volume of a gas. In the case of moist air, parameters specifying the humidity shoud be added. 

If `Unitful` units are used, all parameters of the function (except dimensionless) should have units associated. When no units are provided, all parameters should use SI units (and so does the return type). For further information, checkout the doc for `Psychro` module.

Temperature should be in the range `173.15 K < T < 473.15` and pressure should be below 5 MPa.

# Examples
```julia-repl
julia> density(DryAir, 293.15, 101325.0)
1.2044776705391136

julia> density(DryAir, 20.0u"°C", 1.0u"atm")
1.2044776705391136 kg m^-3

julia> density(Vapor, 293.15)
0.017308548277505224

julia> density(Vapor, 20.0u"°C")
0.017308548277505224 kg m^-3

julia> density(MoistAir, 293.15, RelHum, 0.7, 101325.0)
1.1971436091840653

julia> density(MoistAir, 20.0u"°C", WetBulb, 18.0u"°C", 1.0u"atm")
1.195819185774941 kg m^-3

julia> density(MoistAir, 20.0u"°C", WetBulb, 18.0u"°C", 1.0u"atm", u"lb/inch^3")
4.3201708903793606e-5 in^-3 lb

```
"""
function density(::Type{DryAir}, Tk, P)
    @assert 173.1 < Tk < 473.2 "Out of range error: Temperature should be between 173.15 K and 473.15 K!"
    @assert 0 <= P < 5e6 "Out of range error: Pressure should be below 5×10⁶ Pa"
    1.0/volumeair(Tk, P)
end

"""
    enthalpy(DryAir, T, P[, out_unit=u"m^3/kg"])
    enthalpy(Vapor, T[, out_unit=u"m^3/kg"])
    enthalpy(MoistAir, T, HumidityType, hum, P[, out_unit="m^3/kg"])

Calculates the specific enthalpy of a gas. In the case of moist air, parameters specifying the humidity shoud be added and the the property is based on the mass of dry air not mass of the gas (as is the case for `DryAir` or `Vapor`). 

If `Unitful` units are used, all parameters of the function (except dimensionless) should have units associated. When no units are provided, all parameters should use SI units (and so does the return type). For further information, checkout the doc for `Psychro` module.

Temperature should be in the range `173.15 K < T < 473.15` and pressure should be below 5 MPa.

# Examples
```julia-repl
julia> enthalpy(DryAir, 293.15, 101325.0)
20121.2722813185

julia> enthalpy(DryAir, 20.0u"°C", 1.0u"atm")
20121.2722813185 kg^-1 J

julia> enthalpy(Vapor, 293.15)
2.5374003352412493e6

julia> enthalpy(Vapor, 20.0u"°C")
2.5374003352412493e6 kg^-1 J

julia> enthalpy(MoistAir, 293.15, RelHum, 0.7, 101325.0)
46142.773129687484

julia> enthalpy(MoistAir, 20.0u"°C", WetBulb, 18.0u"°C", 1.0u"atm")
50944.84575377501 kg^-1 J

julia> enthalpy(MoistAir, 20.0u"°C", WetBulb, 18.0u"°C", 1.0u"atm", u"J/lb")
23108.193324739244 J lb^-1

```
"""
function enthalpy(::Type{DryAir}, Tk, P)
    @assert 173.1 < Tk < 473.2 "Out of range error: Temperature should be between 173.15 K and 473.15 K!"
    @assert 0 <= P < 5e6 "Out of range error: Pressure should be below 5×10⁶ Pa"
    enthalpyair(Tk, P)
end

"""
    enthalpym(DryAir, T, P[, out_unit=u"m^3/kg"])
    enthalpym(Vapor, T[, out_unit=u"m^3/kg"])
    enthalpym(MoistAir, T, HumidityType, hum, P[, out_unit="m^3/kg"])

Calculates the molar enthalpy of a gas. In the case of moist air, parameters specifying the humidity shoud be added. 

If `Unitful` units are used, all parameters of the function (except dimensionless) should have units associated. When no units are provided, all parameters should use SI units (and so does the return type). For further information, checkout the doc for `Psychro` module.

Temperature should be in the range `173.15 K < T < 473.15` and pressure should be below 5 MPa.

# Examples
```julia-repl
julia> enthalpym(DryAir, 293.15, 101325.0)
582.7824697199684

julia> enthalpym(DryAir, 20.0u"°C", 1.0u"atm")
582.7824697199684 J mol^-1

julia> enthalpy(Vapor, 293.15)
2.5374003352412493e6

julia> enthalpym(Vapor, 20.0u"°C")
45711.977511464975 J mol^-1

julia> enthalpym(MoistAir, 293.15, RelHum, 0.7, 101325.0)
1314.7744549269437

julia> enthalpym(MoistAir, 20.0u"°C", WetBulb, 18.0u"°C", 1.0u"atm")
1447.2684837312868 J mol^-1

julia> enthalpym(MoistAir, 20.0u"°C", WetBulb, 18.0u"°C", 1.0u"atm", u"kJ/kmol")
1447.2684837312868 kJ kmol^-1

```
"""
function enthalpym(::Type{DryAir}, Tk, P)
    @assert 173.1 < Tk < 473.2 "Out of range error: Temperature should be between 173.15 K and 473.15 K!"
    @assert 0 <= P < 5e6 "Out of range error: Pressure should be below 5×10⁶ Pa"
    molarenthalpyair(Tk, P)
end

"""
    entropy(DryAir, T, P[, out_unit=u"m^3/kg"])
    entropy(Vapor, T[, out_unit=u"m^3/kg"])
    entropy(MoistAir, T, HumidityType, hum, P[, out_unit="m^3/kg"])

Calculates the specific entropy of a gas. In the case of moist air, parameters specifying the humidity shoud be added and the the property is based on the mass of dry air not mass of the gas (as is the case for `DryAir` or `Vapor`). 

If `Unitful` units are used, all parameters of the function (except dimensionless) should have units associated. When no units are provided, all parameters should use SI units (and so does the return type). For further information, checkout the doc for `Psychro` module.

Temperature should be in the range `173.15 K < T < 473.15` and pressure should be below 5 MPa.

# Examples
```julia-repl
julia> entropy(DryAir, 293.15, 101325.0)
71.09262577195514

julia> entropy(DryAir, 20.0u"°C", 1.0u"atm")
71.09262577195514 kg^-1 J K^-1

julia> entropy(Vapor, 293.15)
8665.873003613997

julia> entropy(Vapor, 20.0u"°C")
8665.873003613997 kg^-1 J K^-1

julia> entropy(MoistAir, 293.15, RelHum, 0.7, 101325.0)
166.33538729355692

julia> entropy(MoistAir, 20.0u"°C", WetBulb, 18.0u"°C", 1.0u"atm")
182.97140931650105 kg^-1 J K^-1

julia> entropy(MoistAir, 20.0u"°C", WetBulb, 18.0u"°C", 1.0u"atm", u"kJ/lb/°F")
0.04610801955228433 °F^-1 kJ lb^-1

```
"""
function entropy(::Type{DryAir}, Tk, P)
    @assert 173.1 < Tk < 473.2 "Out of range error: Temperature should be between 173.15 K and 473.15 K!"
    @assert 0 <= P < 5e6 "Out of range error: Pressure should be below 5×10⁶ Pa"
    entropyair(Tk, P)
end

"""
    entropym(DryAir, T, P[, out_unit=u"m^3/kg"])
    entropym(Vapor, T[, out_unit=u"m^3/kg"])
    entropym(MoistAir, T, HumidityType, hum, P[, out_unit="m^3/kg"])

Calculates the molar entropy of a gas. In the case of moist air, parameters specifying the humidity shoud be added. 

If `Unitful` units are used, all parameters of the function (except dimensionless) should have units associated. When no units are provided, all parameters should use SI units (and so does the return type). For further information, checkout the doc for `Psychro` module.

Temperature should be in the range `173.15 K < T < 473.15` and pressure should be below 5 MPa.

# Examples
```julia-repl
julia> entropym(DryAir, 293.15, 101325.0)
2.0590912665460226

julia> entropym(DryAir, 20.0u"°C", 1.0u"atm")
2.0590912665460226 J K^-1 mol^-1

julia> entropy(Vapor, 293.15)
8665.873003613997

julia> entropym(Vapor, 20.0u"°C")
156.11812860454717 J K^-1 mol^-1

julia> entropym(MoistAir, 293.15, RelHum, 0.7, 101325.0)
4.739496639035868

julia> entropym(MoistAir, 20.0u"°C", WetBulb, 18.0u"°C", 1.0u"atm")
5.197949865380578 J K^-1 mol^-1

julia> entropym(MoistAir, 20.0u"°C", WetBulb, 18.0u"°C", 1.0u"atm", u"kJ/kmol/°F")
2.8877499252114327 °F^-1 kJ kmol^-1

```
"""
function entropym(::Type{DryAir}, Tk, P)
    @assert 173.1 < Tk < 473.2 "Out of range error: Temperature should be between 173.15 K and 473.15 K!"
    @assert 0 <= P < 5e6 "Out of range error: Pressure should be below 5×10⁶ Pa"
    molarentropyair(Tk, P)
end

"""
    compressfactor(DryAir, T, P[, out_unit=u"m^3/kg"])
    compressfactor(Vapor, T[, out_unit=u"m^3/kg"])
    compressfactor(MoistAir, T, HumidityType, hum, P[, out_unit="m^3/kg"])

Calculates the compressibility factor Z of a gas. In the case of moist air, parameters specifying the humidity shoud be added. 

If `Unitful` units are used, all parameters of the function (except dimensionless) should have units associated. When no units are provided, all parameters should use SI units (and so does the return type). For further information, checkout the doc for `Psychro` module.

Temperature should be in the range `173.15 K < T < 473.15` and pressure should be below 5 MPa.

# Examples
```julia-repl
julia> compressfactor(DryAir, 293.15, 101325.0)
1.2044776705391136

julia> compressfactor(DryAir, 20.0u"°C", 1.0u"atm")
1.2044776705391136 kg m^-3

julia> compressfactor(Vapor, 293.15)
0.017308548277505224

julia> compressfactor(Vapor, 20.0u"°C")
0.017308548277505224 kg m^-3

julia> compressfactor(MoistAir, 293.15, RelHum, 0.7, 101325.0)
1.1971436091840653

julia> compressfactor(MoistAir, 20.0u"°C", WetBulb, 18.0u"°C", 1.0u"atm")
1.195819185774941 kg m^-3

```
"""
function compressfactor(::Type{DryAir}, Tk, P)
    @assert 173.1 < Tk < 473.2 "Out of range error: Temperature should be between 173.15 K and 473.15 K!"
    @assert 0 <= P < 5e6 "Out of range error: Pressure should be below 5×10⁶ Pa"
    P*molarvolumeair(Tk, P) / (R*Tk)
end



function volume(::Type{Vapor}, Tk)
    @assert 173.1 < Tk < 473.2 "Out of range error: Temperature should be between 173.15 K and 473.15 K!"
    volumevapor(Tk)
end
    
function volumem(::Type{Vapor}, Tk)
    @assert 173.1 < Tk < 473.2 "Out of range error: Temperature should be between 173.15 K and 473.15 K!"
    molarvolvapor(Tk)
end

function density(::Type{Vapor}, Tk)
    @assert 173.1 < Tk < 473.2 "Out of range error: Temperature should be between 173.15 K and 473.15 K!"
    1/volumevapor(Tk)
end

function enthalpy(::Type{Vapor}, Tk)
    @assert 173.1 < Tk < 473.2 "Out of range error: Temperature should be between 173.15 K and 473.15 K!"
    enthalpyvapor(Tk)
end

function enthalpym(::Type{Vapor}, Tk)
    @assert 173.1 < Tk < 473.2 "Out of range error: Temperature should be between 173.15 K and 473.15 K!"
    enthalpyvapor(Tk)*Mv
end

function entropym(::Type{Vapor}, Tk)
    @assert 173.1 < Tk < 473.2 "Out of range error: Temperature should be between 173.15 K and 473.15 K!"
    molarentropyvapor(Tk)
end

function entropy(::Type{Vapor}, Tk)
    @assert 173.1 < Tk < 473.2 "Out of range error: Temperature should be between 173.15 K and 473.15 K!"
    entropyvapor(Tk)
end

function compressfactor(::Type{Vapor}, Tk)
    @assert 173.1 < Tk < 473.2 "Out of range error: Temperature should be between 173.15 K and 473.15 K!"
    Zvapor(Tk)
end



"""
    molarfrac(Tk, HumidityType, x, P)

Calculates the molarfraction of water vapor given a humidity parameter.
The humidity parameter can be one of the following:

 * `MolarFrac` - Identity function basically...
 * `RelHum` - Relative humidity in decimal notation (not percentage!)
 * `DewPoint`- Dew point temperature in K
 * `HumRat` - Humidity ration, kg of vapor / kg of dry air
 * `SpecHum` - Specific Humidity of water vapor = mv / (mv + ma)
 * `WetBulb` - Wet bulb temperature in K- estimated by the adiabatic saturation temperature.
 
As usual all calculations are performed using plain SI units. 
The parameters used by this function are:

 * `Tk` - Temperature in K
 * Humidity type - a type listed above that characterizes the humidity
 * The value of the humidity parameter. 
 * P - Pressure in Pa
"""
function molarfrac(Tk, ::Type{MolarFrac}, xv, P)
    xv
end

function molarfrac(Tk, ::Type{HumRat}, w, P)
    w  / (Mv/Ma + w)
end

function molarfrac(Tk, ::Type{SpecHum}, q, P)
    q * Ma / (Mv + r*(Ma - Mv))
end

function molarfrac(Tk, ::Type{RelHum}, rel, P)
    rel * efactor(Tk, P) * Pws(Tk) / P
end

function molarfrac(Tk, ::Type{DewPoint}, D, P)
    efactor(D,P) * Pws(D) / P
end

function molarfrac(Tk, ::Type{WetBulb}, B, P)
    w = calc_W_from_B(Tk, B, P)
    w / (Mv/Ma + w)
end



function volume(::Type{MoistAir}, Tk, ::Type{T}, y, P) where {T<:PsychroProperty}
    @assert 173.1 < Tk < 473.2 "Out of range error: Temperature should be between 173.15 K and 473.15! K"
    @assert 0 <= P < 5e6 "Out of range error: Pressure should be below 5×10⁶ Pa"
    xv = molarfrac(Tk, T, y, P)
    volumemoist(Tk, P, xv)
end


function volumem(::Type{MoistAir}, Tk, ::Type{T}, y, P) where {T<:PsychroProperty}
    @assert 173.1 < Tk < 473.2 "Out of range error: Temperature should be between 173.15 K and 473.15! K"
    @assert 0 <= P < 5e6 "Out of range error: Pressure should be below 5×10⁶ Pa"

    xv = molarfrac(Tk, T, y, P)
    molarvolumemoist(Tk, P, xv)
end


function density(::Type{MoistAir}, Tk, ::Type{T}, y, P) where {T<:PsychroProperty}
    @assert 173.1 < Tk < 473.2 "Out of range error: Temperature should be between 173.15 K and 473.15 K!"
    @assert 0 <= P < 5e6 "Out of range error: Pressure should be below 5×10⁶ Pa"

    xv = molarfrac(Tk, T, y, P)
    M = xv*Mv + (1-xv)*Ma
    vm = molarvolumemoist(Tk, P, xv)
    return M/vm
end


function enthalpy(::Type{MoistAir}, Tk, ::Type{T}, y, P) where {T<:PsychroProperty}
    @assert 173.1 < Tk < 473.2 "Out of range error: Temperature should be between 173.15 K and 473.15 K!"
    @assert 0 <= P < 5e6 "Out of range error: Pressure should be below 5×10⁶ Pa"

    xv = molarfrac(Tk, T, y, P)
    enthalpymoist(Tk, P, xv)
end

function enthalpym(::Type{MoistAir}, Tk, ::Type{T}, y, P) where {T<:PsychroProperty}
    @assert 173.1 < Tk < 473.2 "Out of range error: Temperature should be between 173.15 K and 473.15 K!"
    @assert 0 <= P < 5e6 "Out of range error: Pressure should be below 5×10⁶ Pa"

    xv = molarfrac(Tk, T, y, P)
    molarenthalpymoist(Tk, P, xv)
end

function entropy(::Type{MoistAir}, Tk, ::Type{T}, y, P) where {T<:PsychroProperty}
    @assert 173.1 < Tk < 473.2 "Out of range error: Temperature should be between 173.15 K and 473.15 K!"
    @assert 0 <= P < 5e6 "Out of range error: Pressure should be below 5×10⁶ Pa"

    xv = molarfrac(Tk, T, y, P)
    entropymoist(Tk, P, xv)
end

function entropym(::Type{MoistAir}, Tk, ::Type{T}, y, P) where {T<:PsychroProperty}
    @assert 173.1 < Tk < 473.2 "Out of range error: Temperature should be between 173.15 K and 473.15 K!"
    @assert 0 <= P < 5e6 "Out of range error: Pressure should be below 5×10⁶ Pa"

    xv = molarfrac(Tk, T, y, P)
    molarentropymoist(Tk, P, xv)
end

function compressfactor(::Type{MoistAir}, Tk, ::Type{T}, y, P) where {T<:PsychroProperty}
    @assert 173.1 < Tk < 473.2 "Out of range error: Temperature should be between 173.15 K and 473.15 K!"
    @assert 0 <= P < 5e6 "Out of range error: Pressure should be below 5×10⁶ Pa"

    xv = molarfrac(Tk, T, y, P)
    Zmoist(Tk, P, xv)
end

function calcdewpoint(Tk, P, xv, EPS=1e-9, MAXITER=100)
    # Use Ideal Gas to 
    D = Tws(xv*P)
    Dnew = D
    i = 0
    err = 0.0
    for i = 1:MAXITER
        f = efactor(D, P)
        Dnew = Tws(xv*P/f)
        err = abs(D-Dnew)
        if err < EPS
            return Dnew
        end
        D = Dnew
    end
    throw(ConvergenceError("Dew point calculation failed to converge!", D, i, err))
    return Dnew
end

"""
# Thermodynamic properties of moist air, dry air and saturated water vapor.

The methods listed below calculate thermodynamic properties of moist air:

    volume(MoistAir, T, HumidityType, humidity, P[, outunit]) 
    volume(MoistAir, T, HumidityType, humidity, P[, outunit]) 
    density(MoistAir, T, HumidityType, humidity, P[, outunit])
    enthalpy(MoistAir, T, HumidityType, humidity, P[, outunit])
    enthalpym(MoistAir, T, HumidityType, humidity, P[, outunit])
    entropy(MoistAir, T, HumidityType, humidity, P[, outunit])
    entropym(MoistAir, T, HumidityType, humidity, P[, outunit])
    compressfactor(MoistAir, T, HumidityType, humidity, P[, outunit])
    dewpoint(MoistAir, T, HumidityType, humidity, P[, outunit]) 
    wetbulb(MoistAir, T, HumidityType, humidity, P[, outunit]) 
    humrat(MoistAir, T, HumidityType, humidity, P) 
    relhum(MoistAir, T, HumidityType, humidity, P) 
    humrat(MoistAir, T, HumidityType, humidity, P) 
    spechum(MoistAir, T, HumidityType, humidity, P) 
    molarfrac(MoistAir, T, HumidityType, humidity, P) 
    
The methods listed above calculate the following thermodynamic properties of moist air:

 * `volume` Specific volume 
 * `volumem` Molar volume
 * `density` Density
 * `enthalpy` Specific enthalpy
 * `enthalpym` Molar enthalpy
 * `entropy` Specific entropy
 * `entropym` Molar entropy
 * `compressfactor` Compressibility factor Z 
 * `dewpoint` Dew point temperature
 * `wetbulb` Adiabatic saturation temperature
 * `humrat` Humidity ratio
 * `relhum` Relative humidity
 * `spechum` Specific humidity
 * `molarfrac` Molar fraction of water vapor

The humidity is specified using two parameters:

 * How the humidity is specified
 * The actual value of humidity

The following types are used to characterize the humidity.

 * `WetBulb` for wet bulb temperature, actually adiabatic saturation temperature
 * `DewPoint` Dew point temperature
 * `RelHum` Relative humidity
 * `HumRat` Humidity ratio (kg of vapor / kg of dry air)
 * `SpecHum` Specific humidity (kg of vapor / kg of moist air)
 * `MolarFrac` Molar fraction of water vapor.

# Examples
```julia-repl
julia> volume(MoistAir, 293.15, WetBulb, 291.15, 101325.0)
0.8464079202783964

julia> volume(MoistAir, 293.15, DewPoint, 291.15, 101325.0)
0.8475219875187474

julia> volume(MoistAir, 293.15, RelHum, 0.7, 101325.0)
0.843889817602806

julia> volume(MoistAir, 20.0u"°C", DewPoint, 60.0u"°F", 1.0u"atm")
0.8449934929585231 kg^-1 m^3

julia> volumem(MoistAir, 293.15, RelHum, 0.5, 93000.0)
0.026199080086890276

julia> volumem(MoistAir, 20.0u"°C", WetBulb, 17.0u"°C", 93u"kPa", u"inch^3/kmol")
1.598733210336603e6 in^3 kmol^-1

julia> density(MoistAir, 20.0u"°C", WetBulb, 17.0u"°C", 93u"kPa")
1.0976075893895811 kg m^-3

julia> density(MoistAir, 20.0u"°C", WetBulb, 17.0u"°C", 93u"kPa", u"lb/inch^3")
3.965358988338535e-5 in^-3 lb

julia> volumem(MoistAir, 20.0u"°C", WetBulb, 17.0u"°C", 93u"kPa", u"inch^3/kmol")
1.598733210336603e6 in^3 kmol^-1

julia> enthalpy(MoistAir, 20.0u"°C", WetBulb, 17.0u"°C", 93u"kPa")
50667.43014746832 kg^-1 J

julia> enthalpym(MoistAir, 20.0u"°C", WetBulb, 17.0u"°C", 93u"kPa")
1439.6551689935861 J mol^-1

julia> compressfactor(MoistAir, -90.0u"°C", RelHum, 0.01, 4.5u"MPa")
0.8552758629097985

julia> wetbulb(MoistAir, 20.0u"°C", WetBulb, 17.0u"°C", 93u"kPa", u"°C")
17.0 °C

julia> dewpoint(MoistAir, 20.0u"°C", WetBulb, 17.0u"°C", 93u"kPa", u"°C")
15.475836053510477 °C

julia> humrat(MoistAir, 20.0u"°C", WetBulb, 17.0u"°C", 93u"kPa")
0.012032930694441925

julia> relhum(MoistAir, 20.0u"°C", WetBulb, 17.0u"°C", 93u"kPa")
0.7517801524436909

julia> spechum(MoistAir, 20.0u"°C", WetBulb, 17.0u"°C", 93u"kPa")
0.011889860823189923


```

"""
function dewpoint(::Type{MoistAir}, Tk, ::Type{T}, y, P) where {T<:PsychroProperty}
    @assert 173.1 < Tk < 473.2 "Out of range error: Temperature should be between 173.15 K and 473.15 K!"
    @assert 0 <= P < 5e6 "Out of range error: Pressure should be below 5×10⁶ Pa"

    if T==DewPoint
        return y
    end
    
    xv = molarfrac(Tk, T, y, P)
    D = calcdewpoint(Tk, P, xv)
    return D
end


function relhum(::Type{MoistAir}, Tk, ::Type{T}, y, P) where {T<:PsychroProperty}
    @assert 173.1 < Tk < 473.2 "Out of range error: Temperature should be between 173.15 K and 473.15 K!"
    @assert 0 <= P < 5e6 "Out of range error: Pressure should be below 5×10⁶ Pa"

    if T==RelHum
        return y
    end
    
    xv = molarfrac(Tk, T, y, P)
    xv * P / (efactor(Tk, P) * Pws(Tk))
end
@doc (@doc dewpoint) relhum





function wetbulb(::Type{MoistAir}, Tk, ::Type{T}, y, P) where {T<:PsychroProperty}
    @assert 173.1 < Tk < 473.2 "Out of range error: Temperature should be between 173.15 K and 473.15 K!"
    @assert 0 <= P < 5e6 "Out of range error: Pressure should be below 5×10⁶ Pa"

    if T==WetBulb
        return y
    end
    xv = molarfrac(Tk, T, y, P)

    B = calcwetbulb(Tk, P, xv)
    return B
    
end
@doc (@doc dewpoint) wetbulb



humrat(xv) = xv / (1-xv) * (Mv/Ma)
function humrat(::Type{MoistAir}, Tk, ::Type{T}, y, P) where {T<:PsychroProperty}
    @assert 173.1 < Tk < 473.2 "Out of range error: Temperature should be between 173.15 K and 473.15 K!"
    @assert 0 <= P < 5e6 "Out of range error: Pressure should be below 5×10⁶ Pa"

    if T==HumRat
        return y
    end
    
    xv = molarfrac(Tk, T, y, P)
    return humrat(xv)
end
@doc (@doc dewpoint) humrat

function molarfrac(::Type{MoistAir}, Tk, ::Type{T}, y, P) where {T<:PsychroProperty}
    @assert 173.1 < Tk < 473.2 "Out of range error: Temperature should be between 173.15 K and 473.15 K!"
    @assert 0 <= P < 5e6 "Out of range error: Pressure should be below 5×10⁶ Pa"

    if T==MolarFrac
        return y
    end
    
    xv = molarfrac(Tk, T, y, P)
    return xv
end
@doc (@doc dewpoint) molarfrac


function spechum(::Type{MoistAir}, Tk, ::Type{T}, y, P) where {T<:PsychroProperty}
    @assert 173.1 < Tk < 473.2 "Out of range error: Temperature should be between 173.15 K and 473.15 K!"
    @assert 0 <= P < 5e6 "Out of range error: Pressure should be below 5×10⁶ Pa"

    if T==SpecHum
        return y
    end
    xv = molarfrac(Tk, T, y, P)
    return xv*Mv / ( (1-xv)*Ma + xv*Mv )
    
end
@doc (@doc dewpoint) spechum


