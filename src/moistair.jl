#=
This file actually implements the user interface of the Psychro library.
=#


"""
    volume(DryAir, T, P)

Specific volume of dry air. Temperature should be in the range `173.15 K < T < 473.15` and pressure should be below 5 MPa.

If `Unitful` units are not used, the output is in `m^3/kg`, temperaute in K and pressure in Pa.
"""
function volume(::Type{DryAir}, Tk, P)
    @assert 173.1 < Tk < 473.2 "Out of range error: Temperature should be between 173.15 K and 473.15 K!"
    @assert 0 <= P < 5e6 "Out of range error: Pressure should be below 5×10⁶ Pa"
    volumeair(Tk, P)
end


"""
    volumem(DryAir, T, P)

Molar volume of dry air. Temperature should be in the range `173.15 K < T < 473.15` and pressure should be below 5 MPa.

If `Unitful` units are not used, the output is in `m^3/mol`, temperaute in K and pressure in Pa.
"""
function volumem(::Type{DryAir}, Tk, P)
    @assert 173.1 < Tk < 473.2 "Out of range error: Temperature should be between 173.15 K and 473.15 K!"
    @assert 0 <= P < 5e6 "Out of range error: Pressure should be below 5×10⁶ Pa"
    molarvolumeair(Tk, P)
end
    
   
"""
    density(DryAir, T, P)

Density of dry air. Temperature should be in the range `173.15 K < T < 473.15` and pressure should be below 5 MPa.

If `Unitful` units are not used, the output is in `kg/m^3`, temperaute in K and pressure in Pa.
"""
function density(::Type{DryAir}, Tk, P)
    @assert 173.1 < Tk < 473.2 "Out of range error: Temperature should be between 173.15 K and 473.15 K!"
    @assert 0 <= P < 5e6 "Out of range error: Pressure should be below 5×10⁶ Pa"
    1.0/volumeair(Tk, P)
end

"""
    enthalpy(DryAir, T, P)

Specific enthalpy of dry air. Temperature should be in the range `173.15 K < T < 473.15` and pressure should be below 5 MPa.

If `Unitful` units are not used, the output is in `J/kg`, temperaute in K and pressure in Pa.
"""
function enthalpy(::Type{DryAir}, Tk, P)
    @assert 173.1 < Tk < 473.2 "Out of range error: Temperature should be between 173.15 K and 473.15 K!"
    @assert 0 <= P < 5e6 "Out of range error: Pressure should be below 5×10⁶ Pa"
    enthalpyair(Tk, P)
end

"""
    enthalpym(DryAir, T, P)

Molar enthalpy of dry air. Temperature should be in the range `173.15 K < T < 473.15` and pressure should be below 5 MPa.

If `Unitful` units are not used, the output is in `J/mol`, temperaute in K and pressure in Pa.
"""
function enthalpym(::Type{DryAir}, Tk, P)
    @assert 173.1 < Tk < 473.2 "Out of range error: Temperature should be between 173.15 K and 473.15 K!"
    @assert 0 <= P < 5e6 "Out of range error: Pressure should be below 5×10⁶ Pa"
    molarenthalpyair(Tk, P)
end

"""
    entropy(DryAir, T, P)

Specific entropy of dry air. Temperature should be in the range `173.15 K < T < 473.15` and pressure should be below 5 MPa.

If `Unitful` units are not used, the output is in `J/kg/K`, temperaute in K and pressure in Pa.
"""
function entropy(::Type{DryAir}, Tk, P)
    @assert 173.1 < Tk < 473.2 "Out of range error: Temperature should be between 173.15 K and 473.15 K!"
    @assert 0 <= P < 5e6 "Out of range error: Pressure should be below 5×10⁶ Pa"
    entropyair(Tk, P)
end

"""
    entropy(DryAir, T, P)

Molar entropy of dry air. Temperature should be in the range `173.15 K < T < 473.15` and pressure should be below 5 MPa.

If `Unitful` units are not used, the output is in `J/mol/K`, temperaute in K and pressure in Pa.
"""
function entropym(::Type{DryAir}, Tk, P)
    @assert 173.1 < Tk < 473.2 "Out of range error: Temperature should be between 173.15 K and 473.15 K!"
    @assert 0 <= P < 5e6 "Out of range error: Pressure should be below 5×10⁶ Pa"
    molarentropyair(Tk, P)
end

"""
    compressfactor(DryAir, T, P)

Compressibility factor Z of dry air. Temperature should be in the range `173.15 K < T < 473.15` and pressure should be below 5 MPa.

If `Unitful` units are not used, the output is dimensionless, temperaute in K and pressure in Pa.
"""
function compressfactor(::Type{DryAir}, Tk, P)
    @assert 173.1 < Tk < 473.2 "Out of range error: Temperature should be between 173.15 K and 473.15 K!"
    @assert 0 <= P < 5e6 "Out of range error: Pressure should be below 5×10⁶ Pa"
    P*molarvolumeair(Tk, P) / (R*Tk)
end



"""
    volume(Vapor, T)

Specific volume of saturated water vapor. Temperature should be in the range `173.15 K < T < 473.15`.

If `Unitful` units are not used, the output is in `m^3/kg`, temperaute in K.
"""
function volume(::Type{Vapor}, Tk)
    @assert 173.1 < Tk < 473.2 "Out of range error: Temperature should be between 173.15 K and 473.15 K!"
    volumevapor(Tk)
end
    
"""
    volumem(Vapor, T)

Molar volume of saturated water vapor. Temperature should be in the range `173.15 K < T < 473.15`.

If `Unitful` units are not used, the output is in `m^3/mol`, temperaute in K.
"""
function volumem(::Type{Vapor}, Tk)
    @assert 173.1 < Tk < 473.2 "Out of range error: Temperature should be between 173.15 K and 473.15 K!"
    molarvolvapor(Tk)
end

"""
    density(Vapor, T)

Density of saturated water vapor. Temperature should be in the range `173.15 K < T < 473.15`.

If `Unitful` units are not used, the output is in `kg/m^3`, temperaute in K.
"""
function density(::Type{Vapor}, Tk)
    @assert 173.1 < Tk < 473.2 "Out of range error: Temperature should be between 173.15 K and 473.15 K!"
    1/volumevapor(Tk)
end

"""
    enthalpy(Vapor, T)

Specific enthalpy of saturated water vapor. Temperature should be in the range `173.15 K < T < 473.15`.

If `Unitful` units are not used, the output is in `J/kg`, temperaute in K.
"""
function enthalpy(::Type{Vapor}, Tk)
    @assert 173.1 < Tk < 473.2 "Out of range error: Temperature should be between 173.15 K and 473.15 K!"
    enthalpyvapor(Tk)
end

"""
    enthalpym(Vapor, T)

Molar enthalpy of saturated water vapor. Temperature should be in the range `173.15 K < T < 473.15`.

If `Unitful` units are not used, the output is in `J/mol`, temperaute in K.
"""
function enthalpym(::Type{Vapor}, Tk)
    @assert 173.1 < Tk < 473.2 "Out of range error: Temperature should be between 173.15 K and 473.15 K!"
    enthalpyvapor(Tk)*Mv
end

"""
    entropy(Vapor, T)

Specific entropy of saturated water vapor. Temperature should be in the range `173.15 K < T < 473.15`.

If `Unitful` units are not used, the output is in `J/kg/K`, temperaute in K.
"""
function entropym(::Type{Vapor}, Tk)
    @assert 173.1 < Tk < 473.2 "Out of range error: Temperature should be between 173.15 K and 473.15 K!"
    molarentropyvapor(Tk)
end

"""
    entropym(Vapor, T)

Molar entropy of saturated water vapor. Temperature should be in the range `173.15 K < T < 473.15`.

If `Unitful` units are not used, the output is in `J/mol/K`, temperaute in K.
"""
function entropy(::Type{Vapor}, Tk)
    @assert 173.1 < Tk < 473.2 "Out of range error: Temperature should be between 173.15 K and 473.15 K!"
    entropyvapor(Tk)
end

"""
    compressfactor(Vapor, T)

Compressibility factor Z of saturated water vapor. Temperature should be in the range `173.15 K < T < 473.15`.

If `Unitful` units are not used, the output is in dimensionless, temperaute in K.
"""
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
@doc (@doc volume) volumem


function density(::Type{MoistAir}, Tk, ::Type{T}, y, P) where {T<:PsychroProperty}
    @assert 173.1 < Tk < 473.2 "Out of range error: Temperature should be between 173.15 K and 473.15 K!"
    @assert 0 <= P < 5e6 "Out of range error: Pressure should be below 5×10⁶ Pa"

    xv = molarfrac(Tk, T, y, P)
    M = xv*Mv + (1-xv)*Ma
    vm = molarvolumemoist(Tk, P, xv)
    return M/vm
end
@doc (@doc volume) density


function enthalpy(::Type{MoistAir}, Tk, ::Type{T}, y, P) where {T<:PsychroProperty}
    @assert 173.1 < Tk < 473.2 "Out of range error: Temperature should be between 173.15 K and 473.15 K!"
    @assert 0 <= P < 5e6 "Out of range error: Pressure should be below 5×10⁶ Pa"

    xv = molarfrac(Tk, T, y, P)
    enthalpymoist(Tk, P, xv)
end
@doc (@doc volume) enthalpy

function enthalpym(::Type{MoistAir}, Tk, ::Type{T}, y, P) where {T<:PsychroProperty}
    @assert 173.1 < Tk < 473.2 "Out of range error: Temperature should be between 173.15 K and 473.15 K!"
    @assert 0 <= P < 5e6 "Out of range error: Pressure should be below 5×10⁶ Pa"

    xv = molarfrac(Tk, T, y, P)
    molarenthalpymoist(Tk, P, xv)
end
@doc (@doc volume) enthalpym

function entropy(::Type{MoistAir}, Tk, ::Type{T}, y, P) where {T<:PsychroProperty}
    @assert 173.1 < Tk < 473.2 "Out of range error: Temperature should be between 173.15 K and 473.15 K!"
    @assert 0 <= P < 5e6 "Out of range error: Pressure should be below 5×10⁶ Pa"

    xv = molarfrac(Tk, T, y, P)
    entropymoist(Tk, P, xv)
end
@doc (@doc volume) entropy

function entropym(::Type{MoistAir}, Tk, ::Type{T}, y, P) where {T<:PsychroProperty}
    @assert 173.1 < Tk < 473.2 "Out of range error: Temperature should be between 173.15 K and 473.15 K!"
    @assert 0 <= P < 5e6 "Out of range error: Pressure should be below 5×10⁶ Pa"

    xv = molarfrac(Tk, T, y, P)
    molarentropymoist(Tk, P, xv)
end
@doc (@doc volume) entropym

function compressfactor(::Type{MoistAir}, Tk, ::Type{T}, y, P) where {T<:PsychroProperty}
    @assert 173.1 < Tk < 473.2 "Out of range error: Temperature should be between 173.15 K and 473.15 K!"
    @assert 0 <= P < 5e6 "Out of range error: Pressure should be below 5×10⁶ Pa"

    xv = molarfrac(Tk, T, y, P)
    Zmoist(Tk, P, xv)
end
@doc (@doc volume) compressfactor

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
@doc (@doc volume) dewpoint


function relhum(::Type{MoistAir}, Tk, ::Type{T}, y, P) where {T<:PsychroProperty}
    @assert 173.1 < Tk < 473.2 "Out of range error: Temperature should be between 173.15 K and 473.15 K!"
    @assert 0 <= P < 5e6 "Out of range error: Pressure should be below 5×10⁶ Pa"

    if T==RelHum
        return y
    end
    
    xv = molarfrac(Tk, T, y, P)
    xv * P / (efactor(Tk, P) * Pws(Tk))
end
@doc (@doc volume) relhum





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
@doc (@doc volume) wetbulb



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

function molarfrac(::Type{MoistAir}, Tk, ::Type{T}, y, P) where {T<:PsychroProperty}
    @assert 173.1 < Tk < 473.2 "Out of range error: Temperature should be between 173.15 K and 473.15 K!"
    @assert 0 <= P < 5e6 "Out of range error: Pressure should be below 5×10⁶ Pa"

    if T==MolarFrac
        return y
    end
    
    xv = molarfrac(Tk, T, y, P)
    return xv
end
@doc (@doc volume) molarfrac


function spechum(::Type{MoistAir}, Tk, ::Type{T}, y, P) where {T<:PsychroProperty}
    @assert 173.1 < Tk < 473.2 "Out of range error: Temperature should be between 173.15 K and 473.15 K!"
    @assert 0 <= P < 5e6 "Out of range error: Pressure should be below 5×10⁶ Pa"

    if T==SpecHum
        return y
    end
    xv = molarfrac(Tk, T, y, P)
    return xv*Mv / ( (1-xv)*Ma + xv*Mv )
    
end
@doc (@doc volume) spechum


