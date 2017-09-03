#=
This file actually implements the user interface of the Psychro library.
=#



volume(::Type{DryAir}, Tk, P) = volumeair(Tk, P)
molarvolume(::Type{DryAir}, Tk, P) = molarvolumeair(Tk, P)
density(::Type{DryAir}, Tk, P) = 1.0/volumeair(Tk, P)
enthalpy(::Type{DryAir}, Tk, P) = enthalpyair(Tk, P)
molarenthalpy(::Type{DryAir}, Tk, P) = molarenthalpyair(Tk, P)
entropy(::Type{DryAir}, Tk, P) = entropyair(Tk, P)
molarentropy(::Type{DryAir}, Tk, P) = molarentropyair(Tk, P)
compressfactor(::Type{DryAir}, Tk, P) = P*molarvolumeair(Tk, P) / (R*Tk)


volume(::Type{Vapor}, Tk) = volumevapor(Tk)
molarvolume(::Type{Vapor}, Tk) = molarvolvapor(Tk)
density(::Type{Vapor}, Tk) = 1/volumevapor(Tk)
enthalpy(::Type{Vapor}, Tk) = enthalpyvapor(Tk)
molarenthalpy(::Type{Vapor}, Tk) = enthalpyvapor(Tk)*Mv
molarentropy(::Type{Vapor}, Tk) = molarentropyvapor(Tk)
entropy(::Type{Vapor}, Tk) = entropyvapor(Tk)
compressfactor(::Type{Vapor}, Tk) = Zvapor(Tk)


"""
    molarfrac(Tk, HumidityType, x, P)

Calculates the molarfraction of water vapor given a humidity parameter.
The humidity parameter can be one of the following:

 * `MolFrac` - Identity function basically...
 * `RelHum` - Relative humidity in decimal notation (not percentage!)
 * `DewPoint`- Dew point temperature in K
 * `HumRat` - Humidity ration, kg of vapor / kg of dry air
 * `MassFrac` - Mass fraction of water vapor = mv / (mv + ma)
 * `WetBulb` - Wet bulb temperature in K- estimated by the adiabatic saturation temperature.
 
As usual all calculations are performed using plain SI units. 
The parameters used by this function are:

 * `Tk` - Temperature in K
 * Humidity type - a type listed above that characterizes the humidity
 * The value of the humidity parameter. 
 * P - Pressure in Pa
"""
function molarfrac(Tk, ::Type{MolFrac}, xv, P)
    xv
end

function molarfrac(Tk, ::Type{HumRat}, w, P)
    w  / (Mv/Ma + w)
end

function molarfrac(Tk, ::Type{RelHum}, rel, P)
    xv = rel * efactor(Tk, P) * Pws(Tk) / P
end

function molarfrac(Tk, ::Type{DewPoint}, D, P)
    xv = efactor(D,P) * Pws(D) / P
end

function molarfrac(Tk, ::Type{WetBulb}, B, P)
    w = calc_W_from_B(Tk, B, P)
    w / (Mv/Ma + w)
end

function volume(::Type{MoistAir}, Tk, ::Type{T}, y, P) where {T<:PsychroProperty}
    xv = molarfrac(Tk, T, y, P)
    volumemoist(Tk, P, xv)
end

function molarvolume(::Type{MoistAir}, Tk, ::Type{T}, y, P) where {T<:PsychroProperty}
    xv = molarfrac(Tk, T, y, P)
    molarvolumemoist(Tk, P, xv)
end

function density(::Type{MoistAir}, Tk, ::Type{T}, y, P) where {T<:PsychroProperty}
    xv = molarfrac(Tk, T, y, P)
    M = xv*Mv + (1-xv)*Ma
    vm = molarvolumemoist(Tk, P, xv)
    return M/vm
end


function enthalpy(::Type{MoistAir}, Tk, ::Type{T}, y, P) where {T<:PsychroProperty}
    xv = molarfrac(Tk, T, y, P)
    enthalpymoist(Tk, P, xv)
end

function molarenthalpy(::Type{MoistAir}, Tk, ::Type{T}, y, P) where {T<:PsychroProperty}
    xv = molarfrac(Tk, T, y, P)
    molarenthalpymoist(Tk, P, xv)
end

function entropy(::Type{MoistAir}, Tk, ::Type{T}, y, P) where {T<:PsychroProperty}
    xv = molarfrac(Tk, T, y, P)
    entropymoist(Tk, P, xv)
end

function molarentropy(::Type{MoistAir}, Tk, ::Type{T}, y, P) where {T<:PsychroProperty}
    xv = molarfrac(Tk, T, y, P)
    molarentropymoist(Tk, P, xv)
end

function compressfactor(::Type{MoistAir}, Tk, ::Type{T}, y, P) where {T<:PsychroProperty}
    xv = molarfrac(Tk, T, y, P)
    Zmoist(Tk, P, xv)
end



function calc_W_from_B(Tk, B, P, EPS=1e-8, MAXITER=100)

    xsv = efactor(B,P) * Pws(B) / P
    w2 = Mv/Ma * xsv/(1-xsv)

    w = (h_a_(B) - h_a_(T) - w2*(h_f_(B) + h_v_(B))) / (h_v(T) - h_f_(B))
    for iter = 1:MAXITER
        f = aux_WB(w, T, B, P)
        df = (aux_WB(w+1e-5*w2, T, B, P) - f) / (1e-5*w2)
        dw = -f/df
        w = w + dw

        if abs(dw) < EPS
            return w
        end

    end

    return w
end

function aux_WB(w, T, B, P)
    xv1 = w / (Mv/Ma+w)
    xv2 = e_factor(B, P) * Pws(B) / P
    w2 = Mv / Ma * xv2 / (1.0-xv2)
    (1.0+w)*h_(T,P,xv1) + (w2-w)*h_f_(B) - (1.0+w2)*h_(B,P,xv2)
end


