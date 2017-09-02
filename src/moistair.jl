#=
This file actually implements the user interface of the Psychro library.
=#

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


