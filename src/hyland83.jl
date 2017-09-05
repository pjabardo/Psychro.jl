#=

 This file implements the functions presented in the paper

 * [1] Wexler, A. and Hyland, R. W., "Formulations for the thermodynamic properties of the saturated phases of H2O from 173.15 K to 473.15 K", ASHRAE Transactions, 1983.

=#    

"""
    volumeice(Tk)

Specific volume of saturated ice. Equation 2 from ref. [1].

 * `Tk` Temperature in K.
 * Output: specific volume in m^3/kg
"""
volumeice(Tk) = 0.1070003e-2 + Tk*(-0.249936e-7 + 0.371611e-9*Tk)


"""
    enthalpyice(Tk)

Specific enthalpy of saturated ice. Equation 3 of ref. [1].

 * `Tk` Temperature in K
 * Output: J/kg

"""
function enthalpyice(Tk)

    -0.647595E6 + Tk*(0.274292e3 + Tk*(0.2910583e1 + Tk*0.1083437e-2)) + 0.107e-2*Pws_s(Tk)
    
end

"""
    densitywater(Tk)

Density of saturated water. Equation 5 from ref. [1].

 * `Tk` Temperature in K.
 * Output: density kg/m^3
"""
densitywater(Tk) = ( -0.2403360201e4 +
                     Tk*(-0.140758895e1 +
                         Tk * (0.1068287657e0 +
                               Tk*(-0.2914492351e-3 + Tk*(0.373497936e-6 - Tk*0.21203787e-9))))) /
                               ( -0.3424442728e1 + 0.1619785e-1*Tk )

"""
    volumewater(Tk)

Specific volume of saturated water. Inverse of equation 5 from ref. [1].

 * `Tk` Temperature in K.
 * Output: specific volume in m^3/kg
"""
volumewater(Tk) = 1.0 / densitywater(Tk)



const L = (-0.11411380e7, 0.41930463e4, -0.8134865e-1,
           0.1451133e-3, -0.1005230e-6, -0.563473e3, -0.036)
const M = (-0.1141837121e7, 0.4194325677e4, -0.6908894163e-1,
           0.105555302e-3, -0.7111382234e-7, 0.6059e-3)

"""
    enthalpywater(Tk)

Specific enthalpy of saturated water. Equations 6-11 of ref. [1].

 * `Tk` Temperature in K
 * Output: J/kg

"""
function enthalpywater(Tk)
    β₀ = Tk * volumewater(273.16) * dPws_l(273.16)
    β = Tk * volumewater(Tk) * dPws_l(Tk) 

    if Tk < 373.125
        α = L[1] + Tk*(L[2] + Tk*(L[3] + Tk*(L[4] + Tk*L[5]))) +
            L[6] * 10^(L[7] * (Tk - 273.16))
    else 
        α = M[1] + Tk*(M[2] + Tk*(M[3] + Tk*(M[4] + Tk*M[5])))
        if Tk > 403.128
            α = α - M[6]*(Tk - 403.128)^3.1
        end
    end
    
    return α + β - β₀
end


"""
    enthalpywi(Tk)

Specific enthalpy of saturated water in condensed phase (ice or water).

 * `Tk` Temperature in K
 * Output: J/kg

"""
function enthalpywi(Tk)
    if Tk >= 273.15
        return enthalpywater(Tk)
    else
        return enthalpyice(Tk)
    end
end



                               


"""
    ```Blin(Tk)```

Virial coefficient B' saturated vapor eq 15 [1]

 * `Tk` Temperature in K
 * Output: B' in Pa^(-1)
"""
Blin(Tk) = 0.70e-8 - 0.147184e-8 * exp(1734.29/Tk) # Pa^(-1)

"""
Second virial coefficient C' saturated vapor eq 16 [1]

 * `Tk` Temperature in K
 * Output: B' in Pa^(-2)
"""
Clin(Tk) = 0.104e-14 - 0.335297e-17*exp(3645.09/Tk) # Pa^(-2)


"""
    ```dBlin(Tk)```

Derivative of virial coefficient dB'/dT saturated vapor eq 15 [1]

 * `Tk` Temperature in K
 * Output: B' in Pa^(-1)/K
"""
dBlin(Tk) = 2.5525974e-6/(Tk*Tk) * exp(1734.29/Tk)


"""
    ```dClin(Tk)```

Derivative of Second virial coefficient dC'/dT saturated vapor eq 16 [1]

 * `Tk` Temperature in K
 * Output: B' in Pa^(-2)/K
"""
dClin(Tk) = 1.2221877e-14/(Tk*Tk) * exp(3645.09/Tk)




const gg = (-0.58002206e4, 0.13914993e1, -0.48640239e-1,  0.41764768e-4, -0.14452093e-7, 0.65459673e1)
const mm = (-0.56745359e4, 0.63925247e1, -0.96778430e-2,  0.62215701e-6,  0.20747825e-8,
           -0.94840240e-12, 0.41635019e1)


"""
    ```Pws_l(Tk)```

Saturation pressure of vapor pressure over liquid water.
This implements equation 17 from [1].

 * `Tk` Temperature in K
 * Output: Pa
"""
function Pws_l(Tk)
    
    exp((gg[1] + Tk*(gg[2] + Tk*(gg[3] + Tk*(gg[4]+Tk*gg[5]))))/Tk + gg[6]*log(Tk))
                    
end



"""
    ```Pws_s(Tk)```

Saturation pressure of vapor pressure over ice.
This implements equation 18 from [1].

 * `Tk` Temperature in K
 * Output: Pa
"""
function Pws_s(Tk)

  
    exp( (mm[1] + Tk*(mm[2] + Tk*(mm[3] + Tk*(mm[4] + Tk*(mm[5] + Tk*mm[6])))))/Tk + mm[7] * log(Tk) )
end


"""
    ```Pws(Tk)```

Saturation pressure of vapor pressure over liquid water or ice.
This function calls either `Pws_l` or `Pws_s`. At a temperature of 
273.16 both expressions are almost exactly the same (6 decimal figures).

 * `Tk` Temperature in K
 * Output: Pa
"""
function Pws(Tk)
    if Tk < 273.15
        Pws_s(Tk)
    elseif Tk > 273.16
        Pws_l(Tk)
    else
        p1 = Pws_s(273.15)
        p2 = Pws_l(273.16)
        p2 * (Tk-273.15)/0.01 + p1 * (273.16-Tk)/0.01
    end
        
end
   

"""
    ```dPws_s(Tk)```

Derivative of saturation pressure of vapor pressure over ice.
This implements the derivative of equation 18 from [1].

 * `Tk` Temperature in K
 * Output: Pa/K
"""
function dPws_s(Tk)
    x1 = Pws_s(Tk)
    x2 = 1.0/Tk*(mm[7] - mm[1]/Tk) + mm[3] + Tk*(2*mm[4] + Tk*(3*mm[5] +4*mm[6]*Tk))
  
    return x1*x2
end

"""
    ```dPws_s(Tk)```

Derivative of saturation pressure of vapor pressure over water.
This implements the derivative of equation 17 from [1].

 * `Tk` Temperature in K
 * Output: Pa/K
"""
function dPws_l(Tk)
    x1 = Pws_l(Tk)
    x2 = 1.0/Tk*(gg[6] - gg[1]/Tk) + gg[3] + Tk*(2*gg[4] + 3*gg[5]*Tk)
    return x1*x2 
end


"""
    ```dPws(Tk)```

Derivative of saturation pressure of vapor pressure over liquid water and ice.
This combines the functions `dPws_l` and `dPws_s`.

 * `Tk` Temperature in K
 * Output: Pa/K
"""
function dPws(Tk)
    if Tk < 273.16
        dPws_s(Tk)
    else Tk 
        dPws_l(Tk)
    end
end
  

const hh = (2.127925e2, 7.305398e0, 1.969953e-1, 1.103701e-2, 1.849307e-3, 5.145087e-6)

"""
    ```Tws(P)```

Calculates the saturation temperature of water vapor. 
This function is the inverse of function `Pws(T)`. 
First an approximation was obtained and then a Newton
iteration is used to obtain more accurate data.

 * `P` Saturation pressure in Pa
 * Output: Saturation temperature in K.
"""
function Tws(P)
    EPS=1e-8
    MAXITER=200

    lnP = log(P)
    T = hh[1] + lnP * (hh[2] + lnP*(hh[3] + lnP*(hh[4] + lnP*hh[5]))) + hh[6] * P

    i = 0
    dT = 0.0
    for i = 0:MAXITER
        f = P - Pws(T)
        df = -dPws(T)
        dT = -f / df
        T += dT

        if abs(dT) < EPS
            return T
        end
    end
    throw(ConvergenceError("Calculation of saturation temperature failed to converge!", T, i, dT))
    return T
end



"""
    enthalpyvapor(Tk)

Specific enthalpy of saturated water vapor. Equation 19 of ref. [1].

Equation 19 ha a problem since the terms using the virial coefficients are actually
molar enthalpies not specific. Dividing these terms by Mv solves any issues and
the equations match the table in the appendix!

 * `Tk` Temperature in K
 * Output: J/kg

"""
function enthalpyvapor(Tk)

    termo1 = 0.199798e7 + Tk*(0.18035706e4 +
                              Tk*(0.36400463e0 +
                                  Tk*(-0.14677622e-2 + Tk*(0.28726608e-5 - Tk*0.17508262e-8))))
    p = Pws(Tk)
    return termo1 - (R*Tk*Tk*p*(dBlin(Tk) + 0.5*dClin(Tk)*p))/Mv
end



"""
    volumevapor(Tk)

Molar volume of saturated vapor. Eq. 21 reference [1].

 * `Tk` Temperature in K
 * Output: m^3/mol
"""
function molarvolvapor(Tk)

    p = Pws(Tk)
    return R*Tk/p * (1.0 + p*(Blin(Tk)  + p*Clin(Tk)))
end

"""
    Zvapor(Tk)

Compressibility factor of saturated vapor. Eq. 21 reference [1].

 * `Tk` Temperature in K
 * Output: m^3/mol
"""
function Zvapor(Tk)
    p = Pws(Tk)
    (1.0 + p*(Blin(Tk)  + p*Clin(Tk)))
end

"""
    volumevapor(Tk)

Specific volume of saturated vapor. Eq. 21 reference [1].

 * `Tk` Temperature in K
 * Output: m^3/kg
"""
function volumevapor(Tk)

    p = Pws(Tk)
    return R*Tk/p * (1.0 + p*(Blin(Tk)  + p*Clin(Tk))) / Mv
end

const ww = (1.92703e3, 0.10959485e1, -0.3892708e-2,
            0.82520238e-5, -0.80622878e-8,
            0.30897984e-11, 1.7920705e3)


"""
    molarentropyvapor(Tk)

Molar entropy of saturated vapor. Eq. 20 reference [1].

 * `Tk` Temperature in K
 * Output: J/mol/K
"""
function molarentropyvapor(Tk)
    p = Pws(Tk)
    s1 = @polyeval(Tk, ww, 6) + ww[7]*log(Tk)
    s2 = -R*log(p) - R*p * ( (Blin(Tk) + Tk*dBlin(Tk)) + p/2 * (Clin(Tk) + Tk*dClin(Tk)) )

    return s1*Mv + s2
end


"""
    entropyvapor(Tk)

Specific entropy of saturated vapor. Eq. 20 reference [1].

 * `Tk` Temperature in K
 * Output: J/kg/K
"""
entropyvapor(Tk) = molarentropyvapor(Tk) / Mv

