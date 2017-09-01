#=

 This file implements the functions presented in the paper

 * [2] Wexler, A. and Hyland, R. W., "Formulations for the thermodynamic properties of dry air from 173.15 K to 473.15 K, and of saturated moist air from 173.15 K to 372.15 K at pressures to 5 MPa

=#    



"""
    ```Bm(Tk,xv)```

First virial coefficient Bm of moist air of saturated vapor eq 2 [2]

 * `Tk` Temperature in K
 * `xv` Molar fraction of water vapor in the moist air
 * Output: Baw in m^3/mol
"""
function Bm(Tk, xv)
    xa = 1 - xv
    xa*xa*Baa(Tk) + 2*xa*xv*Baw(Tk) + xv*xv*Bww(Tk)
end

"""
    ```Cm(Tk,xv)```

Second virial coefficient Cm of moist air of saturated vapor eq 3 [2]

 * `Tk` Temperature in K
 * `xv` Molar fraction of water vapor in the moist air
 * Output: Baw in m^6/mol^2
"""
function Cm(Tk, xv)
    xa = 1-xv
    xa*xa*xa*Caaa(Tk) + 3*xa*xa*xv*Caaw(Tk) + 3*xa*xv*xv*Caww(Tk) + xv*xv*xv*Cwww(Tk)
end



"""
    ```Baa(Tk)```

Virial coefficient Bₐₐ of dry air eq 10 [2]

 * `Tk` Temperature in K
 * Output: Bₐₐ in m^3/mol

"""
Baa(Tk) = 0.349568e-4 + (1.0/Tk)*(-0.668772e-2 + (1.0/Tk) * (-0.210141e1 + 0.924746e2/Tk))

"""
    ```dBaa(Tk)```

Derivative of Bₐₐ(T). Calculated from equation 10 [2].

 * `Tk` Temperature in K
 * Output: Bₐₐ in m^3/mol/K

"""
dBaa(Tk) = 1.0/(Tk*Tk) * (0.668772e-2 + (1.0/Tk)*(0.420282e1 - 2.774238e2/Tk))

"""
    ```Caa(Tk)```

Second virial coefficient Cₐₐₐ of dry air eq 11 [2].

 * `Tk` Temperature in K
 * Output: Bₐₐ in m^6/mol^2

"""
Caaa(Tk) = (0.125975e-8 + 1.0/Tk * (-0.190905e-6 + 0.632467e-4/Tk))

"""
    ```dCaa(Tk)```

Derivative of second virial coefficient Cₐₐₐ of dry air eq 11 [2].

 * `Tk` Temperature in K
 * Output: Bₐₐ in m^6/mol^2/K

"""
dCaaa(Tk) = 1.0/(Tk*Tk) * (0.190905e-6 - 1.264934e-4/Tk) 


"""
    molarvolair(Tk)

Molar volume of dry air. Equation 12 of reference [1].
This function requires iteration to compute the volume.

 * `Tk` Temperature in K
 * `P` Pressure in Pa
 * `EPS`: Acceptable error
 * `MAXITER`: Maixmum number of iterations
 * Output: molar volume in m^3/mol
"""
function molarvolair(Tk, P, EPS=1e-8, MAXITER=100)

    RT=R*Tk
    va = RT/P

    conv = false
    b = Baa(Tk)
    c = Caaa(Tk)
    z0 = 1.0
    z=1.0
    niter = 0
    for i = 1:MAXITER
        z = 1.0 + 1/va*(b + c/va)
        va = RT*z/P
        if abs(z-z0) < EPS
            conv = true
            niter = i
            break
        end
        z0 = z
    end
    return va
end

"""
    volumeair

Specific volume of dry air. Uses `molarvolair`.

 * `Tk` Temperature in K
 * `P` Pressure in Pa
 * `EPS`: Acceptable error
 * `MAXITER`: Maixmum number of iterations
 * Output: molar volume in m^3/kg
"""
volumeair(Tk, P, EPS=1e-8, MAXITER=100) = molarvolair(Tk, P, EPS, MAXITER) / Ma


"""
    enthalpyair(Tk)

Specific enthalpy of dry air. Equation 13 of reference [2]

 * `Tk` Temperature in K
 * `P` Pressure in Pa
 * `EPS`: Acceptable error
 * `MAXITER`: Maixmum number of iterations
 * Output: m^3/kg
"""
function enthalpyair(Tk, P, EPS=1e-8, MAXITER=100)
    h1 = -0.79078691e4 + Tk*(0.28709015e2 +
                             Tk*(0.26431805e-2 +
                                 Tk*(-0.10405863e-4 +
                                     Tk*(0.18660410e-7 - 0.97843331e-11*Tk))))
    va = molarvolair(Tk, P, EPS, MAXITER)

    h2 = R*Tk/va * ( (Baa(Tk) - Tk*dBaa(Tk)) + (Caaa(Tk) - 0.5*Tk*dCaaa(Tk))/va )

    return (h1 + h2)/Ma
    
end


const ℓ =  (-0.16175159e3, 0.52863609e-2, -0.15608795e-4,
            0.24880547e-7, -0.12230416e-10, 0.28709015e2)
            

"""
    entropyair(Tk)


Specific entropy of dry air. Equation 14 of reference [2]

 * `Tk` Temperature in K
 * `P` Pressure in Pa
 * `EPS`: Acceptable error
 * `MAXITER`: Maixmum number of iterations
 * Output: m^3/(kg.K)
"""
function entropyair(Tk, P, EPS=1e-8, MAXITER=100)

    s1 = @evalpoly2(Tk, ℓ, 5) + ℓ[6]*log(Tk) - R*log(P/101325.0)
    va = molarvolair(Tk, P, EPS, MAXITER)
    s2 = R*log(P*va/R/Tk) - R/va * ( (Baa(Tk) + Tk*dBaa(Tk)) + 0.5/va * (Caaa(Tk) + Tk*dCaaa(Tk)))

    return (s1 + s2) / Ma
end


"""
    ```Baw(Tk)```

Cross virial coefficient Baw of saturated vapor eq 15 [2]

 * `Tk` Temperature in K
 * Output: Baw in m^3/mol
"""
Baw(Tk) = 0.32366097e-4 + 1.0/Tk * (-0.141138e-1 + 1.0/Tk * (-0.1244535e1 - 0.2348789e4/(Tk*Tk)))


"""
    ```Baw(Tk)```

Derivative of cross virial coefficient Bww of saturated vapor eq 15 [2]

 * `Tk` Temperature in K
 * Output: dBaw/dT in m^3/mol/K
"""
dBaw(Tk) = 1.0/(Tk*Tk) * (0.141138e-1 + 1.0/Tk * (0.248907e1 + 0.93951568e4 /(Tk*Tk)))

"""
    ```Caaw(Tk)```

Cross virial coefficient Caaw of saturated vapor eq 16 [2]

 * `Tk` Temperature in K
 * Output: Caaw in m^6/mol^2
"""
Caaw(Tk) = 0.482737e-9 + 1.0/Tk * (0.105678e-6 + 1.0/Tk *
                                   (-0.656394e-4 + 1.0/Tk*
                                    (0.294442e-1 - 0.319317e1/Tk)))

"""
    ```dCaaw(Tk)```

Derivative of cross virial coefficient dCaaw/dT of saturated vapor eq 16 [2]

 * `Tk` Temperature in K
 * Output: dCaaw/dT in m^6/mol^2/K
"""
dCaaw(Tk) = 1.0/(Tk*Tk) * (-0.105678e-6 + 1.0/Tk*(1.312788e-4 + 1.0/Tk*(-0.883326e-1 + 1.277268e1/Tk)))


"""
    ```Caww(Tk)```

Cross virial coefficient Caaw of saturated vapor eq 17 [2]

 * `Tk` Temperature in K
 * Output: Caww in m^6/mol^2
"""
Caww(Tk) = -1e-6 * exp( -0.10728876e2 + 1.0/Tk * (0.347802e4 + 1.0/Tk*(-0.383383e6 +  0.33406e8/Tk)))

"""
    ```dCaww(Tk)```

Derivative of cross virial coefficient dCaaw/dT of saturated vapor eq 17 [2]

 * `Tk` Temperature in K
 * Output: dCaww/dT in m^6/mol^2/K
"""
dCaww(Tk) = 1.0/(Tk*Tk) * (-0.347802e4 + 1.0/Tk * (2*0.383383e6 - 3*0.33406e8/Tk)) * Caww(Tk)


"""

    ```e_factor(Tk,P)```

Enhancement factor of moist air. 
Uses an iterative procedure to calculate the Enhancement factor defined 
in equation 18 [2]


"""
function efactor(Tk, P)
    f = 1.0
    EPS = 1e-7
    NMAX = 100
    for iter = 1:NMAX
        xas = (P-f*Pws(Tk)) / P
        fnovo = exp(lnf(Tk,P,xas))

        if abs(fnovo - 1.0) < EPS
            if fnovo < 1
                fnovo = 1.0
            end
            fnovo = 1.0
        end

        f = fnovo
    end
    if fnovo < 1.0
        fnovo = 1.0
    end
    fnovo
end


"""
    ```lnf(Tk, P, xas)```

Auxiliary function used to calculate the enhancement factor.
Actually, this function implements the RHS of eq. 18 of [2].
"""
function lnf(Tk, P, xas)
    if Tk < 273.16
        vc = volumeice(Tk) * Mv
        p = Pws_s(Tk)
        k = 0.0
    else
        vc = volumewater(Tk) * Mv
        p = Pws_l(Tk)
        k = henryk(Tk)
    end
    
    κ = kappa_f(Tk)

    RT = R*Tk
    RT2 = RT*RT
    xas2 = xas*xas
    p2 = p*p
    P2 = P*P
    
    baa = Baa(Tk)
    bww = Bww(Tk)
    baw = Baw(Tk)
    caaa = Caaa(Tk)
    caaw = Caaw(Tk)
    caww = Caww(Tk)
    cwww = Cwww(Tk)
    
    t1 = vc/RT * ( (1 + κ*p)*(P-p) - 0.5*κ*(P2 - p2) )
    t2 = log(1.0 - k*xas*P) + (xas2*P/RT)*baa - (2*xas2*P/RT)*baw
    t3 = -(P-p-xas2*P)/RT*bww + xas2*xas*P*P/RT2 * caaa
    t4 = 3*xas2*(1-2xas)*P2/(2*RT2) * caaw - (3xas2*(1-xas)*P2)/RT2 * caww
    t5 = -( (1+2xas)*(1-xas)^2*P2 - p2 )/(2*RT2) * cwww - (xas2*(1-3xas)*(1-xas)*P2)/RT2 * baa*bww
    t6 = -(2xas2*xas*(2-3xas)*P2) / RT2 * baa*baw + ( 6xas2*(1-xas)^2*P2 )/RT2 * bww*baw
    t7 = -3xas2*xas2*P2/(2RT2)*baa*baa - ( 2xas2*(1-xas)*(1-3xas)*P2 )/RT2 * baw*baw
    t8 = -( p2 - (1+3xas)*(1-xas)^3*P2 ) / (2*RT2) * bww*bww
    

   return t1+t2+t3+t4+t5+t6+t7+t8
   
end





"""
    ```Bww(Tk)```

Virial coefficient Bww of saturated vapor eq 19 [2]

 * `Tk` Temperature in K
 * Output: Bww in m^3/mol
"""
Bww(Tk) = R * Tk * Blin(Tk)

"""
    ```dBww(Tk)```

Derivative of virial coefficient dBww/dT of saturated vapor eq 19 [2]

 * `Tk` Temperature in K
 * Output: dBww/dT in m^3/mol/K
"""
dBww(Tk) = R * (Tk * dBlin(Tk) + Blin(Tk))


"""
    ```Cwww(Tk)```

Second virial coefficient Cwww of saturated vapor eq 20 [2]

 * `Tk` Temperature in K
 * Output: Cwww in m^3/mol
"""
Cwww(Tk) = R*R*Tk*Tk * (Clin(Tk) + Blin(Tk)^2)


"""
    ```dCwww(Tk)```

Derivative of second virial coefficient dCww/dT of saturated vapor eq 20 [2]

 * `Tk` Temperature in K
 * Output: dCwww/dT in m^3/mol/K
"""
dCwww(Tk) = R*Tk*R*Tk * (dClin(Tk) + 2*Blin(Tk) * dBlin(Tk)) + (2*R*R*Tk) * (Clin(Tk) + Blin(Tk)^2)


"""
    kappa_l(Tk)

Isothermal compressibility of saturated liquid water. Equation 21, ref. [2].
Range of applicability 0°C - 150°C but can be extended up to 200°C
More details in reference [4].
 * `Tk` Temperature in K
 * Output: κ in 1/Pa
"""
function kappa_l(Tk)
    Tc = Tk - 273.15

    #@check_range Tc 0 200
    if Tc < 100.0
        k = (50.88496 + Tc*(0.6163813 +
                            Tc*(1.459187e-3 +
                                Tc*(20.08438e-6 + Tc*(-58.47727e-9 + Tc*0.4104110e-9))))) /
                                (1.0 + 0.1967348e-1*Tc)
        
                                    
    else
        k = (50.884917 + Tc*(0.62590623 +
                            Tc*(1.3848668e-3 +
                                Tc*(21.603427e-6 + Tc*(-0.72087667e-7 + Tc*0.46545054e-9))))) /
                                (1.0 + 0.1967348e-1*Tc)
    end

    k*1e-11
end

"""
    kappa_s(Tk)

Isothermal compressibility of ice. Equation 22, ref. [2].

 * `Tk` Temperature in K
 * Output: κ in 1/Pa
"""
kappa_s(Tk) = (8.875 + 0.0165*Tk)*1e-11

kappa_f(Tk) = 
    if Tk < 273.16
        kappa_s(Tk)
    else
        kappa_l(Tk)
    end



"""
    henryk_O2(Tk)

Henry's coefficient for solubility of O2 in water.
Implements equation 23 from ref. [2]. See reference [3]
Valid for temperatures between 273.15 K and 500 K

 * `Tk` Temperature in K
 * Output: k in atm/mol * 10^(-4)
"""
function henryk_O2(Tk)

    τ = 1000.0/Tk
    α = -0.0005943
    β = -0.1470
    γ = -0.05120
    δ = -0.1076
    ɛ = 0.8447
    a2 = α;
    a1 = γ*τ + δ
    a0 = β*τ*τ + ɛ*τ - 1.0

    10^( (-a1-sqrt(a1*a1 - 4.0*a2*a0)) / (2.0*a2))
end


"""
    henryk_N2(Tk)

Henry's coefficient for solubility of N2 in water.
Implements equation 23 from ref. [2]. See reference [3]
Valid for temperatures between 273.15 K and 500 K

 * `Tk` Temperature in K
 * Output: k in atm/mol * 10^(-4)
"""
function henryk_N2(Tk)

    τ = 1000.0/Tk
    α = -0.1021
    β = -0.1482
    γ = -0.019
    δ = -0.03741
    ɛ = 0.851
    a2 = α;
    a1 = γ*τ + δ
    a0 = β*τ*τ + ɛ*τ - 1.0

    10^( (-a1-sqrt(a1*a1 - 4.0*a2*a0)) / (2.0*a2))
end


"""
    henryk(Tk)

Henry's coefficient for solubility of air in water.
Implements equations 24 and 25 from ref. [2]. See reference [3]
Valid for temperatures between 273.15 K and 500 K

 * `Tk` Temperature in K
 * Output: k in atm/mol fract * 10^(-4)
"""
function henryk(Tk)

    kO2 = henryk_O2(Tk)
    kN2 = henryk_N2(Tk)
    xO2 = 0.22
    xN2 = 0.78

    k = 1.0 / (xO2/kO2 + xN2/kN2)

    1e-4/k * 1.0 /101325.0
end




"""
    ```molarvol```

Molar volume of moist air given the molar fraction of water vapor.
Assumes a real gas.

 * `Tk` Temperature in K
 * `P` Pressure in Pa
 * `xv`  Molar fraction of water vapor
 * Output: molar volume in m^3/mol
"""
molarvol(Tk, P, xv) = Z(Tk, P, xv) * R * Tk / P
vM_(Tk, P, xv) = Z(Tk, P, xv) * R * Tk / P

"""
Specific volume of saturated water vapor in m^3/kg
"""
volumevapor1(Tk) = vM_v_(Tk) / Mv

v_v_(Tk) = vM_v_(Tk) / Mv

"""
Density of saturated water vapor in kg/m^3
"""
densityv(Tk) = 1/v_v_(Tk)
r_v_(Tk) = 1/v_v_(Tk)

"""
    ```Z(Tk, P, xv)```

Compressibility factor of moist air.

 * `Tk` Temperature in K
 * `P` Pressure in Pa
 * `xv` molar fraction of water vapor.

"""
function Z(Tk, P, xv)

    xa = 1-xv
    vmi = R*Tk/P
    vm = vmi
    b = Bm(Tk, xv)
    c = Cm(Tk, xv)
    NMAX = 100
    EPS = 1e-8

    for iter = 0:NMAX
        vmn = R*Tk/P * (1 + b/vm + c/(vm*vm))
        erro = abs(vmn-vm)
        vm = vmn

        if erro < EPS
            return vm/vmi
        end
    end

    vm/vmi
end




abstract type AbstractPsychro end

immutable Ashrae <: AbstractPsychro
    Tmin::Float64
    Tmax::Float64
    Pmin::Float64
    Pmax::Float64
    xv::Float64
    W::Float64
    M::Float64
end

function Ashrae(xv)
    w = Mv/Ma * xv/(1.0 - xv)
    m = xv*Mv + (1.0-xv)*Ma

    Ashrae(173.15, 473.15, 0.0, 5.0e6, xv, w, m)
end

Ashrae() = Ashrae(173.15, 473.15, 0.0, 5.0e6, 0.0, 0.0, Ma)


function Ashrae(ch::Char, T, humidity, P)
    a = Ashrae()
    ch = uppercase(ch)
    if ch == 'X'
        xv = humidity
        xsv = e_factor(T, P) * Pws(T) / P
        w = Mv/Ma * xv/(1.0 - xv)
    elseif ch == 'W'
        w = humidity
        xv = w / (Mv/Ma + w)
        xsv = e_factor(T,P) * Pws(T) / P
    elseif ch == 'R'
        rel = humidity
        xv = rel * e_factor(T, P) * Pws(T) / P
        w = Mv/Ma * xv/(1.0 - xv)
    elseif ch == 'D'
        D = humidity
        xv = e_factor(D,P) * Pws(D) / P
    elseif ch == 'B'

        B = humidity
        w = calc_W_from_B(T, B, P)
        xv = w / (Mv/Ma + w)
    elseif ch=='D'
        D = humidity

        xv = e_factor(D,P) * Pws(D) / P
        w = Mv/Ma * xv/(1.0 - xv)

    end

    return Ashrae(xv)
    
        
        
end


function calc_W_from_B(T, B, P)

    NMAX = 100
    xsv = e_factor(B,P) * Pws(B) / P
    w2 = Mv/Ma * xsv/(1-xsv)

    EPS = 1e-8*w2
    w = (h_a_(B) - h_a_(T) - w2*(h_f_(B) + h_v_(B))) / (h_v(T) - h_f_(B))
    for iter = 1:NMAX
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


    
function vM_a_(Tk, P)
    xa = 1.0
    vmi = R*Tk/P
    vm = vmi
    b = Baa(Tk)
    c = Caaa(Tk)
    NMAX = 100
    EPS = 1e-8
    for iter = 1:NMAX
        vmn = R*Tk/P * (1.0 + 1.0/vm*(b + c/vm))
        erro = abs(vmn-vm)
        vm = vmn

        if erro < EPS
            return vm
        end
    end
    return vm
end
   
v_a_(Tk, P) = vM_a_(Tk,P) / Ma
r_a_(Tk, P) = 1.0 / v_a_(Tk,P)

function h_a_(Tk, P)
    b = [-0.79078691e4,
         0.28709015e2,
         0.26431805e-2,
         -0.10405863e-4,
         0.18660410e-7,
	 -0.97843331e-11]

    B = Baa(Tk)
    C = Caaa(Tk)
    dB = dBaa(Tk)
    dC = dCaaa(Tk)
  
    Vm = vM_a_(Tk, P)

    ha = 1000.0*(b[1] + b[2]*Tk + b[3]*Tk*Tk + b[4]*Tk*Tk*Tk +b[5]*Tk^4 + b[6]*Tk^5)
  
    ha = ha + R*Tk * ( (B - Tk*dB)/Vm + (C - 0.5*Tk*dC)/(Vm*Vm)  )
    return ha/Ma;
end


function vM_v_(Tk)
    P = Pws(Tk)
    vmi = R*Tk/P
    vm = vmi
    b = Bww(Tk)
    c = Cww(Tk)
    NMAX = 200
    EPS = 1e-9
    for iter = 1:NMAX
        vmn = R*Tk/P * (1.0 + 1.0/vm*(b + c/vm))
        erro = abs(vmn - vm)
        vm = vmn
        if erro < EPS
            return vm
        end
    end

    return vm
end
    
v_(Tk, P, xv) = vM_(Tk, P, xv) / ((1.0-xv)*Ma + xv*Mv)
r_(Tk, P, xv) = 1.0 / v_(Tk, P, xv)


h_f_(Tk) =
    if Tk < 273.16
        h_s_(Tk)
    else
        h_l_(Tk)
    end

function h_v_(Tk)

    xv = 1.0
    xa = 0.0

    P = Pws(Tk)
    a =[0.63290874e1,
        0.28709015e2,
        0.26431805e-2,
        -0.10405863e-4,
        0.18660410e-7,
	-0.9784331e-11]

    d = [-0.5008e-2,
         0.32491829e2,
         0.65576345e-2,
         -0.26442147e-4,
         0.51751789e-7,
	 -0.31541624e-10]

    B = Bww(Tk)
    C = Cwww(Tk) 

    dB = dBww(Tk)
    dC = dCwww(Tk)
    
    
    hv = 35994.17
    
    termo2 = d[1] + d[2]*Tk + d[3]*Tk*Tk + d[4]*Tk*Tk*Tk +
    d[5]*Tk^4 + d[6]*Tk^5 + hv
    

    # Cálculo do volume molar
    Vm = vM_v_(Tk)
    
    termo3 = (B - Tk*dB)/Vm + (C - 0.5*Tk*dC)/(Vm*Vm)

    hm =  termo2 * 1000.0 + R*Tk*termo3
  
    hm/Mv
end

    

    
        
    

    

    
