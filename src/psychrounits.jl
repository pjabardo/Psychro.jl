#=
Implement the functions using the Unitful package to check units
=#

using Unitful

const uT = u"K"
const uP = u"Pa"
const uV = u"m^3/kg"
const umV = u"m^3/mol"
const uD = u"kg/m^3"
const umD = u"mol/m^3"
const perc = u"cm/m"
const uH = u"J/kg"
const umH = u"J/mol"
const uS = u"J/kg/K"
const umS = u"J/mol/K"


dimless{T,U} = Quantity{T, Unitful.Dimensions{()}, U}
ndim = unit(dimless)


val(u, x::Quantity) = uconvert(u, x).val
val(u::Unitful.FreeUnits{(), Unitful.Dimensions{()}}, x) = uconvert(u, x)

function volume(::Type{DryAir}, Tk::Quantity, P::Quantity, u=uV)
    uconvert(u, volume(DryAir, val(uT, Tk), val(u"Pa", P))*uV)
end

function volumem(::Type{DryAir}, Tk::Quantity, P::Quantity, u=umV)
    uconvert(u, volumem(DryAir, val(uT, Tk), val(uP, P))*umV)
end

function density(::Type{DryAir}, Tk::Quantity, P::Quantity, u=uD)
    uconvert(u, density(DryAir, val(uT, Tk), val(uP, P))*uD)
end

function enthalpy(::Type{DryAir}, Tk::Quantity, P::Quantity, u=uH)
    uconvert(u, enthalpy(DryAir, val(uT, Tk), val(uP, P))*uH)
end

function enthalpym(::Type{DryAir}, Tk::Quantity, P::Quantity, u=umH)
    uconvert(u, enthalpym(DryAir, val(uT, Tk), val(uP, P))*umH)
end

function entropy(::Type{DryAir}, Tk::Quantity, P::Quantity, u=uS)
    uconvert(u, entropy(DryAir, val(uT, Tk), val(uP, P))*uS)
end

function entropym(::Type{DryAir}, Tk::Quantity, P::Quantity, u=umS)
    uconvert(u, entropym(DryAir, val(uT, Tk), val(uP, P))*umS)
end

function compressfactor(::Type{DryAir}, Tk::Quantity, P::Quantity)
    compressfactor(DryAir, val(uT, Tk), val(uP, P))
end


function volume(::Type{Vapor}, Tk::Quantity, u=uV)
    uconvert(u, volume(Vapor, val(uT, Tk))*uV)
end

function volumem(::Type{Vapor}, Tk::Quantity, u=umV)
    uconvert(u, volumem(Vapor, val(uT, Tk))*umV)
end

function density(::Type{Vapor}, Tk::Quantity, u=uD)
    uconvert(u, density(Vapor, val(uT, Tk))*uD)
end

function enthalpy(::Type{Vapor}, Tk::Quantity, u=uH)
    uconvert(u, enthalpy(Vapor, val(uT, Tk))*uH)
end

function enthalpym(::Type{Vapor}, Tk::Quantity, u=umH)
    uconvert(u, enthalpym(Vapor, val(uT, Tk))*umH)
end


function entropy(::Type{Vapor}, Tk::Quantity, u=uS)
    uconvert(u, entropy(Vapor, val(uT, Tk))*uS)
end

function entropym(::Type{Vapor}, Tk::Quantity, u=umS)
    uconvert(u, entropym(Vapor, val(uT, Tk))*umS)
end

function compressfactor(::Type{Vapor}, Tk::Quantity)
    compressfactor(Vapor, val(uT, Tk))
end

function molarfrac(Tk::Quantity, ::Type{MolarFrac}, xv::Number, P::Quantity)
    xv
end
function molarfrac(Tk::Quantity, ::Type{HumRat}, w::Number, P::Quantity)
   w / (Mv/Ma + w)
end
function molarfrac(Tk::Quantity, ::Type{SpecHum}, q::Number, P::Quantity)
    q * Ma / (Mv + q*(Ma - Mv))
end

function molarfrac(Tk::Quantity, ::Type{RelHum}, rel, P::Quantity)
    molarfrac(val(uT, Tk), RelHum, rel, val(uP, P))
end

function molarfrac(Tk::Quantity, ::Type{DewPoint}, D::Quantity, P::Quantity)
    molarfrac(Tk, DewPoint, uconvert(u"K", D).val, val(uP, P))
end

function molarfrac(Tk::Quantity, ::Type{WetBulb}, B::Quantity, P::Quantity)
    molarfrac(val(uT, Tk), WetBulb, uconvert(u"K", B).val, val(uP, P))
end

function volume(::Type{MoistAir}, Tk::Quantity,
                ::Type{T}, y, P::Quantity, u=uV) where {T<:PsychroProperty}

    xv = molarfrac(Tk, T, y, P)
    uconvert(u, volume(MoistAir, val(uT, Tk), MolarFrac, xv, val(uP, P))*uV)
end
                   

function volumem(::Type{MoistAir}, Tk::Quantity,
                 ::Type{T}, y, P::Quantity, u=umV) where {T<:PsychroProperty}

    xv = molarfrac(Tk, T, y, P)
    uconvert(u, volumem(MoistAir, val(uT, Tk), MolarFrac, xv, val(uP, P))*umV)
end

function density(::Type{MoistAir}, Tk::Quantity,
                 ::Type{T}, y, P::Quantity, u=uD) where {T<:PsychroProperty}

    xv = molarfrac(Tk, T, y, P)
    uconvert(u, density(MoistAir, val(uT, Tk), MolarFrac, xv, val(uP, P))*uD)
end

function enthalpy(::Type{MoistAir}, Tk::Quantity,
                  ::Type{T}, y, P::Quantity, u=uH) where {T<:PsychroProperty}

    xv = molarfrac(Tk, T, y, P)
    uconvert(u, enthalpy(MoistAir, val(uT, Tk), MolarFrac, xv, val(uP, P))*uH)
end

function enthalpym(::Type{MoistAir}, Tk::Quantity,
                   ::Type{T}, y, P::Quantity, u=umH) where {T<:PsychroProperty}

    xv = molarfrac(Tk, T, y, P)
    uconvert(u, enthalpym(MoistAir, val(uT, Tk), MolarFrac, xv, val(uP, P))*umH)
end
             
function entropy(::Type{MoistAir}, Tk::Quantity,
                  ::Type{T}, y, P::Quantity, u=uS) where {T<:PsychroProperty}
    
    xv = molarfrac(Tk, T, y, P)
    uconvert(u, entropy(MoistAir, val(uT, Tk), MolarFrac, xv, val(uP, P))*uS)
end

function entropym(::Type{MoistAir}, Tk::Quantity,
                  ::Type{T}, y, P::Quantity, u=umS) where {T<:PsychroProperty}
    
    xv = molarfrac(Tk, T, y, P)
    uconvert(u, entropym(MoistAir, val(uT, Tk), MolarFrac, xv, val(uP, P))*umS)
end


function compressfactor(::Type{MoistAir}, Tk::Quantity,
                        ::Type{T}, y, P::Quantity) where {T<:PsychroProperty}

    xv = molarfrac(Tk, T, y, P)
    compressfactor(MoistAir, val(uT, Tk), MolarFrac, xv, val(uP, P))
end



function dewpoint(::Type{MoistAir}, Tk::Quantity,
                  ::Type{T}, y, P::Quantity, u=u"°C") where {T<:PsychroProperty}

    xv = molarfrac(Tk, T, y, P)
    uconvert(u"°C", dewpoint(MoistAir, val(uT, Tk), MolarFrac, xv, val(uP, P))*uT)
end


function relhum(::Type{MoistAir}, Tk::Quantity,
                ::Type{T}, y, P::Quantity, u=ndim) where {T<:PsychroProperty}

    xv = molarfrac(Tk, T, y, P)
    uconvert(u, relhum(MoistAir, val(uT, Tk), MolarFrac, xv, val(uP, P)))
end


function wetbulb(::Type{MoistAir}, Tk::Quantity,
                 ::Type{T}, y, P::Quantity, u=u"°C") where {T<:PsychroProperty}

    xv = molarfrac(Tk, T, y, P)
    uconvert(u"°C", wetbulb(MoistAir, val(uT, Tk), MolarFrac, xv, val(uP, P))*uT)
end

function humrat(::Type{MoistAir}, Tk::Quantity,
                ::Type{T}, y, P::Quantity, u=ndim) where {T<:PsychroProperty}

    xv = molarfrac(Tk, T, y, P)
    uconvert(u, humrat(MoistAir, val(uT, Tk), MolarFrac, xv, val(uP, P)))
end

function molarfrac(::Type{MoistAir}, Tk::Quantity,
                   ::Type{T}, y, P::Quantity) where {T<:PsychroProperty}

    molarfrac(Tk, T, y, P)
end

function spechum(::Type{MoistAir}, Tk::Quantity,
                ::Type{T}, y, P::Quantity) where {T<:PsychroProperty}

    xv = molarfrac(Tk, T, y, P)
    spechum(MoistAir, val(uT, Tk), MolarFrac, xv, val(uP, P))
end
