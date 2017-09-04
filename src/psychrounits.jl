#=
Implement the functions using the Unitful package to check units
=#

using Unitful

function volume(::Type{DryAir}, Tk::Quantity, P::Quantity)
    volume(DryAir, uconvert(u"K", Tk).val, uconvert(u"Pa", P).val)*1.0u"m^3/kg"
end

function molarvolume(::Type{DryAir}, Tk::Quantity, P::Quantity)
    molarvolume(DryAir, uconvert(u"K", Tk).val, uconvert(u"Pa", P).val)*1.0u"m^3/mol"
end

function density(::Type{DryAir}, Tk::Quantity, P::Quantity)
    density(DryAir, uconvert(u"K", Tk).val, uconvert(u"Pa", P).val)*1.0u"kg/m^3"
end

function enthalpy(::Type{DryAir}, Tk::Quantity, P::Quantity)
    enthalpy(DryAir, uconvert(u"K", Tk).val, uconvert(u"Pa", P).val)*1.0u"J/kg"
end

function molarenthalpy(::Type{DryAir}, Tk::Quantity, P::Quantity)
    molarenthalpy(DryAir, uconvert(u"K", Tk).val, uconvert(u"Pa", P).val)*1.0u"J/mol"
end

function entropy(::Type{DryAir}, Tk::Quantity, P::Quantity)
    entropy(DryAir, uconvert(u"K", Tk).val, uconvert(u"Pa", P).val)*1.0u"J/kg/K"
end

function molarentropy(::Type{DryAir}, Tk::Quantity, P::Quantity)
    molarentropy(DryAir, uconvert(u"K", Tk).val, uconvert(u"Pa", P).val)*1.0u"J/mol/K"
end

function compressfactor(::Type{DryAir}, Tk::Quantity, P::Quantity)
    compressfactor(DryAir, uconvert(u"K", Tk).val, uconvert(u"Pa", P).val)
end


function volume(::Type{Vapor}, Tk::Quantity, P::Quantity)
    volume(Vapor, uconvert(u"K", Tk).val)*1.0u"m^3/kg"
end

function molarvolume(::Type{Vapor}, Tk::Quantity)
    molarvolume(Vapor, uconvert(u"K", Tk).val)*1.0u"m^3/mol"
end

function density(::Type{Vapor}, Tk::Quantity)
    density(Vapor, uconvert(u"K", Tk).val)*1.0u"kg/m^3"
end

function enthalpy(::Type{Vapor}, Tk::Quantity)
    enthalpy(Vapor, uconvert(u"K", Tk).val)*1.0u"J/kg"
end

function molarenthalpy(::Type{Vapor}, Tk::Quantity)
    molarenthalpy(Vapor, uconvert(u"K", Tk).val)*1.0u"J/mol"
end


function entropy(::Type{Vapor}, Tk::Quantity)
    entropy(Vapor, uconvert(u"K", Tk).val)*1.0u"J/kg/K"
end

function molarentropy(::Type{Vapor}, Tk::Quantity)
    molarentropy(Vapor, uconvert(u"K", Tk).val)*1.0u"J/mol/K"
end

function compressfactor(::Type{Vapor}, Tk::Quantity)
    compressfactor(Vapor, uconvert(u"K", Tk).val)
end

function molarfrac(Tk::Quantity, ::Type{MolarFrac}, xv::Number, P::Quantity)
    xv
end
function molarfrac(Tk::Quantity, ::Type{HumRat}, w::Number, P::Quantity)
   w / (Mv/Ma + w)
end
function molarfrac(Tk::Quantity, ::Type{SpecHum}, q::Number, P::Quantity)
    q * Ma / (Mv + r*(Ma - Mv))
end

function molarfrac(Tk::Quantity, ::Type{RelHum}, rel::Number, P::Quantity)
    molarfrac(uconvert(u"K", Tk).val, RelHum, rel, uconvert(u"Pa", P).val)
end

function molarfrac(Tk::Quantity, ::Type{DewPoint}, D::Quantity, P::Quantity)
    molarfrac(Tk, DewPoint, uconvert(u"K", D).val, uconvert(u"Pa", P).val)
end

function molarfrac(Tk::Quantity, ::Type{WetBulb}, B::Quantity, P::Quantity)
    molarfrac(uconvert(u"K", Tk).val, WetBulb, uconvert(u"K", B).val, uconvert(u"Pa", P).val)
end

function volume(::Type{MoistAir}, Tk::Quantity,
                ::Type{T}, y::Quantity, P::Quantity) where {T<:PsychroProperty}

    xv = molarfrac(Tk, T, y, P)
    volume(MoistAir, uconvert(u"K", Tk).val, MolarFrac, xv, uconvert(u"Pa", P).val)*1.0u"m^3/kg"
end


function molarvolume(::Type{MoistAir}, Tk::Quantity,
                ::Type{T}, y::Quantity, P::Quantity) where {T<:PsychroProperty}

    xv = molarfrac(Tk, T, y, P)
    molarvolume(MoistAir, uconvert(u"K", Tk).val, MolarFrac, xv, uconvert(u"Pa", P).val)*1.0u"m^3/mol"
end

function density(::Type{MoistAir}, Tk::Quantity,
                 ::Type{T}, y::Quantity, P::Quantity) where {T<:PsychroProperty}

    xv = molarfrac(Tk, T, y, P)
    density(MoistAir, uconvert(u"K", Tk).val, MolarFrac, xv, uconvert(u"Pa", P).val)*1.0u"kg/m^3"
end

function enthalpy(::Type{MoistAir}, Tk::Quantity,
                ::Type{T}, y::Quantity, P::Quantity) where {T<:PsychroProperty}

    xv = molarfrac(Tk, T, y, P)
    enthalpy(MoistAir, uconvert(u"K", Tk).val, MolarFrac, xv, uconvert(u"Pa", P).val)*1.0u"J/kg"
end

function molarenthalpy(::Type{MoistAir}, Tk::Quantity,
                       ::Type{T}, y::Quantity, P::Quantity) where {T<:PsychroProperty}

    xv = molarfrac(Tk, T, y, P)
    molarenthalpy(MoistAir, uconvert(u"K", Tk).val, MolarFrac, xv, uconvert(u"Pa", P).val)*1.0u"J/mol"
end

function entroppy(::Type{MoistAir}, Tk::Quantity,
                ::Type{T}, y::Quantity, P::Quantity) where {T<:PsychroProperty}

    xv = molarfrac(Tk, T, y, P)
    entropy(MoistAir, uconvert(u"K", Tk).val, MolarFrac, xv, uconvert(u"Pa", P).val)*1.0u"J/kg/K"
end

function molarentroppy(::Type{MoistAir}, Tk::Quantity,
                ::Type{T}, y::Quantity, P::Quantity) where {T<:PsychroProperty}

    xv = molarfrac(Tk, T, y, P)
    molarentropy(MoistAir, uconvert(u"K", Tk).val, MolarFrac, xv, uconvert(u"Pa", P).val)*1.0u"J/mol/K"
end


function compressfactor(::Type{MoistAir}, Tk::Quantity,
                        ::Type{T}, y::Quantity, P::Quantity) where {T<:PsychroProperty}

    xv = molarfrac(Tk, T, y, P)
    compressfactor(MoistAir, uconvert(u"K", Tk).val, MolarFrac, xv, uconvert(u"Pa", P).val)
end



function dewpoint(::Type{MoistAir}, Tk::Quantity,
                  ::Type{T}, y::Quantity, P::Quantity) where {T<:PsychroProperty}

    xv = molarfrac(Tk, T, y, P)
    dewpoint(MoistAir, uconvert(u"K", Tk).val, MolarFrac, xv, uconvert(u"Pa", P).val)*1.0u"K"
end


function relhum(::Type{MoistAir}, Tk::Quantity,
                  ::Type{T}, y::Quantity, P::Quantity) where {T<:PsychroProperty}

    xv = molarfrac(Tk, T, y, P)
    relhum(MoistAir, uconvert(u"K", Tk).val, MolarFrac, xv, uconvert(u"Pa", P).val)
end


function wetbulb(::Type{MoistAir}, Tk::Quantity,
                  ::Type{T}, y::Quantity, P::Quantity) where {T<:PsychroProperty}

    xv = molarfrac(Tk, T, y, P)
    wetbulb(MoistAir, uconvert(u"K", Tk).val, MolarFrac, xv, uconvert(u"Pa", P).val)*1.0u"K"
end

function humrat(::Type{MoistAir}, Tk::Quantity,
                ::Type{T}, y::Quantity, P::Quantity) where {T<:PsychroProperty}

    xv = molarfrac(Tk, T, y, P)
    humrat(MoistAir, uconvert(u"K", Tk).val, MolarFrac, xv, uconvert(u"Pa", P).val)
end

function molarfrac(::Type{MoistAir}, Tk::Quantity,
                   ::Type{T}, y::Quantity, P::Quantity) where {T<:PsychroProperty}

    molarfrac(Tk, T, y, P)
end

function spechum(::Type{MoistAir}, Tk::Quantity,
                ::Type{T}, y::Quantity, P::Quantity) where {T<:PsychroProperty}

    xv = molarfrac(Tk, T, y, P)
    spechum(MoistAir, uconvert(u"K", Tk).val, MolarFrac, xv, uconvert(u"Pa", P).val)
end
