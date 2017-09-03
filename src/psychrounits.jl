#=
Implement the functions using the Unitful package to check units
=#

using Unitful

function volume(::Type{DryAir}, Tk::Q where {Q<:Quantity}, P::Q where {Q<:Quantity})
    volume(DryAir, uconvert(u"K", Tk).val, uconvert(u"Pa", P).val)*1.0u"m^3/kg"
end

function molarvolume(::Type{DryAir}, Tk::Q where {Q<:Quantity}, P::Q where {Q<:Quantity})
    molarvolume(DryAir, uconvert(u"K", Tk).val, uconvert(u"Pa", P).val)*1.0u"m^3/mol"
end

function density(::Type{DryAir}, Tk::Q where {Q<:Quantity}, P::Q where {Q<:Quantity})
    density(DryAir, uconvert(u"K", Tk).val, uconvert(u"Pa", P).val)*1.0u"kg/m^3"
end

function enthalpy(::Type{DryAir}, Tk::Q where {Q<:Quantity}, P::Q where {Q<:Quantity})
    enthalpy(DryAir, uconvert(u"K", Tk).val, uconvert(u"Pa", P).val)*1.0u"J/kg"
end

function molarenthalpy(::Type{DryAir}, Tk::Q where {Q<:Quantity}, P::Q where {Q<:Quantity})
    molarenthalpy(DryAir, uconvert(u"K", Tk).val, uconvert(u"Pa", P).val)*1.0u"J/mol"
end

function entropy(::Type{DryAir}, Tk::Q where {Q<:Quantity}, P::Q where {Q<:Quantity})
    entropy(DryAir, uconvert(u"K", Tk).val, uconvert(u"Pa", P).val)*1.0u"J/kg/K"
end

function molarentropy(::Type{DryAir}, Tk::Q where {Q<:Quantity}, P::Q where {Q<:Quantity})
    molarentropy(DryAir, uconvert(u"K", Tk).val, uconvert(u"Pa", P).val)*1.0u"J/mol/K"
end

function compressfactor(::Type{DryAir}, Tk::Q where {Q<:Quantity}, P::Q where {Q<:Quantity})
    compressfactor(DryAir, uconvert(u"K", Tk).val, uconvert(u"Pa", P).val)
end


function volume(::Type{Vapor}, Tk::Q where {Q<:Quantity}, P::Q where {Q<:Quantity})
    volume(Vapor, uconvert(u"K", Tk).val)*1.0u"m^3/kg"
end

function molarvolume(::Type{Vapor}, Tk::Q where {Q<:Quantity}, P::Q where {Q<:Quantity})
    molarvolume(Vapor, uconvert(u"K", Tk).val)*1.0u"m^3/mol"
end

function density(::Type{Vapor}, Tk::Q where {Q<:Quantity}, P::Q where {Q<:Quantity})
    density(Vapor, uconvert(u"K", Tk).val)*1.0u"kg/m^3"
end

function enthalpy(::Type{Vapor}, Tk::Q where {Q<:Quantity}, P::Q where {Q<:Quantity})
    Tk1 = uconvert(u"K", Tk).val
    enthalpy(Vapor, Tk1.val)*1.0u"J/kg"
end

function molarenthalpy(::Type{Vapor}, Tk::Q where {Q<:Quantity}, P::Q where {Q<:Quantity})
    molarenthalpy(Vapor, uconvert(u"K", Tk).val)*1.0u"J/mol"
end


function entropy(::Type{Vapor}, Tk::Q where {Q<:Quantity}, P::Q where {Q<:Quantity})
    entropy(Vapor, uconvert(u"K", Tk).val)*1.0u"J/kg/K"
end

function molarentropy(::Type{Vapor}, Tk::Q where {Q<:Quantity}, P::Q where {Q<:Quantity})
    molarentropy(Vapor, uconvert(u"K", Tk).val)*1.0u"J/mol/K"
end

function compressfactor(::Type{Vapor}, Tk::Q where {Q<:Quantity}, P::Q where {Q<:Quantity})
    compressfactor(Vapor, uconvert(u"K", Tk).val)
end

