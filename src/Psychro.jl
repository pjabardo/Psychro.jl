

"""
# Thermodynamic properties of moist air

This  module implements functions to calculate the
thermodynamic properties of moist air. It uses
correlations recommended by ASHRAE and described
in the articles listed below.

## References

 * [1] Wexler, A. and Hyland, R. W., "Formulations for the thermodynamic properties of the saturated phases of H2O from 173.15 K to 473.15 K", ASHRAE Transactions, 1983.
 * [2] Wexler, A. and Hyland, R. W., "Formulations for the thermodynamic properties of dry air from 173.15 K to 473.15 K, and of saturated moist air from 173.15 K to 372.15 K at pressures to 5 MPa
 * [3] Himmelblaum D. M., "Solubilities of inert gases in water, 0oC to near the critical point of water", Journal of Chemical and Engineering Data, Vol. 5, No. 1, January 1960.
 * [4] Kell, George S., "Density, thermal expansivity, and compressibility of liquid water from 0oC to 150oC: correlations and tables for atmospheric pressure and saturation reviewed and expressed on 1968 temperature scale", Journal of Chemical and Engineering Data, Vol. 20, No. 1, 1975.
"""

module Psychro

"Molecular weight of air in kg/mol"
const Ma = 0.0289635

"Molecular weight of water in kg/mol"
const Mv = 0.01801528

"Universal gas constant kg m^2 /s^2 /K /mol"
const R = 8.314459848

"Melting point of water"
const T0 = 273.15

abstract type AbstractMaterial end
abstract type AbstractFluid <: AbstractMaterial end

struct DryAir <: AbstractFluid end
struct Vapor <: AbstractFluid end
struct MoistAir <: AbstractFluid end


abstract type ThermodynamicProperty end
abstract type PsychroProperty end

struct Volume <: ThermodynamicProperty end
struct Enthalpy <: ThermodynamicProperty end
struct Entropy <: ThermodynamicProperty end
struct Temperature <: ThermodynamicProperty end
struct Pressure <: ThermodynamicProperty end
struct Density <: ThermodynamicProperty end

struct MolarFrac <: PsychroProperty end

"Specific Humidity mv / (ma+mv)"
struct SpecHum <: PsychroProperty end
struct WetBulb <: PsychroProperty end
struct DewPoint <: PsychroProperty end
struct HumRat <: PsychroProperty end
struct RelHum <: PsychroProperty end



# package code goes here
include("utilities.jl")
include("hyland83.jl")
include("hyland83a.jl")
include("wetbulb.jl")
include("psychrounits.jl")
include("moistair.jl")

export volume, volumem, density, enthalpy, enthalpym, entropy, entropym, compressfactor
export dewpoint, relhum, humrat, wetbulb, molarfrac, spechum
export DryAir, Vapor, MoistAir
export Volume, Enthalpy, Entropy, Temperature, Pressure, Density
export MolarFrac, SpecHum, WetBulb, DewPoint, HumRat, RelHum


end # module
