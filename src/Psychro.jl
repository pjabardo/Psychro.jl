

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
 * [5] ASHRAE, "Psychrometrics: Theory and Practice", ASHRAE, 1996.

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

"Abstract type that represents any kind of matter"
abstract type AbstractMaterial end

"Abstract type that represents any kind of fluid"
abstract type AbstractFluid <: AbstractMaterial end


"""
Type used to model dry air

This uses the formulations presente in the ASHRAE handbook Psychrometrics: Theory and Practice, reference [5].

The methods dealing with this type implement the equations of reference [2].

Some thermodynamic properties are implemented as methods. Every method takes as arguments, the type Vapor` as first argument and the temperature as the second. If `Unitful` package is used, any unit can be used. On the other hand if no units are provided, the packages assumes SI units all through. In particular:

 * Temperature in K
 * Pressure in Pa
 * Length in m
 * Quantity of matter in mol
 * Mass in kg
 * Energy in J


The following methods are implemented:

 * `volume(DryAir, T, P)` specific volume
 * `volumem(DryAir, T, P` molar volume
 * `density(DryAir, T, P)` density
 * `enthalpy(DryAir, T, P)` specific enthalpy
 * `enthalpym(DryAir, T, P)` molar enthalpy
 * `entropy(DryAir, T, P)` specific entropy
 * `entropym(DryAir, T, P)` molar entropy
 * `compressfactor(DryAir, T, P)` Compressibility factor.
 
"""
struct DryAir <: AbstractFluid end


"""
Type used to model saturated water vapor

This uses the formulations presente in the ASHRAE handbook Psychrometrics: Theory and Practice, reference [5].

The methods dealing with this type implement the equations of reference [1].

Some thermodynamic properties are implemented as methods. Every method takes as arguments, the type `DryAir` as first argument and the saturation temperature. If `Unitful` package is used, any unit can be used. On the other hand if no units are provided, the packages assumes SI units all through. In particular:

 * Temperature in K
 * Pressure in Pa
 * Length in m
 * Quantity of matter in mol
 * Mass in kg
 * Energy in J


The following methods are implemented:

 * `volume(Vapor, T)` specific volume
 * `volumem(Vapor, T)` molar volume
 * `density(Vapor, T)` density
 * `enthalpy(Vapor, T)` specific enthalpy
 * `enthalpym(Vapor, T)` molar enthalpy
 * `entropy(Vapor, T)` specific entropy
 * `entropym(Vapor, T)` molar entropy
 * `compressfactor(Vapor, T)` Compressibility factor.
 
"""
struct Vapor <: AbstractFluid end


"""
Type used to model moist air

This uses the formulations presente in the ASHRAE handbook Psychrometrics: Theory and Practice, reference [5]. 

The temperature range is `173.15 < T < 474.15 K` with pressures `P < 5 MPa`.

The methods dealing with this type implement the equations of references [1-2] to calculate some thermodynamics of real moits air.

Some thermodynamic properties are implemented as methods. Every method takes as arguments, the type `MoistAir` as first argument, temperature as second argument, a type declaring how humidity is specified, the value of humidity using and pressure. If `Unitful` package is used, any unit can be used. On the other hand if no units are provided, the packages assumes SI units all through. In particular:

 * Temperature in K
 * Pressure in Pa
 * Length in m
 * Quantity of matter in mol
 * Mass in kg
 * Energy in J

To specify humidity the following types are implemented:

 * `WetBulb` for wet bulb temperature, actually adiabatic saturation temperature
 * `DewPoint` Dew point temperature
 * `RelHum` Relative humidity
 * `HumRat` Humidity ratio (kg of vapor / kg of dry air)
 * `SpecHum` Specific humidity (kg of vapor / kg of moist air)
 * `MolarFrac` Molar fraction of water vapor.

The following methods are implemented:

 * `volume(MoistAir, T, HumType, h, P)` specific volume
 * `volumem(MoistAir, T, HumType, h, P` molar volume
 * `density(MoistAir, T, HumType, h, P)` density
 * `enthalpy(MoistAir, T, HumType, h, P)` specific enthalpy
 * `enthalpym(MoistAir, T, HumType, h, P)` molar enthalpy
 * `entropy(MoistAir, T, HumType, h, P)` specific entropy
 * `entropym(MoistAir, T, HumType, h, P)` molar entropy
 * `compressfactor(MoistAir, T, HumType, h, P)` Compressibility factor.

It should be noted that the specific enthalpy, entropy and volume are calculated with respect to the mass of dry air. Therefore, 1 J/kg is actually 1 J per kg of dry air.
"""
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
