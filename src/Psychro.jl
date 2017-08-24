

"""
# Thermodynamic properties of moist air

This  module implements functions to calculate the
thermodynamic properties of moist air. It uses
correlations recommended by ASHRAE and described
in the articles listed below.

## References

 * [1] Wexler, A. and Hyland, R. W., "Formulations for the thermodynamic properties of the saturated phases of H2O from 173.15 K to 473.15 K", ASHRAE Transactions, 1983.
 * [2] Wexler, A. and Hyland, R. W., "Formulations for the thermodynamic properties of dry air from 173.15 K to 473.15 K, and of saturated moist air from 173.15 K to 372.15 K at pressures to 5 MPa



"""
module Psychro

# package code goes here
include("ashrae.jl")
end # module
