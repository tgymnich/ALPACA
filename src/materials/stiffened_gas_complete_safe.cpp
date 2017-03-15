/*****************************************************************************************
*                                                                                        *
* This file is part of ALPACA                                                            *
*                                                                                        *
******************************************************************************************
*                                                                                        *
*  \\                                                                                    *
*  l '>                                                                                  *
*  | |                                                                                   *
*  | |                                                                                   *
*  | alpaca~                                                                             *
*  ||    ||                                                                              *
*  ''    ''                                                                              *
*                                                                                        *
* ALPACA is a MPI-parallelized C++ code framework to simulate compressible multiphase    *
* flow physics. It allows for advanced high-resolution sharp-interface modeling          *
* empowered with efficient multiresolution compression. The modular code structure       *
* offers a broad flexibility to select among many most-recent numerical methods covering *
* WENO/T-ENO, Riemann solvers (complete/incomplete), strong-stability preserving Runge-  *
* Kutta time integration schemes, level set methods and many more.                       *
*                                                                                        *
* This code is developed by the 'Nanoshock group' at the Chair of Aerodynamics and       *
* Fluid Mechanics, Technical University of Munich.                                       *
*                                                                                        *
******************************************************************************************
*                                                                                        *
* LICENSE                                                                                *
*                                                                                        *
* ALPACA - Adaptive Level-set PArallel Code Alpaca                                       *
* Copyright (C) 2020 Nikolaus A. Adams and contributors (see AUTHORS list)               *
*                                                                                        *
* This program is free software: you can redistribute it and/or modify it under          *
* the terms of the GNU General Public License as published by the Free Software          *
* Foundation version 3.                                                                  *
*                                                                                        *
* This program is distributed in the hope that it will be useful, but WITHOUT ANY        *
* WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A        *
* PARTICULAR PURPOSE. See the GNU General Public License for more details.               *
*                                                                                        *
* You should have received a copy of the GNU General Public License along with           *
* this program (gpl-3.0.txt).  If not, see <https://www.gnu.org/licenses/gpl-3.0.html>   *
*                                                                                        *
******************************************************************************************
*                                                                                        *
* THIRD-PARTY tools                                                                      *
*                                                                                        *
* Please note, several third-party tools are used by ALPACA. These tools are not shipped *
* with ALPACA but available as git submodule (directing to their own repositories).      *
* All used third-party tools are released under open-source licences, see their own      *
* license agreement in 3rdParty/ for further details.                                    *
*                                                                                        *
* 1. tiny_xml           : See LICENSE_TINY_XML.txt for more information.                 *
* 2. expression_toolkit : See LICENSE_EXPRESSION_TOOLKIT.txt for more information.       *
* 3. FakeIt             : See LICENSE_FAKEIT.txt for more information                    *
* 4. Catch2             : See LICENSE_CATCH2.txt for more information                    *
*                                                                                        *
******************************************************************************************
*                                                                                        *
* CONTACT                                                                                *
*                                                                                        *
* nanoshock@aer.mw.tum.de                                                                *
*                                                                                        *
******************************************************************************************
*                                                                                        *
* Munich, July 1st, 2020                                                                 *
*                                                                                        *
*****************************************************************************************/
#include "stiffened_gas_complete_safe.h"
#include "mathematical_functions.h"

/**
 * @brief Constructs a stiffened-gas complete object with material parameters given as input.
 * @param gamma .
 * @param A .
 * @param B .
 * @param C .
 * @param rho0 .
 * @param specific_gas_constant .
 * @param mu_shear Shear viscosity.
 * @param mu_bulk Bulk viscosity.
 */
StiffenedGasCompleteSafe::StiffenedGasCompleteSafe( double const gamma, double const A, double const B, double const C, double const specific_gas_constant, double const thermal_conductivity, double const mu_shear, double const mu_bulk ) :
   gamma_(gamma),
   energy_translation_factor_(A),
   stiffened_pressure_constant_(B),
   thermal_energy_factor_(C),
   specific_gas_constant_(specific_gas_constant),
   thermal_conductivity_(thermal_conductivity),
   mu_shear_(mu_shear),
   mu_bulk_(mu_bulk),
   specific_heat_( 1.0 / (gamma_ - 1.0) * specific_gas_constant_ ) {
   // Empty besides initializer list.
}

/**
 * @brief Computes the pressure from inputs as -gamma * B + (gamma - 1) * (E  - 0.5 * ||v^2||) - rho * A (avoids devision by zero).
 * @param density .
 * @param momentum_x .
 * @param momentum_y .
 * @param momentum_z .
 * @param energy .
 * @return Pressure according to complete stiffened-gas equation of state.
 */
double StiffenedGasCompleteSafe::DoGetPressure( double const density, double const momentum_x, double const momentum_y, double const momentum_z, double const energy ) const {
   double pressure = -gamma_ * stiffened_pressure_constant_ + (gamma_ -1.0) * (energy - 0.5 * DimensionAwareConsistencyManagedSum( momentum_x * momentum_x,
      momentum_y * momentum_y,
      momentum_z * momentum_z) / std::max( density, epsilon_ ) - density * energy_translation_factor_);
   return std::max( pressure ,-stiffened_pressure_constant_ + epsilon_ );
}

/**
 * @brief Computes enthalphy as (E + p) / rho (avoids devision by zero).
 * @param density .
 * @param momentum_x .
 * @param momentum_y .
 * @param momentum_z .
 * @param energy .
 * @return Enthalpy value.
 */
double StiffenedGasCompleteSafe::DoGetEnthalpy( double const density, double const momentum_x, double const momentum_y, double const momentum_z, double const energy ) const {
   return (energy + DoGetPressure(density, momentum_x, momentum_y, momentum_z, energy)) / std::max( density, epsilon_ );
}

/**
 * @brief Computes energy from inputs as (p + gamma * B) / (gamma - 1) + 0.5 * ||v^2|| (avoids devision by zero).
 * @param density .
 * @param momentum_x .
 * @param momentum_y .
 * @param momentum_z .
 * @param pressure .
 * @return Energy according to complete stiffened-gas equation of state.
 */
double StiffenedGasCompleteSafe::DoGetEnergy( double const density, double const momentum_x, double const momentum_y, double const momentum_z, double const pressure ) const {
   return density * energy_translation_factor_ + ( pressure + gamma_ * stiffened_pressure_constant_ ) / ( gamma_ -1.0 )
          + ( 0.5 * DimensionAwareConsistencyManagedSum( momentum_x * momentum_x, momentum_y * momentum_y, momentum_z * momentum_z ) / std::max( density, epsilon_ ) );
}

/**
 * @brief Computes temperature for stiffened-gas EOS (avoids devision by zero).
 * @param density .
 * @param momentum_x .
 * @param momentum_y .
 * @param momentum_z .
 * @param energy .
 * @return Temperature according to complete stiffened-gas EOS.
 */
double StiffenedGasCompleteSafe::DoGetTemperature( double const density, double const momentum_x, double const momentum_y, double const momentum_z, double const energy ) const {
   return (DoGetPressure(density, momentum_x, momentum_y, momentum_z, energy) + stiffened_pressure_constant_) /
      (std::max( density, epsilon_ ) * specific_gas_constant_) + ((gamma_ - 1.0) * thermal_energy_factor_ * std::pow( density, gamma_ - 1.0 )) / specific_gas_constant_;
}

/**
 * @brief Computes Gruneisen coefficent as (gamma-1) for stiffened-gas equation of state.
 * @return Gruneisen coefficent according to complete stiffened-gas EOS.
 */
double StiffenedGasCompleteSafe::DoGetGruneisen() const {
   return (gamma_ - 1.0);
}

/**
 * @brief Computes psi from inputs as (p + gamma * B) / rho.
 * @param pressure .
 * @param one_density (devision by zero is already avoided before) .
 * @return Psi according to complete stiffened-gas equation of state.
 */
double StiffenedGasCompleteSafe::DoGetPsi( double const pressure, double const one_density ) const {
   return (pressure + gamma_ * stiffened_pressure_constant_ ) * one_density;
}

/**
 * @brief Returns gamma.
 * @return gamma.
 */
double StiffenedGasCompleteSafe::DoGetGamma() const {
   return gamma_;
}

/**
 * @brief Returns B.
 * @return B.
 */
double StiffenedGasCompleteSafe::DoGetB() const {
   return stiffened_pressure_constant_;
}

/**
 * @brief Computes speed of sound from inputs as sqrt(gamma * (p + B) / rho).
 * @param density .
 * @param pressure .
 * @return Speed of sound according to complete stiffened-gas equation of state.
 */
double StiffenedGasCompleteSafe::DoGetSpeedOfSound( double const density, double const pressure ) const {
   double speed_of_sound_squared = gamma_ * (pressure + stiffened_pressure_constant_) / std::max( density, epsilon_ );
   return std::sqrt( std::max( speed_of_sound_squared, epsilon_ ) );
}

/**
 * @brief Gives the viscosity of the material.
 * @return first: shear viscosity; second: bulk viscosity
 */
std::vector<double> StiffenedGasCompleteSafe::DoGetViscosity() const {
   return {mu_shear_, mu_bulk_};
}

/**
 * @brief Gives the thermal conductivity of the material.
 * @return Thermal conductivity
 */
double StiffenedGasCompleteSafe::DoGetThermalConductivity() const {
   return thermal_conductivity_;
}

/**
 * @brief Returns the isochoric specific heat of the material.
 * @return Specific heat.
 */
double StiffenedGasCompleteSafe::DoGetSpecificHeat() const {
   return specific_heat_;
}
