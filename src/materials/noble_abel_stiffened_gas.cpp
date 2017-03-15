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
#include "noble_abel_stiffened_gas.h"
#include "mathematical_functions.h"
#include <stdexcept>

/**
 * @brief Constructs a Noble-Abel stiffened-gas object with material parameters given as input.
 * @param mu_shear Shear viscosity.
 * @param mu_bulk Bulk viscosity.
 */
NobleAbelStiffenedGas::NobleAbelStiffenedGas( double const gamma, double const covolume, double const pressure_constant, double const energy_constant,
      double const specific_heat_capacity, double const thermal_conductivity, double const mu_shear, double const mu_bulk ) :
   gamma_(gamma),
   covolume_(covolume),
   pressure_constant_(pressure_constant),
   energy_constant_(energy_constant),
   specific_heat_capacity_(specific_heat_capacity),
   thermal_conductivity_(thermal_conductivity),
   mu_shear_(mu_shear),
   mu_bulk_(mu_bulk) {
   if constexpr( !CC::GruneisenDensityDependent() ) {
      throw std::runtime_error( "To use NobleAbelStiffenedGas you need to activate CC::GruneisenDensityDependent" );
   }
}

/**
 * @brief Computes the pressure from inputs as (gamma - 1) * ((E  - 0.5 * rho * ||v^2||) / rho - e_const) / (1 / rho - covolume) - gamma * p_const.
 * @param density .
 * @param momentum_x .
 * @param momentum_y .
 * @param momentum_z .
 * @param energy .
 * @return Pressure according to NASG equation of state.
 */
double NobleAbelStiffenedGas::DoGetPressure( double const density, double const momentum_x, double const momentum_y, double const momentum_z, double const energy ) const {
   double const one_density = 1.0 / density;
   double const internal_energy = one_density * (energy - one_density * 0.5 * DimensionAwareConsistencyManagedSum( momentum_x * momentum_x, momentum_y * momentum_y, momentum_z * momentum_z ));

   return (gamma_ - 1.0) * (internal_energy - energy_constant_) / (one_density - covolume_) - gamma_ * pressure_constant_;
}

/**
 * @brief Computes enthalphy as (E + p) / rho.
 * @param density .
 * @param momentum_x .
 * @param momentum_y .
 * @param momentum_z .
 * @param energy .
 * @return Enthalpy value.
 */
double NobleAbelStiffenedGas::DoGetEnthalpy( double const density, double const momentum_x, double const momentum_y, double const momentum_z, double const energy ) const {
   return (energy + DoGetPressure( density, momentum_x, momentum_y, momentum_z, energy )) / density;
}

/**
 * @brief Computes energy from inputs as rho * ((p + gamma * p_const) / (gamma - 1) * (1 / rho - covolume) + e_const) + 0.5 * rho * ||v^2||
 * @param density .
 * @param momentum_x .
 * @param momentum_y .
 * @param momentum_z .
 * @param pressure .
 * @return Energy according to NASG equation of state.
 */
double NobleAbelStiffenedGas::DoGetEnergy( double const density, double const momentum_x, double const momentum_y, double const momentum_z, double const pressure ) const {
   double const internal_energy = (pressure + gamma_ * pressure_constant_) / (density * (gamma_ - 1.0)) * (1.0 - covolume_ * density) + energy_constant_;
   return density * internal_energy + 0.5 * DimensionAwareConsistencyManagedSum( momentum_x * momentum_x, momentum_y * momentum_y, momentum_z * momentum_z ) / density;
}

/**
 * @brief Computes temperature for Noble-Abel stiffened-gas EOS.
 * @param density .
 * @param momentum_x .
 * @param momentum_y .
 * @param momentum_z .
 * @param energy .
 * @return Temperature according to NASG EOS.
 */
double NobleAbelStiffenedGas::DoGetTemperature( double const density, double const momentum_x, double const momentum_y, double const momentum_z, double const energy ) const {
   double const pressure = DoGetPressure( density, momentum_x, momentum_y, momentum_z, energy );
   return (pressure + pressure_constant_) / (density * (gamma_ - 1.0) * specific_heat_capacity_) * (1.0 - covolume_ * density);
}

/**
 * @brief Computes Gruneisen coefficent as (gamma - 1) / (1 - covolume * rho)
 * @param density .
 * @return Gruneisen coefficent according to NASG equation of state.
 */
double NobleAbelStiffenedGas::DoGetGruneisen( double const density ) const {
   return (gamma_ - 1.0) / (1.0 - covolume_ * density);
}

/**
 * @brief DO NOT CALL. Throws an error since Gruneisen coefficient for NASG is density dependent. Call NobleAbelStiffenedGas::DoGetGruneisen(double const density) instead.
 */
double NobleAbelStiffenedGas::DoGetGruneisen() const {
   throw std::runtime_error( "NobleAbelStiffenedGas: Gruneisen parameter depends on density!" );
}

/**
 * @brief Computes psi from inputs as (p + gamma * p_const) / (rho * (1 - covolume * rho))
 * @param pressure .
 * @param one_density (devision by zero is already avoided before) .
 * @return Psi according to NASG equation of state.
 */
double NobleAbelStiffenedGas::DoGetPsi( double const pressure, double const one_density ) const {
   return (pressure + gamma_ * pressure_constant_ ) * one_density * one_density / (one_density - covolume_);
}

/**
 * @brief Returns gamma
 * @return gamma.
 */
double NobleAbelStiffenedGas::DoGetGamma() const {
   return gamma_;
}

/**
 * @brief Returns the pressure constant.
 * @return B.
 */
double NobleAbelStiffenedGas::DoGetB() const {
   return pressure_constant_;
}

/**
 * @brief Computes speed of sound from inputs as sqrt(gamma * (p + p_const) / (rho * (1 - covolume * rho)))
 * @param density .
 * @param pressure .
 * @return Speed of sound according to NASG equation of state.
 */
double NobleAbelStiffenedGas::DoGetSpeedOfSound( double const density, double const pressure ) const {
   double const speed_of_sound_squared = gamma_ * (pressure + pressure_constant_) / (density * (1.0 - covolume_ * density));
   return std::sqrt( speed_of_sound_squared );
}

/**
 * @brief Gives the viscosity of the material
 * @return First: shear viscosity; Second: bulk viscosity
 */
std::vector<double> NobleAbelStiffenedGas::DoGetViscosity() const {
   return {mu_shear_, mu_bulk_};
}

/**
 * @brief Gives the thermal conductivity of the material.
 * @return Thermal conductivity
 */
double NobleAbelStiffenedGas::DoGetThermalConductivity() const {
   return thermal_conductivity_;
}

/**
 * @brief Returns the isochoric specific heat of the material.
 * @return Specific heat.
 */
double NobleAbelStiffenedGas::DoGetSpecificHeat() const {
   return specific_heat_capacity_;
}
