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
#include "material_manager.h"

#include "materials/stiffened_gas.h"
#include "materials/stiffened_gas_safe.h"
#include "materials/stiffened_gas_complete_safe.h"
#include "materials/waterlike_fluid.h"
#include "materials/noble_abel_stiffened_gas.h"
#include <stdexcept>

/**
 * @brief Sets up a MaterialManager handling all requests to all materials provided here as input.
 *        Creates the material objects directly in place.
 * @param material_data The raw data of the materials to be created.
 * @param surface_tension_coefficients The surface tension coefficients between two materials. %Currently only implemented for two-material simulations%.
 */
MaterialManager::MaterialManager( std::vector<std::tuple<MaterialName, MaterialName, std::unordered_map<std::string, double>>> const material_data,
   std::vector<double> surface_tension_coefficients ) :
   surface_tension_coefficients_(surface_tension_coefficients) {
#ifdef PERFORMANCE
   material_names_.reserve(material_data.size());
  equations_of_state_.reserve(material_data.size());
#endif
   for(const auto& material : material_data) {
      AddMaterial(material);
   }
}

/**
 * @brief Gives an instance of the material with the given identifier.
 * @param material .
 * @return The Material object (equation of state).
 */
Material const& MaterialManager::GetEquationOfState( MaterialName const material ) const {
#ifndef PERFORMANCE
   return *materials_.at(material);
#else
   std::size_t index = 0;
   for(std::size_t i = 0; i< material_names_.size(); ++i) {
      if(material_names_[i] == material) {
         index = i;
      }
   }
   return *equations_of_state_[index];
#endif
}

/**
 * @brief Proxy function to obtain the pressure for the provided material and inputs.
 * @param material .
 * @param density .
 * @param momentum_x .
 * @param momentum_y .
 * @param momentum_z .
 * @param energy .
 * @return pressure.
 */
double MaterialManager::GetPressure( MaterialName const material, double const density, double const momentum_x, double const momentum_y, double const momentum_z, double const energy ) const {
   Material const& equation_of_state = GetEquationOfState(material);
   return equation_of_state.GetPressure(density, momentum_x, momentum_y, momentum_z, energy);
}

/**
 * @brief Proxy function to obtain the enthalpy for the provided material and inputs.
 * @param material .
 * @param density .
 * @param momentum_x .
 * @param momentum_y .
 * @param momentum_z .
 * @param energy .
 * @return enthalpy.
 */
double MaterialManager::GetEnthalpy( MaterialName const material, double const density, double const momentum_x, double const momentum_y, double const momentum_z, double const energy ) const {
   Material const& equation_of_state = GetEquationOfState(material);
   return equation_of_state.GetEnthalpy(density, momentum_x, momentum_y, momentum_z, energy);
}

/**
 * @brief Proxy function to obtain the energy for the provided material and inputs.
 * @param material .
 * @param density .
 * @param momentum_x .
 * @param momentum_y .
 * @param momentum_z .
 * @param pressure .
 * @return energy .
 */
double MaterialManager::GetEnergy( MaterialName const material, double const density, double const momentum_x, double const momentum_y, double const momentum_z, double const pressure ) const {
   Material const& equation_of_state = GetEquationOfState(material);
   return equation_of_state.GetEnergy(density, momentum_x, momentum_y, momentum_z, pressure);
}

/**
 * @brief Proxy function to obtain the temperature for the provided material and inputs.
 * @param material .
 * @param density .
 * @param momentum_x .
 * @param momentum_y .
 * @param momentum_z .
 * @param energy .
 * @return temperature.
 */
double MaterialManager::GetTemperature( MaterialName const material, double const density, double const momentum_x, double const momentum_y, double const momentum_z, double const energy ) const {
   Material const& equation_of_state = GetEquationOfState(material);
   return equation_of_state.GetTemperature(density, momentum_x, momentum_y, momentum_z, energy);
}

/**
 * @brief Proxy function to obtain the Grueneisen coefficient for a given material according to the generalized equation of state.
 * @param material .
 * @return Grueneisen coefficient .
 */
double MaterialManager::GetGruneisen( MaterialName const material ) const {
#ifndef PERFORMANCE
   if constexpr( CC::GruneisenDensityDependent() ) {
      throw std::runtime_error( "CC::GruneisenDensityDependent is true, you should call MaterialManager::GetGruneisen(material, density) instead" );
   }
#endif
   Material const& equation_of_state = GetEquationOfState(material);
   return equation_of_state.GetGruneisen();
}

/**
 * @brief Proxy function to obtain the Grueneisen coefficient for a given material according to the generalized equation of state.
 * @param material .
 * @param density .
 * @return Grueneisen coefficient .
 */
double MaterialManager::GetGruneisen( MaterialName const material, double const density ) const {
   Material const& equation_of_state = GetEquationOfState( material );
   return equation_of_state.GetGruneisen( density );
}

/**
 * @brief Proxy function to obtain the specific material parameter psi for a given material according to the generalized equation of state.
 * @param material .
 * @param pressure .
 * @param one over density (saves the expensive division operation) .
 * @return Psi for given inputs.
 */
double MaterialManager::GetPsi( MaterialName const material, double const pressure, double const one_density ) const {
   Material const& equation_of_state = GetEquationOfState(material);
   return equation_of_state.GetPsi(pressure, one_density);
}

/**
 * @brief Proxy function to obtain the Grueneisen coefficient for a given material according to the generalized equation of state.
 * @param material .
 * @return Grueneisen coefficient .
 */
double MaterialManager::GetGamma( MaterialName const material ) const {
   Material const& equation_of_state = GetEquationOfState(material);
   return equation_of_state.GetGamma();
}

/**
 * @brief Proxy function to obtain B for a given material.
 * @param material .
 * @return B .
 */
double MaterialManager::GetB( MaterialName const material ) const {
   Material const& equation_of_state = GetEquationOfState(material);
   return equation_of_state.GetB();
}

/**
 * @brief Proxy function to obtain the speed of sound for a given material, density and pressure.
 * @param material .
 * @param density .
 * @param pressure .
 * @return speed of sound for given inputs.
 */
double MaterialManager::GetSpeedOfSound( MaterialName const material, double const density, double const pressure ) const{
   Material const& equation_of_state = GetEquationOfState(material);
   return equation_of_state.GetSpeedOfSound(density, pressure);
}

/**
 * @brief Proxy function to obtain the viscosity of the specified material.
 * @param material .
 * @return Viscosity parameters (multiple parameters, e.g. shear and bulk) according to material. %Constant per Material, may change in future versions%.
 */
std::vector<double> MaterialManager::GetViscosity( MaterialName const material ) const {
   Material const& equation_of_state = GetEquationOfState(material);
   return equation_of_state.GetViscosity();
}

/**
 * @brief Proxy function to obtain the surface tension coefficients between two materials.
 * @param first_material, second_material Unique identifiers of the materials of interest. Should give the same result independent of the order.
 * @return Surface tension coefficient between two materials.
 */
double MaterialManager::GetSurfaceTensionCoefficient( MaterialName const, MaterialName const ) const {
   unsigned int index = 0;
   return surface_tension_coefficients_[index];
}

/**
 * @brief Proxy function to obtain the thermal conductivity of the specified material.
 * @param material .
 * @return Thermal conductivity %Constant per Material, may change in future versions%.
 */
double MaterialManager::GetThermalConductivity( MaterialName const material ) const {
   Material const& equation_of_state = GetEquationOfState(material);
   return equation_of_state.GetThermalConductivity();
}

/**
 * @brief Proxy function to obtain the specific heat of the specified material.
 * @param material .
 * @return Specific heat %Constant per Material, may change in future versions%.
 */
double MaterialManager::GetSpecificHeat( MaterialName const material ) const {
   Material const& equation_of_state = GetEquationOfState(material);
   return equation_of_state.GetSpecificHeat();
}

/**
 * @brief Creates in place and registers a further material to this instance of the MaterialManager.
 * @param data Material data consisting of two material identifiers. The first naming the generic type e.g. "StiffenedGas",
 *        the second the unique identifier e.g. "StiffenedGasOne".
 *        The tuple also holds a variable (depending on material type) number of parameters to be set for the material to be created.
 */
void MaterialManager::AddMaterial( std::tuple<MaterialName, MaterialName, std::unordered_map<std::string, double>> const data ) {

   auto const& map = std::get<2>(data);
   switch (std::get<0>(data)) {
      case MaterialName::StiffenedGas :
#ifndef PERFORMANCE
         materials_.emplace(std::get<1>(data),
            std::make_unique<StiffenedGas const>( map.at("gamma"), map.at("B"), map.at("dynamicShear"), map.at("dynamicBulk"))
         );
#else
         material_names_.push_back(std::get<1>(data));
         equations_of_state_.push_back(
            std::make_unique<StiffenedGas const >( map.at("gamma"), map.at("B"), map.at("dynamicShear"), map.at("dynamicBulk"))
         );
#endif
         break;
      case MaterialName::StiffenedGasSafe :
#ifndef PERFORMANCE
         materials_.emplace(std::get<1>(data),
            std::make_unique<StiffenedGasSafe const>( map.at("gamma"), map.at("B"), map.at("dynamicShear"), map.at("dynamicBulk"))
         );
#else
         material_names_.push_back(std::get<1>(data));
         equations_of_state_.push_back(
            std::make_unique<StiffenedGasSafe const>( map.at("gamma"), map.at("B"), map.at("dynamicShear"), map.at("dynamicBulk"))
         );
#endif
         break;
      case MaterialName::StiffenedGasCompleteSafe :
#ifndef PERFORMANCE
         materials_.emplace(std::get<1>(data),
            std::make_unique<StiffenedGasCompleteSafe const>( map.at("gamma"), map.at("A"), map.at("B"), map.at("C"),
                                                              map.at("specificGasConstant"), map.at("thermalConductivity"),
                                                              map.at("dynamicShear"), map.at("dynamicBulk"))
         );
#else
         material_names_.push_back(std::get<1>(data));
         equations_of_state_.push_back(
            std::make_unique<StiffenedGasCompleteSafe const>( map.at("gamma"), map.at("A"), map.at("B"), map.at("C"),
                                                              map.at("specificGasConstant"), map.at("thermalConductivity"),
                                                              map.at("dynamicShear"), map.at("dynamicBulk"))
         );
#endif
         break;
      case MaterialName::WaterlikeFluid :
#ifndef PERFORMANCE
         materials_.emplace(std::get<1>(data),
            std::make_unique<WaterlikeFluid const>( map.at("gamma"), map.at("A"), map.at("B"), map.at("rho0"), map.at("dynamicShear"), map.at("dynamicBulk"))
         );
#else
         material_names_.push_back(std::get<1>(data));
         equations_of_state_.push_back(
            std::make_unique<WaterlikeFluid const>( map.at("gamma"), map.at("A"), map.at("B"), map.at("rho0"), map.at("dynamicShear"), map.at("dynamicBulk"))
         );
#endif
         break;
      case MaterialName::NobleAbelStiffenedGas :
#ifndef PERFORMANCE
         materials_.emplace(std::get<1>(data),
            std::make_unique<NobleAbelStiffenedGas const>( map.at("gamma"), map.at("covolume"), map.at("pressureConstant"), map.at("energyConstant"),
                                                           map.at("specificHeatCapacity"), map.at("thermalConductivity"),
                                                           map.at("dynamicShear"), map.at("dynamicBulk"))
         );
#else
         material_names_.push_back(std::get<1>(data));
         equations_of_state_.push_back(
            std::make_unique<NobleAbelStiffenedGas const>( map.at("gamma"), map.at("covolume"), map.at("pressureConstant"), map.at("energyConstant"),
                                                           map.at("specificHeatCapacity"), map.at("thermalConductivity"),
                                                           map.at("dynamicShear"), map.at("dynamicBulk"))
         );
#endif
         break;
      default:
         throw std::logic_error("This material has not yet been implemented");
         break;
   }
}
