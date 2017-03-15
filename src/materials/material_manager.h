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
#ifndef MATERIAL_MANAGER_H
#define MATERIAL_MANAGER_H

#include <memory>
#include <tuple>
#include <unordered_map>

#include "materials/material.h"
#include "materials/material_names.h"
#include "mathematical_functions.h"

/**
 * @brief The MaterialManager class provides access to all materials present in the current simulation and forwards the appropriate object to the caller.
 *        It thus acts as a proxy to obtain non-conservative values.
 */
class MaterialManager {
#ifndef PERFORMANCE
   // Sadly, the compiler is not willing to optimize the Hash Map as hard as raw loops (Because of possible exceptions).
   std::unordered_map<MaterialName, std::unique_ptr<const Material>> materials_;
#else
   // So we create our own (very poor) "map" with Black Jack and H...
  std::vector<MaterialName> material_names_;
  std::vector<std::unique_ptr<const Material>> equations_of_state_;
#endif

   std::vector<double> surface_tension_coefficients_;

   void AddMaterial( std::tuple<MaterialName, MaterialName, std::unordered_map<std::string, double>> const data );
   Material const& GetEquationOfState( MaterialName const material ) const;

public:

   MaterialManager() = delete;
   explicit MaterialManager( std::vector<std::tuple<MaterialName, MaterialName, std::unordered_map<std::string, double>>> const material_data, std::vector<double> surface_tension_coefficients );
   ~MaterialManager() = default;
   MaterialManager( MaterialManager const& ) = delete;
   MaterialManager& operator=( MaterialManager const& ) = delete;
   MaterialManager( MaterialManager&& ) = delete;
   MaterialManager& operator=( MaterialManager&& ) = delete;

   double GetPressure   ( MaterialName const material, double const density, double const momentum_x, double const momentum_y, double const momentum_z, double const energy ) const;
   double GetEnthalpy   ( MaterialName const material, double const density, double const momentum_x, double const momentum_y, double const momentum_z, double const energy ) const;
   double GetEnergy     ( MaterialName const material, double const density, double const momentum_x, double const momentum_y, double const momentum_z, double const pressure ) const;
   double GetTemperature( MaterialName const material, double const density, double const momentum_x, double const momentum_y, double const momentum_z, double const energy ) const;

   double GetGruneisen( MaterialName const material ) const;
   double GetGruneisen( MaterialName const material, double const density ) const;
   double GetPsi      ( MaterialName const material, double const pressure, double const one_density ) const;
   double GetSpeedOfSound( MaterialName const material, double const density, double const pressure ) const;

   double GetGamma( MaterialName const material ) const;
   double GetB( MaterialName const material ) const;

   std::vector<double> GetViscosity( MaterialName const material ) const;
   double GetSurfaceTensionCoefficient( MaterialName const first_material, MaterialName const second_material ) const;

   double GetThermalConductivity( MaterialName const material ) const;
   double GetSpecificHeat( MaterialName const material ) const;
};

#endif // MATERIAL_MANAGER_H
