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
#include <catch.hpp>
#include "solvers/source_term_contributions/axisymmetric_viscous_volume_forces.h"
#include "topology/node.h"

namespace {
   /**
    * @brief Creates the material manager that is used during the testing.
    * @param material_name Name of the material to be added.
    * @param shear_viscosity Shear viscosity of the material.
    * @param bulk_viscosity Bulk viscosity of the material.
    * @return Initializes MaterialManager class.
    */
   MaterialManager ReturnMaterialManagerWithViscosities( MaterialName const material_name, double const shear_viscosity, double const bulk_viscosity ) {
      std::unordered_map<std::string, double> const parameter_map = { { "gamma", 1.4 }, { "A", 0.0 }, { "B", 0.0 }, { "C", 0.0 }, { "rho0", 0.0 }, { "specificGasConstant", 0.0 }, { "thermalConductivity", 0.0 },
         { "dynamicShear", shear_viscosity }, { "dynamicBulk", bulk_viscosity } };

      return MaterialManager( { std::make_tuple( material_name, material_name, parameter_map ) }, {} );
   }
}

SCENARIO( "Volume Forces Calculation Correctness", "[1rank]" ) {

   GIVEN( "An axisymmetric viscous volume forces calculator with cell-dependent x- and r-velocity component" ) {
      std::pair<MaterialName const, Block> mat_block( std::piecewise_construct, std::make_tuple( MaterialName::StiffenedGas ), std::make_tuple() );

      for( unsigned int i = 0; i < CC::TCX(); ++i ) {
         for( unsigned int j = 0; j < CC::TCY(); ++j ) {
            for( unsigned int k = 0; k < CC::TCZ(); ++k ) {
               mat_block.second.GetPrimeStateBuffer( PrimeState::VelocityX )[i][j][k] = static_cast<double>( i ) * 3.0;
               mat_block.second.GetPrimeStateBuffer( PrimeState::VelocityY )[i][j][k] = static_cast<double>( j ) * 4.0;
            }
         }
      }

      constexpr double cell_size = 1.0;
      WHEN( "The shear and bulk viscosity are zero" ) {
         constexpr double x_block_coordinate = 0.0;
         MaterialManager const material_manager( ReturnMaterialManagerWithViscosities( MaterialName::StiffenedGas, 0.0, 0.0 ) );
         AxisymmetricViscousVolumeForces const forces_calculator = AxisymmetricViscousVolumeForces( material_manager );
         double volume_forces[FF::ANOE()][CC::ICX()][CC::ICY()][CC::ICZ()];
         for( unsigned int e = 0; e < FF::ANOE(); ++e ) {
            for( unsigned int i = 0; i < CC::ICX(); ++i ) {
               for( unsigned int j = 0; j < CC::ICY(); ++j ) {
                  for( unsigned int k = 0; k < CC::ICZ(); ++k ) {
                     volume_forces[e][i][j][k] = 0.0;
                  }
               } //j
            } //i
         } //equation

         forces_calculator.ComputeForces( mat_block, volume_forces, cell_size, x_block_coordinate );

         THEN( "The volume forces are zero" ) {
            for( unsigned int e = 0; e < FF::ANOE(); ++e ) {
               for( unsigned int i = 0; i < CC::ICX(); ++i ) {
                  for( unsigned int j = 0; j < CC::ICY(); ++j ) {
                     REQUIRE( volume_forces[e][i][j][0] == Approx( 0.0 ) );
                  } //j
               } //i
            } //equation
         }
      }
      WHEN( "Viscosities are twice different" ) {
         constexpr double x_block_coordinate = 0.0;
         double first_volume_forces[FF::ANOE()][CC::ICX()][CC::ICY()][CC::ICZ()];
         double second_volume_forces[FF::ANOE()][CC::ICX()][CC::ICY()][CC::ICZ()];
         for( unsigned int e = 0; e < FF::ANOE(); ++e ) {
            for( unsigned int i = 0; i < CC::ICX(); ++i ) {
               for( unsigned int j = 0; j < CC::ICY(); ++j ) {
                  for( unsigned int k = 0; k < CC::ICZ(); ++k ) {
                     first_volume_forces[e][i][j][k] = 0.0;
                     second_volume_forces[e][i][j][k] = 0.0;
                  }
               } //j
            } //i
         } //equation

         MaterialManager const material_manager_first( ReturnMaterialManagerWithViscosities( MaterialName::StiffenedGas, 3.4, 0.0 ) );
         AxisymmetricViscousVolumeForces const first_forces_calculator = AxisymmetricViscousVolumeForces( material_manager_first );
         first_forces_calculator.ComputeForces( mat_block, first_volume_forces, cell_size, x_block_coordinate);

         MaterialManager const material_manager_second( ReturnMaterialManagerWithViscosities( MaterialName::StiffenedGas, 6.8, 0.0 ) );
         AxisymmetricViscousVolumeForces const second_forces_calculator = AxisymmetricViscousVolumeForces( material_manager_second );
         second_forces_calculator.ComputeForces( mat_block, second_volume_forces, cell_size, x_block_coordinate);

         THEN( "Volume forces are twice different" ) {
             for( unsigned int e = 0; e < FF::ANOE(); ++e ) {
                for( unsigned int i = 0; i < CC::ICX(); ++i ) {
                   for( unsigned int j = 0; j < CC::ICY(); ++j ) {
                       REQUIRE( 2.0 * first_volume_forces[e][i][j][0] == Approx( second_volume_forces[e][i][j][0] ) );
                   } //j
                } //i
             } //equation
         }
      }
      WHEN( "Cell-center coordinates are twice different due to the difference of the block position" ) {
         constexpr double first_x_block_coordinate = 0.0;
         constexpr double second_x_block_coordinate = 4.5;
         MaterialManager const material_manager( ReturnMaterialManagerWithViscosities( MaterialName::StiffenedGas, 5.0, 0.0 ) );
         AxisymmetricViscousVolumeForces const forces_calculator = AxisymmetricViscousVolumeForces( material_manager );
         double first_volume_forces[FF::ANOE()][CC::ICX()][CC::ICY()][CC::ICZ()];
         double second_volume_forces[FF::ANOE()][CC::ICX()][CC::ICY()][CC::ICZ()];
         for( unsigned int e = 0; e < FF::ANOE(); ++e ) {
            for( unsigned int i = 0; i < CC::ICX(); ++i ) {
               for( unsigned int j = 0; j < CC::ICY(); ++j ) {
                  for( unsigned int k = 0; k < CC::ICZ(); ++k ) {
                     first_volume_forces[e][i][j][k] = 0.0;
                     second_volume_forces[e][i][j][k] = 0.0;
                  }
               } //j
            } //i
         } //equation

         forces_calculator.ComputeForces( mat_block, first_volume_forces, cell_size, first_x_block_coordinate );
         forces_calculator.ComputeForces( mat_block, second_volume_forces, cell_size, second_x_block_coordinate );

         THEN( "r-momentum forces are twice different" ) {
            for( unsigned int i = 0; i < CC::ICX(); ++i ) {
               for( unsigned int j = 0; j < CC::ICY(); ++j ) {
                  REQUIRE( first_volume_forces[ETI( Equation::MomentumY )][i][j][0] == Approx( 2.0 * second_volume_forces[ETI( Equation::MomentumY )][i][j][0] ) );
               } //j
            } //i
         }
      }
   }
}