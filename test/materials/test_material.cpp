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
#include "materials/material.h"
#include "materials/stiffened_gas.h"
#include "materials/stiffened_gas_complete_safe.h"
#include "materials/stiffened_gas_safe.h"
#include "materials/waterlike_fluid.h"

SCENARIO( "Materials are constructed and their properties are queried", "[1rank]" ) {

   GIVEN( "A stiffened gas" ) {

      WHEN( "All parameters are set to zero" ) {

         StiffenedGas const material = StiffenedGas( 0.0, 0.0, 0.0, 0.0 );

         THEN( "The temperature is" ) {
            REQUIRE( material.GetTemperature( 0.0, 0.0, 0.0, 0.0, 0.0 ) == -1 );
         }

         THEN( "The gruneisen is" ) {
            REQUIRE( material.GetGruneisen() == -1 );
         }

         THEN( "The psi is" ) {
            REQUIRE( material.GetPsi( 0.0, 0.0 ) == 0 );
         }

         THEN( "The gamma is" ) {
            REQUIRE( material.GetGamma() == 0 );
         }

         THEN( "The background pressure is" ) {
            REQUIRE( material.GetB() == 0 );
         }

         THEN( "The viscosity is" ) {
            std::vector<double> const expected( 2, 0.0 );
            REQUIRE( material.GetViscosity() == expected );
         }

         THEN( "The thermal conductivity is" ) {
            REQUIRE( material.GetThermalConductivity() == 0 );
         }

         THEN( "The specific heat is" ) {
            REQUIRE( material.GetSpecificHeat() == -1 );
         }
      }

      WHEN( "The parameters are reasonable values" ) {

         StiffenedGas const material = StiffenedGas( 6.1, 4.3, 2.0, 3.0 );

         THEN( "The pressure is" ) {
            REQUIRE( material.GetPressure( 1.0, 2.0, 3.0, 4.0, 5.0 ) == Approx( -74.68 ));
         }

         THEN( "The enthalpy is" ) {
            REQUIRE( material.GetEnthalpy( 6.0, 7.0, 8.0, 9.0, 10.0 ) == Approx( -7.9466666667 ));
         }

         THEN( "The energy is" ) {
            REQUIRE( material.GetEnergy( 11.0, 12.0, 13.0, 14.0, 15.0 ) == Approx( 31.2206773619 ));
         }

         THEN( "The temperature is" ) {
            REQUIRE( material.GetTemperature( 16.0, 17.0, 18.0, 19.0, 20.0 ) == -1 );
         }

         THEN( "The gruneisen is" ) {
            REQUIRE( material.GetGruneisen() == 5.1 );
         }

         THEN( "The psi is" ) {
            REQUIRE( material.GetPsi( 21.0, 22.0 ) == 1039.06 );
         }

         THEN( "The gamma is" ) {
            REQUIRE( material.GetGamma() == 6.1 );
         }

         THEN( "The background pressure is" ) {
            REQUIRE( material.GetB() == 4.3 );
         }

         THEN( "The speed of sound is" ) {
            REQUIRE( material.GetSpeedOfSound( 23.0, 24.0 ) == Approx( 2.7396445342 ));
         }

         THEN( "The viscosity is" ) {
            std::vector<double> const expected( { 2.0, 3.0 } );
            REQUIRE( material.GetViscosity() == expected );
         }

         THEN( "The thermal conductivity is" ) {
            REQUIRE( material.GetThermalConductivity() == 0 );
         }

         THEN( "The specific heat is" ) {
            REQUIRE( material.GetSpecificHeat() == -1 );
         }
      }
   }

   GIVEN( "A stiffened gas complete safe" ) {

      WHEN( "All parameters are set to zero" ) {

         StiffenedGasCompleteSafe const material = StiffenedGasCompleteSafe( 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 );

         THEN( "The pressure is" ) {
            REQUIRE( material.GetPressure( 0.0, 0.0, 0.0, 0.0, 0.0 ) == Approx( 0.0 ).margin( 0.000001 ));
         }

         THEN( "The enthalpy is" ) {
            REQUIRE( material.GetEnthalpy( 0.0, 0.0, 0.0, 0.0, 0.0 ) == 1.0 );
         }

         THEN( "The energy is" ) {
            REQUIRE( material.GetEnergy( 0.0, 0.0, 0.0, 0.0, 0.0 ) == 0.0 );
         }

         THEN( "The gruneisen is" ) {
            REQUIRE( material.GetGruneisen() == -1 );
         }

         THEN( "The psi is" ) {
            REQUIRE( material.GetPsi( 0.0, 0.0 ) == 0 );
         }

         THEN( "The gamma is" ) {
            REQUIRE( material.GetGamma() == 0 );
         }

         THEN( "The background pressure is" ) {
            REQUIRE( material.GetB() == 0 );
         }

         THEN( "The speed of sound is" ) {
            REQUIRE( material.GetSpeedOfSound( 0.0, 0.0 ) == Approx( 0.0000000149 ).epsilon( 0.1 ));
         }

         THEN( "The viscosity is" ) {
            std::vector<double> const expected( 2, 0.0 );
            REQUIRE( material.GetViscosity() == expected );
         }

         THEN( "The thermal conductivity is" ) {
            REQUIRE( material.GetThermalConductivity() == 0 );
         }

         THEN( "The specific heat is" ) {
            REQUIRE( material.GetSpecificHeat() == 0.0 );
         }
      }

      WHEN( "The parameters contain reasonable values" ) {

         StiffenedGasCompleteSafe const material = StiffenedGasCompleteSafe( 6.1, 1.0, 4.3, 2.0, 3.0, 4.0, 2.0, 3.0 );

         THEN( "The pressure is" ) {
            REQUIRE( material.GetPressure( 1.0, 2.0, 3.0, 4.0, 5.0 ) == -4.3 );
         }

         THEN( "The enthalpy is" ) {
            REQUIRE( material.GetEnthalpy( 6.0, 7.0, 8.0, 9.0, 10.0 ) == Approx( 0.95 ));
         }

         THEN( "The energy is" ) {
            REQUIRE( material.GetEnergy( 11.0, 12.0, 13.0, 14.0, 15.0 ) == Approx( 42.2206773619 ));
         }

         THEN( "The temperature is" ) {
            REQUIRE( material.GetTemperature( 16.0, 17.0, 18.0, 19.0, 20.0 ) == Approx(  4704254.7119584298 ));
         }

         THEN( "The gruneisen is" ) {
            REQUIRE( material.GetGruneisen() == 5.1 );
         }

         THEN( "The psi is" ) {
            REQUIRE( material.GetPsi( 21.0, 22.0 ) == 1039.06 );
         }

         THEN( "The gamma is" ) {
            REQUIRE( material.GetGamma() == 6.1 );
         }

         THEN( "The background pressure is" ) {
            REQUIRE( material.GetB() == 4.3 );
         }

         THEN( "The speed of sound is" ) {
            REQUIRE( material.GetSpeedOfSound( 23.0, 24.0 ) == Approx( 2.7396445342 ));
         }

         THEN( "The viscosity is" ) {
            std::vector<double> const expected( { 2.0, 3.0 } );
            REQUIRE( material.GetViscosity() == expected );
         }

         THEN( "The thermal conductivity is" ) {
            REQUIRE( material.GetThermalConductivity() == 4.0 );
         }

         THEN( "The specific heat is" ) {
            REQUIRE( material.GetSpecificHeat() == Approx( 0.5882352941 ));
         }
      }
   }

   GIVEN( "A stiffened gas safe" ) {

      WHEN( "All parameters are set to zero" ) {

         StiffenedGasSafe const material = StiffenedGasSafe( 0.0, 0.0, 0.0, 0.0 );

         THEN( "The pressure is" ) {
            REQUIRE( material.GetPressure( 0.0, 0.0, 0.0, 0.0, 0.0 ) == Approx( 0.0 ).margin( 0.000001 ));
         }

         THEN( "The enthalpy is" ) {
            REQUIRE( material.GetEnthalpy( 0.0, 0.0, 0.0, 0.0, 0.0 ) == 1.0 );
         }

         THEN( "The energy is" ) {
            REQUIRE( material.GetEnergy( 0.0, 0.0, 0.0, 0.0, 0.0 ) == 0.0 );
         }

         THEN( "The temperature is" ) {
            REQUIRE( material.GetTemperature( 0.0, 0.0, 0.0, 0.0, 0.0 ) == -1.0 );
         }

         THEN( "The gruneisen is" ) {
            REQUIRE( material.GetGruneisen() == -1 );
         }

         THEN( "The psi is" ) {
            REQUIRE( material.GetPsi( 0.0, 0.0 ) == 0 );
         }

         THEN( "The gamma is" ) {
            REQUIRE( material.GetGamma() == 0 );
         }

         THEN( "The background pressure is" ) {
            REQUIRE( material.GetB() == 0 );
         }

         THEN( "The speed of sound is" ) {
            REQUIRE( material.GetSpeedOfSound( 0.0, 0.0 ) == Approx( 0.0000000149 ).epsilon( 0.1 ));
         }

         THEN( "The viscosity is" ) {
            std::vector<double> const expected( 2, 0.0 );
            REQUIRE( material.GetViscosity() == expected );
         }

         THEN( "The thermal conductivity is" ) {
            REQUIRE( material.GetThermalConductivity() == 0 );
         }

         THEN( "The specific heat is" ) {
            REQUIRE( material.GetSpecificHeat() == -1.0 );
         }
      }

      WHEN( "The parameters contain reasonable values" ) {

         StiffenedGasSafe const material = StiffenedGasSafe( 6.1, 4.3, 2.0, 3.0 );

         THEN( "The pressure is" ) {
            REQUIRE( material.GetPressure( 1.0, 2.0, 3.0, 4.0, 5.0 ) == -4.3 );
         }

         THEN( "The enthalpy is" ) {
            REQUIRE( material.GetEnthalpy( 6.0, 7.0, 8.0, 9.0, 10.0 ) == Approx( 0.95 ));
         }

         THEN( "The energy is" ) {
            REQUIRE( material.GetEnergy( 11.0, 12.0, 13.0, 14.0, 15.0 ) == Approx( 31.2206773619 ));
         }

         THEN( "The temperature is" ) {
            REQUIRE( material.GetTemperature( 16.0, 17.0, 18.0, 19.0, 20.0 ) == -1.0 );
         }

         THEN( "The gruneisen is" ) {
            REQUIRE( material.GetGruneisen() == 5.1 );
         }

         THEN( "The psi is" ) {
            REQUIRE( material.GetPsi( 21.0, 22.0 ) == 1039.06 );
         }

         THEN( "The gamma is" ) {
            REQUIRE( material.GetGamma() == 6.1 );
         }

         THEN( "The background pressure is" ) {
            REQUIRE( material.GetB() == 4.3 );
         }

         THEN( "The speed of sound is" ) {
            REQUIRE( material.GetSpeedOfSound( 23.0, 24.0 ) == Approx( 2.7396445342 ));
         }

         THEN( "The viscosity is" ) {
            std::vector<double> const expected( { 2.0, 3.0 } );
            REQUIRE( material.GetViscosity() == expected );
         }

         THEN( "The thermal conductivity is" ) {
            REQUIRE( material.GetThermalConductivity() == 0.0 );
         }

         THEN( "The specific heat is" ) {
            REQUIRE( material.GetSpecificHeat() == -1.0 );
         }
      }
   }


   GIVEN( "A waterlike fluid" ) {

      WHEN( "All parameters are set to zero" ) {

         WaterlikeFluid const material = WaterlikeFluid( 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 );

         THEN( "The enthalpy is" ) {
            REQUIRE( material.GetEnthalpy( 0.0, 0.0, 0.0, 0.0, 0.0 ) == 0.0 );
         }

         THEN( "The temperature is" ) {
            REQUIRE( material.GetTemperature( 0.0, 0.0, 0.0, 0.0, 0.0 ) == -1.0 );
         }

         THEN( "The gruneisen is" ) {
            REQUIRE( material.GetGruneisen() == 0.0 );
         }

         THEN( "The psi is" ) {
            REQUIRE( material.GetPsi( 0.0, 0.0 ) == 0 );
         }

         THEN( "The gamma is" ) {
            REQUIRE( material.GetGamma() == -1.0 );
         }

         THEN( "The background pressure is" ) {
            REQUIRE( material.GetB() == -1.0 );
         }

         THEN( "The viscosity is" ) {
            std::vector<double> const expected( 2, 0.0 );
            REQUIRE( material.GetViscosity() == expected );
         }

         THEN( "The thermal conductivity is" ) {
            REQUIRE( material.GetThermalConductivity() == 0 );
         }

         THEN( "The specific heat is" ) {
            REQUIRE( material.GetSpecificHeat() == -1.0 );
         }
      }

      WHEN( "The parameters contain reasonable values" ) {

         WaterlikeFluid const material = WaterlikeFluid( 6.1, 1.0, 4.3, 1.0, 2.0, 3.0 );

         THEN( "The pressure is" ) {
            REQUIRE( material.GetPressure( 1.0, 2.0, 3.0, 4.0, 5.0 ) == 1.0 );
         }

         THEN( "The enthalpy is" ) {
            REQUIRE( material.GetEnthalpy( 6.0, 7.0, 8.0, 9.0, 10.0 ) == 0.0 );
         }

         THEN( "The energy is" ) {
            REQUIRE( material.GetEnergy( 11.0, 12.0, 13.0, 14.0, 15.0 ) == Approx( 22.8481283422 ));
         }

         THEN( "The temperature is" ) {
            REQUIRE( material.GetTemperature( 16.0, 17.0, 18.0, 19.0, 20.0 ) == -1.0 );
         }

         THEN( "The gruneisen is" ) {
            REQUIRE( material.GetGruneisen() == 0.0 );
         }

         THEN( "The psi is" ) {
            REQUIRE( material.GetPsi( 21.0, 22.0 ) == 3261.06 );
         }

         THEN( "The gamma is" ) {
            REQUIRE( material.GetGamma() == -1.0 );
         }

         THEN( "The background pressure is" ) {
            REQUIRE( material.GetB() == -1.0 );
         }

         THEN( "The speed of sound is" ) {
            REQUIRE( material.GetSpeedOfSound( 23.0, 24.0 ) == Approx( 2.690805601 ));
         }

         THEN( "The viscosity is" ) {
            std::vector<double> const expected( { 2.0, 3.0 } );
            REQUIRE( material.GetViscosity() == expected );
         }

         THEN( "The thermal conductivity is" ) {
            REQUIRE( material.GetThermalConductivity() == 0.0 );
         }

         THEN( "The specific heat is" ) {
            REQUIRE( material.GetSpecificHeat() == -1.0 );
         }
      }
   }
}
