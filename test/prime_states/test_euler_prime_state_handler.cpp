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

#include "prime_states/euler_prime_state_handler.h"

SCENARIO( "Conservative quantities and prime state values can be converted into each other", "[1rank]" ) {
   std::unordered_map<std::string, double> const material_parameters = {
      {"gamma", 1.4 },
      {"A", 0.0},
      {"B", 0.0},
      {"C", 0.0},
      {"rho0", 0.0},
      {"specificGasConstant", 0.0},
      {"thermalConductivity", 0.0},
      {"dynamicShear", 0.0},
      {"dynamicBulk", 0.0}
   };
   constexpr MaterialName material = MaterialName::UserStiffenedGasOne;
   auto const material_manager = MaterialManager( {std::make_tuple(MaterialName::StiffenedGas, material, material_parameters)}, {1.0} );
   auto const euler_prime_state_handler = EulerPrimeStateHandler( material_manager );

   GIVEN( "A set of prime states" ) {
      std::array<double, FF::ANOP()> prime_states;
      for( unsigned int p = 0; p < FF::ANOP(); ++p ) {
         prime_states[p] = double(p + 1);
      }

      WHEN( "Prime states are converted to conservatives" ) {
         std::array<double, FF::ANOE()> conservatives;
         euler_prime_state_handler.ConvertPrimeStatesToConservatives( material, prime_states, conservatives );

         THEN( "Momenta equal the product of density and velocity" ) {
            for( unsigned int m = 0; m < DTI(CC::DIM()); ++m ) {
               REQUIRE( conservatives[ETI(FF::AME()[m])] == Approx( prime_states[PTI(PrimeState::Density)] * prime_states[PTI(FF::AV()[m])] ) );
            }
         }

         THEN( "Density equals mass" ) {
            REQUIRE( conservatives[ETI(Equation::Mass)] == Approx( prime_states[PTI(PrimeState::Density)] ) );
         }

         // Energy is not suitable to be checked since it depends on the used equation of state
      }

      WHEN( "Prime states are converted to conservatives and converted back to prime states" ) {
         std::array<double, FF::ANOE()> conservatives;
         std::array<double, FF::ANOP()> prime_states_new;
         euler_prime_state_handler.ConvertPrimeStatesToConservatives( material, prime_states, conservatives );
         euler_prime_state_handler.ConvertConservativesToPrimeStates( material, conservatives, prime_states_new );

         THEN( "The resulting prime states equal the original prime states" ) {
            for( unsigned int p = 0; p < FF::ANOP(); ++p ) {
               REQUIRE( prime_states_new[p] == Approx( prime_states[p] ) );
            }
         }
      }
   }
}
