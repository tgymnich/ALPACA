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
#ifndef FIXED_VALUE_BOUNDARY_CONDITION_H
#define FIXED_VALUE_BOUNDARY_CONDITION_H

#include "fluid_boundary_condition.h"
#include "boundary_constants.h"

#include "levelset/multi_phase_manager/two_phase_manager.h"

/**
 * @brief The FixedValueBoundaryCondition class imposes a pre-defined value for each variable in all halo cells.
 * Works only for single-phase and single-level-set simulations!
 */
template<BoundaryLocation LOC>
class FixedValueBoundaryCondition : public FluidBoundaryCondition {

private:
   std::array<std::array<double,FF::ANOE()>, 2> const conservatives_;
   std::array<std::array<double,FF::ANOP()>, 2> const prime_states_;

public:
   FixedValueBoundaryCondition() = delete;
   /**
    * @brief Default constructor. See base class constructor.
    */
   explicit FixedValueBoundaryCondition( std::array<std::array<double,FF::ANOE()>, 2> const conservatives,
                                         std::array<std::array<double,FF::ANOP()>, 2> const prime_states ) :
      conservatives_(conservatives),
      prime_states_(prime_states) {
      // Empty besides initializer list
   }
   ~FixedValueBoundaryCondition() = default;
   FixedValueBoundaryCondition( FixedValueBoundaryCondition const& ) = delete;
   FixedValueBoundaryCondition& operator=( FixedValueBoundaryCondition const& ) = delete;
   FixedValueBoundaryCondition( FixedValueBoundaryCondition&& ) = delete;
   FixedValueBoundaryCondition& operator=( FixedValueBoundaryCondition&& ) = delete;

   /**
    * @brief Imposes predefined values onto the respective halo cells. See base class.
    */
   void UpdateFluidExternal( Node& node, FluidFieldType const field_type ) const override {
      auto start_indices = BoundaryConstants<LOC>::HaloStartIndices();
      auto end_indices = BoundaryConstants<LOC>::HaloEndIndices();

      for( auto& host_mat_block : node.GetPhases() ) {
         unsigned int const n = host_mat_block.first == MaterialSignCapsule::PositiveFluidMaterial() ? 0 : 1;
         switch( field_type ) {
            case FluidFieldType::Conservatives: {
               for( Equation const eq : FF::ASOE() ) {
                  double (&cells)[CC::TCX()][CC::TCY()][CC::TCZ()] = host_mat_block.second.GetRightHandSideBuffer( eq );
                  for( unsigned int i = start_indices[0]; i < end_indices[0]; ++i ) {
                     for( unsigned int j = start_indices[1]; j < end_indices[1]; ++j ) {
                        for( unsigned int k = start_indices[2]; k < end_indices[2]; ++k ) {
                           cells[i][j][k] = conservatives_[n][ETI( eq )];
                        }
                     }
                  }
               }
            }
            break;
#ifndef PERFORMANCE
            case FluidFieldType::PrimeStates: {
               for( PrimeState const ps : FF::ASOP() ) {
                  double (&cells)[CC::TCX()][CC::TCY()][CC::TCZ()] = host_mat_block.second.GetPrimeStateBuffer( ps );
                  for( unsigned int i = start_indices[0]; i < end_indices[0]; ++i ) {
                     for( unsigned int j = start_indices[1]; j < end_indices[1]; ++j ) {
                        for( unsigned int k = start_indices[2]; k < end_indices[2]; ++k ) {
                           cells[i][j][k] = prime_states_[n][PTI( ps )];
                        }
                     }
                  }
               }
            }
            break;
            default: 
               throw std::runtime_error( "Fluid field type not known!" );
#else 
            default: /* FluidFieldType::PrimeStates */ {
               for( PrimeState const ps : FF::ASOP() ) {
                  double (&cells)[CC::TCX()][CC::TCY()][CC::TCZ()] = host_mat_block.second.GetPrimeStateBuffer( ps );
                  for( unsigned int i = start_indices[0]; i < end_indices[0]; ++i ) {
                     for( unsigned int j = start_indices[1]; j < end_indices[1]; ++j ) {
                        for( unsigned int k = start_indices[2]; k < end_indices[2]; ++k ) {
                           cells[i][j][k] = prime_states_[n][PTI( ps )];
                        }
                     }
                  }
               }
            }  
#endif 
         }
      }
   }

   /**
    * @brief Identifies the Location of the BoundaryCondition.
    * @return A BoundaryLocation indicating the position of the BoundaryCondition, i. e. which halo cells are updated by it.
    */
   BoundaryLocation GetLocation() const {
      return LOC;
   }
};

#endif // FIXED_VALUE_BOUNDARY_CONDITION_H
