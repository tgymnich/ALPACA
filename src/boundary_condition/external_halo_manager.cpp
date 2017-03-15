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
#include "boundary_condition/external_halo_manager.h"
#include "boundary_condition/symmetry_boundary_condition.h"
#include "boundary_condition/fixed_value_boundary_condition.h"
#include "boundary_condition/wall_boundary_condition.h"
#include "boundary_condition/zero_gradient_boundary_condition.h"
#include "topology/id_periodic_information.h"

/**
 * @brief Default constructor.
 * @param setup Setup holding all information contained in the input file
 */
ExternalHaloManager::ExternalHaloManager( SimulationSetup const& setup ) :
   fluid_boundary_conditions_( InitializeFluidBoundaryCondition( setup ) ),
   levelset_boundary_conditions_( InitializeLevelsetBoundaryCondition( setup ) ) {
   /** Empty besides initializer list */
}

/**
 * @brief Performs all levelset halo updates in external boundaries.
 * @param node Node for which the halo update is done (indirect return)
 * @param buffer_type Identifier of the levelset buffer that should be updated
 * @param loc Location of the boundary to update.
 */
void ExternalHaloManager::UpdateLevelsetExternal( Node& node, LevelsetBlockBufferType const buffer_type, BoundaryLocation const loc ) const {
   levelset_boundary_conditions_[LTI( loc )]->UpdateLevelsetExternal( node, buffer_type );
}

/**
 * @brief Performs all interface tag halo updates in external boundaries.
 * @param Node for which the halo update is done (indirect return)
 * @param loc Location of the boundary to update.
 */
void ExternalHaloManager::UpdateInterfaceTagExternal( Node& node, BoundaryLocation const loc ) const {
   levelset_boundary_conditions_[LTI( loc )]->UpdateInterfaceTagExternal( node );
}

/**
 * @brief Performs all interface tag halo updates in external boundaries.
 * @param Node for which the halo update is done (indirect return)
 * @param loc Location of the boundary to update.
 */
void ExternalHaloManager::UpdateFluidExternal( Node& node, FluidFieldType const field_type, BoundaryLocation const loc ) const {
   fluid_boundary_conditions_[LTI( loc )]->UpdateFluidExternal( node, field_type );
}

/**
* @brief Initializes the fluid boundary conditions with the given conditions, the conditions for the periodic locations are not created.
* @param setup Setup holding all information contained in the input file
* @return Array with all six pointers to the created fluid boundary conditions
*/
std::array<std::unique_ptr<FluidBoundaryCondition const>, 6> ExternalHaloManager::InitializeFluidBoundaryCondition( SimulationSetup const& setup ) {

   std::array<FluidBoundaryType, 6> const fluid_boundary_types( setup.GetFluidBoundaryConditions() );
   std::array<std::unique_ptr<FluidBoundaryCondition const>, 6> fluid_boundary_conditions;

   // Check periodicity in each direction
   if( !( setup.GetActivePeriodicLocations() & PeriodicBoundariesLocations::EastWest ) ) {
      fluid_boundary_conditions[LTI( BoundaryLocation::East )] = CreateFluidBoundary<BoundaryLocation::East>(
         fluid_boundary_types[LTI( BoundaryLocation::East )], setup );
      fluid_boundary_conditions[LTI( BoundaryLocation::West )] = CreateFluidBoundary<BoundaryLocation::West>(
         fluid_boundary_types[LTI( BoundaryLocation::West )], setup );
   }

   if( !( setup.GetActivePeriodicLocations() & PeriodicBoundariesLocations::NorthSouth ) ) {
      fluid_boundary_conditions[LTI( BoundaryLocation::North)] = CreateFluidBoundary<BoundaryLocation::North>(
         fluid_boundary_types[LTI( BoundaryLocation::North )], setup );
      fluid_boundary_conditions[LTI( BoundaryLocation::South )] = CreateFluidBoundary<BoundaryLocation::South>(
         fluid_boundary_types[LTI( BoundaryLocation::South )], setup );
   }

   if( !( setup.GetActivePeriodicLocations() & PeriodicBoundariesLocations::TopBottom ) ) {
      fluid_boundary_conditions[LTI( BoundaryLocation::Top )] = CreateFluidBoundary<BoundaryLocation::Top>(
         fluid_boundary_types[LTI( BoundaryLocation::Top )], setup );
      fluid_boundary_conditions[LTI( BoundaryLocation::Bottom)] = CreateFluidBoundary<BoundaryLocation::Bottom>(
         fluid_boundary_types[LTI( BoundaryLocation::Bottom )], setup );
   }

   return fluid_boundary_conditions;
}

/**
 * @brief Builds a fluid boundary condition of appropiate type.
 * @tparam LOC Template parameter for the boundary location.
 * @param fluid_boundary_type Boundary type for which the fluid boundary condition should be built
 * @param setup Setup holding all information contained in the input file
 * @return Pointer to the created fluid boundary condition
 */
template<BoundaryLocation LOC>
std::unique_ptr<FluidBoundaryCondition const> ExternalHaloManager::CreateFluidBoundary( FluidBoundaryType const fluid_boundary_type, SimulationSetup const& setup ) {
   //Switch with return do not need break.
   switch( fluid_boundary_type ) {
      case FluidBoundaryType::ZeroGradient :
         return std::make_unique<ZeroGradientBoundaryCondition<LOC> const>();
      case FluidBoundaryType::Symmetry :
         return std::make_unique<SymmetryBoundaryCondition<LOC> const>();
      case FluidBoundaryType::FixedValue :
         return std::make_unique<FixedValueBoundaryCondition<LOC> const>( setup.FixedValueBoundaryConservatives( LOC ), setup.FixedValueBoundaryPrimeStates( LOC ) );
      case FluidBoundaryType::Wall :
         return std::make_unique<WallBoundaryCondition<LOC> const>();
      // no performance flags required here since only used during initialization
      case FluidBoundaryType::Internal :
         throw std::logic_error( "Internals Halos are not managed by the ExternalHaloManager!" );
      case FluidBoundaryType::Periodic :
         throw std::logic_error( "Periodic boundary conditions are not managed by the ExternalHaloManager!" );
      default: 
         throw std::logic_error( "Fluid boundary type not known" );
   }
}

/*
* @brief Initializes the levelset boundary conditions with the given conditions, the conditions for the periodic locations are not created.
* @param setup Setup holding all information contained in the input file
* @return Array with all six pointers to the created levelset boundary conditions
*/
std::array<std::unique_ptr<LevelsetBoundaryCondition const>, 6> ExternalHaloManager::InitializeLevelsetBoundaryCondition( SimulationSetup const& setup ) {

   std::array<LevelSetBoundaryType, 6> const levelset_boundary_types( setup.GetLevelsetBoundaryConditions() );
   std::array<std::unique_ptr<LevelsetBoundaryCondition const>, 6> levelset_boundary_conditions;

   // Check periodicity in each direction
   if( !( setup.GetActivePeriodicLocations() & PeriodicBoundariesLocations::EastWest ) ) {
      levelset_boundary_conditions[LTI( BoundaryLocation::East )] = CreateLevelsetBoundary<BoundaryLocation::East>(
         levelset_boundary_types[LTI( BoundaryLocation::East )] );
      levelset_boundary_conditions[LTI( BoundaryLocation::West )] = CreateLevelsetBoundary<BoundaryLocation::West>(
         levelset_boundary_types[LTI( BoundaryLocation::West )] );
   }

   if( !( setup.GetActivePeriodicLocations() & PeriodicBoundariesLocations::NorthSouth ) ) {
      levelset_boundary_conditions[LTI( BoundaryLocation::North )] = CreateLevelsetBoundary<BoundaryLocation::North>(
         levelset_boundary_types[LTI( BoundaryLocation::North )] );
      levelset_boundary_conditions[LTI( BoundaryLocation::South )] = CreateLevelsetBoundary<BoundaryLocation::South>(
         levelset_boundary_types[LTI( BoundaryLocation::South )] );
   }

   if( !( setup.GetActivePeriodicLocations() & PeriodicBoundariesLocations::TopBottom ) ) {
      levelset_boundary_conditions[LTI( BoundaryLocation::Top )] = CreateLevelsetBoundary<BoundaryLocation::Top>(
         levelset_boundary_types[LTI( BoundaryLocation::Top )] );
      levelset_boundary_conditions[LTI( BoundaryLocation::Bottom )] = CreateLevelsetBoundary<BoundaryLocation::Bottom>(
         levelset_boundary_types[LTI( BoundaryLocation::Bottom )] );
   }

   return levelset_boundary_conditions;
}

/**
 * @brief Builds a fluid boundary condition of appropiate type.
 * @tparam LOC Template parameter for the boundary location.
 * @param levelset_boundary_type Boundary type for which the levelset boundary condition should be built
 * @return Pointer to the created levelset boundary condition
 */
template<BoundaryLocation LOC>
std::unique_ptr<LevelsetBoundaryCondition const> ExternalHaloManager::CreateLevelsetBoundary( LevelSetBoundaryType const levelset_boundary_type ) {

   //Switch with return do not need break.
   switch( levelset_boundary_type ) {
      case LevelSetBoundaryType::ZeroGradient :
         return std::make_unique<ZeroGradientBoundaryCondition<LOC> const>();
      case LevelSetBoundaryType::Symmetry :
         return std::make_unique<SymmetryBoundaryCondition<LOC> const>();
      // no performance flags required here since only used during initialization
      case LevelSetBoundaryType::Internal :
         throw std::logic_error( "Internals Halos are not managed by the ExternalHaloManager!" );
      case LevelSetBoundaryType::Periodic :
         throw std::logic_error( "Periodic boundary conditions are not managed by the ExternalHaloManager!" );
      default: 
         throw std::logic_error( "Levelset boundary type not known" );
   }
}