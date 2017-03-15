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
#include "communication/exchange_types.h"
#include "multiresolution/multiresolution.h"
#include "topology/id_information.h"
#include "levelset/multi_phase_manager/material_sign_capsule.h"
#include "halo_manager.h"

/**
 * @brief Default constructor.
 * @param setup .
 * @param tree .
 * @param external_boundary_manager .
 * @param communication_manager .
 * @param maximum_level .
 */
HaloManager::HaloManager( Tree& tree, ExternalHaloManager const& external_halo_manager, InternalHaloManager& internal_halo_manager,
                          CommunicationManager& communication_manager, unsigned int const maximum_level ) :
   tree_( tree ),
   external_halo_manager_( external_halo_manager ),
   internal_halo_manager_( internal_halo_manager ),
   communication_manager_( communication_manager ),
   maximum_level_( maximum_level )
{
   // Empty besides initializer list.
}

/**
 * @brief Adjusts the values in the halo cells, according to their type (symmetry, internal, ...).
 * @param levels_ascending The levels on which halos of nodes will be modified in ascending order.
 * @param field_type The decider whether a halo update for conservatives or for prime states is done.
 * @param cut_jumps Decider if jump halos should be updated on all specified level. If true: jumps will not be updated on the coarsest level in "upddate_levels".
 */
void HaloManager::FluidHaloUpdate( std::vector<unsigned int> const& levels_ascending, FluidFieldType const field_type, bool const cut_jumps ) {
   std::vector<unsigned int> no_jump_update_levels(levels_ascending);
   /* NH 2017-02-20: It may be that no-jump halos are not to be updated on the coarsest level in the input list. Therefore this level is handled separately.
    * Afterwards a normal Halo update is performed on all remaining levels in the input list.
    */
   if( cut_jumps ) {
      unsigned int no_jump_extra_level = no_jump_update_levels.front();
      no_jump_update_levels.erase( no_jump_update_levels.begin() );
      FluidHaloUpdateOnLevel( no_jump_extra_level, field_type, true );
   }
   for( unsigned int const level : no_jump_update_levels ) {
      FluidHaloUpdateOnLevel( level, field_type, false );
   }
}

/**
 * @brief Adjusts the fluid values in all halo cells, according to their type
 * @param level The level on which halos of nodes will be modified.
 * @param field_type The decider whether a halo update for conservatives or for prime states is done.
 * @param cut_jumps Decider if jump halos should be updated on specified level. If true: jumps will not be updated on the current level.
 */
void HaloManager::FluidHaloUpdateOnLevel( unsigned int const level, FluidFieldType const field_type, bool const cut_jumps ) {
   FluidInternalHaloUpdateOnLevel( level, field_type, cut_jumps );
   FluidExternalHaloUpdateOnLevel( level, field_type );
}

/**
 * @brief Adjusts the values in internal halo cells, according to their type
 * @param level The level on which halos of nodes will be modified.
 * @param field_type The decider whether a halo update for conservatives or for prime states is done.
 * @param cut_jumps Decider if jump halos should be updated on specified level. If true: jumps will not be updated on the current level.
 */
void HaloManager::FluidInternalHaloUpdateOnLevel( unsigned int const level, FluidFieldType const field_type, bool const cut_jumps ) {
   internal_halo_manager_.FluidHaloUpdateOnLevel( level, field_type, cut_jumps );
}

/**
 * @brief Adjusts the fluid values in external halo cells, according to their type
 * @param level The level on which halos of nodes will be modified.
 * @param field_type The decider whether a halo update for conservatives or for prime states is done.
 */
void HaloManager::FluidExternalHaloUpdateOnLevel( unsigned int const level, FluidFieldType const field_type ) {
   for( std::tuple<std::uint64_t, BoundaryLocation> const& boundary : communication_manager_.ExternalBoundaries( level ) ) {
      external_halo_manager_.UpdateFluidExternal( tree_.GetNodeWithId( std::get<0>( boundary ) ), field_type, std::get<1>( boundary ) );
   }
}

/**
 * @brief Adjusts the fluid values in the halo cells on the finest level, according to their type.
 * @param field_type The decider whether a halo update for conservatives or for prime states is done.
 * @param cut_jumps Decider if jump halos should be updated on specified level. If true: jumps will not be updated on the current level.
 * @note The default value for cut_jumps is true.
 */
void HaloManager::FluidHaloUpdateOnLmax( FluidFieldType const field_type, bool const cut_jumps ) {
   FluidHaloUpdateOnLevel( maximum_level_, field_type, cut_jumps );
}

/**
 * @brief Calls an interface tag halo update on Lmax only.
 */
void HaloManager::InterfaceTagHaloUpdateOnLmax() {
   InterfaceTagHaloUpdateOnLevelList( { maximum_level_ } );
}

/**
 * @brief Calls a Levelset halo update of the specified "levelset-tag (levelset + interface tags)" buffer on Lmax (only!).
 * @param type The identifier of the buffer that is to be updated.
 */
void HaloManager::LevelsetHaloUpdateOnLmax( LevelsetBlockBufferType const type ) {
   // perform levelset halo update on Lmax
   LevelsetHaloUpdateOnLevelList( { maximum_level_ }, type );
}

/**
 * @brief Adjusts the values in the interface tag buffer according to their type (symmetry, internal ...).
 * @param updated_levels The levels on which halos of nodes will be modified.
 */
void HaloManager::InterfaceTagHaloUpdateOnLevelList( std::vector<unsigned int> const& updated_levels ) const {
   for( unsigned int const &level : updated_levels ) {
      internal_halo_manager_.InterfaceTagHaloUpdateOnLevel( level );
      // Update of domain boundaries
      for( auto const& domain_boundary : communication_manager_.ExternalBoundaries( level ) ) {
         std::uint64_t const id = std::get<0>( domain_boundary );
         BoundaryLocation const location = std::get<1>( domain_boundary );
         external_halo_manager_.UpdateInterfaceTagExternal( tree_.GetNodeWithId( id ), location );
      }
   } //levels
}

/**
 * @brief Adjusts the values in the stated levelset block buffer according to their type (symmetry, internal ...).
 * @param updated_levels The levels on which halos of nodes will be modified.
 * @param type The identifier of the buffer that is to be updated.
 */
void HaloManager::LevelsetHaloUpdateOnLevelList( std::vector<unsigned int> const updated_levels, LevelsetBlockBufferType const type ) {
   for( auto const& level : updated_levels ) {
      internal_halo_manager_.LevelsetHaloUpdateOnLevel( level, type );
      // Update of domain boundaries
      for( auto const& domain_boundary : communication_manager_.ExternalBoundaries( level ) ) {
         std::uint64_t const id = std::get<0>( domain_boundary );
         BoundaryLocation const location = std::get<1>( domain_boundary );
         external_halo_manager_.UpdateLevelsetExternal( tree_.GetNodeWithId( id ), type, location );
      }
   } //levels

}