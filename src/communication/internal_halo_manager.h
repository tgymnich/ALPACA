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
#ifndef INTERNAL_BOUNDARY_MANAGER_H
#define INTERNAL_BOUNDARY_MANAGER_H

#include "communication/communication_manager.h"
#include "communication/communication_statistics.h"

/**
 * @brief The InternalHaloManager is used within the domain, i.e. a classical Halo. This class exchanges information of neighboring blocks
 *        by filling the halos of the host with data from the domain of the neighbor. In the MR setup two types of internal boundaries are
 *        possible, jump and no-jump; both are treated by this class. In the first case communication with the parent of the host node is
 *        needed in the latter communication with the direct neighbor is needed. Works inter and intra MPI ranks.
 */
class InternalHaloManager {
private:
   Tree& tree_;
   TopologyManager& topology_;
   CommunicationManager& communication_manager_; // Cannot be const (for now NH TODO-19) because of new tagging system.
   unsigned int const number_of_fluids_;

   // Helper function for the local halo filling of special buffers
   template<class T>
   inline void UpdateNoJumpLocal( T (&host_buffer)[CC::TCX()][CC::TCY()][CC::TCZ()],
                                  T const (& partner_buffer)[CC::TCX()][CC::TCY()][CC::TCZ()], BoundaryLocation const loc ) const;

   template<class T>
   inline void ExtendClosestInternalValue( T (& host_buffer)[CC::TCX()][CC::TCY()][CC::TCZ()], BoundaryLocation const loc ) const;

   // Fluid specific functions 
   unsigned int UpdateFluidJumpMpiSend( std::uint64_t id, std::vector<MPI_Request>& requests, std::uint64_t const remote_child_id,
                                        void *send_buffer, BoundaryLocation const loc, FluidFieldType const filed_type );
   void UpdateFluidJumpMpiRecv( std::uint64_t id, std::vector<MPI_Request>& requests, BoundaryLocation const loc,
                                FluidFieldType const field_type );
   void UpdateFluidJumpNoMpi( std::uint64_t id, BoundaryLocation const loc, FluidFieldType const field_type );
   void UpdateFluidHaloCellsMpiSend( std::uint64_t id, std::vector<MPI_Request>& requests, BoundaryLocation const loc,
                                     FluidFieldType const field_type );
   void UpdateFluidHaloCellsMpiRecv( std::uint64_t id, std::vector<MPI_Request>& requests, BoundaryLocation const loc,
                                     FluidFieldType const field_type );
   void UpdateFluidHaloCellsNoMpi( std::uint64_t id, BoundaryLocation const loc, FluidFieldType const field_type );
   void MpiFluidHaloUpdateNoJump( std::vector<MPI_Request>& requests,
                                  std::vector<std::tuple<std::uint64_t, BoundaryLocation, InternalBoundaryType>> const& no_jump_boundaries,
                                  FluidFieldType const field_type );
   void NoMpiFluidHaloUpdate( std::vector<std::tuple<std::uint64_t, BoundaryLocation, InternalBoundaryType>> const& boundaries,
                              FluidFieldType const field_type );
   void MpiFluidHaloUpdateJump( std::vector<MPI_Request>& requests,
                                std::vector<std::tuple<std::uint64_t, BoundaryLocation, InternalBoundaryType>> const& no_jump_boundaries,
                                std::vector<ExchangePlane>& jump_buffer_plane, std::vector<ExchangeStick>& jump_buffer_stick,
                                std::vector<ExchangeCube>& jump_buffer_cube,
                                FluidFieldType const field_type );

   // Levelset specific functions
   void UpdateLevelsetHaloCellsMpiSend( std::uint64_t id, std::vector<MPI_Request>& requests, LevelsetBlockBufferType const buffer_type,
                                        BoundaryLocation const loc );
   void UpdateLevelsetHaloCellsMpiRecv( std::uint64_t id, std::vector<MPI_Request>& requests, LevelsetBlockBufferType const buffer_type,
                                        BoundaryLocation const loc );
   void UpdateLevelsetHaloCellsNoMpi( std::uint64_t id, LevelsetBlockBufferType const buffer_type, BoundaryLocation const loc );
   void NoMpiLevelsetHaloUpdate( std::vector<std::tuple<std::uint64_t, BoundaryLocation, InternalBoundaryType>> const& boundaries,
                                 LevelsetBlockBufferType const buffer_type );
   void MpiLevelsetHaloUpdate( std::vector<std::tuple<std::uint64_t, BoundaryLocation, InternalBoundaryType>> const& boundaries,
                               LevelsetBlockBufferType const buffer_type, std::vector<MPI_Request>& requests );

   // Interface tag specific functions
   void UpdateInterfaceTagHaloCellsMpiSend( std::uint64_t id, std::vector<MPI_Request>& requests, BoundaryLocation const loc );
   void UpdateInterfaceTagHaloCellsMpiRecv( std::uint64_t id, std::vector<MPI_Request>& requests, BoundaryLocation const loc );
   void UpdateInterfaceTagHaloCellsNoMpi( std::uint64_t id, BoundaryLocation const loc );
   void NoMpiInterfaceTagHaloUpdate( std::vector<std::tuple<std::uint64_t, BoundaryLocation, InternalBoundaryType>> const& boundaries );
   void MpiInterfaceTagHaloUpdate( std::vector<std::tuple<std::uint64_t, BoundaryLocation, InternalBoundaryType>> const& boundaries,
                                   std::vector<MPI_Request>& requests );

public:
   InternalHaloManager() = delete;
   explicit InternalHaloManager( Tree& tree, TopologyManager& topology, CommunicationManager& communication_manager, unsigned int const number_of__fluids );
   ~InternalHaloManager() = default;
   InternalHaloManager( InternalHaloManager const&) = delete;
   InternalHaloManager& operator=( InternalHaloManager const& ) = delete;
   InternalHaloManager( InternalHaloManager&& ) = delete;
   InternalHaloManager& operator=( InternalHaloManager&& ) = delete;

   // Function to trigger halo update of specific types
   void FluidHaloUpdateOnLevel( unsigned int const level, FluidFieldType const field_type, bool const cut_jumps );
   void InterfaceTagHaloUpdateOnLevel( unsigned int const level );
   void LevelsetHaloUpdateOnLevel( unsigned int const level, LevelsetBlockBufferType const type );
};

#endif // INTERNAL_BOUNDARY_MANAGER_H
