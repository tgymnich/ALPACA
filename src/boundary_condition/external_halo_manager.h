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
#ifndef EXTERNAL_HALO_MANAGER_H
#define EXTERNAL_HALO_MANAGER_H

#include <array>
#include "simulation_setup.h"
#include "boundary_condition/boundary_specifications.h"
#include "boundary_condition/fluid_boundary_condition.h"
#include "boundary_condition/levelset_boundary_condition.h"
#include "topology/node.h"

/**
 * @brief Container of the external boundaries conditions for fluid and levelsets
 */
class ExternalHaloManager {

private:
   // Arrays holding the boundary conditions for all six natural sides. Always done also for reduced dimensions
   std::array<std::unique_ptr<FluidBoundaryCondition const>, 6> const fluid_boundary_conditions_;
   std::array<std::unique_ptr<LevelsetBoundaryCondition const>, 6> const levelset_boundary_conditions_;

   // Factory functions
   std::array<std::unique_ptr<FluidBoundaryCondition const>, 6> InitializeFluidBoundaryCondition( SimulationSetup const& setup );
   template<BoundaryLocation LOC>
   std::unique_ptr<FluidBoundaryCondition const> CreateFluidBoundary( FluidBoundaryType const fluid_type, SimulationSetup const& setup );

   std::array<std::unique_ptr<LevelsetBoundaryCondition const>, 6> InitializeLevelsetBoundaryCondition( SimulationSetup const& setup );
   template<BoundaryLocation LOC>
   std::unique_ptr<LevelsetBoundaryCondition const> CreateLevelsetBoundary( LevelSetBoundaryType const levelset_type );

public:
   ExternalHaloManager() = delete;
   explicit ExternalHaloManager( SimulationSetup const& setup );
   ~ExternalHaloManager() = default;
   ExternalHaloManager( ExternalHaloManager const& ) = delete;
   ExternalHaloManager& operator=( ExternalHaloManager const& ) = delete;
   ExternalHaloManager( ExternalHaloManager&& ) = delete;
   ExternalHaloManager& operator=( ExternalHaloManager&& ) = delete;

   // Functions performing the halo update on different buffer
   void UpdateLevelsetExternal( Node& node, LevelsetBlockBufferType const buffer_type, BoundaryLocation const loc ) const;
   void UpdateInterfaceTagExternal( Node& node, BoundaryLocation const loc ) const;
   void UpdateFluidExternal( Node& node, FluidFieldType const field_type, BoundaryLocation const loc ) const;
};

#endif /* EXTERNAL_HALO_MANAGER_H */
