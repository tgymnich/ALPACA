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
#include "numerous_phase_manager.h"

/**
 * @brief Default constructor for the NumerousPhaseManager. Calls the default constructor of the base class.
 * @param material_manager Instance of a material manager, which already has been inialized according to the user input.
 * @param halo_manager Instance to a HaloManager which provides MPI-related methods.
 */
NumerousPhaseManager::NumerousPhaseManager( MaterialManager const& material_manager, HaloManager& halo_manager ) :
   MultiPhaseManager( material_manager, halo_manager )
{
   // Empty Constructor, besides call of base class constructor.
}

/**
 * @brief Implements a cut-cell mixing procedure for multi-level set simulations. See also base class.
 * @param nodes See base class.
 * @param stage See base class.
 */
void NumerousPhaseManager::MixImplementation(std::vector<std::reference_wrapper<Node>> const& nodes, unsigned int const stage) const {
   for(Node& node : nodes) {
      cut_cell_mixer_.Mix(node,stage);
   }
}

/**
 * @brief See base class.
 * @param nodes See base class.
 * @param stage See base class.
 * @param is_last_stage See base class.
 */
void NumerousPhaseManager::EnforceWellResolvedDistanceFunctionImplementation(std::vector<std::reference_wrapper<Node>> const& nodes, unsigned int const stage, const bool is_last_stage) const {
   (void) is_last_stage; //Avoid compiler warning
   levelset_reinitializer_.Reinitialize(nodes,stage);
}

/**
 * @brief See base class.
 * @param node See base class.
 * @param stage See base class.
 */
void NumerousPhaseManager::ExtendImplementation(std::vector<std::reference_wrapper<Node>> const& nodes, unsigned int const stage) const {
   ghost_fluid_extender_.Extend(nodes,stage);
}

/**
 * @brief See base class.
 * @param nodes See base class.
 * @param stage See base class.
 */
void NumerousPhaseManager::PropagateLevelsetImplementation(std::vector<std::reference_wrapper<Node>> const& nodes, unsigned int const stage) const {
   (void) nodes; // Avoid compiler warning.
   (void) stage; // Avoid compiler warning.
   throw std::logic_error("Not yet implemented: NumerousPhaseManager::PropagateLevelset");
}

/**
 * @brief See base class.
 * @param nodes See base class.
 */
void NumerousPhaseManager::InitializeVolumeFractionBufferImplementation(std::vector<std::reference_wrapper<Node>> const& nodes) const {
   (void) nodes; // Avoid compiler warning.
   throw std::logic_error("Not yet implemented: NumerousPhaseManager::InitializeVolumeFractionBuffer");
}

/**
 * @brief See base class.
 * @param nodes See base class.
 */
void NumerousPhaseManager::ObtainInterfaceQuantitiesImplementation(std::vector<std::reference_wrapper<Node>> const& nodes, const bool reset_interface_states) const {
   (void) nodes; // Avoid compiler warning.
   (void) reset_interface_states; // Avoid compiler warning.
   throw std::logic_error("Not yet implemented: NumerousPhaseManager::ObtainInterfaceQuantities");
}
