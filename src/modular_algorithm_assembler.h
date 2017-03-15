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
#ifndef MODULAR_ALGORITHM_ASSEMBLER_H
#define MODULAR_ALGORITHM_ASSEMBLER_H

#include "halo_manager.h"
#include "integrator/time_integrator_setup.h"
#include "levelset/multi_phase_manager/multi_phase_manager_setup.h"
#include "prime_states/prime_state_handler_setup.h"
#include "solvers/space_solver.h"
#include "input_output/input_output_manager.h"
#include "communication/communication_manager.h"
#include "multiresolution/multiresolution.h"
#include "multiresolution/averager.h"

using TimeIntegratorConcretization = TimeIntegratorSetup::Concretize<time_integrator>::type;
using MultiPhaseManagerConcretization = MultiPhaseManagerSetup::Concretize<phase_manager>::type;
using PrimeStateHandlerConcretization = PrimeStateHandlerSetup::Concretize<prime_state_handler>::type;

/**
 *  @brief Executes the multiresolution Algorithm according to Kaiser et al. (to appear).
 * Allows different temporal and spatial solvers.
 *  @note This class is the only object which may invoke state changes in Tree [and its nodes, and their blocks] (local data change) or in the
 *        TopologyManager (global data change). The correct order is to execute all local changes first and then propagate the changes to the
 *        TopologyManager.
 */
class ModularAlgorithmAssembler {

   SimulationSetup const& setup_;
   MaterialManager const& material_manager_;

   TimeIntegratorConcretization time_integrator_;

   InputOutputManager& input_output_;

   Tree& tree_;
   TopologyManager& topology_;
   HaloManager& halo_manager_;
   CommunicationManager& communicator_;

   Multiresolution const& multiresolution_;
   Averager averager_;

   MultiPhaseManagerConcretization multi_phase_manager_;  //TODO-19 NH make const

   PrimeStateHandlerConcretization const prime_state_handler_;

   SpaceSolver const space_solver_;

   LogWriter& logger_;

   void CreateNewSimulation();
   void RestartSimulation();

   void Advance();
   void ProvideDebugInformation(const std::string debug_string, const bool plot_this_step, const bool print_this_step, unsigned int& debug_key) const;
   void LogElapsedTimeSinceInProfileRuns( double const start_time, std::string const message );

   void ComputeRightHandSide(const std::vector<unsigned int> levels,const unsigned int stage);
   void SwapBuffers(const std::vector<unsigned int> updated_levels, const unsigned int stage) const;
   void Integrate(const std::vector<unsigned int> updated_levels, const unsigned int stage);
   void JumpFluxAdjustment(const std::vector<unsigned int> finished_levels_descending) const;

   double ComputeTimestepSize() const;

   void ResetAllJumpBuffers() const;
   void ResetJumpConservativeBuffers(const std::vector<unsigned int> levels) const;

   void LoadBalancing(const std::vector<unsigned int> updated_levels, const bool force = false);

   void ImposeInitialCondition(const unsigned int level);

   void UpdateInterfaceTags(const std::vector<unsigned int> levels_with_updated_parents_descending) const;
   void SenseApproachingInterface(const std::vector<unsigned int> levels_ascending, bool refine_if_necessary = true);
   void SenseVanishedInterface(const std::vector<unsigned int> levels_descending);

   void Remesh(const std::vector<unsigned int> levels_to_update_ascending);
   void DetermineRemeshingNodes( std::vector<unsigned int> const parent_levels, std::vector<std::uint64_t>& remove_list,
                                 std::vector<std::uint64_t>& refine_list ) const;

   void RefineNode(const std::uint64_t node_id);

   void UpdateTopology();

   std::vector<unsigned int> GetLevels(const unsigned int timestep) const;

   template<ConservativeBufferType C>
   void ObtainPrimeStatesFromConservatives(const std::vector<unsigned int> updated_levels, const bool skip_levelset_nodes = false) const;

   template<ConservativeBufferType C>
   void DoObtainPrimeStatesFromConservativesForNonLevelsetNodes(Node& node) const;
   template<ConservativeBufferType C>
   void DoObtainPrimeStatesFromConservativesForLevelsetNodes(Node& node) const;

   void LogNodeNumbers() const;

public:
   ModularAlgorithmAssembler() = delete;
   explicit ModularAlgorithmAssembler( Tree& flower, TopologyManager& topology, HaloManager& halo_manager, CommunicationManager& communication,
                                       Multiresolution const& multiresolution, MaterialManager const& material_manager,
                                       SimulationSetup const& setup, InputOutputManager& io );
   ~ModularAlgorithmAssembler() = default;
   ModularAlgorithmAssembler( ModularAlgorithmAssembler const& ) = delete;
   ModularAlgorithmAssembler& operator=( ModularAlgorithmAssembler const& ) = delete;
   ModularAlgorithmAssembler( ModularAlgorithmAssembler&& ) = delete;
   ModularAlgorithmAssembler& operator=( ModularAlgorithmAssembler&& ) = delete;

   void ComputeLoop();

   void Initialization();
};

#endif // MODULAR_ALGORITHM_ASSEMBLER_H