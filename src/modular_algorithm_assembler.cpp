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
#include "modular_algorithm_assembler.h"

#include <bitset>
#include <algorithm>
#include <string>

#include "user_specifications/compile_time_constants.h"
#include "user_specifications/debug_and_profile_setup.h"
#include "user_specifications/riemann_solver_settings.h"
#include "levelset_block.h"
#include "interface_tags/interface_tag_functions.h"
#include "enums/interface_tag_definition.h"
#include "enums/remesh_identifier.h"
#include "enums/restore_mode.h"
#include "topology/id_information.h"
#include "helper_functions.h"
#include "communication/mpi_utilities.h"
#include "multiresolution/multiresolution.h"

namespace {
   void SetTimeInProfileRuns( double& time ) {
      if constexpr( DP::Profile() ) {
         time = MPI_Wtime();
      }
   }
}

/**
 * @brief Default constructor. Creates an instance to advance the simulation with the specified spatial solver and temporal integrator.
 * @param flower The tree instance which gives access to the fluid data.
 * @param topology TopologyManager instance for access of global information.
 * @param communication CommunicationManager instance to update other ranks of local changes.
 * @param material_manager Instance of a material manager, which already has been initialized according to the user input.
 * @param setup Instance to get access to user defined properties, relevant for the simulation.
 * @param input_output Instance of the I/O manager that handles filesystem access and ouput decisions.
 */
ModularAlgorithmAssembler::ModularAlgorithmAssembler( Tree& flower, TopologyManager& topology, HaloManager& halo_manager,
                                                      CommunicationManager& communication, Multiresolution const& multiresolution,
                                                      MaterialManager const& material_manager, SimulationSetup const& setup,
                                                      InputOutputManager& input_output ) :
   setup_( setup ),
   material_manager_( material_manager ),
   time_integrator_( setup_.GetStartTime() ),
   input_output_( input_output ),
   tree_( flower ),
   topology_( topology ),
   halo_manager_( halo_manager ),
   communicator_( communication ),
   multiresolution_( multiresolution ),
   averager_( topology_, communicator_, tree_ ),
   multi_phase_manager_( setup_, material_manager_, halo_manager_ ),
   prime_state_handler_( material_manager ),
   space_solver_( material_manager_, setup_.GetGravity() ),
   logger_( LogWriter::Instance() )
{
   /* Empty besides initializer list*/
}

/**
 * @brief Executes the outermost time loop. I.e. advances the simulation for n macro timesteps,
 *        according to the algorithm of Kaiser et al. (to appear). Also creates outputs if desired and
 *        balances the load between MPI ranks. Information to the outside, e.g. is most meaningfully
 *        created between two macro time step within this function.
 */
void ModularAlgorithmAssembler::ComputeLoop() {

   std::vector<double> loop_times;
   double time_measurement_start;
   double time_measurement_end;
   double current_simulation_time = time_integrator_.CurrentRunTime();

   bool time_step_size_too_small = false;

   /* fast forward if current time is already greater than start time (i.e. simulation was restarted).
    * We have to catch divison by zero in case only the initalization should be performed, i. e. start and end time are set to zero
    */
   double const flush_percentage = setup_.GetEndTime() == setup_.GetStartTime() ? 0.0 : ( current_simulation_time - setup_.GetStartTime() ) / ( setup_.GetEndTime() - setup_.GetStartTime() );
   logger_.FlushAlpaca( flush_percentage , current_simulation_time > setup_.GetStartTime() );

   while(current_simulation_time < setup_.GetEndTime() && ( !time_step_size_too_small )) {
      MPI_Barrier(MPI_COMM_WORLD); //For Time measurement
      time_measurement_start = MPI_Wtime();
      Advance(); // This is the heart of the Simulation, the advancement in Time over the different levels
      ResetAllJumpBuffers();
      MPI_Barrier(MPI_COMM_WORLD); // For Time measurement
      time_measurement_end = MPI_Wtime();
      loop_times.push_back( time_measurement_end - time_measurement_start );
      // Information Logging
      LogNodeNumbers();
      if constexpr(CC::WTL()) {
         input_output_.WriteTimestepFile(time_integrator_.MicroTimestepSizes());
      }
      if( time_integrator_.MicroTimestepSizes().back() < CC::MTS() ) {
         time_step_size_too_small = true;
      }
      time_integrator_.FinishMacroTimestep();
      current_simulation_time = time_integrator_.CurrentRunTime();
      logger_.LogMessage( "Macro timestep done t = " + ToScientificNotationString( setup_.DimensionalizeTime( current_simulation_time ), 9 ) );
      logger_.FlushAlpaca((current_simulation_time - setup_.GetStartTime()) / (setup_.GetEndTime() - setup_.GetStartTime()));
      // writing a restart file has priority over normal output, so call it first
      input_output_.WriteRestartFile(current_simulation_time, time_step_size_too_small);
      if( input_output_.WriteFullOutput(current_simulation_time, time_step_size_too_small) ) {
         // if output has been written this timestep, we also write profiling information
         if constexpr( DP::Profile() ) {logger_.LogMessage( topology_.LeafRankDistribution( MpiUtilities::NumberOfRanks() ) );}
      }
   }

   if constexpr( DP::Profile() ) {
      logger_.LogMessage( SummedCommunicationStatisticsString() );
   }
   logger_.LogMessage(" Total Time Spent in Compute Loop (seconds): " + ToScientificNotationString( std::accumulate( loop_times.begin(), loop_times.end(), 0.0 ), 5 ) );
}

/**
 * @brief Sets-up the starting point of the simulation based on either a restart file or the initial conditions depending on the configuration.
 */
void ModularAlgorithmAssembler::Initialization() {
   switch(setup_.GetRestoreMode()) {
      case RestoreMode::Disabled:
         logger_.LogMessage("Restore disabled: initializing new simulation");
         CreateNewSimulation();
         break;
      case RestoreMode::SoftRestore:
         if( input_output_.CheckIfRestoreFileExists() ) {
            logger_.LogMessage("Soft restore enabled: initializing simulation from restart file");
            RestartSimulation();
         } else {
            logger_.LogMessage("WARNING: restart file specified in the input file does not exist!");
            logger_.LogMessage("Soft restore enabled: initializing new simulation");
            CreateNewSimulation();
         }
         break;
      case RestoreMode::ForcedRestore:
         logger_.LogMessage("Forced restore enabled: initializing simulation from restart file");
         RestartSimulation();
         break;
      default:
         // as this method is only called once at startup, we don't need to make use of the PERFORMANCE flag
         throw std::invalid_argument("This restore mode is not known!");
   }
   logger_.LogMessage("Simulation successfully initialized");

   // Information Logging
   LogNodeNumbers();

   // initial output
   double const run_time = time_integrator_.CurrentRunTime();
   input_output_.WriteFullOutput(run_time, true);
   input_output_.WriteRestartFile(run_time); // Does only trigger writing of restart file if requested by user input (=inputfile).
}

/**
 * @brief Sets-up the starting point of the simulation based on a restart file. I.e. restores topology and conservatives from a previous run of the simulation.
 */
void ModularAlgorithmAssembler::RestartSimulation() {

   double const restart_time = input_output_.RestoreSimulationFromSnapshot();

   // align the simulation with the restart time
   time_integrator_.SetStartTime(restart_time);

   // project interface tags to make sure that cut cell tags are present on lower levels
   std::vector<unsigned int> child_levels_descending( setup_.AllLevels() );
   std::reverse(child_levels_descending.begin(), child_levels_descending.end());
   child_levels_descending.pop_back(); // remove level 0
   UpdateInterfaceTags( child_levels_descending );

   // initialize volume fractions if necessary
   std::vector<std::reference_wrapper<Node>> const nodes_needing_multiphase_treatment = tree_.NodesWithLevelset();
   bool globaly_existing_multi_phase_nodes = MpiUtilities::GloballyReducedBool( !nodes_needing_multiphase_treatment.empty() );
   if( globaly_existing_multi_phase_nodes ) {
      multi_phase_manager_.InitializeVolumeFractionBuffer(nodes_needing_multiphase_treatment);
      multi_phase_manager_.ObtainInterfaceQuantities( nodes_needing_multiphase_treatment );
   }

   // swap buffers
   SwapBuffers( setup_.AllLevels(), 0 );
}

/**
 * @brief Sets-up the starting point of the simulation based on the input file. I.e. inserts the initial conditions and creates the topology.
 * @note In the process the grid is build-up levelwise from coarse to fine. The (analytic) initial condition given by the user is, therein, merely
 * sampled at the cell center coordinates. In the subsequent multi-resolution analysis this can - commonly only for pathological cases - lead to slight
 * inacuracies and might in a worst-case situation alter the mesh-structure compared to (memory-wise infeasable) top-down approach.
 */
void ModularAlgorithmAssembler::CreateNewSimulation() {

   double levelset_temp[CC::TCX()][CC::TCY()][CC::TCZ()];
   std::int8_t initial_interface_tags[CC::TCX()][CC::TCY()][CC::TCZ()];
   std::vector<MaterialName> initial_materials;
   std::vector<std::uint64_t> coarsable_list;
   std::vector<std::uint64_t> globally_coarsable;
   std::vector<std::uint64_t> parents_of_coarsable;
   std::vector<std::uint64_t> refinement_list; //This list is only need as stub in this function.

   int const my_rank = communicator_.MyRankId();
   for( unsigned int level = 0; level <= setup_.GetMaximumLevel(); ++level ) {
      if( level > 0 ) {
         for( const std::uint64_t& node_id : topology_.IdsOnLevelOfRank( level - 1, my_rank ) ) { //We refine the parent to get the level we want to work on
            topology_.RefineNodeWithId( node_id );
         }
         UpdateTopology();
      }
      for( const std::uint64_t& node_id : topology_.IdsOnLevelOfRank( level, my_rank ) ) {
         std::uint64_t parent_id = ParentIdOfNode(node_id);
         if( level == 0 || topology_.IsNodeMultiPhase( parent_id ) ) { //Evaluated left to right, makes it safe on level zero
            // get the materials that initially exists in this node (this is determined on the finest level to avoid losing structures that are unresolved by this node)
            initial_materials = setup_.GetInitialMaterials( node_id );
            if( initial_materials.size() > 1 ) {
               // we have a multi node, thus we need the correct levelset to determine the interface tags
               setup_.GetInitialLevelset(node_id, levelset_temp);
               InterfaceTagFunctions::InitializeInternalInterfaceTags( initial_interface_tags );
               InterfaceTagFunctions::SetInternalCutCellTagsFromLevelset(levelset_temp, initial_interface_tags);
            } else {
               // we have a single node, thus we need only to consider uniform interface tags but no levelset
               const std::int8_t uniform_tag = MaterialSignCapsule::SignOfMaterial( initial_materials.front() ) * ITTI( IT::BulkPhase );
               for( unsigned int i = CC::FICX(); i <= CC::LICX(); ++i ) {
                  for( unsigned int j = CC::FICY(); j <= CC::LICY(); ++j ) {
                     for( unsigned int k = CC::FICZ(); k <= CC::LICZ(); ++k ) {
                        initial_interface_tags[i][j][k] = uniform_tag;
                     }
                  }
               }
            }
            if( level == setup_.GetMaximumLevel() && initial_materials.size() > 1 ) { // multi and leaf
               tree_.CreateNode(node_id, initial_materials, initial_interface_tags, std::make_unique<LevelsetBlock>(levelset_temp));
            } else { // either single or non-leaf multi
               tree_.CreateNode(node_id, initial_materials, initial_interface_tags);
            }
         } else { // parent is single fluid
            initial_materials = topology_.GetFluidsOfNode(parent_id);
            // copying tags of parent is sufficient as they are the same (single fluid)
            tree_.CreateNode(node_id, initial_materials, tree_.GetNodeWithId( parent_id ).GetInterfaceTags());
         }
         for( const MaterialName& material : initial_materials ) {
            topology_.AddFluidToNode(node_id, material);
         }
      }

      // updates weights and fluids in the topology
      UpdateTopology();
      if( level == setup_.GetMaximumLevel() ) {
         halo_manager_.LevelsetHaloUpdateOnLmax( LevelsetBlockBufferType::PhiRightHandSide );
         for( const std::uint64_t& node_id : topology_.IdsOnLevelOfRank( level, my_rank ) ) {
            Node& node = tree_.GetNodeWithId( node_id );
            if( node.HasLevelset() ) {
               InterfaceTagFunctions::SetInternalCutCellTagsFromLevelset(node.GetLevelsetBlock().GetPhiReinitialized(), node.GetInterfaceTags());
            }
         }
      }
      halo_manager_.InterfaceTagHaloUpdateOnLevelList( { level } );
      for( const std::uint64_t& node_id : topology_.IdsOnLevelOfRank( level, my_rank ) ) {
         InterfaceTagFunctions::SetTotalInterfaceTagsFromCutCells( tree_.GetNodeWithId( node_id ).GetInterfaceTags() );
      }
      halo_manager_.InterfaceTagHaloUpdateOnLevelList( { level } );
      SenseApproachingInterface( {level}, false ); // setting refine=false might cause an ill-defined tree afterwards which is corrected in next loop (++level)
      ImposeInitialCondition(level);
      halo_manager_.FluidHaloUpdateOnLevel( level, FluidFieldType::Conservatives );
      if( level == setup_.GetMaximumLevel() ) {
         halo_manager_.LevelsetHaloUpdateOnLmax( LevelsetBlockBufferType::PhiRightHandSide );
      }
      if( level > 1 ) { // Level One may never be coarsened.
         coarsable_list.clear();
         DetermineRemeshingNodes( {level - 1}, coarsable_list, refinement_list ); //Called on parent level

         MpiUtilities::LocalToGlobalData( coarsable_list, MPI_LONG_LONG_INT, MpiUtilities::NumberOfRanks(), globally_coarsable );

         // We get the parents, as coarsening is called on parents. Duplicates do not hurt (TopologyManager handles them).
         parents_of_coarsable = globally_coarsable;
         std::for_each(parents_of_coarsable.begin(), parents_of_coarsable.end(), [](std::uint64_t& to_parent) { to_parent = ParentIdOfNode(to_parent); });
         for( const auto& coarsable_id : coarsable_list ) {
            if( topology_.NodeIsOnRank( coarsable_id, communicator_.MyRankId() ) ) {
               tree_.RemoveNodeWithId( coarsable_id );
            }
         }

         // Update the topology
         for( const std::uint64_t parent_id : parents_of_coarsable ) {
            topology_.CoarseNodeWithId( parent_id );
         }
         if( !parents_of_coarsable.empty() ) {
            communicator_.InvalidateCache();
         }
      }

      /** The swaps around the LoadBalancing are necessary, since after initialization the correct level-set values are in the right-hand side level-set buffer,
       *  but LoadBalancing has to work on the reinitialized level-set buffer. */
      for( Node& levelset_node : tree_.NodesWithLevelset() ) {
         std::swap(levelset_node.GetLevelsetBlock().GetPhiReinitialized(), levelset_node.GetLevelsetBlock().GetPhiRightHandSide());
      }
      LoadBalancing( setup_.AllLevels(), true );
   } // loop: level

   UpdateTopology();

   // project interface tags to make sure that cut cell tags are present on lower levels which underresolve the initial levelset
   std::vector<unsigned int> child_levels_descending( setup_.AllLevels() );
   std::reverse(child_levels_descending.begin(), child_levels_descending.end());
   child_levels_descending.pop_back(); // remove level 0

   /** At this place, the prescribed initial conditions are imposed.
    *  The prescribed level-set field is in phi_rhs
    *  The prescribed initial conditions are in the conservative_rhs buffer
    *  In case the prescribed level-set field is not a perfect distance function, reinitialization has to be done.
    *  Based on the reinitialized level-set buffer, the volume fractions have to be calculated.
    *  Based on the reinitialized level-set field and the volume fractions, the extension can be done. */
   std::vector<std::reference_wrapper<Node>> const nodes_needing_multiphase_treatment = tree_.NodesWithLevelset();
   bool const exist_multi_nodes_global = MpiUtilities::GloballyReducedBool( !nodes_needing_multiphase_treatment.empty() );
   if( exist_multi_nodes_global ) {
      for( const std::uint64_t& node_id : topology_.IdsOnLevelOfRank( setup_.GetMaximumLevel(), my_rank ) ) {
         for( auto& phase : tree_.GetNodeWithId( node_id ).GetPhases() ) {
            for( const Equation eq : FF::ASOE() ) {
               /**
                * During the initialization procedure conservative values are initialized in the conservative_rhs buffer.
                * In the following, prime states have to be calculated based on the initialized conservatives. Prime state calculation requires
                * conservative values in the conservative_avg buffer. Thus, a swap between the conservative_rhs and the conservative_avg buffer is necessary.
                */
               std::swap(phase.second.GetRightHandSideBuffer(eq), phase.second.GetAverageBuffer(eq));
            } // equation
         } // phases
      } // nodes

      ObtainPrimeStatesFromConservatives<ConservativeBufferType::Average>( {setup_.GetMaximumLevel()} );
      multi_phase_manager_.EnforceWellResolvedDistanceFunction(nodes_needing_multiphase_treatment, 0, true);
      multi_phase_manager_.InitializeVolumeFractionBuffer(nodes_needing_multiphase_treatment);
      UpdateInterfaceTags( child_levels_descending );
      multi_phase_manager_.CopyReinitializedLevelsetToRightHandSideLevelset(nodes_needing_multiphase_treatment);
      ObtainPrimeStatesFromConservatives<ConservativeBufferType::Average>( {setup_.GetMaximumLevel()} );
      for( const std::uint64_t& node_id : topology_.IdsOnLevelOfRank( setup_.GetMaximumLevel(), my_rank ) ) {
         for( auto& phase : tree_.GetNodeWithId( node_id ).GetPhases() ) {
            for( const Equation eq : FF::ASOE() ) {
               std::swap(phase.second.GetRightHandSideBuffer(eq), phase.second.GetAverageBuffer(eq));
            } // equation
         } // phases
      } // nodes

      multi_phase_manager_.Extend( nodes_needing_multiphase_treatment, 0 );
   } //if multi nodes exist

   // So far everything happened in the RHS buffer - now we swap for proper start.
   SwapBuffers( setup_.AllLevels(), 0 );

   // Calculate the prime states based on the initialized conservatives and save them in the prime state buffer.
   ObtainPrimeStatesFromConservatives<ConservativeBufferType::Average>(setup_.AllLevels());
   if( exist_multi_nodes_global ) {
      multi_phase_manager_.ObtainInterfaceQuantities( nodes_needing_multiphase_treatment );
      if constexpr(GeneralTwoPhaseSettings::LogConvergenceInformation) {
         logger_.DelayedLogMessage(true, true);
      }
   }
}

/**
 * @brief Advances the simulation in time by one macro time step. To do so 2^level micro time steps need to be executed
 *        each consisting of 1 to n stages, depending on the used time integrator.
 */
void ModularAlgorithmAssembler::Advance() {

   // These variables are only for debugging. Allow fine tuned debug output. see below.
   unsigned int debug_key = 0;
   bool plot_this_step = false;
   bool print_this_step = false;

   // These variables are only for profiling
   double function_timer = 0.0;

   unsigned int const maximum_level = setup_.GetMaximumLevel();

   //number of timesteps to run on the maximum level to run one timestep on level 0
   unsigned int const number_of_timesteps_on_finest_level = 1 << (maximum_level);

   // In the first run all levels need to update
   std::vector<unsigned int> levels_to_update_descending( setup_.AllLevels() );
   std::reverse(levels_to_update_descending.begin(), levels_to_update_descending.end());
   std::vector<unsigned int> levels_to_update_ascending;
   std::vector<unsigned int> levels_with_updated_parents_descending;
   std::vector<std::reference_wrapper<Node>> nodes_needing_multiphase_treatment = tree_.NodesWithLevelset();
   bool exist_multi_nodes_global = MpiUtilities::GloballyReducedBool( !nodes_needing_multiphase_treatment.empty() );

   //Run enough timesteps on the maximum level to run one timestep on level 0
   for(unsigned int timestep = 0; timestep < number_of_timesteps_on_finest_level; ++timestep) {

      SetTimeInProfileRuns( function_timer );
      time_integrator_.AppendMicroTimestep(ComputeTimestepSize());
      LogElapsedTimeSinceInProfileRuns( function_timer, "ComputeTimestepSize                " );
      ProvideDebugInformation("ComputeTimestepSize - Done ", plot_this_step, print_this_step, debug_key);

      for(unsigned int stage = 0; stage < time_integrator_.NumberOfStages(); ++stage) {

         if constexpr( DP::Debug() ) {
            debug_key = 1000 * timestep + 100 * stage;
            // Example how to enable fine grained plotting.
            /*if(time_integrator_.CurrentRunTime() > 0.0 && time_integrator_.CurrentRunTime() < 0.0001) {
                plot_this_step = true;
                print_this_step = true;
            } else {
                plot_this_step = false;
                print_this_step = false;
            }*/
         }
         ProvideDebugInformation("Start of Loop ", false, print_this_step, debug_key);

         // compute rhs on all levels which need to be updated this integer timestep
         SetTimeInProfileRuns( function_timer );
         ComputeRightHandSide(levels_to_update_descending, stage);
         LogElapsedTimeSinceInProfileRuns( function_timer, "ComputeRightHandSide               " );
         ProvideDebugInformation("ComputeRightHandSide - Done ", plot_this_step, print_this_step, debug_key);


         //Flux averaging from levels which run this timestep down to the lowest neighbor level or parent
         SetTimeInProfileRuns( function_timer );
         averager_.AverageFluid( levels_to_update_descending );
         LogElapsedTimeSinceInProfileRuns( function_timer, "AverageFluid                       " );
         ProvideDebugInformation("AverageFluid - Done ", plot_this_step, print_this_step, debug_key);

         SetTimeInProfileRuns( function_timer );
         halo_manager_.FluidHaloUpdate(setup_.AllLevels(), FluidFieldType::Conservatives);
         LogElapsedTimeSinceInProfileRuns( function_timer, "UpdateHalos (all)                  " );
         ProvideDebugInformation( "UpdateHalos( AllLevels ) - Done ", plot_this_step, print_this_step, debug_key );

         if(exist_multi_nodes_global) {
            SetTimeInProfileRuns( function_timer );
            halo_manager_.LevelsetHaloUpdateOnLmax(LevelsetBlockBufferType::PhiRightHandSide);
            LogElapsedTimeSinceInProfileRuns( function_timer, "LevelsetHaloUpdate                 " );
            ProvideDebugInformation("LevelsetHaloUpdate (maximum level) - Done ", plot_this_step, print_this_step, debug_key);
         }

         //Get which levels need to be advanced from here on in this timestep and stage
         levels_to_update_descending = GetLevels( timestep );
         levels_to_update_ascending.clear();
         std::reverse_copy(levels_to_update_descending.begin(), levels_to_update_descending.end(), std::back_inserter(levels_to_update_ascending));
         levels_with_updated_parents_descending = levels_to_update_descending;
         levels_with_updated_parents_descending.pop_back();

         SetTimeInProfileRuns( function_timer );
         Integrate(levels_to_update_descending, stage);
         LogElapsedTimeSinceInProfileRuns( function_timer, "Integrate                          " );
         ProvideDebugInformation( "Integration - Done ", plot_this_step, print_this_step, debug_key );

         if(exist_multi_nodes_global) {
            SetTimeInProfileRuns( function_timer );
            multi_phase_manager_.PropagateLevelset(nodes_needing_multiphase_treatment, stage);
            LogElapsedTimeSinceInProfileRuns( function_timer, "PropagateLevelset                  " );
            ProvideDebugInformation("PropagateLevelset in MultiphaseManager - Done ", plot_this_step, print_this_step, debug_key);
         }

         //Averaging of new mean values - project all levels to be on the safe side
         SetTimeInProfileRuns( function_timer );
         averager_.AverageFluid( levels_with_updated_parents_descending );
         LogElapsedTimeSinceInProfileRuns( function_timer, "ProjectLevel                       " );
         ProvideDebugInformation("ProjectLevel - Done ", plot_this_step, print_this_step, debug_key);

         //to maintain conservation
         if( time_integrator_.IsLastStage(stage) ) {
            // We correct the values at jumps to maintain conservation.
            SetTimeInProfileRuns( function_timer );
            JumpFluxAdjustment(levels_to_update_descending);
            LogElapsedTimeSinceInProfileRuns( function_timer, "AdjustJumpFluxes                   " );
            ProvideDebugInformation("AdjustJumpFluxes - Done ", plot_this_step, print_this_step, debug_key);
         }

         //boundary exchange mean values and jumps on finished levels
         SetTimeInProfileRuns( function_timer );
         halo_manager_.FluidHaloUpdate(levels_to_update_ascending, FluidFieldType::Conservatives, true);
         LogElapsedTimeSinceInProfileRuns( function_timer, "UpdateHalos (cut_jumps)            " );
         ProvideDebugInformation("UpdateHalos( levels_to_update, cut_jump=true ) - Done ", plot_this_step, print_this_step, debug_key );

         if(time_integrator_.IsLastStage(stage)) {
            if(exist_multi_nodes_global) {
               SetTimeInProfileRuns( function_timer );
               SenseVanishedInterface(levels_to_update_descending);
               LogElapsedTimeSinceInProfileRuns( function_timer, "SenseVanishedInterface             " );
               ProvideDebugInformation("SenseVanishedInterface - Done ", plot_this_step, print_this_step, debug_key);
            }

            SetTimeInProfileRuns( function_timer );
            Remesh(levels_to_update_ascending);
            LogElapsedTimeSinceInProfileRuns( function_timer, "Remesh                             " );
            ProvideDebugInformation("Remesh - Done ", plot_this_step, print_this_step, debug_key);

            if(exist_multi_nodes_global) { // TODO-19 JW and NH: Think of removing this if statement in case of a newly created interface (phase change...)
               SetTimeInProfileRuns( function_timer );
               SenseApproachingInterface(levels_to_update_ascending);
               LogElapsedTimeSinceInProfileRuns( function_timer, "SenseApproachingInterface          " );
               ProvideDebugInformation("SenseApproachingInterface - Done ", plot_this_step, print_this_step, debug_key);
            }

            SetTimeInProfileRuns( function_timer );
            LoadBalancing( levels_to_update_descending );
            LogElapsedTimeSinceInProfileRuns( function_timer, "LoadBalancing                      " );

            nodes_needing_multiphase_treatment = tree_.NodesWithLevelset();
            exist_multi_nodes_global = MpiUtilities::GloballyReducedBool( !nodes_needing_multiphase_treatment.empty() );
            ProvideDebugInformation("LoadBalancing - Done ", plot_this_step, print_this_step, debug_key);

            if constexpr( DP::Profile() ) {
               logger_.LogMessage( "Node count after remeshing:" );
               LogNodeNumbers();
            }
         } //last stage

         if(exist_multi_nodes_global) {
            SetTimeInProfileRuns( function_timer );
            multi_phase_manager_.Mix(nodes_needing_multiphase_treatment, stage);
            LogElapsedTimeSinceInProfileRuns( function_timer, "Mixing                             " );
            ProvideDebugInformation("Mixing - Done ", plot_this_step, print_this_step, debug_key);

            bool const is_last_stage = time_integrator_.IsLastStage(stage);
            SetTimeInProfileRuns( function_timer );
            multi_phase_manager_.EnforceWellResolvedDistanceFunction(nodes_needing_multiphase_treatment, stage, is_last_stage );
            LogElapsedTimeSinceInProfileRuns( function_timer, "EnforceWellResolvedDistanceFunction" );
            std::string&& message = is_last_stage ? "EnforceWellResolvedDistanceFunction (possibly with scale separation) - Done " : "EnforceWellResolvedDistanceFunction - Done ";
            ProvideDebugInformation( message, plot_this_step, print_this_step, debug_key );

            SetTimeInProfileRuns( function_timer );
            UpdateInterfaceTags( levels_with_updated_parents_descending );
            LogElapsedTimeSinceInProfileRuns( function_timer, "UpdateInterfaceTags                " );
            ProvideDebugInformation("UpdateInterfaceTags - Done ", plot_this_step, print_this_step, debug_key);

            SetTimeInProfileRuns( function_timer );
            ObtainPrimeStatesFromConservatives<ConservativeBufferType::RightHandSide>( {setup_.GetMaximumLevel()}, true );
            LogElapsedTimeSinceInProfileRuns( function_timer, "ObtainPrimeStatesFromConservatives " );
            ProvideDebugInformation("ObtainPrimeStatesFromConservatives - Done ", plot_this_step, print_this_step, debug_key);

            SetTimeInProfileRuns( function_timer );
            multi_phase_manager_.Extend( nodes_needing_multiphase_treatment, stage );
            LogElapsedTimeSinceInProfileRuns( function_timer, "Extend                             " );
            ProvideDebugInformation("Extend - Done ", plot_this_step, print_this_step, debug_key);
         }

         //SWAP on levels which were integrated this step
         SetTimeInProfileRuns( function_timer );
         SwapBuffers( levels_to_update_descending, stage );
         LogElapsedTimeSinceInProfileRuns( function_timer, "Swap                               " );
         ProvideDebugInformation("SwapOnLevel - Done ", plot_this_step, print_this_step, debug_key);

         // Calculate the prime states based on the integrated conservatives and save them in the prime state buffer.
         SetTimeInProfileRuns( function_timer );
         ObtainPrimeStatesFromConservatives<ConservativeBufferType::Average>(setup_.AllLevels()); //TODO-19 JW and NH: It would suffice to do it on updated_levels if we could track nodes that were load balanced but not integrated.
         LogElapsedTimeSinceInProfileRuns( function_timer, "ObtainPrimeStatesFromConservatives " );
         ProvideDebugInformation("ObtainPrimeStatesFromConservatives - Done ", plot_this_step, print_this_step, debug_key);

         if(exist_multi_nodes_global) {
            SetTimeInProfileRuns( function_timer );
            multi_phase_manager_.ObtainInterfaceQuantities( nodes_needing_multiphase_treatment, time_integrator_.IsLastStage(stage) );
            LogElapsedTimeSinceInProfileRuns( function_timer, "SetInterfaceQuantities             " );
            ProvideDebugInformation("SetInterfaceQuantities - Done ", plot_this_step, print_this_step, debug_key);
            if constexpr(GeneralTwoPhaseSettings::LogConvergenceInformation) {
               logger_.DelayedLogMessage(true, true);
            }
         }
      } // stages

      // Safe stop of the code (after every micro timestep, therefore not in ComputeLoop)
      if( input_output_.CheckIfAbortfileExists() ) {
         logger_.LogMessage("The file 'ABORTFILE' was found in the output folder. Simulation is being terminated");
         throw std::runtime_error("The simulation was aborted by the user! \n");
      }
   }
   communicator_.ResetTagsForPartner();
}

/**
 * @brief Provides debug information for the current state. Depending on user defined compile time constants, information is written to the log-file or
 * debug output is written
 * @param debug_string The string which describes the state for which debug information is provided.
 * @param plot_this_step Decision whether an debug output file is written for the current state.
 * @param print_this_step Decision whether information is print to the terminal for the current state.
 * @param debug_key An integer to allow increasing numbering of the debug output files.
 */
void ModularAlgorithmAssembler::ProvideDebugInformation( std::string const debug_string, bool const plot_this_step, bool const print_this_step,
                                                         unsigned int& debug_key ) const {
   if constexpr( DP::DebugLog() ) {
      if(print_this_step) { logger_.LogMessage(debug_string + std::to_string( debug_key ));}
   }
   if constexpr( DP::DebugOutput() ) {
      if(plot_this_step) {input_output_.WriteDebugFile(debug_key);}
   }
   if constexpr( DP::Debug() ) {debug_key++;}
}

/**
 * @brief Computes the f(u) term in the Runge-Kutta function u^i = u^(i-1) + c_i * dt * f(u).
 *        Stores the result in the right hand side buffers. Computations only done in leaves.
 * @param levels For all leaves on these levels the right hand side will be computed.
 * @param stage The current Runge-Kutta stage.
 */
void ModularAlgorithmAssembler::ComputeRightHandSide( std::vector<unsigned int> const levels, unsigned int const stage ) {

   // Global Lax-Friedrich scheme
   if constexpr( RoeSolverSettings::flux_splitting_scheme == FluxSplitting::GlobalLaxFriedrichs ) {
      // In case of Global Lax Friedrichs Eigenvalues must be collected across ranks and blocks
      double max_eigenvalues[DTI(CC::DIM())][FF::ANOE()];
      // NH Initializing to zero necessary!
      for( unsigned int d = 0; d < DTI(CC::DIM()); ++d ) {
         for( unsigned int e = 0; e < FF::ANOE(); ++e ) {
            max_eigenvalues[d][e] = 0.0;
         }
      }
      double current_eigenvalues[DTI(CC::DIM())][FF::ANOE()];
      for( const auto& level : levels ) {
         for( const Node& node : tree_.LeavesOnLevel(level) ) {
            for( auto& phase : node.GetPhases() ) {
               space_solver_.ComputeMaxEigenvaluesForPhase(phase, current_eigenvalues);
               for( unsigned int d = 0; d < DTI(CC::DIM()); ++d ) {
                  for( unsigned int e = 0; e < FF::ANOE(); ++e ) {
                     max_eigenvalues[d][e] = std::max(max_eigenvalues[d][e], current_eigenvalues[d][e]);
                  }
               }
            } // phases
         } // node
      } // level
      MPI_Allreduce(MPI_IN_PLACE, max_eigenvalues, DTI(CC::DIM()) * FF::ANOE(), MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
      space_solver_.SetFluxFunctionGlobalEigenvalues( max_eigenvalues );
   }
   for( const auto& level : levels ) {
      for( Node& node : tree_.LeavesOnLevel(level) ) {
         time_integrator_.FillInitialBuffer(node, stage);

         // compute fluxes for levelset and fluids (including single phase and interface contributions!)
         space_solver_.UpdateFluxes( node );

         // Integration can only be performed on conservatives, but not on volume averaged conservatives.
         // Thus, we have to transform the volume averaged conservatives, which are currently saved in the average buffer, to conservatives.
         multi_phase_manager_.TransformToConservatives( node );

         // Conservative as well as levelset buffers are prepared for integration
         time_integrator_.PrepareBufferForIntegration(node, stage);
      } // node
   } // level
}

/**
 * @brief Swaps the content of the average and right hand side buffers of all nodes on the specified level. Thereby, it is ensured that
 * the average buffer of the conservatives and the base buffer of the level set contain values based on which the right-hand side values
 * for the next RK-(sub)step can be calculated. In the last RK-stage it is also necessary to copy the values of the reinitialized level-set buffer
 * to the right-hand side level-set buffer before the swap is done. Thereby, it is ensured, that for the first RK-stage the level-set advection
 * right-hand side is calculated based on a reinitialized and potentially scale separated level-set field.
 * @param updated_levels The buffers of the blocks of all nodes on these levels will be swapped.
 * @param stage The current stage of the RK scheme.
 */
void ModularAlgorithmAssembler::SwapBuffers( std::vector<unsigned int> const updated_levels, unsigned int const stage ) const {

   for(const auto& level : updated_levels) {
      for(Node& node : tree_.NodesOnLevel(level)) {
         if(time_integrator_.IsLastStage(stage) && node.HasLevelset()) {

            LevelsetBlock& levelset_block = node.GetLevelsetBlock();
            const double (&phi_reinitialized)[CC::TCX()][CC::TCY()][CC::TCZ()] = levelset_block.GetPhiReinitialized();
            double (&phi_rhs)[CC::TCX()][CC::TCY()][CC::TCZ()] = levelset_block.GetPhiRightHandSide();
            for( unsigned int i = 0; i < CC::TCX(); ++i ) {
               for( unsigned int j = 0; j < CC::TCY(); ++j ) {
                  for( unsigned int k = 0; k < CC::TCZ(); ++k ) {
                     phi_rhs[i][j][k] = phi_reinitialized[i][j][k];
                  } //k
               } //j
            } //i
         } //node with level set and final stage
         time_integrator_.SwapBuffersForNextStage(node);
      } //nodes
   } //levels
}

/**
 * @brief Performs one integration stage with the selected time integrator for all nodes on the specified level.
 * @param updated_levels The levels which should be integrated in time.
 * @param stage The current integrator stage (before the integration is executed).
 */
void ModularAlgorithmAssembler::Integrate( std::vector<unsigned int> const updated_levels, unsigned int const stage ) {

   for( const auto& level : updated_levels ) {
      unsigned int number_of_timesteps = 1 << ( setup_.GetMaximumLevel() - level ); // 2^x

      // We integrate all leaves
      for( Node& node : tree_.LeavesOnLevel( level ) ) {
         time_integrator_.IntegrateNode( node, stage, number_of_timesteps );
      }

      /* We need to integrate the values in jump halos. However, as two jump halos may overlap,
       * we integrate ALL halos values if the node has ANY jump. The wrongly integrated halo cells
       * are later overwritten by the halo update
       */

      // We integrate jump halos on all nodes
      for( auto& [id, node] : tree_.GetLevelContent(level) ) {
         for( const BoundaryLocation& location : CC::HBS() ) {
            if( topology_.FaceIsJump( id, location ) ) {
               std::array<int, 3> start_indices_halo = communicator_.GetStartIndicesHaloRecv( location );
               std::array<int, 3> halo_size = communicator_.GetHaloSize( location );
               time_integrator_.IntegrateJumpHalos( node, stage, number_of_timesteps, start_indices_halo, halo_size );
            }
         }
      }
   } // level
}

/**
 * @brief Updates the interface tags on all levels based on the cut cells on the finest level.
 * @param nodes The nodes on the finest level, for which the level-set field has to be updated.
 * @param levels_with_updated_parents_descending The child levels whose parents were updated in descending order.
 */
void ModularAlgorithmAssembler::UpdateInterfaceTags( std::vector<unsigned int> const levels_with_updated_parents_descending ) const {

   /**
    * Step 0: To project cut-cell tags to parent levels and to subsequently set the narrow-band tags based on then, two level lists are necessary.
    *         1: The child levels whose parents were updated in descending order. This list is used for the interface tag averaging. Note, that this
    *            list must not contain level zero! The list given as input can be used for that.
    *         2: The parent levels which obtained information about cut-cell tag location during the interface-tag averaging.
    *            This list is named parent_levels_with_projected_cut_cell_tags.
    */
   std::vector<unsigned int> parent_levels_with_projected_cut_cell_tags( levels_with_updated_parents_descending );
   if( !levels_with_updated_parents_descending.empty() ) {
      parent_levels_with_projected_cut_cell_tags.push_back( parent_levels_with_projected_cut_cell_tags.back() - 1 );
   }

   /**
    * Step 1: The information about where cut cells are located on the finest level has to be projected to levels whose parent levels were updated.
    *         After that, information about the location of cut cells is available for the finest level and levels whose parents were updated.
    */
   averager_.AverageInterfaceTags( levels_with_updated_parents_descending );
   //Halo Update does not hurt, and we need not only parents but also the child levels
   halo_manager_.InterfaceTagHaloUpdateOnLevelList( setup_.AllLevels() );

   /**
    * Step 2: For all levels where cut cells were newly set, the narrow-band tags are set based on the cut-cell tags.
    */
   for( const unsigned int level : parent_levels_with_projected_cut_cell_tags ) {
      for( const auto node_id : topology_.IdsOnLevelOfRank( level, communicator_.MyRankId() ) ) {
         InterfaceTagFunctions::SetTotalInterfaceTagsFromCutCells( tree_.GetNodeWithId( node_id ).GetInterfaceTags() );
      }
   }

   /**
    * Step 3: A halo update on all levels where narrow-band tags were newly set is necessary in the end.
    */
   halo_manager_.InterfaceTagHaloUpdateOnLevelList( parent_levels_with_projected_cut_cell_tags );
}

/**
 * @brief Senses if an interface is approaching a single-phase node, e.g. by having entered the node's halo, and if necessary
 *        performs the change from single- to multi-phase node. Has to be called on levels in ascending order to make sure
 *        that the change propagates correctly up the tree.
 * @param levels_ascending The levels on which all nodes should be checked in ascending order.
 * @param refine_if_necessary Gives whether a found node should also be refined. If false, it might cause an ill-defined tree.
 */
void ModularAlgorithmAssembler::SenseApproachingInterface( std::vector<unsigned int> const levels_ascending, bool const  refine_if_necessary ) {
   // WARNING: calling this function with refine=false might cause an ill-defined tree, used only during initialization

   // make sure this function is called from lowest (e.g. 0) to highest (e.g. Lmax) level of interest, i.e. in ascending order
   // this guarantees that the appearance of the interface propagates up the tree

   bool levelset_block_created = false;
   bool node_refined = false;
   for( const auto& level : levels_ascending ) {
      for( const std::uint64_t& node_id : topology_.IdsOnLevelOfRank( level, communicator_.MyRankId() ) ) {
         if( !topology_.IsNodeMultiPhase(node_id) ) {
            Node& node = tree_.GetNodeWithId( node_id );
            // TODO-19 TP test total cells, but probably halo is enough for standard case (no phase change in bulk) --> OPTIMIZE?
            if( !InterfaceTagFunctions::TotalInterfaceTagsAreUniform(node.GetInterfaceTags()) ) {
               // not uniform anymore, thus change to multi

               // get additional material
               MaterialName const material_new = (node.GetSinglePhaseMaterial() == MaterialSignCapsule::PositiveFluidMaterial()) ? MaterialSignCapsule::NegativeFluidMaterial()
                                                                                                                                 : MaterialSignCapsule::PositiveFluidMaterial();
               topology_.AddFluidToNode(node_id, material_new);
               if( level == setup_.GetMaximumLevel() ) {
                  // Lmax node with levelset block but without children
                  auto const sign = Signum(node.GetUniformInterfaceTag());
                  node.AddPhase(material_new);
                  node.SetLevelsetBlock( std::make_unique<LevelsetBlock>(sign * CC::LSCOF()) );
                  if( refine_if_necessary ) levelset_block_created = true;
               } else {
                  // node without levelset but with children (if it is a leaf until now)
                  if( refine_if_necessary && topology_.NodeIsLeaf(node_id) ) {
                     // new children have to be generated before the switch to multi because RefineNode works on single fluid nodes only
                     RefineNode(node_id);
                     node_refined = true;
                  }
                  // now change to multi-phase
                  node.AddPhase(material_new);
               }
            }
         }
      }
      // update topology after every level to make sure, that the change can propagate more than one level
      UpdateTopology();
   }
   MPI_Allreduce(MPI_IN_PLACE, &levelset_block_created, 1, MPI_CXX_BOOL, MPI_LOR, MPI_COMM_WORLD);
   MPI_Allreduce(MPI_IN_PLACE, &node_refined, 1, MPI_CXX_BOOL, MPI_LOR, MPI_COMM_WORLD);
   if( levelset_block_created ) {
      halo_manager_.LevelsetHaloUpdateOnLmax( LevelsetBlockBufferType::PhiReinitialized );
   }
   if( node_refined ) {
      std::vector<unsigned int> halo_levels( levels_ascending );
      halo_levels.erase(halo_levels.begin());

      /*** We need to fill the Halo Cells ***/
      halo_manager_.FluidHaloUpdate( halo_levels, FluidFieldType::Conservatives );
   }
}

/**
 * @brief Senses if the interface has completely vanished from multi-phase node and its halo and if necessary
 *        performs the change from multi- to single-phase node. Has to be called on levels in descending order
 *        to make sure that the change propagates correctly down the tree.
 * @param levels_descending The levels on which all nodes should be checked in descending order.
 */
void ModularAlgorithmAssembler::SenseVanishedInterface( std::vector<unsigned int> const levels_descending ) {
   // make sure this function is called from highest (e.g. Lmax) to lowest (e.g. 0) level of interest, i.e. in descending order
   // this guarantees that the disappearance of the interface propagates down the tree

   for( const auto& level : levels_descending ) {
      for( const auto& node_id : topology_.IdsOnLevelOfRank( level, communicator_.MyRankId() ) ) {
         if( topology_.IsNodeMultiPhase(node_id) ) {
            // a multi-phase non-Lmax node has to have children, therefore no additional check for existence of children necessary
            bool all_children_single = true;
            if( level < setup_.GetMaximumLevel() ) {
               for( const auto& child_id : IdsOfChildren(node_id) ) {
                  all_children_single = all_children_single && ( !topology_.IsNodeMultiPhase(child_id) );
               }
            }
            if( all_children_single ) {
               // if all children are single, this node might become single as well
               Node& node = tree_.GetNodeWithId( node_id );
               if( InterfaceTagFunctions::TotalInterfaceTagsAreUniform(node.GetInterfaceTags()) ) {
                  // make single again

                  // get the vanished material (can use an arbitrary interface tag since it's already clear that they are uniform)
                  MaterialName const material_old = (node.GetInterfaceTags()[CC::FICX()][CC::FICY()][CC::FICZ()] < 0) ? MaterialSignCapsule::PositiveFluidMaterial()
                                                                                                                      : MaterialSignCapsule::NegativeFluidMaterial();
                  topology_.RemoveFluidFromNode(node_id, material_old);

                  // remove old material and levelset block (by setting it to nullptr)
                  node.RemovePhase(material_old);
                  node.SetLevelsetBlock();
               }
            }
         }
      }
      // update topology after every level to make sure, that the change can propagate more than one level
      UpdateTopology();
   }
}

/**
 * @brief Ensures conservation at resolution jumps following \cite Roussel2003. Uses the (more precise) fluxes on the finer level to correct the fluxes on the coarser level.
 *        The adjustment consists of 3 steps:
 *        - Averaging jump buffers down one level
 *        - Exchanging jump buffers to leaves on the lower levels
 *        - Resetting the jump buffers
 * @param finished_levels_descending Levels which ran this time instance and thus have correct values in their jump buffers.
 */
void ModularAlgorithmAssembler::JumpFluxAdjustment( std::vector<unsigned int> const finished_levels_descending ) const {

   // We have to cut level zero if present as it cannot send information downwards
   std::vector<unsigned int> levels_averaging_down(finished_levels_descending);
   levels_averaging_down.erase(std::remove(levels_averaging_down.begin(), levels_averaging_down.end(), 0), levels_averaging_down.end());

   // The finest level does not exchange
   std::vector<unsigned int> level_exchanging( finished_levels_descending );
   level_exchanging.erase(level_exchanging.begin());
   const int my_rank = communicator_.MyRankId();

   /*** Sending Down ***/
   // First the parents' jump buffers are filled form the childrens values.
   for( const auto& level : levels_averaging_down ) {
      for( const auto& child_id : topology_.GlobalIdsOnLevel(level) ) {
         const std::uint64_t parent_id = ParentIdOfNode(child_id);
         const int rank_of_child = topology_.GetRankOfNode(child_id);
         const int rank_of_parent = topology_.GetRankOfNode(parent_id);
         if( rank_of_child == my_rank && rank_of_parent == my_rank ) {
            // Non MPI Averaging
            Node& parent = tree_.GetNodeWithId( parent_id );
            const Node& child = tree_.GetNodeWithId( child_id );
            for(const auto material : topology_.GetFluidsOfNode(child_id)) {
               Multiresolution::AverageJumpBuffer( child.GetPhaseByMaterial( material ).GetBoundaryJumpConservatives(), parent.GetPhaseByMaterial( material ).GetBoundaryJumpConservatives(), child_id );
            }
         } else if( rank_of_child == my_rank && rank_of_parent != my_rank ) {
            // MPI_ISend
            const Node& child = tree_.GetNodeWithId( child_id );
            for( const auto material : topology_.GetFluidsOfNode(child_id) ) {
               MPI_Send( &child.GetPhaseByMaterial(material).GetBoundaryJumpConservatives(), CC::SIDES(), communicator_.JumpSurfaceDatatype(), rank_of_parent,
                        communicator_.TagForRank( rank_of_parent ), MPI_COMM_WORLD );
            }
         } else if( rank_of_child != my_rank && rank_of_parent == my_rank ) {
            SurfaceBuffer childs_jump_buffer;

            // MPI_Recv
            Node& parent = tree_.GetNodeWithId( parent_id );
            for(const auto material : topology_.GetFluidsOfNode(child_id)) {
               MPI_Recv(&childs_jump_buffer, CC::SIDES(), communicator_.JumpSurfaceDatatype(),rank_of_child, communicator_.TagForRank( rank_of_child ),MPI_COMM_WORLD,MPI_STATUS_IGNORE);
               Multiresolution::AverageJumpBuffer( childs_jump_buffer,parent.GetPhaseByMaterial( material ).GetBoundaryJumpConservatives(), child_id );
            }
         }
      }
   }

   /*** Exchanging ***/

   /* A node adjusts its values if:
    * 1.) It is leaf
    * 2.) Its neighbor exists
    * 3.) Its neighbor is not a leaf
    */

   std::uint64_t neighbor_id;
   bool neighbor_exists;
   bool neighbor_is_leaf;
   BoundaryLocation neighbor_location;
   std::vector<std::uint64_t> leaf_ids_on_level;
   unsigned int x_start;
   unsigned int x_end;
   unsigned int y_start;
   unsigned int y_end;
   unsigned int z_start;
   unsigned int z_end;
   double one_cell_size = 0.0; //Triggers floating point exception if not changed (intended).

   double direction;
   double coarse_fluxes[6][FF::ANOE()][CC::TCX()][CC::TCY()][CC::TCZ()];
   double fine_fluxes[6][FF::ANOE()][CC::TCX()][CC::TCY()][CC::TCZ()];

   // This initialization is crucial
   for( unsigned int b = 0; b < 6; b++ ) {
      for( const Equation eq : FF::ASOE() ) {
         for( unsigned int i = 0; i < CC::TCX(); ++i ) {
            for( unsigned int j = 0; j < CC::TCY(); ++j ) {
               for( unsigned int k = 0; k < CC::TCZ(); ++k ) {
                  coarse_fluxes[b][ETI(eq)][i][j][k] = 0.0;
                  fine_fluxes[b][ETI(eq)][i][j][k] = 0.0;
               }
            }
         }
      }
   }
   for( const auto& level : level_exchanging ) {
      leaf_ids_on_level = topology_.LeafIdsOnLevel(level);
      for( const auto& leaf_id : leaf_ids_on_level ) {
         for( const MaterialName material : topology_.GetFluidsOfNode(leaf_id) ) {
            for( const auto& location : CC::ANBS() ) {
               neighbor_id = topology_.GetTopologyNeighborId( leaf_id, location );
               neighbor_exists = topology_.NodeExists( neighbor_id );
               if( neighbor_exists ) {
                  neighbor_is_leaf = topology_.NodeIsLeaf(neighbor_id);
               } else {
                  neighbor_is_leaf = false;
               }
               neighbor_location = OppositeDirection(location);
               x_start = CC::FICX();
               x_end   = CC::LICX();
               y_start = CC::FICY();
               y_end   = CC::LICY();
               z_start = CC::FICZ();
               z_end   = CC::LICZ();
               if( neighbor_exists && !neighbor_is_leaf ) {
                  if( topology_.NodeIsOnRank( leaf_id, my_rank ) ) { // I have the Node, I must collect the data and update

                     direction = 0.0;
                     switch (location) {
                        case BoundaryLocation::East:
                           x_start = CC::LICX();
                           direction = -1.0;
                           break;
                        case BoundaryLocation::West:
                           x_end = CC::FICX();
                           direction = 1.0;
                           break;
                        case BoundaryLocation::North:
                           y_start = CC::LICY();
                           direction = -1.0;
                           break;
                        case BoundaryLocation::South:
                           y_end = CC::FICY();
                           direction = 1.0;
                           break;
                        case BoundaryLocation::Top:
                           z_start = CC::LICZ();
                           direction = -1.0;
                           break;
                        case BoundaryLocation::Bottom:
                           z_end = CC::FICZ();
                           direction = 1.0;
                           break;
#ifndef PERFORMANCE
                        default:
                           throw std::invalid_argument(" Why, oh why, did my simulation break?");
#endif
                     }
#ifndef PERFORMANCE
                     if( direction == 0.0 ) { throw std::logic_error("No no no"); }
#endif
                     Node& node = tree_.GetNodeWithId( leaf_id );
                     one_cell_size = 1.0 / node.GetCellSize();
                     Block& block = node.GetPhaseByMaterial(material);
                     double (&jump_buffer)[FF::ANOE()][CC::ICY()][CC::ICZ()] = block.GetBoundaryJumpConservatives( location );
                     unsigned int jump_index_one = 0;
                     unsigned int jump_index_two = 0;

                     // Update Setup One
                     for( const Equation eq : FF::ASOE() ) {
                        jump_index_one = 0;
                        jump_index_two = 0;
                        for( unsigned int i = x_start; i <= x_end; ++i ) {
                           for( unsigned int j = y_start; j <= y_end; ++j ) {
                              for( unsigned int k = z_start; k <= z_end; ++k ) {
                                 coarse_fluxes[LTI(location)][ETI(eq)][i][j][k] = jump_buffer[ETI(eq)][jump_index_one][jump_index_two] * one_cell_size * direction;
                                 jump_index_two++;
                                 // 3D: counter equal to IC, 2D: equal to 1, hence use std::min
                                 if( jump_index_two == std::min(CC::ICX(), std::min(CC::ICY(), CC::ICZ())) ) {
                                    jump_index_one++;
                                    jump_index_two = 0;
                                 }
                              }
                           }
                        }
                     }
                     if( topology_.NodeIsOnRank( neighbor_id, my_rank ) ) {
                        // Non-MPI
                        const double (&neighbor_jump_buffer)[FF::ANOE()][CC::ICY()][CC::ICZ()] = tree_.GetNodeWithId( neighbor_id ).GetPhaseByMaterial(
                                material).GetBoundaryJumpConservatives( neighbor_location );
                        for( unsigned int e = 0; e < FF::ANOE(); ++e ) {
                           for( unsigned int i = 0; i < CC::ICY(); ++i ) {
                              for( unsigned int j = 0; j < CC::ICZ(); ++j ) {
                                 jump_buffer[e][i][j] = neighbor_jump_buffer[e][i][j];
                              }
                           }
                        }
                     } else {
                        // MPI Recv
                        // Overwrites directly into the jump_buffer
                        MPI_Recv(jump_buffer, JumpBufferSendingSize(), MPI_DOUBLE, topology_.GetRankOfNode(neighbor_id), 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                     }

                     // Update Step 2
                     jump_index_one = 0;
                     jump_index_two = 0;
                     for( const Equation eq : FF::ASOE() ) {
                        jump_index_one = 0;
                        jump_index_two = 0;
                        for( unsigned int i = x_start; i <= x_end; ++i ) {
                           for( unsigned int j = y_start; j <= y_end; ++j ) {
                              for( unsigned int k = z_start; k <= z_end; ++k ) {
                                 fine_fluxes[LTI(location)][ETI(eq)][i][j][k] = jump_buffer[ETI(eq)][jump_index_one][jump_index_two] * one_cell_size * direction;
                                 jump_index_two++;
                                 // 3D: counter equal to IC, 2D: equal to 1, hence use std::min
                                 if( jump_index_two == std::min(CC::ICX(), std::min(CC::ICY(), CC::ICZ())) ) {
                                    jump_index_one++;
                                    jump_index_two = 0;
                                 }
                              }
                           }
                        }
                     }
                  } else {
                     if( topology_.NodeIsOnRank( neighbor_id, my_rank ) ) { // I do not have the Node but I do have the Neighbor and must send it
                        // Send MPI
                        double (& neighbor_jump_buffer)[FF::ANOE()][CC::ICY()][CC::ICZ()] = tree_.GetNodeWithId( neighbor_id ).GetPhaseByMaterial(
                                material).GetBoundaryJumpConservatives( neighbor_location );
                        MPI_Send( neighbor_jump_buffer, JumpBufferSendingSize(), MPI_DOUBLE, topology_.GetRankOfNode(leaf_id), 0, MPI_COMM_WORLD );
                     }
                  }
               } // Neighbor qualifies for exchange
            } // location

            // Now add up/subtract everything here
            if( topology_.NodeIsOnRank( leaf_id, my_rank ) ) {
               Node& node = tree_.GetNodeWithId( leaf_id );
               Block& block = node.GetPhaseByMaterial(material);
               for( const Equation eq : FF::ASOE() ) {
                  double (& cells)[CC::TCX()][CC::TCY()][CC::TCZ()] = block.GetRightHandSideBuffer(eq);
                  for( unsigned int i = CC::FICX(); i <= CC::LICX(); ++i ) {
                     for( unsigned int j = CC::FICY(); j <= CC::LICY(); ++j ) {
                        for( unsigned int k = CC::FICZ(); k <= CC::LICZ(); ++k ) {
                           cells[i][j][k] -= ConsistencyManagedSum(coarse_fluxes[0][ETI(eq)][i][j][k] + coarse_fluxes[1][ETI(eq)][i][j][k],
                                                                   coarse_fluxes[2][ETI(eq)][i][j][k] + coarse_fluxes[3][ETI(eq)][i][j][k],
                                                                   coarse_fluxes[4][ETI(eq)][i][j][k] + coarse_fluxes[5][ETI(eq)][i][j][k]);
                           cells[i][j][k] += ConsistencyManagedSum(fine_fluxes[0][ETI(eq)][i][j][k] + fine_fluxes[1][ETI(eq)][i][j][k],
                                                                   fine_fluxes[2][ETI(eq)][i][j][k] + fine_fluxes[3][ETI(eq)][i][j][k],
                                                                   fine_fluxes[4][ETI(eq)][i][j][k] + fine_fluxes[5][ETI(eq)][i][j][k]);
                           for( unsigned int b = 0; b < 6; b++ ) {
                              coarse_fluxes[b][ETI(eq)][i][j][k] = 0.0;
                              fine_fluxes[b][ETI(eq)][i][j][k] = 0.0;
                           }
                        }
                     }
                  }
               }
            } // Now add up/subtract everything here
         } // loop over materials
      } // leaves on id
   } // levels to exchange

   /*** Resetting Buffer ***/
   ResetJumpConservativeBuffers( finished_levels_descending );
}

/**
 * @brief Calculates the prime states from given conservatives. This is done for leaves on the given level list which do not have a
 * level-set block.
 * @tparam c Template parameter that specifies which conservative buffer (average, right-hand side or initial) is used to calculate the prime
 * @param updated_levels A vector containing the levels for whose leaves the prime states are calculated.
 * @param cut_levelset_nodes Indicates whether nodes having a level-set block are skipped or not.
 */
template<ConservativeBufferType C>
void ModularAlgorithmAssembler::ObtainPrimeStatesFromConservatives( std::vector<unsigned int> const updated_levels, bool const skip_levelset_nodes ) const {
   for(const unsigned int &level : updated_levels) {
      for(Node &non_levelset_node : tree_.NonLevelsetLeaves( level )) {
         DoObtainPrimeStatesFromConservativesForNonLevelsetNodes<C>(non_levelset_node);
      } // nodes without levelset
      if( !skip_levelset_nodes && level == setup_.GetMaximumLevel() ) {
         for( Node& node : tree_.NodesWithLevelset() ) {
            DoObtainPrimeStatesFromConservativesForLevelsetNodes<C>(node);
         }
      } // nodes with levelset
   } // levels
}

/**
 * @brief Calculates the prime states from given conservatives. This is done for the given node.
 * @tparam c Template parameter that specifies which conservative buffer (average, right-hand side or initial) is used to calculate the prime.
 * @param node The node for which the prime states are calculated.
 */
template<ConservativeBufferType C>
void ModularAlgorithmAssembler::DoObtainPrimeStatesFromConservativesForNonLevelsetNodes( Node& node ) const {

   for( auto& phase : node.GetPhases() ) {
      PrimeStates& prime_states = phase.second.GetPrimeStateBuffer();
      Conservatives const& conservatives = phase.second.GetConservativeBuffer<C>();
      prime_state_handler_.ConvertConservativesToPrimeStates( phase.first, conservatives, prime_states );
   } // phases
}

/**
 * @brief Calculates the prime states from given conservatives. This is done for the given node.
 * @tparam c Template parameter that specifies which conservative buffer (average, right-hand side or initial) is used to calculate the prime.
 * @param node The node for which the prime states are calculated.
 */
template<ConservativeBufferType C>
void ModularAlgorithmAssembler::DoObtainPrimeStatesFromConservativesForLevelsetNodes( Node& node ) const {
   const std::int8_t (&interface_tags)[CC::TCX()][CC::TCY()][CC::TCZ()] = node.GetInterfaceTags();
   for(auto& phase : node.GetPhases()) {
      PrimeStates& prime_states = phase.second.GetPrimeStateBuffer();
      const Conservatives& conservatives = phase.second.GetConservativeBuffer<C>();
      const MaterialName& material = phase.first;
      const std::int8_t material_sign = MaterialSignCapsule::SignOfMaterial(material);

      for(unsigned int i = 0; i < CC::TCX(); ++i) {
         for(unsigned int j = 0; j < CC::TCY(); ++j) {
            for(unsigned int k = 0; k < CC::TCZ(); ++k) {
               if(interface_tags[i][j][k] * material_sign > 0 || std::abs( interface_tags[i][j][k] ) <= ITTI(IT::ExtensionBand)) {
                  prime_state_handler_.ConvertConservativesToPrimeStates( material, conservatives, prime_states, i, j, k );
               } else {
                  for( PrimeState const p : FF::ASOP() ) {
                     prime_states[p][i][j][k] = 0.0;
                  }
               }
            } // k
         } // j
      } // i
   } //phases
}

/**
 * @brief Determines the maximal allowed size of the next time step (on the finest level).
 * @return Largest non-cfl-violating time step size on the finest level.
 */
double ModularAlgorithmAssembler::ComputeTimestepSize() const {

   std::array<double, DTI(CC::DIM())> velocity_plus_sound;
   double dt = 0.0;
   double sum_of_signalspeeds = 0.0;

   double nu = 0.0;
   double sigma = 0.0;
   double g = 0.0;

   //scaling for viscosity and surface tension limited timestep size - values taken from \cite Sussman2000
   double constexpr nu_timestep_size_constant = 3.0/14.0 ;
   double constexpr sigma_timestep_size_constant = M_PI*8.0;

   //scaling for thermal-diffusivity limited timestep size - value taken from \cite Pieper2016
   double constexpr thermal_diffusivity_dt_constant = 0.1;
   double thermal_diffusivity = 0.0;

   for(Node& node : tree_.Leaves()) {
      for(auto const& [material, block] : node.GetPhases()) {
         const auto material_sign = MaterialSignCapsule::SignOfMaterial(material);
         const PrimeStates& prime_states = block.GetPrimeStateBuffer();

         const std::int8_t (&interface_tags)[CC::TCX()][CC::TCY()][CC::TCZ()] = node.GetInterfaceTags();

         const double thermal_conductivity = material_manager_.GetThermalConductivity(material);
         const double specific_heat        = material_manager_.GetSpecificHeat(material);
         double const thermal_conductivity_over_specific_heat = specific_heat != 0.0 ? thermal_conductivity / specific_heat : 0.0;

         std::array<double, 3> const gravity = setup_.GetGravity();
         if constexpr( CC::GravityIsActive() ) {
            // using DimensionAwareConsistencyManagedSum ignores gravity entries of higher dimensions (e.g. gravity[2] for 2D) even if they are nonzero
            g = std::max( g, std::sqrt( DimensionAwareConsistencyManagedSum( gravity[0] * gravity[0], gravity[1] * gravity[1], gravity[2] * gravity[2] ) ) );
         }

         for( unsigned int i = CC::FICX(); i <= CC::LICX(); ++i ) {
            for( unsigned int j = CC::FICY(); j <= CC::LICY(); ++j ) {
               for( unsigned int k = CC::FICZ(); k <= CC::LICZ(); ++k ) {
                  // only consider cells with the correct sign of the tag
                  // cut cells (sign=0) are not considered!
                  if( interface_tags[i][j][k] * material_sign > 0 || std::abs( interface_tags[i][j][k] ) == ITTI(IT::NewCutCell) ) {
                     //only required when Euler equations are solved
                     if constexpr( CC::InviscidExchangeActive() ) {
                        double const c = material_manager_.GetSpeedOfSound( material, prime_states[PrimeState::Density][i][j][k], prime_states[PrimeState::Pressure][i][j][k] );

                        for( unsigned int d = 0; d < DTI(CC::DIM()); ++d ) {
                           velocity_plus_sound[d] = std::abs( prime_states[FF::AV()[d]][i][j][k] ) + c;
                        }

                        sum_of_signalspeeds = std::max( sum_of_signalspeeds, ConsistencyManagedSum(velocity_plus_sound) );
                     }

                     double const one_density = 1.0 / prime_states[PrimeState::Density][i][j][k];

                     //only required for viscous cases
                     if constexpr( CC::ViscosityIsActive() ) {
                        nu = std::max(nu, material_manager_.GetViscosity(material)[0] * one_density);
                     }

                     //only required if surface tension is active
                     if constexpr( CC::CapillaryForcesActive() ) {
                        sigma = std::max(sigma, material_manager_.GetSurfaceTensionCoefficient(MaterialSignCapsule::PositiveFluidMaterial(), MaterialSignCapsule::NegativeFluidMaterial()) * one_density);
                     }

                     //only required if heat conduction is considered
                     if constexpr( CC::HeatConductionActive() ){
                        thermal_diffusivity = std::max(thermal_diffusivity, thermal_conductivity_over_specific_heat * one_density );
                     }

                  }
               } // k
            } // j
         } // i

      } // block
   } // node

   /* We need the smallest possible cell size (and not the smallest currently present one) as the timestep on the other levels is deduced from it later in the algorithm.
    * An example is if Lmax is "missing" at the time of this computation the timesteps on the other levels would be too large by a factor of 2.
    */
   double const smallest_cell_size = setup_.SmallestPossibleCellSize();
   dt = sum_of_signalspeeds / smallest_cell_size;

   if constexpr( CC::ViscosityIsActive() ) {
      dt = std::max(dt, (nu / (nu_timestep_size_constant * smallest_cell_size * smallest_cell_size)));
   }

   if constexpr( CC::CapillaryForcesActive() ) {
      dt = std::max(dt, std::sqrt(sigma_timestep_size_constant * sigma) / (std::pow(smallest_cell_size,1.5)));
   }

   if constexpr( CC::HeatConductionActive() ){
      dt = std::max(dt, thermal_diffusivity / ( smallest_cell_size * smallest_cell_size * thermal_diffusivity_dt_constant ) );
   }

   //described in \cite Sussman 2000
   if constexpr( CC::GravityIsActive() ){
      dt = std::max( dt, 0.5* ( sum_of_signalspeeds + std::sqrt( sum_of_signalspeeds * sum_of_signalspeeds + 4.0 * g * smallest_cell_size ) ) / smallest_cell_size );
   }

   double local_dt_on_finest_level;

   //NH Check against double is okay in this case (and this case only!)
   if(dt == 0.0) { // The rank does not have any nodes, thus it could not compute a dt
      local_dt_on_finest_level = std::numeric_limits<double>::max();
   } else {
      local_dt_on_finest_level = setup_.GetCflNumber() / dt;
   }

   //limit the time-step size in the last macro time step to the exact end time
   if constexpr(CC::LET()){
      const std::vector<double> micro_time_steps = time_integrator_.MicroTimestepSizes();
      const double current_run_time = std::accumulate(micro_time_steps.cbegin(),micro_time_steps.cend(),time_integrator_.CurrentRunTime());

      if(current_run_time + local_dt_on_finest_level > setup_.GetEndTime()){
         local_dt_on_finest_level = setup_.GetEndTime() - current_run_time;

         //just for safety - make sure that the limiter does not lead to negative time-step sizes for the last output
         if(local_dt_on_finest_level < 0.0) {
            throw std::logic_error("Time-step limiter gives negative time-step size!");
         }
      }
   }

   double global_min_dt;
   //NH 2016-10-28 dt_in_finest_level needs to be the GLOBAL minimum of the computed values
   MPI_Allreduce(&local_dt_on_finest_level, &global_min_dt, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);

   logger_.LogMessage( "Timestep = " + ToScientificNotationString( setup_.DimensionalizeTime( global_min_dt ), 9 ) );

   return global_min_dt;
}

/**
 * @brief Sets all values in all jump buffers on all levels to 0.0.
 */
void ModularAlgorithmAssembler::ResetAllJumpBuffers() const {

   for(auto& level : tree_.FullNodeList()) {
      for(auto& id_node: level) {
         for(auto& phase : id_node.second.GetPhases()) {
            for(const auto& location : CC::ANBS()) {
               phase.second.ResetJumpConservatives( location );
               phase.second.ResetJumpFluxes( location );
            } //locations
         } //phases
      } //nodes
   } //levels
}

/**
 * @brief For all nodes on the provided levels resets, i.e. set values to 0.0, the buffers for the conservative values at jump faces.
 * @param levels The level on which buffers are reset.
 */
void ModularAlgorithmAssembler::ResetJumpConservativeBuffers( std::vector<unsigned int> const levels ) const {

   for(const auto& level : levels) {
      for(Node& node : tree_.NodesOnLevel(level) ) {
         for(auto& phase : node.GetPhases()) {
            for(auto& location : CC::ANBS()) {
               phase.second.ResetJumpConservatives( location );
            } //locations
         } //phases
      } //node
   } //levels
}

/**
 * @brief Checks if DoLoadBalancing is necessary and eventually executes it.
 * @param updated_levels_descending Gives the list of levels which have integrated values in the righthand side buffers. For these sending the RHS Buffer is enough.
 * @param force Enforces execution of load balancing independent of other indicators.
 */
void ModularAlgorithmAssembler::LoadBalancing( std::vector<unsigned int> const updated_levels_descending, bool const force ) {
   if( topology_.IsLoadBalancingNecessary() || force ) {
      //id - Current Rank - Future Rank
      const std::vector<std::tuple<const std::uint64_t, const int, const int>> ids_rank_map = topology_.GetLoadBalancedTopology( MpiUtilities::NumberOfRanks() );
      // ^ Changes the rank assignment in the Topology.
      communicator_.InvalidateCache();

      const MPI_Datatype conservatives_datatype = communicator_.ConservativesDatatype();
      const MPI_Datatype boundary_jump_datatype = communicator_.JumpSurfaceDatatype();

      std::vector<MPI_Request> requests;

      int const my_rank_id = MpiUtilities::MyRankId();

      for( auto const& [id, current_rank, future_rank] : ids_rank_map) { //We traverse current topology
         /*If the node has not been updated, i.e. integrated values in RHS buffer, we need to handle the AVG buffer as well*/
         const bool send_averages = std::find( updated_levels_descending.begin(), updated_levels_descending.end(), LevelOfNode( id ) ) == updated_levels_descending.end();

         if( current_rank == my_rank_id ) {
            Node const& node = tree_.GetNodeWithId( id );
            for( MaterialName const material : topology_.GetFluidsOfNode( id ) ) {
               Block const& block = node.GetPhaseByMaterial( material );

               communicator_.Send( &block.GetRightHandSideBuffer(), FF::ANOE(), conservatives_datatype, future_rank, requests );
               /*If the node has not been updated we need to send the AVG buffer as well*/
               if(send_averages) {
                  communicator_.Send( &block.GetAverageBuffer(), FF::ANOE(), conservatives_datatype, future_rank, requests );
                  communicator_.Send( &block.GetInitialBuffer(), FF::ANOE(), conservatives_datatype, future_rank, requests );
               }
               communicator_.Send( &block.GetBoundaryJumpFluxes(), CC::SIDES(), boundary_jump_datatype, future_rank, requests );
               communicator_.Send( &block.GetBoundaryJumpConservatives(), CC::SIDES(), boundary_jump_datatype, future_rank, requests );
            }
            if( topology_.IsNodeMultiPhase( id ) ) {
               communicator_.Send( node.GetInterfaceTags(), FullBlockSendingSize(), MPI_INT8_T, future_rank, requests );
               if( node.HasLevelset() ) { // Lmax node

                  for( MaterialName const material : topology_.GetFluidsOfNode( id ) ) {
                     communicator_.Send( &node.GetPhaseByMaterial( material ).GetPrimeStateBuffer(), FF::ANOP(), conservatives_datatype, future_rank, requests );
                  }
                  communicator_.Send( node.GetLevelsetBlock().GetPhiReinitialized(), FullBlockSendingSize(), MPI_DOUBLE, future_rank, requests );
                  /**
                   * It is not necessary to send the phi_initial buffer, since the levelset is integrated each micro timestep.
                   */
                  communicator_.Send( node.GetLevelsetBlock().GetVolumeFraction(), FullBlockSendingSize(), MPI_DOUBLE, future_rank, requests );
                  communicator_.Send( node.GetLevelsetBlock().GetInterfaceQuantityBuffer( InterfaceQuantity::Velocity ), FullBlockSendingSize(), MPI_DOUBLE, future_rank, requests );
               }
            }
            if constexpr( DP::Profile() ) {
               CommunicationStatistics::balance_send_++;
            }
         } else if( future_rank == my_rank_id ) { // The node is currently NOT ours, but will be in the future
            // Create Node first, then post asynchronous Recv.
            Node& new_node = tree_.CreateNode( id, topology_.GetFluidsOfNode( id ) );

            for( MaterialName const material : topology_.GetFluidsOfNode( id ) ) {
               Block& block = new_node.GetPhaseByMaterial( material );
               communicator_.Recv( &block.GetRightHandSideBuffer(), FF::ANOE(), conservatives_datatype, current_rank, requests );
               /*If the node has not been updated, we need to receive (and apply) the AVG buffer as well*/
               if(send_averages) {
                  communicator_.Recv( &block.GetAverageBuffer(), FF::ANOE(), conservatives_datatype, current_rank, requests );
                  communicator_.Recv( &block.GetInitialBuffer(), FF::ANOE(), conservatives_datatype, current_rank, requests );
               }
               communicator_.Recv( &block.GetBoundaryJumpFluxes(), CC::SIDES(), boundary_jump_datatype, current_rank, requests );
               communicator_.Recv( &block.GetBoundaryJumpConservatives(), CC::SIDES(), boundary_jump_datatype, current_rank, requests );
            }
            if( topology_.IsNodeMultiPhase( id ) ) {
               communicator_.Recv( &new_node.GetInterfaceTags(), FullBlockSendingSize(), MPI_INT8_T, current_rank, requests );
               if( LevelOfNode( id ) == setup_.GetMaximumLevel() ) {
                  for( MaterialName const material : topology_.GetFluidsOfNode( id ) ) {
                     communicator_.Recv(&new_node.GetPhaseByMaterial( material ).GetPrimeStateBuffer(), FF::ANOP(), conservatives_datatype, current_rank, requests );
                  }
                  // We have not yet created a LS field in our recieving Node. At this point it is clear it need one, so we create it with dummys and receive the correct values.
                  new_node.SetLevelsetBlock( std::make_unique<LevelsetBlock>( 0.0 ) );
                  communicator_.Recv( new_node.GetLevelsetBlock().GetPhiReinitialized(), FullBlockSendingSize(), MPI_DOUBLE, current_rank, requests );
                  communicator_.Recv( new_node.GetLevelsetBlock().GetVolumeFraction(), FullBlockSendingSize(), MPI_DOUBLE, current_rank, requests );
                  communicator_.Recv( new_node.GetLevelsetBlock().GetInterfaceQuantityBuffer( InterfaceQuantity::Velocity ), FullBlockSendingSize(), MPI_DOUBLE, current_rank, requests );
               }
            } else {
               std::int8_t uniform_tag = MaterialSignCapsule::SignOfMaterial( topology_.GetFluidsOfNode( id ).back() ) * ITTI( IT::BulkPhase );
               std::int8_t(&new_tags)[CC::TCX()][CC::TCY()][CC::TCZ()] = new_node.GetInterfaceTags();
               for(unsigned int i = 0; i < CC::TCX(); ++i) {
                  for(unsigned int j = 0; j < CC::TCY(); ++j) {
                     for(unsigned int k = 0; k < CC::TCZ(); ++k) {
                        new_tags[i][j][k] = uniform_tag;
                     }
                  }
               }
            }
            if constexpr( DP::Profile() ) {
               CommunicationStatistics::balance_recv_++;
            }
         }
      }
      MPI_Waitall( requests.size(), requests.data(), MPI_STATUSES_IGNORE );
      requests.clear();

      //remove nodes only after Data was received by partner
      for( unsigned int i = 0; i < ids_rank_map.size(); i++ ) { //We traverse current topology
         const std::uint64_t id = std::get<0>( ids_rank_map[i] );
         const int current_rank = std::get<1>( ids_rank_map[i] );
         if( current_rank == my_rank_id ) {
            tree_.RemoveNodeWithId( id );
         }
      }

      logger_.LogMessage( "Load Balancing (" + std::to_string( ids_rank_map.size() ) + ")" );
   }
}

/**
 * @brief Sets the values in the internal cells on the given level to match the user-input initial condition
 * @param level The level on which the initial condition values are applied to.
 */
void ModularAlgorithmAssembler::ImposeInitialCondition( unsigned int const level ) {

   double initial_prime_states[FF::ANOP()][CC::ICX()][CC::ICY()][CC::ICZ()];

   for( auto& [id, node]: tree_.GetLevelContent( level ) ) {
      for( const MaterialName material : topology_.GetFluidsOfNode( id ) ) {

         const std::int8_t (&interface_tags)[CC::TCX()][CC::TCY()][CC::TCZ()] = node.GetInterfaceTags();

         Block& block = node.GetPhaseByMaterial( material );

         Conservatives& conservatives = block.GetRightHandSideBuffer();

         setup_.GetInitialPrimeStates( id, material, initial_prime_states );

         const auto material_sign = MaterialSignCapsule::SignOfMaterial( material );

         for( unsigned int i = CC::FICX(); i <= CC::LICX(); ++i ) {
            for( unsigned int j = CC::FICY(); j <= CC::LICY(); ++j ) {
               for( unsigned int k = CC::FICZ(); k <= CC::LICZ(); ++k ) {
                  // only assign in interface/narrowband cells and in bulk cells of same sign
                  // the following evaluates to true only in bulk cells of DIFFERENT sign
                  if( material_sign * interface_tags[i][j][k] > 0 || std::abs( interface_tags[i][j][k] ) < ITTI(IT::BulkPhase ) ) {
                     std::array<double, FF::ANOP()> prime_states_cell;
                     for( unsigned int p = 0; p < FF::ANOP(); ++p ) {
                        prime_states_cell[p] = initial_prime_states[p][i - CC::FICX()][j - CC::FICY()][k - CC::FICZ()];
                     }
                     auto&& conservatives_cell = conservatives.GetCellView( i, j, k );
                     prime_state_handler_.ConvertPrimeStatesToConservatives( material, prime_states_cell, conservatives_cell );
                  } else {
                     for( Equation const e : FF::ASOE() ) {
                        conservatives[e][i][j][k] = 0.0;
                     }
                  }
               } //k
            } //j
         } //i
      } //phases
   } //nodes
}

void ModularAlgorithmAssembler::UpdateTopology() {
   if(topology_.UpdateTopology()) {
      //Caches for Halo Updates are only invalid if nodes have been refined, coarsened or moved. Changes in the number of fluids are no Problem.
      communicator_.InvalidateCache();
   }
}

/**
 * @brief Changes the refinement depth in the domain. I.e. coarsens and/or refines blocks on certain levels.
 *        All levels in the input list may be refined but only those which also find their parent level in
 *        the input list may be coarsened, in order to maintain a consistent buffer state.
 * @param levels_to_update_ascending List of the levels to be considered for coarsening/refinement.
 */
void ModularAlgorithmAssembler::Remesh( std::vector<unsigned int> const levels_to_update_ascending ) {

   /*** We create run the Wavelet analysis on finished levels ***/
   std::vector<unsigned int> parent_levels( levels_to_update_ascending );
   int const number_of_ranks = MpiUtilities::NumberOfRanks();

   // The maximum level is not a parent hence it is removed form the parent list.
   parent_levels.erase(std::remove_if(parent_levels.begin(), parent_levels.end(),
                                      [&](const unsigned int level) { return level >= topology_.GetCurrentMaximumLevel(); }), parent_levels.end());
   std::vector<std::uint64_t> nodes_to_be_coarsened;
   std::vector<std::uint64_t> nodes_needing_refinement;
   DetermineRemeshingNodes( parent_levels, nodes_to_be_coarsened, nodes_needing_refinement );

   /* First we deal with the refinement. We keep the nodes to be coarsened until after the halo update, which is need in the refinement process,
    * to reduce the possible number of jumphalos.
    */

   // We need to cut all non-leaves and all Lmax-leaves from the refinement list.
   nodes_needing_refinement.erase(std::remove_if(nodes_needing_refinement.begin(), nodes_needing_refinement.end(),
                                                 [&](const std::uint64_t id) { return ( !topology_.NodeIsLeaf(id) || LevelOfNode(id) == setup_.GetMaximumLevel() ); }),
                                  nodes_needing_refinement.end());

   // Global Distibution of the refine list
   std::vector<std::uint64_t> global_refine_list;
   MpiUtilities::LocalToGlobalData(nodes_needing_refinement, MPI_LONG_LONG_INT, number_of_ranks, global_refine_list);

   // Duplicates can not exist - no check needed
   for( const std::uint64_t leaf_id : global_refine_list ) {
      if( topology_.NodeIsOnRank( leaf_id, communicator_.MyRankId() ) ) {
         RefineNode(leaf_id);
      }
   }
   UpdateTopology();
   std::vector<unsigned int> halo_levels( levels_to_update_ascending );
   halo_levels.erase(halo_levels.begin()); // The lowest level did not change during refinement

   /*** We need to fill the Halo Cells ***/
   halo_manager_.FluidHaloUpdate( halo_levels, FluidFieldType::Conservatives );

   /* Second we deal with coarsening*/

   // Global Distribution of the coarse list
   std::vector<std::uint64_t> global_remove_list;
   MpiUtilities::LocalToGlobalData( nodes_to_be_coarsened, MPI_LONG_LONG_INT, number_of_ranks, global_remove_list );

   //Gives ALL local nodes which need to be deleted
   std::vector<std::uint64_t> local_cut;
   std::copy_if(global_remove_list.begin(), global_remove_list.end(), std::back_inserter(local_cut),
                [&](const std::uint64_t id) { return topology_.NodeIsOnRank( id, communicator_.MyRankId() ); });

   // Coarsening operation is called on parents, hence we convert the coarse list to a list of the respective parents.
   // NH TODO-19 this cumbersome step should go away
   std::vector<std::uint64_t> parents_of_coarsened = global_remove_list;
   std::for_each(parents_of_coarsened.begin(), parents_of_coarsened.end(), [](std::uint64_t& to_parent) { to_parent = ParentIdOfNode(to_parent); });

   // Housekeeping - Due to ancestors in list duplicates may arise - these need to be cut
   // Sort - Erase - Unique - Idiom
   std::sort(parents_of_coarsened.begin(), parents_of_coarsened.end());
   parents_of_coarsened.erase(std::unique(parents_of_coarsened.begin(), parents_of_coarsened.end()), parents_of_coarsened.end());

   // Level zero parents are not allowed, as this means a coarsening of level 1.
   parents_of_coarsened.erase(std::remove_if(parents_of_coarsened.begin(), parents_of_coarsened.end(), [&](const std::uint64_t parent_id) { return LevelOfNode(parent_id) == 0; }),
                              parents_of_coarsened.end());

   // Updating the topology (light data)
   for( const std::uint64_t parent_id : parents_of_coarsened ) {
      topology_.CoarseNodeWithId( parent_id );
   }
   if( !parents_of_coarsened.empty() ) {
      communicator_.InvalidateCache();
   }

   // Updating the tree (hard data)

   //Sort - Erase - Unique - Idiom
   std::sort(local_cut.begin(), local_cut.end());
   local_cut.erase(std::unique(local_cut.begin(), local_cut.end()), local_cut.end());

   // Level One may not be coarsened, ever.
   local_cut.erase(std::remove_if(local_cut.begin(), local_cut.end(), [&](const std::uint64_t id) { return (LevelOfNode(id) == 1); }), local_cut.end());
   for( const auto& id_to_be_removed : local_cut ) {
      tree_.RemoveNodeWithId( id_to_be_removed );
   }
}

/**
 * @brief Gives the ids of nodes which may be coarsened or need refinement according to the wavelet-analysis of \cite Harten 1993
 * @param parent_levels The levels of the parents, i.e. Children of these parents might be coarsened.
 * @param remove_list A list of all ids of nodes which may be coarsened (indirect return parameter).
 * @param refine_list A list of all ids of nodes which must be refined (indirect return parameter).
 */
void ModularAlgorithmAssembler::DetermineRemeshingNodes( std::vector<unsigned int> const parent_levels, std::vector<std::uint64_t>& remove_list,
                                                         std::vector<std::uint64_t>& refine_list ) const {

   std::vector<RemeshIdentifier> remesh_list;
   int const my_rank = communicator_.MyRankId();
   /**
    *  We need to check whether or not a node may be coarsened stay or even be refined. This is done by an analysis between the node and its parent.
    *  We only coarse or refine leaves and siblings may only be coarsened together.
    *  In two-phase simulations further checks are needed as multi nodes may only be leaves if they reside on Lmax.
    */
   MPI_Datatype const conservatives_struct_ = communicator_.ConservativesDatatype();
   for( const auto& level_of_parent : parent_levels ) {
      for( const auto& parent_id : topology_.GlobalIdsOnLevel(level_of_parent) ) {
         bool const parent_on_my_rank = topology_.NodeIsOnRank( parent_id, my_rank );
         std::vector<std::uint64_t> const children = IdsOfChildren(parent_id);
         for( const auto& child_id : children ) {
            if( topology_.NodeExists( child_id ) ) {
               if( !topology_.IsNodeMultiPhase(child_id) ) { // only single nodes may be coarsened or refined
                  if( topology_.NodeIsLeaf(child_id) ) { // We only check leaves (for now, TODO-19 NH).
                     bool const child_on_my_rank = topology_.NodeIsOnRank( child_id, my_rank );
                     if( parent_on_my_rank ) {
                        if( child_on_my_rank ) { // We hold parent and Child -> NO MPI

                           remesh_list.emplace_back( multiresolution_.ChildNeedsRemeshing<CC::NFWA()>( tree_.GetNodeWithId( parent_id ).GetPhaseByMaterial( topology_.SingleFluidOfNode( child_id ) ),
                                                                                                       tree_.GetNodeWithId( child_id ).GetSinglePhase(),
                                                                                                       child_id ) );

                        } else { // We hold parent but not the Child (which exists) -> MPI_Receive

                           // get material of child via topology (since it is single-phase taking front() of the material vector is valid)
                           Block received_child_block = Block();
                           int const sender_rank = topology_.GetRankOfNode(child_id);
                           MPI_Recv(&received_child_block.GetRightHandSideBuffer(), FF::ANOE(), conservatives_struct_, sender_rank, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

                           remesh_list.emplace_back( multiresolution_.ChildNeedsRemeshing<CC::NFWA()>( tree_.GetNodeWithId( parent_id ).GetPhaseByMaterial( topology_.SingleFluidOfNode( child_id ) ),
                                                                                                       received_child_block, child_id ) );
                        }
                     } else {
                        if( child_on_my_rank ) { // We do NOT hold the parent, but do hold the Child -> MPI Send
                           int receiver_rank = topology_.GetRankOfNode(parent_id);
                           const Block& send_child = tree_.GetNodeWithId( child_id ).GetSinglePhase();
                           MPI_Send( &send_child.GetRightHandSideBuffer(), FF::ANOE(), conservatives_struct_, receiver_rank, 0, MPI_COMM_WORLD );
                        }
                     }
                  } else { // IF: IsLeaf
                     //Add a meshing-lock, i. e. neutral
                     if( parent_on_my_rank ) { // Only the parent needs to fill this list.
                        remesh_list.emplace_back( RemeshIdentifier::Neutral );
                     }
                  } // ELSE: IsLeaf

               } else {// Child is multi-phase
                  // In case the child is multi-phase we don't do anything
                  if( parent_on_my_rank ) {
                     remesh_list.emplace_back( RemeshIdentifier::Neutral );
                  }
               } // else : child is multi-phase
            } // child node exists
         }//children

         //Now we have checked all siblings
#ifndef PERFORMANCE
         if( !remesh_list.empty() && remesh_list.size() != children.size() ) {
            throw std::logic_error("This must not happen");
         }
#endif
         for( unsigned int i = 0; i < remesh_list.size(); ++i ) {
            if( remesh_list[i] == RemeshIdentifier::Refine ) {
               refine_list.emplace_back( children[i] );
            }
         }
         // siblings may only be coarsened together. List can be empty if children do not exist.
         if( !remesh_list.empty() && !topology_.IsNodeMultiPhase(parent_id) &&
             std::all_of(remesh_list.begin(), remesh_list.end(), [](const RemeshIdentifier condition) { return condition == RemeshIdentifier::Coarse; }) ) {
            remove_list.insert(remove_list.end(), children.begin(), children.end());
         }
         remesh_list.clear();
      } // parents
   } // level_of_parent
}

/**
 * @brief Triggers the MPI consistent refinement of the node with the given id.
 * @param id Node identifier of the node to be refined.
 */
void ModularAlgorithmAssembler::RefineNode( std::uint64_t const id ) {
   topology_.RefineNodeWithId( id );
   std::vector<std::uint64_t> const ids_of_children = tree_.RefineNode(id);
   const Node& parent = tree_.GetNodeWithId( id );
   for( const auto& child_id : ids_of_children ) {
      Node& child = tree_.GetNodeWithId( child_id );
      // only single fluid nodes are supposed to be refined, hence using the single phase block is valid
      // only predicts internal cells
      for(const Equation eq : FF::ASOE()) {
         Multiresolution::Prediction(parent.GetSinglePhase().GetRightHandSideBuffer(eq), child.GetSinglePhase().GetRightHandSideBuffer(eq), child_id,
            CC::FICX(), CC::ICX(),CC::FICY(), CC::ICY(), CC::FICZ(), CC::ICZ());
      }
      // add material
      topology_.AddFluidToNode(child_id, parent.GetSinglePhaseMaterial());
   }
}

/**
 * @brief Gives a list of levels which need to be updated in the next stage depending on the chosen time-integration scheme. List is in descending order.
 * @param timestep The current micro time step.
 * @return List of levels to be updated in descending order.
 */
std::vector<unsigned int> ModularAlgorithmAssembler::GetLevels( unsigned int const timestep ) const {

   std::bitset<16> height_counter((timestep + 1) ^ (timestep));
   std::size_t height = height_counter.count();
   std::vector<unsigned int> levels_to_update(height);
   std::iota(levels_to_update.begin(), levels_to_update.end(), (setup_.GetMaximumLevel() - height + 1));
   std::reverse(levels_to_update.begin(), levels_to_update.end());
   return levels_to_update;
}

/**
 * @brief Logs the number of nodes currently active in the domain. Depending on compile-time settings the level of details varies.
 */
void ModularAlgorithmAssembler::LogNodeNumbers() const {
   auto&& [number_of_nodes, number_of_leaves] = topology_.NodeAndLeafCount();
   logger_.LogMessage( "Global number of nodes : " + std::to_string( number_of_nodes ) );
   logger_.LogMessage( "Global number of leaves: " + std::to_string( number_of_leaves ) );
   if constexpr( GeneralTwoPhaseSettings::LogMultiPhaseNodeCount ) {
      logger_.LogMessage( "Global number of multi-phase nodes: " + std::to_string( topology_.MultiPhaseNodeCount() ) );
   }
   if constexpr( GeneralTwoPhaseSettings::LogLevelsetLeafCount ) {
      unsigned int global_levelset_leaves = tree_.NodesWithLevelset().size();
      MPI_Allreduce( MPI_IN_PLACE, &global_levelset_leaves, 1, MPI_UNSIGNED, MPI_SUM, MPI_COMM_WORLD );
      logger_.LogMessage( "Global number of levelset leaves  : " + std::to_string( global_levelset_leaves ) );
   }
}

void ModularAlgorithmAssembler::LogElapsedTimeSinceInProfileRuns( double const start_time, std::string const function_name ) {
   if( DP::Profile() ){
      logger_.LogMessage( function_name + " - elapsed time: " + ToScientificNotationString( MPI_Wtime() - start_time ) );
   }
}
