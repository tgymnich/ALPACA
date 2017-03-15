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
#include "input_output_manager.h"

#include <cstdio> // needed for file deletion
#include <fstream>
#include <mpi.h>
#include <unistd.h>

#include "user_specifications/debug_and_profile_setup.h"
#include "xdmf_output_writer.h"
#include "communication/mpi_utilities.h"

/**
 * @brief Default constructor using the inputs from the operating system.
 * @param setup Instance to get access to user defined properties, relevant for the simulation.
 * @param topology TopologyManager instance for access of global node information.
 * @param flower The tree instance which gives access to the local fluid data.
 */
InputOutputManager::InputOutputManager(SimulationSetup const& setup, TopologyManager& topology, Tree& flower) :
   simulation_name_(RemoveFilePath( RemoveFileExtension( setup.InputFileName() ) )),
   output_folder_name_( SetupOutputFolder() ),
   setup_(setup),
   restart_manager_(flower, topology, setup, output_folder_name_),
   logger_(LogWriter::Instance()),
   output_writer_( CreateOutputWriter(topology, flower, setup, output_folder_name_, simulation_name_)),
   output_timestamps_(setup.OutputTimes()),
   snapshot_timestamps_(setup.SnapshotTimestamps()),
   symlink_latest_snapshot_name_(output_folder_name_ + RestartManager::RestartSubfolderName() + "/latest_snapshot"),
   wall_time_of_last_snapshot_(std::chrono::system_clock::now())
{
  if( MpiUtilities::MyRankId() == 0 ) {
    // copy inputfile into directory
    std::ifstream input_stream(setup.InputFileName(), std::ios::binary);
    std::ofstream output_stream(output_folder_name_ + "/" + RemoveFilePath(setup.InputFileName()), std::ios::binary);
    output_stream << input_stream.rdbuf();
    output_stream.flush();
    output_stream.close();
  }

  // logging
  logger_.SetLogfileName(output_folder_name_ + "/" + simulation_name_ + ".log");
  logger_.LogMessage("Simulation Name : " + simulation_name_);
  logger_.LogMessage("Output Folder   : " + output_folder_name_);
}

/**
 * @brief This routine creates the output folders for the simulation.
 * @return The name of the output folder.
 * @note Only the master-rank "zero" creates the folder.
 */
std::string InputOutputManager::SetupOutputFolder() const {

  std::string output_folder_name = AddUnusedNumberToPath( simulation_name_ );
  MPI_Barrier(MPI_COMM_WORLD); // This Barrier is needed, otherwise we get inconsistent folder names across the ranks.
  if( MpiUtilities::MyRankId() == 0 ) {
    // create folder and subfolders
    CreateFolder(output_folder_name);
    CreateFolder(output_folder_name + OutputWriter::OutputSubfolderName());
    CreateFolder(output_folder_name + RestartManager::RestartSubfolderName());
    if constexpr( DP::DebugOutput() ) {
      CreateFolder(output_folder_name + OutputWriter::DebugSubfolderName());
    }

#ifndef PERFORMANCE
    // create an .gitignore file in the output folder
    std::ofstream output_stream(output_folder_name + "/.gitignore", std::ios::out);
    output_stream << "*";
    output_stream.flush();
    output_stream.close();
#endif
  }

  return output_folder_name;
}

/**
 * @brief Create the output writer instance.
 * @param topology TopologyManager instance for access of global node information.
 * @param flower The tree instance which gives access to the local fluid data.
 * @param setup Instance to get access to user defined properties, relevant for the simulation.
 * @param output_folder_name Path to the output folder where output files are written into.
 * @param simulation_name Name of the simulation that is computed.
 * @return Pointer to output writer instance.
 */
std::unique_ptr<OutputWriter const> InputOutputManager::CreateOutputWriter(TopologyManager const& topology, Tree const& flower,
  SimulationSetup const& setup, std::string const output_folder_name, std::string const simulation_name ) const {
  // setup output writer
  switch(setup.GetOutputFormat()) {
    case OutputType::Xdmf:
      return std::make_unique<XdmfOutputWriter const>(flower, topology, setup, output_folder_name, simulation_name);
    default:
      throw std::invalid_argument("Selected output Type is not possible");
  }
}

/**
 * @brief Writes all micro timestep sizes performed during the last macro timestep.
 * @param timesteps_on_finest_level Sizes of the micro timesteps.
 */
void InputOutputManager::WriteTimestepFile(std::vector<double> const& timesteps_on_finest_level) const {
  if( MpiUtilities::MyRankId() == 0 ) {
    std::string filename = output_folder_name_ + "/time.txt";
    std::ofstream output_stream(filename, std::ios::app);
    for(const auto& timestep : timesteps_on_finest_level) {
      output_stream << timestep << "\n";
    }

    output_stream.flush();
    output_stream.close();
  }
}

/**
 * @brief Writes the full output (standard simulation output and, if applicable, debug output) if necessary at the current timestep.
 *        If the force_output flag is set, output is written in any case.
 * @param timestep The current timestep.
 * @param force_output A flag indicating whether output should be forced.
 * @return Return whether output was written or not.
 */
bool InputOutputManager::WriteFullOutput(double const timestep, bool const force_output) {
  if(force_output || output_timestamps_.front() <= timestep) {
    // erase the timestamps that are handled now
    output_timestamps_.erase(output_timestamps_.begin(), std::find_if_not(output_timestamps_.begin(),
       output_timestamps_.end(), [&timestep](double const timestamp){return timestamp <= timestep;}));

    // call the output functions
    double const dimensionalized_time = setup_.DimensionalizeTime(timestep);
    output_writer_->WriteOutputFile(dimensionalized_time);
    if constexpr( DP::DebugOutput() ) {
      output_writer_->WriteDebugFile(dimensionalized_time);
    }

    logger_.Flush();

    return true;
  }
  return false;
}

/**
 * @brief Writes debug output with the given debug key.
 * @param debug_key The debug key indicating at which algorithm substep the output is triggered.
 */
void InputOutputManager::WriteDebugFile(unsigned int const debug_key) const {
  output_writer_->WriteDebugFile(debug_key);
}

/**
 * @brief Writes a restart snapshot if necessary at the current timestep. Old snapshots are deleted depending on user configuration.
 * @param timestep The current timestep.
 * @param force_output A flag indicating whether output should be forced.
 */
void InputOutputManager::WriteRestartFile(double const timestep, bool const force_output) {
  bool snapshot_interval_triggered = false;
  // only consider interval-based snapshots if the interval is greater zero
  if(setup_.SnapshotInterval() > 0) {
    if( MpiUtilities::MyRankId() == 0 ) {
      std::chrono::time_point<std::chrono::system_clock> const current_wall_time = std::chrono::system_clock::now();
      int const wall_seconds_since_snapshot = std::chrono::duration_cast<std::chrono::seconds>(current_wall_time - wall_time_of_last_snapshot_).count();
      if(setup_.SnapshotInterval() <= wall_seconds_since_snapshot) {
        wall_time_of_last_snapshot_ = current_wall_time;
        snapshot_interval_triggered = true;
      }
    }
    // distribute restart decision among all ranks
    MPI_Bcast(&snapshot_interval_triggered, 1, MPI_CXX_BOOL, 0, MPI_COMM_WORLD);
  }

  bool const snapshot_timestamp_triggered = snapshot_timestamps_.front() <= timestep;
  if(snapshot_interval_triggered || snapshot_timestamp_triggered || force_output) {
    // erase the timestamps that are handled now
    snapshot_timestamps_.erase(snapshot_timestamps_.begin(), std::find_if_not(snapshot_timestamps_.begin(),
       snapshot_timestamps_.end(), [&timestep](double const timestamp){return timestamp <= timestep;}));

    // write the actual restart file and store its name
    std::string snapshot_filename = restart_manager_.WriteSnapshotFile(timestep);

    // handle filesystem access only on rank zero
    if( MpiUtilities::MyRankId() == 0 ) {
      // update symbolic link to latest snapshot file
      std::remove(symlink_latest_snapshot_name_.c_str());
      [[maybe_unused]] int const result_io = symlink(RemoveFilePath(snapshot_filename).c_str(), symlink_latest_snapshot_name_.c_str());

      // only consider non-timestamp snapshots for deletion
      if(!snapshot_timestamp_triggered) {
        if(snapshot_files_written_.size() == setup_.SnapshotsToKeep()) {
          // remove the oldest file
          std::remove(snapshot_files_written_.front().c_str());
          snapshot_files_written_.erase(snapshot_files_written_.begin());
        }
        // add the newest file to the list
        snapshot_files_written_.push_back(snapshot_filename);
      }
    }
  }
}

/**
 * @brief Restores the simulation from the user configured snapshot and returns the simulation time of the snapshot.
 * @return The time of the simulation snapshot.
 */
double InputOutputManager::RestoreSimulationFromSnapshot() {
  if(!CheckIfRestoreFileExists()) {
    throw std::runtime_error("Restart file '" + setup_.RestoreFileName() + "' does not exist!");
  }
  // call the manager that actually restores the simulation
  double const restart_time = restart_manager_.RestoreSimulation();

  // adapt timestamp lists to new start time
  output_timestamps_.erase(output_timestamps_.begin(), std::find_if_not(output_timestamps_.begin(),
     output_timestamps_.end(), [&restart_time](double const timestamp){return timestamp <= restart_time;}));
  snapshot_timestamps_.erase(snapshot_timestamps_.begin(), std::find_if_not(snapshot_timestamps_.begin(),
     snapshot_timestamps_.end(), [&restart_time](double const timestamp){return timestamp <= restart_time;}));

  return restart_time;
}
