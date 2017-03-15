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
#ifndef INPUT_OUTPUT_MANAGER_H
#define INPUT_OUTPUT_MANAGER_H

#include <chrono>
#include <memory>

#include "helper_functions.h"
#include "log_writer.h"
#include "output_writer.h"
#include "restart_manager.h"
#include "topology/topology_manager.h"
#include "topology/tree.h"

/**
 * @brief The InputOutputManager class handles creation of and access to a unique output folder
 *        and delegates all output calls. It decides whether simulation output or restart snapshots
 *        have to be written based on user configuration and calls the respective routines.
 */
class InputOutputManager {

   // path parameters
   std::string const simulation_name_;
   std::string const output_folder_name_;

   SimulationSetup const& setup_;

   RestartManager const restart_manager_;

   //Logger must not be const (otherwise no logbook cannot be appended
   LogWriter& logger_;

   // simulation output
   std::unique_ptr<OutputWriter const> const output_writer_;
   std::vector<double> output_timestamps_;

   // restart output
   std::vector<double> snapshot_timestamps_;
   std::vector<std::string> snapshot_files_written_;
   std::string const symlink_latest_snapshot_name_;
   // the wall clock time at the last restart output (only used on rank zero!)
   std::chrono::time_point<std::chrono::system_clock> wall_time_of_last_snapshot_;

   // Local factory functions
   std::string SetupOutputFolder() const;
   std::unique_ptr<OutputWriter const> CreateOutputWriter(TopologyManager const& topology, Tree const& flower,
      SimulationSetup const& setup, std::string const output_folder_name, std::string const simulation_name ) const;

public:
   InputOutputManager() = delete;
   explicit InputOutputManager( SimulationSetup const& setup, TopologyManager& topology, Tree& flower);
   ~InputOutputManager() = default;
   InputOutputManager( InputOutputManager const& ) = delete;
   InputOutputManager& operator=( InputOutputManager const& ) = delete;
   InputOutputManager( InputOutputManager&& ) = delete;
   InputOutputManager& operator=( InputOutputManager&& ) = delete;

   /**
   * @brief Checks whether the file "ABORTFILE" exists in the output folder indicating that the simulation should be aborted.
   * @return True if "ABORTFILE" exists, false otherwise.
   */
   inline bool CheckIfAbortfileExists() const {return CheckIfPathExists(OutputFolderName() + "/ABORTFILE");}

   /**
    * @brief Checks whether the restore file specified in the input file exists.
    * @return True is the file exists, false otherwise.
    */
   inline bool CheckIfRestoreFileExists() const {return CheckIfPathExists(setup_.RestoreFileName());}

   /**
    * @brief Gives the name of the output folder used in this run.
    * @return The path to the folder.
    */
   inline std::string OutputFolderName() const {return output_folder_name_;}

   // Output and Restart functions
   void WriteTimestepFile(std::vector<double> const& timesteps_on_finest_level) const;

   bool WriteFullOutput(double const timestep, const bool force_output = false);
   void WriteDebugFile(unsigned int const debug_key) const;

   void WriteRestartFile(double const timestep, bool const force_output = false);
   double RestoreSimulationFromSnapshot();
};

#endif // INPUT_OUTPUT_MANAGER_H
