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
#ifndef SIMULATION_SETUP_H
#define SIMULATION_SETUP_H

#include <vector>
#include <array>
#include <tuple>

#include "boundary_condition/boundary_specifications.h"
#include "enums/dimension_definition.h"
#include "enums/restore_mode.h"
#include "materials/material_names.h"
#include "input_output/inputfile_parser.h"
#include "initial_condition.h"
#include "unit_handler.h"
#include "input_output/output_types.h"
#include "log_writer.h"
#include "user_specifications/compile_time_constants.h"
#include "fluid_fields_definitions.h"

/**
 * @brief The SimulationSetup class gives the interface between the user and the simulation kernel. The user input via the xml file is processed
 *        and the given data converted into quantities usable by the kernel. E.g. inputs are given in terms of primary states and are converted to
 *        conservative states.
 */
class SimulationSetup {

   const std::string input_file_name_;

   //Delegating Objects
   const InputFileParser input_file_parser_;
   const InitialCondition initial_condition_;
   const UnitHandler unit_handler_;

   // IO-Parameters
   const OutputType output_format_;
   const std::vector<double> output_times_;

   const RestoreMode restore_mode_;
   const std::string restore_file_name_;
   const int snapshot_interval_;
   const unsigned int snapshots_to_keep_;
   const std::vector<double> snapshot_timestamps_;

   //time parameters
   const double start_time_;
   const double end_time_;
   const double time_naming_factor_;

   //domain parameters
   const unsigned int maximum_level_;
   const double block_size_on_level_zero_;
   const unsigned int x_number_of_level_zero_blocks_;
   const unsigned int y_number_of_level_zero_blocks_;
   const unsigned int z_number_of_level_zero_blocks_;
   const std::array<FluidBoundaryType,6> fluid_boundary_conditions_;
   const std::array<LevelSetBoundaryType,6> levelset_boundary_conditions_;
   const unsigned int active_periodic_locations_;

   // MR parameters
   double const user_reference_epsilon_;
   unsigned int const user_reference_level_;

   const double cfl_number_;

   //fluid parameters
   const std::array<double, 3> gravity_;

   const std::vector<double> surface_tension_coefficients_;

   const std::vector<std::pair<MaterialName, std::unordered_map<std::string, double>>>  user_input_fluids_;
   const std::vector<MaterialName> treated_fluids_;


   //fixed value BC parameters - has to be defined here after the fluid parameters
   const std::array<std::array<std::array<double,FF::ANOE()>, 2>, 6> fixed_boundary_conservatives_;
   const std::array<std::array<std::array<double,FF::ANOP()>, 2>, 6> fixed_boundary_prime_states_;


   //Logger must not be const (otherwise no logbook cannot be appended
   LogWriter& logger_;

   std::vector<double> DetermineOutputTimes(const InputFileParser& parser) const;

   std::vector<double> DetermineSnapshotTimes(const InputFileParser& parser) const;

   //functions for conversion of fixed value BC input parameters into full set of prime states and conservatives
   std::array<std::array<std::array<double, FF::ANOE()>, 2>, 6> ComputeFixedValueBoundaryConservatives(const InputFileParser& parser) const;
   std::array<std::array<std::array<double, FF::ANOP()>, 2>, 6> ComputeFixedValueBoundaryPrimeStates() const;
   std::array<std::array<double, FF::ANOE()>, 2> ConvertInputToConservatives(const std::array<double, FF::ANOP()> values) const;

   double DimensionalizeRho(const double rho) const;
   double DimensionalizeMomentum(const double momentum) const;
   double DimensionalizeEnergy(const double energy) const;

public:
   SimulationSetup() = delete;
   explicit SimulationSetup(const std::string input_file_name);
   ~SimulationSetup() = default;
   SimulationSetup( SimulationSetup const&) = delete;
   SimulationSetup& operator=( SimulationSetup const& ) = delete;
   SimulationSetup( SimulationSetup&& ) = delete;
   SimulationSetup& operator=( SimulationSetup&& ) = delete;

   /**
    * @brief Gives the name of the input file.
    * @return The path to the input file.
    */
   inline std::string InputFileName() const {return input_file_name_;}

   std::array<FluidBoundaryType,6> GetFluidBoundaryConditions() const;
   std::array<LevelSetBoundaryType,6> GetLevelsetBoundaryConditions() const;

   std::array<std::array<double,FF::ANOE()>, 2> FixedValueBoundaryConservatives(const BoundaryLocation location) const;
   std::array<std::array<double,FF::ANOP()>, 2> FixedValueBoundaryPrimeStates(const BoundaryLocation location) const;

   std::vector<std::tuple<MaterialName, MaterialName, std::unordered_map<std::string, double>>> FluidDataForMaterialManager() const;

   std::vector<unsigned int> AllLevels() const;

   void GetInitialPrimeStates( const std::uint64_t node_id, const MaterialName material, double (&initial_values)[FF::ANOP()][CC::ICX()][CC::ICY()][CC::ICZ()]) const;
   void GetInitialLevelset(const std::uint64_t node_id, double (&initial_levelset)[CC::TCX()][CC::TCY()][CC::TCZ()]) const;
   std::vector<MaterialName> GetInitialMaterials(const std::uint64_t node_id) const;

   double NonDimensionalizeValue( double const value, UnitType const dimension ) const;
   double DimensionalizeValue( double const value, UnitType const dimension ) const;
   double DimensionalizeLength(const double length) const;
   double DimensionalizeTime(const double time) const;

   double SmallestPossibleCellSize() const;

   /**
    * @brief Gives a list of all Fluids present in this simulation.
    * @return List of fluids.
    */
   inline std::vector<MaterialName> AllFluids() const {return treated_fluids_;}

   /**
    * @brief Give the number of fluids present in this simulation.
    * @return number of fluids.
    */
   inline unsigned int NumberOfFluids() const {return treated_fluids_.size();}

   /**
    * @brief Gives the time at which the simulation should start.
    * @return Dimensionless time at which the simulation should start.
    */
   inline double GetStartTime() const {return start_time_;}

   /**
    * @brief Gives the time at which the simulation should terminate.
    * @return Dimensionless time at which the simulation should end.
    */
   inline double GetEndTime() const {return end_time_;}

   /**
    * @brief Gives the dimensional multiplication factor for the output file timestamps.
    * @return Dimensional multiplication factor for output time files.
    */
   inline double GetTimeNamingFactor() const {return time_naming_factor_;}

   /**
    * @brief Gives the times at which output should be created.
    * @return List of output times.
    */
   inline std::vector<double> OutputTimes() const {return output_times_;}

   /**
    * @brief Gives whether restart is enabled.
    * @return True if the simulation should be restarted from a file.
    */
   inline RestoreMode GetRestoreMode() const {return restore_mode_;}

   /**
    * @brief Gives the name of the restart file.
    * @return The path to the restart file.
    */
   inline std::string RestoreFileName() const {return restore_file_name_;}

   inline std::vector<double> SnapshotTimestamps() const {return snapshot_timestamps_;}
   inline int SnapshotInterval() const {return snapshot_interval_;}
   inline unsigned int SnapshotsToKeep() const {return snapshots_to_keep_;}

   /**
    * @brief Gives the Courant–Friedrichs–Lewy number.
    * @return CFL number.
    */
   inline double GetCflNumber() const {return cfl_number_;}
   /**
    * @brief Gives the size of a block on level zero.
    * @return Block Size on level 0.
    */
   inline double LevelZeroBlockSize() const {return block_size_on_level_zero_;}
   /**
    * @brief Gives the user epsilon_ref.
    * @return .
    */
   inline double GetUserReferenceEpsilon() const { return user_reference_epsilon_; }

   /**
    * @brief Gives the user level of reference for the multiresolution threshold.
    * @return .
    */
   inline unsigned int GetUserReferenceLevel() const { return user_reference_level_; }

   /**
    * @brief Gives the output format to be used.
    * @return Output format.
    */
   inline OutputType GetOutputFormat() const {return output_format_;}
   /**
    * @brief Gives the maximum level used in this simulation.
    * @return Maximum level.
    */
   inline unsigned int GetMaximumLevel() const {return maximum_level_;}
   /**
    * @brief Gives the number of blocks on level zero in X-Direction.
    * @return Number of blocks in X-Direction on Level zero.
    */
   inline unsigned int GetLevelZeroBlocksX() const {return x_number_of_level_zero_blocks_;}
   /**
    * @brief Gives the number of blocks on level zero in Y-Direction.
    * @return Number of blocks in Y-Direction on Level zero.
    */
   inline unsigned int GetLevelZeroBlocksY() const {return y_number_of_level_zero_blocks_;}
   /**
    * @brief Gives the number of blocks on level zero in Z-Direction.
    * @return Number of blocks in Z-Direction on Level zero.
    */
   inline unsigned int GetLevelZeroBlocksZ() const {return z_number_of_level_zero_blocks_;}
   /**
    * @brief Gives the bitwise active periodic conditions
    * @return bitwise active periodic conditions
    */
   inline unsigned int GetActivePeriodicLocations() const {return active_periodic_locations_;}
   /**
    * @brief Gives the gravity as array 0: X-Component, 1: Y-Component, 2: Z-Component.
    * @return Gravity.
    */
   inline std::array<double, 3> GetGravity() const {return gravity_;}

   /**
    * @brief Gives the surface-tension-coefficients vector.
    * @return Surface-tension-coefficients vector.
    */
   inline std::vector<double> GetSurfaceTensionCoefficients() const {return surface_tension_coefficients_;}

   static std::vector<MaterialName> MapUserInputToMaterialType(const std::vector<std::pair<MaterialName, std::unordered_map<std::string, double>>>& input);

   /**
   * @brief Gives the velocity reference
   * @return velocity reference
   */
   inline double GetVelocityReference() const {
      return unit_handler_.GetVelocityReference();
   }

   /**
    * @brief Gives the density reference
    * @return density reference
    */
   inline double GetDensityReference() const {
      return unit_handler_.GetDensityReference();
   }

   /**
    * @brief Gives the length reference
    * @return length reference
    */
   inline double GetLengthReference() const {
      return unit_handler_.GetLengthReference();
   }

   /**
    * @brief Gives the temperature reference
    * @return temperature reference
    */
   inline double GetTemperatureReference() const {
      return unit_handler_.GetTemperatureReference();
   }
};

#endif // SIMULATION_SETUP_H
