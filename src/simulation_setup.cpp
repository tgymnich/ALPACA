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
#include "simulation_setup.h"

#include <algorithm>
#include <numeric>
#include <mpi.h>
#include <type_traits>
#include "materials/material_manager.h"
#include "topology/id_information.h"
#include "helper_functions.h"
#include "prime_states/prime_state_handler_setup.h"

/**
 * @brief Default constructor to create a simulation setup from an input file.
 * @param input_file_name The path to the inputfile.
 */
SimulationSetup::SimulationSetup(const std::string input_file_name) :
   input_file_name_(input_file_name),
   input_file_parser_(input_file_name),
   initial_condition_(input_file_parser_),
   unit_handler_(input_file_parser_),
   output_format_(input_file_parser_.ReadOutputFormat()),
   output_times_(DetermineOutputTimes(input_file_parser_)),
   restore_mode_(input_file_parser_.ReadRestoreMode()),
   restore_file_name_(input_file_parser_.ReadRestoreFileName()),
   snapshot_interval_(input_file_parser_.ReadSnapshotInterval()),
   snapshots_to_keep_(input_file_parser_.ReadSnapshotsToKeep()),
   snapshot_timestamps_(DetermineSnapshotTimes(input_file_parser_)),
   start_time_(unit_handler_.NonDimensionalizeTime(input_file_parser_.ReadStartTime())),
   end_time_(unit_handler_.NonDimensionalizeTime(input_file_parser_.ReadEndTime())),
   time_naming_factor_(input_file_parser_.ReadTimeNamingFactor()),
   maximum_level_(input_file_parser_.ReadMaximumLevel()),
   block_size_on_level_zero_(unit_handler_.NonDimensionalizeLength(input_file_parser_.ReadBlockSize())),
   x_number_of_level_zero_blocks_(input_file_parser_.ReadNumberOfBlocksX()),
   y_number_of_level_zero_blocks_(CC::DIM() != Dimension::One   ? input_file_parser_.ReadNumberOfBlocksY() : 1),
   z_number_of_level_zero_blocks_(CC::DIM() == Dimension::Three ? input_file_parser_.ReadNumberOfBlocksZ() : 1),
   fluid_boundary_conditions_(input_file_parser_.ReadFluidBoundaryConditions()),
   levelset_boundary_conditions_(input_file_parser_.ReadLevelSetBoundaryConditions()),
   active_periodic_locations_(input_file_parser_.ReadActivePeriodicLocations()),
   user_reference_epsilon_( input_file_parser_.ReadReferenceMultiresolutionEpsilon() ),
   user_reference_level_( input_file_parser_.ReadReferenceMultiresolutionLevel() ),
   cfl_number_(input_file_parser_.ReadCFLnumber()),
   gravity_(unit_handler_.NonDimensionalizeGravity(input_file_parser_.ReadGravity())),
   surface_tension_coefficients_(unit_handler_.NonDimensionalizeSurfaceTensionCoefficients(input_file_parser_.ReadSurfaceTensionCoefficients())),
   user_input_fluids_ (unit_handler_.NonDimensionalizeMaterialProperties(input_file_parser_.ReadParameterOfAllFluids())),
   treated_fluids_(MapUserInputToMaterialType(user_input_fluids_)),
   fixed_boundary_conservatives_(ComputeFixedValueBoundaryConservatives(input_file_parser_)),
   fixed_boundary_prime_states_(ComputeFixedValueBoundaryPrimeStates()),
   logger_(LogWriter::Instance())
{
   // Sanity Checks - Do not claim completeness.
   int rank = -1;
   MPI_Comm_rank(MPI_COMM_WORLD,&rank);
   if( rank == 0 ) {
      if( end_time_ < 0 ) {throw std::invalid_argument("End time must not be negative");}
      if( maximum_level_ > CC::AMNL() - 1u ) {throw std::invalid_argument("Maximum level must not exceed 13");}
      if( x_number_of_level_zero_blocks_ + y_number_of_level_zero_blocks_ + z_number_of_level_zero_blocks_ == 0 ) {
         throw std::invalid_argument("At least one Block needs to be present on level zero");
      }
      if( block_size_on_level_zero_ <= 0 ) {throw std::invalid_argument("Block sizes must be greater zero");}
      if( cfl_number_ < 0 ) {throw std::invalid_argument("CFL-Number must be larger zero");}
      if( x_number_of_level_zero_blocks_ > 128 || y_number_of_level_zero_blocks_ > 128 || z_number_of_level_zero_blocks_ > 128 ) {
         throw std::invalid_argument("Block number on level zero 'Block Ratio' must not exceed 128");
      }
   }

   logger_.LogMessage("CFL number         : " + std::to_string(cfl_number_));
   logger_.LogMessage("Starting time      : 0.0",true,true);
   logger_.LogMessage("End time           : " + std::to_string(unit_handler_.DimensionalizeTime(end_time_)));
   logger_.LogMessage("Time naming factor : " + std::to_string(time_naming_factor_));
   logger_.AddBreakLine(true);

   //NH 2017-02-09: Dirty Hack with casts ...
   std::string fluid_boundaries_string = "Fluid Boundaries      : " + std::to_string(static_cast<typename std::underlying_type<FluidBoundaryType>::type>(fluid_boundary_conditions_[0])) + ","
      + std::to_string(static_cast<typename std::underlying_type<FluidBoundaryType>::type>(fluid_boundary_conditions_[1]));
   std::string levelset_boundaries_string = "Level-Set Boundaries  : " + std::to_string(static_cast<typename std::underlying_type<LevelSetBoundaryType>::type>(levelset_boundary_conditions_[0])) + ","
      + std::to_string(static_cast<typename std::underlying_type<LevelSetBoundaryType>::type>(levelset_boundary_conditions_[1]));
   std::string domain_size_string = "Domain Size                  : " + std::to_string(x_number_of_level_zero_blocks_ * unit_handler_.DimensionalizeLength(block_size_on_level_zero_));
   std::string blocks_on_level_zero = "Number of blocks on Level 0  : " + std::to_string(x_number_of_level_zero_blocks_);
   std::string resolution_on_level_zero_string = "Resolution on Level Zero     : " + std::to_string(x_number_of_level_zero_blocks_ * CC::ICX());
   std::string resolution_on_level_max_string  = "Resolution on Maximum Level  : " + std::to_string(x_number_of_level_zero_blocks_ * CC::ICX() * (1 << maximum_level_));

   if constexpr( CC::DIM() != Dimension::One ) {
      fluid_boundaries_string.append("," + std::to_string(static_cast<typename std::underlying_type<FluidBoundaryType>::type>(fluid_boundary_conditions_[2])) +
                                        "," + std::to_string(static_cast<typename std::underlying_type<FluidBoundaryType>::type>(fluid_boundary_conditions_[3]))  );
      levelset_boundaries_string.append("," + std::to_string(static_cast<typename std::underlying_type<LevelSetBoundaryType>::type>(levelset_boundary_conditions_[2])) +
                                           "," + std::to_string(static_cast<typename std::underlying_type<LevelSetBoundaryType>::type>(levelset_boundary_conditions_[3]))  );
      domain_size_string.append("x" + std::to_string(y_number_of_level_zero_blocks_ * unit_handler_.DimensionalizeLength(block_size_on_level_zero_)) );
      blocks_on_level_zero.append("x" + std::to_string(y_number_of_level_zero_blocks_));
      resolution_on_level_zero_string.append("x" + std::to_string(y_number_of_level_zero_blocks_ * CC::ICY()));
      resolution_on_level_max_string.append("x" + std::to_string(y_number_of_level_zero_blocks_ * CC::ICY() * (1 << maximum_level_)));
   }

   if constexpr( CC::DIM() == Dimension::Three ) {
      fluid_boundaries_string.append( "," + std::to_string(static_cast<typename std::underlying_type<FluidBoundaryType>::type>(fluid_boundary_conditions_[4])) +
                                         "," + std::to_string(static_cast<typename std::underlying_type<FluidBoundaryType>::type>(fluid_boundary_conditions_[5]))  );
      levelset_boundaries_string.append( "," + std::to_string(static_cast<typename std::underlying_type<LevelSetBoundaryType>::type>(levelset_boundary_conditions_[4])) +
                                            "," + std::to_string(static_cast<typename std::underlying_type<LevelSetBoundaryType>::type>(levelset_boundary_conditions_[5]))  );
      domain_size_string.append("x" + std::to_string(z_number_of_level_zero_blocks_ * unit_handler_.DimensionalizeLength(block_size_on_level_zero_)) );
      blocks_on_level_zero.append("x" + std::to_string(z_number_of_level_zero_blocks_));
      resolution_on_level_zero_string.append("x" + std::to_string(z_number_of_level_zero_blocks_ * CC::ICZ()));
      resolution_on_level_max_string.append("x" + std::to_string(z_number_of_level_zero_blocks_ * CC::ICZ() * (1 << maximum_level_)));
   }
   blocks_on_level_zero.append(" Blocks ");
   resolution_on_level_zero_string.append(" Internal Cells");
   resolution_on_level_max_string.append(" Internal Cells");

   logger_.LogMessage(domain_size_string);
   logger_.LogMessage("Internal cells per Block     : " + std::to_string(CC::ICX()));
   logger_.LogMessage("Maximum Level                : " + std::to_string(maximum_level_));

   logger_.LogMessage(resolution_on_level_zero_string);
   logger_.LogMessage( "Cell Size on Level Zero      : " + ToScientificNotationString( unit_handler_.DimensionalizeLength( block_size_on_level_zero_ ) / ( CC::ICX() ), 9) );
   logger_.LogMessage(resolution_on_level_max_string);
   //choose direction with maximum cell number
   logger_.LogMessage( "Cell Size on Maximum Level   : " + ToScientificNotationString( unit_handler_.DimensionalizeLength( block_size_on_level_zero_ ) / ( CC::ICX() * ( 1 << maximum_level_ ) ), 9 ) );

   logger_.AddBreakLine(true);

   logger_.LogMessage(fluid_boundaries_string);
   logger_.LogMessage(levelset_boundaries_string);
   //NH 2017-02-09 Dirty Hack with casting ...
   std::string           fluids;
   std::string base_name("Fluid");
   int i = 0;
   for( const auto& fluid : FluidDataForMaterialManager() ) {
      ++i;
      fluids = base_name + std::to_string(i);
      logger_.LogMessage(fluids + "                      : " + std::to_string(static_cast<int>(std::get<0>(fluid))) + " -> " + std::to_string(static_cast<int>(std::get<1>(fluid))));
      for( auto const& [parameter, value] : std::get<2>(fluid)) {
         logger_.LogMessage(parameter + std::string( 28 - parameter.size(), ' ' ) + ": " + ToScientificNotationString( value, 9 ) );
      }
   }
   //has to be adjusted for multilevelset
   logger_.LogMessage( "surface tension coefficient : "
                     + ToScientificNotationString( unit_handler_.DimensionalizeSurfaceTensionCoefficients( surface_tension_coefficients_ )[0], 9 ) );
   // Gruneisen calculation
   if constexpr( CC::GruneisenDensityDependent() ) {
      logger_.LogMessage( "Gruneisen calculated density dependent!" );
   } else {
      logger_.LogMessage( "Gruneisen is material constant!" );
   }
   logger_.AddBreakLine(true);
#ifdef HILBERT
   logger_.LogMessage("Load balancing      : Hilbert-Curve");
#else
   logger_.LogMessage("Load Balancing      : Z-Curve");
#endif

   int number_of_ranks = -1;
   MPI_Comm_size(MPI_COMM_WORLD,&number_of_ranks);
   logger_.LogMessage("Number of MPI ranks : " + std::to_string(number_of_ranks));
   logger_.AddBreakLine(true);
}

/**
 * @brief Gives the output times to be used during the simulation in treated form.
 * @param parser Inputfile parser to get the user inputs.
 * @return List of time instances at which output should be written.
 */
std::vector<double> SimulationSetup::DetermineOutputTimes(const InputFileParser& parser) const {

   const std::string output_type = parser.ReadOutputTimesType();
   std::vector<double> output_times;

   const double start_time = unit_handler_.NonDimensionalizeTime(parser.ReadStartTime());
   const double end_time   = unit_handler_.NonDimensionalizeTime(parser.ReadEndTime());

   if( !(output_type.compare("Interval")) ) {
      const double period = unit_handler_.NonDimensionalizeTime(parser.ReadOutputPeriod());
      if( period <= 0.0 ) {
         throw std::invalid_argument("Output Period must not be smaller equal zero");
      }
      // fill output_times with timestamps generated from given interval
      size_t number_of_timesteps = size_t(std::ceil(std::nextafter((end_time - start_time) / period, 0)));
      if( number_of_timesteps != 0 ) {
         number_of_timesteps -= 1;
      }
      output_times.resize(number_of_timesteps); // make sure end_time is not included
      double time = start_time;
      std::generate(output_times.begin(), output_times.end(), [&time,&period]() {return (time += period);});
   } else if( !output_type.compare("Timestamps") ) {
      output_times = parser.ReadOutputTimestamps();
      // local reference for lambda capture
      const UnitHandler& unit_handler = unit_handler_;
      std::for_each(output_times.begin(), output_times.end(), [&unit_handler](double& timestamp) {timestamp = unit_handler.NonDimensionalizeTime(timestamp);});
      output_times.erase(std::remove_if( output_times.begin(), output_times.end(), [&start_time](const double& timestamp) {return timestamp <= start_time;} ), output_times.end());
      output_times.erase(std::remove_if( output_times.begin(), output_times.end(),   [&end_time](const double& timestamp) {return timestamp >= end_time;  } ), output_times.end());
   } else {
      throw std::invalid_argument("Outputtype does not exist");
   }

   // now add the exact end_time as final timestamp
   output_times.push_back(end_time);

   // final checks
   if( std::any_of(output_times.begin(), output_times.end(), [&start_time](const double& timestamp) {return timestamp < start_time;}) ) {
      throw std::invalid_argument("Output timestamps must be larger than start_time");
   }
   if( std::any_of(output_times.begin(), output_times.end(), [&end_time](const double& timestamp) {return timestamp > end_time;}) ) {
      throw std::invalid_argument("Output timestamps must be smaller than end_time");
   }

   if( !std::is_sorted(output_times.begin(), output_times.end()) ) {
      throw std::invalid_argument("Output timestamps must be ascending");
   }

   return output_times;
}

/**
 * @brief Gives the snapshot timestamps to be used during the simulation in treated form.
 * @param parser Inputfile parser to get the user inputs.
 * @return List of time instances at which restart snapshots should be written.
 */
std::vector<double> SimulationSetup::DetermineSnapshotTimes(const InputFileParser& parser) const {

   std::vector<double> snapshot_times;

   const double start_time = unit_handler_.NonDimensionalizeTime(parser.ReadStartTime());
   const double end_time   = unit_handler_.NonDimensionalizeTime(parser.ReadEndTime());

   snapshot_times = parser.ReadSnapshotTimestamps();
   // local reference for lambda capture
   const UnitHandler& unit_handler = unit_handler_;
   std::for_each(snapshot_times.begin(), snapshot_times.end(), [&unit_handler](double& timestamp) {timestamp = unit_handler.NonDimensionalizeTime(timestamp);});
   snapshot_times.erase(std::remove_if( snapshot_times.begin(), snapshot_times.end(), [&start_time](const double& timestamp) {return timestamp < start_time;} ), snapshot_times.end());
   snapshot_times.erase(std::remove_if( snapshot_times.begin(), snapshot_times.end(),   [&end_time](const double& timestamp) {return timestamp >= end_time;  } ), snapshot_times.end());

   snapshot_times.push_back(end_time);

   // final checks
   if( std::any_of(snapshot_times.begin(), snapshot_times.end(), [&start_time](const double& timestamp) {return timestamp < start_time;}) ) {
      throw std::invalid_argument("Snapshot timestamps must be larger than start_time");
   }
   if( std::any_of(snapshot_times.begin(), snapshot_times.end(), [&end_time](const double& timestamp) {return timestamp > end_time;}) ) {
      throw std::invalid_argument("Snapshot timestamps must be smaller than end_time");
   }

   if( !std::is_sorted(snapshot_times.begin(), snapshot_times.end()) ) {
      throw std::invalid_argument("Snapshot timestamps must be ascending");
   }

   return snapshot_times;
}

/**
 * @brief Reads the prime states given in the inputfile for the fixed value BCs
 * and converts them into conservative states.
 * @param parser Inputfile parser to get the user inputs.
 * @return Array with conservative states in potentially all external boundary conditions.
 */
std::array<std::array<std::array<double, FF::ANOE()>, 2>, 6> SimulationSetup::ComputeFixedValueBoundaryConservatives(const InputFileParser& parser) const {

   //initialize output and temporary buffers; uses maximum number of boundary conditions (6) in 3D as this
   //makes functions below for input-file parser and conversion to primes simple
   std::array<std::array<std::array<double, FF::ANOE()>, 2>, 6> fixed_boundary_values;
   std::array<double, FF::ANOP()> temp;

   //initialize fixed values at boundary with 0
   for( unsigned int i= 0; i<6; ++i ) {
      for( unsigned int j= 0; j<2; ++j ) {
         for( unsigned int k= 0; k<FF::ANOE(); ++k ) {
            fixed_boundary_values[i][j][k] = 0.0;
         }
      }
   }

   //read in fixed value boundary conditions (prime states) from inputfile and convert into conservatives
   if( fluid_boundary_conditions_[0] == FluidBoundaryType::FixedValue ) {
      temp = parser.ReadFixedValueBoundaryCondition(BoundaryLocation::East);
      fixed_boundary_values[0] = ConvertInputToConservatives(temp);
   }

   if( fluid_boundary_conditions_[1] == FluidBoundaryType::FixedValue ) {
      temp = parser.ReadFixedValueBoundaryCondition(BoundaryLocation::West);
      fixed_boundary_values[1] = ConvertInputToConservatives(temp);
   }

   if constexpr( CC::DIM() != Dimension::One ) {
      if( fluid_boundary_conditions_[2] == FluidBoundaryType::FixedValue ) {
         temp = parser.ReadFixedValueBoundaryCondition(BoundaryLocation::North);
         fixed_boundary_values[2] = ConvertInputToConservatives(temp);
      }

      if( fluid_boundary_conditions_[3] == FluidBoundaryType::FixedValue ) {
         temp = parser.ReadFixedValueBoundaryCondition(BoundaryLocation::South);
         fixed_boundary_values[3] = ConvertInputToConservatives(temp);
      }
   }


   if constexpr( CC::DIM() == Dimension::Three ) {
      if( fluid_boundary_conditions_[4] == FluidBoundaryType::FixedValue ) {
         temp = parser.ReadFixedValueBoundaryCondition(BoundaryLocation::Top);
         fixed_boundary_values[4] = ConvertInputToConservatives(temp);
      }

      if( fluid_boundary_conditions_[5] == FluidBoundaryType::FixedValue ) {
         temp = parser.ReadFixedValueBoundaryCondition(BoundaryLocation::Bottom);
         fixed_boundary_values[5] = ConvertInputToConservatives(temp);
      }
   }


   return fixed_boundary_values;
}


/**
 * @brief Gets prime states from inputfile and converts them into conservatives
 * required for the BCs.
 * @param values Array of fixed value BC prime states.
 * @return conservatives.
 */
std::array<std::array<double, FF::ANOE()>, 2> SimulationSetup::ConvertInputToConservatives(const std::array<double, FF::ANOP()> values) const {

   using PrimeStateHandlerConcretization = PrimeStateHandlerSetup::Concretize<prime_state_handler>::type;

   std::array<std::array<double, FF::ANOE()>, 2> conservatives;
   std::array<double, FF::ANOP()> prime_states_nondim;

   // non-dimensionalize prime states first
   for( PrimeState const p : FF::ASOP() ) {
      prime_states_nondim[PTI(p)] = NonDimensionalizeValue(values[PTI(p)], FF::FieldUnit( p ));
   }

   //temporary material manager here necessary to compute fixed value conservatives. Only use here, generally not a good idea!!
   MaterialManager const temporary_material_manager(FluidDataForMaterialManager(), GetSurfaceTensionCoefficients());
   PrimeStateHandlerConcretization const prime_state_handler(temporary_material_manager);

   //set conservatives from inputfile prime states
   for( unsigned int n = 0; n < treated_fluids_.size(); ++n ) {
      prime_state_handler.ConvertPrimeStatesToConservatives(treated_fluids_[n], prime_states_nondim, conservatives[n]);
   }

   return conservatives;
}

/**
 * @brief Converts the conservatives computed from the (reduced) prime states given in the inputfile
 *        for the fixed value BCs into a full set of prime states.
 * @return Array with prime states in potentially all external boundary conditions.
 */
std::array<std::array<std::array<double, FF::ANOP()>, 2>, 6> SimulationSetup::ComputeFixedValueBoundaryPrimeStates() const {
   // approach: computed the full set of prime states from the conservatives that were calculated from
   // the (potentially partial) set of prime states given in the input file.

   using PrimeStateHandlerConcretization = PrimeStateHandlerSetup::Concretize<prime_state_handler>::type;

   std::array<std::array<std::array<double, FF::ANOP()>, 2>, 6> fixed_boundary_values;

   //initialize fixed values at boundary with 0
   for( unsigned int i= 0; i<6; ++i ) {
      for( unsigned int j= 0; j<2; ++j ) {
         for( unsigned int k= 0; k<FF::ANOP(); ++k ) {
            fixed_boundary_values[i][j][k] = 0.0;
         }
      }
   }

   //temporary material manager here necessary to compute fixed value conservatives. Only use here, generally not a good idea!!
   MaterialManager const temporary_material_manager(FluidDataForMaterialManager(), GetSurfaceTensionCoefficients());
   PrimeStateHandlerConcretization const prime_state_handler(temporary_material_manager);

   //set active set of prime states from inputfile prime states
   for( unsigned int s = 0; s < CC::SIDES(); ++s ) {
      if( fluid_boundary_conditions_[s] == FluidBoundaryType::FixedValue ) {
         for( unsigned int n = 0; n < treated_fluids_.size(); ++n ) {
            prime_state_handler.ConvertConservativesToPrimeStates(treated_fluids_[n], fixed_boundary_conservatives_[s][n], fixed_boundary_values[s][n]);
         }
      }
   }

   return fixed_boundary_values;
 }

/**
 * @brief Gives the type of the external fluid boundaries
 * @return The fluid boundary identifer at the edge of the domain.
 */
std::array<FluidBoundaryType, 6> SimulationSetup::GetFluidBoundaryConditions() const {
    return fluid_boundary_conditions_;
}

/**
 * @brief Gives the type of the external levelset boundary conditions.
 * @return The levelset boundary identifer at the edge of the domain.
 */
std::array<LevelSetBoundaryType, 6> SimulationSetup::GetLevelsetBoundaryConditions() const{
    return levelset_boundary_conditions_;
}

/**
 * @brief Gives the type of the external boundary condition at the specified location.
 * @param location Direction of the edge of the domain.
 * @return The boundary identifer at the edge of the domain.
 * @note Not safe. Wrong input leads buffer overflow, hence to segmentation faults.
 */
std::array<std::array<double,FF::ANOE()>, 2> SimulationSetup::FixedValueBoundaryConservatives(const BoundaryLocation location) const {
   return fixed_boundary_conservatives_[LTI(location)];
}

/**
 * @brief Gives the type of the external boundary condition at the specified location.
 * @param location Direction of the edge of the domain.
 * @return The boundary identifer at the edge of the domain.
 * @note Not safe. Wrong input leads buffer overflow, hence to segmentation faults.
 */
std::array<std::array<double,FF::ANOP()>, 2> SimulationSetup::FixedValueBoundaryPrimeStates(const BoundaryLocation location) const {
   return fixed_boundary_prime_states_[LTI(location)];
}

/**
 * @brief Converts the fluid input from the user into input used by the MaterialManager class.
 * @return Tuple holding the Equation of State identifer, the unique material identifier and the fluid data.
 */
std::vector<std::tuple<MaterialName, MaterialName, std::unordered_map<std::string, double>>> SimulationSetup::FluidDataForMaterialManager() const {

   std::vector<std::tuple<MaterialName, MaterialName, std::unordered_map<std::string, double>>> result;
   if( user_input_fluids_.size() != treated_fluids_.size() ) {
      throw std::invalid_argument("Fluid type vectors do not have the same size");
   }

   for( unsigned int i = 0; i < user_input_fluids_.size(); ++i ) {
      result.emplace_back(user_input_fluids_[i].first, treated_fluids_[i], user_input_fluids_[i].second);
   }

   return result;
}

/**
 * @brief Gives a list of all levels which may be present in this simulation run.
 * @return Vector holding the levels in ascending order.
 */
std::vector<unsigned int> SimulationSetup::AllLevels() const {

   std::vector<unsigned int> all_levels(maximum_level_+1); //Level zero need to be counted as well
   std::iota(all_levels.begin(),all_levels.end(),0);
   return all_levels;
}

/**
 * @brief Gives the initial density at the provided location for the given material.
 * @param node_id The id of the node to be initialized.
 * @param material The material in the block to be filled with the returned data.
 * @param initial_values Reference to array holding the resulting density. Indirect return value.
 */
void SimulationSetup::GetInitialPrimeStates(const std::uint64_t node_id, const MaterialName material,
   double (&initial_values)[FF::ANOP()][CC::ICX()][CC::ICY()][CC::ICZ()]) const {
   //block size needs dimension here since initial conditions given in inputfiles have dimensions also
   const double block_size = unit_handler_.DimensionalizeLength( block_size_on_level_zero_ );
   std::array<double,3> const origin = DomainCoordinatesOfId( node_id, DomainSizeOfId(node_id, block_size) );
   double const cell_size = (block_size / double(CC::ICX())) / double(1 << (LevelOfNode(node_id))); // cell_size on level zero divided by 2^level
   initial_condition_.GetInitialPrimeStates( origin, cell_size, material, initial_values );

   // non-dimensionalize obtained prime states
   for( unsigned int i = 0; i < CC::ICX(); ++i ) {
      for( unsigned int j = 0; j < CC::ICY(); ++j ) {
         for( unsigned int k = 0; k < CC::ICZ(); ++k ) {
            for( PrimeState const p : FF::ASOP() ) {
               initial_values[PTI(p)][i][j][k] = NonDimensionalizeValue(initial_values[PTI(p)][i][j][k], FF::FieldUnit( p ));
            }
         }
      }
   }
}

/**
 * @brief Gives the initial levelset for the node with the given id.
 * @param node_id The id of the node to be initialized.
 * @param initial_levelset Indirect return parameter for the determined levelset.
 */
void SimulationSetup::GetInitialLevelset(const std::uint64_t node_id, double (&initial_levelset)[CC::TCX()][CC::TCY()][CC::TCZ()]) const {
   //block size needs dimension here since initial conditions given in inputfiles have dimensions also
   const double block_size = unit_handler_.DimensionalizeLength(block_size_on_level_zero_);
   std::array<double,3> origin = DomainCoordinatesOfId(node_id, DomainSizeOfId(node_id, block_size));
   double cell_size = (block_size / double(CC::ICX())) / double(1 << (LevelOfNode(node_id))); // cell_size on level zero divided by 2^level
   initial_condition_.GetInitialLevelset(origin,cell_size,initial_levelset);
}

/**
 * @brief Gives the initial materials for the node with the given id.
 * @param node_id The id of the node to be initialized.
 * @return A list of the materials present in the considered node.
 */
std::vector<MaterialName> SimulationSetup::GetInitialMaterials(const std::uint64_t node_id) const {
   //block size needs dimension here since initial conditions given in inputfiles have dimensions also
   const double block_size = unit_handler_.DimensionalizeLength(block_size_on_level_zero_);
   std::array<double, 3> origin = DomainCoordinatesOfId(node_id, DomainSizeOfId(node_id, block_size));
   unsigned int level_factor = (1 << (maximum_level_ - LevelOfNode(node_id))); // bit shift is of type "(unsigned?) int"

   const std::vector<bool> material_is_contained = initial_condition_.GetInitialMaterials(origin, unit_handler_.DimensionalizeLength(SmallestPossibleCellSize()), level_factor);

   std::vector<MaterialName> initial_materials;
   for( unsigned int i = 0; i < material_is_contained.size(); ++i ) {
      if( material_is_contained[i] ) {
         initial_materials.push_back(treated_fluids_[i]);
      }
   }
   return initial_materials;
}

/*
 * @brief Translates the non-dimensional value to a value with units.
 * @param value Non-dimensional value.
 * @param dimension Decider which dimensionalization routine has to be used to obtain the correct unit.
 * @return Conservative in unit representation.
 */
double SimulationSetup::NonDimensionalizeValue( double const value, UnitType const dimension ) const {
   switch(dimension) {
      case UnitType::Unitless:
         return value;
      case UnitType::Length:
         return unit_handler_.NonDimensionalizeLength(value);
      case UnitType::Time:
         return unit_handler_.NonDimensionalizeTime(value);
      case UnitType::Density:
         return unit_handler_.NonDimensionalizeDensity(value);
      case UnitType::Velocity:
         return unit_handler_.NonDimensionalizeVelocity(value);
      case UnitType::Momentum:
         return unit_handler_.NonDimensionalizeMomentum(value);
      case UnitType::Pressure:
         return unit_handler_.NonDimensionalizePressure(value);
      case UnitType::Temperature:
         return unit_handler_.NonDimensionalizeTemperature(value);
      case UnitType::Energy:
         return unit_handler_.NonDimensionalizeEnergy(value);
      default:
         throw std::logic_error("UnitType " + std::to_string(static_cast<int>(dimension)) + " is unknown and cannot be non-dimensionalized"); //suppresses Compiler Warning "control reaches end of non-void function [-Wreturn-type]"
   }
}

/*
 * @brief Translates the non-dimensional value to a value with units.
 * @param value Non-dimensional value.
 * @param dimension Decider which dimensionalization routine has to be used to obtain the correct unit.
 * @return Conservative in unit representation.
 */
double SimulationSetup::DimensionalizeValue( double const value, UnitType const dimension ) const {
   switch(dimension) {
      case UnitType::Unitless:
         return value;
      case UnitType::Length:
         return unit_handler_.DimensionalizeLength(value);
      case UnitType::Time:
         return unit_handler_.DimensionalizeTime(value);
      case UnitType::Density:
         return unit_handler_.DimensionalizeDensity(value);
      case UnitType::Velocity:
         return unit_handler_.DimensionalizeVelocity(value);
      case UnitType::Momentum:
         return unit_handler_.DimensionalizeMomentum(value);
      case UnitType::Pressure:
         return unit_handler_.DimensionalizePressure(value);
      case UnitType::Temperature:
         return unit_handler_.DimensionalizeTemperature(value);
      case UnitType::Energy:
         return unit_handler_.DimensionalizeEnergy(value);
      default:
         throw std::logic_error("UnitType " + std::to_string(static_cast<int>(dimension)) + " is unknown and cannot be dimensionalized"); //suppresses Compiler Warning "control reaches end of non-void function [-Wreturn-type]"
   }
}

/**
 * @brief Translates the non-dimensional length/size value to a value with units.
 * @param rho Unitless length.
 * @return Length/size in unit representation.
 */
double SimulationSetup::DimensionalizeLength(const double length) const {
   return unit_handler_.DimensionalizeLength(length);
}

/**
 * @brief Translates the non-dimensional time value to a value with unitz.
 * @param time Unitless time.
 * @return Time in unit representation.
 */
double SimulationSetup::DimensionalizeTime(const double time) const {
   return unit_handler_.DimensionalizeTime(time);
}

/**
 * @brief Gives the smallest cell size possible in the current simulation (number of internal cells and maximum level).
 * @return .
 */
double SimulationSetup::SmallestPossibleCellSize() const {
   double level_zero_cell_size = block_size_on_level_zero_ / double(CC::ICX()); //ICX is always filled.
   double factor =  (1 << maximum_level_); // bit shift is of type "(unsigned?) int" then converison to double happens.
   return level_zero_cell_size/factor;
}

/**
 * @brief Converts the User Input Material Type, e.g. generic equation of state to a unique material
 *        identifier. Does not convert final material classes, e.g. "Water".
 * @param input The material data as retrieved from the inputfile.
 * @return Vector holding one material identifier for each input.
 */
std::vector<MaterialName> SimulationSetup::MapUserInputToMaterialType(const std::vector<std::pair<MaterialName, std::unordered_map<std::string, double>>>& input) {

   std::vector<MaterialName> result;
   MaterialName stiffened_gas_counter = MaterialName::UserStiffenedGasOne;
   MaterialName stiffened_gas_limited_counter = MaterialName::UserStiffenedGasSafeOne;
   MaterialName stiffened_gas_complete_safe_counter = MaterialName::UserStiffenedGasCompleteSafeOne;
   MaterialName waterlike_fluid_counter = MaterialName::UserWaterlikeFluidOne;
   MaterialName noble_abel_stiffened_gas_counter = static_cast<MaterialName>( MTI(MaterialName::NobleAbelStiffenedGas) + 1 );
   MaterialName material;

   for( const auto& fluid : input ) {
      switch(fluid.first) {
         case MaterialName::StiffenedGas : {
            switch (stiffened_gas_counter) {
               case MaterialName::UserStiffenedGasOne: {
                  material = stiffened_gas_counter;
                  stiffened_gas_counter = MaterialName::UserStiffenedGasTwo;
               }
                  break;
               case MaterialName::UserStiffenedGasTwo: {
                  material = stiffened_gas_counter;
                  stiffened_gas_counter = MaterialName::UserStiffenedGasThree;
               }
                  break;
               case MaterialName::UserStiffenedGasThree: {
                  material = stiffened_gas_counter;
                  stiffened_gas_counter = MaterialName::StiffenedGasOutOfBounds;
               }
                  break;
               default: {
                  throw std::invalid_argument("Out of Stiffened Gasses - Please use less User Defined Ones");
               }
                  break;
            }
         }
            break;
         case MaterialName::StiffenedGasSafe : {
            switch (stiffened_gas_limited_counter) {
               case MaterialName::UserStiffenedGasSafeOne: {
                  material = stiffened_gas_limited_counter;
                  stiffened_gas_limited_counter = MaterialName::UserStiffenedGasSafeTwo;
               }
                  break;
               case MaterialName::UserStiffenedGasSafeTwo: {
                  material = stiffened_gas_limited_counter;
                  stiffened_gas_limited_counter = MaterialName::UserStiffenedGasSafeThree;
               }
                  break;
               case MaterialName::UserStiffenedGasSafeThree: {
                  material = stiffened_gas_limited_counter;
                  stiffened_gas_limited_counter = MaterialName::StiffenedGasSafeOutOfBounds;
               }
                  break;
               default: {
                  throw std::invalid_argument("Out of limited Stiffened Gasses - Please use less User Defined Ones");
               }
                  break;
            }
         }
            break;
         case MaterialName::StiffenedGasCompleteSafe : {
            switch (stiffened_gas_complete_safe_counter) {
               case MaterialName::UserStiffenedGasCompleteSafeOne: {
                  material = stiffened_gas_complete_safe_counter;
                  stiffened_gas_complete_safe_counter = MaterialName::UserStiffenedGasCompleteSafeTwo;
               }
                  break;
               case MaterialName::UserStiffenedGasCompleteSafeTwo: {
                  material = stiffened_gas_complete_safe_counter;
                  stiffened_gas_complete_safe_counter = MaterialName::UserStiffenedGasCompleteSafeThree;
               }
                  break;
               case MaterialName::UserStiffenedGasCompleteSafeThree: {
                  material = stiffened_gas_complete_safe_counter;
                  stiffened_gas_complete_safe_counter = MaterialName::StiffenedGasCompleteSafeOutOfBounds;
               }
                  break;
               default: {
                  throw std::invalid_argument("Out of complete Stiffened Gasses - Please use less User Defined Ones");
               }
                  break;
            }
         }
            break;
         case MaterialName::WaterlikeFluid : {
            switch (waterlike_fluid_counter) {
               case MaterialName::UserWaterlikeFluidOne: {
                  material = waterlike_fluid_counter;
                  waterlike_fluid_counter = MaterialName::UserWaterlikeFluidTwo;
               }
                  break;
               case MaterialName::UserWaterlikeFluidTwo: {
                  material = waterlike_fluid_counter;
                  waterlike_fluid_counter = MaterialName::UserWaterlikeFluidThree;
               }
                  break;
               case MaterialName::UserWaterlikeFluidThree: {
                  material = waterlike_fluid_counter;
                  waterlike_fluid_counter = MaterialName::WaterlikeFluidOutOfBounds;
               }
                  break;
               default: {
                  throw std::invalid_argument("Out of Waterlike Fluids - Please use less User Defined Ones");
               }
                  break;
            }
         }
            break;
         case MaterialName::NobleAbelStiffenedGas : {
            material = noble_abel_stiffened_gas_counter;
            noble_abel_stiffened_gas_counter = static_cast<MaterialName>( MTI(noble_abel_stiffened_gas_counter) + 1 );
            if( noble_abel_stiffened_gas_counter == MaterialName::NobleAbelStiffenedGasOutOfBounds ) {
               throw std::invalid_argument("Out of Noble-Abel stiffened gases - Please use less User Defined Ones");
            }
         }
            break;
         default: {
            throw std::logic_error("The provided Material Type Input does not exist");
            break;
         }

      }
      result.emplace_back(material);
   }

   return result;
}
