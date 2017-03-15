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
#include <mpi.h>
#include <fenv.h> // Floating-Point raising exceptions.

#include "log_writer.h"
#include "input_output/input_output_manager.h"
#include "simulation_setup.h"
#include "materials/material_manager.h"
#include "topology/topology_manager.h"
#include "topology/tree.h"
#include "communication/communication_manager.h"
#include "modular_algorithm_assembler.h"
#include "multiresolution/multiresolution.h"
#include "halo_manager.h"
#include "boundary_condition/external_halo_manager.h"

/**
 * @brief Starting function of ALPACA, called from the operating system.
 * @param argc Argument count set by your operating system
 * @param argv Input arguments, MPI settings, openMP settings, etc.
 * @return Zero if program finished correctly.
 * @note Please note, due to the usage of openMP #pragmas valgrind and such tools might produce inaccurate measurements.
 * This passes valgrind tests with disabled openMP without errors
 */
int main( int argc, char* argv[] ) {

   MPI_Init( &argc, &argv );
   //Triggers signals on floating point errors, i.e. prohibits quiet NaNs and alike
   feenableexcept( FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW );

   //NH Seperate Scope for MPI.
   {
      LogWriter& logger = LogWriter::Instance();

      // determine the name of the input file (default: inputfile.xml)
      std::string const input_file( argc > 1 ? argv[1] : "inputfile.xml" );

      // determine the name of the executable and write it to the logger
      std::string const  executable_name( argv[0] );
      logger.LogMessage( "Using executable: " + executable_name );

      const SimulationSetup setup( input_file );


      MaterialManager material_manager( setup.FluidDataForMaterialManager(), setup.GetSurfaceTensionCoefficients() );

      TopologyManager topology_manager( setup.GetMaximumLevel(), setup.GetLevelZeroBlocksX(), setup.GetLevelZeroBlocksY(), setup.GetLevelZeroBlocksZ(),
                                        setup.GetActivePeriodicLocations());

      Tree flower( topology_manager, setup.GetMaximumLevel(), setup.LevelZeroBlockSize() );

      ExternalHaloManager external_boundary_manager( setup );
      CommunicationManager communication_manager( topology_manager, setup.GetMaximumLevel() );

      InputOutputManager input_output( setup, topology_manager, flower );
      InternalHaloManager internal_boundary_manager( flower,topology_manager, communication_manager, setup.NumberOfFluids() );
      HaloManager halo_manager( flower, external_boundary_manager, internal_boundary_manager, communication_manager, setup.GetMaximumLevel() );


      Multiresolution multiresolution = InstantiateMultiresolution( setup.GetMaximumLevel(), setup.GetUserReferenceLevel(), setup.GetUserReferenceEpsilon() );

      logger.FlushWelcomeMessage();
      logger.Flush();

      ModularAlgorithmAssembler mr_based_algorithm = ModularAlgorithmAssembler( flower, topology_manager, halo_manager, communication_manager,
                                                                                multiresolution, material_manager, setup, input_output );

      mr_based_algorithm.Initialization();

      mr_based_algorithm.ComputeLoop();

      logger.Flush();
   }

   MPI_Finalize();

   return 0;
}
