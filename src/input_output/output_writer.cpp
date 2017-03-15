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
#include "output_writer.h"
#include "helper_functions.h"

/**
 * @brief Creates an object to get the simulation data from the RAM to the hard disk.
 * @param flower The tree to read out the data from.
 * @param topology The topology to get information about the MPI status of the simulation.
 * @param setup The setup provides information about user-settings of the simulation.
 * @param output_folder Path to the output folder of the current simulation.
 */
OutputWriter::OutputWriter(Tree const& flower, TopologyManager const& topology, SimulationSetup const& setup, std::string const& output_folder) :
   // Start initializer list
   tree_(flower),
   topology_(topology),
   setup_(setup),
   output_folder_(output_folder),
   logger_(LogWriter::Instance()) {
   /** Empty besides initializer list */
}

/**
 * @brief Writes an output file at the given time.
 * @param output_time The (dimensionalized) time at which the ouput is written.
 */
void OutputWriter::WriteOutputFile(double const output_time) const {
  DoWriteOutputFile(output_time);
  logger_.LogMessage( "Output file written at t = " + ToScientificNotationString( output_time, 9 ), true, true );
}

/**
 * @brief Writes a debug file at the given time.
 * @param output_time The (dimensionalized) time at which the ouput is written.
 */
void OutputWriter::WriteDebugFile(double const debug_time) const {
  DoWriteDebugFile(debug_time);
  logger_.LogMessage( "Debug file written at t = " + ToScientificNotationString( debug_time, 9 ), true, true );
}
