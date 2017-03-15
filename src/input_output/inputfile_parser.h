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
#ifndef INPUTFILE_PARSER_H
#define INPUTFILE_PARSER_H

#include <string>
#include <array>
#include <unordered_map>
#include <vector>
#include <tinyxml2.h>
#include "boundary_condition/boundary_specifications.h"
#include "enums/restore_mode.h"
#include "materials/material_names.h"
#include "input_output/output_types.h"
#include "fluid_fields_definitions.h"

/**
 * @brief The InputFileParser class provides functions to read out a xml-inputfile. The input is parsed and converted into machine
 *        readable formats (strings, enums, doubles, etc.). The input is neither checked for plausibility nor converted, e.g. into conservative form.
 * @note Uses the TinyXML library by Lee Thomason (www.grinninglizard.com). See respective files for License and Copyright information.
 */
class InputFileParser{
   // The opened xml-file
   tinyxml2::XMLDocument inputfile_;

   // Functions to read single values from xml file
   double ReadDouble( tinyxml2::XMLElement const *node, std::string const tag, double const default_value ) const;
   int ReadInt( tinyxml2::XMLElement const *node, std::string const tag, int const default_value ) const;
   std::string ReadString( tinyxml2::XMLElement const *node ) const;

   // Functions to read specific types from input files
   std::vector<double> ReadTimestamps( tinyxml2::XMLElement const *parent_node ) const;

   FluidBoundaryType SelectFluidBoundaryCondition( std::string const boundary_string ) const;
   LevelSetBoundaryType SelectLevelSetBoundaryCondition( std::string const boundary_string ) const;
   std::string ReadInitialConditionString( std::string const name ) const;

   MaterialName ReadMaterialName( std::string const fluid ) const;
   MaterialName SelectMaterial( std::string const material_name ) const;
   std::unordered_map<std::string, double> ReadMaterialProperties( std::string const fluid ) const;

public:
   InputFileParser() = delete;
   explicit InputFileParser( std::string const& filename );
   ~InputFileParser() = default;
   InputFileParser( InputFileParser const& ) = delete;
   InputFileParser& operator=( InputFileParser const& ) = delete;
   InputFileParser( InputFileParser&& ) = delete;
   InputFileParser& operator=( InputFileParser&& ) = delete;

   double ReadBlockSize() const;
   int ReadNumberOfBlocksX() const;
   int ReadNumberOfBlocksY() const;
   int ReadNumberOfBlocksZ() const;
   std::array<FluidBoundaryType, 6> ReadFluidBoundaryConditions() const;
   std::array<LevelSetBoundaryType, 6> ReadLevelSetBoundaryConditions() const;
   std::array<double,FF::ANOP()> ReadFixedValueBoundaryCondition( const BoundaryLocation location ) const;
   unsigned int ReadActivePeriodicLocations() const;

   std::vector<std::string> ReadInitialConditionOfFluids() const;
   std::vector<std::string> ReadInitialConditionOfLevelSet() const;

   int ReadNumberOfFluids() const;
   std::vector<std::pair<MaterialName, std::unordered_map<std::string, double>>> ReadParameterOfAllFluids() const;
   std::vector<double> ReadSurfaceTensionCoefficients() const;

   std::array<double, 3> ReadGravity() const;

   int ReadMaximumLevel() const;
   int ReadReferenceMultiresolutionLevel() const;
   double ReadReferenceMultiresolutionEpsilon() const;

   RestoreMode ReadRestoreMode() const;
   std::string ReadRestoreFileName() const;
   int ReadSnapshotInterval() const;
   int ReadSnapshotsToKeep() const;
   std::vector<double> ReadSnapshotTimestamps() const;

   double ReadStartTime() const;
   double ReadEndTime() const;
   double ReadCFLnumber() const;

   double ReadReferenceLength() const ;
   double ReadReferenceVelocity() const;
   double ReadReferenceDensity() const;
   double ReadReferenceTemperature() const;

   double ReadMixingThreshold() const;

   bool ReadEnableOutput() const;
   OutputType ReadOutputFormat() const;
   std::string ReadOutputTimesType() const;
   double ReadOutputPeriod() const;
   std::vector<double> ReadOutputTimestamps() const;
   double ReadTimeNamingFactor() const;

};
#endif // INPUTFILE_PARSER_H
