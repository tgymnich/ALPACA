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
#include "inputfile_parser.h"

#include <stdexcept>
#include <iostream>
#include <algorithm>
#include "topology/id_periodic_information.h"

/**
 * @brief Standard constructor opening and registering the xml file for later parsings.
 * @param filename XML-File including its path to be opened and read-out.
 */
InputFileParser::InputFileParser(std::string const& filename){

   // Open Config file
   inputfile_.LoadFile( filename.c_str());
   if(inputfile_.FirstChildElement() == NULL){
      throw std::logic_error("Error parsing the XML inputfile file");
   }

}

/**
 * @brief Reads out a numeric value from an XML node, treats and converts it into a double value.
 * @param node The XML node holding the desired information.
 * @param tag Name of the variable to be read-out
 * @param default_value Emergency return value.
 * @return The read-out and converted value.
 */
double InputFileParser::ReadDouble( tinyxml2::XMLElement const *node, std::string const tag, double const default_value ) const {

   double value;
   if(node->QueryDoubleText( &value) != tinyxml2::XML_SUCCESS){
      std::cerr << "Error while reading user defined input argument" << tag << "' set to default value: " << default_value << std::endl;
      return default_value;
   } else {
      return (double)value;
   }

}

/**
 * @brief Reads out a numeric value from an XML node, treats and converts it into an int value.
 * @param node The XML node holding the desired information.
 * @param tag Name of the variable to be read-out.
 * @param default_value Emergency return value.
 * @return The read-out and converted value.
 */
int InputFileParser::ReadInt( tinyxml2::XMLElement const *node, std::string const tag, int const default_value ) const {

   int value;
   if(node->QueryIntText(&value) != tinyxml2::XML_SUCCESS){
      std::cerr << "Error while reading argument (int)" << std::endl;
      std::cerr << "'" << tag << "' set to default value: " << default_value << std::endl;
      return default_value;
   } else {
      return value;
   }

}

/**
 * @brief Reads out a string value from an XML node.
 * @param node The XML node holding the desired information.
 * @return Found String.
 */
std::string InputFileParser::ReadString(tinyxml2::XMLElement const *node) const{

   std::string value = node->GetText();

   if(value.empty() || !value.compare("")){
      std::string const nodename = node->Name();
      throw std::logic_error("ERROR in inputfile file in Node " + nodename);
   }
   return value;
}

/**
 * @brief Reads the full time stamp list for the standard output.
 * @param parent_node The parent node, where all time stamps can be found.
 * @return Vector holding all timestamps (order is the same as the time stamp indices).
 *
 * @note the time stamps are specified as child nodes under the parent node in ascending order, i.e.
 *       <parentNode>
 *          <t1> </t1>
 *          <t2> </t2>
 *          ...
 *       </parentNode>
 */
std::vector<double> InputFileParser::ReadTimestamps(tinyxml2::XMLElement const *parent_node) const {

   std::string timestamp_name;
   std::string base_name("ts");

   std::vector<double> timestamps;

   bool keep_reading = true;
   int index = 1;
   while(keep_reading) {
      timestamp_name = base_name + std::to_string(index++);
      tinyxml2::XMLElement const* node = parent_node->FirstChildElement(timestamp_name.c_str());
      if(node == nullptr) {
         keep_reading = false;
      } else {
         timestamps.emplace_back(ReadDouble(node,timestamp_name.c_str(),0.0));
      }
   }

   return timestamps;
}

/**
 * @brief Reads out the initial condition of the simulation. Gives the initial condition as string (equation form).
 *        No conversion to machine-readbale inputs is performed.
 * @param name The name of the object (first/second/n-th fluid, phi, ..) whose initial condition is to be retrieved.
 * @return Initial Conditions in equation form.
 */
std::string InputFileParser::ReadInitialConditionString(std::string const name) const {
   tinyxml2::XMLElement const *node = inputfile_.FirstChildElement()->FirstChildElement("domain")->FirstChildElement("initialConditions")->FirstChildElement(name.c_str());
   return ReadString(node);
}


/**
 * @brief Reads-out the types of the materials, i. e. Equation of State.
 * @param fluid The name of the considered fluid as in XML file (fluid1, fluid2, ..., fluidN).
 * @return The material identifer enum associated with the equation of state.
 */
MaterialName InputFileParser::ReadMaterialName(std::string const fluid) const {

   tinyxml2::XMLElement const* node = inputfile_.FirstChildElement()->FirstChildElement("fluids")->FirstChildElement(fluid.c_str())->FirstChildElement("type");
   std::string  material_name = ReadString(node);
   return SelectMaterial(material_name);
}

/**
 * @brief Associates a fluid boundary type to the string describing the boundary condition.
 * @param boundary_string The boundary description as found in the XML file.
 * @return Machine readable fluid boundary definition (type identifier).
 */
FluidBoundaryType InputFileParser::SelectFluidBoundaryCondition(std::string const boundary_string) const {
   // NH 2016-10-24: Switch Statements do not work with String

   //make entire string upper case and erase whitespace
   std::string boundary_string_upper_case = boundary_string;
   boundary_string_upper_case.erase(std::remove_if(boundary_string_upper_case.begin(), boundary_string_upper_case.end(), ::isspace), boundary_string_upper_case.end());
   std::transform(boundary_string_upper_case.begin(), boundary_string_upper_case.end(),boundary_string_upper_case.begin(), ::toupper);

   if(boundary_string_upper_case == "ZEROGRADIENT") {return FluidBoundaryType::ZeroGradient;}
   else if(boundary_string_upper_case == "SYMMETRY") {return FluidBoundaryType::Symmetry;}
   else if(boundary_string_upper_case == "FIXEDVALUE") {return FluidBoundaryType::FixedValue;}
   else if(boundary_string_upper_case == "WALL") {return FluidBoundaryType::Wall;}
   else if(boundary_string_upper_case == "PERIODIC") {return FluidBoundaryType::Periodic;}
   else { throw std::invalid_argument("Fluid boundary condition in input file is undefined / not yet implemented");}
}

/**
 * @brief Associates a level-set boundary type to the string describing the boundary condition.
 * @param boundary_string The boundary description as found in the XML file.
 * @return Machine readable level-set boundary definition (type identifier).
 */
LevelSetBoundaryType InputFileParser::SelectLevelSetBoundaryCondition(std::string const boundary_string) const {

   //make entire string upper case and erase whitespace
   std::string boundary_string_upper_case = boundary_string;
   boundary_string_upper_case.erase(std::remove_if(boundary_string_upper_case.begin(), boundary_string_upper_case.end(), ::isspace), boundary_string_upper_case.end());
   std::transform(boundary_string_upper_case.begin(), boundary_string_upper_case.end(),boundary_string_upper_case.begin(), ::toupper);

   if (boundary_string_upper_case == "ZEROGRADIENT") {return LevelSetBoundaryType::ZeroGradient;}
   else if(boundary_string_upper_case == "SYMMETRY") {return LevelSetBoundaryType::Symmetry;}
   else if(boundary_string_upper_case == "PERIODIC") {return LevelSetBoundaryType::Periodic;}
   else { throw std::invalid_argument("Level-set boundary condition in input file is undefined / not yet implemented");}
}

/**
 * @brief Associates a material identifer to the string describing the material.
 * @param material_name The material description as found in the XML file.
 * @return Machine readable material (type identifier).
 */
MaterialName InputFileParser::SelectMaterial(std::string const material_name) const {

   //make entire string upper case and erase whitespace
   std::string material_name_upper_case = material_name;
   material_name_upper_case.erase(std::remove_if(material_name_upper_case.begin(), material_name_upper_case.end(), ::isspace), material_name_upper_case.end());
   std::transform(material_name_upper_case.begin(), material_name_upper_case.end(),material_name_upper_case.begin(), ::toupper);

   if (material_name_upper_case == "STIFFENEDGAS") {return MaterialName::StiffenedGas;}
   else if(material_name_upper_case == "STIFFENEDGASSAFE") {return MaterialName::StiffenedGasSafe;}
   else if(material_name_upper_case == "STIFFENEDGASCOMPLETESAFE") {return MaterialName::StiffenedGasCompleteSafe;}
   else if(material_name_upper_case == "WATERLIKEFLUID") {return MaterialName::WaterlikeFluid;}
   else if(material_name_upper_case == "NOBLEABELSTIFFENEDGAS") {return MaterialName::NobleAbelStiffenedGas;}
   else { throw std::invalid_argument("Material name in input file is undefined / not yet implemented");}
}

/**
 * @brief Reads out the fluid properties, i.e. material parameters.
 * @param fluid The name of the considered fluid as in the xml file (fluid1, fluid2, ..., fluidN).
 * @return The obtained values, usually:
 *         0: gamma Default: 0.0,
 *         1: A Default: 0.0,
 *         2: B Default: 0.0,
 *         3: C Default: 0.0,
 *         4: rho0 Default: 0.0,
 *         5: specific gas constant Default 0.0,
 *         6: thermal conductivity Default 0.0,
 *         7: dynamic shear viscosity Default 0.0,
 *         8: dynamic bulk viscosity Default 0.0,
 */
std::unordered_map<std::string, double> InputFileParser::ReadMaterialProperties( std::string const fluid ) const {

   std::unordered_map<std::string, double> material_properties;
   tinyxml2::XMLElement const* node = inputfile_.FirstChildElement()->FirstChildElement( "fluids" )->FirstChildElement( fluid.c_str() )->FirstChildElement();
   while( node != nullptr ) {
      std::string const name = node->Name();
      if( name == "viscosity" ) {
         material_properties["dynamicShear"] = ReadDouble( node->FirstChildElement( "dynamicShear" ), "dynamicShear", 0.0 );
         material_properties["dynamicBulk"] = ReadDouble( node->FirstChildElement( "dynamicBulk" ), "dynamicBulk", 0.0 );
      } else if( name != "type" ) {
         material_properties[name] = ReadDouble( node, name, 0.0 );
      }
      node = node->NextSiblingElement();
   }

   return material_properties;
}

/**
 * @brief Reads out the size of blocks on level zero.
 * @return Size of the blocks. Default: 1.
 */
double InputFileParser::ReadBlockSize() const{
   tinyxml2::XMLElement const* node = inputfile_.FirstChildElement()->FirstChildElement("domain")->FirstChildElement("blockSize");
   return ReadDouble(node, "Block Size", 1);
}

/**
 * @brief Reads out how many blocks on level zero will be present in x-direction.
 * @return Number of blocks in x-direction. Default: 1.
 */
int InputFileParser::ReadNumberOfBlocksX() const {
   tinyxml2::XMLElement const* node = inputfile_.FirstChildElement()->FirstChildElement("domain")->FirstChildElement("blockRatio")->FirstChildElement("x");
   return ReadInt(node, "Block Ratio x", 1);
}

/**
 * @brief Reads out how many blocks on level zero will be present in y-direction.
 * @return Number of blocks in y-direction. Default: 1.
 */
int InputFileParser::ReadNumberOfBlocksY() const {
   tinyxml2::XMLElement const* node = inputfile_.FirstChildElement()->FirstChildElement("domain")->FirstChildElement("blockRatio")->FirstChildElement("y");
   return ReadInt(node, "Block Ratio y", 1);
}

/**
 * @brief Reads out how many blocks on level zero will be present in z-direction.
 * @return Number of blocks in z-direction. Default: 1.
 */
int InputFileParser::ReadNumberOfBlocksZ() const{
   tinyxml2::XMLElement const* node = inputfile_.FirstChildElement()->FirstChildElement("domain")->FirstChildElement("blockRatio")->FirstChildElement("z");
   return ReadInt(node, "Block Ratio z", 1);
}

/**
 * @brief Reads out the type of the external fluid boundary conditions and converts them into a machine-readable FluidBoundaryType enum.
 * @return Boundary types at the respective edges of the domain. 0: West, 1: East, 2: South, 3: North, 4: Bottom, 5: Top.
 */
std::array<FluidBoundaryType, 6> InputFileParser::ReadFluidBoundaryConditions() const {

   std::array<FluidBoundaryType, 6> fluid_boundaries;

   tinyxml2::XMLElement const* const node0 = inputfile_.FirstChildElement()->FirstChildElement( "domain" )->FirstChildElement( "boundaryConditions" )->FirstChildElement( "fluid" )->FirstChildElement( "east" );
   fluid_boundaries[LTI( BoundaryLocation::East )]   = SelectFluidBoundaryCondition( ReadString( node0 ) );
   tinyxml2::XMLElement const* const node1 = inputfile_.FirstChildElement()->FirstChildElement( "domain" )->FirstChildElement( "boundaryConditions" )->FirstChildElement( "fluid" )->FirstChildElement( "west" );
   fluid_boundaries[LTI( BoundaryLocation::West )]   =  SelectFluidBoundaryCondition( ReadString( node1 ) );
   tinyxml2::XMLElement const* const node2 = inputfile_.FirstChildElement()->FirstChildElement( "domain" )->FirstChildElement( "boundaryConditions" )->FirstChildElement( "fluid" )->FirstChildElement( "north" );
   fluid_boundaries[LTI( BoundaryLocation::North )]  = SelectFluidBoundaryCondition( ReadString( node2 ) );
   tinyxml2::XMLElement const* const node3 = inputfile_.FirstChildElement()->FirstChildElement( "domain" )->FirstChildElement( "boundaryConditions" )->FirstChildElement( "fluid" )->FirstChildElement( "south" );
   fluid_boundaries[LTI( BoundaryLocation::South )]  = SelectFluidBoundaryCondition( ReadString( node3 ) );
   tinyxml2::XMLElement const* const node4 = inputfile_.FirstChildElement()->FirstChildElement( "domain" )->FirstChildElement( "boundaryConditions" )->FirstChildElement( "fluid" )->FirstChildElement( "top" );
   fluid_boundaries[LTI( BoundaryLocation::Top )]    =  SelectFluidBoundaryCondition( ReadString( node4 ) );
   tinyxml2::XMLElement const* const node5 = inputfile_.FirstChildElement()->FirstChildElement( "domain" )->FirstChildElement( "boundaryConditions" )->FirstChildElement( "fluid" )->FirstChildElement( "bottom" );
   fluid_boundaries[LTI( BoundaryLocation::Bottom )] = SelectFluidBoundaryCondition( ReadString( node5 ) );

   return fluid_boundaries;
}

/**
 * @brief Reads out the type of the external levelset boundary conditions and converts them into a machine-readable LevelSetBoundaryType enum.
 * @return Boundary types at the respective edges of the domain. 0: West, 1: East, 2: South, 3: North, 4: Bottom, 5: Top.
 */
std::array<LevelSetBoundaryType, 6> InputFileParser::ReadLevelSetBoundaryConditions() const {

   std::array<LevelSetBoundaryType, 6> levelset_boundaries;

   tinyxml2::XMLElement const* const node0 = inputfile_.FirstChildElement()->FirstChildElement( "domain" )->FirstChildElement( "boundaryConditions" )->FirstChildElement( "levelSet" )->FirstChildElement( "east" );
   levelset_boundaries[LTI( BoundaryLocation::East )]   = SelectLevelSetBoundaryCondition( ReadString( node0 ) );
   tinyxml2::XMLElement const* const node1 = inputfile_.FirstChildElement()->FirstChildElement( "domain" )->FirstChildElement( "boundaryConditions" )->FirstChildElement( "levelSet" )->FirstChildElement( "west" );
   levelset_boundaries[LTI( BoundaryLocation::West )]   = SelectLevelSetBoundaryCondition( ReadString( node1 ) );
   tinyxml2::XMLElement const* const node2 = inputfile_.FirstChildElement()->FirstChildElement( "domain" )->FirstChildElement( "boundaryConditions" )->FirstChildElement( "levelSet" )->FirstChildElement( "north" );
   levelset_boundaries[LTI( BoundaryLocation::North )]  = SelectLevelSetBoundaryCondition( ReadString( node2 ) );
   tinyxml2::XMLElement const* const node3 = inputfile_.FirstChildElement()->FirstChildElement( "domain" )->FirstChildElement( "boundaryConditions" )->FirstChildElement( "levelSet" )->FirstChildElement( "south" );
   levelset_boundaries[LTI( BoundaryLocation::South )]  = SelectLevelSetBoundaryCondition( ReadString( node3 ) );
   tinyxml2::XMLElement const* const node4 = inputfile_.FirstChildElement()->FirstChildElement( "domain" )->FirstChildElement( "boundaryConditions" )->FirstChildElement( "levelSet" )->FirstChildElement( "top" );
   levelset_boundaries[LTI( BoundaryLocation::Top )]    = SelectLevelSetBoundaryCondition( ReadString( node4 ) );
   tinyxml2::XMLElement const* const node5 = inputfile_.FirstChildElement()->FirstChildElement( "domain" )->FirstChildElement( "boundaryConditions" )->FirstChildElement( "levelSet" )->FirstChildElement( "bottom" );
   levelset_boundaries[LTI( BoundaryLocation::Bottom )] = SelectLevelSetBoundaryCondition( ReadString( node5 ) );

   return levelset_boundaries;
}

/**
 * @brief Reads out prime states specified in the fixed value boundary condition.
 * @return Boundary locations at the respective edges of the domain. 0: West, 1: East, 2: South, 3: North, 4: Bottom, 5: Top.
 */
std::array<double,FF::ANOP()> InputFileParser::ReadFixedValueBoundaryCondition(BoundaryLocation const location) const {

   std::array<double,FF::ANOP()> fixed_values;
   std::string location_name;

   //select name of BC that needs to be read in from inputfile
   if(location == BoundaryLocation::East){
      location_name = "valuesEast";
   } else if(location == BoundaryLocation::West){
      location_name = "valuesWest";
   } else if(location == BoundaryLocation::North){
      location_name = "valuesNorth";
   } else if(location == BoundaryLocation::South){
      location_name = "valuesSouth";
   } else if(location == BoundaryLocation::Top){
      location_name = "valuesTop";
   } else if(location == BoundaryLocation::Bottom){
      location_name = "valuesBottom";
   }

   //read values from inputfile
   for( PrimeState const prime : FF::ASOP() ) {
      // ATTENTION: here we use the OUTPUT name of the prime state, not the input name (TODO-19: maybe unify those two names)
      std::string const prime_name(FF::FieldOutputName( prime ));
      // check both input and output name for empty (input name gives the intention, output name the actual used name, see comment above)
      if( !prime_name.empty() && !FF::FieldInputName( prime ).empty() ) {
         tinyxml2::XMLElement const* node = inputfile_.FirstChildElement()->FirstChildElement("domain")->FirstChildElement("boundaryConditions")->FirstChildElement("fluid")->FirstChildElement(location_name.c_str())->FirstChildElement(prime_name.c_str());
         fixed_values[PTI(prime)] = ReadDouble(node, "fixedValue_" + prime_name, 0.0);
      } else {
         fixed_values[PTI(prime)] = 0.0;
      }
   }

   return fixed_values;
}

/**
 * @brief Reads the plane boundary locations that are periodic BCs.
 * @return three bit identifier that gives all active periodic regions:
 *         Bit 1: East-West periodic, Bit 2: North-South periodic, Bit 3: Top-Bottom periodic.
 */
unsigned int InputFileParser::ReadActivePeriodicLocations() const {

   std::array<FluidBoundaryType, 6> const fluid_boundary_conditions( ReadFluidBoundaryConditions() );
   std::array<LevelSetBoundaryType, 6> const levelset_boundary_conditions( ReadLevelSetBoundaryConditions() );
   unsigned int active_periodic_locations = 0;

   // The if-else statements contain the following checks (same holds for y- and z-direction):
   // 1. Check if any levelset or fluid boundary has periodicity in east or west.
   // 2. If so, Check if the opposite side is also periodic.
   // 3. If not, throw error.
   if( ( fluid_boundary_conditions[LTI( BoundaryLocation::East )] == FluidBoundaryType::Periodic ) ||
       ( fluid_boundary_conditions[LTI( BoundaryLocation::West )] == FluidBoundaryType::Periodic ) ||
       ( levelset_boundary_conditions[LTI( BoundaryLocation::East )] == LevelSetBoundaryType::Periodic ) ||
       ( levelset_boundary_conditions[LTI( BoundaryLocation::West )] == LevelSetBoundaryType::Periodic ) ) {
      if( ( fluid_boundary_conditions[LTI( BoundaryLocation::East )] == FluidBoundaryType::Periodic ) &&
          ( fluid_boundary_conditions[LTI( BoundaryLocation::West )] == FluidBoundaryType::Periodic ) &&
          ( levelset_boundary_conditions[LTI( BoundaryLocation::East )] == LevelSetBoundaryType::Periodic ) &&
          ( levelset_boundary_conditions[LTI( BoundaryLocation::West )] == LevelSetBoundaryType::Periodic ) ) {

         active_periodic_locations |= PeriodicBoundariesLocations::EastWest;
      } else {
        throw std::invalid_argument("Incorrect use of EastWest periodic condition, both east and west boundaries from the fluid and levelset must be periodic");
      }
   }

   if( ( fluid_boundary_conditions[LTI( BoundaryLocation::North )] == FluidBoundaryType::Periodic ) ||
       ( fluid_boundary_conditions[LTI( BoundaryLocation::South )] == FluidBoundaryType::Periodic ) ||
       ( levelset_boundary_conditions[LTI( BoundaryLocation::North )] == LevelSetBoundaryType::Periodic ) ||
       ( levelset_boundary_conditions[LTI( BoundaryLocation::South )] == LevelSetBoundaryType::Periodic ) ) {
      if( ( fluid_boundary_conditions[LTI( BoundaryLocation::North )] == FluidBoundaryType::Periodic ) &&
          ( fluid_boundary_conditions[LTI( BoundaryLocation::South )] == FluidBoundaryType::Periodic ) &&
          ( levelset_boundary_conditions[LTI( BoundaryLocation::North )] == LevelSetBoundaryType::Periodic ) &&
          ( levelset_boundary_conditions[LTI( BoundaryLocation::South )] == LevelSetBoundaryType::Periodic ) ) {

         active_periodic_locations |= PeriodicBoundariesLocations::NorthSouth;
      }else{
         throw std::invalid_argument("Incorrect use of NorthSouth periodic condition, both north and south boundaries from the fluid and levelset must be periodic");
      }
   }

   if( ( fluid_boundary_conditions[LTI( BoundaryLocation::Top )] == FluidBoundaryType::Periodic ) ||
       ( fluid_boundary_conditions[LTI( BoundaryLocation::Bottom )] == FluidBoundaryType::Periodic ) ||
       ( levelset_boundary_conditions[LTI( BoundaryLocation::Top )] == LevelSetBoundaryType::Periodic ) ||
       ( levelset_boundary_conditions[LTI( BoundaryLocation::Bottom )] == LevelSetBoundaryType::Periodic ) ) {
      if( ( fluid_boundary_conditions[LTI( BoundaryLocation::Top )] == FluidBoundaryType::Periodic ) &&
          ( fluid_boundary_conditions[LTI( BoundaryLocation::Bottom )] == FluidBoundaryType::Periodic ) &&
          ( levelset_boundary_conditions[LTI( BoundaryLocation::Top )] == LevelSetBoundaryType::Periodic ) &&
          ( levelset_boundary_conditions[LTI( BoundaryLocation::Bottom )] == LevelSetBoundaryType::Periodic ) ) {
         active_periodic_locations |= PeriodicBoundariesLocations::TopBottom;
      } else {
         throw std::invalid_argument("Incorrect use of TopBottom periodic condition, both top and bottom boundaries from the fluid and levelset must be periodic");
      }
   }
   return active_periodic_locations;
}

/**
 * @brief Reads out the initial condition of all fluids specifed in the inputfile.
 * @return Equation form of the initial condition of all fluids.
 */
std::vector<std::string>  InputFileParser::ReadInitialConditionOfFluids() const {

   std::vector<std::string> fluid_initialisation;
   std::string fluid;
   std::string base_name("fluid");

   for(int i = 1; i <= ReadNumberOfFluids(); ++i){
      fluid = base_name + std::to_string(i);
      fluid_initialisation.push_back(ReadInitialConditionString(fluid));
   }

   return fluid_initialisation;
}

/**
 * @brief Reads out the initial condition of the (all) level sets for a two-phase (multi-phase) simulation.
 * @return Equation form of the inital levelset(s).
 */
std::vector<std::string>  InputFileParser::ReadInitialConditionOfLevelSet() const {

   std::vector<std::string> phi_expressions;
   std::string phi;
   std::string name_basic("levelSet");

   phi_expressions.push_back(ReadInitialConditionString("levelSet1"));

   for(int i = 2; i < ReadNumberOfFluids(); ++i){
      phi = name_basic + std::to_string(i);
      phi_expressions.push_back(ReadInitialConditionString(phi));
   }

   return phi_expressions;
}

/**
 * @brief Reads out the number of fluids present in the simulation.
 * @return Number of fluids. Default: 1.
 */
int InputFileParser::ReadNumberOfFluids() const {
   tinyxml2::XMLElement const* node = inputfile_.FirstChildElement()->FirstChildElement("fluids")->FirstChildElement("numberOfFluids");
   return ReadInt(node, "numberOfFluids", 1);
}

/**
 * @brief Reads out the material parameters for all fluids defined in the inputfile.
 * @return Vector of pairs defining a fluid by its material type, i.e. equation of state, and its properties, e.g. gamma, viscosity, etc.
 */
std::vector<std::pair<MaterialName, std::unordered_map<std::string, double>>>  InputFileParser::ReadParameterOfAllFluids() const {

   std::string           fluid;
   std::string base_name("fluid");
   std::vector<std::pair<MaterialName, std::unordered_map<std::string, double>>>  fluids;

   for(int i = 1; i <= ReadNumberOfFluids(); ++i){
      fluid = base_name + std::to_string(i);

      // add pair of fluid data to vector
      fluids.emplace_back( ReadMaterialName(fluid), ReadMaterialProperties(fluid) );
   }

   return fluids;
}

/**
 * @brief Reads out the surface tension coefficients. After it will be adapted for multiphase, the order of the elements should be as follows:
 *        [0] = sigma_12, [1] = sigma_13, [2] = sigma_14 and so forth.
 * @return Surface tension coefficients. Default: 0.0.
 */
std::vector<double> InputFileParser::ReadSurfaceTensionCoefficients() const {
   std::vector<double> surface_tension_coefficients;
   tinyxml2::XMLElement const* node = inputfile_.FirstChildElement()->FirstChildElement("fluids")->FirstChildElement("surfaceTensionCoefficients");
   surface_tension_coefficients.push_back(ReadDouble(node, "surfaceTensionCoefficients", 0.0));
   return surface_tension_coefficients;
}

/**
 * @brief Reads out the gravitational pull as three dimensional array, one for the pull in each cartesian direction.
 * @return Gravity in x-, y- and z-direction. Default: 0 in all directions.
 */
std::array<double, 3> InputFileParser::ReadGravity() const {
   std::array<double, 3> gravity;
   tinyxml2::XMLElement const* node0 = inputfile_.FirstChildElement()->FirstChildElement("sourceTerms")->FirstChildElement("gravity")->FirstChildElement("x");
   gravity[0] = ReadDouble(node0, "gravity x", 0);
   tinyxml2::XMLElement const* node1 = inputfile_.FirstChildElement()->FirstChildElement("sourceTerms")->FirstChildElement("gravity")->FirstChildElement("y");
   gravity[1] = ReadDouble(node1, "gravity y", 0);
   tinyxml2::XMLElement const* node2 = inputfile_.FirstChildElement()->FirstChildElement("sourceTerms")->FirstChildElement("gravity")->FirstChildElement("z");
   gravity[2] = ReadDouble(node2, "gravity z", 0);
   return gravity;
}

/**
 * @brief Reads out the maximum level (of refinement of the simulation).
 * @return Finest level. Default: 0.
 */
int InputFileParser::ReadMaximumLevel() const {
   tinyxml2::XMLElement const* node = inputfile_.FirstChildElement()->FirstChildElement("multiResolution")->FirstChildElement("maximumLevel");
   return ReadInt(node, "maximum level", 0);
}

/**
 * @brief Reads-out epsilon_ref, which defines the threshold for coarsening of levels.
 *        Epsilon_ref is thereby not the finally used threshold directly, but is modified according to Roussel et al. (2003)- "A conservative
 *        fully adaptive multiresolution algorithm for parabolic PDEs".
 * @return Epsilon_ref Reference epsilon. Default: 0.05.
 */
double InputFileParser::ReadReferenceMultiresolutionEpsilon() const {
   tinyxml2::XMLElement const *node = inputfile_.FirstChildElement()->FirstChildElement("multiResolution")->FirstChildElement("refinementCriterion")->FirstChildElement("epsilonReference");
   return ReadDouble(node, "epsilonReference", 0.01);
}

/**
 * @brief Reads out the level on which the reference epsilon (see ReadReferenceMultiresolutionEpsilon()) is enforced.
 * @return Level on which reference epsilon is to be enforced. Default: 1.
 */
int InputFileParser::ReadReferenceMultiresolutionLevel() const {
   tinyxml2::XMLElement const* node = inputfile_.FirstChildElement()->FirstChildElement("multiResolution")->FirstChildElement("refinementCriterion")->FirstChildElement("levelOfEpsilonReference");
   return ReadInt(node, "levelOfEpsilonReference", 1);
}

/**
 * @brief Reads-out whether the simulation is started from a restart file. 0 means "no", 1 means "yes, if file exists", 2 means "yes in any case".
 * @return Identifier of the mode used for restororing. Default: 0 => Disabled.
 */
RestoreMode InputFileParser::ReadRestoreMode() const {
   tinyxml2::XMLElement const* node = inputfile_.FirstChildElement()->FirstChildElement("restart")->FirstChildElement("restoreMode");
   return static_cast<RestoreMode>(ReadInt(node, "restoreMode", 0));
}

/**
 * @brief Reads-out the name of the restart file. Should be "false" if no restart file is to be used.
 * @return Filename string.
 * @note The absolute path needs to be specified, otherwise it used the inputfile folder to find restart file.
 */
std::string InputFileParser::ReadRestoreFileName() const {
   tinyxml2::XMLElement const *node = inputfile_.FirstChildElement()->FirstChildElement("restart")->FirstChildElement("restoreFileName");
   std::string filename = ReadString(node);
   filename.erase(std::remove_if(filename.begin(), filename.end(), ::isspace), filename.end());
   return filename;
}

/**
 * @brief Reads-out the wall clock interval between restart snapshots. Unit is wall seconds.
 * @return The restart snapshot interval in wall seconds.
 */
int InputFileParser::ReadSnapshotInterval() const {
   tinyxml2::XMLElement const* node = inputfile_.FirstChildElement()->FirstChildElement("restart")->FirstChildElement("snapshotInterval");
   return ReadInt(node, "snapshotInterval", 0);
}

/**
 * @brief Reads-out the number of interval-based restart snapshots to keep. Once the number is exceeded, the oldest snapshots are deleted.
 *        Timestamp-based snapshots are always kept.
 * @return The number of snapshots to keep.
 */
int InputFileParser::ReadSnapshotsToKeep() const {
   tinyxml2::XMLElement const* node = inputfile_.FirstChildElement()->FirstChildElement("restart")->FirstChildElement("snapshotsToKeep");
   return ReadInt(node, "snapshotsToKeep", 1);
}

/**
 * @brief Reads-out the timestamps at which restart snapshots should be written. Unit is simulation time.
 * @return A list containing the snapshot timestamps.
 */
std::vector<double> InputFileParser::ReadSnapshotTimestamps() const {
   return ReadTimestamps(inputfile_.FirstChildElement()->FirstChildElement("restart")->FirstChildElement("snapshotTimestamps"));
}


/**
 * @brief Reads-out the time at which the simulation is supposed to start.
 * @return Start time. Default: 0.0.
 */
double InputFileParser::ReadStartTime() const {
   tinyxml2::XMLElement const* node = inputfile_.FirstChildElement()->FirstChildElement("calculation")->FirstChildElement("timeControl")->FirstChildElement("startTime");
   return ReadDouble(node, "startTime", 0.0);
}

/**
 * @brief Reads out the time at which the simulation is supposed to end. I. e. EndTime - StartTime = SimulatedTime.
 * @return End time. Default: 0.0.
 */
double InputFileParser::ReadEndTime() const {
   tinyxml2::XMLElement const* node = inputfile_.FirstChildElement()->FirstChildElement("calculation")->FirstChildElement("timeControl")->FirstChildElement("endTime");
   return ReadDouble(node, "endTime", 0.0);
}

/**
 * @brief Read-out the Courant–Friedrichs–Lewy number.
 * @return CFL number CFL number. Default 0.6.
 */
double InputFileParser::ReadCFLnumber() const {
   tinyxml2::XMLElement const* node = inputfile_.FirstChildElement()->FirstChildElement("calculation")->FirstChildElement("CFLNumber");
   return ReadDouble(node, "CFLNumber", 0.6);
}

/**
 * @brief Reads-out the reference length to (non-)dimensionalize the (inputs) outputs.
 * @return Reference length. Default: 1.0.
 */
double InputFileParser::ReadReferenceLength() const {
   tinyxml2::XMLElement const* node = inputfile_.FirstChildElement()->FirstChildElement("calculation")->FirstChildElement("referenceParameter")->FirstChildElement("lengthReference");
   return ReadDouble(node, "lengthReference", 1.0);
}

/**
 * @brief Reads out the reference velocity to (non-)dimensionalize the (inputs) outputs.
 * @return Reference velocity. Default: 1.0.
 */
double InputFileParser::ReadReferenceVelocity() const {
   tinyxml2::XMLElement const* node = inputfile_.FirstChildElement()->FirstChildElement("calculation")->FirstChildElement("referenceParameter")->FirstChildElement("velocityReference");
   return ReadDouble(node, "velocityReference", 1.0);
}

/**
 * @brief Reads-out the reference density to (non-)dimensionalize the (inputs) outputs.
 * @return Reference density. Default: 1.0.
 */
double InputFileParser::ReadReferenceDensity() const {
   tinyxml2::XMLElement const* node = inputfile_.FirstChildElement()->FirstChildElement("calculation")->FirstChildElement("referenceParameter")->FirstChildElement("densityReference");
   return ReadDouble(node, "densityReference", 1.0);
}

/**
 * @brief Reads-out the reference density to (non-)dimensionalize the (inputs) outputs.
 * @return Reference density. Default: 1.0.
 */
double InputFileParser::ReadReferenceTemperature() const {
   tinyxml2::XMLElement const* node = inputfile_.FirstChildElement()->FirstChildElement("calculation")->FirstChildElement("referenceParameter")->FirstChildElement("temperatureReference");
   return ReadDouble(node, "temperatureReference", 1.0);
}


/**
 * @brief Reads out the cell-volume threshold below which cut-cells are mixed.
 * @return Cell-volume threshold for cut-cell mixing. Default: 0.6.
 */
double InputFileParser::ReadMixingThreshold() const {
   tinyxml2::XMLElement const* node = inputfile_.FirstChildElement()->FirstChildElement("calculation")->FirstChildElement("levelSet")->FirstChildElement("mixingThreshold");
   return ReadDouble(node, "mixingThreshold", 0.6);
}

/**
 * @brief Reads out if output should be created at all during the simulation.
 * @return True if output is desired. False otherwise. Default: 1.
 */
bool InputFileParser::ReadEnableOutput() const {
   tinyxml2::XMLElement const* node = inputfile_.FirstChildElement()->FirstChildElement("output")->FirstChildElement("enableOutput");
   return static_cast<bool>(ReadInt(node, "enableOutput", 1));
}

/**
 * @brief Reads out the format in which the output should be written.
 * @return The output format in machine-readable format (type identifier).
 */
OutputType InputFileParser::ReadOutputFormat() const {
   tinyxml2::XMLElement const *node = inputfile_.FirstChildElement()->FirstChildElement("output")->FirstChildElement("outputFileType");
   std::string type = ReadString(node);
   // NH 2017-02-08: Switch Statements do not work with String
   if(type == "XDMF") {return OutputType::Xdmf;}
   else { throw std::logic_error( "Unknown output type" ); }
}

/**
 * @brief Reads out the way output times are determined. Inputs should be 'Interval' or 'Timestamps'.
 * @return Output times type.
 */
std::string InputFileParser::ReadOutputTimesType() const {
   tinyxml2::XMLElement const* node = inputfile_.FirstChildElement()->FirstChildElement("output")->FirstChildElement("outputTimesType");
   return ReadString(node);
}

/**
 * @brief Reads out the output period. The period defines the time interval in between two outputs.
 * @return Interval after which the next output is to be written. Default: 1.0.
 */
double InputFileParser::ReadOutputPeriod() const {
   tinyxml2::XMLElement const *node = inputfile_.FirstChildElement()->FirstChildElement("output")->FirstChildElement("outputPeriod");
   return ReadDouble(node, "outputPeriod", 1.0);
}

/**
 * @brief Reads out the timestamps to be used for outputs.
 * @return Timestamps.
 */
std::vector<double> InputFileParser::ReadOutputTimestamps() const {
   return ReadTimestamps(inputfile_.FirstChildElement()->FirstChildElement("output")->FirstChildElement("timestamps"));
}

/**
 * @brief Reads out the detail to which the timestep is printed, i. e. number of decimals shown.
 * @return Number with the corresponding amount of decimals. Default: 1, i.e. output time = simulation time.
 */
double InputFileParser::ReadTimeNamingFactor() const {
   tinyxml2::XMLElement const* node = inputfile_.FirstChildElement()->FirstChildElement("output")->FirstChildElement("timeNamingFactor");
   return ReadDouble(node, "timeNamingFactor", 1.0);
}
