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
#include "initial_condition.h"
#include "simulation_setup.h"
#include "user_expression.h"

#include <stdexcept>

const std::string InitialCondition::variable_name_phi_ = "phi";

/**
 * @brief Default constructor requiring a parse-able inputfile %Currently only XML file%
 * @param parser An user inputfile parser instance.
 */
InitialCondition::InitialCondition(const InputFileParser& parser) :
   phi_expression_(parser.ReadInitialConditionOfLevelSet()),
   number_of_fluids_(parser.ReadNumberOfFluids()),
   fluid_initialisation_(parser.ReadInitialConditionOfFluids()),
   fluids_(SimulationSetup::MapUserInputToMaterialType(parser.ReadParameterOfAllFluids()))
{
}

/**
 * @brief Compiles the given expression such that it is capable to give the desired fluid variables.
 * @param expression The expression in original text form.
 * @param x,y,z Reference to spatial variable.
 * @return The compiled expression ready to be evaluated.
 */
UserExpression InitialCondition::CreateInputExpression( std::string const expression, std::vector<std::string> const& variables_out, double &x, double &y, double &z ) const {
   std::vector<std::tuple<std::string,double&>> variables_in;
   variables_in.push_back( std::make_tuple(std::string("x"), std::ref(x)) );
   variables_in.push_back( std::make_tuple(std::string("y"), std::ref(y)) );
   variables_in.push_back( std::make_tuple(std::string("z"), std::ref(z)) );

   return UserExpression( expression, variables_in, variables_out );
}

/**
 * @brief Evaluates the user input for the initial density at the provided X-/Y- and Z-coordinates.
 * @param origin The origin coordinates of the node.
 * @param cell_size The size of cells in the node.
 * @param material Identifier to fill the cells with the state associated with its material.
 */
void InitialCondition::GetInitialPrimeStates( const std::array<double, 3> origin, const double cell_size, const MaterialName material,
   double (&initial_values)[FF::ANOP()][CC::ICX()][CC::ICY()][CC::ICZ()] ) const {

   std::string initial_condition_input;

   for( unsigned int ii = 0; ii < number_of_fluids_; ++ii ){
      if( fluids_[ii] == material ){
         initial_condition_input = fluid_initialisation_[ii];
      }
   }

   std::vector<std::string> variable_names;
   for( PrimeState const p : FF::ASOP() ) {
      // TP: Do not check for empty names here! We need them later to assign 0.0 to the respective prime states and exprtk can handle the empty names
      variable_names.emplace_back( FF::FieldInputName( p ) );
   }

   // TP here we need non-const variables as UserExpression stores references to them in order to reflect changes
   double running_x;
   double running_y;
   double running_z;

   UserExpression const input_expression = CreateInputExpression( initial_condition_input, variable_names, running_x, running_y, running_z );

   for(unsigned int i = 0; i < CC::ICX(); ++i) {
      running_x = origin[0] + (double(i) + 0.5) * cell_size;
      for(unsigned int j = 0; j < CC::ICY(); ++j) {
         running_y = CC::DIM() != Dimension::One ? origin[1] + (double(j) + 0.5) * cell_size : 0.0;
         for(unsigned int k = 0; k < CC::ICZ(); ++k) {
            running_z = CC::DIM() == Dimension::Three ? origin[2] + (double(k) + 0.5) * cell_size : 0.0;
            for( unsigned int p = 0; p < FF::ANOP(); ++p ) {
               if( variable_names[p].empty() ) {
                  initial_values[p][i][j][k] = 0.0;
               } else {
                  initial_values[p][i][j][k] = input_expression.GetValue( variable_names[p] );
               }
            }
         }
      }
   }
}

/**
 * @brief Evaluates the user input for the initial levelset for a node with the given coordinates and expansion.
 * @param origin The origin coordinates of the node.
 * @param cell_size The size of cells in the node.
 * @param initial_values Indirect return parameter for the determined levelset values.
 * @note Works on total cells
 */
void InitialCondition::GetInitialLevelset(const std::array<double, 3> origin, const double cell_size, double (&initial_values)[CC::TCX()][CC::TCY()][CC::TCZ()]) const{

   std::string number_phi = phi_expression_[0];

   // TP here we need non-const variables as UserExpression stores references to them in order to reflect changes
   double running_x;
   double running_y;
   double running_z;
   const double one_cell_size = 1.0 / cell_size;

   UserExpression const phi_expr = CreateInputExpression(number_phi, {variable_name_phi_}, running_x, running_y, running_z);

   for(unsigned int i = 0; i < CC::TCX(); ++i) {
      running_x = origin[0] + (double(i) - double(CC::FICX()) + 0.5) * cell_size;
      for(unsigned int j = 0; j < CC::TCY(); ++j) {
         running_y = CC::DIM() != Dimension::One ? origin[1] + (double(j) - double(CC::FICY()) + 0.5) * cell_size : 0.0;
         for(unsigned int k = 0; k < CC::TCZ(); ++k) {
            running_z = CC::DIM() == Dimension::Three ? origin[2] + (double(k) - double(CC::FICZ()) + 0.5) * cell_size : 0.0;
            initial_values[i][j][k] = phi_expr.GetValue(variable_name_phi_) * one_cell_size;
         }
      }
   }

}

/**
 * @brief Evaluates the user input for the initial materials for a node with the given coordinates and expansion.
 * @param origin The origin coordinates of the node.
 * @param smallest_cell_size The smallest possible size of cells in the simulation, i.e. on the maximum level.
 * @param level_factor The refinement factor between the considered level and the finest level, e.g. 4 if level 2 is considered and the maximum is level 4.
 * @return A boolean value for each material indicating whether the material is present in the considered region.
 * @note Works on total cells.
 */
std::vector<bool> InitialCondition::GetInitialMaterials(const std::array<double, 3> origin, const double smallest_cell_size, const unsigned int level_factor) const {

   std::vector<bool> material_is_contained(fluids_.size());

   if(fluids_.size() == 1) {
      material_is_contained[0] = true;
      return material_is_contained;
   } else {

      // TP here we need non-const variables as UserExpression stores references to them in order to reflect changes
      double running_x;
      double running_y;
      double running_z;

      if(fluids_.size() == 2) {
         // if we only have two materials (positive, negative) we only consider the sign of the single level-set function
         // first material is the positive, second one the negative

         std::string number_phi = phi_expression_[0];
         UserExpression const phi_expr = CreateInputExpression(number_phi, {variable_name_phi_}, running_x, running_y, running_z);

         material_is_contained[0] = false; // positive material
         material_is_contained[1] = false; // negative material

         const unsigned int level_factor_y = CC::DIM() != Dimension::One   ? level_factor : 1;
         const unsigned int level_factor_z = CC::DIM() == Dimension::Three ? level_factor : 1;

         // in the following loop we run over all cells on the finest level that would cover the same region as the node of interest including its halo region
         double phi;
         for(unsigned int i = 0; i < level_factor*CC::TCX(); ++i) {
            running_x = origin[0] + (double(i) - double(level_factor*CC::FICX()) + 0.5) * smallest_cell_size;
            for(unsigned int j = 0; j < level_factor_y*CC::TCY(); ++j) {
               running_y = CC::DIM() != Dimension::One ? origin[1] + (double(j) - double(level_factor_y*CC::FICY()) + 0.5) * smallest_cell_size : 0.0;
               for(unsigned int k = 0; k < level_factor_z*CC::TCZ(); ++k) {
                  running_z = CC::DIM() == Dimension::Three ? origin[2] + (double(k) - double(level_factor_z*CC::FICZ()) + 0.5) * smallest_cell_size : 0.0;
                  phi = phi_expr.GetValue(variable_name_phi_);
                  if(phi < 0.0) {
                     material_is_contained[1] = true;
                  } else {
                     material_is_contained[0] = true;
                  }
                  // if both materials are contained, we can already return
                  if(material_is_contained[0] && material_is_contained[1]) return material_is_contained;
               }
            }
         }
      } else { // real multi-fluid
         // for multi fluids we have one level-set function for each fluid, thus there is a one-to-one mapping
         throw std::logic_error("InitialCondition::GetInitialMaterials not implemented yet for multi-fluid");
      }
   }

   return material_is_contained;
}
