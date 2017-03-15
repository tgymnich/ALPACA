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
#include "heat_fluxes.h"
#include "stencils/differentiation_utilities.h"
#include "stencils/spatial_derivative_stencils/fourth_order_cell_face.h"

/**
 * @brief The default constructor for the class.
 * @param material_manager The material manager.
 */
HeatFluxes::HeatFluxes(MaterialManager const& material_manager) :
    material_manager_(material_manager)
{
   // Empty constructor besides initializer list
}

/**
 * @brief Computes the fluxes and adds them to buffers for fluxes in x-, y- and z-direction.
 * @param mat_block The material block pair of the fluid for which the fluxes are calculated.
 * @param heat_fluxes_x The heat fluxes in x-direction (indirect return parameter).
 * @param heat_fluxes_y The heat fluxes in y-direction (indirect return parameter).
 * @param heat_fluxes_z The heat fluxes in z-direction (indirect return parameter).
 * @param cell_size The cell size of the block.
 */
void HeatFluxes::ComputeFluxes(const std::pair<const MaterialName, Block> &mat_block,
   double (&heat_fluxes_x)[FF::ANOE()][CC::ICX() + 1][CC::ICY() + 1][CC::ICZ() + 1],
   double (&heat_fluxes_y)[FF::ANOE()][CC::ICX() + 1][CC::ICY() + 1][CC::ICZ() + 1],
   double (&heat_fluxes_z)[FF::ANOE()][CC::ICX() + 1][CC::ICY() + 1][CC::ICZ() + 1],
   double const cell_size) const {

   using TemperatureGradientDerivativeStencilCellFaceConcretization = DerivativeStencilSetup::Concretize<temperature_gradient_derivative_stencil_cell_face>::type;

   double const thermal_conductivity = material_manager_.GetThermalConductivity(mat_block.first);

   double const (&temperature)[CC::TCX()][CC::TCY()][CC::TCZ()] = mat_block.second.GetPrimeStateBuffer(PrimeState::Temperature);

   // Offset for flux-arrays
   constexpr int offset_x = CC::FICX() - 1;
   constexpr int offset_y = CC::DIM() != Dimension::One   ? CC::FICY() - 1 : -1;
   constexpr int offset_z = CC::DIM() == Dimension::Three ? CC::FICZ() - 1 : -1;

   // Compute changes due to heat transfer - compute second order derivative of temperature with fourth-order central-difference stencil
   for(unsigned int i = CC::FICX()-1; i <= CC::LICX(); ++i){
      for(unsigned int j = CC::FICY(); j <= CC::LICY(); ++j){
         for(unsigned int k = CC::FICZ(); k <= CC::LICZ(); ++k){
            heat_fluxes_x[ETI(Equation::Energy)][i-offset_x][j-offset_y][k-offset_z] += -thermal_conductivity * DifferentiationUtilities::ComputeDerivative<TemperatureGradientDerivativeStencilCellFaceConcretization, Direction::X, double>(temperature, i, j, k, cell_size);
         }
      }
   }

   if constexpr(CC::DIM() != Dimension::One) {
      for(unsigned int i = CC::FICX(); i <= CC::LICX(); ++i){
         for(unsigned int j = CC::FICY()-1; j <= CC::LICY(); ++j){
            for(unsigned int k = CC::FICZ(); k <= CC::LICZ(); ++k){
               heat_fluxes_y[ETI(Equation::Energy)][i-offset_x][j-offset_y][k-offset_z] += -thermal_conductivity * DifferentiationUtilities::ComputeDerivative<TemperatureGradientDerivativeStencilCellFaceConcretization, Direction::Y, double>(temperature, i, j, k, cell_size);
            }
         }
      }
   }

   if constexpr(CC::DIM() == Dimension::Three) {
      for(unsigned int i = CC::FICX(); i <= CC::LICX(); ++i){
         for(unsigned int j = CC::FICY(); j <= CC::LICY(); ++j){
            for(unsigned int k = CC::FICZ()-1; k <= CC::LICZ(); ++k){
               heat_fluxes_z[ETI(Equation::Energy)][i-offset_x][j-offset_y][k-offset_z] += -thermal_conductivity * DifferentiationUtilities::ComputeDerivative<TemperatureGradientDerivativeStencilCellFaceConcretization, Direction::Z, double>(temperature, i, j, k, cell_size);
            }
         }
      }
   }
}
