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
#include "axisymmetric_viscous_volume_forces.h"
#include "levelset/multi_phase_manager/material_sign_capsule.h"
#include "index_transformations.h"
#include "stencils/differentiation_utilities.h"

using ViscousFluxesDerivativeStencil = DerivativeStencilSetup::Concretize<viscous_fluxes_derivative_stencil>::type;
using ViscousFluxesReconstructionStencil = ReconstructionStencilSetup::Concretize<viscous_fluxes_reconstruction_stencil>::type;

/**
 * @brief The constructor for the AxisymmetricViscousVolumeForces class.
 * @param material_manager The material manager which contains information about the viscosities.
 */
AxisymmetricViscousVolumeForces::AxisymmetricViscousVolumeForces( MaterialManager const& material_manager ) :
   material_manager_( material_manager )
{
   // Empty constructor besides initializer list
}

/**
 * @brief Computes axisymmetric increments for cell averages due to viscosity.
 * @param mat_block A pair containing a block and its material.
 * @param axisymmetric_viscous_volume_forces Reference to array of volume forces increments to be filled here (indirect return parameter).
 * @param cell_size The cell size.
 */
void AxisymmetricViscousVolumeForces::ComputeForces( std::pair<MaterialName const, Block> const& mat_block,
                                                     double (&axisymmetric_viscous_volume_forces)[FF::ANOE()][CC::ICX()][CC::ICY()][CC::ICZ()],
                                                     double const cell_size, double const x_block_coordinate ) const {

   double const (&u)[CC::TCX()][CC::TCY()][CC::TCZ()] = mat_block.second.GetPrimeStateBuffer( PrimeState::VelocityX );
   double const (&v)[CC::TCX()][CC::TCY()][CC::TCZ()] = mat_block.second.GetPrimeStateBuffer( PrimeState::VelocityY );
   double const (&w)[CC::TCX()][CC::TCY()][CC::TCZ()] = u;

   double velocity_gradient[CC::TCX()][CC::TCY()][1][dim_][dim_];
   for( unsigned int i = 0; i < CC::TCX(); ++i ) {
      for( unsigned int j = 0; j < CC::TCY(); ++j ) {
         for( unsigned int r = 0; r < dim_; ++r ) {
            for( unsigned int s = 0; s < dim_; ++s ) {
               velocity_gradient[i][j][0][r][s] = 0.0;
            }
         }
      } //j
   } //i

   ComputeVelocityGradient( u, v, w, cell_size, velocity_gradient );

   std::vector<double> const viscosity = material_manager_.GetViscosity( mat_block.first );
   double const mu_1 = viscosity[0];
   double const mu_2 = viscosity[1] - 2.0 * viscosity[0] / 3.0;

   for(unsigned int i = 0; i < CC::ICX(); ++i) {
      double const one_radius = 1.0 / ( x_block_coordinate + ( static_cast<double>( i ) + 0.5 ) * cell_size );
      for(unsigned int j = 0; j < CC::ICY(); ++j) {

         std::array<unsigned int, 3> const indices = { BIT::I2TX( i ), BIT::I2TY( j ), 0 };

         // Add up to volume forces
         axisymmetric_viscous_volume_forces[ETI( Equation::Mass  )  ][i][j][0] += 0.0;
         axisymmetric_viscous_volume_forces[ETI( Equation::Energy ) ][i][j][0] +=
            ( 2.0 * mu_1 + mu_2 ) * u[indices[0]][indices[1]][indices[2]] * u[indices[0]][indices[1]][indices[2]] * one_radius * one_radius +
              2.0 * mu_2          * u[indices[0]][indices[1]][indices[2]] * one_radius * ( velocity_gradient[indices[0]][indices[1]][indices[2]][0][0] + velocity_gradient[indices[0]][indices[1]][indices[2]][1][1] );

         axisymmetric_viscous_volume_forces[ETI( Equation::MomentumX )][i][j][0] +=
            ( 2.0 * mu_1 + mu_2 ) * ( velocity_gradient[indices[0]][indices[1]][indices[2]][0][0] - u[indices[0]][indices[1]][indices[2]] * one_radius ) * one_radius;
         axisymmetric_viscous_volume_forces[ETI( Equation::MomentumY )][i][j][0] +=
            ( mu_1 * velocity_gradient[indices[0]][indices[1]][indices[2]][1][0] + ( mu_1 + mu_2 ) * velocity_gradient[indices[0]][indices[1]][indices[2]][0][1] ) * one_radius;
      }
   }
}

/**
 * @brief Computes the gradient of the velocity field at cell centers.
 * @param u The velocity field in x direction.
 * @param v The velocity field in y direction.
 * @param w The velocity field in z direction.
 * @param cell_size The cell size.
 * @param velocity_gradient The velocity gradient field as indirect return parameter.
 */
void AxisymmetricViscousVolumeForces::ComputeVelocityGradient( double const (&u)[CC::TCX()][CC::TCY()][CC::TCZ()],
                                                               double const (&v)[CC::TCX()][CC::TCY()][CC::TCZ()],
                                                               double const (&w)[CC::TCX()][CC::TCY()][CC::TCZ()], double const cell_size,
                                                               double (&velocity_gradient)[CC::TCX()][CC::TCY()][1][dim_][dim_] ) const {
   /**
    * @brief Offsets in order to also calculate first derivatives in halo cells. This is necessary for the reconstruction to cell faces.
    */
   constexpr unsigned int offset_x = ViscousFluxesDerivativeStencil::DownstreamStencilSize();
   constexpr unsigned int offset_y = ViscousFluxesDerivativeStencil::DownstreamStencilSize();

   for( unsigned int i = 0 + offset_x; i < CC::TCX() - offset_x; ++i ) {
      for( unsigned int j = 0 + offset_y; j < CC::TCY() - offset_y; ++j ) {
         std::array< std::array<double, 3>, 3> const single_gradient = DifferentiationUtilities::VectorGradient<ViscousFluxesDerivativeStencil, double, Dimension::Two>( u, v, w, i, j, 0, cell_size );
         for( unsigned int r = 0; r < dim_; ++r ) {
            for( unsigned int s = 0; s < dim_; ++s ) {
               velocity_gradient[i][j][0][r][s] = single_gradient[r][s];
            }
         }
      }
   }
}
