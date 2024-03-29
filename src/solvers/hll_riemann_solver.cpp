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
#include "hll_riemann_solver.h"

#include "mathematical_functions.h"
#include "solvers/hll_signal_speed_calculator.h"
#include "stencils/differentiation_utilities.h"

/**
 * @brief Standard constructor using an already existing MaterialManager and EigenDecomposition object. See base class.
 * @param material_manager .
 * @param eigendecomposition_calculator .
 */
HllRiemannSolver::HllRiemannSolver( MaterialManager const& material_manager, EigenDecomposition const& eigendecomposition_calculator ) :
   RiemannSolver( material_manager, eigendecomposition_calculator )
{
   /* Empty besides initializer list*/
}

/**
 * @brief Solving the right hand side of the underlying system of equations. Using local characteristic decomposition
 * in combination with spatial reconstruction of cell averaged values (finite volume approach) and flux determination by HLL procedure.
 * See base class.
 */
void HllRiemannSolver::UpdateImplementation( std::pair<MaterialName const, Block> const& mat_block, double const cell_size,
   double (&fluxes_x)[FF::ANOE()][CC::ICX()+1][CC::ICY()+1][CC::ICZ()+1],
   double (&fluxes_y)[FF::ANOE()][CC::ICX()+1][CC::ICY()+1][CC::ICZ()+1],
   double (&fluxes_z)[FF::ANOE()][CC::ICX()+1][CC::ICY()+1][CC::ICZ()+1] ) const {

   double  roe_eigenvectors_left[CC::ICX()+1][CC::ICY()+1][CC::ICZ()+1][FF::ANOE()][FF::ANOE()];
   double roe_eigenvectors_right[CC::ICX()+1][CC::ICY()+1][CC::ICZ()+1][FF::ANOE()][FF::ANOE()];
   double        roe_eigenvalues[CC::ICX()+1][CC::ICY()+1][CC::ICZ()+1][FF::ANOE()];

   // Resetting the eigenvectors, eigenvalues is necessary for two-phase simulations.
   for( unsigned int i = 0; i < CC::ICX() + 1; ++i ) {
      for( unsigned int j = 0; j < CC::ICY() + 1; ++j ) {
         for( unsigned int k = 0; k < CC::ICZ() + 1; ++k ) {
            for( unsigned int e = 0; e < FF::ANOE(); ++e ) {
               for( unsigned int f = 0; f < FF::ANOE(); ++f ) {
                  roe_eigenvectors_left[i][j][k][e][f] = 0.0;
                  roe_eigenvectors_right[i][j][k][e][f] = 0.0;
               }
               roe_eigenvalues[i][j][k][e] = 0.0;
            }
         }
      }
   }

   eigendecomposition_calculator_.ComputeRoeEigendecomposition<Direction::X>( mat_block, roe_eigenvectors_left, roe_eigenvectors_right, roe_eigenvalues );
   ComputeFluxes<Direction::X>( mat_block, fluxes_x, roe_eigenvectors_left, roe_eigenvectors_right, cell_size );

   if constexpr( CC::DIM() != Dimension::One ) {
      eigendecomposition_calculator_.ComputeRoeEigendecomposition<Direction::Y>( mat_block, roe_eigenvectors_left, roe_eigenvectors_right, roe_eigenvalues );
      ComputeFluxes<Direction::Y>( mat_block, fluxes_y, roe_eigenvectors_left, roe_eigenvectors_right, cell_size );
   }

   if constexpr( CC::DIM() == Dimension::Three ) {
      eigendecomposition_calculator_.ComputeRoeEigendecomposition<Direction::Z>( mat_block, roe_eigenvectors_left, roe_eigenvectors_right, roe_eigenvalues );
      ComputeFluxes<Direction::Z>( mat_block, fluxes_z, roe_eigenvectors_left, roe_eigenvectors_right, cell_size );
   }
}

/**
 * @brief Computes the cell face fluxes with the set stencil using local characteristic decomposition & HLL procedure
 * @param mat_block The block and material information of the phase under consideration.
 * @param fluxes Reference to an array which is filled with the computed fluxes (indirect return parameter).
 * @param roe_eigenvectors_left .
 * @param roe_eigenvectors_right .
 * @param cell_size .
 * @tparam DIR Indicates which spatial direction is to be computed.
 * @note Hotpath function.
 */
template<Direction DIR>
void HllRiemannSolver::ComputeFluxes( std::pair<MaterialName const, Block> const& mat_block,
   double (&fluxes)[FF::ANOE()][CC::ICX()+1][CC::ICY()+1][CC::ICZ()+1],
   double const  (&Roe_eigenvectors_left)[CC::ICX()+1][CC::ICY()+1][CC::ICZ()+1][FF::ANOE()][FF::ANOE()],
   double const (&Roe_eigenvectors_right)[CC::ICX()+1][CC::ICY()+1][CC::ICZ()+1][FF::ANOE()][FF::ANOE()],
   double const cell_size ) const {

   using ReconstructionStencilConcretization = ReconstructionStencilSetup::Concretize<reconstruction_stencil>::type;

   // declaration of applied vectors
   std::vector<double> state_face_left( FF::ANOE() );             //variable vector containing interpolated states of left patch of cell face i/j/k+1/2
   std::vector<double> state_face_right( FF::ANOE() );            //variable vector containing interpolated states of right patch of cell face i/j/k+1/2
   std::vector<double> u_characteristic( ReconstructionStencilConcretization::StencilSize() );    //temp storage for characteristic decomposition
   std::vector<double> characteristic_average_plus( FF::ANOE() ); // state_face_left  in characteristic space
   std::vector<double> characteristic_average_minus( FF::ANOE() );// state_face_right in characteristic space
   std::vector<double> flux_left( FF::ANOE() );                   // F(state_face_left)
   std::vector<double> flux_right( FF::ANOE() );                  // F(state_face_right)

   unsigned int constexpr principal_momentum = momentum_order_[DTI(DIR)][0];

   constexpr unsigned int x_start = DIR == Direction::X ? CC::FICX()-1 : CC::FICX();
   constexpr unsigned int y_start = DIR == Direction::Y ? CC::FICY()-1 : CC::FICY();
   constexpr unsigned int z_start = DIR == Direction::Z ? CC::FICZ()-1 : CC::FICZ();

   constexpr unsigned int x_reconstruction_offset  = DIR == Direction::X ? 1 : 0;
   constexpr unsigned int y_reconstruction_offset  = DIR == Direction::Y ? 1 : 0;
   constexpr unsigned int z_reconstruction_offset  = DIR == Direction::Z ? 1 : 0;

   constexpr unsigned int x_end = CC::LICX();
   constexpr unsigned int y_end = CC::LICY();
   constexpr unsigned int z_end = CC::LICZ();

   constexpr int total_to_internal_offset_x = CC::FICX()-1 ;
   constexpr int total_to_internal_offset_y = CC::DIM() != Dimension::One   ? static_cast<int>( CC::FICY() ) - 1 : -1;
   constexpr int total_to_internal_offset_z = CC::DIM() == Dimension::Three ? static_cast<int>( CC::FICZ() ) - 1 : -1;

   // Access the pair's elements directly.
   auto const& [material, block] = mat_block;

   // To check for invalid cells due to ghost fluid method
   double const B = material_manager_.GetB( material );

   // Required for Toro signal speeds
   double const gamma = material_manager_.GetGamma( material );

   for( unsigned int i = x_start; i <= x_end; ++i ) {
      for( unsigned int j = y_start; j <= y_end; ++j ) {
         for( unsigned int k = z_start; k <= z_end; ++k ) {
            // Shifted indices to match block index system and roe-ev index system
            int const i_index = i-total_to_internal_offset_x;
            int const j_index = j-total_to_internal_offset_y;
            int const k_index = k-total_to_internal_offset_z;

            // Reconstruct fluxes at face i+1/2 by characteristic decomposition using a selected reconstruction procedure
            for( unsigned int n = 0; n < FF::ANOE(); ++n ) { // n is index of characteristic field (eigenvalue, eigenvector)
               // Characteristic decomposition
               for( unsigned int m = 0; m < ReconstructionStencilConcretization::StencilSize(); ++m ) {
                  u_characteristic[m] = 0.0;
                  // Compute characteristics for U
                  for( const unsigned int l : conservative_equation_summation_sequence_[DTI(DIR)] ) { // l is index of conservative equation, iterated in symmetry-preserving sequence
                     u_characteristic[m] += Roe_eigenvectors_left[i_index][j_index][k_index][n][l] *
                        block.GetAverageBuffer( FF::ASOE()[l] )[i + x_reconstruction_offset * ( m - ReconstructionStencilConcretization::DownstreamStencilSize())]
                                                               [j + y_reconstruction_offset * ( m - ReconstructionStencilConcretization::DownstreamStencilSize())]
                                                               [k + z_reconstruction_offset * ( m - ReconstructionStencilConcretization::DownstreamStencilSize())];
                  } // L-Loop
               } // M-Loop

               // Apply reconstruction scheme to characteristic values
               characteristic_average_minus[n]  = DifferentiationUtilities::ApplyStencil<ReconstructionStencilConcretization, StencilProperty::UpwindLeft,  double>( u_characteristic, cell_size );
               characteristic_average_plus[n]   = DifferentiationUtilities::ApplyStencil<ReconstructionStencilConcretization, StencilProperty::UpwindRight, double>( u_characteristic, cell_size );
            } // N-Loop

            // Back-transformation into physical space
            for( unsigned int l = 0; l < FF::ANOE(); ++l ) {
               state_face_left[l]  = 0.0;
               state_face_right[l] = 0.0;
               for( unsigned int const n : characteristic_field_summation_sequence_[DTI( DIR )] ) {
                  state_face_left[l]  += characteristic_average_minus[n] * Roe_eigenvectors_right[i_index][j_index][k_index][l][n];
                  state_face_right[l] += characteristic_average_plus[n]  * Roe_eigenvectors_right[i_index][j_index][k_index][l][n];
               } // N-Loop
               // Non-linear contributions have to be added together to maintain full symmetry
               state_face_left [l] += ( characteristic_average_minus[0]            * Roe_eigenvectors_right[i_index][j_index][k_index][l][0] +
                                        characteristic_average_minus[FF::ANOE()-1] * Roe_eigenvectors_right[i_index][j_index][k_index][l][FF::ANOE()-1] );
               state_face_right[l] += ( characteristic_average_plus[0]             * Roe_eigenvectors_right[i_index][j_index][k_index][l][0] +
                                        characteristic_average_plus[FF::ANOE()-1]  * Roe_eigenvectors_right[i_index][j_index][k_index][l][FF::ANOE()-1] );
            } // L-Loop

            // Check for invalid cells due to ghost fluid method
            if( state_face_left[ETI( Equation::Mass )] <= std::numeric_limits<double>::epsilon() || state_face_right[ETI( Equation::Mass )] <= std::numeric_limits<double>::epsilon() ) continue;

            // Compute pressure, velocity and speed of sound for both cells for reconstructed values
            double const pressure_left  = material_manager_.GetPressure( material,
                                                                         state_face_left[ETI( Equation::Mass )],
                                                                         state_face_left[ETI( Equation::MomentumX )],
                                         CC::DIM()  > Dimension::One   ? state_face_left[ETI( Equation::MomentumY )] : 0.0,
                                         CC::DIM() == Dimension::Three ? state_face_left[ETI( Equation::MomentumZ )] : 0.0,
                                                                         state_face_left[ETI( Equation::Energy )] );
            double const pressure_right = material_manager_.GetPressure( material,
                                                                         state_face_right[ETI( Equation::Mass )],
                                                                         state_face_right[ETI( Equation::MomentumX )],
                                         CC::DIM()  > Dimension::One   ? state_face_right[ETI( Equation::MomentumY )] : 0.0,
                                         CC::DIM() == Dimension::Three ? state_face_right[ETI( Equation::MomentumZ )] : 0.0,
                                                                         state_face_right[ETI( Equation::Energy )] );

            // Check for invalid cells due to ghost fluid method
            if( pressure_left <= -B || pressure_right <= -B ) continue;

            double const velocity_left  = state_face_left [principal_momentum] / state_face_left [ETI( Equation::Mass )];
            double const velocity_right = state_face_right[principal_momentum] / state_face_right[ETI( Equation::Mass )];
            double const speed_of_sound_left  = material_manager_.GetSpeedOfSound( material, state_face_left[ ETI( Equation::Mass )], pressure_left );
            double const speed_of_sound_right = material_manager_.GetSpeedOfSound( material, state_face_right[ETI( Equation::Mass )], pressure_right );

            // Calculation of signal speeds
            auto const [wave_speed_left_simple, wave_speed_right_simple] = CalculateSignalSpeed (state_face_left [ETI( Equation::Mass )], state_face_right [ETI( Equation::Mass )],
                                                                                                 velocity_left, velocity_right,
                                                                                                 pressure_left, pressure_right,
                                                                                                 speed_of_sound_left, speed_of_sound_right,
                                                                                                 gamma);

            double const wave_speed_left  = std::min( wave_speed_left_simple, 0.0 );
            double const wave_speed_right = std::max( wave_speed_right_simple, 0.0 );

            // Calculation of left and right flux
            flux_left[ETI( Equation::Mass )]      = state_face_left[principal_momentum];
            flux_left[principal_momentum]         = ( ( state_face_left[principal_momentum] * state_face_left[principal_momentum] ) / state_face_left[ETI( Equation::Mass )] ) + pressure_left;
            flux_left[ETI( Equation::Energy )]    = velocity_left * (state_face_left[ETI( Equation::Energy )] + pressure_left);
            flux_right[ETI( Equation::Mass )]     = state_face_right[principal_momentum];
            flux_right[principal_momentum]        = ( ( state_face_right[principal_momentum] * state_face_right[principal_momentum] ) / state_face_right[ETI( Equation::Mass )] ) + pressure_right;
            flux_right[ETI( Equation::Energy )]   = velocity_right * (state_face_right[ETI( Equation::Energy )] + pressure_right);

            // Momenta besides principal momentum (start iterating at 1)
            for( unsigned int d = 1; d < DTI( CC::DIM() ); ++d ) {
               flux_left[ momentum_order_[DTI( DIR )][d]] = velocity_left  *  state_face_left[momentum_order_[DTI( DIR )][d]];
               flux_right[momentum_order_[DTI( DIR )][d]] = velocity_right * state_face_right[momentum_order_[DTI( DIR )][d]];
            }

            for( unsigned int n = 0; n < FF::ANOE(); ++n ) {
               fluxes[n][i_index][j_index][k_index] += ( ( ( wave_speed_right * flux_left[n] ) - ( wave_speed_left * flux_right[n] ) ) +
                                                         ( ( wave_speed_right * wave_speed_left ) * ( state_face_right[n] - state_face_left[n] ) ) )
                                                       / ( wave_speed_right - wave_speed_left );
            }
         } // k
      } // j
   } // i
}
