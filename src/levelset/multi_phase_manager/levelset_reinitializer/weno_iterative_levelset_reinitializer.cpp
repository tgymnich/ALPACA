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
#include "weno_iterative_levelset_reinitializer.h"
#include "user_specifications/two_phase_constants.h"
#include "stencils/differentiation_utilities.h"

/**
 * @brief Computes the advection velocity used for the iterative reinitialization. The advection velocity is a smoothened signum of the levelset.
 * Smoothening is necessary for small absolute level-set values to not reinitialize there.
 * @param levelset The levelset-value.
 * @param godunov_hamiltonian The godunov hamiltonian, i.e. the absolute value of the level-set gradient.
 * @return The advection velocity.
 */
double ComputeAdvectionVelocity(double const levelset, double const godunov_hamiltonian) {
   constexpr double epsilon = 2.0e-5;
   return levelset / std::sqrt(levelset * levelset + (godunov_hamiltonian * godunov_hamiltonian) * (epsilon * epsilon));
}

/**
 * @brief The default constructor for a WenoIterativeLevelsetReinitializer. Calls the default constructor of the base class.
 * @param halo_manager See base class.
 */
WenoIterativeLevelsetReinitializer::WenoIterativeLevelsetReinitializer( HaloManager& halo_manager ) :
   IterativeLevelsetReinitializerBase( halo_manager )
{
   // Empty Constructor, besides call of base class constructor.
}

/**
 * @brief Reinitializes the levelset field of a node as described in \cite Min2010 .
 * @param node The levelset block which has to be reinitialized.
 * @return The residuum for the current node.
 */
double WenoIterativeLevelsetReinitializer::ReinitializeSingleNodeImplementation(Node& node) const {

   using ReconstructionStencilConcretization = ReconstructionStencilSetup::Concretize<reconstruction_stencil>::type;

   std::int8_t const  (&interface_tags)[CC::TCX()][CC::TCY()][CC::TCZ()] = node.GetInterfaceTags();
   LevelsetBlock& levelset_block = node.GetLevelsetBlock();
   double (&levelset_orig)[CC::TCX()][CC::TCY()][CC::TCZ()] = levelset_block.GetPhiReinitialized();
   double const (&levelset_0_orig)[CC::TCX()][CC::TCY()][CC::TCZ()] = levelset_block.GetPhiRightHandSide();

   double reinitialization_rhs[CC::TCX()][CC::TCY()][CC::TCZ()];

   double residuum = 0.0;

   std::vector<double> interpolation_array;

   double derivatives[DTI(CC::DIM())][2];
   for(unsigned int d = 0; d < DTI(CC::DIM()); ++d) {
      for(unsigned int e = 0; e < 2; ++e) {
         derivatives[d][e] = 0.0;
      }
   }

   for(unsigned int i = 0; i < CC::TCX(); ++i) {
      for(unsigned int j = 0; j < CC::TCY(); ++j) {
         for(unsigned int k = 0; k < CC::TCZ(); ++k) {
            reinitialization_rhs[i][j][k] = 0.0;
         } // k
      } // j
   } // i

   for(unsigned int i = CC::FICX(); i <= CC::LICX(); ++i) {
      for(unsigned int j = CC::FICY(); j <= CC::LICY(); ++j) {
         for(unsigned int k = CC::FICZ(); k <= CC::LICZ(); ++k) {
            if(std::abs(interface_tags[i][j][k]) > ITTI(IT::NewCutCell) || reinitialize_cut_cells_) {

               // We normalize the level-set field on the cell size, thus 1.0 is necessary as cell size for stencil evaluation.
               interpolation_array.clear();
               for(unsigned int ii = i - (ReconstructionStencilConcretization::DownstreamStencilSize() + 1); ii <= i + ReconstructionStencilConcretization::DownstreamStencilSize(); ii++) {
                  interpolation_array.push_back(levelset_orig[ii+1][j][k] - levelset_orig[ii][j][k]);
               }

               derivatives[0][0] = DifferentiationUtilities::ApplyStencil<ReconstructionStencilConcretization, StencilProperty::UpwindLeft, double>(interpolation_array, 1.0);
               derivatives[0][1] = DifferentiationUtilities::ApplyStencil<ReconstructionStencilConcretization, StencilProperty::UpwindRight, double>(interpolation_array, 1.0);


               if constexpr(CC::DIM() != Dimension::One) {

                  interpolation_array.clear();
                  for(unsigned int jj = j - (ReconstructionStencilConcretization::DownstreamStencilSize() + 1); jj <= j + ReconstructionStencilConcretization::DownstreamStencilSize(); jj++) {
                     interpolation_array.push_back(levelset_orig[i][jj+1][k] - levelset_orig[i][jj][k]);
                  }

                  derivatives[1][0] = DifferentiationUtilities::ApplyStencil<ReconstructionStencilConcretization, StencilProperty::UpwindLeft, double>(interpolation_array, 1.0);
                  derivatives[1][1] = DifferentiationUtilities::ApplyStencil<ReconstructionStencilConcretization, StencilProperty::UpwindRight, double>(interpolation_array, 1.0);
               }


               if constexpr(CC::DIM() == Dimension::Three) {

                  interpolation_array.clear();
                  for(unsigned int kk = k - (ReconstructionStencilConcretization::DownstreamStencilSize() + 1); kk <= k + ReconstructionStencilConcretization::DownstreamStencilSize(); kk++) {
                     interpolation_array.push_back(levelset_orig[i][j][kk+1] - levelset_orig[i][j][kk]);
                  }

                  derivatives[2][0] = DifferentiationUtilities::ApplyStencil<ReconstructionStencilConcretization, StencilProperty::UpwindLeft, double>(interpolation_array, 1.0);
                  derivatives[2][1] = DifferentiationUtilities::ApplyStencil<ReconstructionStencilConcretization, StencilProperty::UpwindRight, double>(interpolation_array, 1.0);
               }

               double const old_phi_sign = Signum(levelset_0_orig[i][j][k]);
               double const godunov_hamiltonian = GetGodunovHamiltonian(derivatives, old_phi_sign);
               double const advection_velocity = ComputeAdvectionVelocity(levelset_0_orig[i][j][k], godunov_hamiltonian);
               double const increment = ReinitializationConstants::Dtau * advection_velocity * (1.0 - godunov_hamiltonian );
               if(ReinitializationConstants::TrackConvergence && std::abs(interface_tags[i][j][k]) < ITTI(IT::ReinitializationBand)) {
                  residuum = std::max(residuum, std::abs(increment));
               }

               reinitialization_rhs[i][j][k] = increment;
            }
         } //k
      } //j
   } //i

   //add to level-set and also the residuum
   for(unsigned int i = CC::FICX(); i <= CC::LICX(); ++i) {
      for(unsigned int j = CC::FICY(); j <= CC::LICY(); ++j) {
         for(unsigned int k = CC::FICZ(); k <= CC::LICZ(); ++k) {
            levelset_orig[i][j][k] += reinitialization_rhs[i][j][k];
         } //k
      } //j
   } //i

   return residuum;
}


