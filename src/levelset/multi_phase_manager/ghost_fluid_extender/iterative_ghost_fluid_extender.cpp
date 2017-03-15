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
#include "iterative_ghost_fluid_extender.h"

#include "user_specifications/two_phase_constants.h"
#include "levelset/multi_phase_manager/two_phase_manager.h"

/**
 * @brief The default constructor for the IterativeGhostFluidExtender. Calls the default constructor of the base class.
 * @param material_manager Instance of a material manager, which already has been inialized according to the user input.
 * @param halo_manager Instance to a HaloManager which provides MPI-related methods.
 */
IterativeGhostFluidExtender::IterativeGhostFluidExtender( MaterialManager const& material_manager, HaloManager& halo_manager ) :
   GhostFluidExtender( material_manager, halo_manager )
{
   // Empty Constructor, besides call of base class constructor.
}

namespace {

constexpr unsigned int number_of_convergence_tracking_quantities = FF::ANOP() + 1;

/**
 * @brief Determines the maximum prime states in the cells where we extend. This is necessary to have a global normalization constant for convergence tracking.
 * @param node The node for which the maximum prime-states are determined.
 * @param convergence_tracking_quantities An array holding information about the convergence status of the iterative extension method.
 */
void DetermineMaximumValueOfQuantitiesToExtend(Node const& node, double (&convergence_tracking_quantities)[2][number_of_convergence_tracking_quantities]){

   std::int8_t const (&interface_tags)[CC::TCX()][CC::TCY()][CC::TCZ()] = node.GetInterfaceTags();
   double const (&phi_reinitialized)[CC::TCX()][CC::TCY()][CC::TCZ()] = node.GetLevelsetBlock().GetPhiReinitialized();
   double const (&volume_fraction)[CC::TCX()][CC::TCY()][CC::TCZ()] = node.GetLevelsetBlock().GetVolumeFraction();

   for(auto const& phase : node.GetPhases()) {

      std::int8_t const material_sign = MaterialSignCapsule::SignOfMaterial(phase.first);
      unsigned int const material_index = phase.first == MaterialSignCapsule::PositiveFluidMaterial() ? 0 : 1;

      double const reference_volume_fraction = (material_sign > 0) ? 0.0 : 1.0;
      double const material_sign_double = double(material_sign);

      for(const PrimeState p : FF::ASOP()) {
         double const (&prime_state)[CC::TCX()][CC::TCY()][CC::TCZ()] = phase.second.GetPrimeStateBuffer(p);
         for(unsigned int i = CC::FICX(); i <= CC::LICX(); ++i) {
            for(unsigned int j = CC::FICY(); j <= CC::LICY(); ++j) {
               for(unsigned int k = CC::FICZ(); k <= CC::LICZ(); ++k) {
                  double const cell_volume_fraction = reference_volume_fraction + material_sign_double * volume_fraction[i][j][k];
                  if(std::abs(interface_tags[i][j][k]) <= ITTI(IT::ReinitializationBand) && (material_sign * phi_reinitialized[i][j][k] < 0.0 || cell_volume_fraction < CC::MITH())) {
                     convergence_tracking_quantities[material_index][PTI(p)] = std::max(convergence_tracking_quantities[material_index][PTI(p)], std::abs(prime_state[i][j][k]));
                  }
               } // k
            } // j
         } // i
      } // prime states
   } // phases
}

/**
 * @brief Extend to cut-cell neighbours and extension band iteratively.
 * @param node The node which is extended.
 * @param convergence_tracking_quantities An array holding information about the convergence status of the iterative extension method.
 */
void IterativeExtension(Node& node, double (&convergence_tracking_quantities)[2][number_of_convergence_tracking_quantities]){



   double const (&phi_reinitialized)[CC::TCX()][CC::TCY()][CC::TCZ()] = node.GetLevelsetBlock().GetPhiReinitialized();
   std::int8_t const (&interface_tags)[CC::TCX()][CC::TCY()][CC::TCZ()] = node.GetInterfaceTags();
   double const (&volume_fraction)[CC::TCX()][CC::TCY()][CC::TCZ()] = node.GetLevelsetBlock().GetVolumeFraction();

   constexpr double epsilon = std::numeric_limits<double>::epsilon();

   for(auto& phase : node.GetPhases()) {

      std::int8_t const material_sign = MaterialSignCapsule::SignOfMaterial(phase.first);
      unsigned int const material_index = phase.first == MaterialSignCapsule::PositiveFluidMaterial() ? 0 : 1;

      double const reference_volume_fraction = (material_sign > 0) ? 0.0 : 1.0;
      double const material_sign_double = double(material_sign);

      double extension_rhs[FF::ANOP()][CC::TCX()][CC::TCY()][CC::TCZ()];
      double one_normalization_constant[FF::ANOP()];
      for(unsigned int p = 0; p < FF::ANOP(); ++p) {
         one_normalization_constant[p] = 1.0 / (std::max(epsilon, convergence_tracking_quantities[material_index][p]));
      }

      /**
       * Setting the extension_rhs buffer to zero is crucial!
       */
      for(const PrimeState p : FF::ASOP()) {
         for(unsigned int i = 0; i < CC::TCX(); ++i) {
            for(unsigned int j = 0; j < CC::TCY(); ++j) {
               for(unsigned int k = 0; k < CC::TCZ(); ++k) {
                  extension_rhs[PTI(p)][i][j][k] = 0.0;
               } // k
            } // j
         } // i
      } // prime states

      std::array<double, DTI(CC::DIM())> rhs_contributions;

      unsigned int derivative_indices[DTI(CC::DIM())][2];
      for(unsigned int d = 0; d < DTI(CC::DIM()); ++d) {
         for(unsigned int l = 0; l < 2; ++l) {
            derivative_indices[d][l] = 0;
         }
      }

      // Loop through internal block - finally, fill extension band and cut-cell neighbours
      for(unsigned int i = CC::FICX(); i <= CC::LICX(); ++i) {
         for(unsigned int j = CC::FICY(); j <= CC::LICY(); ++j) {
            for(unsigned int k = CC::FICZ(); k <= CC::LICZ(); ++k) {
               double const cell_volume_fraction = reference_volume_fraction + material_sign_double * volume_fraction[i][j][k];
               /**
                * JW: We also extend in the reinitialization band in order to have better convergence behaviour (Reduce influence of implicitly imposed boundary
                * conditions at the end of the narrow band). The convergence criteria is only checked for the extension band.
                */
               if(std::abs(interface_tags[i][j][k]) <= ITTI(IT::ExtensionBand) && (material_sign * phi_reinitialized[i][j][k] < 0.0 || cell_volume_fraction < CC::MITH())) {

                  std::array<double, 3> const normal = GetNormal(phi_reinitialized, i, j, k, material_sign);

                  if(normal[0] > 0.0) {
                     derivative_indices[0][0] = i + 1;
                     derivative_indices[0][1] = i;
                  } else {
                     derivative_indices[0][0] = i;
                     derivative_indices[0][1] = i - 1;
                  }

                  if constexpr(CC::DIM() != Dimension::One) {
                     if(normal[1] > 0.0) {
                        derivative_indices[1][0] = j + 1;
                        derivative_indices[1][1] = j;
                     } else {
                        derivative_indices[1][0] = j;
                        derivative_indices[1][1] = j - 1;
                     }
                  }

                  if constexpr(CC::DIM() == Dimension::Three) {
                     if(normal[2] > 0.0) {
                        derivative_indices[2][0] = k + 1;
                        derivative_indices[2][1] = k;
                     } else {
                        derivative_indices[2][0] = k;
                        derivative_indices[2][1] = k - 1;
                     }
                  }

                  //calculate gradients
                  for(const PrimeState p : FF::ASOP()) {
                     double const (&cell)[CC::TCX()][CC::TCY()][CC::TCZ()] = phase.second.GetPrimeStateBuffer(p);

                     rhs_contributions[0] = cell[derivative_indices[0][0]][j][k] - cell[derivative_indices[0][1]][j][k];
                     rhs_contributions[0] *= normal[0];

                     if constexpr(CC::DIM() != Dimension::One) {
                        rhs_contributions[1] = cell[i][derivative_indices[1][0]][k] - cell[i][derivative_indices[1][1]][k];
                        rhs_contributions[1] *= normal[1];
                     }

                     if constexpr(CC::DIM() == Dimension::Three) {
                        rhs_contributions[2] = cell[i][j][derivative_indices[2][0]] - cell[i][j][derivative_indices[2][1]];
                        rhs_contributions[2] *= normal[2];
                     }

                     extension_rhs[PTI(p)][i][j][k] = ConsistencyManagedSum(rhs_contributions) * ExtensionConstants::Dtau;

                     if(ExtensionConstants::TrackConvergence && std::abs(interface_tags[i][j][k]) <= ITTI(IT::ExtensionBand)) {
                        convergence_tracking_quantities[material_index][FF::ANOP()] = std::max(convergence_tracking_quantities[material_index][FF::ANOP()], std::abs(extension_rhs[PTI(p)][i][j][k] * one_normalization_constant[PTI(p)]));
                     }

                  } // prime states
               } // cells to extend
            } // k
         } // j
      } // i

      for(const PrimeState p : FF::ASOP()) {
         double (&prime_state)[CC::TCX()][CC::TCY()][CC::TCZ()] = phase.second.GetPrimeStateBuffer(p);
         for(unsigned int i = CC::FICX(); i <= CC::LICX(); ++i) {
            for(unsigned int j = CC::FICY(); j <= CC::LICY(); ++j) {
               for(unsigned int k = CC::FICZ(); k <= CC::LICZ(); ++k) {
                  prime_state[i][j][k] += extension_rhs[PTI(p)][i][j][k];
               } // k
            } // j
         } // i
      } // prime states

   }
}

}

/**
 * @brief Iteratively solves the extension equation. For description of functionality also see base class.
 * @param nodes The nodes for which the extension equation is solved.
 * @param stage The current stage of the Runge-Kutta method.
 */
void IterativeGhostFluidExtender::ExtendImplementation(std::vector<std::reference_wrapper<Node>> const& nodes, unsigned int const stage) const {
#ifndef PERFORMANCE
   (void) stage; // Avoid compiler warning
#endif

   constexpr unsigned int number_of_convergence_tracking_quantities = FF::ANOP() + 1;

   //The 2 is hardcoded on purpose. It corresponds to the number of fluids. An issue about that is already in the git.
   double convergence_tracking_quantities[2][number_of_convergence_tracking_quantities];
   for(unsigned int p = 0; p < 2; ++p) {
      for(unsigned int i = 0; i <= FF::ANOP(); ++i) {
         convergence_tracking_quantities[p][i] = 0.0;
      }
   }
   for(unsigned int iteration_number = 0; iteration_number < ExtensionConstants::MaximumNumberOfIterations; ++iteration_number) {

      if constexpr(ExtensionConstants::TrackConvergence) {
         for(unsigned int p = 0; p < 2; ++p) {
            for(unsigned int i = 0; i < FF::ANOP(); ++i) {
               convergence_tracking_quantities[p][i] = 0.0;
            }
         }
         for(Node const& node : nodes) {
            DetermineMaximumValueOfQuantitiesToExtend(node, convergence_tracking_quantities);
         }

         MPI_Allreduce(MPI_IN_PLACE, &convergence_tracking_quantities, number_of_convergence_tracking_quantities * 2, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

         if(convergence_tracking_quantities[0][FF::ANOP()] < ExtensionConstants::MaximumResiduum && convergence_tracking_quantities[1][FF::ANOP()] < ExtensionConstants::MaximumResiduum && iteration_number != 0) {
            if constexpr(GeneralTwoPhaseSettings::LogConvergenceInformation) {
               logger_.AppendDelayedLog("Ext: " + std::to_string(static_cast<int>(iteration_number)) + " ");
            }
            break;
         } else if(iteration_number == ExtensionConstants::MaximumNumberOfIterations - 1) {
            if constexpr(GeneralTwoPhaseSettings::LogConvergenceInformation) {
               logger_.AppendDelayedLog("Ext: nc   !!!   ");
            }
         }
         convergence_tracking_quantities[0][FF::ANOP()] = 0.0;
         convergence_tracking_quantities[1][FF::ANOP()] = 0.0;
      }

      //iterative extension on prime state buffer
      for(Node& node : nodes) {
         IterativeExtension(node, convergence_tracking_quantities);
      } //nodes

      halo_manager_.FluidHaloUpdateOnLmax(FluidFieldType::PrimeStates);
   }
}
