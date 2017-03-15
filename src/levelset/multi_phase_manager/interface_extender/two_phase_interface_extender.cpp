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
#include "two_phase_interface_extender.h"


#include <vector>
#include <functional>

#include "user_specifications/two_phase_constants.h"
#include "enums/interface_tag_definition.h"
#include "mathematical_functions.h"

/**
 * @brief Default constructor of the TwoPhaseInterfaceExtender class.
 * @param halo_manager Instance of a HaloManager which provides MPI-related methods.
 */
TwoPhaseInterfaceExtender::TwoPhaseInterfaceExtender( HaloManager& halo_manager ) : InterfaceExtender( halo_manager )
{
   // Empty Constructor, besides call of base class constructor.
}

/**
 * @brief Extends interface quantities into the narrow band.
 * @param nodes The nodes for which scale separation is done.
 */
void TwoPhaseInterfaceExtender::ExtendInterfaceQuantitiesImplementation(std::vector<std::reference_wrapper<Node>> const& nodes) const {

   std::vector<double> convergence_tracking_quantities(number_of_convergence_tracking_quantities_,0.0);
   for(unsigned int iteration_number = 0; iteration_number < InterfaceQuantityExtensionConstants::MaximumNumberOfIterations; ++iteration_number) {
      if constexpr(InterfaceQuantityExtensionConstants::TrackConvergence) {
         for(unsigned int i = 0; i < FF::NOIQTE(); ++i) {
            convergence_tracking_quantities[i] = 0.0;
         }
         for(auto const& node : nodes) {
            DetermineMaximumInterfaceQuantityValues(node, convergence_tracking_quantities);
         }

         MPI_Allreduce(MPI_IN_PLACE, convergence_tracking_quantities.data(), number_of_convergence_tracking_quantities_, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

         if(convergence_tracking_quantities[FF::NOIQTE()] < InterfaceQuantityExtensionConstants::MaximumResiduum && iteration_number != 0) {
            if constexpr(GeneralTwoPhaseSettings::LogConvergenceInformation) {
               logger_.AppendDelayedLog("IntExt: " + std::to_string(static_cast<int>(iteration_number)) + " ");
            }
            break;
         } else if(iteration_number == InterfaceQuantityExtensionConstants::MaximumNumberOfIterations - 1) {
            if constexpr(GeneralTwoPhaseSettings::LogConvergenceInformation) {
               logger_.AppendDelayedLog("IntExt: nc   !!!   ");
            }
         }
         convergence_tracking_quantities[FF::NOIQTE()] = 0.0;
      }

      for(auto const& node : nodes) {
         IterativeInterfaceExtension(node, convergence_tracking_quantities);
      }

      for(const InterfaceQuantity iq : FF::IQTE()) {
         halo_manager_.LevelsetHaloUpdateOnLmax( MapInterfaceQuantityTypeToLevelsetBlockBufferType( iq ) );
      }
   }
}

/**
 * @brief Extend interface velocity into narrow band.
 * @param node The node which contains the phase which is extended.
 * @param convergence_tracking_quantities An array holding information about the convergence status of the iterative extension method.
 */
void TwoPhaseInterfaceExtender::IterativeInterfaceExtension(Node& node, std::vector<double>& convergence_tracking_quantities) const {

   double const   (&phi_reinitialized)[CC::TCX()][CC::TCY()][CC::TCZ()] = node.GetLevelsetBlock().GetPhiReinitialized();
   std::int8_t const (&interface_tags)[CC::TCX()][CC::TCY()][CC::TCZ()] = node.GetInterfaceTags();

   std::array<double, DTI(CC::DIM())> rhs_contributions;

   double interface_quantity_change[CC::TCX()][CC::TCY()][CC::TCZ()];
   unsigned int interface_quantity_index = 0;

   for(InterfaceQuantity iq : FF::IQTE()) {
      double const one_normalization_constant = 1.0 / (convergence_tracking_quantities[interface_quantity_index] + epsilon_);
      double (&interface_quantity)[CC::TCX()][CC::TCY()][CC::TCZ()] = node.GetLevelsetBlock().GetInterfaceQuantityBuffer(iq);

      for(unsigned int i = 0; i < CC::TCX(); ++i){
         for(unsigned int j = 0; j < CC::TCY(); ++j){
            for(unsigned int k = 0; k < CC::TCZ(); ++k){
               interface_quantity_change[i][j][k] = 0.0;
            } //k
         } //j
      } //i

      //Loop through all narrow-band cells which are no cut-cells to compute increment for iterative reinitialization
      for(unsigned int i = CC::FICX(); i <= CC::LICX(); ++i) {
         for(unsigned int j = CC::FICY(); j <= CC::LICY(); ++j) {
            for(unsigned int k = CC::FICZ(); k <= CC::LICZ(); ++k) {
               /**
                 * JW: We also extend in the reinitialization band in order to have better convergence behaviour (Reduce influence of implicitly imposed boundary
                 * conditions at the end of the narrow band). The convergence criteria is only checked for the extension band.
                 */
               if(std::abs(interface_tags[i][j][k]) <= ITTI(IT::ReinitializationBand) && std::abs(interface_tags[i][j][k]) > ITTI(IT::NewCutCell)) {

                  // compute normal
                  std::array<double, 3> const normal = GetNormal(phi_reinitialized,i,j,k);

                  double const ls_sign = Signum(phi_reinitialized[i][j][k]);

                  // compute extension terms and store in vector
                  // i-direction
                  // direction of the gradient from interface towards the fluid in both directions
                  int const n_i = (normal[0] * ls_sign < 0.0) - (normal[0] * ls_sign > 0.0);
                  rhs_contributions[0] = (interface_quantity[i+n_i][j][k] - interface_quantity[i][j][k] );
                  rhs_contributions[0] *= std::abs( normal[0] );

                  // j-direction
                  if( CC::DIM() != Dimension::One ){
                     int const n_j = (normal[1] * ls_sign < 0.0) - (normal[1] * ls_sign > 0.0);
                     rhs_contributions[1] = (interface_quantity[i]    [j+n_j][k] - interface_quantity[i][j][k]);
                     rhs_contributions[1] *= std::abs( normal[1] );
                  }

                  // k-direction
                  if( CC::DIM() == Dimension::Three ){
                     int const n_k = (normal[2] * ls_sign < 0.0) - (normal[2] * ls_sign > 0.0);
                     rhs_contributions[2] = (interface_quantity[i]    [j][k+n_k] - interface_quantity[i][j][k]);
                     rhs_contributions[2] *= std::abs( normal[2] );
                  }

                  interface_quantity_change[i][j][k] = ConsistencyManagedSum(rhs_contributions) * InterfaceQuantityExtensionConstants::Dtau;
                  if(InterfaceQuantityExtensionConstants::TrackConvergence && std::abs(interface_tags[i][j][k]) <= ITTI(IT::ExtensionBand)) {
                     convergence_tracking_quantities[FF::NOIQTE()] = std::max(convergence_tracking_quantities[FF::NOIQTE()], std::abs(interface_quantity_change[i][j][k] * one_normalization_constant ));
                  }
               } //if cut cell
            } //k
         } //j
      } //i

      for(unsigned int i = CC::FICX(); i <= CC::LICX(); ++i) {
         for(unsigned int j = CC::FICY(); j <= CC::LICY(); ++j) {
            for(unsigned int k = CC::FICZ(); k <= CC::LICZ(); ++k) {
               interface_quantity[i][j][k] += interface_quantity_change[i][j][k];
            } //k
         } //j
      } //i

      interface_quantity_index += 1;
   } //interface quantity
}

/**
 * @brief Determines the maximum interface quantities for a node. Necessary to track the convergence during the
 * extrapolation of the interface quantities.
 * @param node The node for which the maximum number of interface quantities is determined.
 * @param convergence_tracking_quantities The indirect return parameter.
 */
void TwoPhaseInterfaceExtender::DetermineMaximumInterfaceQuantityValues(Node const& node, std::vector<double>& convergence_tracking_quantities) const {

   std::int8_t const (&interface_tags)[CC::TCX()][CC::TCY()][CC::TCZ()] = node.GetInterfaceTags();
   unsigned int interface_quantity_index = 0;

   for(InterfaceQuantity iq : FF::IQTE()) {
      double const (&interface_quantity)[CC::TCX()][CC::TCY()][CC::TCZ()] = node.GetLevelsetBlock().GetInterfaceQuantityBuffer(iq);
      for(unsigned int i = CC::FICX(); i <= CC::LICX(); ++i) {
         for(unsigned int j = CC::FICY(); j <= CC::LICY(); ++j) {
            for(unsigned int k = CC::FICZ(); k <= CC::LICZ(); ++k) {
               if(std::abs(interface_tags[i][j][k]) < ITTI(IT::BulkPhase) && std::abs(interface_tags[i][j][k]) > ITTI(IT::NewCutCell)) {
                  convergence_tracking_quantities[interface_quantity_index] = std::max(convergence_tracking_quantities[interface_quantity_index], std::abs(interface_quantity[i][j][k]));
               }
            } //k
         } //j
      } //i
      interface_quantity_index += 1;
   } //interface quantity
}
