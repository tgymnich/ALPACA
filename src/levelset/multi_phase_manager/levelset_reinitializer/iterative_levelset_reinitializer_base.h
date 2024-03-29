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
#ifndef ITERATIVE_LEVELSET_REINITIALIZER_BASE_H
#define ITERATIVE_LEVELSET_REINITIALIZER_BASE_H

#include "levelset_reinitializer.h"
#include "user_specifications/two_phase_constants.h"
#include "mathematical_functions.h"
#include "enums/interface_tag_definition.h"

/**
 * @brief The class IterativeLevelsetReinitializerBase ensures the (signed-)distance property of a level-set field.
 * @tparam Typename as template parameter due to CRTP.
 */
template<typename DerivedIterativeLevelsetReinitializer>
class IterativeLevelsetReinitializerBase : public LevelsetReinitializer<DerivedIterativeLevelsetReinitializer> {

   friend LevelsetReinitializer<DerivedIterativeLevelsetReinitializer>;
   friend DerivedIterativeLevelsetReinitializer;

   using LevelsetReinitializer<DerivedIterativeLevelsetReinitializer>::halo_manager_;
   using LevelsetReinitializer<DerivedIterativeLevelsetReinitializer>::logger_;

   /**
    * @brief The default constructor for a LevelsetReinitializer object.
    * @param halo_manager Instance to a HaloManager which provides MPI-related methods.
    */
   explicit IterativeLevelsetReinitializerBase( HaloManager& halo_manager ) :
      LevelsetReinitializer<DerivedIterativeLevelsetReinitializer>( halo_manager )
   {
      // Empty Constructor, besides initializer list.
   }

   /**
    * @brief Calculate the Godunov Hamiltonian (GH) as described in \cite Min2010.
    * @param derivatives A reference to the single components (level-set derivatives) for which the GH has to be calculated.
    * @param old_levelset_sign The sign of the original level-set value.
    * @return The GH as a double value.
    */
   double GetGodunovHamiltonian(double (&derivatives)[DTI(CC::DIM())][2], double const old_levelset_sign) const {
      std::array<double, DTI(CC::DIM())> godunov_hamiltonian_contributions;

      for(unsigned int d = 0; d < DTI(CC::DIM()); ++d) {
         derivatives[d][0] = std::max( 0.0 , old_levelset_sign * derivatives[d][0] );
         derivatives[d][1] = std::min( 0.0 , old_levelset_sign * derivatives[d][1] );
         godunov_hamiltonian_contributions[d] = std::max(derivatives[d][0]*derivatives[d][0] , derivatives[d][1]*derivatives[d][1]);
      }

      return std::sqrt(ConsistencyManagedSum(godunov_hamiltonian_contributions));
   }

   /**
    * @brief Sets the cut-off in the levelset field of a node .
    * @param node The node which has to be reinitialized.
    */
   void CutOffSingleNode(Node& node) const {

      std::int8_t const  (&interface_tags)[CC::TCX()][CC::TCY()][CC::TCZ()] = node.GetInterfaceTags();
      LevelsetBlock& levelset_block = node.GetLevelsetBlock();
      double (&levelset)[CC::TCX()][CC::TCY()][CC::TCZ()] = levelset_block.GetPhiReinitialized();

      //cells which are outside the reinitialization band are set to cut-off value
      double const cutoff = CC::LSCOF();
      for(unsigned int i = 0; i < CC::TCX(); ++i){
         for(unsigned int j = 0; j < CC::TCY(); ++j){
            for(unsigned int k = 0; k < CC::TCZ(); ++k){
               if(std::abs(interface_tags[i][j][k]) > ITTI(IT::ReinitializationBand)) {
                  levelset[i][j][k] = Signum(levelset[i][j][k]) * cutoff;
               }
            } //k
         } //j
      } //i
   }

   /**
    * @brief Reinitializes a single-level set field as described in \cite Sussman1994.
    * @param nodes See base class.
    * @param stage See base class.
    */
   void ReinitializeImplementation(std::vector<std::reference_wrapper<Node>> const& nodes, unsigned int const stage) const {
#ifndef PERFORMANCE
      (void) stage; // Avoid compiler warning
#endif

      for(Node& node : nodes) {
         LevelsetBlock& levelset_block = node.GetLevelsetBlock();
         double const (&phi_reinitialized)[CC::TCX()][CC::TCY()][CC::TCZ()] = levelset_block.GetPhiReinitialized();
         double (&phi_rhs)[CC::TCX()][CC::TCY()][CC::TCZ()] = levelset_block.GetPhiRightHandSide();

         for(unsigned int i = 0; i < CC::TCX(); ++i) {
            for(unsigned int j = 0; j < CC::TCY(); ++j) {
               for(unsigned int k = 0; k < CC::TCZ(); ++k) {
                  phi_rhs[i][j][k] = phi_reinitialized[i][j][k];
               } //k
            } //j
         } //i
      }

      double residuum = 0.0;
      for(unsigned int iteration_number = 0; iteration_number < ReinitializationConstants::MaximumNumberOfIterations; ++iteration_number){


         if constexpr(ReinitializationConstants::TrackConvergence) {
            residuum = 0.0;
         }
         for(auto& node : nodes) {
            residuum = std::max(residuum, static_cast<DerivedIterativeLevelsetReinitializer const&>(*this).ReinitializeSingleNodeImplementation(node));
         }

         //halo update
         halo_manager_.LevelsetHaloUpdateOnLmax( LevelsetBlockBufferType::PhiReinitialized );
         if constexpr(ReinitializationConstants::TrackConvergence) {
            MPI_Allreduce(MPI_IN_PLACE, &residuum, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

            if(residuum < ReinitializationConstants::MaximumResiduum) {
               if constexpr(GeneralTwoPhaseSettings::LogConvergenceInformation) {
                  logger_.AppendDelayedLog("Reinit: " + std::to_string(static_cast<int>(iteration_number)) + " ");
               }
               break;
            } else if(iteration_number == ReinitializationConstants::MaximumNumberOfIterations - 1) {
               if constexpr(GeneralTwoPhaseSettings::LogConvergenceInformation) {
                  logger_.AppendDelayedLog("Reinit: nc   !!!   ");
               }
            }
         }
      }
   }

public:
   IterativeLevelsetReinitializerBase() = delete;
   virtual ~IterativeLevelsetReinitializerBase() = default;
   IterativeLevelsetReinitializerBase( IterativeLevelsetReinitializerBase const& ) = delete;
   IterativeLevelsetReinitializerBase& operator=( IterativeLevelsetReinitializerBase const& ) = delete;
   IterativeLevelsetReinitializerBase( IterativeLevelsetReinitializerBase&& ) = delete;
   IterativeLevelsetReinitializerBase& operator=( IterativeLevelsetReinitializerBase&& ) = delete;
};


#endif //ITERATIVE_LEVELSET_REINITIALIZER_BASE_H
