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
#ifndef LEVELSET_BLOCK
#define LEVELSET_BLOCK

#include "user_specifications/compile_time_constants.h"
#include "fluid_field_buffers.h"


/**
 * @brief Identifier of the differen levelset buffer.
 */
enum class PhiType : unsigned short {Base, RightHandSide, Reinitialized};

/**
 * @brief Identifier of all buffers in a levelset block.
 */
enum class LevelsetBlockBufferType : unsigned short {PhiBase, PhiRightHandSide, PhiReinitialized, VolumeFraction, InterfaceVelocity, InterfacePressurePositive, InterfacePressureNegative};

constexpr LevelsetBlockBufferType MapInterfaceQuantityTypeToLevelsetBlockBufferType( InterfaceQuantity const iq ) {
   switch( iq ) {
      case InterfaceQuantity::Velocity: {
         return LevelsetBlockBufferType::InterfaceVelocity;
      }
      case InterfaceQuantity::PressurePositive: {
         return LevelsetBlockBufferType::InterfacePressurePositive;
      }
      // Last possibility
      default: {
        return LevelsetBlockBufferType::InterfacePressureNegative;
      }
   }
}

/**
 * @brief The LevelsetBlock class holds the levelset data aka the scalar levelset and the matieral indicators for one node. Does NOT manipulate the data itself,
 *        but provides access to the data.
 */
class LevelsetBlock {

   double                 phi_[CC::TCX()][CC::TCY()][CC::TCZ()];
   double             phi_rhs_[CC::TCX()][CC::TCY()][CC::TCZ()];
   double   phi_reinitialized_[CC::TCX()][CC::TCY()][CC::TCZ()];
   double         phi_initial_[CC::TCX()][CC::TCY()][CC::TCZ()];
   double     volume_fraction_[CC::TCX()][CC::TCY()][CC::TCZ()];

   InterfaceQuantities interface_quantities_;

public:
   LevelsetBlock() = delete;
   explicit LevelsetBlock(const double initial_phi);
   explicit LevelsetBlock(const double (&phi_new)[CC::TCX()][CC::TCY()][CC::TCZ()]);
   ~LevelsetBlock() = default;
   LevelsetBlock( LevelsetBlock const& ) = delete;
   LevelsetBlock& operator=( LevelsetBlock const& ) = delete;
   LevelsetBlock( LevelsetBlock&& ) = delete;
   LevelsetBlock& operator=( LevelsetBlock&& ) = delete;

   auto GetPhi() -> double (&)[CC::TCX()][CC::TCY()][CC::TCZ()];
   auto GetPhiRightHandSide() -> double (&)[CC::TCX()][CC::TCY()][CC::TCZ()];
   auto GetPhiReinitialized() -> double (&)[CC::TCX()][CC::TCY()][CC::TCZ()];
   auto GetPhiInitial() -> double (&)[CC::TCX()][CC::TCY()][CC::TCZ()];
   auto GetVolumeFraction() -> double (&)[CC::TCX()][CC::TCY()][CC::TCZ()];

   auto GetPhi() const -> double const (&)[CC::TCX()][CC::TCY()][CC::TCZ()];
   auto GetPhiRightHandSide() const -> double const (&)[CC::TCX()][CC::TCY()][CC::TCZ()];
   auto GetPhiReinitialized() const -> double const (&)[CC::TCX()][CC::TCY()][CC::TCZ()];
   auto GetPhiInitial() const -> double const (&)[CC::TCX()][CC::TCY()][CC::TCZ()];
   auto GetVolumeFraction() const -> double const (&)[CC::TCX()][CC::TCY()][CC::TCZ()];

   auto GetPhiBuffer(const PhiType type) -> double (&)[CC::TCX()][CC::TCY()][CC::TCZ()];
   auto GetPhiBuffer(const PhiType type) const -> double const (&)[CC::TCX()][CC::TCY()][CC::TCZ()];

   auto GetBuffer(const LevelsetBlockBufferType type) -> double (&)[CC::TCX()][CC::TCY()][CC::TCZ()];
   auto GetBuffer(const LevelsetBlockBufferType type) const -> double const (&)[CC::TCX()][CC::TCY()][CC::TCZ()];

   auto GetInterfaceQuantityBuffer(const InterfaceQuantity iq) -> double (&)[CC::TCX()][CC::TCY()][CC::TCZ()];
   auto GetInterfaceQuantityBuffer(const InterfaceQuantity iq) const -> double const (&)[CC::TCX()][CC::TCY()][CC::TCZ()];

};

#endif // LEVELSET_BLOCK
