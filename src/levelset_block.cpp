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
#include "levelset_block.h"
#include <stdexcept>

/**
 * @brief Constructor to create a levelset according to an already computed levelset field.
 * @param phi_new Reference to Array holding the levelset field to be in this levelset block.
 */
LevelsetBlock::LevelsetBlock(const double (&phi_new)[CC::TCX()][CC::TCY()][CC::TCZ()]) {
   for(unsigned int i = 0; i < CC::TCX(); ++i) {
      for(unsigned int j = 0; j < CC::TCY(); ++j) {
         for(unsigned int k = 0; k < CC::TCZ(); ++k) {
            phi_reinitialized_[i][j][k] = phi_new[i][j][k];
            phi_rhs_[i][j][k] = phi_new[i][j][k];
            phi_[i][j][k] = 0.0;
            phi_initial_[i][j][k] = 0.0;
         }
      }
   }
   for(unsigned int i = 0; i < CC::TCX(); ++i) {
      for(unsigned int j = 0; j < CC::TCY(); ++j) {
         for(unsigned int k = 0; k < CC::TCZ(); ++k) {
            volume_fraction_[i][j][k] = 0.0;
         }
      }
   }


   for(const InterfaceQuantity iq : FF::ASIQ()) {
      double (&cells)[CC::TCX()][CC::TCY()][CC::TCZ()] = GetInterfaceQuantityBuffer(iq);
      for(unsigned int i = 0; i < CC::TCX(); ++i){
         for(unsigned int j = 0; j < CC::TCY(); ++j){
            for(unsigned int k = 0; k < CC::TCZ(); ++k){
               cells[i][j][k] = 0.0;
            } // k
         } // j
      } // i
   } // interface quantities
}

/**
 * @brief Constructor to create an initial homogenous levelset field.
 * @param initial_phi The value to be imposed as levelset function.
 */
LevelsetBlock::LevelsetBlock(const double initial_phi) {

   const double uniform_volume_fraction = initial_phi > 0 ? 1.0 : 0.0;
   for(unsigned int i = 0; i < CC::TCX(); ++i) {
      for(unsigned int j = 0; j < CC::TCY(); ++j) {
         for(unsigned int k = 0; k < CC::TCZ(); ++k) {
            phi_reinitialized_[i][j][k] = initial_phi;
            phi_rhs_[i][j][k] = initial_phi;
            phi_[i][j][k] = 0.0;
         }
      }
   }
   for(unsigned int i = 0; i < CC::TCX(); ++i) {
      for(unsigned int j = 0; j < CC::TCY(); ++j) {
         for(unsigned int k = 0; k < CC::TCZ(); ++k) {
            volume_fraction_[i][j][k] = uniform_volume_fraction;
         }
      }
   }

   for(const InterfaceQuantity iq : FF::ASIQ()) {
      double (&cells)[CC::TCX()][CC::TCY()][CC::TCZ()] = GetInterfaceQuantityBuffer(iq);
      for(unsigned int i = 0; i < CC::TCX(); ++i){
         for(unsigned int j = 0; j < CC::TCY(); ++j){
            for(unsigned int k = 0; k < CC::TCZ(); ++k){
               cells[i][j][k] = 0.0;
            }
         }
      }
   }
}

/**
 * @brief Gives the phi values in the normal buffer.
 */
auto LevelsetBlock::GetPhi() -> double (&)[CC::TCX()][CC::TCY()][CC::TCZ()] {
   return phi_;
}

/**
 * @brief Gives the phi values in the right hand side buffer.
 */
auto LevelsetBlock::GetPhiRightHandSide() -> double (&)[CC::TCX()][CC::TCY()][CC::TCZ()] {
   return phi_rhs_;
}

/**
 * @brief Gives the phi reinitialized buffer.
 */
auto LevelsetBlock::GetPhiReinitialized() -> double (&)[CC::TCX()][CC::TCY()][CC::TCZ()] {
   return phi_reinitialized_;
}

/**
 * @brief Gives the internal volume fraction buffer.
 */
auto LevelsetBlock::GetVolumeFraction() -> double (&)[CC::TCX()][CC::TCY()][CC::TCZ()] {
   return volume_fraction_;
}

/**
 * @brief Const overload of GetPhi(). See there for details.
 */
auto LevelsetBlock::GetPhi() const -> double const (&)[CC::TCX()][CC::TCY()][CC::TCZ()] {
   return phi_;
}

/**
 * @brief Const overload of GetPhiRightHandSide(). See there for details.
 */
auto LevelsetBlock::GetPhiRightHandSide() const -> double const (&)[CC::TCX()][CC::TCY()][CC::TCZ()] {
   return phi_rhs_;
}

/**
 * @brief Const overlaod of GetPhiReinitialized(). See there for details.
 */
auto LevelsetBlock::GetPhiReinitialized() const -> double const (&)[CC::TCX()][CC::TCY()][CC::TCZ()] {
   return phi_reinitialized_;
}

/**
 * @brief Const overlaod of GetVolumeFraction(). See there for details.
 */
auto LevelsetBlock::GetVolumeFraction() const -> double const (&)[CC::TCX()][CC::TCY()][CC::TCZ()] {
   return volume_fraction_;
}

/**
 * @brief Gives a reference to the corresponding initial phi buffer.
 * @param equation Decider which buffer is to be returned.
 * @return Reference to the array, that is the requested buffer.
 */
auto LevelsetBlock::GetPhiInitial() -> double (&)[CC::TCX()][CC::TCY()][CC::TCZ()] {
   return phi_initial_;
}

/**
 * @brief Const overload.
 */
auto LevelsetBlock::GetPhiInitial() const -> const double (&)[CC::TCX()][CC::TCY()][CC::TCZ()] {
   return phi_initial_;
}

/**
 * @brief Gives the requested phi buffer.
 * @param type .
 */
auto LevelsetBlock::GetPhiBuffer(const PhiType type) -> double (&)[CC::TCX()][CC::TCY()][CC::TCZ()] {
   switch(type) {
      case PhiType::Base:
         return phi_;
      case PhiType::RightHandSide:
         return phi_rhs_;
#ifdef PERFORMANCE
      default:
        return phi_reinitialized_;
#else
      case PhiType::Reinitialized:
         return phi_reinitialized_;
      default:
         throw std::invalid_argument("Requested phi buffer does not exsitst (impossible error)");
         break;
#endif
   }
}

/**
 * @brief Const overload.
 */
auto LevelsetBlock::GetPhiBuffer(const PhiType type) const -> double const (&)[CC::TCX()][CC::TCY()][CC::TCZ()] {
   switch(type) {
      case PhiType::Base:
         return phi_;
      case PhiType::RightHandSide:
         return phi_rhs_;
#ifdef PERFORMANCE
      default:
        return phi_reinitialized_;
#else
      case PhiType::Reinitialized:
         return phi_reinitialized_;
      default:
         throw std::invalid_argument("Requested phi buffer does not exist (impossible error)");
         break;
#endif
   }
}

/**
 * @brief Gives the requested buffer of the levelset block.
 * @param type .
 */
auto LevelsetBlock::GetBuffer(const LevelsetBlockBufferType type) -> double (&)[CC::TCX()][CC::TCY()][CC::TCZ()] {
   switch(type) {
      case LevelsetBlockBufferType::PhiBase:
         return phi_;
      case LevelsetBlockBufferType::PhiRightHandSide:
         return phi_rhs_;
      case LevelsetBlockBufferType::PhiReinitialized:
         return phi_reinitialized_;
      case LevelsetBlockBufferType::InterfaceVelocity:
         return interface_quantities_[InterfaceQuantity::Velocity];
      case LevelsetBlockBufferType::InterfacePressurePositive:
         return interface_quantities_[InterfaceQuantity::PressurePositive];

      case LevelsetBlockBufferType::InterfacePressureNegative:
         return  interface_quantities_[InterfaceQuantity::PressureNegative];

#ifdef PERFORMANCE
      default:
        return volume_fraction_;
#else
      case LevelsetBlockBufferType::VolumeFraction:
         return volume_fraction_;
      default:
         throw std::invalid_argument("Requested buffer does not exist (impossible error");
         break;
#endif
   }
}

/**
 * @brief Const overload.
 */
auto LevelsetBlock::GetBuffer(const LevelsetBlockBufferType type) const -> double const (&)[CC::TCX()][CC::TCY()][CC::TCZ()] {
   switch(type) {
      case LevelsetBlockBufferType::PhiBase:
         return phi_;
      case LevelsetBlockBufferType::PhiRightHandSide:
         return phi_rhs_;
      case LevelsetBlockBufferType::PhiReinitialized:
         return phi_reinitialized_;
      case LevelsetBlockBufferType::InterfaceVelocity:
         return interface_quantities_[InterfaceQuantity::Velocity];
      case LevelsetBlockBufferType::InterfacePressurePositive:
         return interface_quantities_[InterfaceQuantity::PressurePositive];
      case LevelsetBlockBufferType::InterfacePressureNegative:
         return  interface_quantities_[InterfaceQuantity::PressureNegative];
#ifdef PERFORMANCE
      default:
        return volume_fraction_;
#else
      case LevelsetBlockBufferType::VolumeFraction:
         return volume_fraction_;
      default:
         throw std::invalid_argument("Requested buffer does not exist (impossible error)");
         break;
#endif
   }
}

/**
 * @brief Gives a reference to the corresponding interface quantity buffer.
 * @param iq Decider which buffer is to be returned.
 * @return Reference to the array, that is the requested buffer.
 */
auto LevelsetBlock::GetInterfaceQuantityBuffer(const InterfaceQuantity iq) -> double (&)[CC::TCX()][CC::TCY()][CC::TCZ()] {
   return interface_quantities_[iq];
}

/**
 * @brief Const overload.
 */
auto LevelsetBlock::GetInterfaceQuantityBuffer(const InterfaceQuantity iq) const -> const double (&)[CC::TCX()][CC::TCY()][CC::TCZ()] {
   return interface_quantities_[iq];
}
