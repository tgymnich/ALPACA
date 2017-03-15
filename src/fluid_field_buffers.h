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
#ifndef FLUID_FIELD_BUFFER_H
#define FLUID_FIELD_BUFFER_H

#include "user_specifications/compile_time_constants.h"
#include "fluid_fields_definitions.h"

/**
 * @brief Bundles buffers for fluid fields of a certain type to have them contiguous in memory and allows accessing each field separately.
 * @tparam N Number of fields.
 * @tparam FieldEnum Enumeration type allowing to access the fluid fields.
 * @tparam int(*const FieldToIndex)(FieldEnum) Function converting the field enumeration to an index in the range [0;N).
 */
template<std::size_t N, typename FieldEnum, unsigned int(*const FieldToIndex)(FieldEnum)>
struct FluidFieldBuffer {
   std::array<double[CC::TCX()][CC::TCY()][CC::TCZ()], N> Fields;

   /**
    * @brief Access the buffer corresponding to field f.
    */
   auto operator[]( FieldEnum const f ) -> double(&)[CC::TCX()][CC::TCY()][CC::TCZ()] {
      return Fields[FieldToIndex(f)];
   }

   /**
    * @brief Access the buffer corresponding to field f. Const overload.
    */
   auto operator[]( FieldEnum const f ) const -> double const(&)[CC::TCX()][CC::TCY()][CC::TCZ()] {
      return Fields[FieldToIndex(f)];
   }

   /**
    * @brief Access the buffer at the given index.
    */
   auto operator[]( unsigned short const index ) -> double(&)[CC::TCX()][CC::TCY()][CC::TCZ()] {
      return Fields[index];
   }

   /**
    * @brief Access the buffer at the given index. Const overload.
    */
   auto operator[]( unsigned short const index ) const -> double const(&)[CC::TCX()][CC::TCY()][CC::TCZ()] {
      return Fields[index];
   }

   /**
    * @brief Allows to get and set the field values in a certain fluid cell.
    */
   struct CellView {
      FluidFieldBuffer& buffer_;
      unsigned int const i, j, k;

      double& operator[]( FieldEnum const f ) {
         return buffer_[f][i][j][k];
      }

      double& operator[]( unsigned short const index ) {
         return buffer_[index][i][j][k];
      }
   };

   /**
    * @brief Allows to get the field values in a certain fluid cell. Read-only!
    */
   struct CellViewConst {
      FluidFieldBuffer const& buffer_;
      unsigned int const i, j, k;

      double operator[]( FieldEnum const f ) const {
         return buffer_[f][i][j][k];
      }

      double operator[]( unsigned short const index ) const {
         return buffer_[index][i][j][k];
      }
   };

   /**
    * @brief Get a non-modifiable object representing the fluid cell at the given indices.
    * @param i, j, k Fluid cell index.
    * @return A constant view of the field values in the cell.
    */
   auto GetCellView( unsigned int const i, unsigned int const j, unsigned int const k ) const {
      return CellViewConst{ *this, i, j, k };
   }

   /**
    * @brief Get a modifiable object representing the fluid cell at the given indices.
    * @param i, j, k Fluid cell index.
    * @return A handle to the field values in the cell.
    */
   auto GetCellView( unsigned int const i, unsigned int const j, unsigned int const k ) {
      return CellView{ *this, i, j, k };
   }
};

/**
 * @brief Bundles the conservative values to have them contiguous in memory.
 */
using Conservatives = FluidFieldBuffer<FF::ANOE(), Equation, ETI>;
// Check Memory Layout at compile time for safe MPI sending (Ensures Compiler did not pad the struct)
static_assert( sizeof(Conservatives) == FF::ANOE()*CC::TCX()*CC::TCY()*CC::TCZ()*sizeof(double), "Conservative Struct is not contiguous in Memory" );

/**
 * @brief Bundles the prime state values to have them contiguous in memory.
 */
using PrimeStates = FluidFieldBuffer<FF::ANOP(), PrimeState, PTI>;
// Check Memory Layout at compile time for safe MPI sending (Ensures Compiler did not pad the using)
static_assert( sizeof(PrimeStates) == FF::ANOP()*CC::TCX()*CC::TCY()*CC::TCZ()*sizeof(double), "Prime State Struct is not contiguous in Memory" );

/**
 * @brief Bundles the interface quantities to have them contiguous in memory.
 */
using InterfaceQuantities = FluidFieldBuffer<FF::ANIQ(), InterfaceQuantity, IQTI>;
static_assert( sizeof(InterfaceQuantities) == FF::ANIQ()*CC::TCX()*CC::TCY()*CC::TCZ()*sizeof(double), "InterfaceQuantities Struct is not contiguous in Memory" );

#endif // FLUID_FIELD_BUFFER_H
