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
#ifndef SPATIAL_RECONSTRUCTION_STENCIL_H
#define SPATIAL_RECONSTRUCTION_STENCIL_H

#include <vector>
#include <limits>

#include "enums/direction_definition.h"
#include "stencils/stencil_properties.h"
#include "user_specifications/compile_time_constants.h"

/**
 * @brief Helper functions to evaluate stencils.
 */
namespace {
/**
 * @brief Get the stencil offset and stencil sign based on whether the stencil should be evaluated UpwindLeft, UpwindRight or Central.
 * @param derivative_property Indicates whether the stencil should be evaluated UpwindLeft, UpwindRight or Central.
 * @return The stencil_offset and stencil_sign as an array.
 */
constexpr std::array<int, 2> GetStencilParameters( StencilProperty const derivative_property ) {
   switch( derivative_property ) {
      case StencilProperty::UpwindLeft: return {0, 1};
      case StencilProperty::UpwindRight: return {1, -1};
      default: return {0, 0};
   }
}

/**
 * @brief Returns an array indicating in which direction the stencil should be evaluated.
 * @param direction The direction in which the stencil should be evaluated.
 * @return The array indicating the direction in which the stencil should be evaluated.
 */
constexpr std::array<unsigned short, 3> GetDirectionIndicator( Direction const direction ) {
   switch( direction ) {
      case Direction::X: return {1, 0, 0};
      case Direction::Y: return {0, 1, 0};
      default: return {0, 0, 1};
   }
}
}

/**
 * @brief The SpatialReconstructionStencil class provides an interface for computational Stencils needed throughout the simulation.
 */
template<typename DerivedStencil>
class Stencil {

   friend DerivedStencil;

   explicit Stencil() = default;

   /**
    * @brief Returns the part of a complete buffer which is necessary to evaluate a stencil.
    * @tparam D The direction in which the stencil should be evaluated.
    * @param buffer The buffer on which the stencil should be applied.
    * @param i The index in x-direction.
    * @param j The index in y-direction.
    * @param k The index in z-direction.
    * @return The part of the buffer necessary to evaluate the stencil as a vector.
    */
   template<Direction D>
   std::vector<double> GetVectorFromBuffer( double const (&buffer)[CC::TCX()][CC::TCY()][CC::TCZ()], unsigned int i, unsigned int j, unsigned int k ) const {

      constexpr std::array<unsigned short, 3> direction_indicator = GetDirectionIndicator(D);

      std::vector<double> values;
      values.reserve(StencilSize());
      for( unsigned int s = 0; s < StencilSize(); ++s ) {
         values.push_back(buffer[i + direction_indicator[0] * (-DownstreamStencilSize() + s)][j + direction_indicator[1] * (-DownstreamStencilSize() + s)][k + direction_indicator[2] * (-DownstreamStencilSize() + s)]);
      }
      return values;
   }

protected:
   static constexpr double epsilon_ = std::numeric_limits<double>::epsilon();

public:
   ~Stencil() = default;
   Stencil( Stencil const& ) = delete;
   Stencil& operator=( Stencil const& ) = delete;
   Stencil( Stencil&& ) = delete;
   Stencil& operator=( Stencil&& ) = delete;

   /**
    * @brief Applies the SpatialReconstructionStencil to the provided Array
    * @tparam P The manner in which the stencil is applied (UpwindLeft, UpwindRight or Central).
    * @param array The array on which to apply the spatial reconstruction stencil.
    * @param cell_size The cell size of the corresponding block.
    * @return Value at the position of interest.
    * @note Hotpath function.
    */
   template<StencilProperty P>
   double Apply( const std::vector<double>& array, double const cell_size ) const {
      constexpr std::array<int, 2> stencil_properties = GetStencilParameters( P );
      return static_cast<DerivedStencil const&>(*this).ApplyImplementation( array, stencil_properties[0], stencil_properties[1], cell_size );
   }

   /**
    * @brief Applies the stencil on a given buffer.
    * @tparam P The manner in which the stencil is applied (UpwindLeft, UpwindRight or Central).
    * @tparam D The direction in which the buffer should be applied.
    * @param buffer The buffer on which the stencil is applied.
    * @param i The index in x-direction.
    * @param j The index in y-direction.
    * @param k The index in z-direction.
    * @param cell_size The cell size of the block to which the buffer belongs.
    * @return The result of the stencil.
    */
   template<StencilProperty P, Direction D>
   double Apply( double const (&buffer)[CC::TCX()][CC::TCY()][CC::TCZ()], unsigned int i, unsigned int j, unsigned int k, double const cell_size ) const {
      const std::vector<double> array = GetVectorFromBuffer<D>( buffer, i, j, k );
      constexpr std::array<int, 2> stencil_properties = GetStencilParameters( P );
      return static_cast<DerivedStencil const&>(*this).ApplyImplementation( array, stencil_properties[0], stencil_properties[1], cell_size );
   }

   /**
    * @brief Gives the number of cells needed for a single stencil evaluation.
    * @return Size of the stencil, i.e. number of data cells the stencil works on.
    */
   static constexpr unsigned int StencilSize() {
      return DerivedStencil::stencil_size_;
   }

   /**
    * @brief Gives the size of the stencil in down stream direction.
    * @return Size of the stencil arm reaching down stream, i.e. number of data cells that lay down stream the stencil works on.
    */
   static constexpr unsigned int DownstreamStencilSize() {
      return DerivedStencil::downstream_stencil_size_;
   }

   /**
    * @brief Return the type of a stencil.
    * @return The type of a stencil.
    */
   static constexpr StencilType GetStencilType() {
      return DerivedStencil::stencil_type_;
   }
};

#endif // SPATIAL_RECONSTRUCTION_STENCIL_H
