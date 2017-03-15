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
#ifndef DIFFERENTIATION_UTILITIES_H
#define DIFFERENTIATION_UTILITIES_H

#include <math.h>
#include "mathematical_functions.h"

#include "user_specifications/stencil_setup.h"
#include "user_specifications/compile_time_constants.h"
#include "spatial_reconstruction_stencils/reconstruction_stencil_setup.h"
#include "spatial_derivative_stencils/derivative_stencil_setup.h"

static_assert( ReconstructionStencilSetup::Concretize<reconstruction_stencil>::type::DownstreamStencilSize() < CC::HS(), "Halo size too small for RECONSTRUCTION_STENCIL. Increase the halo size in compile_time_constants.h!" );
static_assert( DerivativeStencilSetup::Concretize<derivative_stencil>::type::DownstreamStencilSize() < CC::HS(), "Halo size too small for DERIVATIVE_STENCIL. Increase the halo size in compile_time_constants.h!" );

/**
 * @brief This namespace provides functionality to calculate differential operators on a field using a given stencil using a specified dimension.
 */
namespace DifferentiationUtilities {

/**
 * @brief      Calculates the gradient at a specified location of a scalar field.
 *
 * @param      buffer     The buffer describing the scalar field.
 * @param[in]  i          The index in x-direction.
 * @param[in]  j          The index in y-direction.
 * @param[in]  k          The index in z-direction.
 * @param[in]  cell_size  The cell size.
 *
 * @tparam     S          The stencil to be used.
 * @tparam     T          The type of the scalar field (for example int or double).
 * @tparam     DIM        The dimension to be used.
 *
 * @return     The gradient as an array.
 */
template<typename S, typename T, Dimension DIM>
std::array<T, 3> Gradient( T const (&buffer)[CC::TCX()][CC::TCY()][CC::TCZ()],
   unsigned int const i, unsigned int const j, unsigned int const k,
   double const cell_size ) {
   const S&& stencil = S();
   return {stencil.template Apply<StencilProperty::Central, Direction::X>( buffer, i, j, k, cell_size ),
           DIM != Dimension::One ? stencil.template Apply<StencilProperty::Central, Direction::Y>( buffer, i, j, k, cell_size ) : 0.0,
           DIM == Dimension::Three ? stencil.template Apply<StencilProperty::Central, Direction::Z>( buffer, i, j, k, cell_size ): 0.0};
}

/**
 * @brief      Calculates the gradient at a specified location of a scalar field.
 *
 * @param      buffer     The buffer describing the scalar field.
 * @param[in]  i          The index in x-direction.
 * @param[in]  j          The index in y-direction.
 * @param[in]  k          The index in z-direction.
 * @param[in]  cell_size  The cell size.
 *
 * @tparam     S          The stencil to be used.
 * @tparam     T          The type of the scalar field (for example int or double).
 *
 * @return     The gradient as an array.
 */
template<typename S, typename T>
std::array<T, 3> ComputeGradient( T const (&buffer)[CC::TCX()][CC::TCY()][CC::TCZ()],
   unsigned int const i, unsigned int const j, unsigned int const k,
   double const cell_size ) {
   return Gradient<S,T,CC::DIM()>( buffer, i, j, k, cell_size );
}

/**
 * @brief      Calculates the derivative at a specified location of a scalar field in a specified direction.
 *
 * @param      buffer     The buffer describing the scalar field.
 * @param[in]  i          The index in x-direction.
 * @param[in]  j          The index in y-direction.
 * @param[in]  k          The index in z-direction.
 * @param[in]  cell_size  The cell size.
 *
 * @tparam     S          The stencil to be used.
 * @tparam     D          The direction in which the derivative should be calculated.
 * @tparam     T          The type of the scalar field (for example int or double).
 *
 * @return     The derivative.
 */
template<typename S, Direction D, typename T>
T ComputeDerivative( T const (&buffer)[CC::TCX()][CC::TCY()][CC::TCZ()],
   unsigned int const i, unsigned int const j, unsigned int const k,
   double const cell_size ) {
   const S&& stencil = S();
   return stencil.template Apply<StencilProperty::Central, D>( buffer, i, j, k, cell_size );
}

/**
 * @brief      Calculates the vector gradient, i.e. the Jacobian of a vector field at a specified location using a specified stencil
 * and a specified dimension number.
 *
 * @param      buffer_x   The x-component of the vector field.
 * @param      buffer_y   The y-component of the vector field.
 * @param      buffer_z   The z-component of the vector field.
 * @param[in]  i          The index in x-direction.
 * @param[in]  j          The index in y-direction.
 * @param[in]  k          The index in z-direction.
 * @param[in]  cell_size  The cell size.
 *
 * @tparam     S          The stencil to be used.
 * @tparam     T          The type of the scalar field (for example int or double).
 * @tparam     DIM        The dimension number to be used.
 *
 * @return     The Jacobian as a two-dimensional array.
 */
template<typename S, typename T, Dimension DIM>
std::array< std::array<T, 3>, 3> VectorGradient( T const (&buffer_x)[CC::TCX()][CC::TCY()][CC::TCZ()],
   T const (&buffer_y)[CC::TCX()][CC::TCY()][CC::TCZ()],
   T const (&buffer_z)[CC::TCX()][CC::TCY()][CC::TCZ()],
   unsigned int const i, unsigned int const j, unsigned int const k,
   double const cell_size ) {
   return { Gradient<S, T, DIM>( buffer_x, i, j, k, cell_size )
          , Gradient<S, T, DIM>( buffer_y, i, j, k, cell_size )
          , Gradient<S, T, DIM>( buffer_z, i, j, k, cell_size ) };
}

/**
 * @brief      Calculates the vector gradient, i.e. the Jacobian of a vector field at a specified location using a specified stencil.
 *
 * @param      buffer_x   The x-component of the vector field.
 * @param      buffer_y   The y-component of the vector field.
 * @param      buffer_z   The z-component of the vector field.
 * @param[in]  i          The index in x-direction.
 * @param[in]  j          The index in y-direction.
 * @param[in]  k          The index in z-direction.
 * @param[in]  cell_size  The cell size.
 *
 * @tparam     S          The stencil to be used.
 * @tparam     T          The type of the scalar field (for example int or double).
 *
 * @return     The Jacobian as a two-dimensional array.
 */
template<typename S, typename T>
std::array< std::array<T, 3>, 3> ComputeVectorGradient( T const (&buffer_x)[CC::TCX()][CC::TCY()][CC::TCZ()],
   T const (&buffer_y)[CC::TCX()][CC::TCY()][CC::TCZ()],
   T const (&buffer_z)[CC::TCX()][CC::TCY()][CC::TCZ()],
   unsigned int const i, unsigned int const j, unsigned int const k,
   double const cell_size ) {
   return { ComputeGradient<S,T>( buffer_x, i, j, k, cell_size )
          , ComputeGradient<S,T>( buffer_y, i, j, k, cell_size )
          , ComputeGradient<S,T>( buffer_z, i, j, k, cell_size )};
}

/**
 * @brief Applies a stencil to a given vector of values.
 * @tparam S The stencil to be used.
 * @tparam P The stencil property to be used (upwind, central...).
 * @tparam T The type of the values of the vector to which the stencil is applied.
 * @param array The vector to which the stencil is applied.
 * @param cell_size The cell size.
 * @return The result of the stencil evaluation.
 */
template<typename S, StencilProperty P, typename T>
T ApplyStencil( const std::vector<double>& array, double const cell_size ) {
   const S&& stencil = S();
   return stencil.template Apply<P>( array, cell_size );
}

/**
 * @brief Applies a stencil to a given vector of values.
 * @tparam S The stencil to be used.
 * @tparam T The type of the values of the vector to which the stencil is applied.
 * @param array The vector to which the stencil is applied.
 * @param upwind_decision Indicates whether left or right upwinding has to be done (if upwind_decision > 0.0 -> upwind left, else vice versa)
 * @param cell_size The cell size.
 * @return The result of the stencil evaluation.
 */
template<typename S, typename T>
T ApplyStencilUpwind( const std::vector<double>& array, double const upwind_decision, double const cell_size ) {
   const S&& stencil = S();
   if( upwind_decision >= 0 ) {
      return stencil.template Apply<StencilProperty::UpwindLeft>( array, cell_size );
   } else {
      return stencil.template Apply<StencilProperty::UpwindRight>( array, cell_size );
   }
}

/**
 * @brief Applies a stencil on a field.
 * @tparam S The stencil to be used.
 * @tparam P The stencil property to be used (upwind, central...).
 * @tparam D
 * @tparam T The type of the values of the vector to which the stencil is applied.
 * @param buffer The buffer to which the stencil is applied.
 * @param i The index in x-direction.
 * @param j The index in y-direction.
 * @param k The index in z-direction.
 * @param cell_size The cell size.
 * @return The result of the stencil evaluation.
 */
template<typename S, StencilProperty P, Direction D, typename T>
T ApplyStencil( T const (&buffer)[CC::TCX()][CC::TCY()][CC::TCZ()], unsigned int const i, unsigned int const j, unsigned int const k, double const cell_size ) {
   const S&& stencil = S();
   return stencil.template Apply<P, D>( buffer, i, j, k, cell_size );
}

}

#endif // DIFFERENTIATION_UTILITIES_H
