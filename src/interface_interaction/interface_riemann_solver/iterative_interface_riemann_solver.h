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
#ifndef ITERATIVE_INTERFACE_RIEMANN_SOLVER_H
#define ITERATIVE_INTERFACE_RIEMANN_SOLVER_H

#include <cmath>
#include <limits>

#include "interface_riemann_solver.h"
#include "user_specifications/two_phase_constants.h"

/**
 * Constants used during the iterative solution of the Riemann problem.
 */
namespace IterationConstants {
/**
 * @brief      Computes the constant A.
 *
 * @param[in]  gamma    The gamma of the material.
 * @param[in]  density  The density.
 *
 * @return     The constant A.
 */
constexpr double A( double const gamma, double const density ) {
   return 2.0 / ( ( gamma + 1.0 ) * density );
}

/**
 * @brief      Computes the constant B.
 *
 * @param[in]  gamma  The gamma of the material.
 * @param[in]  p_ref  The left/right pressure of the Riemann problem.
 *
 * @return     The result.
 */
constexpr double B( double const gamma, double const p_ref ) {
   return p_ref * ( ( gamma - 1.0 ) / ( gamma + 1.0 ) );
}

/**
 * @brief      Computes the constant C.
 *
 * @param[in]  gamma           The gamma of the material.
 * @param[in]  speed_of_sound  The speed of sound.
 *
 * @return     The result.
 */
constexpr double C( double const gamma, double const speed_of_sound ) {
   return 2.0 * speed_of_sound / ( gamma - 1.0 );
}

/**
 * @brief      Computes the constant D.
 *
 * @param[in]  gamma  The ratio of specific heats of the material.
 *
 * @return     The result.
 */
constexpr double D( double const gamma ) {
   return ( gamma - 1.0 ) / ( 2.0 * gamma );
}
}

namespace IterationUtilities {

/**
 * @brief      Computes a value necessary for the iterative solution procedure.
 *
 * @param[in]  p           The pressure.
 * @param[in]  material_B  The material constant B.
 *
 * @return     The result.
 */
constexpr double MaterialPressureFunction( double const p, double const material_B ) {
   return p + material_B;
}

/**
 * @brief      Returns the current residuum of the iterative solution procedure.
 *
 * @param[in]  old_value  The old value.
 * @param[in]  new_value  The new value.
 *
 * @return     The tolerance.
 */
inline double GetTolerance( double const old_value, double const new_value ) {
   return std::abs( new_value - old_value ) / std::max( std::abs( old_value ), std::numeric_limits<double>::epsilon() );
}
}

namespace {

/**
 * @brief Evaluates the function whose zero point is calculated iteratively.
 * @param function_left The evaluated relation for the left state.
 * @param function_right The evaluated relation for the right state.
 * @param velocity_difference The velocity difference between the left and right state.
 * @return The function value.
 */
constexpr double RootFunction ( double const function_left, double const function_right, double const velocity_difference ) {
   return function_left + function_right + velocity_difference;
}

/**
 * @brief Evaluates the derivative of the function whose zero point is calculated iteratively.
 * @param derivative_left The derivative of the relation of the left side.
 * @param derivative_right The derivative of the relation of the right side.
 * @return The derivative.
 */
constexpr double DerivativeOfRootFunction ( double const derivative_left, double const derivative_right ) {
   return derivative_left + derivative_right;
}
}


/**
 * @brief Implements an iterative solution procedure for a Riemann problem. In ALPACA, it is used to solve the interface Riemann problem.
 * @tparam DerivedIterativeInterfaceRiemannSolver Template parameter necessary for the CRTP.
 */
template<typename DerivedIterativeInterfaceRiemannSolver>
class IterativeInterfaceRiemannSolver : public InterfaceRiemannSolver<DerivedIterativeInterfaceRiemannSolver> {

   friend InterfaceRiemannSolver<DerivedIterativeInterfaceRiemannSolver>;
   friend DerivedIterativeInterfaceRiemannSolver;

   using InterfaceRiemannSolver<DerivedIterativeInterfaceRiemannSolver>::material_manager_;

   /**
    * @brief Constructor for the IterativeInterfaceRiemannSolver.
    * @param material_manager See base class.
    */
   explicit IterativeInterfaceRiemannSolver( MaterialManager const& material_manager ) : InterfaceRiemannSolver<DerivedIterativeInterfaceRiemannSolver>( material_manager ) {
      // Empty besides call of base class constructor.
   }

   /**
    * @brief      Calculates the shock/rarefaction relation and its derivative.
    *
    * @param[in]  initial_root           The interface pressure of the current iteration.
    * @param[in]  p                      The left/right pressure defining the interface Riemann problem.
    * @param[in]  pressure_function      A constant dependent on the EoS.
    * @param[in]  one_pressure_function  The inverse of the pressure function.
    * @param[in]  pressure_constant      The pressure constant B.
    * @param[in]  A                      A constant. See IterationConstants::A() for details.
    * @param[in]  B                      A constant. See IterationConstants::B() for details.
    * @param[in]  C                      A constant. See IterationConstants::C() for details.
    * @param[in]  D                      A constant. See IterationConstants::D() for details.
    *
    * @return     The shock/rarefaction relation and its derivative.
    */
   std::array<double, 2> ObtainFunctionAndDerivative( double const initial_root, double const p,
      double const pressure_function, double const one_pressure_function, double const pressure_constant, double const A, double const B, double const C, double const D ) const {
      return static_cast<DerivedIterativeInterfaceRiemannSolver const&>( *this ).ObtainFunctionAndDerivativeImplementation( initial_root, p,
         pressure_function, one_pressure_function, pressure_constant, A, B, C, D );
   }

   /**
    * @brief Solves the Riemann problem at the interface iteratively.
    * @param rho_left Density of the left fluid.
    * @param p_left Pressure of the left fluid.
    * @param velocity_normal_left Velocity normal to the interface of the left fluid.
    * @param material_left Material of the left fluid.
    * @param rho_right Density of the right fluid.
    * @param p_right Pressure of the right fluid.
    * @param velocity_normal_right Velocity normal to the interface of the right fluid.
    * @param material_right Material of the right fluid.
    * @param delta_p Pressure jump due to capillarity.
    * @return An array that contains following information in the given order: interface_velocity, interface_pressure_positive, interface_pressure_negative.
    */
   std::array<double, 3> SolveInterfaceRiemannProblemImplementation( double const rho_left, double const p_left, double const velocity_normal_left, MaterialName const material_left,
      double const rho_right, double const p_right, double const velocity_normal_right, MaterialName const material_right,
      double const delta_p ) const {
      /**
       * For the iteration procedure several constants can be computed in advance. Those constants can be found in \cite Toro2009 chapter 4.3.
       */
      double const speed_of_sound_left = material_manager_.GetSpeedOfSound( material_left, rho_left, p_left );
      double const speed_of_sound_right = material_manager_.GetSpeedOfSound( material_right, rho_right, p_right );

      double const impedance_left = rho_left * speed_of_sound_left;
      double const impedance_right = rho_right * speed_of_sound_right;

      double const inverse_impedance_sum = 1.0 / std::max( ( impedance_left + impedance_right ), std::numeric_limits<double>::epsilon() );

      double const gamma_left = material_manager_.GetGamma( material_left );
      double const gamma_right = material_manager_.GetGamma( material_right );

      double const pressure_constant_left = material_manager_.GetB( material_left );
      double const pressure_constant_right = material_manager_.GetB( material_right );

      double const pressure_function_left = IterationUtilities::MaterialPressureFunction( p_left, pressure_constant_left );
      double const one_pressure_function_left = 1.0 / std::max( pressure_function_left, std::numeric_limits<double>::epsilon() );
      double const pressure_function_right = IterationUtilities::MaterialPressureFunction( p_right, pressure_constant_right );
      double const one_pressure_function_right = 1.0 / std::max( pressure_function_right, std::numeric_limits<double>::epsilon() );

      double const velocity_difference = velocity_normal_right - velocity_normal_left;
      double const velocity_sum = velocity_normal_right + velocity_normal_left;

      double const A_left = IterationConstants::A( gamma_left, rho_left );
      double const B_left = IterationConstants::B( gamma_left, pressure_function_left );
      double const C_left = IterationConstants::C( gamma_left, speed_of_sound_left );
      double const D_left = IterationConstants::D( gamma_left );

      double const A_right = IterationConstants::A( gamma_right, rho_right );
      double const B_right = IterationConstants::B( gamma_right, pressure_function_right );
      double const C_right = IterationConstants::C( gamma_right, speed_of_sound_right );
      double const D_right = IterationConstants::D( gamma_right );

      double initial_root = ( impedance_left  *  p_right           +
         impedance_right * ( p_left - delta_p ) +
         impedance_left  *  impedance_right   * ( velocity_normal_left - velocity_normal_right ) ) * inverse_impedance_sum;

      bool converged = false;
      double interface_velocity = 0.0;
      double interface_pressure_positive = 0.0;
      double interface_pressure_negative = 0.0;

      /**
       * Iterative procedure to compute the interface states.
       */
      for( unsigned it = 0; it < IterativeInterfaceRiemannSolverConstants::MaximumNumberOfIterations; ++it ) {

         const std::array<double, 2> relations_left = ObtainFunctionAndDerivative(initial_root, p_left, pressure_function_left, one_pressure_function_left, pressure_constant_left, A_left, B_left, C_left, D_left);
         const std::array<double, 2> relations_right = ObtainFunctionAndDerivative(initial_root, p_right, pressure_function_right, one_pressure_function_right, pressure_constant_right, A_right, B_right, C_right, D_right);

         double const derivative_of_root_function = DerivativeOfRootFunction( relations_left[1], relations_right[1] );

         if( derivative_of_root_function != 0.0 ) {

            double const next_root = initial_root - RootFunction( relations_left[0], relations_right[0], velocity_difference ) / std::max( derivative_of_root_function, std::numeric_limits<double>::epsilon() );

            if( IterationUtilities::GetTolerance( initial_root, next_root ) < IterativeInterfaceRiemannSolverConstants::MaximumResiduum ) {
               converged = true;

               interface_pressure_positive = next_root;
               interface_pressure_negative = interface_pressure_positive - delta_p; //interface pressure right
               interface_velocity = 0.5 * velocity_sum + 0.5 * (relations_right[0] - relations_left[0]);

               break;
            }
            initial_root = next_root;
         } else {
            break;
         }
      }

      /**
       * In case the iterative procedure did not converge within the specified maximum number of iterations, we take the linearized solution.
       */
      if( !converged ) {
         interface_velocity = ( impedance_left *  velocity_normal_left  +
            impedance_right * velocity_normal_right +
            p_left - p_right - delta_p ) * inverse_impedance_sum;

         interface_pressure_positive = ( impedance_left  *  p_right           +
            impedance_right * ( p_left - delta_p ) +
            impedance_left  *  impedance_right   * ( velocity_normal_left - velocity_normal_right ) ) * inverse_impedance_sum;

         interface_pressure_negative = ( impedance_left  * ( p_right + delta_p ) +
            impedance_right *  p_left             +
            impedance_left  *  impedance_right    * ( velocity_normal_left - velocity_normal_right ) ) * inverse_impedance_sum;
      }

      return {interface_velocity, interface_pressure_positive, interface_pressure_negative};
   }

public:
   IterativeInterfaceRiemannSolver() = delete;
   ~IterativeInterfaceRiemannSolver() = default;
   IterativeInterfaceRiemannSolver( IterativeInterfaceRiemannSolver const& ) = delete;
   IterativeInterfaceRiemannSolver& operator=( IterativeInterfaceRiemannSolver const& ) = delete;
   IterativeInterfaceRiemannSolver( IterativeInterfaceRiemannSolver&& ) = delete;
   IterativeInterfaceRiemannSolver& operator=( IterativeInterfaceRiemannSolver&& ) = delete;
};

/**
 * Implements the shock relations as described in \cite Toro2009 chapter 4.3
 */
namespace ShockRelations {

/**
 * @brief      Computes the shock relation.
 *
 * @param[in]  p      The pressure for which the function should be evaluated.
 * @param[in]  p_ref  The left/right pressure defining the Riemann problem.
 * @param[in]  A      The constant A.
 * @param[in]  B      The constant B.
 *
 * @return     The resulting function value.
 */
inline double Function( double p, double const p_ref, double const A, double const B ) {
   p = std::max( p, std::numeric_limits<double>::epsilon() );
   return ( p - p_ref ) * std::sqrt( A / ( p + B ) );
}

/**
 * @brief      Computes the derivative of the shock relation.
 *
 * @param[in]  p      The pressure for which the function should be evaluated.
 * @param[in]  p_ref  The left/right pressure defining the Riemann problem.
 * @param[in]  A      The constant A.
 * @param[in]  B      The constant B.
 *
 * @return     The derivative of the function evaluation.
 */
inline double Derivative( double p, double const p_ref, double const A, double const B ) {
   p = std::max( p, std::numeric_limits<double>::epsilon() );
   return std::sqrt( A / ( B + p ) ) * ( 1.0 - ( p - p_ref ) / ( 2.0 * ( B + p ) ) );
}
}

/**
 * Implements the rarefaction relations as described in \cite Toro2009 chapter 4.3
 */
namespace RarefactionRelations {

/**
 * @brief      Computes the rarefaction relation.
 *
 * @param[in]  p          The pressure for which the function should be evaluated.
 * @param[in]  one_p_ref  The left/right pressure of the Riemann problem.
 * @param[in]  C          The constant C.
 * @param[in]  D          The constant D.
 *
 * @return     The result.
 */
inline double Function( double p, double const one_p_ref, double const C, double const D ) {
   p = std::max( p, std::numeric_limits<double>::epsilon() );
   return C * ( std::pow( p * one_p_ref, D ) - 1.0 );
}

/**
 * @brief      The derivative of the rarefaction relation.
 *
 * @param[in]  p          The pressure for which the function should be evaluated.
 * @param[in]  one_p_ref  The left/right pressure of the Riemann problem.
 * @param[in]  C          The constant C.
 * @param[in]  D          The constant D.
 *
 * @return     The result.
 */

inline double Derivative( double p, double const one_p_ref, double const C, double const D ) {
   p = std::max( p, std::numeric_limits<double>::epsilon() );
   return C * D * std::pow( one_p_ref, D ) * std::pow( p, ( D - 1.0 ) );
}
}


#endif //ITERATIVE_INTERFACE_RIEMANN_SOLVER_H
