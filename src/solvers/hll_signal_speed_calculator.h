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
#ifndef HLL_SIGNAL_SPEED_CALCULATOR_H
#define HLL_SIGNAL_SPEED_CALCULATOR_H

#include <cmath>
#include "user_specifications/riemann_solver_settings.h"

/**
 * @brief Computes the signal speed estimates for Hll-type solvers according to Einfeldt \cite Einfeldt88, Davis \cite Davis88, Toro \cite Toro94, Arithmetic-averaged \cite Coralic14
 * @param density_left .
 * @param density_right .
 * @param velocity_left .
 * @param velocity_right .
 * @param pressure_left .
 * @param pressure_right .
 * @param speed_of_sound_left .
 * @param speed_of_sound_right .
 * @param gamma .
 */
inline std::pair<double, double> CalculateSignalSpeed ( double const density_left, double const density_right,
                                                        double const velocity_left, double const velocity_right,
                                                        double const pressure_left, double const pressure_right,
                                                        double const speed_of_sound_left, double const speed_of_sound_right,
                                                        double const gamma ) {

   if constexpr( HllSolverSettings::signal_speed_selection == SignalSpeed::Arithmetic ) {
      double const speed_of_sound_face_average = 0.5 * ( speed_of_sound_left + speed_of_sound_right );
      double const velocity_face_average       = 0.5 * ( velocity_left + velocity_right );

      double const wave_speed_left_simple   = std::min( velocity_face_average - speed_of_sound_face_average , velocity_left  - speed_of_sound_left );
      double const wave_speed_right_simple  = std::max( velocity_face_average + speed_of_sound_face_average , velocity_right + speed_of_sound_right );

      return std::make_pair( wave_speed_left_simple, wave_speed_right_simple );
    }

   if constexpr( HllSolverSettings::signal_speed_selection == SignalSpeed::Davis ) {
      double const wave_speed_left_simple   = std::min( velocity_left - speed_of_sound_left , velocity_right - speed_of_sound_right );
      double const wave_speed_right_simple  = std::max( velocity_left + speed_of_sound_left , velocity_right + speed_of_sound_right );

      return std::make_pair( wave_speed_left_simple, wave_speed_right_simple );
   }

   if constexpr( HllSolverSettings::signal_speed_selection == SignalSpeed::Einfeldt ) {
      double const sqrt_density_left   = std::sqrt( density_left );
      double const sqrt_density_right  = std::sqrt( density_right );
      double const rho_div = 1.0 / ( sqrt_density_left + sqrt_density_right );

      double const u_roe_avg = ( ( velocity_left * sqrt_density_left ) + ( velocity_right * sqrt_density_right ) ) * rho_div;
      double const speed_of_sound_einfeldt = std::sqrt( ( ( speed_of_sound_left * speed_of_sound_left * sqrt_density_left ) + ( speed_of_sound_right * speed_of_sound_right * sqrt_density_right ) ) * rho_div
                                            + ( 0.5 * sqrt_density_left * sqrt_density_right * rho_div * rho_div ) * ( velocity_right-velocity_left ) * ( velocity_right - velocity_left ) );
      double const wave_speed_left_simple   = std::min( u_roe_avg - speed_of_sound_einfeldt, velocity_left - speed_of_sound_left );
      double const wave_speed_right_simple  = std::max( u_roe_avg + speed_of_sound_einfeldt, velocity_right + speed_of_sound_right );

      return std::make_pair( wave_speed_left_simple, wave_speed_right_simple );
   }

   if constexpr( HllSolverSettings::signal_speed_selection == SignalSpeed::Toro ) {
      double const z = ( gamma - 1.0 ) / ( 2.0 * gamma );
      double const p_star = std::pow ( ( ( speed_of_sound_left + speed_of_sound_right ) - 0.5 * ( velocity_right - velocity_left ) * ( gamma - 1.0 ) )
                                       / ( speed_of_sound_left / std::pow( pressure_left, z ) + speed_of_sound_right / std::pow( pressure_right,z ) ), ( 1.0 / z ) );
      double q_l = 1.0;
      double q_r = 1.0;

      if( p_star > pressure_left ) {
         q_l = std::sqrt( 1.0 + ( gamma + 1.0 ) / ( 2.0 * gamma ) * ( p_star / pressure_left - 1.0 ) );
      }

      if( p_star > pressure_right ) {
         q_r = std::sqrt( 1.0 + ( gamma + 1.0 ) / ( 2.0 * gamma ) * ( p_star / pressure_right - 1.0 ) );
      }

      double const wave_speed_left_simple   = velocity_left  - speed_of_sound_left  * q_l;
      double const wave_speed_right_simple  = velocity_right + speed_of_sound_right * q_r;

      return std::make_pair( wave_speed_left_simple, wave_speed_right_simple );
   }
}

#endif // HLL_SIGNAL_SPEED_CALCULATOR_H