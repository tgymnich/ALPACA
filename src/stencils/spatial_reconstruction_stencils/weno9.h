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
#ifndef WENO9_H
#define WENO9_H

#include "stencils/stencil.h"

/**
 * @brief Discretization of the SpatialReconstructionStencil class to compute WENO9 fluxes according to \cite Balsara2000.
 */
class WENO9 : public Stencil<WENO9> {

   friend Stencil;

   static constexpr StencilType stencil_type_ = StencilType::Reconstruction;

   // Coefficients for WENO9 scheme
   static constexpr double coef_smoothness_0_01_ =  22658.0;
   static constexpr double coef_smoothness_0_02_ = -208501.0;
   static constexpr double coef_smoothness_0_03_ =  364863.0;
   static constexpr double coef_smoothness_0_04_ = -288007.0;
   static constexpr double coef_smoothness_0_05_ =  86329.0;
   static constexpr double coef_smoothness_0_06_ =  482963.0;
   static constexpr double coef_smoothness_0_07_ = -1704396.0;
   static constexpr double coef_smoothness_0_08_ =  1358458.0;
   static constexpr double coef_smoothness_0_09_ = -411487.0;
   static constexpr double coef_smoothness_0_10_ =  1521393.0;
   static constexpr double coef_smoothness_0_11_ = -2462076.0;
   static constexpr double coef_smoothness_0_12_ =  758823.0;
   static constexpr double coef_smoothness_0_13_ =  1020563.0;
   static constexpr double coef_smoothness_0_14_ = -649501.0;
   static constexpr double coef_smoothness_0_15_ =  107918.0;

   static constexpr double coef_smoothness_1_01_ =  6908.0;
   static constexpr double coef_smoothness_1_02_ = -60871.0;
   static constexpr double coef_smoothness_1_03_ =  99213.0;
   static constexpr double coef_smoothness_1_04_ = -70237.0;
   static constexpr double coef_smoothness_1_05_ =  18079.0;
   static constexpr double coef_smoothness_1_06_ =  138563.0;
   static constexpr double coef_smoothness_1_07_ = -464976.0;
   static constexpr double coef_smoothness_1_08_ =  337018.0;
   static constexpr double coef_smoothness_1_09_ = -88297.0;
   static constexpr double coef_smoothness_1_10_ =  406293.0;
   static constexpr double coef_smoothness_1_11_ = -611976.0;
   static constexpr double coef_smoothness_1_12_ =  165153.0;
   static constexpr double coef_smoothness_1_13_ =  242723.0;
   static constexpr double coef_smoothness_1_14_ = -140251.0;
   static constexpr double coef_smoothness_1_15_ =  22658.0;

   static constexpr double coef_smoothness_2_01_ =  6908.0;
   static constexpr double coef_smoothness_2_02_ = -51001.0;
   static constexpr double coef_smoothness_2_03_ =  67923.0;
   static constexpr double coef_smoothness_2_04_ = -38947.0;
   static constexpr double coef_smoothness_2_05_ =  8209.0;
   static constexpr double coef_smoothness_2_06_ =  104963.0;
   static constexpr double coef_smoothness_2_07_ = -299076.0;
   static constexpr double coef_smoothness_2_08_ =  179098.0;
   static constexpr double coef_smoothness_2_09_ = -38947.0;
   static constexpr double coef_smoothness_2_10_ =  231153.0;
   static constexpr double coef_smoothness_2_11_ = -299076.0;
   static constexpr double coef_smoothness_2_12_ =  67923.0;
   static constexpr double coef_smoothness_2_13_ =  104963.0;
   static constexpr double coef_smoothness_2_14_ = -51001.0;
   static constexpr double coef_smoothness_2_15_ =  6908.0;

   static constexpr double coef_smoothness_3_01_ =  22658.0;
   static constexpr double coef_smoothness_3_02_ = -140251.0;
   static constexpr double coef_smoothness_3_03_ =  165153.0;
   static constexpr double coef_smoothness_3_04_ = -88297.0;
   static constexpr double coef_smoothness_3_05_ =  18079.0;
   static constexpr double coef_smoothness_3_06_ =  242723.0;
   static constexpr double coef_smoothness_3_07_ = -611976.0;
   static constexpr double coef_smoothness_3_08_ =  337018.0;
   static constexpr double coef_smoothness_3_09_ = -70237.0;
   static constexpr double coef_smoothness_3_10_ =  406293.0;
   static constexpr double coef_smoothness_3_11_ = -464976.0;
   static constexpr double coef_smoothness_3_12_ =  99213.0;
   static constexpr double coef_smoothness_3_13_ =  138563.0;
   static constexpr double coef_smoothness_3_14_ = -60871.0;
   static constexpr double coef_smoothness_3_15_ =  6908.0;

   static constexpr double coef_smoothness_4_01_ =  107918.0;
   static constexpr double coef_smoothness_4_02_ = -649501.0;
   static constexpr double coef_smoothness_4_03_ =  758823.0;
   static constexpr double coef_smoothness_4_04_ = -411487.0;
   static constexpr double coef_smoothness_4_05_ =  86329.0;
   static constexpr double coef_smoothness_4_06_ =  1020563.0;
   static constexpr double coef_smoothness_4_07_ = -2462076.0;
   static constexpr double coef_smoothness_4_08_ =  1358458.0;
   static constexpr double coef_smoothness_4_09_ = -288007.0;
   static constexpr double coef_smoothness_4_10_ =  1521393.0;
   static constexpr double coef_smoothness_4_11_ = -1704396.0;
   static constexpr double coef_smoothness_4_12_ =  364863.0;
   static constexpr double coef_smoothness_4_13_ =  482963.0;
   static constexpr double coef_smoothness_4_14_ = -208501.0;
   static constexpr double coef_smoothness_4_15_ =  22658.0;

   static constexpr double coef_weights_1_ =  1.0/126.0;
   static constexpr double coef_weights_2_ = 10.0/63.0;
   static constexpr double coef_weights_3_ = 10.0/21.0;
   static constexpr double coef_weights_4_ = 20.0/63.0;
   static constexpr double coef_weights_5_ =  5.0/126.0;

   static constexpr double coef_stencils_1_  =  12.0/60.0;
   static constexpr double coef_stencils_2_  = -63.0/60.0;
   static constexpr double coef_stencils_3_  = 137.0/60.0;
   static constexpr double coef_stencils_4_  =-163.0/60.0;
   static constexpr double coef_stencils_5_  = 137.0/60.0;

   static constexpr double coef_stencils_6_  =  -3.0/60.0;
   static constexpr double coef_stencils_7_  =  17.0/60.0;
   static constexpr double coef_stencils_8_  = -43.0/60.0;
   static constexpr double coef_stencils_9_  =  77.0/60.0;
   static constexpr double coef_stencils_10_ =  12.0/60.0;

   static constexpr double coef_stencils_11_ =   2.0/60.0;
   static constexpr double coef_stencils_12_ = -13.0/60.0;
   static constexpr double coef_stencils_13_ =  47.0/60.0;
   static constexpr double coef_stencils_14_ =  27.0/60.0;
   static constexpr double coef_stencils_15_ =  -3.0/60.0;

   static constexpr double coef_stencils_16_ =  -3.0/60.0;
   static constexpr double coef_stencils_17_ =  27.0/60.0;
   static constexpr double coef_stencils_18_ =  47.0/60.0;
   static constexpr double coef_stencils_19_ = -13.0/60.0;
   static constexpr double coef_stencils_20_ =   2.0/60.0;

   static constexpr double coef_stencils_21_ =  12.0/60.0;
   static constexpr double coef_stencils_22_ =  77.0/60.0;
   static constexpr double coef_stencils_23_ = -43.0/60.0;
   static constexpr double coef_stencils_24_ =  17.0/60.0;
   static constexpr double coef_stencils_25_ =  -3.0/60.0;

   // Small values to avoid division by 0, but also to adjust dissipation
   static constexpr double epsilon_weno9_ = 1.0e-10;

   static constexpr unsigned int stencil_size_            = 10;
   static constexpr unsigned int downstream_stencil_size_ = 4;

   double ApplyImplementation( std::vector<double> const& array, int const stencil_offset, int const stencil_sign, double const cell_size ) const;

public:
   explicit WENO9() = default;
   ~WENO9() = default;
   WENO9( WENO9 const& ) = delete;
   WENO9& operator=( WENO9 const& ) = delete;
   WENO9( WENO9&& ) = delete;
   WENO9& operator=( WENO9&& ) = delete;
};

#endif // STENCIL_WENO9_H
