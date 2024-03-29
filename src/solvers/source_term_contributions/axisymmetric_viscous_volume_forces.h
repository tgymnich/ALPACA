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
#ifndef AXISYMMETRIC_VISCOUS_VOLUME_FORCES_H
#define AXISYMMETRIC_VISCOUS_VOLUME_FORCES_H

#include <vector>

#include "materials/material_manager.h"
#include "block.h"

/**
 * @brief This class calculates the viscous contribution to axisymmetric forces as described in \cite Meng2016b.
 */
class AxisymmetricViscousVolumeForces {

private:
   MaterialManager const& material_manager_;
   static constexpr unsigned int dim_ = 2; //only sane configuration; enables unit testing

   void ComputeVelocityGradient( double const (&u)[CC::TCX()][CC::TCY()][CC::TCZ()]
                               , double const (&v)[CC::TCX()][CC::TCY()][CC::TCZ()]
                               , double const (&w)[CC::TCX()][CC::TCY()][CC::TCZ()], double const cell_size
                               , double (&velocity_gradient)[CC::TCX()][CC::TCY()][1][dim_][dim_] ) const;

public:
   AxisymmetricViscousVolumeForces() = delete;
   explicit AxisymmetricViscousVolumeForces( MaterialManager const& material_manager );
   ~AxisymmetricViscousVolumeForces() = default;
   AxisymmetricViscousVolumeForces( AxisymmetricViscousVolumeForces const& ) = delete;
   AxisymmetricViscousVolumeForces( AxisymmetricViscousVolumeForces&& ) = delete;
   AxisymmetricViscousVolumeForces& operator=( AxisymmetricViscousVolumeForces const& ) = delete;
   AxisymmetricViscousVolumeForces& operator=( AxisymmetricViscousVolumeForces&& ) = delete;

   void ComputeForces( std::pair<const MaterialName, Block> const& mat_block, double (&axisymmetric_viscous_volume_forces)[FF::ANOE()][CC::ICX()][CC::ICY()][CC::ICZ()],
                       double const cell_size, double const x_block_coordinate ) const;
};

#endif //AXISYMMETRIC_VISCOUS_VOLUME_FORCES_H
