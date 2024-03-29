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
#ifndef VISCOUS_FLUXES_H
#define VISCOUS_FLUXES_H

#include <vector>

#include "materials/material_manager.h"
#include "block.h"

/**
 * @brief This class calculates the viscous source terms and adds them to a flux buffer.
 */
class ViscousFluxes {

private:
   MaterialManager const& material_manager_;

   void ComputeVelocityGradient( double const (&u)[CC::TCX()][CC::TCY()][CC::TCZ()], double const (&v)[CC::TCX()][CC::TCY()][CC::TCZ()], double const (&w)[CC::TCX()][CC::TCY()][CC::TCZ()],
      double const cell_size, double (&velocity_gradient)[CC::TCX()][CC::TCY()][CC::TCZ()][DTI(CC::DIM())][DTI(CC::DIM())] ) const;

   void ComputeVelocityGradientAtCellFaces( double const (&velocity_gradient)[CC::TCX()][CC::TCY()][CC::TCZ()][DTI(CC::DIM())][DTI(CC::DIM())], double const cell_size,
      double const (&u)[CC::TCX()][CC::TCY()][CC::TCZ()], double const (&v)[CC::TCX()][CC::TCY()][CC::TCZ()], double const (&w)[CC::TCX()][CC::TCY()][CC::TCZ()],
      double (&velocity_gradient_at_cell_faces)[CC::ICX() + 1][CC::ICY() + 1][CC::ICZ() + 1][DTI(CC::DIM())][DTI(CC::DIM())][DTI(CC::DIM())] ) const;

   void ComputeVelocityAtCellFaces( double const (&u)[CC::TCX()][CC::TCY()][CC::TCZ()], double const (&v)[CC::TCX()][CC::TCY()][CC::TCZ()], double const (&w)[CC::TCX()][CC::TCY()][CC::TCZ()],
      double const cell_size, double (&velocity_at_cell_faces)[CC::ICX()+1][CC::ICY()+1][CC::ICZ()+1][DTI(CC::DIM())][DTI(CC::DIM())]) const;

   void ComputeTau( double const (&velocity_gradient_at_cell_faces)[CC::ICX()+1][CC::ICY()+1][CC::ICZ()+1][DTI(CC::DIM())][DTI(CC::DIM())][DTI(CC::DIM())],
      std::vector<double> const viscosity, double (&tau)[CC::ICX()+1][CC::ICY()+1][CC::ICZ()+1][DTI(CC::DIM())][DTI(CC::DIM())]) const;

public:
   ViscousFluxes() = delete;
   explicit ViscousFluxes( MaterialManager const& material_manager );
   ~ViscousFluxes() = default;
   ViscousFluxes( ViscousFluxes const& ) = delete;
   ViscousFluxes( ViscousFluxes&& ) = delete;
   ViscousFluxes& operator=( ViscousFluxes const& ) = delete;
   ViscousFluxes& operator=( ViscousFluxes&& ) = delete;

   void ComputeFluxes( std::pair<MaterialName const, Block> const& mat_block,
      double (&dissipative_flux_x)[FF::ANOE()][CC::ICX() + 1][CC::ICY() + 1][CC::ICZ() + 1],
      double (&dissipative_flux_y)[FF::ANOE()][CC::ICX() + 1][CC::ICY() + 1][CC::ICZ() + 1],
      double (&dissipative_flux_z)[FF::ANOE()][CC::ICX() + 1][CC::ICY() + 1][CC::ICZ() + 1], double const cell_size ) const;
};

#endif //VISCOUS_FLUXES_H
