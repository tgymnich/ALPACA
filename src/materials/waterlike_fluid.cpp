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
#include "waterlike_fluid.h"
#include "mathematical_functions.h"
#include <cmath>

/**
 * @brief Constructs a waterlike fluid object with material parameters given as input.
 * @param gamma .
 * @param A .
 * @param B .
 * @param C .
 * @param rho0 .
 * @param specific_gas_constant .
 * @param mu_shear Shear viscosity.
 * @param mu_bulk Bulk viscosity.
 */
WaterlikeFluid::WaterlikeFluid( double const gamma, double const A, double const B, double const rho0, double const mu_shear, double const mu_bulk ) :
   gamma_(gamma),
   A_(A),
   B_(B),
   rho0_(rho0),
   mu_shear_(mu_shear),
   mu_bulk_(mu_bulk) {
   /* Empty besides initializer list*/
}

/**
 * @brief Computes pressure from inputs as A - B + B * (rho / rho0)^gamma.
 * @param density .
 * @param momentum_x (not applied for Tait) .
 * @param momentum_y (not applied for Tait) .
 * @param momentum_z (not applied for Tait) .
 * @param energy (not applied for Tait) .
 * @return Pressure according to Tait's equation of state.
 */
double WaterlikeFluid::DoGetPressure( double const density, double const, double const, double const, double const ) const {
   return A_ - B_ + B_ * std::pow( density / rho0_, gamma_ );
}

/**
 * @brief Gives enthalpy for the given inputs. (Not available for classic Tait).
 * @param density (not applied for Tait) .
 * @param momentum_x (not applied for Tait) .
 * @param momentum_y (not applied for Tait) .
 * @param momentum_z (not applied for Tait) .
 * @param energy (not applied for Tait) .
 * @return Zero. This is according to Tait's equation of state correct.
 */
double WaterlikeFluid::DoGetEnthalpy(double const, double const, double const, double const, double const ) const {
   return 0.0;
}

/**
 * @brief Computes energy according to 1/(gamma-1) * (p + B - A) + B - A + 1/2 * rho *|v|^2.
 * @param density .
 * @param momentum_x .
 * @param momentum_y .
 * @param momentum_z .
 * @param pressure .
 * @return Energy according to Tait's equation of state.
 */
double WaterlikeFluid::DoGetEnergy( double const density, double const momentum_x, double const momentum_y, double const momentum_z, double const pressure ) const {
   return 1.0 / ( 1.0 - gamma_ ) * ( pressure + B_ - A_ ) + B_ - A_ + ( 0.5 * DimensionAwareConsistencyManagedSum( momentum_x * momentum_x,
                                                                                                                   momentum_y * momentum_y,
                                                                                                                   momentum_z * momentum_z ) / density);
}

/**
 * @brief Computes Gruneisen coefficent as (gamma-1) for stiffened-gas equation of state.
 * @return Gruneisen coefficent according to Tait's equation of state.
 */
double WaterlikeFluid::DoGetGruneisen() const {
   return 0.0;
}

/**
 * @brief Computes psi from inputs as gamma * (p + B - A) / rho.
 * @param pressure .
 * @param one_density .
 * @return Psi according to Tait's equation of state.
 */
double WaterlikeFluid::DoGetPsi( double const pressure, double const one_density ) const {
   return gamma_ * ( pressure + B_ - A_ ) * one_density ;
}

/**
 * @brief Computes speed of sound from inputs as sqrt(gamma * (p + B - A) / rho).
 * @param density .
 * @param pressure .
 * @return Speed of sound according to Tait's equation of state.
 */
double WaterlikeFluid::DoGetSpeedOfSound( double const density, double const pressure ) const {
   return std::sqrt(gamma_ * (pressure + B_ - A_) / density);
}

/**
 * @brief Gives the viscosity of the material.
 * @return first: shear viscosity; second: bulk viscosity
 */
std::vector<double> WaterlikeFluid::DoGetViscosity() const {
   return {mu_shear_, mu_bulk_};
}
