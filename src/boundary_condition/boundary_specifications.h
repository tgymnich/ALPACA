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
#ifndef BOUNDARY_SPECIFICATIONS_H
#define BOUNDARY_SPECIFICATIONS_H

#include <type_traits>

/**
 * @brief Unique identifier to separate the (implemented) types of fluid boundary condition
 */
enum class FluidBoundaryType { Internal, ZeroGradient, Symmetry, FixedValue, Wall, Periodic };

/**
 * @brief Unique identifier to separate the (implemented) types of level-set boundary condition
 */
enum class LevelSetBoundaryType { Internal, ZeroGradient, Symmetry, Periodic };

// ATTENTION DO NOT GIVE THEM OTHER INDICIES. WILL BREAK THE ALGORITHM IN MULTIPLE PLACES!
/**
 * @brief Unique Identifier for the Boundary Location (Compass orientations plus Top/Bottom plus Diagonals)
 */
enum class BoundaryLocation : unsigned int {
   // normal Types used in 1D, for planes, and for DomainBoundaries
      East = 0, West = 1, North = 2, South = 3, Top = 4, Bottom = 5,
   // Diagonals in 2 Dimensions for Sticks in 3D and cubes in 2D
      BottomNorth, BottomSouth, TopNorth,  TopSouth, // x-Axis Sticks
      BottomEast,  BottomWest,  TopEast,   TopWest,  // y-Axis Sticks
      NorthEast,   NorthWest,   SouthEast, SouthWest,// z-Axis Sticks
   // Diagonals in 3 Dimensions for cubes in 3D
      EastNorthTop, EastNorthBottom, EastSouthTop, EastSouthBottom, //East-Side
      WestNorthTop, WestNorthBottom, WestSouthTop, WestSouthBottom, //West-Side
   //DO NOT EXTEND WITHOUT ADJUSTING THE OppositeDirection FUNCTION
};

/**
 * @brief Converts a BoundaryLocation identifier to a (C++11 standard compliant, i. e. positive) array index. "LTI = Location to Index"
 * @param l The location identifier.
 * @return Index to be used in Arrays.
 */
static inline constexpr std::underlying_type<BoundaryLocation>::type LTI(const BoundaryLocation l) { return static_cast<typename std::underlying_type<BoundaryLocation>::type>( l ); }

/**
 * @brief Gives the inverse direction. Convenience function.
 * @param location The direction that is to be inverted.
 * @return The inverse direction.
 */
constexpr BoundaryLocation OppositeDirection(const BoundaryLocation location) {
   switch (location) {
      case BoundaryLocation::East: //Planes
         return BoundaryLocation::West;
      case BoundaryLocation::West:
         return BoundaryLocation::East;
      case BoundaryLocation::North:
         return BoundaryLocation::South;
      case BoundaryLocation::South:
         return BoundaryLocation::North;
      case BoundaryLocation::Top:
         return BoundaryLocation::Bottom;
      case BoundaryLocation::Bottom:
         return BoundaryLocation::Top;

      case BoundaryLocation::BottomNorth: //Sticks
         return BoundaryLocation::TopSouth;
      case BoundaryLocation::BottomSouth:
         return BoundaryLocation::TopNorth;
      case BoundaryLocation::TopNorth:
         return BoundaryLocation::BottomSouth;
      case BoundaryLocation::TopSouth:
         return BoundaryLocation::BottomNorth;

      case BoundaryLocation::BottomEast:
         return BoundaryLocation::TopWest;
      case BoundaryLocation::BottomWest:
         return BoundaryLocation::TopEast;
      case BoundaryLocation::TopEast:
         return BoundaryLocation::BottomWest;
      case BoundaryLocation::TopWest:
         return BoundaryLocation::BottomEast;

      case BoundaryLocation::NorthEast:
         return BoundaryLocation::SouthWest;
      case BoundaryLocation::NorthWest:
         return BoundaryLocation::SouthEast;
      case BoundaryLocation::SouthEast:
         return BoundaryLocation::NorthWest;
      case BoundaryLocation::SouthWest:
         return BoundaryLocation::NorthEast;

      case BoundaryLocation::EastNorthTop: //Cubes
         return BoundaryLocation::WestSouthBottom;
      case BoundaryLocation::EastNorthBottom:
         return BoundaryLocation::WestSouthTop;
      case BoundaryLocation::EastSouthTop:
         return BoundaryLocation::WestNorthBottom;
      case BoundaryLocation::EastSouthBottom:
         return BoundaryLocation::WestNorthTop;
      case BoundaryLocation::WestNorthTop:
         return BoundaryLocation::EastSouthBottom;
      case BoundaryLocation::WestNorthBottom:
         return BoundaryLocation::EastSouthTop;
      case BoundaryLocation::WestSouthTop:
         return BoundaryLocation::EastNorthBottom;
      default: // Only remaning case: BoundaryLocation::WestSouthBottom:
         return BoundaryLocation::EastNorthTop;
   }
}

#endif // BOUNDARY_SPECIFICATIONS_H
