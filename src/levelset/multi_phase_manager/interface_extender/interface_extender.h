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
#ifndef INTERFACE_EXTENDER_H
#define INTERFACE_EXTENDER_H


#include "user_specifications/numerical_setup.h"
#include "levelset/geometry/geometry_calculator_marching_cubes.h"
#include "halo_manager.h"

/**
 * @brief Base class to perform two-way extrapolation for quantities defined at the interface.
 * @tparam DerivedInterfaceExtender Typename as template parameter due to CRTP.
 */
template<typename DerivedInterfaceExtender>
class InterfaceExtender {

   friend DerivedInterfaceExtender;

   HaloManager& halo_manager_;
   LogWriter& logger_;

   /**
    * @brief Default constructor of the InterfaceExtender class.
    * @param halo_manager Instance to a HaloManager which provides MPI-related methods.
    */
   explicit InterfaceExtender( HaloManager& halo_manager ) : halo_manager_( halo_manager ), logger_( LogWriter::Instance() ) {
      // Empty besides initializer list.
   }

public:
   InterfaceExtender() = delete;
   ~InterfaceExtender() = default;
   InterfaceExtender( InterfaceExtender const& ) = delete;
   InterfaceExtender& operator=( InterfaceExtender const& ) = delete;
   InterfaceExtender( InterfaceExtender&& ) = delete;
   InterfaceExtender& operator=( InterfaceExtender&& ) = delete;

   /**
    * @brief Performs an extension of interface quantities.
    * @param nodes The nodes for which scale separation should be done.
    */
   void ExtendInterfaceQuantities(std::vector<std::reference_wrapper<Node>> const& nodes) const {
      static_cast<DerivedInterfaceExtender const&>(*this).ExtendInterfaceQuantitiesImplementation(nodes);
   }
};


#endif //INTERFACE_EXTENDER_H
