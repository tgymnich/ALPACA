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
#ifndef TWO_PHASE_INTERFACE_EXTENDER_H
#define TWO_PHASE_INTERFACE_EXTENDER_H


#include "interface_extender.h"

/**
 * @brief Implements the two-way extrapolation of interface quantities for two phases.
 */
class TwoPhaseInterfaceExtender : public InterfaceExtender<TwoPhaseInterfaceExtender> {

   friend InterfaceExtender;

private:

   /**
    * The number of quantities that are necessary to track convergence. Those are the maximum values of the quantities to extend and the residuum.
    */
   static constexpr unsigned int number_of_convergence_tracking_quantities_ = FF::NOIQTE() + 1;

   double const epsilon_ = std::numeric_limits<double>::epsilon();

   void IterativeInterfaceExtension(Node& node, std::vector<double>& convergence_tracking_quantities) const;
   void DetermineMaximumInterfaceQuantityValues(Node const& node, std::vector<double>& convergence_tracking_quantities) const;

   void ExtendInterfaceQuantitiesImplementation(std::vector<std::reference_wrapper<Node>> const& nodes) const;

public:
   TwoPhaseInterfaceExtender() = delete;
   explicit TwoPhaseInterfaceExtender( HaloManager& halo_manager );
   ~TwoPhaseInterfaceExtender() = default;
   TwoPhaseInterfaceExtender( TwoPhaseInterfaceExtender const& ) = delete;
   TwoPhaseInterfaceExtender& operator=( TwoPhaseInterfaceExtender const& ) = delete;
   TwoPhaseInterfaceExtender( TwoPhaseInterfaceExtender&& ) = delete;
   TwoPhaseInterfaceExtender& operator=( TwoPhaseInterfaceExtender&& ) = delete;
};


#endif //TWO_PHASE_INTERFACE_EXTENDER_H
