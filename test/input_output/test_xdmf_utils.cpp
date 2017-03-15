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
#include <catch.hpp>
#include "input_output/xdmf_utilities.h"

using namespace XdmfUtilities;

SCENARIO( "DataItems can be correctly populated", "[1rank]" ) {
   GIVEN( "An hdf5 filename (test.h5) + 42 global cells + a dataset named ds" ) {
      std::string const filename = "test.h5";
      constexpr unsigned int number_of_cells = 42;
      std::string const dataset = "ds";
      WHEN( "Created with defaults" ) {
         REQUIRE( DataItemString( filename, dataset, number_of_cells ) == "<DataItem Format=\"HDF\" NumberType=\"Float\" Precision=\"8\" Dimensions=\"42\">test.h5:/ds</DataItem>" );
      }
      WHEN( "Created with 5 indentions" ) {
         REQUIRE( DataItemString( filename, dataset, number_of_cells, 5 ) == "     <DataItem Format=\"HDF\" NumberType=\"Float\" Precision=\"8\" Dimensions=\"42\">test.h5:/ds</DataItem>" );
      }
   }
}

SCENARIO( "Attributes can be properly created", "[1rank]" ) {
   GIVEN( "An attribute name (eman)" ) {
      std::string const name = "eman";
      WHEN( "Created with empty string" ) {
         REQUIRE( AttributeString( name, "" ) == "<Attribute Name=\"eman\" Center=\"Cell\">\n\n</Attribute>" );
      }
      WHEN( "Created with a string (blubb)" ) {
         std::string const blubb = "blubb";
         REQUIRE( AttributeString( name, blubb ) == "<Attribute Name=\"eman\" Center=\"Cell\">\nblubb\n</Attribute>" );
      }
      WHEN( "Created with indention" ) {
         REQUIRE( AttributeString( name, "", 3 ) == "   <Attribute Name=\"eman\" Center=\"Cell\">\n\n   </Attribute>" );
      }
   }
   GIVEN( "A different attribute name (aemn)" ) {
      std::string const name = "aemn";
      WHEN( "Created with empty string" ) {
         REQUIRE( AttributeString( name, "" ) == "<Attribute Name=\"aemn\" Center=\"Cell\">\n\n</Attribute>" );
      }
   }
}
