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

#include "helper_functions.h"

SCENARIO( "File extensions can be changed", "[1rank]" ) {
   GIVEN( "A filename with absolute path and .txt extension" ) {
      std::string const txt_filename_with_path = "/scratch/directory/file.txt";
      THEN( "The file extension changes to .oop" ) {
         REQUIRE( ChangeFileExtension( txt_filename_with_path, ".oop" ) == "/scratch/directory/file.oop" );
      }
   }
   GIVEN( "A filename with absolute paths containing dots" ) {
      std::string const txt_filename_with_path = "/scratch/directory/fi.le.txt";
      THEN( "The file extension should change to .oop and the string between the dots is not changed" ) {
         REQUIRE( ChangeFileExtension( txt_filename_with_path, ".oop" ) == "/scratch/directory/fi.le.oop" );
      }
   }
   GIVEN( "A filename with long extension" ) {
      std::string const txt_filename_with_path = "/scratch/directory/file.longext";
      THEN( "The file extension changes to .oop" ) {
         REQUIRE( ChangeFileExtension( txt_filename_with_path, ".oop" ) == "/scratch/directory/file.oop" );
      }
   }
}

SCENARIO( "File extensions can be removed", "[1rank]" ) {
   GIVEN( "A filename with absolute path but no extension" ) {
      std::string const txt_filename_with_path = "/scratch/directory/file";
      THEN( "The new filename is the original one" ) {
         REQUIRE( RemoveFileExtension( txt_filename_with_path ) == "/scratch/directory/file" );
      }
   }
   GIVEN( "A filename with absolute path and .txt extension" ) {
      std::string const txt_filename_with_path = "/scratch/directory/file.txt";
      THEN( "The new filename does not contain the .txt extension" ) {
         REQUIRE( RemoveFileExtension( txt_filename_with_path ) == "/scratch/directory/file" );
      }
   }
   GIVEN( "A filename with absolute paths containing dots" ) {
      std::string const txt_filename_with_path = "/scratch/directory/fi.le.txt";
      THEN( "The new filename does not contain the .txt extension, but the intermediate dot" ) {
         REQUIRE( RemoveFileExtension( txt_filename_with_path ) == "/scratch/directory/fi.le" );
      }
   }
}

SCENARIO( "File paths can be removed", "[1rank]" ) {
   GIVEN( "A filename with absolute path but no extension" ) {
      std::string const txt_filename_with_path = "/scratch/directory/file";
      THEN( "Eversthing up to the last /-character is be removed" ) {
         REQUIRE( RemoveFilePath( txt_filename_with_path ) == "file" );
      }
   }
   GIVEN( "A filename with absolute path .txt extension" ) {
      std::string const txt_filename_with_path = "/scratch/directory/file.txt";
      THEN( "Eversthing up to the last /-character is be removed" ) {
         REQUIRE( RemoveFilePath( txt_filename_with_path ) == "file.txt" );
      }
   }
   GIVEN( "A filename with absolute paths containing dots" ) {
      std::string const txt_filename_with_path = "/scratch/directory/fi.le.txt";
      THEN( "Eversthing up to the last /-character is be removed" ) {
         REQUIRE( RemoveFilePath( txt_filename_with_path ) == "fi.le.txt" );
      }
   }
}

SCENARIO( "Numeric values can be converted into scientific notation with given number of digits", "[1rank]" ) {
   GIVEN( "Numeric value with 5 pre-digits and 5 digits (all different)" ) {
      constexpr double value = 12345.06789;
      WHEN( "No digit is desired" ) {
         THEN( "The string should be equal to 1e+04" ) {
            REQUIRE( ToScientificNotationString( value, 0 ) == "1e+04" );
         }
      }
      WHEN( "4 digits are desired" ) {
         THEN( "The string should be equal to 1.2345e+04" ) {
            REQUIRE( ToScientificNotationString( value, 4 ) == "1.2345e+04" );
         }
      }
      WHEN( "8 digits are desired" ) {
         THEN( "The string should be equal to 1.23450679e+04" ) {
            REQUIRE( ToScientificNotationString( value, 8 ) == "1.23450679e+04" );
         }
      }
      WHEN( "12 digits are desired" ) {
         THEN( "The string should be equal to 1.234506789000e+04" ) {
            REQUIRE( ToScientificNotationString( value, 12 ) == "1.234506789000e+04" );
         }
      }
   }
}