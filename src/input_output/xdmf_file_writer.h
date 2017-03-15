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
#ifndef XDMF_FILE_WRITER_H
#define XDMF_FILE_WRITER_H

#include <string>
#include <vector>

/**
 * @brief The XdmfFileWriter class creates an xdmf file to reference the output data of the full simulation, i.e. all time steps in a single file.
 *        The data is only written by MPI rank 0.
 */
class XdmfFileWriter {

   std::string const series_filename_;
   int const mpi_rank_;
   unsigned int const number_of_fluids_;
   bool const print_time_series_;

   //NH We use private static members to cheat the const-ness. Making the class a Singelton is therefore up for discussion.
   static std::vector<double> output_times_;

   void AppendTimestep( std::string const hdf5_filename, unsigned int const global_number_of_cells,
                        unsigned int const global_number_of_vertices ) const;
   void CreateTimestepFile( std::string const& hdf5_file_with_path, unsigned int const global_number_of_cells,
                            unsigned int const global_number_of_vertices, double const time ) const;

public:
   XdmfFileWriter() = delete;
   explicit XdmfFileWriter(std::string const filename, int const mpi_rank, bool const print_time_series = true, unsigned int const number_of_fluids = 1);
   ~XdmfFileWriter();
   XdmfFileWriter( XdmfFileWriter const& ) = delete;
   XdmfFileWriter& operator=( XdmfFileWriter const& ) = delete;
   XdmfFileWriter( XdmfFileWriter&& ) = delete;
   XdmfFileWriter& operator=( XdmfFileWriter&& ) = delete;

   void WriteXdmfForTimestep(double const timestep, std::string const hdf5_filename, unsigned int const global_number_of_cells, unsigned int const global_number_of_vertices) const;
   void WriteDebugXdmfFile(double const timestep, std::string const hdf5_filename, unsigned int const global_number_of_cells, unsigned int const global_number_of_vertices) const;
};

#endif // XDMF_FILE_WRITER_H
