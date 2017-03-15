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
#ifndef HDF5_UTILS_H
#define HDF5_UTILS_H

#include <hdf5.h>
#include <string>
#include <vector>

namespace Hdf5Utilities {

   /**
    * @brief Writes the given values in according to the given specifications into an hdf5 file into the specified hdf-group.
    * @param group Hdf5 group into which the dataset is to be placed.
    * @param name Name of the dataset.
    * @param filespace Hdf5 filespace used when writing the dataset.
    * @param property_list_create_dataset Properties used when creating the dataset.
    * @param property_list_write_dataset Properties used when when writing the dataset.
    * @param write_count Vector of the number of elements to be written into the dataset (per dataset dimension).
    * @param write_start Vector indicating how many elements the given values shall be offsetted with respect to the start of the dataset.
    * @param values Vector of the values to be written into the file
    * @param datatype The Hdf5 datatype.
    * @note * Does not perform sanity checks, inputs need to be proper, in particular the created files/group/memspaces and the data
    *        dimension and offset must match.
    */
   template<typename T>
   void WriteValuesIntoHdf5Group( hid_t const group, std::string const& name, hid_t const& filespace, hid_t const& property_list_create_dataset,
                                  hid_t const& property_list_write_dataset, std::vector<hsize_t> const& write_count,
                                  std::vector<hsize_t> const& write_start, std::vector<T> const& values,
                                  hid_t const& datatype = H5T_NATIVE_DOUBLE ) {

      hid_t const dataset_to_write_into = H5Dcreate2( group, name.c_str(), datatype, filespace, H5P_DEFAULT, property_list_create_dataset,
                                                      H5P_DEFAULT );
      hid_t const memspace = H5Screate_simple( 1, write_count.data(), NULL );
      hid_t const slap = H5Dget_space( dataset_to_write_into );
      H5Sselect_hyperslab( slap, H5S_SELECT_SET, write_start.data(), 0, write_count.data(), NULL );
      H5Dwrite( dataset_to_write_into, H5T_NATIVE_DOUBLE, memspace, slap, property_list_write_dataset, values.data() );
   
      H5Sclose( memspace );
      H5Sclose( slap );
      H5Dclose( dataset_to_write_into );
   }
}

#endif // HDF5_UTILS_H
