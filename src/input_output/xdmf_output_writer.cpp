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
#include "xdmf_output_writer.h"

#include <algorithm> //lower_bound, sort

#include <hdf5.h>
#include <numeric>
#include "topology/id_information.h"
#include "block.h"
#include "levelset/multi_phase_manager/two_phase_manager.h"
#include "user_specifications/compile_time_constants.h"
#include "user_specifications/debug_and_profile_setup.h"
#include "mathematical_functions.h"
#include "enums/interface_tag_definition.h"
#include "enums/vertex_filter_type.h"
#include "hdf5_utilities.h"
#include "communication/mpi_utilities.h"

/**
 * @brief Creates an object to get the simulation data from the RAM to the hard disk.
 * @param flower The tree to read out the local data from.
 * @param topology The topology to get information about the global nodes of the simulation on all ranks.
 * @param setup The setup provides information about user-settings of the simulation.
 * @param output_folder Path to the output folder of the current simulation.
 * @param simulation_name Name of the simulation. Used for the XDMF file.
 */
XdmfOutputWriter::XdmfOutputWriter(Tree const& flower, TopologyManager const& topology, SimulationSetup const& setup,
   std::string const& output_folder, std::string const& simulation_name) :
   // Start initializer list
   OutputWriter(flower, topology, setup, output_folder),
   xdmf_file_name_(output_folder + OutputSubfolderName() + "/" + simulation_name + ".xdmf"),
   xdmf_file_writer_( xdmf_file_name_, MpiUtilities::MyRankId(), true, setup_.NumberOfFluids() ),
   mesh_generator_( topology_ ) {
   /*+ Empty besides initializer list */
}

/**
 * @brief Triggers the output of the simulation results. Based on user Input the correct type of output is created.
 * @param output_time The time at which the simulation is currently at.
 *
 * @note All leaves are written without the halo cells
 */
void XdmfOutputWriter::DoWriteOutputFile( double const output_time ) const {
   // Initial logging
   double write_output_start_time = MPI_Wtime();
   if constexpr( DP::Profile() ) {
      std::string message = "You enabled the Output Vertex Filter: " + VertexFilterTypeToString( CC::VERTEX_FILTER() );
      logger_.LogMessage(message);
   }

   std::string const hdf_filename = output_folder_ + OutputSubfolderName() + "/data_" + std::to_string(output_time*setup_.GetTimeNamingFactor()) + ".h5";

   int const rank = MpiUtilities::MyRankId();
   long long unsigned int const offset = topology_.LeafOffsetOfRank( rank );

   std::pair<unsigned int, unsigned int> node_leaf_count = topology_.NodeAndLeafCount();
   unsigned int const number_of_leaves_globally = node_leaf_count.second;

   auto const& leaves = tree_.Leaves();
   std::size_t const number_of_leaves_locally = leaves.size();

   /** Open the hdf5 file */
   // hdf5 file where the content is written into
   hid_t property_list_create_file = H5Pcreate(H5P_FILE_ACCESS);
   H5Pset_fapl_mpio(property_list_create_file,MPI_COMM_WORLD,MPI_INFO_NULL);
   hid_t file = H5Fcreate(hdf_filename.c_str(),H5F_ACC_TRUNC,H5P_DEFAULT,property_list_create_file);

   /** Write complete mesh topology information into the hdf5 file */
   /** Open mesh topology group */
   hid_t property_list_create_domain_dataset = H5Pcreate(H5P_DATASET_CREATE);
   hid_t property_list_write_domain_dataset = H5Pcreate(H5P_DATASET_XFER);
   H5Pset_dxpl_mpio(property_list_write_domain_dataset,H5FD_MPIO_INDEPENDENT);
   // Put a group for the mesh information into the hdf5 file
   hid_t group_domain = H5Gcreate(file, "domain", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

   /** Vertex IDs */
   // Get the vertex IDs from the mesh generator
   std::vector<unsigned long long int> cell_vertex_ids = mesh_generator_.VertexIds( number_of_leaves_locally );
   // Specify parameter of the dataset the vertex IDs are written into
   std::vector<hsize_t> global_dimensions_cell_vertex_ids = {number_of_leaves_globally * NumberOfCellsPerBlock(), 8};
   std::vector<hsize_t> vertex_ids_write_count = {cell_vertex_ids.size() / 8, 8}; // Multiple of 8 is ensured.
   std::vector<hsize_t> vertex_ids_write_start = {offset * NumberOfCellsPerBlock(), 0};
   // Open the dataset with the given dimensions
   hid_t filespace_vertex_ids  = H5Screate_simple(global_dimensions_cell_vertex_ids.size()   ,global_dimensions_cell_vertex_ids.data()   ,NULL);
   hid_t dataset_vertex_ids = H5Dcreate2(group_domain,"cell_vertices",H5T_NATIVE_ULLONG,filespace_vertex_ids,H5P_DEFAULT,property_list_create_domain_dataset,H5P_DEFAULT);
   hid_t memspace_vertex_ids = H5Screate_simple(2,vertex_ids_write_count.data(),NULL);
   hid_t slap_vertex_ids = H5Dget_space(dataset_vertex_ids);
   // Write the local data into the correct hyperslab position of global dataset
   H5Sselect_hyperslab(slap_vertex_ids,H5S_SELECT_SET,vertex_ids_write_start.data(),0,vertex_ids_write_count.data(),NULL);
   H5Dwrite(dataset_vertex_ids,H5T_NATIVE_ULLONG,memspace_vertex_ids,slap_vertex_ids,property_list_write_domain_dataset,cell_vertex_ids.data());
   // Release the dataspace occupied by the dataset
   H5Sclose(filespace_vertex_ids);
   H5Sclose(memspace_vertex_ids);
   H5Sclose(slap_vertex_ids);
   H5Dclose(dataset_vertex_ids);

   /** Vertex coordinates */
   // Get the vertex coordinates from the mesh generator
   std::vector<double> vertex_coordinates = mesh_generator_.VertexCoordinates<CC::VERTEX_FILTER()>( setup_.LevelZeroBlockSize() );
   // Dimensionalizing
   std::transform( std::begin( vertex_coordinates ), std::end( vertex_coordinates ), std::begin( vertex_coordinates ),
                   [=]( auto const& element ) { return setup_.DimensionalizeLength( element ); } );
   // Specify parameter of the dataset the vertex coordinates are written into
   std::vector<hsize_t> global_dimensions_vertex_coordinates = mesh_generator_.GlobalVertexCoordinatesDimensions<CC::VERTEX_FILTER()>();
   auto const write_start_and_count = mesh_generator_.WriteStartAndCountCoordinates<CC::VERTEX_FILTER()>( number_of_leaves_locally );
   std::vector<hsize_t> coordinates_write_start = std::get<0>( write_start_and_count );
   std::vector<hsize_t> coordinates_write_count = std::get<1>( write_start_and_count );
   // Open the dataset with the given dimensions
   hid_t filespace_coordinates = H5Screate_simple(global_dimensions_vertex_coordinates.size(),global_dimensions_vertex_coordinates.data(),NULL);
   hid_t dataset_coordinates = H5Dcreate2(group_domain,"vertex_coordinates",H5T_NATIVE_DOUBLE,filespace_coordinates,H5P_DEFAULT,property_list_create_domain_dataset,H5P_DEFAULT);
   hid_t memspace_coordinates = H5Screate_simple(2,coordinates_write_count.data(),NULL);
   hid_t slap_coordinates = H5Dget_space(dataset_coordinates);
   // Write the local data into the correct hyperslab position of global dataset
   H5Sselect_hyperslab(slap_coordinates,H5S_SELECT_SET,coordinates_write_start.data(),0,coordinates_write_count.data(),NULL);
   H5Dwrite(dataset_coordinates,H5T_NATIVE_DOUBLE,memspace_coordinates,slap_coordinates,property_list_write_domain_dataset,vertex_coordinates.data());
   // Release the dataspace occupied by the dataset
   H5Pclose(property_list_create_domain_dataset);
   H5Pclose(property_list_write_domain_dataset);
   H5Sclose(filespace_coordinates);
   H5Sclose(memspace_coordinates);
   H5Sclose(slap_coordinates);
   // Close and release memory of the topology group
   H5Dclose(dataset_coordinates);
   H5Gclose(group_domain);
   cell_vertex_ids.clear();
   vertex_coordinates.clear();

   /** Write cell fields one by one into the hdf5 file */
   // Specify vector and counter where the node data is written into. Hdf5 requires scalar fields beeing present in a
   // one-dimensional buffer (not 3D array).
   std::vector<double> simulation_data( number_of_leaves_locally * NumberOfCellsPerBlock() );
   unsigned long long int simulation_data_counter = 0;

   /** Open simulation data group */
   hid_t group_simulation = H5Gcreate2(file,"simulation",H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
   // Specify dimensions for each data field
   std::vector<hsize_t> global_dimensions_simulation = { number_of_leaves_globally * NumberOfCellsPerBlock() };
   std::vector<hsize_t> const simulation_write_count = {simulation_data.size()};
   std::vector<hsize_t> const simulation_start = {offset * NumberOfCellsPerBlock()};
   hid_t filespace_simulation = H5Screate_simple(global_dimensions_simulation.size(),global_dimensions_simulation.data(),NULL);
   // Specify properties for writing the data into the group (parallel writing)
   hid_t property_list_create_dataset = H5Pcreate(H5P_DATASET_CREATE);
   hid_t property_list_write_dataset = H5Pcreate(H5P_DATASET_XFER);
   H5Pset_dxpl_mpio(property_list_write_dataset,H5FD_MPIO_INDEPENDENT);

   // Start writing each field separately
   std::string dataset_name;

   /** Writing Fluid fields (prime states) */
   for(PrimeState const& prime_state : FF::ASOP()) {
      dataset_name = FF::FieldOutputName( prime_state );
      if( dataset_name.empty() ) continue; // skip prime states for which output is disabled (by empty name)
      // Loop through all local leaves
      for(Node const& node : leaves) {
         for(unsigned int k = CC::FICZ(); k <= CC::LICZ(); ++k) {
            for(unsigned int j = CC::FICY(); j <= CC::LICY(); ++j ) {
               for(unsigned int i = CC::FICX(); i <= CC::LICX(); ++i) {
                  MaterialName material = MaterialSignCapsule::PositiveFluidMaterial();  // Initially set to positive for simplicity.
                  std::int8_t const (&interface_tags)[CC::TCX()][CC::TCY()][CC::TCZ()] = node.GetInterfaceTags();
                  // Depending on the interface tags in the current cell, take the positive or negative material
                  if(std::abs(interface_tags[i][j][k]) <= ITTI(IT::NewCutCell)) {
                     LevelsetBlock const& levelset = node.GetLevelsetBlock();
                     double const (&phi)[CC::TCX()][CC::TCY()][CC::TCZ()] = levelset.GetPhi();
                     if(phi[i][j][k] < 0.0) {
                        material = MaterialSignCapsule::NegativeFluidMaterial();
                     }
                     //DO stuff with levelset (must be available as we go over leaves with tag Cutcell)
                  } else if(interface_tags[i][j][k] < 0) {
                     // Search for negative fluid in all_blocks and print it.
                     material = MaterialSignCapsule::NegativeFluidMaterial();
                  } // All other cases imply positive fluid, which is the inital one.

                  /* NH Fetching the buffer anew for every cell is inefficent, but WAY less cumbersome than the alternative
                   * of pre-evaluating where to place the value of which fuild.
                   */
                  Block const& real_fluid = node.GetPhaseByMaterial(material); // throws error if material does not exist in the node
                  double const (&prime_values)[CC::TCX()][CC::TCY()][CC::TCZ()] = real_fluid.GetPrimeStateBuffer(prime_state);

                  // get re-dimensionalized value
                  simulation_data[simulation_data_counter++] = setup_.DimensionalizeValue(prime_values[i][j][k], FF::FieldUnit( prime_state ));
               }
            }
         }
      }
      simulation_data_counter = 0;

      // After each prime state write the data into the hdf5 file at the appropriate hyperslab positions
      Hdf5Utilities::WriteValuesIntoHdf5Group( group_simulation, dataset_name, filespace_simulation, property_list_create_dataset,
                                               property_list_write_dataset, simulation_write_count, simulation_start, simulation_data );

   }

   /** Writing Levelset field (phi) */
   dataset_name = "levelset";
   // Loop through all local leaves
   for(Node const& node : leaves) {
      if(node.HasLevelset()) {
         LevelsetBlock const& levelset_block = node.GetLevelsetBlock();
         double const (&phi)[CC::TCX()][CC::TCY()][CC::TCZ()] = levelset_block.GetPhi();
         for(unsigned int k = CC::FICZ(); k <= CC::LICZ(); ++k) {
            for(unsigned int j = CC::FICY(); j <= CC::LICY(); ++j) {
               for(unsigned int i = CC::FICX(); i <= CC::LICX(); ++i) {
                  simulation_data[simulation_data_counter++] = phi[i][j][k];
               }
            }
         }
      } else {
         std::int8_t const (&interface_tags)[CC::TCX()][CC::TCY()][CC::TCZ()] = node.GetInterfaceTags();
         // We just have to insert cut-off with the correct sign in all cells
         double cut_off = Signum(interface_tags[0][0][0])*CC::LSCOF();
         for(unsigned int k = CC::FICZ(); k <= CC::LICZ(); ++k) {
            for(unsigned int j = CC::FICY(); j <= CC::LICY(); ++j) {
               for(unsigned int i = CC::FICX(); i <= CC::LICX(); ++i) {
                  simulation_data[simulation_data_counter++] = cut_off;
               }
            }
         }
      }
   }
   simulation_data_counter = 0;

   // Write full levelset data into the hdf5 file at the appropriate hyperslab positions
   Hdf5Utilities::WriteValuesIntoHdf5Group( group_simulation, dataset_name, filespace_simulation, property_list_create_dataset,
                                            property_list_write_dataset, simulation_write_count, simulation_start, simulation_data );

   /** Writing Interface quantity fields */
   for( InterfaceQuantity const iq : FF::ASIQ() ) {
      dataset_name = FF::FieldOutputName( iq );
      if( dataset_name.empty() ) continue; // skip interface quantities for which output is disabled (by empty name)
      for(Node const& node : leaves) {
         if(node.HasLevelset()) {
            LevelsetBlock const& levelset_block = node.GetLevelsetBlock();
            double const (&interface_quantity)[CC::TCX()][CC::TCY()][CC::TCZ()] = levelset_block.GetInterfaceQuantityBuffer(iq);
            for(unsigned int k = CC::FICZ(); k <= CC::LICZ(); ++k) {
               for(unsigned int j = CC::FICY(); j <= CC::LICY(); ++j) {
                  for(unsigned int i = CC::FICX(); i <= CC::LICX(); ++i) {
                     // get re-dimensionalized value
                     simulation_data[simulation_data_counter++] = setup_.DimensionalizeValue(interface_quantity[i][j][k], FF::FieldUnit( iq ));
                  }
               }
            }
         } else {
            for(unsigned int k = CC::FICZ(); k <= CC::LICZ(); ++k) {
               for(unsigned int j = CC::FICY(); j <= CC::LICY(); ++j) {
                  for(unsigned int i = CC::FICX(); i <= CC::LICX(); ++i) {
                     simulation_data[simulation_data_counter++] = 0.0;
                  }
               }
            }
         }
      }
      simulation_data_counter = 0;

      // Write full interface quantity data into the hdf5 file at the appropriate hyperslab positions
      Hdf5Utilities::WriteValuesIntoHdf5Group( group_simulation, dataset_name, filespace_simulation, property_list_create_dataset,
                                               property_list_write_dataset, simulation_write_count, simulation_start, simulation_data );
   }

   /** Writing partition (rank) data */
   dataset_name = "partition";
   double const rank_in_double_format = double( rank );
   //only number of nodes is important since the same for all local nodes
   for(unsigned int n = 0; n < leaves.size(); ++n) {
      for(unsigned int k = 0; k < CC::ICZ(); ++k) {
         for(unsigned int j = 0; j < CC::ICY(); ++j ) {
            for(unsigned int i = 0; i < CC::ICX(); ++i) {
               simulation_data[simulation_data_counter++] = rank_in_double_format;
            }
         }
      }
   } //nodes

   simulation_data_counter = 0;

   // Write full partition data into the hdf5 file at the appropriate hyperslab positions
   Hdf5Utilities::WriteValuesIntoHdf5Group( group_simulation, dataset_name, filespace_simulation, property_list_create_dataset,
                                            property_list_write_dataset, simulation_write_count, simulation_start, simulation_data );

   /** Releasing group data */
   H5Pclose(property_list_create_dataset);
   H5Pclose(property_list_write_dataset);
   H5Sclose(filespace_simulation);
   H5Gclose(group_simulation);

   /** Releasing file data */
   H5Pclose(property_list_create_file);
   H5Fclose(file);

   /** Write the Xdmf file for the created Hdf5 file */
   xdmf_file_writer_.WriteXdmfForTimestep( output_time, hdf_filename, number_of_leaves_globally * NumberOfCellsPerBlock(), global_dimensions_vertex_coordinates[0]);

   /** Final debugging messages */
   if( DP::Profile() ) {
      double const seconds_elapsed = MPI_Wtime() - write_output_start_time;
      std::vector<double> all_times;
      int const number_of_ranks = MpiUtilities::NumberOfRanks();
      all_times.resize( number_of_ranks );
      //communicate all individual times
      MPI_Allgather( &seconds_elapsed, 1, MPI_DOUBLE, all_times.data(), 1, MPI_DOUBLE, MPI_COMM_WORLD );
      std::vector<double>::iterator max_time = max_element(all_times.begin(), all_times.end());
      std::vector<double>::iterator min_time = min_element(all_times.begin(), all_times.end());
      double const avg_time = accumulate(all_times.begin(), all_times.end(), 0.0) / number_of_ranks;
      std::string message = "Writing output: avg/max(rank)/min(rank) ";
      message.append(std::to_string(avg_time));
      message.append("/");
      message.append(std::to_string(*max_time));
      message.append("(");
      message.append(std::to_string(distance(all_times.begin(), max_time)));
      message.append(")/");
      message.append(std::to_string(*min_time));
      message.append("(");
      message.append(std::to_string(distance(all_times.begin(), min_time)));
      message.append(") seconds elapsed.");
      logger_.LogMessage(message);
   }

}

/**
 * @brief Triggers the creation of an outputfile in the format specified in the inputfile.
 * @param debug_time The 'time' of this debug output. Commonly only an integer number.
 *
 * @note Compared to the standard output, in the debug output all nodes on all levels are written with the halo cells (not only leaves).
 */
void XdmfOutputWriter::DoWriteDebugFile(double const debug_time) const {

   // Same logic as in standard XDMF output. See there for details.
   std::string const hdf_filename = output_folder_ + DebugSubfolderName() + "/debug_" + std::to_string(debug_time*setup_.GetTimeNamingFactor()) + ".h5";
   hid_t property_list_create_file = H5Pcreate(H5P_FILE_ACCESS);
   H5Pset_fapl_mpio(property_list_create_file,MPI_COMM_WORLD,MPI_INFO_NULL);
   hid_t file = H5Fcreate(hdf_filename.c_str(),H5F_ACC_TRUNC,H5P_DEFAULT,property_list_create_file);

   // We determine an offset for the numbering of the vertices.
   std::vector<std::pair<unsigned int,unsigned int>> const rank_node_map = topology_.NodesAndLeavesPerRank();
   long long unsigned const my_rank = MpiUtilities::MyRankId();

   long long unsigned int const offset = std::accumulate(rank_node_map.begin(),rank_node_map.begin() + my_rank, 0ll,
      [](unsigned int const& a, const std::pair<unsigned int,unsigned int>& b){return a + b.first;});

   std::pair<unsigned int, unsigned int> const node_leaf_count = topology_.NodeAndLeafCount();
   // Specify dimensions of coordinates and IDs
   unsigned int const number_of_nodes_globally = node_leaf_count.first;
   unsigned int const number_of_nodes_locally = rank_node_map[my_rank].first;
   unsigned int const number_of_cells_per_block = CC::TCX()*CC::TCY()*CC::TCZ();
   unsigned int const number_of_vertices_per_block = (CC::TCX()+1)*(CC::TCY()+1)*(CC::TCZ()+1);
   // Define vector holding the vertex IDs and cooridnates
   std::vector<double> vertex_coordinates(number_of_nodes_locally * number_of_vertices_per_block * 3); //We have #Nodes with #Vertecies in 3 dimensions.
   unsigned long long int vertex_coordinates_counter = 0;
   std::vector<unsigned long long int> cell_vertex_ids(number_of_nodes_locally*number_of_cells_per_block*8); //We have #Nodes with #Cells, cells are build by connecting 8 vertecies.
   unsigned long long int cell_vertex_id_counter = 0;

   /** Comptue the vertex coordinates and IDs */
   unsigned int const j_ids_skew = CC::TCX()+1;
   unsigned int const k_ids_skew = (CC::TCX()+1)*(CC::TCY()+1);

   unsigned int leaves_counter = 0;
   double z_offset = 0;
   unsigned long long shift;

   for(auto const& level : tree_.FullNodeList()) {
      for(auto const& [id, node]: level) {
         double const block_size  = setup_.DimensionalizeLength(node.GetBlockSize());
         std::array<double, 3> const coordinates = DomainCoordinatesOfId(id, block_size);

         /*** Coordinates ***/
         for(unsigned int k = 0; k <= CC::TCZ(); ++k) {
            for(unsigned int j = 0; j <= CC::TCY(); ++j ) {
               for(unsigned int i = 0; i <= CC::TCX(); ++i) {
                  //NH TCX is correct for all dimensions. The "+1" gives a space between any two blocks, which is nice for debugging visualization.
                  vertex_coordinates[vertex_coordinates_counter    ] =            coordinates[0] + (double(i)/double(CC::TCX()+1)) * block_size;
                  vertex_coordinates[vertex_coordinates_counter + 1] =            coordinates[1] + (double(j)/double(CC::TCX()+1)) * block_size;
                  vertex_coordinates[vertex_coordinates_counter + 2] = z_offset + coordinates[2] + (double(k)/double(CC::TCX()+1)) * block_size;
                  // update the counter skipping the 3 elements
                  vertex_coordinates_counter += 3;
               }
            }
         }

         /*** Vertex Ids ***/
         for(unsigned int k = 0; k < CC::TCZ(); ++k) {
            for(unsigned int j = 0; j < CC::TCY(); ++j ) {
               for(unsigned int i = 0; i < CC::TCX(); ++i) {
                  shift = leaves_counter * number_of_vertices_per_block + offset*number_of_vertices_per_block;
                  cell_vertex_ids[cell_vertex_id_counter    ] =  i    +  j    * j_ids_skew + k     * k_ids_skew + shift;
                  cell_vertex_ids[cell_vertex_id_counter + 1] = (i+1) +  j    * j_ids_skew + k     * k_ids_skew + shift;
                  cell_vertex_ids[cell_vertex_id_counter + 2] = (i+1) + (j+1) * j_ids_skew + k     * k_ids_skew + shift;
                  cell_vertex_ids[cell_vertex_id_counter + 3] =  i    + (j+1) * j_ids_skew + k     * k_ids_skew + shift;
                  cell_vertex_ids[cell_vertex_id_counter + 4] =  i    +  j    * j_ids_skew + (k+1) * k_ids_skew + shift;
                  cell_vertex_ids[cell_vertex_id_counter + 5] = (i+1) +  j    * j_ids_skew + (k+1) * k_ids_skew + shift;
                  cell_vertex_ids[cell_vertex_id_counter + 6] = (i+1) + (j+1) * j_ids_skew + (k+1) * k_ids_skew + shift;
                  cell_vertex_ids[cell_vertex_id_counter + 7] =  i    + (j+1) * j_ids_skew + (k+1) * k_ids_skew + shift;
                  cell_vertex_id_counter += 8;
               }
            }
         }

         leaves_counter++;
      } // Nodes
      z_offset += (setup_.GetLevelZeroBlocksZ() + 1) * (setup_.LevelZeroBlockSize());
   }

   /*** Writing Mesh Data and Freeing Memory ***/
   // IDs
   hid_t group_domain = H5Gcreate2(file,"domain",H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
   hid_t property_list_create_domain_dataset = H5Pcreate(H5P_DATASET_CREATE);
   hid_t property_list_write_domain_dataset = H5Pcreate(H5P_DATASET_XFER);
   H5Pset_dxpl_mpio(property_list_write_domain_dataset,H5FD_MPIO_INDEPENDENT);

   std::vector<hsize_t> global_dimensions_cell_vertex_ids    = {number_of_nodes_globally*number_of_cells_per_block   ,8};
   std::vector<hsize_t> vertex_ids_write_count = {cell_vertex_ids.size() / 8, 8}; // Multiple of 8 is ensured (see above) and due to XDMF we artificially split the vector.
   std::vector<hsize_t> vertex_ids_write_start = {offset*number_of_cells_per_block, 0};

   hid_t filespace_vertex_ids  = H5Screate_simple(global_dimensions_cell_vertex_ids.size()   ,global_dimensions_cell_vertex_ids.data()   ,NULL);
   hid_t dataset_vertex_ids = H5Dcreate2(group_domain,"cell_vertices",H5T_NATIVE_ULLONG,filespace_vertex_ids ,H5P_DEFAULT,property_list_create_domain_dataset,H5P_DEFAULT);
   hid_t memspace_vertex_ids = H5Screate_simple(2,vertex_ids_write_count.data(),NULL);
   hid_t slap_vertex_ids = H5Dget_space(dataset_vertex_ids);

   H5Sselect_hyperslab(slap_vertex_ids,H5S_SELECT_SET,vertex_ids_write_start.data(),0,vertex_ids_write_count.data(),NULL);
   H5Dwrite(dataset_vertex_ids,H5T_NATIVE_ULLONG,memspace_vertex_ids,slap_vertex_ids,property_list_write_domain_dataset,cell_vertex_ids.data());

   H5Sclose(slap_vertex_ids);
   H5Sclose(memspace_vertex_ids);
   H5Dclose(dataset_vertex_ids);
   H5Sclose(filespace_vertex_ids);

   // Coordinates
   std::vector<hsize_t> global_dimensions_vertex_coordinates = {number_of_nodes_globally*number_of_vertices_per_block,3};
   std::vector<hsize_t> coordinates_write_count = {vertex_coordinates.size()/3, 3}; // Multiple of 3 is ensured (see above).
   std::vector<hsize_t> coordinates_write_start = {offset*number_of_vertices_per_block, 0};

   hid_t filespace_coordinates = H5Screate_simple(global_dimensions_vertex_coordinates.size(),global_dimensions_vertex_coordinates.data(),NULL);
   hid_t dataset_coordinates = H5Dcreate2(group_domain,"vertex_coordinates",H5T_NATIVE_DOUBLE,filespace_coordinates,H5P_DEFAULT,property_list_create_domain_dataset,H5P_DEFAULT);
   hid_t memspace_coordinates = H5Screate_simple(2,coordinates_write_count.data(),NULL);
   hid_t slap_coordinates = H5Dget_space(dataset_coordinates);

   H5Sselect_hyperslab(slap_coordinates,H5S_SELECT_SET,coordinates_write_start.data(),0,coordinates_write_count.data(),NULL);
   H5Dwrite(dataset_coordinates,H5T_NATIVE_DOUBLE,memspace_coordinates,slap_coordinates,property_list_write_domain_dataset,vertex_coordinates.data());

   H5Sclose(slap_coordinates);
   H5Sclose(memspace_coordinates);
   H5Dclose(dataset_coordinates);
   H5Sclose(filespace_coordinates);

   // Close group and free memory
   H5Pclose(property_list_create_domain_dataset);
   H5Pclose(property_list_write_domain_dataset);
   H5Gclose(group_domain);
   cell_vertex_ids.clear();
   vertex_coordinates.clear();

   /*** Gathering Simulation (and Partition Data) and Writing to File **/
   std::vector<double> simulation_data(number_of_nodes_locally*number_of_cells_per_block);
   unsigned long long int simulation_data_counter = 0;

   //Preparing the HDF file
   hid_t group_simulation = H5Gcreate2(file,"simulation",H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);

   std::vector<hsize_t> global_dimensions_simulation = {number_of_nodes_globally*number_of_cells_per_block};
   std::vector<hsize_t> simulation_write_count = {simulation_data.size()};
   std::vector<hsize_t> simulation_start = {offset*number_of_cells_per_block};

   hid_t filespace_simulation = H5Screate_simple(global_dimensions_simulation.size(),global_dimensions_simulation.data(),NULL);
   hid_t property_list_create_dataset = H5Pcreate(H5P_DATASET_CREATE);
   hid_t property_list_write_dataset = H5Pcreate(H5P_DATASET_XFER);
   H5Pset_dxpl_mpio(property_list_write_dataset,H5FD_MPIO_INDEPENDENT);

   double dummy_value = -1; // We use this value to output non-existing buffers.

   std::string dataset_name;

   unsigned int material_count = 0;

   hid_t dataset_simulation;
   hid_t memspace_simulation;
   hid_t slap_simulation;

   /** Writing Fluid fields (conservatives, prime states) */
   for(MaterialName const& material : setup_.AllFluids()) {
      for(PrimeState const& prime_state : FF::ASOP()) {
         std::string const prime_name(FF::FieldOutputName( prime_state ));
         if( prime_name.empty() ) continue; // skip prime states with disabled output
         dataset_name = "fluid_" + std::to_string(material_count) + "_" + prime_name;
         for(auto const& level : tree_.FullNodeList()) {
            for(auto const& id_node : level) {
               Node const& node = id_node.second;
               // Use dummy value if material does not exist on node otherwise write real data
               if(node.ContainsMaterial(material)) {
                  Block const& phase = id_node.second.GetPhaseByMaterial(material);
                  double const (&prime_state_values)[CC::TCX()][CC::TCY()][CC::TCZ()] = phase.GetPrimeStateBuffer(prime_state);
                  for(unsigned int k = 0; k < CC::TCZ(); ++k) {
                     for(unsigned int j = 0; j < CC::TCY(); ++j ) {
                        for(unsigned int i = 0; i < CC::TCX(); ++i) {
                           simulation_data[simulation_data_counter++] = setup_.DimensionalizeValue(prime_state_values[i][j][k], FF::FieldUnit( prime_state ));
                        }
                     }
                  }
               } else {
                  for(unsigned int k = 0; k < CC::TCZ(); ++k) {
                     for(unsigned int j = 0; j < CC::TCY(); ++j ) {
                        for(unsigned int i = 0; i < CC::TCX(); ++i) {
                           simulation_data[simulation_data_counter++] = dummy_value;
                        }
                     }
                  }
               }
            } //nodes
         } //level
         simulation_data_counter = 0;

         // Write the data into the hdf5 file
         dataset_simulation = H5Dcreate2(group_simulation,dataset_name.c_str(),H5T_NATIVE_DOUBLE,filespace_simulation,H5P_DEFAULT,property_list_create_dataset,H5P_DEFAULT);
         memspace_simulation = H5Screate_simple(1,simulation_write_count.data(),NULL);
         slap_simulation = H5Dget_space(dataset_simulation);
         H5Sselect_hyperslab(slap_simulation,H5S_SELECT_SET,simulation_start.data(),0,simulation_write_count.data(),NULL);
         H5Dwrite(dataset_simulation,H5T_NATIVE_DOUBLE,memspace_simulation,slap_simulation,property_list_write_dataset,simulation_data.data());

         H5Sclose(memspace_simulation);
         H5Sclose(slap_simulation);
         H5Dclose(dataset_simulation);
      } // Prime states

      for(Equation const& eq : FF::ASOE()) {
         std::string const equation_name(FF::FieldOutputName( eq ));
         if( equation_name.empty() ) continue; // skip prime states with disabled output
         dataset_name = "fluid_" + std::to_string(material_count) + "_" + equation_name + "_avg";
         for(auto const& level : tree_.FullNodeList()) {
            for(auto const& id_node : level) {
               Node const& node = id_node.second;
               if(node.ContainsMaterial(material)) {
                  Block const& phase = id_node.second.GetPhaseByMaterial(material);
                  double const (&avg_values)[CC::TCX()][CC::TCY()][CC::TCZ()] = phase.GetAverageBuffer(eq);
                  for(unsigned int k = 0; k < CC::TCZ(); ++k) {
                     for(unsigned int j = 0; j < CC::TCY(); ++j ) {
                        for(unsigned int i = 0; i < CC::TCX(); ++i) {
                           simulation_data[simulation_data_counter++] = setup_.DimensionalizeValue(avg_values[i][j][k], FF::FieldUnit( eq ));
                        }
                     }
                  }
               } else {
                  for(unsigned int k = 0; k < CC::TCZ(); ++k) {
                     for(unsigned int j = 0; j < CC::TCY(); ++j ) {
                        for(unsigned int i = 0; i < CC::TCX(); ++i) {
                           simulation_data[simulation_data_counter++] = dummy_value;
                        }
                     }
                  }
               }
            } //nodes
         } //level
         simulation_data_counter = 0;

         // Write the data into the hdf5 file
         dataset_simulation = H5Dcreate2(group_simulation,dataset_name.c_str(),H5T_NATIVE_DOUBLE,filespace_simulation,H5P_DEFAULT,property_list_create_dataset,H5P_DEFAULT);
         memspace_simulation = H5Screate_simple(1,simulation_write_count.data(),NULL);
         slap_simulation = H5Dget_space(dataset_simulation);
         H5Sselect_hyperslab(slap_simulation,H5S_SELECT_SET,simulation_start.data(),0,simulation_write_count.data(),NULL);
         H5Dwrite(dataset_simulation,H5T_NATIVE_DOUBLE,memspace_simulation,slap_simulation,property_list_write_dataset,simulation_data.data());

         H5Sclose(memspace_simulation);
         H5Sclose(slap_simulation);
         H5Dclose(dataset_simulation);

      } // Avg Conservatives

      for(Equation const& eq : FF::ASOE()) {
         std::string const equation_name(FF::FieldOutputName( eq ));
         if( equation_name.empty() ) continue; // skip prime states with disabled output
         dataset_name = "fluid_" + std::to_string(material_count) + "_" + equation_name + "_rhs";
         for(auto const& level : tree_.FullNodeList()) {
            for(auto const& id_node : level) {
               Node const& node = id_node.second;
               if(node.ContainsMaterial(material)) {
                  Block const& phase = node.GetPhaseByMaterial(material);
                  double const (&rhs_values)[CC::TCX()][CC::TCY()][CC::TCZ()] = phase.GetRightHandSideBuffer(eq);
                  for(unsigned int k = 0; k < CC::TCZ(); ++k) {
                     for(unsigned int j = 0; j < CC::TCY(); ++j ) {
                        for(unsigned int i = 0; i < CC::TCX(); ++i) {
                           simulation_data[simulation_data_counter++] = setup_.DimensionalizeValue(rhs_values[i][j][k], FF::FieldUnit( eq ));
                        }
                     }
                  }
               } else {
                  for(unsigned int k = 0; k < CC::TCZ(); ++k) {
                     for(unsigned int j = 0; j < CC::TCY(); ++j ) {
                        for(unsigned int i = 0; i < CC::TCX(); ++i) {
                           simulation_data[simulation_data_counter++] = dummy_value;
                        }
                     }
                  }
               }
            } //nodes
         } //level
         simulation_data_counter = 0;

         // Write the data into the hdf5 file
         dataset_simulation = H5Dcreate2(group_simulation,dataset_name.c_str(),H5T_NATIVE_DOUBLE,filespace_simulation,H5P_DEFAULT,property_list_create_dataset,H5P_DEFAULT);
         memspace_simulation = H5Screate_simple(1,simulation_write_count.data(),NULL);
         slap_simulation = H5Dget_space(dataset_simulation);
         H5Sselect_hyperslab(slap_simulation,H5S_SELECT_SET,simulation_start.data(),0,simulation_write_count.data(),NULL);
         H5Dwrite(dataset_simulation,H5T_NATIVE_DOUBLE,memspace_simulation,slap_simulation,property_list_write_dataset,simulation_data.data());

         H5Sclose(memspace_simulation);
         H5Sclose(slap_simulation);
         H5Dclose(dataset_simulation);

      } // Rhs Conservatives

      material_count++;
   } // material

   /** Writing interface tags */
   dataset_name = "interface_tags";
   for(auto const& level : tree_.FullNodeList()) {
      for(auto const& id_node : level) {
         std::int8_t const (&interface_tags)[CC::TCX()][CC::TCY()][CC::TCZ()] = id_node.second.GetInterfaceTags();
         for(unsigned int k = 0; k < CC::TCZ(); ++k) {
            for(unsigned int j = 0; j < CC::TCY(); ++j ) {
               for(unsigned int i = 0; i < CC::TCX(); ++i) {
                  simulation_data[simulation_data_counter++] = interface_tags[i][j][k];
               }
            }
         }
      }
   }
   simulation_data_counter = 0;

   // Write the data into the hdf5 file
   dataset_simulation = H5Dcreate2(group_simulation,dataset_name.c_str(),H5T_NATIVE_DOUBLE,filespace_simulation,H5P_DEFAULT,property_list_create_dataset,H5P_DEFAULT);
   memspace_simulation = H5Screate_simple(1,simulation_write_count.data(),NULL);
   slap_simulation = H5Dget_space(dataset_simulation);
   H5Sselect_hyperslab(slap_simulation,H5S_SELECT_SET,simulation_start.data(),0,simulation_write_count.data(),NULL);
   H5Dwrite(dataset_simulation,H5T_NATIVE_DOUBLE,memspace_simulation,slap_simulation,property_list_write_dataset,simulation_data.data());
   H5Sclose(memspace_simulation);
   H5Sclose(slap_simulation);
   H5Dclose(dataset_simulation);

   /** Writing levelset data (phi and volume fraction) */
   double const debug_levelset_cutoff  = CC::LSCOF()*ITTI(IT::BulkPhase);
   double const one_bulk_phase = 1.0 / ITTI(IT::BulkPhase);

   dataset_name = "levelset";
   for(auto const& level : tree_.FullNodeList()) {
      for(auto const& id_node : level) {
         Node const& node = id_node.second;
         if(node.HasLevelset()) {
            LevelsetBlock const& levelset_block = node.GetLevelsetBlock();
            double const        (&phi)[CC::TCX()][CC::TCY()][CC::TCZ()] = levelset_block.GetPhi();
            for(unsigned int k = 0; k < CC::TCZ(); ++k) {
               for(unsigned int j = 0; j < CC::TCY(); ++j ) {
                  for(unsigned int i = 0; i < CC::TCX(); ++i) {
                     simulation_data[simulation_data_counter++] = phi[i][j][k];
                  }
               }
            }
         } else {
            std::int8_t const (&interface_tags)[CC::TCX()][CC::TCY()][CC::TCZ()] = node.GetInterfaceTags();
            for(unsigned int k = 0; k < CC::TCZ(); ++k) {
               for(unsigned int j = 0; j < CC::TCY(); ++j ) {
                  for(unsigned int i = 0; i < CC::TCX(); ++i) {
                     simulation_data[simulation_data_counter++] = (interface_tags[i][j][k] * one_bulk_phase) * debug_levelset_cutoff;
                  }
               }
            }
         }
      }
   } // phi base
   simulation_data_counter = 0;

   // Write the data into the hdf5 file
   dataset_simulation = H5Dcreate2(group_simulation,dataset_name.c_str(),H5T_NATIVE_DOUBLE,filespace_simulation,H5P_DEFAULT,property_list_create_dataset,H5P_DEFAULT);
   memspace_simulation = H5Screate_simple(1,simulation_write_count.data(),NULL);
   slap_simulation = H5Dget_space(dataset_simulation);
   H5Sselect_hyperslab(slap_simulation,H5S_SELECT_SET,simulation_start.data(),0,simulation_write_count.data(),NULL);
   H5Dwrite(dataset_simulation,H5T_NATIVE_DOUBLE,memspace_simulation,slap_simulation,property_list_write_dataset,simulation_data.data());
   H5Sclose(memspace_simulation);
   H5Sclose(slap_simulation);
   H5Dclose(dataset_simulation);

   dataset_name = "levelset_rhs";
   for(auto const& level : tree_.FullNodeList()) {
      for(auto const& id_node : level) {
         Node const& node = id_node.second;
         if(node.HasLevelset()) {
            LevelsetBlock const& levelset_block = node.GetLevelsetBlock();
            double const    (&phi_rhs)[CC::TCX()][CC::TCY()][CC::TCZ()] = levelset_block.GetPhiRightHandSide();
            for(unsigned int k = 0; k < CC::TCZ(); ++k) {
               for(unsigned int j = 0; j < CC::TCY(); ++j ) {
                  for(unsigned int i = 0; i < CC::TCX(); ++i) {
                     simulation_data[simulation_data_counter++] = phi_rhs[i][j][k];
                  }
               }
            }
         } else {
            std::int8_t const (&interface_tags)[CC::TCX()][CC::TCY()][CC::TCZ()] = node.GetInterfaceTags();
            for(unsigned int k = 0; k < CC::TCZ(); ++k) {
               for(unsigned int j = 0; j < CC::TCY(); ++j ) {
                  for(unsigned int i = 0; i < CC::TCX(); ++i) {
                     simulation_data[simulation_data_counter++] = (interface_tags[i][j][k] * one_bulk_phase) * debug_levelset_cutoff;
                  }
               }
            }
         }
      }
   } // phi rhs
   simulation_data_counter = 0;

   // Write the data into the hdf5 file
   dataset_simulation = H5Dcreate2(group_simulation,dataset_name.c_str(),H5T_NATIVE_DOUBLE,filespace_simulation,H5P_DEFAULT,property_list_create_dataset,H5P_DEFAULT);
   memspace_simulation = H5Screate_simple(1,simulation_write_count.data(),NULL);
   slap_simulation = H5Dget_space(dataset_simulation);
   H5Sselect_hyperslab(slap_simulation,H5S_SELECT_SET,simulation_start.data(),0,simulation_write_count.data(),NULL);
   H5Dwrite(dataset_simulation,H5T_NATIVE_DOUBLE,memspace_simulation,slap_simulation,property_list_write_dataset,simulation_data.data());
   H5Sclose(memspace_simulation);
   H5Sclose(slap_simulation);
   H5Dclose(dataset_simulation);

   dataset_name = "levelset_reinitialized";
   for(auto const& level : tree_.FullNodeList()) {
      for(auto const& id_node : level) {
         Node const& node = id_node.second;
         if(node.HasLevelset()) {
            LevelsetBlock const& levelset_block = node.GetLevelsetBlock();
            double const (&phi_reinit)[CC::TCX()][CC::TCY()][CC::TCZ()] = levelset_block.GetPhiReinitialized();
            for(unsigned int k = 0; k < CC::TCZ(); ++k) {
               for(unsigned int j = 0; j < CC::TCY(); ++j ) {
                  for(unsigned int i = 0; i < CC::TCX(); ++i) {
                     simulation_data[simulation_data_counter++] = phi_reinit[i][j][k];
                  }
               }
            }
         } else {
            std::int8_t const (&interface_tags)[CC::TCX()][CC::TCY()][CC::TCZ()] = node.GetInterfaceTags();
            for(unsigned int k = 0; k < CC::TCZ(); ++k) {
               for(unsigned int j = 0; j < CC::TCY(); ++j ) {
                  for(unsigned int i = 0; i < CC::TCX(); ++i) {
                     simulation_data[simulation_data_counter++] = (interface_tags[i][j][k] * one_bulk_phase) * debug_levelset_cutoff;
                  }
               }
            }
         }
      }
   } // phi reinit
   simulation_data_counter = 0;

   // Write the data into the hdf5 file
   dataset_simulation = H5Dcreate2(group_simulation,dataset_name.c_str(),H5T_NATIVE_DOUBLE,filespace_simulation,H5P_DEFAULT,property_list_create_dataset,H5P_DEFAULT);
   memspace_simulation = H5Screate_simple(1,simulation_write_count.data(),NULL);
   slap_simulation = H5Dget_space(dataset_simulation);
   H5Sselect_hyperslab(slap_simulation,H5S_SELECT_SET,simulation_start.data(),0,simulation_write_count.data(),NULL);
   H5Dwrite(dataset_simulation,H5T_NATIVE_DOUBLE,memspace_simulation,slap_simulation,property_list_write_dataset,simulation_data.data());
   H5Sclose(memspace_simulation);
   H5Sclose(slap_simulation);
   H5Dclose(dataset_simulation);

   // Volume fraction
   dataset_name = "volume_fraction";
   for(auto const& level : tree_.FullNodeList()) {
      for(auto const& id_node : level) {
         Node const& node = id_node.second;
         if(node.HasLevelset()) {
            LevelsetBlock const& levelset_block = node.GetLevelsetBlock();
            double const (&volume_fraction)[CC::TCX()][CC::TCY()][CC::TCZ()] = levelset_block.GetVolumeFraction();
            for(unsigned int k = 0; k < CC::TCZ(); ++k) {
               for(unsigned int j = 0; j < CC::TCY(); ++j ) {
                  for(unsigned int i = 0; i < CC::TCX(); ++i) {
                     simulation_data[simulation_data_counter++] = volume_fraction[i][j][k];
                  }
               }
            }
         } else {
            std::int8_t const (&interface_tags)[CC::TCX()][CC::TCY()][CC::TCZ()] = node.GetInterfaceTags();
            for(unsigned int k = 0; k < CC::TCZ(); ++k) {
               for(unsigned int j = 0; j < CC::TCY(); ++j ) {
                  for(unsigned int i = 0; i < CC::TCX(); ++i) {
                     simulation_data[simulation_data_counter++] = (interface_tags[i][j][k] * one_bulk_phase) * 1.0;
                  }
               }
            }
         }
      }
   }
   simulation_data_counter = 0;

   // Write the data into the hdf5 file
   dataset_simulation = H5Dcreate2(group_simulation,dataset_name.c_str(),H5T_NATIVE_DOUBLE,filespace_simulation,H5P_DEFAULT,property_list_create_dataset,H5P_DEFAULT);
   memspace_simulation = H5Screate_simple(1,simulation_write_count.data(),NULL);
   slap_simulation = H5Dget_space(dataset_simulation);
   H5Sselect_hyperslab(slap_simulation,H5S_SELECT_SET,simulation_start.data(),0,simulation_write_count.data(),NULL);
   H5Dwrite(dataset_simulation,H5T_NATIVE_DOUBLE,memspace_simulation,slap_simulation,property_list_write_dataset,simulation_data.data());
   H5Sclose(memspace_simulation);
   H5Sclose(slap_simulation);
   H5Dclose(dataset_simulation);

   /** Writing interface quantities */
   for( InterfaceQuantity const iq : FF::ASIQ() ) {
      dataset_name = FF::FieldOutputName( iq );
      if( dataset_name.empty() ) continue; // skip interface quantities with disabled output
      for(auto const& level : tree_.FullNodeList()) {
         for(auto const& id_node : level) {
            Node const& node = id_node.second;
            if(node.HasLevelset()) {
               LevelsetBlock const& levelset_block = node.GetLevelsetBlock();
               double const (&interface_quantity)[CC::TCX()][CC::TCY()][CC::TCZ()] = levelset_block.GetInterfaceQuantityBuffer(iq);
               for(unsigned int k = 0; k < CC::TCZ(); ++k) {
                  for(unsigned int j = 0; j < CC::TCY(); ++j ) {
                     for(unsigned int i = 0; i < CC::TCX(); ++i) {
                        simulation_data[simulation_data_counter++] = setup_.DimensionalizeValue(interface_quantity[i][j][k], FF::FieldUnit( iq ));
                     }
                  }
               }
            } else {
               for(unsigned int k = 0; k < CC::TCZ(); ++k) {
                  for(unsigned int j = 0; j < CC::TCY(); ++j ) {
                     for(unsigned int i = 0; i < CC::TCX(); ++i) {
                        simulation_data[simulation_data_counter++] = 0.0;
                     }
                  }
               }
            }
         }
      }
      simulation_data_counter = 0;

      // Write the data into the hdf5 file
      dataset_simulation = H5Dcreate2(group_simulation,dataset_name.c_str(),H5T_NATIVE_DOUBLE,filespace_simulation,H5P_DEFAULT,property_list_create_dataset,H5P_DEFAULT);
      memspace_simulation = H5Screate_simple(1,simulation_write_count.data(),NULL);
      slap_simulation = H5Dget_space(dataset_simulation);
      H5Sselect_hyperslab(slap_simulation,H5S_SELECT_SET,simulation_start.data(),0,simulation_write_count.data(),NULL);
      H5Dwrite(dataset_simulation,H5T_NATIVE_DOUBLE,memspace_simulation,slap_simulation,property_list_write_dataset,simulation_data.data());
      H5Sclose(memspace_simulation);
      H5Sclose(slap_simulation);
      H5Dclose(dataset_simulation);
   }

   /** Writing partition*/
   dataset_name = "partition";
   double const rank = double( MpiUtilities::MyRankId() );
   for(auto const& level : tree_.FullNodeList()) {
      for(unsigned int n=0; n < level.size(); ++n) {
         for(unsigned int k = 0; k < CC::TCZ(); ++k) {
            for(unsigned int j = 0; j < CC::TCY(); ++j ) {
               for(unsigned int i = 0; i < CC::TCX(); ++i) {
                  simulation_data[simulation_data_counter++] = rank;
               }
            }
         }
      } //nodes
   } //level

   simulation_data_counter = 0;

   // Write the data into the hdf5 file
   dataset_simulation = H5Dcreate2(group_simulation,dataset_name.c_str(),H5T_NATIVE_DOUBLE,filespace_simulation,H5P_DEFAULT,property_list_create_dataset,H5P_DEFAULT);
   memspace_simulation = H5Screate_simple(1,simulation_write_count.data(),NULL);
   slap_simulation = H5Dget_space(dataset_simulation);
   H5Sselect_hyperslab(slap_simulation,H5S_SELECT_SET,simulation_start.data(),0,simulation_write_count.data(),NULL);
   H5Dwrite(dataset_simulation,H5T_NATIVE_DOUBLE,memspace_simulation,slap_simulation,property_list_write_dataset,simulation_data.data());
   H5Sclose(memspace_simulation);
   H5Sclose(slap_simulation);
   H5Dclose(dataset_simulation);

   /** Closing group */
   H5Pclose(property_list_create_dataset);
   H5Pclose(property_list_write_dataset);
   H5Sclose(filespace_simulation);
   H5Gclose(group_simulation);
   /** Closing file */
   H5Pclose(property_list_create_file);
   H5Fclose(file);

   /** Write the xdmf file for the debug hdf5 file */
   xdmf_file_writer_.WriteDebugXdmfFile(debug_time,hdf_filename,number_of_nodes_globally * number_of_cells_per_block, number_of_nodes_globally* number_of_vertices_per_block);
}
