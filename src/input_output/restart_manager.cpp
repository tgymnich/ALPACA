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
#include "restart_manager.h"

#include <hdf5.h>
#include <fstream>
#include <numeric>
#include <vector>

#include "user_specifications/compile_time_constants.h"
#include "block.h"
#include "levelset/multi_phase_manager/two_phase_manager.h"
#include "enums/interface_tag_definition.h"
#include "communication/mpi_utilities.h"

/**
 * @brief Constructs the manager for writing and reading restart snapshots.
 * @param flower The tree to read and write the simulation data. It is only modified if the simulation is restored from a snapshot.
 * @param topology The topology to get and set information about the global structure of the simulation. It is only modified if the simulation is restored from a snapshot.
 * @param setup The setup provides information about user-settings of the simulation.
 * @param output_folder Path to the output folder of the current simulation.
 */
RestartManager::RestartManager(Tree& flower, TopologyManager& topology, SimulationSetup const& setup, std::string const& output_folder) :
   tree_(flower),
   topology_(topology),
   setup_(setup),
   output_folder_(output_folder),
   logger_(LogWriter::Instance()) {
   /** Empty besides initializer list */
}

/**
 * @brief Restores the simulation topology and tree from the restart file.
 * @return The time at which the restart snapshot file was written.
 */
double RestartManager::RestoreSimulation() const {

  // prepare and open the HDF5 file
  hid_t property_list_open_file = H5Pcreate(H5P_FILE_ACCESS);
  H5Pset_fapl_mpio(property_list_open_file, MPI_COMM_WORLD, MPI_INFO_NULL);
  hid_t property_list_read_dataset = H5Pcreate(H5P_DATASET_XFER);
  H5Pset_dxpl_mpio(property_list_read_dataset, H5FD_MPIO_INDEPENDENT);

  hid_t file = H5Fopen(setup_.RestoreFileName().c_str(), H5F_ACC_RDONLY, property_list_open_file);

  // METADATA
  hid_t group = H5Gopen(file,"Metadata",H5P_DEFAULT);

  // sanity check: match metadata of HDF5 against executable
  unsigned int dimensions;
  hid_t metadata_attribute = H5Aopen(group, "Dimensions", H5P_DEFAULT);
  H5Aread(metadata_attribute, H5T_NATIVE_UINT, &dimensions);
  H5Aclose(metadata_attribute);
  if(dimensions != DTI(CC::DIM())) {
    throw std::runtime_error("Compile-time constants of the executable and restart file do not match! (dimensions)");
  }

  unsigned int internal_cells;
  metadata_attribute = H5Aopen(group, "InternalCells", H5P_DEFAULT);
  H5Aread(metadata_attribute, H5T_NATIVE_UINT, &internal_cells);
  H5Aclose(metadata_attribute);
  if(internal_cells != CC::ICX()) {
    throw std::runtime_error("Compile-time constants of the executable and restart file do not match! (internal cells)");
  }

  unsigned int halo_size;
  metadata_attribute = H5Aopen(group, "HaloSize", H5P_DEFAULT);
  H5Aread(metadata_attribute, H5T_NATIVE_UINT, &halo_size);
  H5Aclose(metadata_attribute);
  if(halo_size != CC::HS()) {
    throw std::runtime_error("Compile-time constants of the executable and restart file do not match! (halo size)");
  }

  unsigned int maximum_level_file;
  metadata_attribute = H5Aopen(group, "MaximumLevel", H5P_DEFAULT);
  H5Aread(metadata_attribute, H5T_NATIVE_UINT, &maximum_level_file);
  H5Aclose(metadata_attribute);
  if(maximum_level_file != setup_.GetMaximumLevel()) {
    throw std::runtime_error("Maximum level of input and restart file do not match!");
  }

  //we prevent physically wrong simulations by preventing non-matching reference parameters of restart and inputfile. Might be changed in the future as one could only initialize some case, write a restart file
  //and then intentionally change the reference parameters to e.g. simulate different bubble sizes. As this is an unlikely usecase, we prevent it for now.
  double length_reference;
  metadata_attribute = H5Aopen(group, "LengthReference", H5P_DEFAULT);
  H5Aread(metadata_attribute, H5T_NATIVE_DOUBLE, &length_reference);
  H5Aclose(metadata_attribute);
  if(length_reference != setup_.GetLengthReference()) {
    throw std::runtime_error("Length reference of input and restart file do not match!");
  }

  double velocity_reference;
  metadata_attribute = H5Aopen(group, "VelocityReference", H5P_DEFAULT);
  H5Aread(metadata_attribute, H5T_NATIVE_DOUBLE, &velocity_reference);
  H5Aclose(metadata_attribute);
  if(velocity_reference!= setup_.GetVelocityReference()) {
    throw std::runtime_error("Velocity reference of input and restart file do not match!");
  }

  double density_reference;
  metadata_attribute = H5Aopen(group, "DensityReference", H5P_DEFAULT);
  H5Aread(metadata_attribute, H5T_NATIVE_DOUBLE, &density_reference);
  H5Aclose(metadata_attribute);
  if(density_reference != setup_.GetDensityReference()) {
    throw std::runtime_error("Density reference of input and restart file do not match!");
  }

  double temperature_reference;
  metadata_attribute = H5Aopen(group, "TemperatureReference", H5P_DEFAULT);
  H5Aread(metadata_attribute, H5T_NATIVE_DOUBLE, &temperature_reference);
  H5Aclose(metadata_attribute);
  if(temperature_reference != setup_.GetTemperatureReference()) {
    throw std::runtime_error("Temperature reference of input and restart file do not match!");
  }

  // perhaps add more checks in the future

  // read in the time of the restart file
  double restart_time;
  metadata_attribute = H5Aopen(group, "Time", H5P_DEFAULT);
  H5Aread(metadata_attribute, H5T_NATIVE_DOUBLE, &restart_time);
  H5Aclose(metadata_attribute);

  H5Gclose(group);

  // TOPOLOGY
  // read everything needed to reconstruct the topology tree

  group = H5Gopen2(file, "Topology", H5P_DEFAULT);

  hid_t dataset = H5Dopen2(group, "Ids", H5P_DEFAULT);
  hid_t slap = H5Dget_space(dataset);
  // number of nodes equals the number of data entries in the "Ids" group
  unsigned int const node_count = H5Sget_simple_extent_npoints(slap);
  std::vector<std::uint64_t> ids(node_count);
  H5Dread(dataset, H5T_NATIVE_ULLONG, H5S_ALL, H5S_ALL, property_list_read_dataset, ids.data());
  H5Sclose(slap);
  H5Dclose(dataset);

  dataset = H5Dopen2(group, "NumberOfPhases", H5P_DEFAULT);
  slap = H5Dget_space(dataset);
  std::vector<unsigned short> number_of_phases(node_count);
  H5Dread(dataset, H5T_NATIVE_USHORT, H5S_ALL, H5S_ALL, property_list_read_dataset, number_of_phases.data());
  H5Sclose(slap);
  H5Dclose(dataset);

  dataset = H5Dopen2(group, "Materials", H5P_DEFAULT);
  slap = H5Dget_space(dataset);
  std::vector<unsigned short> materials(std::accumulate(number_of_phases.begin(), number_of_phases.end(), 0));
  H5Dread(dataset, H5T_NATIVE_USHORT, H5S_ALL, H5S_ALL, property_list_read_dataset, materials.data());
  H5Sclose(slap);
  H5Dclose(dataset);

  dataset = H5Dopen2(group, "NumberOfLevelsets", H5P_DEFAULT);
  slap = H5Dget_space(dataset);
  std::vector<unsigned short> number_of_levelsets(node_count);
  H5Dread(dataset, H5T_NATIVE_USHORT, H5S_ALL, H5S_ALL, property_list_read_dataset, number_of_levelsets.data());
  H5Sclose(slap);
  H5Dclose(dataset);

  H5Gclose(group);

  // let the topology manager restore the topology of the restart file and get the ids of the nodes I have to handle
  auto const indices_local_ids = topology_.RestoreTopology(ids, number_of_phases, materials);

  // NODEDATA
  // read the actual node-wise data, i.e. conservatives, primes, interface tags, level sets

  // hyperslap definition for conservatives+primes dataset
  constexpr unsigned int rank_conservatives_primes = 5;
  constexpr hsize_t count_conservatives_primes[rank_conservatives_primes] = {1, 1, 1, 1, 1};
  constexpr hsize_t block_conservatives_primes[5] = { 1, 1, CC::TCX(), CC::TCY(), CC::TCZ() };
  hsize_t offset_conservatives_primes[rank_conservatives_primes] = {0};

  // common hyperslap definition for levelset and interfacetags datasets
  constexpr unsigned int rank_levelsets_interfacetags = 4;
  constexpr hsize_t count_levelsets_interfacetags[rank_levelsets_interfacetags] = {1, 1, 1, 1};
  constexpr hsize_t block_levelsets_interfacetags[rank_levelsets_interfacetags] = { 1, CC::TCX(), CC::TCY(), CC::TCZ() };
  hsize_t offset_levelsets_interfacetags[rank_levelsets_interfacetags] = {0};

  // define memspaces (dimensions equal to one block of hyperslap)
  hid_t memspace_conservatives_primes = H5Screate_simple(rank_conservatives_primes,block_conservatives_primes,NULL);
  hid_t memspace_levelsets_interfacetags = H5Screate_simple(rank_levelsets_interfacetags,block_levelsets_interfacetags,NULL);

  // open resources (group, datasets, hyperslaps)
  group = H5Gopen2(file, "NodeData", H5P_DEFAULT);
  hid_t dataset_conservatives_primes = H5Dopen2(group, "ConservativesPrimes", H5P_DEFAULT);
  hid_t slap_conservatives_primes = H5Dget_space(dataset_conservatives_primes);
  hid_t dataset_levelsets = H5Dopen2(group, "Levelsets", H5P_DEFAULT);
  hid_t slap_levelsets = H5Dget_space(dataset_levelsets);
  hid_t dataset_interfacetags = H5Dopen2(group, "InterfaceTags", H5P_DEFAULT);
  hid_t slap_interfacetags = H5Dget_space(dataset_interfacetags);

  std::int8_t interface_tags[CC::TCX()][CC::TCY()][CC::TCZ()];
  std::vector<MaterialName> node_materials;

  // iterate all nodes that ended up on this rank and read in the respective data
  for(auto const index_node : indices_local_ids) {

    unsigned int const offset_block = std::accumulate(number_of_phases.begin(), number_of_phases.begin()+index_node, 0);
    // iterate all materials of this block

    node_materials.clear();
    for(unsigned int index_material = 0; index_material < number_of_phases[index_node]; ++index_material) {
      node_materials.push_back(static_cast<MaterialName>(materials[offset_block+index_material]));
    }

    std::unique_ptr<LevelsetBlock> levelset_block;
    if(number_of_levelsets[index_node] == 1) {
      offset_levelsets_interfacetags[0] = std::accumulate(number_of_levelsets.begin(), number_of_levelsets.begin()+index_node, 0);

      // read levelset
      double phi[CC::TCX()][CC::TCY()][CC::TCZ()];
      H5Sselect_hyperslab(slap_levelsets, H5S_SELECT_SET, offset_levelsets_interfacetags, NULL, count_levelsets_interfacetags, block_levelsets_interfacetags);
      H5Dread(dataset_levelsets, H5T_NATIVE_DOUBLE, memspace_levelsets_interfacetags, slap_levelsets, property_list_read_dataset, phi);
      levelset_block = std::make_unique<LevelsetBlock>(phi);

      // read interface tags
      H5Sselect_hyperslab(slap_interfacetags, H5S_SELECT_SET, offset_levelsets_interfacetags, NULL, count_levelsets_interfacetags, block_levelsets_interfacetags);
      H5Dread(dataset_interfacetags, H5T_NATIVE_CHAR, memspace_levelsets_interfacetags, slap_interfacetags, property_list_read_dataset, interface_tags);
    } else {
      if(number_of_phases[index_node] == 1) {
        // for a single-phase node the interface tags are uniform
        std::int8_t const uniform_tag = MaterialSignCapsule::SignOfMaterial(node_materials.front())*ITTI(IT::BulkPhase);
        for(unsigned int i = 0; i < CC::TCX(); ++i){
          for(unsigned int j = 0; j < CC::TCY(); ++j){
            for(unsigned int k = 0; k < CC::TCZ(); ++k){
              interface_tags[i][j][k] = uniform_tag;
            }
          }
        }
      }
    }

    Node& new_node = tree_.CreateNode(ids[index_node],node_materials,interface_tags,std::move(levelset_block));

    for(unsigned int index_material = 0; index_material < number_of_phases[index_node]; ++index_material) {
      MaterialName const material = static_cast<MaterialName>(materials[offset_block+index_material]);
      Block& block = new_node.GetPhaseByMaterial(material);

      offset_conservatives_primes[0] = offset_block + index_material;
      for(Equation const& equation : FF::ASOE()) {
        // read conservatives
        offset_conservatives_primes[1] = ETI(equation);
        H5Sselect_hyperslab(slap_conservatives_primes, H5S_SELECT_SET, offset_conservatives_primes, NULL, count_conservatives_primes, block_conservatives_primes);
        H5Dread(dataset_conservatives_primes, H5T_NATIVE_DOUBLE, memspace_conservatives_primes, slap_conservatives_primes, property_list_read_dataset, block.GetRightHandSideBuffer(equation));
      }
      for(PrimeState const& prime : FF::ASOP()) {
        // read prime state
        offset_conservatives_primes[1] = FF::ANOE() + PTI(prime);
        H5Sselect_hyperslab(slap_conservatives_primes, H5S_SELECT_SET, offset_conservatives_primes, NULL, count_conservatives_primes, block_conservatives_primes);
        H5Dread(dataset_conservatives_primes, H5T_NATIVE_DOUBLE, memspace_conservatives_primes, slap_conservatives_primes, property_list_read_dataset, block.GetPrimeStateBuffer(prime));
      }
    }
  }

  // close resources
  H5Sclose(slap_conservatives_primes);
  H5Sclose(slap_levelsets);
  H5Sclose(slap_interfacetags);
  H5Dclose(dataset_conservatives_primes);
  H5Dclose(dataset_levelsets);
  H5Dclose(dataset_interfacetags);
  H5Sclose(memspace_conservatives_primes);
  H5Sclose(memspace_levelsets_interfacetags);
  H5Gclose(group);

  // FINISH FILE
  H5Pclose(property_list_open_file);
  H5Pclose(property_list_read_dataset);
  H5Fclose(file);


  // return the time of the restart file
  return restart_time;
}

/**
 * @brief Writes a restart file containing the topology and all relevant node data at the current time step.
 * @param timestep The current time step which is written to the restart snapshot file.
 */
std::string RestartManager::WriteSnapshotFile(double const timestep) const {

  double const dimensionalized_time = setup_.DimensionalizeTime(timestep)*setup_.GetTimeNamingFactor();
  std::string const file_name = output_folder_ + RestartSubfolderName() + "/restart_" + std::to_string(dimensionalized_time) + ".h5";
  //Specify the different datasets and groups that are written
  std::string const group_name_nodedata = "NodeData";
  std::string const group_name_topology = "Topology";
  std::string const dataset_name_ids = "Ids";
  std::string const dataset_name_phasecount = "NumberOfPhases";
  std::string const dataset_name_materials = "Materials";
  std::string const dataset_name_levelsetcount = "NumberOfLevelsets";
  std::string const dataset_name_conservatives_primes = "ConservativesPrimes";
  std::string const dataset_name_levelsets = "Levelsets";
  std::string const dataset_name_interfacetags = "InterfaceTags";

  // NH: The HDF5 C Interface is used as it is better documented (more readable) and the h5c++ compiler is not configured in the institute
  hid_t property_list_create_file = H5Pcreate(H5P_FILE_ACCESS);
  H5Pset_fapl_mpio(property_list_create_file, MPI_COMM_WORLD, MPI_INFO_NULL);

  hid_t file = H5Fcreate(file_name.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, property_list_create_file);
  hid_t group_nodedata = H5Gcreate2(file, group_name_nodedata.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  hid_t group_topology = H5Gcreate2(file, group_name_topology.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  // determine global sizes and local offsets of the data

  std::vector<std::pair<unsigned int, unsigned int>> const nodes_blocks_per_rank = topology_.NodesAndBlocksPerRank();
  std::pair<unsigned int, unsigned int> const nodes_blocks_global = topology_.NodeAndBlockCount();

  unsigned int const my_rank = MpiUtilities::MyRankId(); //implicit cast intentionally.

  unsigned int const local_nodes_offset = std::accumulate(nodes_blocks_per_rank.begin(), nodes_blocks_per_rank.begin() + my_rank, 0u,
     [](unsigned int const& sum, std::pair<unsigned int, unsigned int> const& entry){return sum + entry.first;});
  unsigned int const global_nodes_count = nodes_blocks_global.first;

  unsigned int const local_blocks_offset = std::accumulate(nodes_blocks_per_rank.begin(), nodes_blocks_per_rank.begin() + my_rank, 0u,
     [](unsigned int const& sum, std::pair<unsigned int, unsigned int> const& entry){return sum + entry.second;});
  unsigned int const global_blocks_count = nodes_blocks_global.second;

  // count the local levelsets in the Lmax nodes. So far nodes can have only one levelset. So this way of counting is fine.
  unsigned int const local_levelsets_count = tree_.NodesWithLevelset().size();

  std::vector<unsigned int> levelsets_per_rank( MpiUtilities::NumberOfRanks() );
  MPI_Allgather(&local_levelsets_count, 1, MPI_UNSIGNED, levelsets_per_rank.data(), 1, MPI_UNSIGNED, MPI_COMM_WORLD);

  unsigned int const local_levelsets_offset = std::accumulate(levelsets_per_rank.begin(), levelsets_per_rank.begin() + my_rank, 0u);
  unsigned int const global_levelsets_count = std::accumulate(levelsets_per_rank.begin() + my_rank, levelsets_per_rank.end(), local_levelsets_offset);


  // rank of the different data sets in the HDF5 file. "Rank" is meant analogous to arrays, hence the number of indices
  // ********************************************************************************************
  // data sets with one entry per block or node
  constexpr unsigned int rank_nodes_blocks = 1; // dimension is: node/block
  // data sets containing simulation data (conservative+prime states) per cell per node
  constexpr unsigned int rank_conservatives_primes = 5; // dimensions are: node, equation/prime, x, y, z
  // data sets containing one entry (phi value or interface tag) per cell per node
  constexpr unsigned int rank_levelsets_interfacetags = 4; // dimensions are: node, x, y, z

  // corresponding dimensions for each data rank
  hsize_t const global_dimensions_nodes[rank_nodes_blocks] = { global_nodes_count }; //We have N nodes so N entries.
  hsize_t const global_dimensions_blocks[rank_nodes_blocks] = { global_blocks_count };
  hsize_t const global_dimensions_conservatives_primes[rank_conservatives_primes] = { global_blocks_count, FF::ANOE()+FF::ANOP(), CC::TCX(), CC::TCY(), CC::TCZ() };
  hsize_t const global_dimensions_levelsets_interfacetags[rank_levelsets_interfacetags] = { global_levelsets_count, CC::TCX(), CC::TCY(), CC::TCZ() };

  // HDF5 hyperslap dimensions
  // * hyperslap chunk = one whole node/block of simulation data)
  // * hyperslap block = one data entry per simulation cell (e.g. one conservative variable for each cell)
  // * hyperslap count = number of hyperslap blocks to read/write per call (always one in our case)
  // * hyperslap offset = at which hyperslap block to start reading/writing
  // * hyperslap stride = how many elements to move in each dimension to select the next block (not used in our case,
  // *                    thus set to NULL which defaults to stride 1 for every dimension)
  // ***************************************************************************

  // for datasets that have one entry per node/block (in Topology group)
  constexpr hsize_t chunk_nodes_blocks[rank_nodes_blocks] = { 1 }; // one entry per node/block
  constexpr hsize_t count_nodes_blocks[rank_nodes_blocks] = { 1 };
  constexpr hsize_t block_nodes_blocks[rank_nodes_blocks] = { 1 };
  hsize_t  offset_nodes[rank_nodes_blocks] = { local_nodes_offset };
  hsize_t offset_blocks[rank_nodes_blocks] = { local_blocks_offset };

  // for datasets that have N entries per block and simulation cell where N is the number of equations + prime states
  constexpr hsize_t  chunk_conservatives_primes[rank_conservatives_primes] = { 1, FF::ANOE()+FF::ANOP(), CC::TCX(), CC::TCY(), CC::TCZ() }; // only one block (first index)
  constexpr hsize_t  count_conservatives_primes[rank_conservatives_primes] = { 1, 1, 1, 1, 1 };
  constexpr hsize_t  block_conservatives_primes[rank_conservatives_primes] = { 1, 1, CC::TCX(), CC::TCY(), CC::TCZ() }; // only one equation
  hsize_t offset_conservatives_primes[rank_conservatives_primes] = { local_blocks_offset, 0, 0, 0, 0 };

  // for datasets that have one entry per levelset-block and simulation cell
  constexpr hsize_t  chunk_levelsets_interfacetags[rank_levelsets_interfacetags] = { 1, CC::TCX(), CC::TCY(), CC::TCZ() }; // only one levelset block
  constexpr hsize_t  count_levelsets_interfacetags[rank_levelsets_interfacetags] = { 1, 1, 1, 1 };
  constexpr hsize_t  block_levelsets_interfacetags[rank_levelsets_interfacetags] = { 1, CC::TCX(), CC::TCY(), CC::TCZ() };
  hsize_t offset_levelsets_interfacetags[rank_levelsets_interfacetags] = { local_levelsets_offset, 0, 0, 0 };


  hid_t property_list_create_dataset_nodes_blocks = H5Pcreate(H5P_DATASET_CREATE);
  hid_t property_list_create_dataset_conservatives_primes = H5Pcreate(H5P_DATASET_CREATE);
  hid_t property_list_create_dataset_levelsets_interfacetags = H5Pcreate(H5P_DATASET_CREATE);

  H5Pset_chunk(property_list_create_dataset_nodes_blocks, rank_nodes_blocks, chunk_nodes_blocks);
  H5Pset_chunk(property_list_create_dataset_conservatives_primes, rank_conservatives_primes, chunk_conservatives_primes);
  H5Pset_chunk(property_list_create_dataset_levelsets_interfacetags, rank_levelsets_interfacetags, chunk_levelsets_interfacetags);

  hid_t filespace_ids = H5Screate_simple(rank_nodes_blocks, global_dimensions_nodes, NULL);
  hid_t filespace_phasecount = H5Screate_simple(rank_nodes_blocks, global_dimensions_nodes, NULL);
  hid_t filespace_materials = H5Screate_simple(rank_nodes_blocks, global_dimensions_blocks, NULL);
  hid_t filespace_levelsetcount = H5Screate_simple(rank_nodes_blocks, global_dimensions_nodes, NULL);
  hid_t filespace_conservatives_primes = H5Screate_simple(rank_conservatives_primes, global_dimensions_conservatives_primes, NULL);
  hid_t filespace_levelsets = H5Screate_simple(rank_levelsets_interfacetags, global_dimensions_levelsets_interfacetags, NULL);
  hid_t filespace_interfacetags = H5Screate_simple(rank_levelsets_interfacetags, global_dimensions_levelsets_interfacetags, NULL);

  hid_t dataset_ids = H5Dcreate2(group_topology,dataset_name_ids.c_str(), H5T_NATIVE_ULLONG, filespace_ids, H5P_DEFAULT, property_list_create_dataset_nodes_blocks, H5P_DEFAULT);
  hid_t dataset_phasecount = H5Dcreate2(group_topology,dataset_name_phasecount.c_str(), H5T_NATIVE_USHORT, filespace_phasecount, H5P_DEFAULT, property_list_create_dataset_nodes_blocks, H5P_DEFAULT);
  hid_t dataset_materials = H5Dcreate2(group_topology,dataset_name_materials.c_str(), H5T_NATIVE_USHORT, filespace_materials, H5P_DEFAULT, property_list_create_dataset_nodes_blocks, H5P_DEFAULT);
  hid_t dataset_levelsetcount = H5Dcreate2(group_topology,dataset_name_levelsetcount.c_str(), H5T_NATIVE_USHORT, filespace_levelsetcount, H5P_DEFAULT, property_list_create_dataset_nodes_blocks, H5P_DEFAULT);
  hid_t dataset_conservatives_primes = H5Dcreate2(group_nodedata,dataset_name_conservatives_primes.c_str(), H5T_NATIVE_DOUBLE,filespace_conservatives_primes, H5P_DEFAULT, property_list_create_dataset_conservatives_primes, H5P_DEFAULT);
  hid_t dataset_levelsets = H5Dcreate2(group_nodedata,dataset_name_levelsets.c_str(), H5T_NATIVE_DOUBLE, filespace_levelsets, H5P_DEFAULT, property_list_create_dataset_levelsets_interfacetags, H5P_DEFAULT);
  hid_t dataset_interfacetags = H5Dcreate2(group_nodedata,dataset_name_interfacetags.c_str(), H5T_NATIVE_CHAR, filespace_interfacetags, H5P_DEFAULT, property_list_create_dataset_levelsets_interfacetags, H5P_DEFAULT);

  hid_t property_list_write_dataset = H5Pcreate(H5P_DATASET_XFER);
  H5Pset_dxpl_mpio(property_list_write_dataset, H5FD_MPIO_INDEPENDENT);

  // create memspaces (dimensions are equal to the dimensions of a hyperslap block)
  hid_t memspace_nodes_blocks = H5Screate_simple(rank_nodes_blocks, block_nodes_blocks, NULL);
  hid_t memspace_conservatives_primes = H5Screate_simple(rank_conservatives_primes, block_conservatives_primes, NULL);
  hid_t memspace_levelsets_interfacetags = H5Screate_simple(rank_levelsets_interfacetags, block_levelsets_interfacetags, NULL);

  hid_t slap_ids = H5Dget_space(dataset_ids);
  hid_t slap_phasecount = H5Dget_space(dataset_phasecount);
  hid_t slap_materials = H5Dget_space(dataset_materials);
  hid_t slap_levelsetcount = H5Dget_space(dataset_levelsetcount);
  hid_t slap_conservatives_primes = H5Dget_space(dataset_conservatives_primes);
  hid_t slap_levelsets = H5Dget_space(dataset_levelsets);
  hid_t slap_interfacetags = H5Dget_space(dataset_interfacetags);

  for(auto const& level : tree_.FullNodeList()) {
    for(auto const& [id, node]: level) {
      unsigned short const number_of_phases = node.GetPhases().size();
      // for now we only have to deal with one or none levelset
      unsigned short const number_of_levelsets = node.HasLevelset() ? 1 : 0;

      // write data that is one single entry per node
      // write Id
      H5Sselect_hyperslab(slap_ids, H5S_SELECT_SET, offset_nodes, NULL, count_nodes_blocks, block_nodes_blocks);
      H5Dwrite(dataset_ids, H5T_NATIVE_ULLONG, memspace_nodes_blocks, slap_ids, property_list_write_dataset, &id);
      // write NumberOfPhases
      H5Sselect_hyperslab(slap_phasecount, H5S_SELECT_SET, offset_nodes, NULL, count_nodes_blocks, block_nodes_blocks);
      H5Dwrite(dataset_phasecount, H5T_NATIVE_USHORT, memspace_nodes_blocks, slap_phasecount, property_list_write_dataset, &number_of_phases);
      // write NumberOfLevelsets
      H5Sselect_hyperslab(slap_levelsetcount, H5S_SELECT_SET, offset_nodes, NULL, count_nodes_blocks, block_nodes_blocks);
      H5Dwrite(dataset_levelsetcount, H5T_NATIVE_USHORT, memspace_nodes_blocks, slap_levelsetcount, property_list_write_dataset, &number_of_levelsets);

      // all single entry per node data is written, increase offset:
      offset_nodes[0]++;

      // loop over materials
      for(auto const& mat_block : node.GetPhases()) {
        unsigned short const material = MTI(mat_block.first);

        // write materials (one entry per block)
        H5Sselect_hyperslab(slap_materials, H5S_SELECT_SET, offset_blocks, NULL, count_nodes_blocks, block_nodes_blocks);
        H5Dwrite(dataset_materials, H5T_NATIVE_USHORT, memspace_nodes_blocks, slap_materials, property_list_write_dataset, &material);
        offset_blocks[0]++;

        // write conservative and prime state data
        for(Equation const& equation : FF::ASOE()) {
          offset_conservatives_primes[1] = ETI(equation);
          // write conservative variable
          H5Sselect_hyperslab(slap_conservatives_primes, H5S_SELECT_SET, offset_conservatives_primes, NULL, count_conservatives_primes, block_conservatives_primes);
          H5Dwrite(dataset_conservatives_primes, H5T_NATIVE_DOUBLE, memspace_conservatives_primes, slap_conservatives_primes, property_list_write_dataset, mat_block.second.GetAverageBuffer(equation));
        }
        for(PrimeState const& prime : FF::ASOP()) {
          offset_conservatives_primes[1] = FF::ANOE() + PTI(prime);
          // write prime state variable
          H5Sselect_hyperslab(slap_conservatives_primes, H5S_SELECT_SET, offset_conservatives_primes, NULL, count_conservatives_primes, block_conservatives_primes);
          H5Dwrite(dataset_conservatives_primes, H5T_NATIVE_DOUBLE, memspace_conservatives_primes, slap_conservatives_primes, property_list_write_dataset, mat_block.second.GetPrimeStateBuffer(prime));
        }
        offset_conservatives_primes[0]++; // next block

      } // end material loop

      if(node.HasLevelset()) {
        // write levelset(s)
        H5Sselect_hyperslab(slap_levelsets, H5S_SELECT_SET, offset_levelsets_interfacetags, NULL, count_levelsets_interfacetags, block_levelsets_interfacetags);
        H5Dwrite(dataset_levelsets, H5T_NATIVE_DOUBLE, memspace_levelsets_interfacetags, slap_levelsets, property_list_write_dataset, node.GetLevelsetBlock().GetPhi());
        // write interface tags
        H5Sselect_hyperslab(slap_interfacetags, H5S_SELECT_SET, offset_levelsets_interfacetags, NULL, count_levelsets_interfacetags, block_levelsets_interfacetags);
        H5Dwrite(dataset_interfacetags, H5T_NATIVE_CHAR, memspace_levelsets_interfacetags, slap_interfacetags, property_list_write_dataset, node.GetInterfaceTags());

        offset_levelsets_interfacetags[0]++; // next levelset-block
      }
    }
  }

  // close resources
  H5Pclose(property_list_create_dataset_conservatives_primes);
  H5Pclose(property_list_create_dataset_levelsets_interfacetags);
  H5Pclose(property_list_create_dataset_nodes_blocks);
  H5Pclose(property_list_write_dataset);
  H5Sclose(slap_conservatives_primes);
  H5Sclose(slap_levelsets);
  H5Sclose(slap_interfacetags);
  H5Sclose(slap_ids);
  H5Sclose(slap_phasecount);
  H5Sclose(slap_materials);
  H5Sclose(slap_levelsetcount);
  H5Sclose(filespace_conservatives_primes);
  H5Sclose(filespace_levelsets);
  H5Sclose(filespace_interfacetags);
  H5Sclose(filespace_ids);
  H5Sclose(filespace_phasecount);
  H5Sclose(filespace_materials);
  H5Sclose(filespace_levelsetcount);
  H5Sclose(memspace_nodes_blocks);
  H5Sclose(memspace_conservatives_primes);
  H5Sclose(memspace_levelsets_interfacetags);
  H5Dclose(dataset_conservatives_primes);
  H5Dclose(dataset_levelsets);
  H5Dclose(dataset_interfacetags);
  H5Dclose(dataset_ids);
  H5Dclose(dataset_phasecount);
  H5Dclose(dataset_materials);
  H5Dclose(dataset_levelsetcount);
  H5Gclose(group_nodedata);
  H5Gclose(group_topology);

  /* Writing Meta data:
   * Time of restart file
   * Number of Dimensions
   * Number of Internal Cells
   * Halo Size
   * Maximum Level

   // perhaps the following
   // * CFL Number
   // * Number of Fluids
   // * MR Epsilon
   // * Lmax
   // * Gravity
   // * Boundary Types
   // * Material Parameters
   */
  hid_t metadata_group = H5Gcreate2(file, "Metadata", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  // Time
  hid_t metadata_attribute_type = H5Screate(H5S_SCALAR);
  hid_t metadata_attribute = H5Acreate2(metadata_group, "Time", H5T_NATIVE_DOUBLE, metadata_attribute_type, H5P_DEFAULT, H5P_DEFAULT);
  H5Awrite(metadata_attribute, H5T_NATIVE_DOUBLE, &timestep);
  H5Sclose(metadata_attribute_type);
  H5Aclose(metadata_attribute);
  // Number of Dimensions
  metadata_attribute_type = H5Screate(H5S_SCALAR);
  metadata_attribute = H5Acreate2(metadata_group, "Dimensions", H5T_NATIVE_USHORT, metadata_attribute_type, H5P_DEFAULT, H5P_DEFAULT);
  constexpr unsigned int number_of_dimensions = DTI(CC::DIM());
  H5Awrite(metadata_attribute, H5T_NATIVE_UINT, &number_of_dimensions);
  H5Sclose(metadata_attribute_type);
  H5Aclose(metadata_attribute);
  // Number of Internal Cells
  metadata_attribute_type = H5Screate(H5S_SCALAR);
  metadata_attribute = H5Acreate2(metadata_group, "InternalCells", H5T_NATIVE_UINT, metadata_attribute_type, H5P_DEFAULT, H5P_DEFAULT);
  constexpr unsigned int number_of_cells = CC::ICX(); //X is always "full"
  H5Awrite(metadata_attribute, H5T_NATIVE_UINT, &number_of_cells);
  H5Sclose(metadata_attribute_type);
  H5Aclose(metadata_attribute);
  // Number of Halo Cells
  metadata_attribute_type = H5Screate(H5S_SCALAR);
  metadata_attribute = H5Acreate2(metadata_group, "HaloSize", H5T_NATIVE_UINT, metadata_attribute_type, H5P_DEFAULT, H5P_DEFAULT);
  constexpr unsigned int halo_size = CC::HS();
  H5Awrite(metadata_attribute, H5T_NATIVE_UINT, &halo_size);
  H5Sclose(metadata_attribute_type);
  H5Aclose(metadata_attribute);
  // Maximum level
  metadata_attribute_type = H5Screate(H5S_SCALAR);
  metadata_attribute = H5Acreate2(metadata_group, "MaximumLevel", H5T_NATIVE_UINT, metadata_attribute_type, H5P_DEFAULT, H5P_DEFAULT);
  unsigned int const maximum_level = setup_.GetMaximumLevel();
  H5Awrite(metadata_attribute, H5T_NATIVE_UINT, &maximum_level);
  H5Sclose(metadata_attribute_type);
  H5Aclose(metadata_attribute);
  // Length reference
  metadata_attribute_type = H5Screate(H5S_SCALAR);
  metadata_attribute = H5Acreate2(metadata_group, "LengthReference", H5T_NATIVE_DOUBLE, metadata_attribute_type, H5P_DEFAULT, H5P_DEFAULT);
  double const length_reference = setup_.GetLengthReference();
  H5Awrite(metadata_attribute, H5T_NATIVE_DOUBLE, &length_reference);
  H5Sclose(metadata_attribute_type);
  H5Aclose(metadata_attribute);
  // Velocity reference
  metadata_attribute_type = H5Screate(H5S_SCALAR);
  metadata_attribute = H5Acreate2(metadata_group, "VelocityReference", H5T_NATIVE_DOUBLE, metadata_attribute_type, H5P_DEFAULT, H5P_DEFAULT);
  double const velocity_reference = setup_.GetVelocityReference();
  H5Awrite(metadata_attribute, H5T_NATIVE_DOUBLE, &velocity_reference);
  H5Sclose(metadata_attribute_type);
  H5Aclose(metadata_attribute);
  // Density reference
  metadata_attribute_type = H5Screate(H5S_SCALAR);
  metadata_attribute = H5Acreate2(metadata_group, "DensityReference", H5T_NATIVE_DOUBLE, metadata_attribute_type, H5P_DEFAULT, H5P_DEFAULT);
  double const density_reference = setup_.GetDensityReference();
  H5Awrite(metadata_attribute, H5T_NATIVE_DOUBLE, &density_reference);
  H5Sclose(metadata_attribute_type);
  H5Aclose(metadata_attribute);
  // Temperature reference
  metadata_attribute_type = H5Screate(H5S_SCALAR);
  metadata_attribute = H5Acreate2(metadata_group, "TemperatureReference", H5T_NATIVE_DOUBLE, metadata_attribute_type, H5P_DEFAULT, H5P_DEFAULT);
  double const temperature_reference = setup_.GetTemperatureReference();
  H5Awrite(metadata_attribute, H5T_NATIVE_DOUBLE, &temperature_reference);
  H5Sclose(metadata_attribute_type);
  H5Aclose(metadata_attribute);
  // Add other variables here (e.g. number of equations) if a more throrough metadata check is prefered

  H5Gclose(metadata_group);

  H5Pclose(property_list_create_file);
  H5Fclose(file);

  logger_.LogMessage("Restart File " + file_name + " written");
  return file_name;
}
