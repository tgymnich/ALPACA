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
#ifndef OUTPUT_MESH_GENERATOR_H
#define OUTPUT_MESH_GENERATOR_H

#include <hdf5.h>
#include <vector>
#include <cstdint>
#include "topology/topology_manager.h"
#include "user_specifications/compile_time_constants.h"
#include "communication/mpi_utilities.h"
#include "topology/id_information.h"

/**
 * @brief Helper functions for the mesh generation.
 */
namespace {

   constexpr unsigned int NumberOfVerticesPerBlock() { return ( CC::ICX() + 1 ) * ( CC::ICY() + 1 ) * ( CC::ICZ() + 1 ); }
   constexpr unsigned int NumberOfCellsPerBlock() { return CC::ICX()*CC::ICY()*CC::ICZ(); }

   /**
    * @brief Gives the number of vertices in each dimension for a given maximum level.
    * @param maximum_level Level for which the number of vertices should be obtained.
    * @param level_zero_nodes_xyz Number of nodes specified on level zero (per dimension).
    * @return Number of vertices in each dimension.
    */
   std::array<unsigned long long int, 3> VertexCountPerDimension( unsigned int const maximum_level,
                                                                  std::array<unsigned int, 3> const level_zero_nodes_xyz ) {
      unsigned long long int resolution = 1 << maximum_level;
      unsigned long long int vertex_count_x = level_zero_nodes_xyz[0] * resolution * CC::ICX() + 1;
      unsigned long long int vertex_count_y = level_zero_nodes_xyz[1] * resolution * CC::ICY() + 1;
      unsigned long long int vertex_count_z = level_zero_nodes_xyz[2] * resolution * CC::ICZ() + 1;
      return {vertex_count_x, vertex_count_y, vertex_count_z};
   }

   /**
    * @brief Total number of all vertices (x + y + z-direction) for a given maximum level.
    * @param maximum_level Level for which the number of vertices should be obtained.
    * @param level_zero_nodes_xyz Number of nodes specified on level zero (per dimension).
    * @return Total number of vertices summed for all dimensions.
    */
   unsigned long long int TotalVertexCountOnFinestLevel( unsigned int const maximum_level,
                                                         std::array<unsigned int, 3> const level_zero_nodes_xyz ) {
      auto&& count_per_dimension = VertexCountPerDimension( maximum_level, level_zero_nodes_xyz );
      return count_per_dimension[0] * count_per_dimension[1] * count_per_dimension[2];
   }

   /**
    * @brief Gives the local (rank-wise) start index and number of elements that should be written in the overall (output) vector for
    *        vertex IDs and coordinates.
    * @param total_count Total number of elements that should be written (summed for all ranks).
    * @param rank Rank identifier.
    * @param number_of_ranks Number of ranks used for the simulation.
    * @return Rank-wise number of elements and start index.
    */
   std::array<unsigned long long int, 2> GetLocalCountAndStart( unsigned long long int const total_count, int const rank,
                                                                int const number_of_ranks ) {
      int const remainder = total_count % number_of_ranks;
      unsigned long long int local_count = total_count / number_of_ranks;
      unsigned long long int const start_index = local_count * rank + ( rank < remainder ? rank : remainder );
      if(  rank  < remainder ) { local_count++; }
      return {local_count, start_index};
   }

   /**
    * @brief Gives the cell size of a node for a given size.
    * @param node_size Size of the node.
    * @return cell size .
    */
   double CellSizeForBlockSize( double const node_size = 1.0 ) {
      // ICX() is used since always present
      return node_size / CC::ICX();
   }
}

/**
 * @brief The OutputMeshGenerator provides the functionality to compute the vertex IDs and coordinates for the standard output. Three 
 *        versions are provided, Off, Mpi and FinestLevel.
 *        Off: Here the mesh is generated with the current topology nodes on each rank separately. Therefore, it contains double created vertices 
 *             that are not filtered in a post-step. Only the coordinates for the used IDs need to be stored.
 *        Mpi: The same as for the Off version, but here a post-filtering step is done that removes double created vertices. The 
 *             coordinates for the non-used IDs are still present in the coordinates list.
 *        FinestLevel: Here, the current topology nodes are projected onto the finest level (CC::AMNL()) and the IDs are taken from the that finest level. 
 *                     Therefore, no filtering step is required. With this version an increased amount of data is stored into the output file 
 *                     since always all coordinates from the finest level need to be stored even though they are not used depending on the 
 *                     current level present for the output time step.
 * @note The mesh geenrator class is designed to provide the correct format for coordinates and IDs required for the HDF5 file format.
 */
class OutputMeshGenerator {

   TopologyManager const& topology_;

public:
   OutputMeshGenerator() = delete;
   explicit OutputMeshGenerator( TopologyManager const& topology ) : topology_( topology ) {}
   ~OutputMeshGenerator() = default;
   OutputMeshGenerator( OutputMeshGenerator const& ) = delete;
   OutputMeshGenerator& operator=( OutputMeshGenerator const& ) = delete;
   OutputMeshGenerator( OutputMeshGenerator&& ) = delete;
   OutputMeshGenerator& operator=( OutputMeshGenerator&& ) = delete;

   /**
    * @brief Gives the global dimensions (summed for all ranks) of the vertex coordinates that need to be written:
    *        In general: Number of vertices, 3 (x,y,z coordinate for each vertex).
    * @return The global dimensions for vertex coordinates.
    * @tparam VF Identifier of vertex filter type (Off, Mpi, FinestLevel).
    */
   template<VertexFilterType VF>
   std::vector<hsize_t> GlobalVertexCoordinatesDimensions() const {
      hsize_t value;
      if constexpr( VF == VertexFilterType::FinestLevel ) {
         value = TotalVertexCountOnFinestLevel( topology_.GetCurrentMaximumLevel(), topology_.GetLevelZeroBlockRatioXyz() );
      } else {
         value = std::get<1>( topology_.NodeAndLeafCount() ) * NumberOfVerticesPerBlock();
      }

      return {value, 3};
   }

   /**
    * @brief Gives the local dimensions and start index (rank-wise) of the vertex coordinates that need to be written:
    * @param number_of_leaves_locally Total number of leave nodes present on the current rank.
    * @param rank Rank identifier.
    * @param number_of_ranks Total number of ranks used in the simulation.
    * @return The start indices and local dimensions for the vertex coordinates (rank-wise).
    * @tparam VF Identifier of vertex filter type (Off, Mpi, FinestLevel)
    */
   template<VertexFilterType VF>
   std::pair<std::vector<hsize_t>, std::vector<hsize_t>> WriteStartAndCountCoordinates( [[maybe_unused]] std::size_t const number_of_leaves_locally,
                                                                                        [[maybe_unused]] int const rank = MpiUtilities::MyRankId(),
                                                                                        [[maybe_unused]] int const number_of_ranks = MpiUtilities::NumberOfRanks() ) const {
      std::vector<hsize_t> coordinates_write_count;
      std::vector<hsize_t> coordinates_write_start;

      if constexpr( VF == VertexFilterType::FinestLevel ) {
         unsigned long long int const total_vertex_count = TotalVertexCountOnFinestLevel( topology_.GetCurrentMaximumLevel(),
                                                                                          topology_.GetLevelZeroBlockRatioXyz() );

         std::array<unsigned long long int, 2>&& local_count_and_start = GetLocalCountAndStart( total_vertex_count, rank, number_of_ranks );

         coordinates_write_count = { local_count_and_start[0], 3 };
         coordinates_write_start = { local_count_and_start[1], 0 };
      } else {
         coordinates_write_count = { hsize_t( number_of_leaves_locally ) * hsize_t( NumberOfVerticesPerBlock() ), hsize_t ( 3 ) };
         coordinates_write_start = { hsize_t( topology_.LeafOffsetOfRank( rank ) ) * hsize_t ( NumberOfVerticesPerBlock() ), hsize_t ( 0 ) };
      }

      return { coordinates_write_start, coordinates_write_count };
   }

   /**
    * @brief Computes the vertex coordinates for the current rank
    * @param node_size_on_level_zero Size of one node on level zero
    * @param rank Rank identifier.
    * @param number_of_ranks Total number of ranks used in the simulation. 
    * @return Vector with all coordinates in the already correct order for HDf5 files. 
    * @tparam VF Identifier of vertex filter type (Off, Mpi, FinestLevel) 
    * 
    * @note The HDF5 file requires the coordinates to be in the order:
    *       vertex1 | vertex2 | vertex3 |     | vertexN
    *        x y z  |  x y z  |  x y z  | ... |  x y z
    * @note  The index position must align with the vertex ID specification for the hexahedrons (see function VertexIds() )

    */
   template<VertexFilterType VF>
   std::vector<double> VertexCoordinates( double const node_size_on_level_zero = 1.0, [[maybe_unused]] int const rank = MpiUtilities::MyRankId(),
                                          [[maybe_unused]] int const number_of_ranks = MpiUtilities::NumberOfRanks() ) const {

      std::vector<double> vertex_coordinates;

      if constexpr( VF == VertexFilterType::FinestLevel ) {
         unsigned int const maximum_level = topology_.GetCurrentMaximumLevel();
         unsigned long long int const total_vertex_count = TotalVertexCountOnFinestLevel( maximum_level, topology_.GetLevelZeroBlockRatioXyz() );

         std::array<unsigned long long int, 2>&& local_count_and_start = GetLocalCountAndStart( total_vertex_count, rank, number_of_ranks );

         //generate part of vertex coordinates
         unsigned long long int const resolution = 1 << maximum_level;
         vertex_coordinates = std::vector<double>(local_count_and_start[0] * 3);
         double const vertex_step_x = node_size_on_level_zero / resolution / CC::ICX();
         double const vertex_step_y = node_size_on_level_zero / resolution / CC::ICY();
         double const vertex_step_z = node_size_on_level_zero / resolution / CC::ICZ();

         //generate vertex coordinates
         std::array<unsigned long long int, 3> const vertex_counts = VertexCountPerDimension( maximum_level,
                                                                                              topology_.GetLevelZeroBlockRatioXyz() );
         for( unsigned long long int local_vertex_index = 0; local_vertex_index < local_count_and_start[0]; local_vertex_index++ ) {
            unsigned long long int  const global_vertex_index = local_vertex_index + local_count_and_start[1];
            //invert 3D to 1D array indexing; global_vertex_index = x + X * y + X * Y * z
            unsigned long int const x =   global_vertex_index %  vertex_counts[0];
            unsigned long int const y = ( global_vertex_index /  vertex_counts[0] ) % vertex_counts[1];
            unsigned long int const z =   global_vertex_index / (vertex_counts[0]   * vertex_counts[1] );
            vertex_coordinates[local_vertex_index * 3    ] = x * vertex_step_x;
            vertex_coordinates[local_vertex_index * 3 + 1] = y * vertex_step_y;
            vertex_coordinates[local_vertex_index * 3 + 2] = z * vertex_step_z;
         }
      } else {
         // We have #Nodes with #Cells, cells are build by connecting 8 vertices.
         std::vector<std::uint64_t> local_leaf_ids = topology_.LocalLeafIds();
         vertex_coordinates.resize( local_leaf_ids.size() * NumberOfVerticesPerBlock() * 3 );
         std::size_t vertex_coordinates_counter = 0;
         for( auto const& id : local_leaf_ids ) {
            double const block_size = DomainSizeOfId( id, node_size_on_level_zero );
            double const cell_size = CellSizeForBlockSize( block_size );
            std::array<double, 3> node_origin = DomainCoordinatesOfId( id, block_size );
            for( unsigned int k = 0; k <= CC::ICZ(); ++k ) {
               for( unsigned int j = 0; j <= CC::ICY(); ++j ) {
                  for( unsigned int i = 0; i <= CC::ICX(); ++i ) {
                     vertex_coordinates[vertex_coordinates_counter    ] = node_origin[0] + double(i) * cell_size;
                     vertex_coordinates[vertex_coordinates_counter + 1] = node_origin[1] + double(j) * cell_size;
                     vertex_coordinates[vertex_coordinates_counter + 2] = node_origin[2] + double(k) * cell_size;
                     vertex_coordinates_counter += 3;
                  }
               }
            }
         }
      }
      return vertex_coordinates;
   }

   // Function to generate vertex IDs
   std::vector<unsigned long long int> VertexIds( std::size_t const number_of_leaves_locally,
                                                  int const number_of_ranks = MpiUtilities::NumberOfRanks() ) const;

};

#endif // OUTPUT_MESH_GENERATOR_H
