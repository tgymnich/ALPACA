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
#include <functional>
#include "output_mesh_generator.h"

/**
 * @brief Declare all functions that are required for filtering the double created vertices with an Mpi routine.
 * @note Only required for VertexFilterType::Mpi
 */
namespace VertexFilter {

   enum class Direction {X, Y, Z};

   /**
    * @brief Initializes and sorts lists of send/receive pairs of leaves for each rank.
    * @param send_list Output list for sending leaves.
    * @param recv_list Output list for receiving leaves.
    * @param local_leaf_ids List of all leaves of calling MPI rank.
    * @param direction Axis along which to send/receive.
    * @param topology_ The topology to get information about the MPI status of the simulation.
    */
    void InitializeSendRecvRankLists(std::vector<std::vector<std::pair<std::uint64_t, std::uint64_t>>>& send_list,
                                     std::vector<std::vector<std::pair<std::uint64_t, std::uint64_t>>>& recv_list,
                                     std::vector<std::uint64_t> const& local_leaf_ids,
                                     Direction const direction,
                                     TopologyManager const& topology_) {

      int max_rank = 0;
      MPI_Comm_size(MPI_COMM_WORLD, &max_rank);

      //init array with empty vectors
      for(int r = 0; r < max_rank; r++) {
         send_list[r] = std::vector<std::pair<std::uint64_t, std::uint64_t>>();
         recv_list[r] = std::vector<std::pair<std::uint64_t, std::uint64_t>>();
      }

      //init send recv direction
      BoundaryLocation send_direction;
      BoundaryLocation recv_direction;
      switch( direction ) {
         case Direction::X: {
            send_direction = BoundaryLocation::West;
            recv_direction = BoundaryLocation::East;
         }
         break;
         case Direction::Y: {
            send_direction = BoundaryLocation::South;
            recv_direction = BoundaryLocation::North;
         }
         break;
         default: /* Direction::Z */ {
            send_direction = BoundaryLocation::Bottom;
            recv_direction = BoundaryLocation::Top;
         }
      }

      //iterate over all leaves to get target ranks
      for( std::uint64_t const& node_id: local_leaf_ids ) {
         //set up send list
         std::vector<std::uint64_t> send_id_list;
         send_id_list = topology_.GetNeighboringLeaves(node_id, send_direction);
         //add every neighbor
         for(auto const& target_id: send_id_list) {
            int const target_rank = topology_.GetRankOfNode(target_id);
            std::pair<std::uint64_t, std::uint64_t> target_source_id(target_id, node_id);
            send_list[target_rank].push_back(target_source_id);
         }
         //set up recv list
         std::vector<uint64_t> recv_id_list;
         recv_id_list = topology_.GetNeighboringLeaves(node_id, recv_direction);
         //add every neighbor
         for(auto const& source_id: recv_id_list) {
            int const source_rank = topology_.GetRankOfNode(source_id);
            std::pair<std::uint64_t, std::uint64_t> target_source_id(node_id, source_id);
            recv_list[source_rank].push_back(target_source_id);
         }
      }
      //sort ids of each rank -> can use position of rank as mpi tag
      for(int r = 0; r < max_rank; r++) {
         std::sort(send_list[r].begin(), send_list[r].end());
         std::sort(recv_list[r].begin(), recv_list[r].end());
      }
      return;
   }

   /**
    * @brief Calculate two dimensional offset of topology node with regards to root at specific tree depth.
    * @param neighbor_id ID of node.
    * @param neighbor_level_diff Number of depth levels to consider.
    * @param direction Axis defining two dimensional offset in 3D space.
    * @return A 2D-Tuple corresponding to the 2D offset.
    */
   std::pair<unsigned int, unsigned int> GetNeighborBlockOffset( std::uint64_t const neighbor_id, int const neighbor_level_diff, Direction const direction ) {
      unsigned int node_offset_1 = 0;
      unsigned int node_offset_2 = 0;
      int power_two = 1;
      //set up location function dependend on direction
      std::function<bool(std::uint64_t)> sibling_func_1;
      std::function<bool(std::uint64_t)> sibling_func_2;
      switch( direction ) {
         case Direction::X: {
            sibling_func_1 = NorthInSiblingPack;
            sibling_func_2 = TopInSiblingPack;
         }
         break;
         case Direction::Y: {
            sibling_func_1 = EastInSiblingPack;
            sibling_func_2 = TopInSiblingPack;
         }
         break;
         default: /* Direction::Z */ {
            sibling_func_1 = EastInSiblingPack;
            sibling_func_2 = NorthInSiblingPack;
         }
      }
      std::uint64_t i_node_id = neighbor_id;
      //calc offsets
      for(int l = 0; l < neighbor_level_diff; l++) {
         if(sibling_func_1(i_node_id)) node_offset_1 += power_two;
         if(sibling_func_2(i_node_id)) node_offset_2 += power_two;
         i_node_id = ParentIdOfNode(i_node_id);
         power_two <<= 1;
      }
      return std::pair<unsigned int, unsigned int>(node_offset_1, node_offset_2);
   }

   /**
    * @brief Load a face of vertex-IDs out of hexahedron topology defining array.
    * @param cell_vertex_ids Array to load IDs from.
    * @param buffer Outputbuffer for face.
    * @param id_jump_a Stride in cell_vertex_ids in first dimension.
    * @param id_jump_b Stride in cell_vertex_ids in second dimension.
    * @param size_a Number of points in first dimension.
    * @param size_b Number of points in second dimension.
    * @param center Offset of center point in one hexahedron.
    * @param edge_a Offset of point shifted into first dimension by one in one hexahedron.
    * @param edge_b Offset of point shifted into second dimension by one in one hexahedron.
    * @param corner Offset of corner point in one hexahedron.
    */
   void LoadBoundaryBuffer(unsigned long long int const* cell_vertex_ids, unsigned long long int* buffer, unsigned int const id_jump_a, unsigned int const id_jump_b,
                           unsigned int const size_a, unsigned int const  size_b, int const center, int const edge_a, int const edge_b, int const corner) {
      for(unsigned int a = 0; a < size_a; a++) {
         for(unsigned int b = 0; b < size_b; b++) {
            unsigned long long int const vertex_index = id_jump_a * a + id_jump_b * b;
            unsigned int const buffer_index = a + b * (size_a+1);
            buffer[buffer_index] = cell_vertex_ids[vertex_index + center];
            if(a == size_a-1) {
               buffer[buffer_index + 1] = cell_vertex_ids[vertex_index + edge_a];
            }
            if(b == size_b-1) {
               buffer[buffer_index + size_a + 1] = cell_vertex_ids[vertex_index + edge_b];
               if(a == size_a-1) {
                  buffer[buffer_index + size_a + 2] = cell_vertex_ids[vertex_index + corner];
               }
            }
         }
      }
   }

   /**
    * @brief Sends specified buffer to all nodes out of a list using MPI.
    * @param node_id ID of node to send from.
    * @param send_list_id_leveldiff List of node IDs paired with the tree leveldifference to node_id to send data to.
    * @param send_buffer Buffer containing data to send.
    * @param direction Axis defining direction of communication.
    * @param send_list_per_rank List of total communications of MPI rank to define a MPI tag.
    * @param mpi_send_request_list All MPI request get appended to this list.
    * @param topology The topology to get information about the MPI status of the simulation.
    */
   void SendFromNode(std::uint64_t const node_id, const std::vector<std::uint64_t>& send_id_list,
                        unsigned long long int const* send_buffer, Direction const direction,
                        std::vector<std::vector<std::pair<std::uint64_t, std::uint64_t>>> const& send_list_per_rank,
                        std::vector<MPI_Request>& mpi_send_request_list,
                        TopologyManager const& topology) {

      unsigned int cell_count_axis_1 = 1;
      unsigned int cell_count_axis_2 = 1;
      //init block size
      switch( direction ) {
         case Direction::X: {
            cell_count_axis_1 = CC::ICY();
            cell_count_axis_2 = CC::ICZ();
         }
         break;
         case Direction::Y: {
            cell_count_axis_1 = CC::ICX();
            cell_count_axis_2 = CC::ICZ();
         }
         break;
         default: /* Direction::Z */ {
            cell_count_axis_1 = CC::ICX();
            cell_count_axis_2 = CC::ICY();
         }
      }

      //iterate over sendlist
      for(auto const& target_id: send_id_list) {
         int const target_level_diff = LevelOfNode(target_id) - LevelOfNode(node_id);

         //compute send buffer size
         unsigned int const abs_level_diff = abs(target_level_diff);
         unsigned int const buf_size_1 = (cell_count_axis_1 >> abs_level_diff) + 1;
         unsigned int const buf_size_2 = (cell_count_axis_2 >> abs_level_diff) + 1;
         //offset of send buffer for target_id
         unsigned int local_send_buffer_offset = 0;

         // if buffer size is one there are whole neighbor nodes to be skipped
         std::pair<unsigned int, unsigned int> const node_offset = GetNeighborBlockOffset(target_level_diff > 0 ? target_id : node_id, abs_level_diff, direction);
         if(buf_size_1 == 1) {
            unsigned int const resolution_difference = (1 << abs_level_diff) / cell_count_axis_1;
            if(resolution_difference > 1 && node_offset.first % resolution_difference == 0) { //lower boundary aligns
               if(target_level_diff > 0) local_send_buffer_offset += node_offset.first / resolution_difference;
            } else if (resolution_difference > 1 && (node_offset.first + 1) % resolution_difference == 0) { //high boundary aligns
               if(target_level_diff < 0 && (node_offset.first + 1) != resolution_difference * cell_count_axis_1) {
                  continue;
               }
               if(target_level_diff > 0) {
                  local_send_buffer_offset += (node_offset.first + 1) / resolution_difference;
               }
               else {
                  local_send_buffer_offset += cell_count_axis_1;
               }
            } else {
               continue; //no boundary aligns -> skip send to target node
            }
         }
         if(buf_size_2==1) {
            unsigned int const resolution_difference= (1 << abs_level_diff) / cell_count_axis_2;
            if(resolution_difference > 1 && node_offset.second % resolution_difference == 0) {//lower boundary aligns
               if(target_level_diff > 0) {
                  local_send_buffer_offset += node_offset.second / resolution_difference * (cell_count_axis_1+1);
               }
            } else if (resolution_difference > 1 && (node_offset.second + 1) % resolution_difference == 0) {  //high boundary aligns
               if(target_level_diff < 0 && (node_offset.second + 1) != resolution_difference * cell_count_axis_2) {
                  continue;
               }
               if(target_level_diff > 0) {
                  local_send_buffer_offset += (node_offset.second + 1) / resolution_difference * (cell_count_axis_1 + 1);
               } else {
                  local_send_buffer_offset += (cell_count_axis_1 + 1) * cell_count_axis_2;
               }
            } else {
               continue; //no boundary aligns -> skip send to target node
            }
         }

         //create mpi type to access block in send buffer
         MPI_Datatype mpi_block_type;
         if(target_level_diff > 0) {
            //case neighbor is smaller
            local_send_buffer_offset += (buf_size_1-1) * node_offset.first + (cell_count_axis_1 + 1) * (buf_size_2 - 1) * node_offset.second; //here buffersizes are used as a scaling of node_offset
            MPI_Type_vector(buf_size_2, buf_size_1, cell_count_axis_1 + 1, MPI_UNSIGNED_LONG_LONG, &mpi_block_type);
         } else {
            //case: node is smaller or same size as neighbor
            int const stride = 1 << abs_level_diff;
            MPI_Datatype mpi_stride_type;
            //stack two vector types to get strides in two dimensions
            MPI_Type_vector(buf_size_1, 1, stride, MPI_UNSIGNED_LONG_LONG, &mpi_stride_type);
            MPI_Type_vector(buf_size_2, 1, stride, mpi_stride_type, &mpi_block_type);
            MPI_Type_free(&mpi_stride_type);
         }
         MPI_Type_commit(&mpi_block_type);

         //send buffer
         int const target_rank = topology.GetRankOfNode(target_id);
         std::pair<std::uint64_t, std::uint64_t> const target_source_id(target_id, node_id);
         std::vector<std::pair<std::uint64_t, std::uint64_t>> const send_list_target_rank = send_list_per_rank[target_rank];
         //look up rank in send_list; use its position as tag
         int const tag = std::distance(send_list_target_rank.begin(), lower_bound(send_list_target_rank.begin(), send_list_target_rank.end(), target_source_id)) % MpiUtilities::MpiTagUb();

         int rank;
         MPI_Comm_rank(MPI_COMM_WORLD, &rank);

         mpi_send_request_list.push_back(MPI_Request());
         MPI_Isend(send_buffer + local_send_buffer_offset, 1, mpi_block_type, target_rank, tag, MPI_COMM_WORLD, &mpi_send_request_list.back());
         MPI_Type_free(&mpi_block_type);
      }
   }

   /**
    * @brief Receive from a list of nodes into specified buffer.
    * @param node_id ID of receiving node.
    * @param recv_list_id_leveldiff List of node IDs paired with the tree leveldifference to node_id to receive data from.
    * @param boundary_buffer Buffer to write received data into.
    * @param direction Axis defining direction of communication.
    * @param recv_list_per_rank List of total communications of MPI rank to define a MPI tag.
    * @param topology The topology to get information about the MPI status of the simulation.
    */
   void RecvToNode(std::uint64_t const node_id, std::vector<std::uint64_t> const& recv_id_list,
                   unsigned long long int* boundary_buffer, Direction const direction,
                  std::vector<std::vector<std::pair<std::uint64_t, std::uint64_t>>> const& recv_list_per_rank,
                   TopologyManager const& topology) {

      unsigned int cell_count_axis_1 = 1;
      unsigned int cell_count_axis_2 = 1;
      switch( direction ) {
         case Direction::X: {
            cell_count_axis_1 = CC::ICY();
            cell_count_axis_2 = CC::ICZ();
         }
         break;
         case Direction::Y: {
            cell_count_axis_1 = CC::ICX();
            cell_count_axis_2 = CC::ICZ();
         }
         break;
         default: /* Direction::Z */ {
            cell_count_axis_1 = CC::ICX();
            cell_count_axis_2 = CC::ICY();
         }
      }

      for(auto const& source_id: recv_id_list) {
         //get node information
         int const source_level_diff = LevelOfNode(source_id) - LevelOfNode(node_id);
         int const source_rank = topology.GetRankOfNode(source_id);

         //compute recv buffer size
         unsigned int const abs_level_diff = abs(source_level_diff);
         unsigned int const buf_size_1 = (cell_count_axis_1 >> abs_level_diff) + 1;
         unsigned int const buf_size_2 = (cell_count_axis_2 >> abs_level_diff) + 1;

         unsigned int buffer_insert_offset = 0;

         std::pair<unsigned int, unsigned int> const node_offset = GetNeighborBlockOffset(source_level_diff > 0 ? source_id : node_id, abs_level_diff, direction);
         // if buffer size is one there are whole neighbor nodes to be skipped
         if(buf_size_1 == 1) {
            unsigned int const resolution_difference = (1 << abs_level_diff) / cell_count_axis_1;
            if(resolution_difference > 1 && node_offset.first % resolution_difference == 0) { //lower boundary aligns
               if(source_level_diff > 0) {
                  buffer_insert_offset += node_offset.first / resolution_difference;
               }
            } else if (resolution_difference > 1 && (node_offset.first + 1) % resolution_difference == 0) { //high boundary aligns
               if(source_level_diff > 0 && (node_offset.first + 1) != resolution_difference * cell_count_axis_1) {
                  continue;
               }
               if(source_level_diff > 0) {
                  buffer_insert_offset += (node_offset.first + 1) / resolution_difference;
               } else {
                  buffer_insert_offset += cell_count_axis_1;
               }
            } else {
               continue; //no boundary aligns -> skip send to target node
            }
         }
         if(buf_size_2==1) { //same as for direction 1
            unsigned int const resolution_difference= (1 << abs_level_diff) / cell_count_axis_2;
            if(resolution_difference > 1 && node_offset.second % resolution_difference == 0) {
               if(source_level_diff > 0) {
                  buffer_insert_offset += node_offset.second / resolution_difference * (cell_count_axis_1+1);
               }
            } else if (resolution_difference > 1 && (node_offset.second + 1) % resolution_difference == 0) {
               if(source_level_diff > 0 && (node_offset.second + 1) != resolution_difference * cell_count_axis_2) {
                  continue;
               }
               if(source_level_diff > 0) {
                  buffer_insert_offset += (node_offset.second + 1) / resolution_difference * (cell_count_axis_1 + 1);
               } else {
                  buffer_insert_offset += (cell_count_axis_1 + 1) * cell_count_axis_2;
               }
            } else {
               continue;
            }
         }

         //insert using MPI Types
         MPI_Datatype mpi_insert_type;
         if(source_level_diff > 0) {
            //case neighbor is smaller -> calculate buffer offset
            buffer_insert_offset += node_offset.second * (buf_size_2-1) * (cell_count_axis_1 + 1) + node_offset.first * (buf_size_1 - 1);
            MPI_Type_vector(buf_size_2, buf_size_1, cell_count_axis_1 + 1, MPI_UNSIGNED_LONG_LONG, &mpi_insert_type);
         } else {
            //case: node is smaller or same size as neighbor
            int const stride = 1 << abs_level_diff;
            MPI_Datatype mpi_stride_type;
            //stack two vector types to get strides in two dimensions
            MPI_Type_vector(buf_size_1, 1, stride, MPI_UNSIGNED_LONG_LONG, &mpi_stride_type);
            MPI_Type_vector(buf_size_2, 1, stride, mpi_stride_type, &mpi_insert_type);
            MPI_Type_free(&mpi_stride_type);
         }
         MPI_Type_commit(&mpi_insert_type);
         std::vector<std::pair<std::uint64_t, std::uint64_t>> const recv_list_source_rank = recv_list_per_rank[source_rank];
         std::pair<std::uint64_t, std::uint64_t> const target_source_id(node_id, source_id);
         //look up location in recv list for tag
         int const tag = distance(recv_list_source_rank.begin(), lower_bound(recv_list_source_rank.begin(), recv_list_source_rank.end(), target_source_id)) % MpiUtilities::MpiTagUb();
         //recv data
         int rank;
         MPI_Comm_rank(MPI_COMM_WORLD, &rank);
         MPI_Recv(boundary_buffer + buffer_insert_offset, 1, mpi_insert_type, source_rank, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
         MPI_Type_free(&mpi_insert_type);
      }
   }
} //namespace vertex_filter

/**
 * @brief Computes the vertex IDs for all hexahedrons for the current rank.
 * @param number_of_leaves_locally Number of leaves present on the current rank.
 * @param number_of_ranks Total number of ranks used in the simulation. 
 * @return Vector with all vertex IDs in the already correct order for HDf5 files. 
 * 
 * @note The HDF5 file requires the vertexIDs to be in the order (i,j,k) is west-south-bottom corner of hexahedron (order in mathermatical positive sense)
 * 
 *                                                    hex1                                               |     |                                            hexN
 *       ID_i_j_k ID_i+1_j_k ID_i+1_j+1_k ID_i_j+1_k ID_i_j_k+1 ID_i+1_j_k+1 ID_i+1_j+1_k+1 ID_i_j+1_k+1 | ... | ID_i_j_k ID_i+1_j_k ID_i+1_j+1_k ID_i_j+1_k ID_i_j_k+1 ID_i+1_j_k+1 ID_i+1_j+1_k+1 ID_i_j+1_k+1 
 * 
 * @note The coordinates for each vertex ID are taken by acessing the coordinates vector with the ID as vector index.
 */
std::vector<unsigned long long int> OutputMeshGenerator::VertexIds( std::size_t const number_of_leaves_locally, int const number_of_ranks ) const {

   std::vector<unsigned long long int> vertex_ids( number_of_leaves_locally * NumberOfCellsPerBlock() * 8 );

   if constexpr( CC::VERTEX_FILTER() == VertexFilterType::FinestLevel ) {

      unsigned int maximum_level = topology_.GetCurrentMaximumLevel();
      std::array<unsigned long long int, 3>&& vertex_count_per_dimension = VertexCountPerDimension( maximum_level, topology_.GetLevelZeroBlockRatioXyz() );
      unsigned long long int const vertex_count_x = vertex_count_per_dimension[0];
      unsigned long long int const vertex_count_y = vertex_count_per_dimension[1];

      //iterate over all leafs to generate vertex ids
      unsigned long long int vertex_id_counter = 0;
      unsigned long long int const resolution = 1 << maximum_level;
      for( std::uint64_t const& id: topology_.LocalLeafIds() ) {
         unsigned long int const point_distance = resolution >> LevelOfNode(id);
         //Find index of nodes domain by coordinates function; Cast from double to long! -> for indices greater 2^53 this is incorrect
         std::array<double, 3> const node_origin = DomainCoordinatesOfId(id, 1);
         for(unsigned int k = 0; k < CC::ICZ(); ++k) {
            for(unsigned int j = 0; j < CC::ICY(); ++j) {
               for(unsigned int i = 0; i < CC::ICX(); ++i) {
                  unsigned long int vertex_index = 0;
                  unsigned long int const index_x = (node_origin[0] * CC::ICX() + i) * point_distance;
                  unsigned long int const index_y = (node_origin[1] * CC::ICY() + j )* point_distance;
                  unsigned long int const index_z = (node_origin[2] * CC::ICZ() + k) * point_distance;

                  vertex_index = index_x + index_y * vertex_count_x + index_z * vertex_count_x * vertex_count_y;
                  vertex_ids[vertex_id_counter] = vertex_index;

                  vertex_index = (index_x + point_distance) + index_y * vertex_count_x + index_z * vertex_count_x * vertex_count_y;
                  vertex_ids[vertex_id_counter + 1] = vertex_index;

                  vertex_index = (index_x + point_distance) + (index_y + point_distance) * vertex_count_x + index_z * vertex_count_x * vertex_count_y;
                  vertex_ids[vertex_id_counter + 2] = vertex_index;

                  vertex_index = index_x + (index_y + point_distance) * vertex_count_x + index_z * vertex_count_x * vertex_count_y;
                  vertex_ids[vertex_id_counter + 3] = vertex_index;

                  vertex_index = index_x + index_y * vertex_count_x + (index_z + point_distance) * vertex_count_x * vertex_count_y;
                  vertex_ids[vertex_id_counter + 4] = vertex_index;

                  vertex_index = (index_x + point_distance) + index_y * vertex_count_x + (index_z + point_distance) * vertex_count_x * vertex_count_y;
                  vertex_ids[vertex_id_counter + 5] = vertex_index;

                  vertex_index = (index_x + point_distance) + (index_y + point_distance) * vertex_count_x + (index_z + point_distance) * vertex_count_x * vertex_count_y;
                  vertex_ids[vertex_id_counter + 6] = vertex_index;

                  vertex_index = index_x + (index_y + point_distance) * vertex_count_x + (index_z + point_distance) * vertex_count_x * vertex_count_y;
                  vertex_ids[vertex_id_counter + 7] = vertex_index;

                  vertex_id_counter += 8;
               }
            }
         }
      }//for node id
   } else {
      //if Vertex Filter Finest Level not activated generate all vertex coordinates and vertex ids locally for every node

      vertex_ids.resize( number_of_leaves_locally * NumberOfCellsPerBlock() * 8 );
      unsigned long long int vertex_id_counter = 0;
      constexpr unsigned int j_ids_skew = ( CC::ICX() + 1 );
      constexpr unsigned int k_ids_skew = ( CC::ICX() + 1 ) * ( CC::ICY() + 1 );

      unsigned int leaves_counter = 0;
      unsigned int const offset = topology_.LeafOffsetOfRank( MpiUtilities::MyRankId() );

      //Store offset in cell_vertex_ids for each node; needed for vertex filter
      std::vector<unsigned long long int> node_cell_vertex_offset;
      node_cell_vertex_offset.reserve(number_of_leaves_locally);

      for( [[maybe_unused]] auto const& node_id : topology_.LocalLeafIds() ) {
         //store cell_vetrex offset for each node
         node_cell_vertex_offset.push_back(vertex_id_counter);
         for( unsigned int k = 0; k < CC::ICZ(); ++k ) {
            for( unsigned int j = 0; j < CC::ICY(); ++j ) {
               for( unsigned int i = 0; i < CC::ICX(); ++i ) {
                  unsigned long long const shift = leaves_counter * NumberOfVerticesPerBlock() + offset * NumberOfVerticesPerBlock();
                  vertex_ids[vertex_id_counter    ] =  i    +  j    * j_ids_skew +  k      * k_ids_skew + shift;
                  vertex_ids[vertex_id_counter + 1] = (i+1) +  j    * j_ids_skew +  k      * k_ids_skew + shift;
                  vertex_ids[vertex_id_counter + 2] = (i+1) + (j+1) * j_ids_skew +  k      * k_ids_skew + shift;
                  vertex_ids[vertex_id_counter + 3] =  i    + (j+1) * j_ids_skew +  k      * k_ids_skew + shift;
                  vertex_ids[vertex_id_counter + 4] =  i    +  j    * j_ids_skew + (k + 1) * k_ids_skew + shift;
                  vertex_ids[vertex_id_counter + 5] = (i+1) +  j    * j_ids_skew + (k + 1) * k_ids_skew + shift;
                  vertex_ids[vertex_id_counter + 6] = (i+1) + (j+1) * j_ids_skew + (k + 1) * k_ids_skew + shift;
                  vertex_ids[vertex_id_counter + 7] =  i    + (j+1) * j_ids_skew + (k + 1) * k_ids_skew + shift;
                  vertex_id_counter += 8;
               }
            }
        }
        leaves_counter++;
      }

      // Filter operation for the Mpi vertex filter type
      if constexpr( CC::VERTEX_FILTER() == VertexFilterType::Mpi ) {
         //Approach fixe node values afterwards
         //separate synchronization for x, y and z axis

         //create list of send vectors to each other rank; this is introduced to deal with the creation of tags
         std::vector<std::vector<std::pair<std::uint64_t, std::uint64_t>>> send_list_per_rank( number_of_ranks );
         std::vector<std::vector<std::pair<std::uint64_t, std::uint64_t>>> recv_list_per_rank( number_of_ranks );
         std::vector<MPI_Request> mpi_send_request_list; //vector to track all send requests
         std::vector<std::uint64_t> const local_leaf_ids = topology_.LocalLeafIds(); //get local leaves
         unsigned int const local_leaf_count = local_leaf_ids.size();

         //x-Direction (send West to East)
         VertexFilter::InitializeSendRecvRankLists(send_list_per_rank, recv_list_per_rank, local_leaf_ids, VertexFilter::Direction::X, topology_);
         unsigned int const x_boundary_size = (CC::ICZ()+1)*(CC::ICY()+1);
         //init send buffer for sends; exact size unknown -> size is upper bound
         std::vector<unsigned long long int> send_buffers(local_leaf_count * x_boundary_size);

         //send data
         for(unsigned int leaf_index = 0; leaf_index< local_leaf_count; leaf_index++) {
            std::uint64_t node_id = local_leaf_ids[leaf_index];
            unsigned int send_buffer_offset = leaf_index * x_boundary_size;
            unsigned long long int const vertex_offset = node_cell_vertex_offset[leaf_index];
            unsigned int const vertex_id_jump_y = 8 * CC::ICX();
            unsigned int const vertex_id_jump_z = 8 * CC::ICX() * CC::ICY();
            //load inner points
            VertexFilter::LoadBoundaryBuffer(vertex_ids.data() + vertex_offset, send_buffers.data() + send_buffer_offset, vertex_id_jump_y, vertex_id_jump_z, CC::ICY(), CC::ICZ(), 0, 3, 4, 7);
            //send data to all ranks necessary
            VertexFilter::SendFromNode(node_id, topology_.GetNeighboringLeaves(node_id, BoundaryLocation::West), send_buffers.data() + send_buffer_offset, VertexFilter::Direction::X, send_list_per_rank, mpi_send_request_list, topology_);
         }

         //recv data
         for(unsigned int leaf_index = 0; leaf_index< local_leaf_ids.size(); leaf_index++) {
            std::uint64_t const node_id = local_leaf_ids[leaf_index];
            std::vector<unsigned long long int> east_boundary_points(x_boundary_size);
            unsigned long long int const vertex_offset = node_cell_vertex_offset[leaf_index] + 8 * (CC::ICX()-1);
            unsigned int const vertex_id_jump_y = 8 * CC::ICX();
            unsigned int const vertex_id_jump_z = 8 * CC::ICX() * CC::ICY();
            //init east vetrices
            VertexFilter::LoadBoundaryBuffer( vertex_ids.data() + vertex_offset, east_boundary_points.data(), vertex_id_jump_y, vertex_id_jump_z, CC::ICY(), CC::ICZ(), 1, 2, 5, 6 );
            //receive data for node
            VertexFilter::RecvToNode(node_id, topology_.GetNeighboringLeaves(node_id, BoundaryLocation::East), east_boundary_points.data(), VertexFilter::Direction::X, recv_list_per_rank, topology_ );

            //insert east boundary into cell vertices; iterate over every boundary block
            for(unsigned int k = 0; k < CC::ICZ(); k++) {
               for(unsigned int j = 0; j < CC::ICY(); j++) {
                  unsigned long long int const vertex_offset = node_cell_vertex_offset[leaf_index];
                  unsigned long long int const block_index = 8 * (CC::ICX()-1) + 8 * j * CC::ICX() + 8 * k * CC::ICX() * CC::ICY() + vertex_offset;
                  vertex_ids[block_index + 1] = east_boundary_points[j     + k     * (CC::ICY()+1)];
                  vertex_ids[block_index + 2] = east_boundary_points[(j+1) + k     * (CC::ICY()+1)];
                  vertex_ids[block_index + 5] = east_boundary_points[j     + (k+1) * (CC::ICY()+1)];
                  vertex_ids[block_index + 6] = east_boundary_points[(j+1) + (k+1) * (CC::ICY()+1)];
               }
            }
         }

         //wait for all messages from x direction
         MPI_Waitall(mpi_send_request_list.size(), mpi_send_request_list.data(), MPI_STATUSES_IGNORE);

         //y direction (send South to North)
         unsigned int const y_boundary_size = (CC::ICZ()+1)*(CC::ICX()+1);

         if constexpr(CC::DIM() != Dimension::One) { //in one-dimensional case no communication along the y-axis is needed

            VertexFilter::InitializeSendRecvRankLists(send_list_per_rank, recv_list_per_rank, local_leaf_ids, VertexFilter::Direction::Y, topology_);
            mpi_send_request_list.clear();
            if(x_boundary_size != y_boundary_size) //only init new sendbuffer if size differs from existing one
               send_buffers.resize(local_leaf_count * y_boundary_size);
            for(unsigned int leaf_index = 0; leaf_index< local_leaf_count; leaf_index++) {
               std::uint64_t const node_id = local_leaf_ids[leaf_index];
               unsigned int const send_buffer_offset = leaf_index * y_boundary_size;
               unsigned long long int const vertex_offset = node_cell_vertex_offset[leaf_index];
               unsigned int const vertex_id_jump_x = 8;
               unsigned int const vertex_id_jump_z = 8 * CC::ICX() * CC::ICY();
               //load into buffer
               VertexFilter::LoadBoundaryBuffer(vertex_ids.data() + vertex_offset, send_buffers.data() + send_buffer_offset, vertex_id_jump_x, vertex_id_jump_z, CC::ICX(), CC::ICZ(), 0, 1, 4, 5);

               //send data to all ranks necessary
               VertexFilter::SendFromNode(node_id, topology_.GetNeighboringLeaves(node_id, BoundaryLocation::South), send_buffers.data() + send_buffer_offset, VertexFilter::Direction::Y, send_list_per_rank, mpi_send_request_list, topology_);
            }

            //recv all data
            for(unsigned int leaf_index = 0; leaf_index< local_leaf_ids.size(); leaf_index++) {
               std::uint64_t const node_id = local_leaf_ids[leaf_index];
               //init north vetrices
               std::vector<unsigned long long int> north_boundary_points(y_boundary_size);
               unsigned long long int const vertex_offset = node_cell_vertex_offset[leaf_index] + 8 * CC::ICX() * (CC::ICY()-1);
               unsigned int const vertex_id_jump_x = 8;
               unsigned int const vertex_id_jump_z = 8 * CC::ICX() * CC::ICY();
               VertexFilter::LoadBoundaryBuffer( vertex_ids.data() + vertex_offset, north_boundary_points.data(), vertex_id_jump_x, vertex_id_jump_z,CC::ICX(), CC::ICZ(), 3, 2, 7, 6 );
               VertexFilter::RecvToNode(node_id, topology_.GetNeighboringLeaves(node_id, BoundaryLocation::North), north_boundary_points.data(), VertexFilter::Direction::Y, recv_list_per_rank, topology_);
               //insert north boundary into cell vertices; iterate over every boundary block
               for(unsigned int k = 0; k < CC::ICZ(); k++) {
                  for(unsigned int i = 0; i < CC::ICX(); i++) {
                     unsigned long long int const vertex_offset = node_cell_vertex_offset[leaf_index];
                     unsigned long long int const block_index = 8 * CC::ICX() * CC::ICY() * k + i * 8 + 8 * CC::ICX() * (CC::ICY()-1) + vertex_offset;
                     vertex_ids[block_index + 3] = north_boundary_points[i     + k     * (CC::ICX() + 1)];
                     vertex_ids[block_index + 2] = north_boundary_points[(i+1) + k     * (CC::ICX() + 1)];
                     vertex_ids[block_index + 7] = north_boundary_points[i     + (k+1) * (CC::ICX() + 1)];
                     vertex_ids[block_index + 6] = north_boundary_points[(i+1) + (k+1) * (CC::ICX() + 1)];
                  }
               }
            }

            MPI_Waitall(mpi_send_request_list.size(), mpi_send_request_list.data(), MPI_STATUSES_IGNORE);
         } // if Dim != 1

          //z-direction (send Bottom to Top)
         if constexpr(CC::DIM() == Dimension::Three) {

            VertexFilter::InitializeSendRecvRankLists(send_list_per_rank, recv_list_per_rank, local_leaf_ids, VertexFilter::Direction::Z, topology_);
            mpi_send_request_list.clear();
            unsigned int const z_boundary_size = (CC::ICX()+1)*(CC::ICY()+1);
            if(y_boundary_size != z_boundary_size)
               send_buffers.resize(local_leaf_count * z_boundary_size);
            for(unsigned int leaf_index = 0; leaf_index< local_leaf_count; leaf_index++) {
               std::uint64_t const node_id = local_leaf_ids[leaf_index];
               //Do MPI communication of boundaries
               unsigned int const send_buffer_offset = leaf_index * z_boundary_size;
               //load inner points
               unsigned long long int const vertex_offset = node_cell_vertex_offset[leaf_index];
               unsigned int const vertex_id_jump_x = 8;
               unsigned int const vertex_id_jump_y = 8 * CC::ICX();
               VertexFilter::LoadBoundaryBuffer( vertex_ids.data() + vertex_offset, send_buffers.data() + send_buffer_offset, vertex_id_jump_x, vertex_id_jump_y, CC::ICX(), CC::ICY(), 0, 1, 3, 2);
               //send data to all ranks necessary
               VertexFilter::SendFromNode(node_id, topology_.GetNeighboringLeaves(node_id, BoundaryLocation::Bottom), send_buffers.data() + send_buffer_offset, VertexFilter::Direction::Z, send_list_per_rank, mpi_send_request_list, topology_);
            }

            //recv all data
            for(unsigned int leaf_index = 0; leaf_index< local_leaf_ids.size(); leaf_index++) {
               std::uint64_t const node_id = local_leaf_ids[leaf_index];
               //init boundary buffer
               std::vector<unsigned long long int> top_boundary_points(z_boundary_size);

               unsigned long long int const vertex_offset = node_cell_vertex_offset[leaf_index] + 8 * CC::ICX() * CC::ICY() * (CC::ICZ()-1);
               unsigned int const vertex_id_jump_x = 8;
               unsigned int const vertex_id_jump_y = 8 * CC::ICX();
               VertexFilter::LoadBoundaryBuffer(vertex_ids.data() + vertex_offset, top_boundary_points.data(), vertex_id_jump_x, vertex_id_jump_y, CC::ICX(), CC::ICY(), 4, 5, 7, 6);
               VertexFilter::RecvToNode(node_id, topology_.GetNeighboringLeaves(node_id, BoundaryLocation::Top), top_boundary_points.data(), VertexFilter::Direction::Z, recv_list_per_rank, topology_);
               //insert boundary buffer into cell vertices; iterate over every boundary block
               for(unsigned int j = 0; j < CC::ICY(); j++) {
                  for(unsigned int i = 0; i < CC::ICX(); i++) {
                     unsigned long long int const vertex_offset = node_cell_vertex_offset[leaf_index];
                     unsigned long long int const block_index = 8 * (i + j * CC::ICX() + CC::ICX() * CC::ICY() * (CC::ICZ()-1)) + vertex_offset;
                     vertex_ids[block_index + 4] = top_boundary_points[i     +  j    * (CC::ICX()+1)];
                     vertex_ids[block_index + 5] = top_boundary_points[(i+1) +  j    * (CC::ICX()+1)];
                     vertex_ids[block_index + 7] = top_boundary_points[i     + (j+1) * (CC::ICX()+1)];
                     vertex_ids[block_index + 6] = top_boundary_points[(i+1) + (j+1) * (CC::ICX()+1)];
                  }
               }
            }
            //wait for communication in z direction to be completed
            MPI_Waitall(mpi_send_request_list.size(), mpi_send_request_list.data(), MPI_STATUSES_IGNORE);
         } // if Dim == 3

      } // Vertex Filter Mpi
   } // Not Filter Finest Level
   return vertex_ids;
}
