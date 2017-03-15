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
#include "topology/topology_manager.h"
#include "input_output/output_mesh_generator.h"

namespace {

   /**
    * @brief Helper function to generate load balanced topology with single node on level zero and one refinement step.
    * @param topology Topology to be created (indirect return).
    * @param number_of_ranks Number of ranks that are used.
    */
   void RefineFirstNodeInTopologyAndLoadBalance( TopologyManager& topology, int const number_of_ranks ) {
      REQUIRE( topology.NodeAndLeafCount() == std::pair<unsigned int, unsigned int>( 1, 1 ) );
      topology.RefineNodeWithId( 0x1400000 );
      topology.UpdateTopology();
      // For load balancing to work, the nodes must have a weight = fluids inside them
      topology.AddFluidToNode( 0x1400000, MaterialName::StiffenedGas );
      topology.AddFluidToNode( 0xA000000, MaterialName::StiffenedGas );
      topology.AddFluidToNode( 0xA000001, MaterialName::StiffenedGas );
      topology.AddFluidToNode( 0xA000002, MaterialName::StiffenedGas );
      topology.AddFluidToNode( 0xA000003, MaterialName::StiffenedGas );
      topology.AddFluidToNode( 0xA000004, MaterialName::StiffenedGas );
      topology.AddFluidToNode( 0xA000005, MaterialName::StiffenedGas );
      topology.AddFluidToNode( 0xA000006, MaterialName::StiffenedGas );
      topology.AddFluidToNode( 0xA000007, MaterialName::StiffenedGas );
      topology.UpdateTopology();
      topology.GetLoadBalancedTopology( number_of_ranks );
   }

   /**
    * @brief Helper function mirroring partly GetLocalCountAndStart() in output_mesh_generator.h. Get local number of vertices.
    * @param total_vertex_count Total number of vertices on all ranks.
    * @param rank Rank identifier.
    * @param number_of_ranks Number of ranks considered.
    * @return Local number of vertices (rank-wise).
    */
   constexpr unsigned long long LocalVertexCount( unsigned long long const total_vertex_count, int const rank, int const number_of_ranks ) {
      unsigned long long local_count = total_vertex_count / number_of_ranks;
      int const remainder = total_vertex_count % number_of_ranks;
      if( rank < remainder ) { local_count++; }
      return local_count;
   }
}

SCENARIO( "Dimension of global vertex coordinates vector is correctly computed", "[1rank]" ) {

   GIVEN( "An output-mesh generator with an underlying topology with Lmax being one" ) {
      TopologyManager topology = TopologyManager( 1 );
      OutputMeshGenerator const mesh_generator = OutputMeshGenerator( topology );

      WHEN( "The topology consists of just one node and Vertex Filter is off" ) {
         // Default for the topology.
         constexpr VertexFilterType vertex_filter = VertexFilterType::Off;

         THEN( "The global vertex coordinates dimension vector has two entries with the first being equal the number of vertices and the second being three" ) {
            auto const global_vertex_coordinates_dimensions = mesh_generator.GlobalVertexCoordinatesDimensions<vertex_filter>();
            REQUIRE( global_vertex_coordinates_dimensions.size() == 2 );
            REQUIRE( global_vertex_coordinates_dimensions[0] == NumberOfVerticesPerBlock() );
            REQUIRE( global_vertex_coordinates_dimensions[1] == 3 );
         }
      }
      WHEN( "The topology consists of just one node and Vertex Filter is Mpi" ) {
         // Default for the topology.
         constexpr VertexFilterType vertex_filter = VertexFilterType::Mpi;

         THEN( "The global vertex coordinates dimension vector has two entries with the first being equal the number of vertices and the second being three" ) {
            auto const global_vertex_coordinates_dimensions = mesh_generator.GlobalVertexCoordinatesDimensions<vertex_filter>();
            REQUIRE( global_vertex_coordinates_dimensions.size() == 2 );
            REQUIRE( global_vertex_coordinates_dimensions[0] == NumberOfVerticesPerBlock() );
            REQUIRE( global_vertex_coordinates_dimensions[1] == 3 );
         }
      }
      WHEN( "The topology consists of just one node and Vertex Filter is Finest Level" ) {
         // Default for the topology.
         constexpr VertexFilterType vertex_filter = VertexFilterType::FinestLevel;

         THEN( "The global vertex coordinates dimension vector has two entries with the first being equal the number of vertices and the second being three" ) {
            auto const global_vertex_coordinates_dimensions = mesh_generator.GlobalVertexCoordinatesDimensions<vertex_filter>();
            REQUIRE( global_vertex_coordinates_dimensions.size() == 2 );
            REQUIRE( global_vertex_coordinates_dimensions[0] == NumberOfVerticesPerBlock() );
            REQUIRE( global_vertex_coordinates_dimensions[1] == 3 );
         }
      }
      WHEN( "The node in the topology is refined and Vertex Filter is Off" ) {
         topology.RefineNodeWithId( 0x1400000 );
         topology.UpdateTopology();
         REQUIRE( topology.NodeAndLeafCount() == std::pair<unsigned int, unsigned int>( 9, 8 ) );
         constexpr VertexFilterType vertex_filter = VertexFilterType::Off;

         THEN( "The global vertex coordinates dimension vector has two entries with the first being equal the number of vertices time the number of leafs (8) and the second being three" ) {
            auto const global_vertex_coordinates_dimensions = mesh_generator.GlobalVertexCoordinatesDimensions<vertex_filter>();
            REQUIRE( global_vertex_coordinates_dimensions.size() == 2 );
            REQUIRE( global_vertex_coordinates_dimensions[0] == NumberOfVerticesPerBlock() * 8 );
            REQUIRE( global_vertex_coordinates_dimensions[1] == 3 );
         }
      }
      WHEN( "The node in the topology is refined and Vertex Filter is Mpi" ) {
         topology.RefineNodeWithId( 0x1400000 );
         topology.UpdateTopology();
         REQUIRE( topology.NodeAndLeafCount() == std::pair<unsigned int, unsigned int>( 9, 8 ) );
         constexpr VertexFilterType vertex_filter = VertexFilterType::Mpi;

         THEN( "The global vertex coordinates dimension vector has two entries with the first being equal the number of vertices time the number of leafs (8) and the second being three" ) {
            auto const global_vertex_coordinates_dimensions = mesh_generator.GlobalVertexCoordinatesDimensions<vertex_filter>();
            REQUIRE( global_vertex_coordinates_dimensions.size() == 2 );
            REQUIRE( global_vertex_coordinates_dimensions[0] == NumberOfVerticesPerBlock() * 8 );
            REQUIRE( global_vertex_coordinates_dimensions[1] == 3 );
         }
      }
      WHEN( "The node in the topology is refined and Vertex Filter is FinestLevel" ) {
         topology.RefineNodeWithId( 0x1400000 );
         topology.UpdateTopology();
         REQUIRE( topology.NodeAndLeafCount() == std::pair<unsigned int, unsigned int>( 9, 8 ) );
         constexpr VertexFilterType vertex_filter = VertexFilterType::FinestLevel;

         THEN( "The global vertex coordinates dimension vector has two entries with the first being is equal (ICX*2+1)*(ICY*2+1)*(ICZ*2+1) and the second being three" ) {
            auto const global_vertex_coordinates_dimensions = mesh_generator.GlobalVertexCoordinatesDimensions<vertex_filter>();
            REQUIRE( global_vertex_coordinates_dimensions.size() == 2 );
            REQUIRE( global_vertex_coordinates_dimensions[0] == ( CC::ICX() * 2 + 1 ) * ( CC::ICY() * 2 + 1 ) * ( CC::ICZ() * 2 + 1 ) );
            REQUIRE( global_vertex_coordinates_dimensions[1] == 3 );
         }
      }
   }
}

SCENARIO( "Vertex Coordinates write count and start are correctly determined", "[1rank]" ) {

   GIVEN( "An output-mesh generator with an underlying topology with Lmax being one" ) {
      TopologyManager topology = TopologyManager( 1 );
      OutputMeshGenerator const mesh_generator = OutputMeshGenerator( topology );

      WHEN( "The topology consists of just one node and Vertex Filter is Mpi" ) {
         // Default for the topology.
         constexpr VertexFilterType vertex_filter = VertexFilterType::Mpi;

         THEN( "The write count vector has two entries with the first equal the number of vertieces in a block and the second being three" ) {
            auto const start_and_write = mesh_generator.WriteStartAndCountCoordinates<vertex_filter>( std::get<1>( topology.NodeAndLeafCount() ) );
            REQUIRE( std::get<1>( start_and_write ).size() == 2 );
            REQUIRE( std::get<1>( start_and_write )[0] == NumberOfVerticesPerBlock() );
            REQUIRE( std::get<1>( start_and_write )[1] == 3 );
         }
         THEN( "The start vector has two entries both beeing zero on rank 0" ) {
            auto const start_and_write = mesh_generator.WriteStartAndCountCoordinates<vertex_filter>( std::get<1>( topology.NodeAndLeafCount() ) );
            REQUIRE( std::get<0>( start_and_write ).size() == 2 );
            REQUIRE( std::get<0>( start_and_write )[0] == 0 );
            REQUIRE( std::get<0>( start_and_write )[1] == 0 );
         }
      }
      WHEN( "The topology consists of just one node and Vertex Filter is Off" ) {
         // Default for the topology.
         constexpr VertexFilterType vertex_filter = VertexFilterType::Off;

         THEN( "The write count vector has two entries with the first equal the number of vertieces in a block and the second being three" ) {
            auto const start_and_write = mesh_generator.WriteStartAndCountCoordinates<vertex_filter>( std::get<1>( topology.NodeAndLeafCount() ) );
            REQUIRE( std::get<1>( start_and_write ).size() == 2 );
            REQUIRE( std::get<1>( start_and_write )[0] == NumberOfVerticesPerBlock() );
            REQUIRE( std::get<1>( start_and_write )[1] == 3 );
         }
         THEN( "The start vector has two entries both beeing zero on rank 0" ) {
            auto const start_and_write = mesh_generator.WriteStartAndCountCoordinates<vertex_filter>( std::get<1>( topology.NodeAndLeafCount() ) );
            REQUIRE( std::get<0>( start_and_write ).size() == 2 );
            REQUIRE( std::get<0>( start_and_write )[0] == 0 );
            REQUIRE( std::get<0>( start_and_write )[1] == 0 );
         }
      }
      WHEN( "The topology consists of just one node and Vertex Filter is FinestLevel" ) {
         // Default for the topology.
         constexpr VertexFilterType vertex_filter = VertexFilterType::FinestLevel;

         THEN( "The write count vector has two entries with the first equal the total vertex count on the finest level and the second being three" ) {
            auto const start_and_write = mesh_generator.WriteStartAndCountCoordinates<vertex_filter>( std::get<1>( topology.NodeAndLeafCount() ) );
            REQUIRE( std::get<1>( start_and_write ).size() == 2 );
            REQUIRE( std::get<1>( start_and_write )[0] == TotalVertexCountOnFinestLevel( topology.GetCurrentMaximumLevel(), topology.GetLevelZeroBlockRatioXyz() ) );
            REQUIRE( std::get<1>( start_and_write )[1] == 3 );
         }
         THEN( "The start vector has two entries both beeing zero on rank 0" ) {
            auto const start_and_write = mesh_generator.WriteStartAndCountCoordinates<vertex_filter>( std::get<1>( topology.NodeAndLeafCount() ) );
            REQUIRE( std::get<0>( start_and_write ).size() == 2 );
            REQUIRE( std::get<0>( start_and_write )[0] == 0 );
            REQUIRE( std::get<0>( start_and_write )[1] == 0 );
         }
      }
      WHEN( "The node in the topology is refined and load balanced onto 3 ranks and Vertex Filter is Off" ) {
         RefineFirstNodeInTopologyAndLoadBalance( topology, 3 );
         auto const nodes_and_leaves_per_rank = topology.NodesAndLeavesPerRank();

         REQUIRE( topology.NodeAndLeafCount() == std::pair<unsigned int, unsigned int>( 9, 8 ) );
         constexpr VertexFilterType vertex_filter = VertexFilterType::Off;
         THEN( "The write count vector has two entries with the first equal the number of vertices in a block times the number of leaves on rank 0 and the second being three" ) {
            auto const start_and_write = mesh_generator.WriteStartAndCountCoordinates<vertex_filter>( std::get<1>( nodes_and_leaves_per_rank[0] ) );
            REQUIRE( std::get<1>( start_and_write ).size() == 2 );
            REQUIRE( std::get<1>( start_and_write )[0] == NumberOfVerticesPerBlock() * std::get<1>( nodes_and_leaves_per_rank[0] ) );
            REQUIRE( std::get<1>( start_and_write )[1] == 3 );
         }
         THEN( "The write count vector has two entries with the first equal the number of vertices in a block times the number of leaves on rank 1 and the second being three" ) {
            auto const start_and_write = mesh_generator.WriteStartAndCountCoordinates<vertex_filter>( std::get<1>( nodes_and_leaves_per_rank[1] ), 1 );
            REQUIRE( std::get<1>( start_and_write ).size() == 2 );
            REQUIRE( std::get<1>( start_and_write )[0] == NumberOfVerticesPerBlock() * std::get<1>( nodes_and_leaves_per_rank[1] ) );
            REQUIRE( std::get<1>( start_and_write )[1] == 3 );
         }
         THEN( "The write count vector has two entries with the first equal the number of vertices in a block times the number of leaves on rank 2second being three" ) {
            auto const start_and_write = mesh_generator.WriteStartAndCountCoordinates<vertex_filter>( std::get<1>( nodes_and_leaves_per_rank[2] ), 2 );
            REQUIRE( std::get<1>( start_and_write ).size() == 2 );
            REQUIRE( std::get<1>( start_and_write )[0] == NumberOfVerticesPerBlock() * std::get<1>( nodes_and_leaves_per_rank[2] ) );
            REQUIRE( std::get<1>( start_and_write )[1] == 3 );
         }

         THEN( "The start vector has two entries both beeing zero on rank 0" ) {
            auto const start_and_write = mesh_generator.WriteStartAndCountCoordinates<vertex_filter>( std::get<1>( topology.NodeAndLeafCount() ) );
            REQUIRE( std::get<0>( start_and_write ).size() == 2 );
            REQUIRE( std::get<0>( start_and_write )[0] == 0 );
            REQUIRE( std::get<0>( start_and_write )[1] == 0 );
         }
         THEN( "The start vector has two entries #leaves on rank 0 * Number of Vertices per Block and 0 on rank 1" ) {
            auto const start_and_write = mesh_generator.WriteStartAndCountCoordinates<vertex_filter>( std::get<1>( topology.NodeAndLeafCount() ), 1 );
            REQUIRE( std::get<0>( start_and_write ).size() == 2 );
            REQUIRE( std::get<0>( start_and_write )[0] == NumberOfVerticesPerBlock() * std::get<1>( nodes_and_leaves_per_rank[0] ) );
            REQUIRE( std::get<0>( start_and_write )[1] == 0 );
         }
         THEN( "The start vector has two entries ( #leaves on rank 0 + #leaves on rank 1 ) * Number of Vertices per Block - 0 on rank 2" ) {
            auto const start_and_write = mesh_generator.WriteStartAndCountCoordinates<vertex_filter>( std::get<1>( topology.NodeAndLeafCount() ), 2 );
            REQUIRE( std::get<0>( start_and_write ).size() == 2 );
            REQUIRE( std::get<0>( start_and_write )[0] == NumberOfVerticesPerBlock() * ( std::get<1>( nodes_and_leaves_per_rank[0] ) + std::get<1>( nodes_and_leaves_per_rank[1] ) ) );
            REQUIRE( std::get<0>( start_and_write )[1] == 0 );
         }
      }

      WHEN( "The node in the topology is refined and load balanced onto 3 ranks and Vertex Filter is Mpi" ) {
         RefineFirstNodeInTopologyAndLoadBalance( topology, 3 );
         auto const nodes_and_leaves_per_rank = topology.NodesAndLeavesPerRank();

         REQUIRE( topology.NodeAndLeafCount() == std::pair<unsigned int, unsigned int>( 9, 8 ) );
         constexpr VertexFilterType vertex_filter = VertexFilterType::Mpi;
         THEN( "The write count vector has two entries with the first equal the number of vertices in a block times the number of leaves on rank 0 and the second being three" ) {
            auto const start_and_write = mesh_generator.WriteStartAndCountCoordinates<vertex_filter>( std::get<1>( nodes_and_leaves_per_rank[0] ) );
            REQUIRE( std::get<1>( start_and_write ).size() == 2 );
            REQUIRE( std::get<1>( start_and_write )[0] == NumberOfVerticesPerBlock() * std::get<1>( nodes_and_leaves_per_rank[0] ) );
            REQUIRE( std::get<1>( start_and_write )[1] == 3 );
         }
         THEN( "The write count vector has two entries with the first equal the number of vertices in a block times the number of leaves on rank 1 and the second being three" ) {
            auto const start_and_write = mesh_generator.WriteStartAndCountCoordinates<vertex_filter>( std::get<1>( nodes_and_leaves_per_rank[1] ), 1 );
            REQUIRE( std::get<1>( start_and_write ).size() == 2 );
            REQUIRE( std::get<1>( start_and_write )[0] == NumberOfVerticesPerBlock() * std::get<1>( nodes_and_leaves_per_rank[1] ) );
            REQUIRE( std::get<1>( start_and_write )[1] == 3 );
         }
         THEN( "The write count vector has two entries with the first equal the number of vertices in a block times the number of leaves on rank 2second being three" ) {
            auto const start_and_write = mesh_generator.WriteStartAndCountCoordinates<vertex_filter>( std::get<1>( nodes_and_leaves_per_rank[2] ), 2 );
            REQUIRE( std::get<1>( start_and_write ).size() == 2 );
            REQUIRE( std::get<1>( start_and_write )[0] == NumberOfVerticesPerBlock() * std::get<1>( nodes_and_leaves_per_rank[2] ) );
            REQUIRE( std::get<1>( start_and_write )[1] == 3 );
         }

         THEN( "The start vector has two entries both beeing zero on rank 0" ) {
            auto const start_and_write = mesh_generator.WriteStartAndCountCoordinates<vertex_filter>( std::get<1>( topology.NodeAndLeafCount() ) );
            REQUIRE( std::get<0>( start_and_write ).size() == 2 );
            REQUIRE( std::get<0>( start_and_write )[0] == 0 );
            REQUIRE( std::get<0>( start_and_write )[1] == 0 );
         }
         THEN( "The start vector has two entries #leaves on rank 0 * Number of Vertices per Block and 0 on rank 1" ) {
            auto const start_and_write = mesh_generator.WriteStartAndCountCoordinates<vertex_filter>( std::get<1>( topology.NodeAndLeafCount() ), 1 );
            REQUIRE( std::get<0>( start_and_write ).size() == 2 );
            REQUIRE( std::get<0>( start_and_write )[0] == NumberOfVerticesPerBlock() * std::get<1>( nodes_and_leaves_per_rank[0] ) );
            REQUIRE( std::get<0>( start_and_write )[1] == 0 );
         }
         THEN( "The start vector has two entries ( #leaves on rank 0 + #leaves on rank 1 ) * Number of Vertices per Block - 0 on rank 2" ) {
            auto const start_and_write = mesh_generator.WriteStartAndCountCoordinates<vertex_filter>( std::get<1>( topology.NodeAndLeafCount() ), 2 );
            REQUIRE( std::get<0>( start_and_write ).size() == 2 );
            REQUIRE( std::get<0>( start_and_write )[0] == NumberOfVerticesPerBlock() * ( std::get<1>( nodes_and_leaves_per_rank[0] ) + std::get<1>( nodes_and_leaves_per_rank[1] ) ) );
            REQUIRE( std::get<0>( start_and_write )[1] == 0 );
         }
      }
      WHEN( "The node in the topology is refined and load balanced onto 3 ranks and Vertex Filter is FinestLevel" ) {
         constexpr int number_of_ranks = 3;
         RefineFirstNodeInTopologyAndLoadBalance( topology, number_of_ranks );
         auto const nodes_and_leaves_per_rank = topology.NodesAndLeavesPerRank();

         REQUIRE( topology.NodeAndLeafCount() == std::pair<unsigned int, unsigned int>( 9, 8 ) );
         constexpr VertexFilterType vertex_filter = VertexFilterType::FinestLevel;
         THEN( "The write count vector has two entries with the first equal the total vertex count divided by the number of ranks on rank 0 and the second being three" ) {
            constexpr int rank = 0;
            auto const start_and_write = mesh_generator.WriteStartAndCountCoordinates<vertex_filter>( std::get<1>( nodes_and_leaves_per_rank[0] ), rank, number_of_ranks );
            REQUIRE( std::get<1>( start_and_write ).size() == 2 );
            auto const local_count = LocalVertexCount( TotalVertexCountOnFinestLevel( topology.GetCurrentMaximumLevel(), topology.GetLevelZeroBlockRatioXyz() ), rank,  number_of_ranks );
            REQUIRE( std::get<1>( start_and_write )[0] ==  local_count );
            REQUIRE( std::get<1>( start_and_write )[1] == 3 );
         }
         THEN( "The write count vector has two entries with the first equal the total vertex count divided by the number of ranks on rank 1 and the second being three" ) {
            constexpr int rank = 1;
            auto const start_and_write = mesh_generator.WriteStartAndCountCoordinates<vertex_filter>( std::get<1>( nodes_and_leaves_per_rank[1] ), rank, number_of_ranks );
            REQUIRE( std::get<1>( start_and_write ).size() == 2 );
            auto const local_count = LocalVertexCount( TotalVertexCountOnFinestLevel( topology.GetCurrentMaximumLevel(), topology.GetLevelZeroBlockRatioXyz() ), rank,  number_of_ranks );
            REQUIRE( std::get<1>( start_and_write )[0] == local_count );
            REQUIRE( std::get<1>( start_and_write )[1] == 3 );
         }
         THEN( "The write count vector has two entries with the first equal the total vertex count divided by the number of ranks on rank 2 and the second being three" ) {
            constexpr int rank = 2;
            auto const start_and_write = mesh_generator.WriteStartAndCountCoordinates<vertex_filter>( std::get<1>( nodes_and_leaves_per_rank[2] ), rank, number_of_ranks );
            REQUIRE( std::get<1>( start_and_write ).size() == 2 );
            auto const local_count = LocalVertexCount( TotalVertexCountOnFinestLevel( topology.GetCurrentMaximumLevel(), topology.GetLevelZeroBlockRatioXyz() ), rank,  number_of_ranks );
            REQUIRE( std::get<1>( start_and_write )[0] == local_count );
            REQUIRE( std::get<1>( start_and_write )[1] == 3 );
         }

         THEN( "The start vector has two entries both beeing zero on rank 0" ) {
            auto const start_and_write = mesh_generator.WriteStartAndCountCoordinates<vertex_filter>( std::get<1>( topology.NodeAndLeafCount() ), 0, number_of_ranks );
            REQUIRE( std::get<0>( start_and_write ).size() == 2 );
            REQUIRE( std::get<0>( start_and_write )[0] == 0 );
            REQUIRE( std::get<0>( start_and_write )[1] == 0 );
         }
         THEN( "The start vector has two entries #leaves on rank 0 * Number of Vertices per Block and 0 on rank 1" ) {
            auto const start_and_write = mesh_generator.WriteStartAndCountCoordinates<vertex_filter>( std::get<1>( topology.NodeAndLeafCount() ), 1, number_of_ranks );
            REQUIRE( std::get<0>( start_and_write ).size() == 2 );
            auto const local_count_rank_zero = LocalVertexCount( TotalVertexCountOnFinestLevel( topology.GetCurrentMaximumLevel(), topology.GetLevelZeroBlockRatioXyz() ), 0, number_of_ranks );
            REQUIRE( std::get<0>( start_and_write )[0] == local_count_rank_zero );
            REQUIRE( std::get<0>( start_and_write )[1] == 0 );
         }
         THEN( "The start vector has two entries ( #leaves on rank 0 + #leaves on rank 1 ) * Number of Vertices per Block - 0 on rank 2" ) {
            auto const start_and_write = mesh_generator.WriteStartAndCountCoordinates<vertex_filter>( std::get<1>( topology.NodeAndLeafCount() ), 2, number_of_ranks );
            REQUIRE( std::get<0>( start_and_write ).size() == 2 );
            auto const local_count_rank_zero = LocalVertexCount( TotalVertexCountOnFinestLevel( topology.GetCurrentMaximumLevel(), topology.GetLevelZeroBlockRatioXyz() ), 0, number_of_ranks );
            auto const local_count_rank_one = LocalVertexCount( TotalVertexCountOnFinestLevel( topology.GetCurrentMaximumLevel(), topology.GetLevelZeroBlockRatioXyz() ), 1, number_of_ranks );
            REQUIRE( std::get<0>( start_and_write )[0] == local_count_rank_zero + local_count_rank_one );
            REQUIRE( std::get<0>( start_and_write )[1] == 0 );
         }
      }
   }
}

SCENARIO( "Vertex Coordinates are correctly written out", "[1rank]" ) {

   GIVEN( "An output-mesh generator with an underlying topology with Lmax being one" ) {
      TopologyManager topology = TopologyManager( 1 );
      OutputMeshGenerator const mesh_generator = OutputMeshGenerator( topology );

      WHEN( "The topology consist of just one node and the Vertex Filter is off" ) {
         auto const coordinates = mesh_generator.VertexCoordinates<VertexFilterType::Off>();

         THEN( "The coordinates vector has 3 * #of Leaves * #vertices per block entries and the first three all read zero and the last is one" ) {
            REQUIRE( coordinates.size() == NumberOfVerticesPerBlock() * 3 );
            REQUIRE( coordinates[0] == Approx( 0.0 ) );
            REQUIRE( coordinates[0] == Approx( 0.0 ) );
            REQUIRE( coordinates[0] == Approx( 0.0 ) );
            REQUIRE( coordinates.back() == Approx( 1.0 ) );
         }
      }

      WHEN( "The topology consist of just one node and the Vertex Filter is off" ) {
         auto const coordinates = mesh_generator.VertexCoordinates<VertexFilterType::Mpi>();

         THEN( "The coordinates vector has 3 * #of Leaves * #vertices per block entries and the first three all read zero and the last is one" ) {
            REQUIRE( coordinates.size() == NumberOfVerticesPerBlock() * 3 );
            REQUIRE( coordinates[0] == Approx( 0.0 ) );
            REQUIRE( coordinates[0] == Approx( 0.0 ) );
            REQUIRE( coordinates[0] == Approx( 0.0 ) );
            REQUIRE( coordinates.back() == Approx( 1.0 ) );
         }
      }

      WHEN( "The topology consist of just one node and the Vertex Filter is off" ) {
         auto const coordinates = mesh_generator.VertexCoordinates<VertexFilterType::FinestLevel>();

         THEN( "The coordinates vector has 3 * #of Leaves * #vertices per block entries and the first three all read zero and the last is one" ) {
            REQUIRE( coordinates.size() == NumberOfVerticesPerBlock() * 3 );
            REQUIRE( coordinates[0] == Approx( 0.0 ) );
            REQUIRE( coordinates[0] == Approx( 0.0 ) );
            REQUIRE( coordinates[0] == Approx( 0.0 ) );
            REQUIRE( coordinates.back() == Approx( 1.0 ) );
         }
      }
   }

   GIVEN( "An output-mesh generator with an underlying topology with Lmax being one and 1x1x2 blocks on L0" ) {
      TopologyManager topology = TopologyManager( 1, 1, 1, 2 );
      OutputMeshGenerator const mesh_generator = OutputMeshGenerator( topology );
      REQUIRE( topology.NodeAndLeafCount() == std::pair<unsigned int, unsigned int>( 2, 2 ) );

      WHEN( "The vertex filter is Off" ) {
         auto const coordinates = mesh_generator.VertexCoordinates<VertexFilterType::Off>();

         THEN( "The coordinates vector has 3 * #of Leaves * #vertices per block entries and the first three all read zero and the last one is 2" ) {
            REQUIRE( coordinates.size() == 2 * NumberOfVerticesPerBlock() * 3 );
            REQUIRE( coordinates[0] == Approx( 0.0 ) );
            REQUIRE( coordinates[0] == Approx( 0.0 ) );
            REQUIRE( coordinates[0] == Approx( 0.0 ) );
            REQUIRE( coordinates.back() == Approx( 2.0 ) );
         }
      }

      WHEN( "The vertex filter is Mpi" ) {
         auto const coordinates = mesh_generator.VertexCoordinates<VertexFilterType::Mpi>();

         THEN( "The coordinates vector has 3 * #of Leaves * #vertices per block entries and the first three all read zero and the last one is 2" ) {
            REQUIRE( coordinates.size() == 2 * NumberOfVerticesPerBlock() * 3 );
            REQUIRE( coordinates[0] == Approx( 0.0 ) );
            REQUIRE( coordinates[0] == Approx( 0.0 ) );
            REQUIRE( coordinates[0] == Approx( 0.0 ) );
            REQUIRE( coordinates.back() == Approx( 2.0 ) );
         }
      }

      WHEN( "The vertex filter is Finest Level" ) {
         auto const coordinates = mesh_generator.VertexCoordinates<VertexFilterType::FinestLevel>();

         THEN( "The coordinates vector has 3 * Sum ( #Block_i on L0 * IC + 1) entries and the first three all read zero and the last one is 2" ) {
            auto const level_zero_blocks_xyz = topology.GetLevelZeroBlockRatioXyz();
            unsigned int const required_size = 3 * ( ( CC::ICX() * level_zero_blocks_xyz[0] ) + 1 )
                                                 * ( ( CC::ICY() * level_zero_blocks_xyz[1] ) + 1 )
                                                 * ( ( CC::ICZ() * level_zero_blocks_xyz[2] ) + 1 );

            REQUIRE( coordinates.size() == required_size );
            REQUIRE( coordinates[0] == Approx( 0.0 ) );
            REQUIRE( coordinates[0] == Approx( 0.0 ) );
            REQUIRE( coordinates[0] == Approx( 0.0 ) );
            REQUIRE( coordinates.back() == Approx( 2.0 ) );
         }
      }
   }
}