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
#ifndef TOPOLOGY_NODE_H
#define TOPOLOGY_NODE_H

#include <vector>
#include <tuple>
#include <cstdint>
#include "materials/material_names.h"

#include "space_filling_curves.h"

/**
 * @brief The TopologyNode class organizes the light weight global (over MPI ranks) node information in a tree structure. Allowing the TopologyManager efficient searches.
 */
class TopologyNode {

   std::uint64_t const unique_id_;
   int current_rank_;
   int future_rank_;
   std::vector<TopologyNode> children_;
   bool is_leaf_;
   std::vector<MaterialName> materials_;

   void Refine();

   void AddFluid(MaterialName const material);
   void RemoveFluid(MaterialName const material);
   std::vector<MaterialName> GetFluids() const;
   MaterialName GetSingleFluid() const;

   TopologyNode& GetChildWithId(std::uint64_t const id);
   const TopologyNode& GetChildWithId(std::uint64_t const id) const;

   unsigned int Weight() const;

public:
   TopologyNode() = delete;
   explicit TopologyNode(std::uint64_t const id, int const rank = -1);
   ~TopologyNode() = default;
   TopologyNode( TopologyNode const& ) = delete;
   TopologyNode& operator=( TopologyNode const& ) = delete;
   TopologyNode( TopologyNode&& ) = default; //Needed for usage in vector.
   TopologyNode& operator=( TopologyNode&& ) = delete;

   void Refine(std::uint64_t const id);
   void Coarse(std::uint64_t const id);

   int GetRank(std::uint64_t const id) const;
   void NodeLeafCount(std::pair<unsigned int, unsigned int>& nodes_and_leaves) const;
   void RankWiseNodeLeafCount(std::vector<std::pair<unsigned int, unsigned int>>& nodes_leaves_per_rank) const;
   void NodeBlockCount(std::pair<unsigned int, unsigned int>& nodes_and_leaves) const;
   unsigned int MultiPhaseNodeCount() const;
   void RankWiseNodeBlockCount(std::vector<std::pair<unsigned int, unsigned int>>& nodes_blocks_per_rank) const;
   unsigned int GetDepth() const;

   void IdsOnLevel(unsigned int const level,std::vector<std::uint64_t>& ids) const;
   void LocalIdsOnLevel(unsigned int const level,std::vector<std::uint64_t>& ids, int const local_rank) const;

   unsigned int LocalLeaves(std::vector<std::uint64_t>& local_leaves, int const rank) const;
   unsigned int GetLeafIds(std::vector<std::uint64_t>& leaves) const;
   unsigned int LocalLeavesOnLevel(std::vector<std::uint64_t>& local_leaves, int const rank, unsigned int const level, unsigned int const current_level = 0) const;
   unsigned int GetLeafIdsOnLevel(std::vector<std::uint64_t>& leaves, unsigned int const level, unsigned int const current_level = 0) const;

   void SetTargetRankForLeaf(std::vector<std::vector<std::tuple<unsigned int,int>>>& count_rank_map, unsigned int const level = 0);
   void SetTargetRankForLeaf(std::vector<std::vector<std::tuple<unsigned int,int>>>& count_rank_map, const HilbertPosition position,unsigned int const level = 0);

   void ListUnbalancedNodes(std::vector<std::tuple<std::uint64_t const, int const, int const> >& ids_current_future_rank_map);

   int BalanceTargetRanks();
   void SetCurrentRankOfLeaf(std::uint64_t const id, int const rank);

   bool NodeExists(std::uint64_t const id) const;
   bool NodeIsLeaf(std::uint64_t const id) const;

   void ChildWeight(std::vector<unsigned int>& weights_on_level, unsigned int const current_level = 0) const;

   void AddFluid(std::uint64_t const id, MaterialName const material);
   void RemoveFluid(std::uint64_t const id, MaterialName const material);
   std::vector<MaterialName> GetFluids(std::uint64_t const id) const;
   MaterialName GetSingleFluid(std::uint64_t const id) const;

   /**
    * @brief Gives the id of this topology node.
    * @return id.
    */
   inline std::uint64_t Id() const {return unique_id_;}

   /**
    * @brief Gives the current rank of the TopologyNode.
    * @return Rank.
    */
   inline int GetRank() const {return current_rank_;}

   /**
    * @brief Specification of == operator for the id of topology nodes
    * @param rhs right hand side value of == operator
    * @return True if Ids of two nodes are equal, False otherwise 
    */
   inline bool operator==(std::uint64_t const rhs) {return (rhs == unique_id_);}
   inline bool operator==(const TopologyNode& rhs) {return (rhs.Id() == unique_id_);}
};

#endif // TOPOLOGY_NODE_H
