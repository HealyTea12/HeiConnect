#pragma once

#include <algorithm>
#include <map>
#include <utility>
#include <vector>

#include "HeiConnect/conn_aug/connectivity_augmentation.hpp"

namespace HeiConnect::conn_aug
{
    class WheelCon
    {
    public:
        template<typename InstanceT>
        AugmentationSolution<typename InstanceT::WeightType> solve(const InstanceT& instance) const
        {
            using NodeId = typename InstanceT::NodeIdType;
            using Weight = typename InstanceT::WeightType;

            const auto& graph = instance.base_graph();
            AugmentationSolution<Weight> solution;
            if (graph.num_vertices() == 0 || graph.num_edges() == 0)
            {
                solution.status = SolveStatus::Feasible;
                return solution;
            }

            auto cactus_order = std::vector<NodeId>{};
            graph.cactus_for_each_hamiltonian_vertex(
                NodeId{0},
                [&](NodeId u) { cactus_order.push_back(u); });

            auto appearances = std::vector<size_t>(graph.num_vertices(), 0);
            for (const NodeId u : cactus_order)
            {
                ++appearances[u];
            }

            auto leaves = std::vector<NodeId>{};
            for (const NodeId u : cactus_order)
            {
                // A vertex appears once exactly when its cactus weighted degree is two.
                if (appearances[u] == 1)
                {
                    leaves.push_back(u);
                }
            }

            if (leaves.size() % 2 != 0)
            {
                solution.status = SolveStatus::Error;
                solution.message = "WheelCon requires an even number of cactus leaves";
                return solution;
            }

            struct LinkChoice
            {
                size_t id;
                Weight cost;
            };
            auto link_choices = std::map<std::pair<NodeId, NodeId>, LinkChoice>{};
            size_t link_id = 0;
            for (const auto& link : instance.links())
            {
                auto endpoints = std::minmax(link.u, link.v);
                const auto key = std::pair<NodeId, NodeId>{endpoints.first, endpoints.second};
                const auto found = link_choices.find(key);
                if (found == link_choices.end() || link.cost < found->second.cost)
                {
                    link_choices.insert_or_assign(key, LinkChoice{link_id, link.cost});
                }
                ++link_id;
            }

            const size_t offset = leaves.size() / 2;
            solution.selected_link_ids.reserve(offset);
            for (size_t i = 0; i < offset; ++i)
            {
                auto endpoints = std::minmax(leaves[i], leaves[i + offset]);
                const auto found = link_choices.find({endpoints.first, endpoints.second});
                if (found == link_choices.end())
                {
                    solution.status = SolveStatus::Error;
                    solution.selected_link_ids.clear();
                    solution.objective_value = Weight{};
                    solution.message = "WheelCon matching link is missing from the instance";
                    return solution;
                }
                solution.selected_link_ids.push_back(found->second.id);
                solution.objective_value += found->second.cost;
            }

            solution.status = SolveStatus::Feasible;
            return solution;
        }
    };
}
