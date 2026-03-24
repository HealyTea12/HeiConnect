#pragma once
#include <vector>
#include <concepts>

template <typename node_T>
    requires std::integral<node_T>
struct CactusMinCuts
{
public:
    CactusMinCuts(std::vector<node_T> tree_cuts, std::vector<std::pair<node_T, node_T>> cycle_cuts)
        : TREE_CUTS(std::move(tree_cuts)), CYCLE_CUTS(std::move(cycle_cuts)) {}
    size_t get_n_min_cuts() const
    {
        return TREE_CUTS.size() + CYCLE_CUTS.size();
    }

    /// @brief the child node of the tree edge that defines the cut
    const std::vector<node_T> TREE_CUTS;
    /// @brief the pair of child nodes of the cycle edges that define the cut, stored in traversal order
    const std::vector<std::pair<node_T, node_T>> CYCLE_CUTS;
};