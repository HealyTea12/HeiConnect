#pragma once

#include "HeiConnect/set_cover/common.hpp"

#include <ranges>

// TODO: set covers are assuming size_t indices everywhere, should improve with templates
template<class link_node_T, class link_edge_T, class link_weight_T>
    requires std::integral<link_node_T> && std::integral<link_edge_T>
class SetCoverOracle
{
public:
    using TIn = uint32_t;
    using TOut = uint32_t;
    using ElementID = size_t;
    using SetID = size_t;
    using SetCost = link_weight_T;
    using ull = unsigned long long;
    struct Link
    {
        link_node_T u;
        link_node_T v;
        SetCost weight;
    };

public:
    SetCoverOracle(
        std::vector<TIn> tin,
        std::vector<TOut> tout,
        CactusMinCuts<link_node_T> min_cuts,
        std::vector<Link> m_links) :
        m_tIn(std::move(tin)),
        m_tOut(std::move(tout)),
        m_minCuts(std::move(min_cuts)),
        m_links(std::move(m_links)) {

        };

private:
    std::vector<TIn> m_tIn;
    std::vector<TOut> m_tOut;
    CactusMinCuts<link_node_T> m_minCuts;
    std::vector<Link> m_links;

public:
    size_t get_num_sets() const noexcept
    {
        return m_links.size();
    }
    size_t get_num_elements() const noexcept
    {
        return m_minCuts.get_n_min_cuts();
    }
    SetCost get_set_cost(SetID set_index) const
    {
        return m_links[set_index].weight;
    }
    /* NOT SUPPORTED
    auto set_elements(size_t set_index) const
    {
        (void)set_index;
        return std::ranges::empty_view<size_t>{};
    }
    */
    size_t get_num_tree_cuts() const
    {
        return m_minCuts.TREE_CUTS.size();
    }
    size_t get_num_cycle_cuts() const
    {
        return m_minCuts.CYCLE_CUTS.size();
    }
    bool covers_tree(SetID set_index, ElementID element_index) const
    {
        auto [u, v, w] = m_links[set_index];
        auto cut = m_minCuts.TREE_CUTS[element_index];
        bool is_ancestor_cut = isAncestor(cut, u) != isAncestor(cut, v);
        return is_ancestor_cut;
    }
    bool covers_cycle(SetID set_index, ElementID element_index) const
    {
        auto [u, v, w] = m_links[set_index];
        auto cut = m_minCuts.CYCLE_CUTS[element_index];
        bool u_belongs = isAncestor(cut.first, u) && !isAncestor(cut.second, u);
        bool v_belongs = isAncestor(cut.first, v) && !isAncestor(cut.second, v);
        return u_belongs != v_belongs;
    }
    void forEachElement(SetID set_index, const std::function<void(ElementID)>& func) const
    {
        auto [u, v, w] = m_links[set_index];
        (void)w;
        for (size_t i = 0; i < m_minCuts.TREE_CUTS.size(); i++)
        {
            if (covers_tree(set_index, i))
                func(i);
        }
        for (size_t i = 0; i < m_minCuts.CYCLE_CUTS.size(); i++)
        {
            if (covers_cycle(set_index, i))
                func(m_minCuts.TREE_CUTS.size() + i);
        }
    }
    void forEachSet(ElementID element_index, const std::function<void(SetID)>& func) const
    {
        if (element_index < m_minCuts.TREE_CUTS.size())
        {
            for (size_t set_index = 0; set_index < m_links.size(); set_index++)
            {
                if (covers_tree(set_index, element_index))
                    func(set_index);
            }
        }
        else
        {
            size_t cycle_element_index = element_index - m_minCuts.TREE_CUTS.size();
            for (size_t set_index = 0; set_index < m_links.size(); set_index++)
            {
                if (covers_cycle(set_index, cycle_element_index))
                    func(set_index);
            }
        }
    }

    template<typename NewSetCost>
    SetCoverOracle<link_node_T, link_edge_T, NewSetCost> discretize_costs(size_t num_bins) const
    {
        using NewSetCoverOracle = SetCoverOracle<link_node_T, link_edge_T, NewSetCost>;
        auto weights_range = m_links | std::views::transform([](const Link& link) { return link.weight; });
        auto dc = discretize_weights<NewSetCost>(weights_range, num_bins);
        std::vector<typename NewSetCoverOracle::Link> new_links;
        new_links.reserve(m_links.size());
        for (size_t i = 0; i < m_links.size(); ++i)
        {
            const auto& link = m_links[i];
            new_links.emplace_back(Link{link.u, link.v, static_cast<NewSetCost>(dc[i])});
        }
        return NewSetCoverOracle(m_tIn, m_tOut, m_minCuts, new_links);
    }

private:
    bool isAncestor(size_t anc, size_t node) const
    {
        return m_tIn[anc] <= m_tIn[node] && m_tOut[anc] >= m_tOut[node];
    }
};
