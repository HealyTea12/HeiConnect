#pragma once

#include <algorithm>
#include <concepts>
#include <cstddef>
#include <ranges>
#include <type_traits>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace HeiConnect::conn_aug
{

    template <typename T, typename NodeId>
    concept NeighborValue = std::convertible_to<T, NodeId>;

    template <typename GraphT, typename NodeId, typename EdgeId, typename Weight>
    concept ImmutableGraph = std::integral<NodeId> && std::integral<EdgeId> && std::floating_point<Weight> && requires(const GraphT &g, const NodeId u, const NodeId v, const EdgeId e) {
        { g.num_vertices() } -> std::convertible_to<NodeId>;
        { g.num_edges() } -> std::convertible_to<EdgeId>;
        { g.is_edge(u, v) } -> std::convertible_to<bool>;
        { g.neighbors(u) } -> std::ranges::input_range;
        requires NeighborValue<std::ranges::range_value_t<decltype(g.neighbors(u))>, NodeId>;
        { g.edge_weight(e) } -> std::convertible_to<Weight>;
    };

    template <std::integral NodeId = std::size_t, std::floating_point Weight = double>
    struct Link
    {
        NodeId u{};
        NodeId v{};
        Weight cost{};

        bool is_self_loop() const
        {
            return u == v;
        }

        void normalize_undirected()
        {
            if (v < u)
            {
                std::swap(u, v);
            }
        }
    };

    template <typename LinkRangeT, typename LinkT>
    concept LinkRange = std::ranges::input_range<LinkRangeT> &&
                        std::convertible_to<std::ranges::range_value_t<LinkRangeT>, LinkT>;

    template <typename GraphT,
              std::integral NodeId = std::size_t,
              std::integral EdgeId = std::size_t,
              std::floating_point Weight = double,
              typename LinkRangeT = std::vector<Link<NodeId, Weight>>>
        requires ImmutableGraph<GraphT, NodeId, EdgeId, Weight> && LinkRange<LinkRangeT, Link<NodeId, Weight>>
    class Instance
    {
    public:
        using GraphType = GraphT;
        using NodeIdType = NodeId;
        using EdgeIdType = EdgeId;
        using WeightType = Weight;
        using LinkType = Link<NodeId, Weight>;
        using LinkRangeType = LinkRangeT;

        Instance(GraphT base_graph,
                 LinkRangeT candidate_links,
                 const std::string instance_name = "")
            : m_graph(std::move(base_graph)),
              m_links(std::move(candidate_links)),
              m_name(std::move(instance_name))
        {
        }

        const GraphT &base_graph() const
        {
            return m_graph;
        }

        const LinkRangeT &links() const
        {
            return m_links;
        }

        const std::string &name() const
        {
            return m_name;
        }

        size_t num_candidate_links() const
        {
            if constexpr (std::ranges::sized_range<LinkRangeT>)
            {
                return static_cast<std::size_t>(std::ranges::size(m_links));
            }
            return static_cast<std::size_t>(std::ranges::distance(m_links));
        }

        void allow_links_on_existing_edges(bool allow = true)
        {
            m_allow_links_on_existing_edges = allow;
        }

        bool allows_links_on_existing_edges() const
        {
            return m_allow_links_on_existing_edges;
        }

        void validate() const
        {
            NodeId n = m_graph.num_vertices();
            std::size_t i = 0;
            for (const auto &link : m_links)
            {
                if (link.is_self_loop())
                {
                    throw std::invalid_argument("Self-loop candidate link at index " + std::to_string(i));
                }
                if (link.u >= n || link.v >= n)
                {
                    throw std::out_of_range("Candidate link endpoint out of graph range at index " + std::to_string(i));
                }
                if (!m_allow_links_on_existing_edges && m_graph.is_edge(static_cast<size_t>(link.u), static_cast<std::size_t>(link.v)))
                {
                    throw std::invalid_argument("Candidate link duplicates base graph edge at index " + std::to_string(i));
                }
                ++i;
            }
        }

        void normalize_links()
            requires std::ranges::range<LinkRangeT> && std::is_lvalue_reference_v<std::ranges::range_reference_t<LinkRangeT>>
        {
            for (auto &link : m_links)
            {
                link.normalize_undirected();
            }
        }

    private:
        GraphT m_graph;
        LinkRangeT m_links;
        std::string m_name;

        bool m_allow_links_on_existing_edges = false;
    };

    enum class SolveStatus
    {
        Optimal,
        Feasible,
        Infeasible,
        TimeLimit,
        Error,
        Unknown // maybe not necessary
    };

    template <std::floating_point Weight = double>
    struct AugmentationSolution
    {
        SolveStatus status = SolveStatus::Unknown;
        std::vector<size_t> selected_link_ids;
        Weight objective_value{};
        std::string message;
    };

    template <typename SolverT, typename InstanceT, typename ResultT>
    concept ConnectivityAugmentationSolver = requires(SolverT solver, const InstanceT &instance) {
        { solver.solve(instance) } -> std::same_as<ResultT>;
    };
}
