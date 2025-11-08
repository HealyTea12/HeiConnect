#pragma once
#include <vector>
#include <ranges>
#include <unordered_map>
#include <unordered_set>
#include <concepts>

template <typename R, typename V>
concept RangeOf = std::ranges::range<R> && std::same_as<std::ranges::range_value_t<R>, V>;

template <typename T>
concept MutableGraph = requires(T g, T::NodeID u, T::NodeID v, T::EdgeID e) {
    { g.add_node(u) } -> std::same_as<void>;
    { g.add_edge(u, v) } -> std::same_as<void>;
    { g.num_nodes() } -> std::same_as<typename T::NodeID>;
    { g.num_edges() } -> std::same_as<typename T::EdgeID>;
    { g.contract_edge(u, v) } -> std::same_as<void>;
    { g.get_neighbors(u) } -> RangeOf<typename T::NodeID>;
};

typedef struct CRFGraph
{
    std::vector<size_t> vertices;
    std::vector<size_t> edges;

    CRFGraph(const std::vector<size_t> &vertices,
             const std::vector<size_t> &edges)
        : vertices(vertices), edges(edges) {}

} CRFGraph;

typedef struct WeightedCRFGraph
{
    CRFGraph graph;
    std::vector<double> weights;

} WeightedCRFGraph;

// A simple mutable directed graph implementation using an adjacency list
template <typename NodeID, typename EdgeID>
class MutableGraphAdjacency
{
public:
    // Construct graph with n nodes [0, n-1]
    MutableGraphAdjacency(NodeID n) : m_num_nodes(n), m_num_edges(0)
    {
        m_adjacency_list = std::unordered_map<NodeID, std::vector<NodeID>>{};
        m_adjacency_list.reserve(n);
        for (NodeID u = 0; u < n; ++u)
        {
            m_adjacency_list[u] = std::vector<NodeID>{};
        }
    }
    void add_edge(NodeID u, NodeID v)
    {
        m_adjacency_list[u].emplace_back(v);
        m_num_edges++;
    }
    void add_node(NodeID u)
    {
        m_adjacency_list[u] = std::vector<NodeID>{};
        m_num_nodes++;
    }
    NodeID num_nodes()
    {
        return this->m_num_nodes;
    }
    EdgeID num_edges()
    {
        return this->m_num_edges;
    }

    void contract_edge(NodeID u, NodeID v)
    {
        for (auto &neighbor : m_adjacency_list[v])
        {
            m_adjacency_list[u].emplace_back(neighbor);
        }
        drop_node(v);
    }

    void drop_node(NodeID u)
    {
        m_adjacency_list.erase(u);
        m_num_nodes--;
    }

private:
    std::unordered_map<NodeID, std::vector<NodeID>> m_adjacency_list;
    NodeID m_num_nodes;
    EdgeID m_num_edges;
};

// static_assert(MutableGraph<MutableGraphAdjacency<size_t>>);

WeightedCRFGraph
generate_links(const CRFGraph &graph, double (*edge_weight)(size_t u, size_t v));

WeightedCRFGraph generate_links(const WeightedCRFGraph &graph, double (*edge_weight)(size_t u, size_t v));

// Adds edges (for instance from the set cover solution) to the mutable graph
// static_assert(MutableGraph<mutable_graph>);
WeightedCRFGraph add_links(const WeightedCRFGraph &graph,
                           const WeightedCRFGraph &link_graph,
                           const std::unordered_set<size_t> &selected_edges);

WeightedCRFGraph make_bidirectional(const WeightedCRFGraph &graph);