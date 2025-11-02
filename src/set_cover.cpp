#include "set_cover.hpp"
#include "graph.hpp"
#include <omp.h>

namespace solver
{
    template <typename node_T, typename edge_T>
    std::vector<graph::CactusCut> get_all_cuts_cactus(
        const std::vector<edge_T> &Vertices,
        const std::vector<node_T> &Edges,
        const std::vector<double> &weight)
        requires std::is_integral<node_T>::value && std::is_integral<edge_T>::value
    {
        std::vector<graph::CactusCut> cuts = {};
        // edges that are cycles
        std::vector<std::pair<node_T, node_T>> cycle_edges{};
        // Parallel BFS for detecting cycles in cactus graph
        node_T root = 0;
        std::vector<node_T> parent{Vertices.size(), -1};
        std::vector<int> distance{Vertices.size(), -1};
        std::vector<double> weight_to_parent{Vertices.size(), -1.};
        parent[root] = root;
        distance[root] = 0;
        std::vector<std::vector<node_T>> levels{};

        std::vector<node_T> current_frontier{Vertices.size(), -1}; // current frontier
        std::vector<node_T> next_frontier{Vertices.size(), -1};    // next frontier
        int size_current_frontier = 0;
        int size_next_frontier = 0;
        current_frontier[size_current_frontier++] = root;

#pragma omp parallel firstprivate(S, T)
        {
            std::vector<node_T> local_current_frontier{Vertices.size(), -1};
            auto local_count = 0ull;
            while (size_current_frontier)
            {
                local_count = 0;
#pragma omp for nowait
                for (int i{0}; i < size_current_frontier; i++)
                {
                    int u = current_frontier[i];
                    for (int j{Vertices[u]}; j < Vertices[u + 1]; j++)
                    {
                        v = Edges[j];
                        if (d[v] < 0)
                        { // not visited
                            distance[v] = distance[u] + 1;
                            parent[v] = u;
                            weight_to_parent[v] = weight[j];
                            local_current_frontier[local_count++] = v;
                        }
                        else
                        { // visited => implies cycle
                            parent[v] = u;
                            distance[v] = distance[u] + 1;
                            weight_to_parent[v] = weight[j];
#pragma omp critical
                            {
                                cycle_edges.emplace_back({u, v}); // probably need critical
                            }
                        }
                    }
                }
            }

            // calculate the offset: from what point on in the global frontier the local frontier should be copied to
            int offset;
#pragma omp critical
            {
                offset = size_next_frontier;
                size_next_frontier += local_count;
            }

            // move the local frontier into the global frontier
            for (int i{0}; i < local_count; i++)
            {
                next_frontier[i + offset] = local_current_frontier[i];
            }
            // Swap next and current frontier (there is a difference between single and master)
#pragma omp barrier
#pragma omp single
            {
                levels.push_back(current_frontier);
                std::swap(next_frontier, current_frontier);
                size_current_frontier = size_next_frontier;
                size_next_frontier = 0;
            }
        }

        // Calculate cycles since it is a cactus, each cycle will have one entry point in the BFS.
        // If there was another entry point, then it would break the cactus property.
        // So to find the cycles, we go backwards from both ends of the cycle edges until we get to a common one.

        // The point of calculating the cycles is that the min cuts are either two edges of a cycle, or a single
        // tree edge. So here we store the edges we have already looked at.
        double min_cut = min_cut_cactus(Vertices, Edges, weight);
        std::unordered_set<std::pair<node_T, node_T>> all_cycle_edges;
#pragma omp parallel for
        for (size_t i{0}; i < cycle_edges.size(); i++)
        {
            std::vector<std::pair<node_T, node_T>> cycle{};
            node_T u = cycle_edges[i].first;
            node_T v = cycle_edges[i].second;
            // make sure u is always the 'top' one
            if (distance[u] > distance[v])
            {
                std::swap(u, v);
                all_cycle_edges.emplace({v, parent[v]});
                cycle.emplace_back({v, parent[v]});
                v = parent[v];
            }
            all_cycle_edges.emplace({u, v});
            cycle.emplace_back({u, v});
            while (u != v)
            {
                all_cycle_edges.emplace({u, parent[u]});
                cycle.emplace_back({u, parent[u]});
                all_cycle_edges.emplace({v, parent[v]});
                cycle.emplace_back({v, parent[v]});
                u = parent[u];
                v = parent[v];
            }
            // find the cuts in this cycle: all combinations of two edges.
            for (size_t j{0}; j < cycle.size() - 1; j++) // ce = cycle edge
            {
                node_T u = cycle[j].first;
                node_T v = cycle[j].second;
                for (size_t k{j + 1}; k < cycle.size(); k++)
                {
                    node_T x = cycle[k].first;
                    node_T y = cycle[k].second;
                    if (
                        weight_between_nodes(Edges, weight, cycle[j].first, cycle[j].second) +
                            weight_between_nodes(Edges, weight, cycle[k].first, cycle[k].second) ==
                        min_cut)
                    {
#pragma omp critical
                        {
                            cuts.emplace_back(
                                Cut{{u, v},
                                    {x, y}});
                        }
                    }
                }
            }
        }
// now add all tree edges that are min cuts
#pragma omp parallel for
        for (size_t i{0}; i < Vertices.size() - 1; i++)
        {
            for (size_t j{Vertices[i]}; j < Vertices[i + 1]; j++)
            {
                node_T u = i;
                node_T v = Edges[j];
                if (all_cycle_edges.find({u, v}) == all_cycle_edges.end() && // not a cycle edge
                    all_cycle_edges.find({v, u}) == all_cycle_edges.end() &&
                    weight[j] == min_cut)
                {
#pragma omp critical
                    {
                        cuts.emplace_back(
                            Cut{{u, v}});
                    }
                }
            }
        }
        return cuts;
    }

    struct SetCover
    {
        std::vector<size_t> a;
        std::vector<size_t> b;
        std::vector<double> costs;
    };

    template <typename node_T, typename edge_T>
    SetCover generate_set_cover_memory(
        const std::vector<graph::CactusCut> &cuts,
        const std::vector<edge_T> &Vertices,
        const std::vector<node_T> &Edges,
        const std::vector<node_T> &children,
        const std::vector<graph::Edge> &links)
    {
        std::vector<size_t> a{cuts.size() + 1, 0};
        std::vector<size_t> b{};
        std::vector<std::unordered_set<node_T>> cut_nodes{cuts.size()};
#pragma omp parallel for
        for (size_t i{0}; i < cuts.size(); i++)
        {
            auto cut = cuts[i];
            if (!cut.e2.has_value()) // tree cut
            {
                auto edge = cut.e1;
                if (parent[edge.first] == edge.second)
                {
                    if (children[edge.first].size() < Vertices.size() / 2)
                    {
                        cut_nodes[i] = bfs_ignore_edges(Vertices, Edges, edge.first, edge);
                    }
                    else
                    {
                        cut_nodes[i] = bfs_ignore_edges(Vertices, Edges, edge.second, edge);
                    }
                }
                else
                {
                    if (children[edge.second].size() < Vertices.size() / 2)
                    {
                        cut_nodes[i] = bfs_ignore_edges(Vertices, Edges, edge.second, edge);
                    }
                    else
                    {
                        cut_nodes[i] = bfs_ignore_edges(Vertices, Edges, edge.first, edge);
                    }
                }
            }
            else // cycle cut
            {
                auto edge1 = cut.e1;
                auto edge2 = cut.e2.value();
                cut_nodes[i] = bfs_ignore_edges(
                    Vertices, Edges, edge1.first, edge1, edge2);
            }
        }
        for (size_t i{0}; i < links.size(); i++)
        {
            auto u = links[i].first;
            auto v = links[i].second;
            for (size_t j{0}; j < cuts.size(); j++)
            {
                if (cut_nodes[j].find(u) == cut_nodes[j].end() && cut_nodes[j].find(v) != cut_nodes[j].end())
                {
                    // link covers cut
                    a[i]++;
                    b.push_back(j);
                }
            }
        }
        std::vector<double> costs{links.size()};
        for (size_t i{0}; i < links.size(); i++)
        {
            costs[i] = links[i].weight;
        }
        return SetCover(a, b, costs);
    }

    // under construction
    template <typename node_T, typename edge_T>
    SetCover generate_set_cover(
        const std::vector<graph::CactusCut> &cuts,
        const std::vector<edge_T> &Vertices,
        const std::vector<node_T> &Edges,
        const std::vector<node_T> &children,
        const std::vector<graph::Edge> &links)
    {
        std::vector<size_t> a{cuts.size() + 1, 0};
        std::vector<size_t> b{};
#pragma omp parallel for
        for (size_t i{0}; i < cuts.size(); i++)
        {
            auto cut = cuts[i];
            std::unordered_set<node_T> cut_nodes{};
            if (!cut.e2.has_value()) // tree cut
            {
                auto edge = cut.e1;
                if (parent[edge.first] == edge.second)
                {
                    if (children[edge.first].size() < Vertices.size() / 2)
                    {
                        cut_nodes = bfs_ignore_edges(Vertices, Edges, edge.first, edge);
                    }
                    else
                    {
                        cut_nodes = bfs_ignore_edges(Vertices, Edges, edge.second, edge);
                    }
                }
                else
                {
                    if (children[edge.second].size() < Vertices.size() / 2)
                    {
                        cut_nodes = bfs_ignore_edges(Vertices, Edges, edge.second, edge);
                    }
                    else
                    {
                        cut_nodes = bfs_ignore_edges(Vertices, Edges, edge.first, edge);
                    }
                }
            }
            else // cycle cut
            // TODO: take into consideration which side is smaller
            {
                auto edge1 = cut.e1;
                auto edge2 = cut.e2.value();
                cut_nodes = bfs_ignore_edges(
                    Vertices, Edges, edge1.first, edge1, edge2)
            }
        }
        return SetCover(a, b, costs);
    }

    std::list<graph::Edge>
    set_cover(const graph::Cactus &g, const std::vector<graph::Edge> &links)
    {
        std::vector<graph::Edge> solution;
        std::list<graph::CactusCut> cactus_cuts = g.get_min_cuts();

        std::vector<unsigned long> edges = {links.size()};
        std::vector<unsigned long> cuts = {};
        for (auto &l : links)
        {
            for (auto i = 0ull; i < cactus_cuts.size(); i++)
            {
                if (cactus_cuts[i].is_edge_cut(l.first, l.second))
                {
                    cuts.push_back(c);
                }
            }
        }
    }
}
}
}