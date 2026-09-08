#pragma once

#include <bit>
#include <concepts>
#include <cstdint>
#include <iostream>
#include <limits>
#include <queue>
#include <tuple>
#include <type_traits>
#include <vector>

#include <omp.h>
#include <immintrin.h>

#include "HeiConnect/bfs.hpp"
#include "HeiConnect/data_structures/immutable_graph.hpp"
#include "HeiConnect/sc_reduction/cactus_min_cuts.hpp"
#include "HeiConnect/sc_reduction/transform_single_core.hpp"
#include "HeiConnect/sc_reduction/transform_single_partition_matrix_utils.hpp"

using ull = unsigned long long;

namespace HeiConnect_details
{

    template <typename node_T, typename edge_T, typename weight_T>
        requires std::unsigned_integral<node_T> &&
                 std::unsigned_integral<edge_T>
    int calculate_number_min_cuts(
        const std::vector<edge_T> &vertices,
        const std::vector<node_T> &edges,
        const std::vector<weight_T> &weights,
        const std::vector<std::vector<std::pair<node_T, node_T>>> &cycles,
        const weight_T min_cut)
    {
        int n_min_cuts{0};
        for (node_T u{}; u < vertices.size() - 1; u++)
        {
            for (edge_T e{vertices[u]}; e < vertices[u + 1]; e++)
            {
                node_T v = edges[e];
                if (u < v)
                {
                    if (weights[e] == min_cut)
                    {
                        n_min_cuts++;
                    }
                }
            }
        }
        for (const auto &cyc : cycles)
        {
            for (size_t i{}; i < cyc.size() - 1; i++)
            {
                for (size_t j = i + 1; j < cyc.size(); j++)
                {
                    auto edge1 = cyc[i];
                    auto edge2 = cyc[j];
                    auto edge1_idx = get_edge_index(vertices, edges, edge1.first, edge1.second);
                    auto edge2_idx = get_edge_index(vertices, edges, edge2.first, edge2.second);
                    weight_T cycle_cut = weights[edge1_idx] + weights[edge2_idx];
                    if (cycle_cut == min_cut)
                    {
                        n_min_cuts++;
                    }
                }
            }
        }
        return n_min_cuts;
    }

    template <typename node_T, typename edge_T, typename weight_T>
        requires std::integral<node_T> && std::integral<edge_T>
    std::vector<ull> calculate_min_cut_partitions_ull( // maybe uint64_t
        const std::vector<edge_T> &vertices,
        const std::vector<node_T> &edges,
        const std::vector<weight_T> &weights,
        const std::vector<std::pair<node_T, node_T>> &cycle_edges_vec,
        const std::vector<bool> &cycle_edges,
        const std::vector<std::vector<std::pair<node_T, node_T>>> &cycles,
        weight_T min_cut,
        size_t n_min_cuts)
    {
        constexpr const size_t B = 8 * sizeof(ull);
        size_t n_vertices = vertices.size() - 1;
        const size_t N_COLS = n_min_cuts / B + (n_min_cuts % B != 0);
        std::vector<ull> min_cuts = std::vector<ull>(n_vertices * N_COLS, 0ULL);
        size_t current_min_cut_idx = 0;
        for (node_T u{}; u < vertices.size() - 1; u++)
        {
            for (edge_T i{vertices[u]}; i < vertices[u + 1]; i++)
            {
                node_T v = edges[i];
                if (u < v && !cycle_edges[i] && weights[i] <= min_cut)
                {
                    set_bit(min_cuts, u * N_COLS * B + current_min_cut_idx);
                    bfs_single_threaded(
                        vertices,
                        edges,
                        u,
                        [](node_T from, node_T to) noexcept {},
                        [&min_cuts, current_min_cut_idx, N_COLS](node_T from, node_T to) noexcept
                        {
                            set_bit(min_cuts, to * N_COLS * 8 * sizeof(ull) + current_min_cut_idx);
                        },
                        [u, v](node_T from, node_T to) noexcept
                        {
                            if (from == u && to == v)
                                return true;
                            return false;
                        });
                    current_min_cut_idx++;
                }
            }
        }
        for (const auto &cycle : cycles)
        {
            for (size_t i{}; i < cycle.size() - 1; i++)
            {
                for (size_t j = i + 1; j < cycle.size(); j++)
                {
                    auto edge1 = cycle[i];
                    auto edge2 = cycle[j];
                    auto edge1_idx = get_edge_index(vertices, edges, edge1.first, edge1.second);
                    auto edge2_idx = get_edge_index(vertices, edges, edge2.first, edge2.second);
                    double cycle_cut = weights[edge1_idx] + weights[edge2_idx];
                    if (cycle_cut == min_cut)
                    {
                        set_bit(min_cuts, edge1.first * N_COLS * 8 * sizeof(ull) + current_min_cut_idx);
                        bfs_single_threaded(
                            vertices,
                            edges,
                            edge1.first,
                            [](node_T from, node_T to) noexcept {},
                            [&min_cuts, current_min_cut_idx, N_COLS](node_T from, node_T to) noexcept
                            {
                                set_bit(min_cuts, to * N_COLS * 8 * sizeof(ull) + current_min_cut_idx);
                            },
                            [&edge1, &edge2](node_T from, node_T to) noexcept
                            {
                                if ((from == edge1.first && to == edge1.second) ||
                                    (from == edge2.first && to == edge2.second) ||
                                    (from == edge2.second && to == edge2.first))
                                    return true;
                                return false;
                            });
                        current_min_cut_idx++;
                    }
                }
            }
        }
        return min_cuts;
    }

    template <typename node_T, typename edge_T, typename weight_T>
        requires std::integral<node_T> && std::integral<edge_T>
    std::vector<unsigned char> calculate_min_cut_partitions_char(
        const std::vector<edge_T> &vertices,
        const std::vector<node_T> &edges,
        const std::vector<weight_T> &weights,
        const std::vector<std::pair<node_T, node_T>> &cycle_edges_vec,
        const std::vector<bool> &cycle_edges,
        const std::vector<std::vector<std::pair<node_T, node_T>>> &cycles,
        weight_T min_cut,
        size_t n_min_cuts)
    {
        std::vector<unsigned char> min_cuts{};
        size_t n_vertices = vertices.size() - 1;
        min_cuts.resize(n_vertices * n_min_cuts, 0);
        size_t current_min_cut_idx = 0;
        for (node_T u{}; u < vertices.size() - 1; u++)
        {
            for (edge_T i{vertices[u]}; i < vertices[u + 1]; i++)
            {
                node_T v = edges[i];
                if (u < v && !cycle_edges[i])
                {
                    if (weights[i] == min_cut)
                    {
                        min_cuts[u * n_min_cuts + current_min_cut_idx] = 1;
                        bfs_single_threaded(
                            vertices,
                            edges,
                            u,
                            [](node_T from, node_T to) noexcept {},
                            [&min_cuts, current_min_cut_idx, n_min_cuts](node_T from, node_T to) noexcept
                            {
                                min_cuts[to * n_min_cuts + current_min_cut_idx] = 1;
                            },
                            [u, v](node_T from, node_T to) noexcept
                            {
                                if (from == u && to == v)
                                    return true;
                                return false;
                            });
                        current_min_cut_idx++;
                    }
                }
            }
        }
        for (const auto &cycle : cycles)
        {
            for (size_t i{}; i < cycle.size() - 1; i++)
            {
                for (size_t j = i + 1; j < cycle.size(); j++)
                {
                    auto edge1 = cycle[i];
                    auto edge2 = cycle[j];
                    auto edge1_idx = get_edge_index(vertices, edges, edge1.first, edge1.second);
                    auto edge2_idx = get_edge_index(vertices, edges, edge2.first, edge2.second);
                    double cycle_cut = weights[edge1_idx] + weights[edge2_idx];
                    if (cycle_cut == min_cut)
                    {
                        min_cuts[edge1.first * n_min_cuts + current_min_cut_idx] = 1;
                        bfs_single_threaded(
                            vertices,
                            edges,
                            edge1.first,
                            [](node_T from, node_T to) noexcept {},
                            [&min_cuts, current_min_cut_idx, n_min_cuts](node_T from, node_T to) noexcept
                            {
                                min_cuts[to * n_min_cuts + current_min_cut_idx] = 1;
                            },
                            [&edge1, &edge2](node_T from, node_T to) noexcept
                            {
                                if ((from == edge1.first && to == edge1.second) ||
                                    (from == edge2.first && to == edge2.second) ||
                                    (from == edge2.second && to == edge2.first))
                                    return true;
                                return false;
                            });
                        current_min_cut_idx++;
                    }
                }
            }
        }
        return min_cuts;
    }

    // node_T/edge_T is a generic type that indexes nodes/edges e.g unsigned int
    template <typename node_T, typename weight_T>
    struct Edge
    {
        node_T u;
        node_T v;
        weight_T weight;
    };

    template <typename node_T, typename edge_T>
        requires std::integral<node_T> && std::integral<edge_T>
    std::tuple<std::vector<ull>, ull> generate_min_cut_matrix(
        const std::vector<edge_T> &vertices,
        const std::vector<node_T> &edges,
        const std::vector<double> &weights,
        const std::vector<size_t> &link_vertices,
        const std::vector<size_t> &link_edges,
        const std::vector<double> &link_weights)
    {
        // for checking during the second BFS
        auto cycle_edges = std::vector<bool>(edges.size(), false);
        std::vector<std::pair<node_T, node_T>> cycle_edge_vec{};
        std::vector<int> distances = std::vector<int>(vertices.size(), 0);
        std::vector<bool> visited_nodes = std::vector<bool>(vertices.size(), false);
        std::vector<node_T> parent = std::vector<node_T>(vertices.size(), 0);

        double start = omp_get_wtime();
        // BFS to find cycles and rooted tree
        std::queue<std::pair<node_T, int>> q;
        node_T root = static_cast<node_T>(0);
        parent[root] = root;
        q.push({root, 0});
        visited_nodes[root] = true;
        while (!q.empty())
        {
            auto curr = q.front();
            node_T current_node = curr.first;
            auto depth = curr.second;
            q.pop();
            for (size_t i{vertices[current_node]}; i < vertices[current_node + 1]; i++)
            {
                node_T neighbor = edges[i];
                if (!visited_nodes[neighbor])
                {
                    visited_nodes[neighbor] = true;
                    parent[neighbor] = current_node;
                    distances[neighbor] = depth + 1;
                    q.push({neighbor, depth + 1});
                }
                else
                {
                    if (distances[neighbor] < depth)
                        continue; // skip back edges to ancestors
                    if (distances[neighbor] == depth && neighbor < current_node)
                        continue; // count same-depth undirected edges only once
                    cycle_edge_vec.emplace_back(current_node, neighbor);
                }
            }
        }
        double end = omp_get_wtime();
        std::cout << "BFS to find cycles took " << (end - start) << " seconds." << std::endl;
        start = omp_get_wtime();
        find_cycles_and_root_tree(
            vertices,
            edges,
            weights,
            link_vertices,
            link_edges,
            link_weights);
        end = omp_get_wtime();
        std::cout << "Alternative finding cycles and rooting tree took " << (end - start) << " seconds." << std::endl;
        // for each cycle edge, reconstruct the cycle
        // could rethink data struct
        start = omp_get_wtime();
        auto cycles = std::vector<std::vector<std::pair<node_T, node_T>>>(cycle_edge_vec.size());
        std::cout << cycles.size() << " cycles found." << std::endl;
        assert(cycle_edge_vec.size() == cycles.size());
        for (size_t i{}; i < cycle_edge_vec.size(); i++)
        {
            auto &cycle = cycles[i];
            cycle.emplace_back(cycle_edge_vec[i]);
            node_T u = cycle_edge_vec[i].first;
            node_T v = cycle_edge_vec[i].second;
            cycle_edges[get_edge_index(vertices, edges, u, v)] = true;
            cycle_edges[get_edge_index(vertices, edges, v, u)] = true;
            if (distances[u] < distances[v])
            {
                cycle_edges[get_edge_index(vertices, edges, v, parent[v])] = true;
                cycle_edges[get_edge_index(vertices, edges, parent[v], v)] = true;
                cycle.emplace_back(v, parent[v]);
                v = parent[v];
            }
            while (u != v)
            {
                cycle.emplace_back(u, parent[u]);
                cycle.emplace_back(v, parent[v]);
                // could store the indices, so you don't have to look them up
                // could also set only one direction, since we know in what direction it will ge accessed
                cycle_edges[get_edge_index(vertices, edges, u, parent[u])] = true;
                cycle_edges[get_edge_index(vertices, edges, parent[u], u)] = true;
                cycle_edges[get_edge_index(vertices, edges, v, parent[v])] = true;
                cycle_edges[get_edge_index(vertices, edges, parent[v], v)] = true;
                u = parent[u];
                v = parent[v];
            }
        }
        end = omp_get_wtime();
        for (auto &cyc : cycles)
        {
            std::cout << cyc.size() << " cycle edges." << std::endl;
        }
        std::cout << "Cycle reconstruction took " << (end - start) << " seconds." << std::endl;
        // iterate over all edges to find min cuts
        // auto min_cut = global_mincut_simple({{vertices, edges},
        //                                     weights}); // should substitute with cactus min cut
        // this part is awful
        start = omp_get_wtime();
        auto min_cut = std::numeric_limits<double>::max();
        for (node_T u{}; u < vertices.size() - 1; u++)
        {
            for (edge_T e{vertices[u]}; e < vertices[u + 1]; e++)
            {
                if (!cycle_edges[e] && weights[e] < min_cut)
                {
                    min_cut = weights[e];
                }
            }
        }
        for (const auto &cycle : cycles)
        {
            for (auto i = 0; i < cycle.size() - 1; i++)
            {
                auto edge1 = cycle[i];
                auto edge1_idx = get_edge_index(vertices, edges, edge1.first, edge1.second);
                for (auto j = i + 1; j < cycle.size(); j++)
                {
                    auto edge2 = cycle[j];
                    auto edge_idx = get_edge_index(vertices, edges, edge2.first, edge2.second);
                    double cycle_cut = weights[edge1_idx] + weights[edge_idx];
                    if (cycle_cut < min_cut)
                    {
                        min_cut = cycle_cut;
                    }
                }
            }
        }
        end = omp_get_wtime();
        std::cout << "Min cut: " << min_cut << std::endl;
        std::cout << "Cactus min cut computation took " << (end - start) << " seconds." << std::endl;

        size_t n_min_cuts = 0;
        start = omp_get_wtime();
        n_min_cuts = calculate_number_min_cuts(
            vertices,
            edges,
            weights,
            cycles,
            min_cut);
        end = omp_get_wtime();
        std::cout << "Number of min cuts: " << n_min_cuts << std::endl;
        std::cout << "Counting number of min cuts took " << (end - start) << " seconds." << std::endl;
        // find and partition min cuts
        start = omp_get_wtime();
        auto min_cuts = calculate_min_cut_partitions_ull(
            vertices,
            edges,
            weights,
            cycle_edge_vec,
            cycle_edges,
            cycles,
            min_cut,
            n_min_cuts);
        end = omp_get_wtime();
        std::cout << "Min cut vertex partitioning took " << (end - start) << " seconds." << std::endl;
        return {min_cuts, n_min_cuts};
    }

    template <typename node_T, typename edge_T>
        requires std::integral<node_T> && std::integral<edge_T>
    std::vector<ull> __calculate_min_cut_partitions_ull_ancestry(
        const std::vector<edge_T> &vertices,
        const CactusMinCuts<node_T> &cactus_min_cuts,
        const std::vector<timer_type> &tin,
        const std::vector<timer_type> &tout)
    {
        size_t n_min_cuts = cactus_min_cuts.get_n_min_cuts();
        size_t n_vertices = vertices.size() - 1;
        const unsigned long B = 8 * sizeof(ull);
        const unsigned long N_COLS = (n_min_cuts + B - 1) / B;
        std::vector<ull> min_cuts = std::vector<ull>(n_vertices * N_COLS, 0ULL);
        for (node_T u{0}; u < n_vertices; u++)
        {

            for (size_t cut_idx{0}; cut_idx < cactus_min_cuts.TREE_CUTS.size(); cut_idx++)
            {
                // THIS CAN AND SHOULD BE VECTORISED
                const node_T &cut = cactus_min_cuts.TREE_CUTS[cut_idx];
                if (isAncestor(tin, tout, cut, u))
                {
                    size_t word_idx = cut_idx / B;
                    size_t bit_idx = cut_idx % B;
                    min_cuts[u * N_COLS + word_idx] |= (1ULL << bit_idx);
                }
            }
            // in construction
            for (size_t cut_idx{0}; cut_idx < cactus_min_cuts.CYCLE_CUTS.size(); cut_idx++)
            {
                const auto &[cut1, cut2] = cactus_min_cuts.CYCLE_CUTS[cut_idx];
                if (isAncestor(tin, tout, cut1, u) && !isAncestor(tin, tout, cut2, u))
                {
                    size_t overall_cut_idx = cactus_min_cuts.TREE_CUTS.size() + cut_idx;
                    size_t word_idx = overall_cut_idx / B;
                    size_t bit_idx = overall_cut_idx % B;
                    min_cuts[u * N_COLS + word_idx] |= (1ULL << bit_idx);
                }
            }
        }
        return min_cuts;
    }

    // Constructs the partition matrix ull by checking ancestry relationships in the DFS tree
    // Instead of doing a BFS for each min cut
    template <class node_T, class edge_T, class weight_T>
        requires std::integral<node_T> && std::integral<edge_T>
    std::vector<ull> calculate_min_cut_partitions_ull_ancestry(
        const std::vector<edge_T> &vertices,
        const std::vector<node_T> &edges,
        const std::vector<weight_T> &weights)
    {
        auto [tin, tout, cycles, is_cycle_edge, parent] = dfs_tin_tout_cycles(vertices, edges);
        weight_T min_cut = calculate_cactus_min_cut(vertices, edges, weights, cycles, is_cycle_edge);
        auto cactus_min_cuts = calculate_all_cactus_min_cuts(vertices, edges, weights, cycles, is_cycle_edge, min_cut, parent, static_cast<node_T>(0));
        return __calculate_min_cut_partitions_ull_ancestry(vertices, cactus_min_cuts, tin, tout);
    }

    template <typename node_T, typename edge_T>
        requires std::integral<node_T> && std::integral<edge_T>
    std::vector<ull> __calculate_min_cut_partitions_ull_ancestry_vec(
        const std::vector<edge_T> &vertices,
        const CactusMinCuts<node_T> &cactus_min_cuts,
        const std::vector<timer_type> &tin,
        const std::vector<timer_type> &tout)
    {
        size_t n_min_cuts = cactus_min_cuts.get_n_min_cuts();
        size_t n_vertices = vertices.size() - 1;
        const unsigned long B = 8 * sizeof(ull);
        const unsigned long N_COLS = (n_min_cuts + B - 1) / B;
        std::vector<ull> min_cuts = std::vector<ull>(n_vertices * N_COLS, 0ULL);
        for (node_T u{0}; u < n_vertices; u++)
        {
            size_t cut_idx{0};
            // #if defined(__AVX2__)
            if constexpr (std::is_same_v<node_T, uint32_t> && std::is_same_v<timer_type, uint32_t>)
            {
                const __m256i tin_u = _mm256_set1_epi32(static_cast<int>(tin[u]));
                const __m256i tout_u = _mm256_set1_epi32(static_cast<int>(tout[u]));
                // flip sign bit to change signed to unsigned comparison
                const __m256i sign_mask = _mm256_set1_epi32(static_cast<int>(0x80000000u));
                const __m256i tin_u_x = _mm256_xor_si256(tin_u, sign_mask);
                const __m256i tout_u_x = _mm256_xor_si256(tout_u, sign_mask);
                for (; cut_idx + 16 <= cactus_min_cuts.TREE_CUTS.size(); cut_idx += 16)
                {
                    const __m256i cuts = _mm256_loadu_si256(
                        reinterpret_cast<const __m256i *>(&cactus_min_cuts.TREE_CUTS[cut_idx]));
                    const __m256i cuts_high = _mm256_loadu_si256(
                        reinterpret_cast<const __m256i *>(&cactus_min_cuts.TREE_CUTS[cut_idx + 8]));
                    const __m256i tin_cut = _mm256_i32gather_epi32(
                        reinterpret_cast<const int *>(tin.data()),
                        cuts,
                        4);
                    const __m256i tout_cut = _mm256_i32gather_epi32(
                        reinterpret_cast<const int *>(tout.data()),
                        cuts,
                        4);
                    const __m256i tin_cut_high = _mm256_i32gather_epi32(
                        reinterpret_cast<const int *>(tin.data()),
                        cuts_high,
                        4);
                    const __m256i tout_cut_high = _mm256_i32gather_epi32(
                        reinterpret_cast<const int *>(tout.data()),
                        cuts_high,
                        4);

                    const __m256i tin_cut_x = _mm256_xor_si256(tin_cut, sign_mask);
                    const __m256i tout_cut_x = _mm256_xor_si256(tout_cut, sign_mask);
                    const __m256i tin_cut_high_x = _mm256_xor_si256(tin_cut_high, sign_mask);
                    const __m256i tout_cut_high_x = _mm256_xor_si256(tout_cut_high, sign_mask);

                    const __m256i tin_ok = _mm256_or_si256(
                        _mm256_cmpgt_epi32(tin_u_x, tin_cut_x),
                        _mm256_cmpeq_epi32(tin_u_x, tin_cut_x));
                    const __m256i tout_ok = _mm256_or_si256(
                        _mm256_cmpgt_epi32(tout_cut_x, tout_u_x),
                        _mm256_cmpeq_epi32(tout_cut_x, tout_u_x));
                    const __m256i tin_ok_high = _mm256_or_si256(
                        _mm256_cmpgt_epi32(tin_u_x, tin_cut_high_x),
                        _mm256_cmpeq_epi32(tin_u_x, tin_cut_high_x));
                    const __m256i tout_ok_high = _mm256_or_si256(
                        _mm256_cmpgt_epi32(tout_cut_high_x, tout_u_x),
                        _mm256_cmpeq_epi32(tout_cut_high_x, tout_u_x));

                    uint32_t lane_mask = static_cast<uint32_t>(
                        _mm256_movemask_ps(_mm256_castsi256_ps(_mm256_and_si256(tin_ok, tout_ok))));
                    uint32_t lane_mask_high = static_cast<uint32_t>(
                        _mm256_movemask_ps(_mm256_castsi256_ps(_mm256_and_si256(tin_ok_high, tout_ok_high))));
                    uint64_t combined_mask = (static_cast<uint64_t>(lane_mask_high) << 32) | lane_mask;
                    min_cuts[u * N_COLS + cut_idx / B] = combined_mask;
                }
            }
            // #endif

            for (; cut_idx < cactus_min_cuts.TREE_CUTS.size(); cut_idx++)
            {
                const node_T &cut = cactus_min_cuts.TREE_CUTS[cut_idx];
                if (isAncestor(tin, tout, cut, u))
                {
                    size_t word_idx = cut_idx / B;
                    size_t bit_idx = cut_idx % B;
                    min_cuts[u * N_COLS + word_idx] |= (1ULL << bit_idx);
                }
            }
            // in construction
            for (size_t cut_idx{0}; cut_idx < cactus_min_cuts.CYCLE_CUTS.size(); cut_idx++)
            {
                const auto &[cut1, cut2] = cactus_min_cuts.CYCLE_CUTS[cut_idx];
                const bool in_cut1_subtree = isAncestor(tin, tout, cut1, u);
                const bool in_cut2_subtree = isAncestor(tin, tout, cut2, u);
                if (in_cut1_subtree != in_cut2_subtree)
                {
                    size_t overall_cut_idx = cactus_min_cuts.TREE_CUTS.size() + cut_idx;
                    size_t word_idx = overall_cut_idx / B;
                    size_t bit_idx = overall_cut_idx % B;
                    min_cuts[u * N_COLS + word_idx] |= (1ULL << bit_idx);
                }
            }
        }
        return min_cuts;
    }

    // Constructs the partition matrix ull by checking ancestry relationships in the DFS tree
    // Instead of doing a BFS for each min cut
    template <class node_T, class edge_T, class weight_T>
        requires std::integral<node_T> && std::integral<edge_T>
    std::vector<ull> calculate_min_cut_partitions_ull_ancestry_vec(
        const std::vector<edge_T> &vertices,
        const std::vector<node_T> &edges,
        const std::vector<weight_T> &weights)
    {
        auto [tin, tout, cycles, is_cycle_edge, parent] = dfs_tin_tout_cycles(vertices, edges);
        weight_T min_cut = calculate_cactus_min_cut(vertices, edges, weights, cycles, is_cycle_edge);
        auto cactus_min_cuts = calculate_all_cactus_min_cuts(vertices, edges, weights, cycles, is_cycle_edge, min_cut, parent, static_cast<node_T>(0));
        return __calculate_min_cut_partitions_ull_ancestry_vec(vertices, cactus_min_cuts, tin, tout);
    }

};
