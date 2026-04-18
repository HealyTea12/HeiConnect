#pragma once

#include <concepts>
#include <iostream>
#include <vector>

#include <omp.h>

#include "HeiConnect/set_cover/set_cover.hpp"
#include "HeiConnect/set_cover/transform_single_adjancency_matrix.hpp"
#include "HeiConnect/set_cover/transform_single_core.hpp"
#include "HeiConnect/set_cover/transform_single_csr.hpp"
#include "HeiConnect/set_cover/transform_single_oracle_ancestry.hpp"
#include "HeiConnect/set_cover/transform_single_partition_matrix.hpp"
#include "HeiConnect/set_cover/transform_single_partition_matrix_utils.hpp"

/*
 * This is the public api for the reduction from connectivity augmentation to set cover.
 */

template <typename node_T, typename edge_T>
    requires std::integral<node_T> && std::integral<edge_T>
SetCover construct_set_cover(
    const std::vector<edge_T> &vertices,
    const std::vector<node_T> &edges,
    const std::vector<double> &weights,
    const std::vector<size_t> &link_vertices,
    const std::vector<size_t> &link_edges,
    const std::vector<double> &link_weights)
{
    // generate min cut matrix
    auto [min_cuts, n_min_cuts] = HeiConnect_details::generate_min_cut_matrix(
        vertices,
        edges,
        weights,
        link_vertices,
        link_edges,
        link_weights);

    // for each link, determine which cuts it crosses
    double start = omp_get_wtime();
    auto set_cover = HeiConnect_details::construct_set_cover_csr_ull(
        vertices,
        edges,
        weights,
        link_vertices,
        link_edges,
        link_weights,
        min_cuts,
        n_min_cuts);
    double end = omp_get_wtime();
    std::cout << "Constructing set cover in CSR form took " << (end - start) << " seconds." << std::endl;
    return set_cover;
}

template <typename node_T, typename edge_T>
    requires std::integral<node_T> && std::integral<edge_T>
SetCoverBit construct_set_cover_bit_matrix(
    const std::vector<edge_T> &vertices,
    const std::vector<node_T> &edges,
    const std::vector<double> &weights,
    const std::vector<size_t> &link_vertices,
    const std::vector<size_t> &link_edges,
    const std::vector<double> &link_weights)
{
    // generate min cut matrix
    auto [min_cuts, n_min_cuts] = HeiConnect_details::generate_min_cut_matrix(
        vertices,
        edges,
        weights,
        link_vertices,
        link_edges,
        link_weights);

    // for each link, determine which cuts it crosses
    double start = omp_get_wtime();
    auto set_cover = HeiConnect_details::construct_set_cover_bit_matrix_ull(
        vertices,
        edges,
        weights,
        link_vertices,
        link_edges,
        link_weights,
        min_cuts,
        n_min_cuts);
    double end = omp_get_wtime();
    std::cout << "Constructing set cover in BIT MATRIX form took " << (end - start) << " seconds." << std::endl;
    return set_cover;
}

// Constructs "set cover" by only calculating the min cut partitions into a matrix
template <typename node_T, typename edge_T, class link_node_T, class link_edge_T>
    requires std::integral<node_T> && std::integral<edge_T>
SetCoverPseudo<link_node_T, link_edge_T> construct_set_cover_pseudo(
    const std::vector<edge_T> &vertices,
    const std::vector<node_T> &edges,
    const std::vector<double> &weights,
    const std::vector<link_edge_T> &link_vertices,
    const std::vector<link_node_T> &link_edges,
    const std::vector<double> &link_weights)
{
    auto [min_cuts, n_min_cuts] = HeiConnect_details::generate_min_cut_matrix(
        vertices,
        edges,
        weights,
        link_vertices,
        link_edges,
        link_weights);
    return {min_cuts, n_min_cuts, link_vertices, link_edges, link_weights};
}

template <class node_T, class edge_T, class weight_T, class link_node_T, class link_edge_T, class link_weight_T>
    requires std::integral<node_T> && std::integral<edge_T> && std::integral<link_node_T> && std::integral<link_edge_T>
SetCoverOracle<link_node_T, link_edge_T, link_weight_T> construct_set_cover_oracle(
    std::vector<edge_T> &vertices,
    std::vector<node_T> &edges,
    std::vector<weight_T> &weights,
    std::vector<link_node_T> &link_vertices,
    std::vector<link_edge_T> &link_edges,
    std::vector<link_weight_T> &link_weights)
{
    auto [tin, tout, cycles, is_cycle_edge, parent] = HeiConnect_details::dfs_tin_tout_cycles(vertices, edges);
    weight_T min_cut = HeiConnect_details::calculate_cactus_min_cut(vertices, edges, weights, cycles, is_cycle_edge);
    auto cactus_min_cuts = HeiConnect_details::calculate_all_cactus_min_cuts(
        vertices, edges, weights, cycles, is_cycle_edge, min_cut, parent, static_cast<node_T>(0));
    return SetCoverOracle<link_node_T, link_edge_T, link_weight_T>{tin, tout, cactus_min_cuts, link_vertices, link_edges, link_weights};
}

template <typename node_T, typename edge_T, typename weight_T, typename link_node_T, typename link_edge_T, typename link_weight_T>
    requires std::integral<node_T> && std::integral<edge_T> && std::integral<link_node_T> && std::integral<link_edge_T>
SetCoverPseudo<link_node_T, link_edge_T> construct_set_cover_pseudo_ancestry(
    const std::vector<edge_T> &vertices,
    const std::vector<node_T> &edges,
    const std::vector<weight_T> &weights,
    const std::vector<link_node_T> &link_vertices,
    const std::vector<link_edge_T> &link_edges,
    const std::vector<link_weight_T> &link_weights)
{
    auto [tin, tout, cycles, is_cycle_edge, parent] = HeiConnect_details::dfs_tin_tout_cycles(vertices, edges);
    weight_T min_cut = HeiConnect_details::calculate_cactus_min_cut(vertices, edges, weights, cycles, is_cycle_edge);
    auto cactus_min_cuts = HeiConnect_details::calculate_all_cactus_min_cuts(vertices, edges, weights, cycles, is_cycle_edge, min_cut, parent, static_cast<node_T>(0));
    auto min_cuts = HeiConnect_details::__calculate_min_cut_partitions_ull_ancestry(
        vertices,
        cactus_min_cuts,
        tin,
        tout);
    return {min_cuts, cactus_min_cuts.get_n_min_cuts(), link_vertices, link_edges, link_weights};
}

// Construct "set cover" by only calculating the parition matrix using ancestry relations instead of DFS/BFS
template <typename node_T, typename edge_T, typename weight_T, typename link_node_T, typename link_edge_T, typename link_weight_T>
    requires std::integral<node_T> && std::integral<edge_T> && std::integral<link_node_T> && std::integral<link_edge_T>
SetCoverPseudo<link_node_T, link_edge_T> construct_set_cover_pseudo_ancestry_vec(
    const std::vector<edge_T> &vertices,
    const std::vector<node_T> &edges,
    const std::vector<weight_T> &weights,
    const std::vector<link_node_T> &link_vertices,
    const std::vector<link_edge_T> &link_edges,
    const std::vector<link_weight_T> &link_weights)
{
    auto [tin, tout, cycles, is_cycle_edge, parent] = HeiConnect_details::dfs_tin_tout_cycles(vertices, edges);
    weight_T min_cut = HeiConnect_details::calculate_cactus_min_cut(vertices, edges, weights, cycles, is_cycle_edge);
    auto cactus_min_cuts = HeiConnect_details::calculate_all_cactus_min_cuts(vertices, edges, weights, cycles, is_cycle_edge, min_cut, parent, static_cast<node_T>(0));
    auto min_cuts = HeiConnect_details::__calculate_min_cut_partitions_ull_ancestry_vec(
        vertices,
        cactus_min_cuts,
        tin,
        tout);
    return {min_cuts, cactus_min_cuts.get_n_min_cuts(), link_vertices, link_edges, link_weights};
}

// TODO: this will return a new type of set cover that is a adapted for better cycle handling
// but it is not implemented yet, so for now it just outputs void
template <typename node_T, typename edge_T, typename weight_T, typename link_node_T, typename link_edge_T, typename link_weight_T>
    requires std::integral<node_T> && std::integral<edge_T> && std::integral<link_node_T> && std::integral<link_edge_T>
SetCoverCyc<link_node_T, link_edge_T, link_weight_T> construct_set_cover_cyc_pseudo_ancestry_vec(
    const std::vector<edge_T> &vertices,
    const std::vector<node_T> &edges,
    const std::vector<weight_T> &weights,
    const std::vector<link_node_T> &link_vertices,
    const std::vector<link_edge_T> &link_edges,
    const std::vector<link_weight_T> &link_weights)
{
    WeightedCRFGraph<> original_graph{{vertices, edges}, weights};
    auto [tin, tout, cycles, is_cycle_edge, parent] = HeiConnect_details::dfs_tin_tout_cycles(vertices, edges);
    weight_T min_cut = HeiConnect_details::calculate_cactus_min_cut(vertices, edges, weights, cycles, is_cycle_edge);
    std::vector<node_T> tree_cactus_min_cuts = HeiConnect_details::calculate_tree_cactus_min_cuts(vertices, edges, weights, is_cycle_edge, min_cut, parent, static_cast<node_T>(0));
    CactusMinCuts<node_T> cactus_min_cuts{tree_cactus_min_cuts, {}};
    auto partition_matrix = HeiConnect_details::__calculate_min_cut_partitions_ull_ancestry(
        vertices,
        cactus_min_cuts,
        tin,
        tout);
    auto start = omp_get_wtime();
    auto [block_tree, cycle_positions] = original_graph.cactus_generate_block_tree(0);
    auto end = omp_get_wtime();
    std::cout << "Generating block tree took " << (end - start) << " seconds." << std::endl;
    auto cycle_crosses = std::vector<std::vector<CycleCross<cycle_pos_T, cycle_id_T>>>(link_edges.size());
    /*For each link, project it into each cycle to figure out where they cross it.
     *We start by doing this with DFS, but we can improve it later
     */
    const node_T bt_n = static_cast<node_T>(block_tree.num_vertices());
    const node_T bt_root = static_cast<node_T>(0);
    std::vector<node_T> parent_bt(bt_n, bt_n);
    std::vector<size_t> depth_bt(bt_n, 0);
    std::vector<node_T> tree_stack;
    tree_stack.push_back(bt_root);
    parent_bt[bt_root] = bt_root;

    while (!tree_stack.empty())
    {
        node_T current_node = tree_stack.back();
        tree_stack.pop_back();
        for (edge_T e{block_tree.graph.vertices[current_node]}; e < block_tree.graph.vertices[current_node + 1]; e++)
        {
            node_T neighbor = block_tree.graph.edges[e];
            if (neighbor == parent_bt[current_node])
                continue;
            parent_bt[neighbor] = current_node;
            depth_bt[neighbor] = depth_bt[current_node] + 1;
            tree_stack.push_back(neighbor);
        }
    }

    start = omp_get_wtime();
    double path_time = 0.0;
    double emit_time = 0.0;
    for (node_T u{0}; u < link_vertices.size() - 1; u++)
    {
        for (edge_T e{link_vertices[u]}; e < link_vertices[u + 1]; e++)
        {
            node_T v = link_edges[e];

            auto section_start = omp_get_wtime();
            std::vector<node_T> path_u;
            std::vector<node_T> path_v;
            path_u.reserve(depth_bt[u] + 1);
            path_v.reserve(depth_bt[v] + 1);

            node_T a = u;
            node_T b = v;
            while (depth_bt[a] > depth_bt[b])
            {
                path_u.emplace_back(a);
                a = parent_bt[a];
            }
            while (depth_bt[b] > depth_bt[a])
            {
                path_v.emplace_back(b);
                b = parent_bt[b];
            }
            while (a != b)
            {
                path_u.emplace_back(a);
                path_v.emplace_back(b);
                a = parent_bt[a];
                b = parent_bt[b];
            }
            path_u.emplace_back(a);
            std::reverse(path_v.begin(), path_v.end());
            path_u.insert(path_u.end(), path_v.begin(), path_v.end());
            path_time += omp_get_wtime() - section_start;

            section_start = omp_get_wtime();
            for (size_t i{0}; i + 2 < path_u.size(); ++i)
            {
                const node_T left = path_u[i];
                const node_T middle = path_u[i + 1];
                const node_T right = path_u[i + 2];
                if (middle > original_graph.num_vertices() - 1) // cycle node
                {
                    cycle_id_T cid = static_cast<cycle_id_T>(middle - original_graph.num_vertices());

                    auto a_pos = cycle_positions[cid][left];
                    auto b_pos = cycle_positions[cid][right];
                    if (a_pos > b_pos)
                    {
                        std::swap(a_pos, b_pos);
                    }
                    cycle_crosses[e].emplace_back(CycleCross<cycle_pos_T, cycle_id_T>{
                        cid,
                        a_pos,
                        b_pos});
                }
            }
            emit_time += omp_get_wtime() - section_start;
        }
    }
    end = omp_get_wtime();
    std::cout << "Calculating cycle crosses took " << (end - start) << " seconds." << std::endl;
    std::cout << "  path reconstruction: " << path_time << " seconds." << std::endl;
    std::cout << "  emit cycle crosses: " << emit_time << " seconds." << std::endl;
    // auto cycle_coverages = std::vector<std::vector<uint32_t>>(cycles.size());
    for (size_t cycle_idx{0}; cycle_idx < cycles.size(); ++cycle_idx)
    {
        const auto cycle_len = cycles[cycle_idx].size();
        // cycle_coverages[cycle_idx] = std::vector<uint32_t>(cycle_len * cycle_len, 1);
    }
    size_t n_cycle_min_cuts = 0;
    for (size_t c{0}; c < cycles.size(); c++)
    {
        const auto cycle = cycles[c];
        const size_t cycle_len = cycle.size();
        for (size_t i{0}; i < cycle_len; i++)
        {
            node_T u = cycle[i];
            node_T v = cycle[(i + 1) % cycle_len];
            edge_T e = get_edge_index(vertices, edges, u, v);
            for (size_t j{i + 1}; j < cycle_len; j++)
            {
                node_T u2 = cycle[j];
                node_T v2 = cycle[(j + 1) % cycle_len];
                edge_T e2 = get_edge_index(vertices, edges, u2, v2);
                if (weights[e] + weights[e2] <= min_cut)
                {
                    n_cycle_min_cuts++;
                }
            }
        }
    }

    std::vector<size_t> cycle_sizes(cycles.size());
    for (size_t c{0}; c < cycles.size(); c++)
    {
        cycle_sizes[c] = cycles[c].size();
    }

    return SetCoverCyc<link_node_T, link_edge_T, link_weight_T>{
        partition_matrix,
        tree_cactus_min_cuts.size(),
        n_cycle_min_cuts,
        // cycle_coverages,
        cycle_crosses,
        cycle_positions,
        cycle_sizes,
        link_vertices,
        link_edges,
        link_weights};
}