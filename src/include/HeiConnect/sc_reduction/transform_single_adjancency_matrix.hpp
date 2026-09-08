#pragma once

#include <concepts>
#include <cstddef>
#include <vector>

#include "HeiConnect/set_cover/set_cover.hpp"

using ull = unsigned long long;

namespace HeiConnect_details
{

    // haven't bothered making a char version of this one since ull is faster
    template <class node_T, class edge_T>
        requires std::integral<node_T> && std::integral<edge_T>
    SetCoverBit<> construct_set_cover_bit_matrix_ull(
        const std::vector<edge_T> &vertices,
        const std::vector<node_T> &edges,
        const std::vector<double> &weights,
        const std::vector<size_t> &link_vertices,
        const std::vector<size_t> &link_edges,
        const std::vector<double> &link_weights,
        const std::vector<ull> &min_cuts,
        size_t n_min_cuts)
    {
        const unsigned long B = 8 * sizeof(ull);
        const unsigned long N_ROWS = link_edges.size();
        const unsigned long N_COLS = n_min_cuts / B + (n_min_cuts % B != 0);
        std::vector<ull> set_cover = std::vector<ull>(N_ROWS * N_COLS, 0ULL);
        size_t current_link_idx = 0;
        for (size_t u{}; u < link_vertices.size() - 1; u++)
        {
            for (size_t e{link_vertices[u]}; e < link_vertices[u + 1]; e++)
            {
                node_T v = link_edges[e];
                for (auto w = 0; w < N_COLS; ++w)
                {
                    ull wu = min_cuts[u * N_COLS + w];
                    ull wv = min_cuts[v * N_COLS + w];
                    ull x = wu ^ wv;
                    set_cover[current_link_idx * N_COLS + w] = x;
                }
                current_link_idx++;
            }
        }
        return SetCoverBit<>{std::move(set_cover), N_ROWS, n_min_cuts, std::move(link_weights)};
    }
};
