#pragma once

#include <bit>
#include <concepts>
#include <cstddef>
#include <iostream>
#include <vector>

#include <omp.h>

#include "HeiConnect/set_cover/set_cover.hpp"

using ull = unsigned long long;

namespace HeiConnect_details
{
    template <typename node_T, typename link_edge_T>
        requires std::unsigned_integral<node_T> && std::unsigned_integral<link_edge_T>
    size_t calculate_size_b_set_char(
        const std::vector<unsigned char> &min_cuts,
        const std::vector<link_edge_T> &link_vertices,
        const std::vector<node_T> &link_edges)
    {
        size_t size_b_set = 0;
        size_t n_min_cuts = min_cuts.size() / (link_vertices.size() - 1);
        for (node_T u{}; u < link_vertices.size() - 1; u++)
        {
            for (link_edge_T e{link_vertices[u]}; e < link_vertices[u + 1]; e++)
            {
                node_T v = link_edges[e];
                for (size_t k{}; k < n_min_cuts; k++)
                {
                    if (min_cuts[u * n_min_cuts + k] ^ min_cuts[v * n_min_cuts + k])
                    {
                        size_b_set++;
                    }
                }
            }
        }
        return size_b_set;
    }

    template <typename node_T, typename link_edge_T>
        requires std::unsigned_integral<node_T> && std::unsigned_integral<link_edge_T>
    size_t calculate_size_b_set_ull(
        const std::vector<ull> &min_cuts,
        const std::vector<link_edge_T> &link_vertices,
        const std::vector<node_T> &link_edges)
    {
        size_t size_b_set = 0;
        constexpr const size_t B = 8 * sizeof(ull);
        size_t words_per_row = min_cuts.size() / (link_vertices.size() - 1);
        for (node_T u{}; u < link_vertices.size() - 1; u++)
        {
            for (link_edge_T e{link_vertices[u]}; e < link_vertices[u + 1]; e++)
            {
                node_T v = link_edges[e];
                for (ull word_i = 0; word_i < words_per_row; ++word_i)
                {
                    ull wu = min_cuts[u * words_per_row + word_i];
                    ull wv = min_cuts[v * words_per_row + word_i];
                    ull x = wu ^ wv;
                    size_t count = std::popcount(x);
                    size_b_set += count;
                }
            }
        }
        return size_b_set;
    }

    template <class node_T, class edge_T>
        requires std::integral<node_T> && std::integral<edge_T>
    SetCover<> construct_set_cover_csr_b_size_ull( // with ull min_cuts
        const std::vector<edge_T> &vertices,
        const std::vector<node_T> &edges,
        const std::vector<double> &weights,
        const std::vector<size_t> &link_vertices,
        const std::vector<size_t> &link_edges,
        const std::vector<double> &link_weights,
        const std::vector<ull> &min_cuts,
        size_t n_min_cuts)
    {
        std::cout << n_min_cuts << " min cuts found." << std::endl;
        auto start = omp_get_wtime();
        size_t size_of_b = calculate_size_b_set_ull(min_cuts, link_vertices, link_edges);
        std::vector<size_t> a = std::vector<size_t>(link_edges.size() + 1, static_cast<size_t>(0));
        std::vector<size_t> b = std::vector<size_t>(size_of_b, 0ULL);
        // b.reserve(link_edges.size() * n_min_cuts / 2);
        size_t n_vertices = vertices.size() - 1;
        constexpr const size_t B = 8 * sizeof(ull);
        auto end = omp_get_wtime();
        std::cout << "Preprocessing time before constructing set cover: " << (end - start) << " seconds." << std::endl;
        start = omp_get_wtime();
        size_t current_idx = 0ULL;
        size_t count = 0ULL;
        for (size_t u{}; u < link_vertices.size() - 1; u++)
        {
            for (size_t e{link_vertices[u]}; e < link_vertices[u + 1]; e++)
            {
                // auto link = Edge{u, link_edges[e], link_weights[e]};
                node_T v = link_edges[e];
                const size_t words_per_vertex = min_cuts.size() / (link_vertices.size() - 1);
                for (size_t word_i = 0; word_i < words_per_vertex; ++word_i)
                {
                    ull wu = min_cuts[u * words_per_vertex + word_i];
                    ull wv = min_cuts[v * words_per_vertex + word_i];
                    ull x = wu ^ wv;
                    while (x)
                    {
                        int bit = std::countr_zero(x);
                        ull mask = (1ULL << bit);
                        size_t k = word_i * B + bit;
                        b[current_idx++] = k;
                        x ^= mask;
                        count++;
                    }
                }
                a[e + 1] = count;
            }
        }
        end = omp_get_wtime();
        std::cout << "Number of links: " << link_edges.size() << std::endl;
        std::cout << "Total time constructing set cover: " << (end - start) << " seconds." << std::endl;
        std::cout << "Average time processing link: " << (end - start) / link_edges.size() << " seconds." << std::endl;
        std::cout << "Size of b set: " << b.size() << std::endl;
        return SetCover<>{std::move(a), std::move(b), std::move(link_weights)};
    }

    template <class node_T, class edge_T>
        requires std::integral<node_T> && std::integral<edge_T>
    SetCover<> construct_set_cover_csr_char( // with char min_cuts
        const std::vector<edge_T> &vertices,
        const std::vector<node_T> &edges,
        const std::vector<double> &weights,
        const std::vector<size_t> &link_vertices,
        const std::vector<size_t> &link_edges,
        const std::vector<double> &link_weights,
        const std::vector<unsigned char> &min_cuts,
        size_t n_min_cuts)
    {
        std::cout << n_min_cuts << " min cuts found." << std::endl;
        auto start = omp_get_wtime();
        std::vector<size_t> a = std::vector<size_t>(link_edges.size() + 1, static_cast<size_t>(0));
        std::vector<size_t> b{};
        // b.reserve(link_edges.size() * n_min_cuts / 2);
        size_t n_vertices = vertices.size() - 1;
        constexpr size_t B = 8 * sizeof(ull);
        auto end = omp_get_wtime();
        std::cout << "Preprocessing time before constructing set cover: " << (end - start) << " seconds." << std::endl;
        start = omp_get_wtime();
        for (size_t u{0}; u < link_vertices.size() - 1; u++)
        {
            for (size_t e{link_vertices[u]}; e < link_vertices[u + 1]; e++)
            {
                // auto link = Edge{u, link_edges[e], link_weights[e]};
                node_T v = link_edges[e];
                for (size_t k{}; k < n_min_cuts; k++)
                {
                    if (min_cuts[u * n_min_cuts + k] ^ min_cuts[v * n_min_cuts + k])
                    {
                        b.emplace_back(k);
                    }
                }
                a[e + 1] = b.size();
            }
        }
        end = omp_get_wtime();
        std::cout << "Number of links: " << link_edges.size() << std::endl;
        std::cout << "Total time constructing set cover: " << (end - start) << " seconds." << std::endl;
        std::cout << "Average time processing link: " << (end - start) / link_edges.size() << " seconds." << std::endl;
        return SetCover<>{std::move(a), std::move(b), std::move(link_weights)};
    }

    template <class node_T, class edge_T>
        requires std::integral<node_T> && std::integral<edge_T>
    SetCover<> construct_set_cover_csr_ull( // with ull min_cuts
        const std::vector<edge_T> &vertices,
        const std::vector<node_T> &edges,
        const std::vector<double> &weights,
        const std::vector<size_t> &link_vertices,
        const std::vector<size_t> &link_edges,
        const std::vector<double> &link_weights,
        const std::vector<ull> &min_cuts,
        size_t n_min_cuts)
    {
        auto start = omp_get_wtime();
        std::vector<size_t> a = std::vector<size_t>(link_edges.size() + 1, static_cast<size_t>(0));
        std::vector<size_t> b = std::vector<size_t>(0ULL);
        b.reserve(link_edges.size() * n_min_cuts / 2);
        size_t n_vertices = vertices.size() - 1;
        constexpr const size_t B = 8 * sizeof(ull);
        auto end = omp_get_wtime();
        std::cout << "Preprocessing time before constructing set cover: " << (end - start) << " seconds." << std::endl;
        start = omp_get_wtime();
        for (size_t u{}; u < link_vertices.size() - 1; u++)
        {
            for (size_t e{link_vertices[u]}; e < link_vertices[u + 1]; e++)
            {
                // auto link = Edge{u, link_edges[e], link_weights[e]};
                node_T v = link_edges[e];
                const size_t words_per_vertex = min_cuts.size() / (link_vertices.size() - 1);
                for (size_t word_i = 0; word_i < words_per_vertex; ++word_i)
                {
                    ull wu = min_cuts[u * words_per_vertex + word_i];
                    ull wv = min_cuts[v * words_per_vertex + word_i];
                    ull x = wu ^ wv;
                    while (x)
                    {
                        int bit = std::countr_zero(x);
                        ull mask = (1ULL << bit);
                        size_t k = word_i * B + bit;
                        b.emplace_back(k);
                        x ^= mask;
                    }
                }
                a[e + 1] = b.size();
            }
        }
        end = omp_get_wtime();
        std::cout << "Number of links: " << link_edges.size() << std::endl;
        std::cout << "Total time constructing set cover: " << (end - start) << " seconds." << std::endl;
        std::cout << "Average time processing link: " << (end - start) / link_edges.size() << " seconds." << std::endl;
        std::cout << "Size of b set: " << b.size() << std::endl;
        return SetCover<>{std::move(a), std::move(b), std::move(link_weights), n_min_cuts};
    }
};