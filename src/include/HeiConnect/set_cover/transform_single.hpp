#pragma once
#include <concepts>
#include <vector>
#include <queue>
#include <unordered_set>
#include <cassert>
#include <omp.h>

#include "HeiConnect/min_cut/simple_mincut.hpp"
#include "HeiConnect/set_cover/set_cover.hpp"
#include "HeiConnect/bfs.hpp"

using ull = unsigned long long;

// require T to be an unsigned integral type
template <class T>
    requires std::unsigned_integral<T>
void inline set_bit(std::vector<T> &bit_vector, size_t index) noexcept
{
    bit_vector[index / (8 * sizeof(T))] |= (static_cast<T>(1) << (index % (8 * sizeof(T))));
}

template <class T>
    requires std::unsigned_integral<T>
bool inline constexpr get_bit(const std::vector<T> &bit_vector, size_t index)
{
    // constexpr unsigned SHIFT = std::countr_zero(8 * sizeof(T));                    // C++20
    // constexpr T MASK = (static_cast<T>(1) << (8 * sizeof(T))) - static_cast<T>(1); // not used directly
    // size_t word = index >> SHIFT;
    // unsigned offset = index & MASK;
    // return (bit_vector[word] >> offset) & T(1);
    constexpr unsigned B = 8 * sizeof(T);
    return bit_vector[index / B] & (static_cast<T>(1) << (index % B));
}

// returns the edge idx between u and v.
// Assumes that there is an edge between u and v.
// If there is no such edge, assert false.
// Could consider adding an optional return type instead.
// But, in this context, if the edge doesn't exist, then there is a bug in the code.
template <typename NodeID, typename EdgeID>
inline EdgeID get_edge_index(
    const std::vector<EdgeID> &vertices,
    const std::vector<NodeID> &edges,
    NodeID u,
    NodeID v) noexcept
{
    for (auto i{vertices[u]}; i < vertices[u + 1]; i++)
    {
        if (edges[i] == v)
        {
            return i;
        }
    }
    return std::numeric_limits<EdgeID>::max(); // max index is reserved for invalid
}

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
    std::vector<ull> min_cuts{};
    size_t n_vertices = vertices.size() - 1;
    constexpr const size_t B = 8 * sizeof(ull);
    // const size_t N_COLS = (n_vertices - 1) / (8 * sizeof(ull)) + 1;    // 64 vertices -> (64 -1)/64 + 1
    // const size_t MAX_NUM_MIN_CUTS = n_vertices * (n_vertices - 1) / 2; // n choose 2 in case of cycle, the number of rows we need
    const size_t N_COLS = (n_min_cuts - 1) / B + 1;
    min_cuts.resize(n_vertices * N_COLS, 0); // 1 extra
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
SetCover construct_set_cover_csr_b_size_ull( // with ull min_cuts
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
    return SetCover{std::move(a), std::move(b), std::move(link_weights)};
}

template <class node_T, class edge_T>
    requires std::integral<node_T> && std::integral<edge_T>
SetCover construct_set_cover_csr_ull( // with ull min_cuts
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
    return SetCover{std::move(a), std::move(b), std::move(link_weights), n_min_cuts};
}

// haven't bothered making a char version of this one since ull is faster
template <class node_T, class edge_T>
    requires std::integral<node_T> && std::integral<edge_T>
SetCoverBit construct_set_cover_bit_matrix_ull(
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
    const unsigned long N_COLS = (n_min_cuts - 1) / B + 1;
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
    return SetCoverBit{std::move(set_cover), N_ROWS, n_min_cuts, std::move(link_weights)};
}

template <class node_T, class edge_T>
    requires std::integral<node_T> && std::integral<edge_T>
SetCover construct_set_cover_csr_char( // with char min_cuts
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
    return SetCover{std::move(a), std::move(b), std::move(link_weights)};
}

namespace details
{
    template <typename node_T, typename edge_T>
        requires std::unsigned_integral<node_T> && std::unsigned_integral<edge_T>
    struct Info
    {
        node_T parent;
        edge_T parent_edge_index; // from parent to v
        int distance;
    };
};

template <typename node_T, typename edge_T>
    requires std::integral<node_T> && std::integral<edge_T>
static std::tuple<
    std::vector<std::pair<node_T, node_T>>,
    std::vector<details::Info<node_T, edge_T>>>
find_cycles_and_root_tree(
    const std::vector<edge_T> &vertices,
    const std::vector<node_T> &edges,
    const std::vector<double> &weights,
    const std::vector<size_t> &link_vertices,
    const std::vector<size_t> &link_edges,
    const std::vector<double> &link_weights)
{
    std::vector<details::Info<node_T, edge_T>> info = std::vector<details::Info<node_T, edge_T>>(vertices.size() - 1);
    std::vector<bool> visited_nodes = std::vector<bool>(vertices.size() - 1, false);
    // store parent as (parent_node, edge_index) edge_index is for the edge going from v to parent(v)
    // otherwise, we would have to do a linear search over the neighbours to find the edge index
    // when reconstructing the cycles
    std::queue<std::pair<node_T, int>> q;
    node_T root = static_cast<node_T>(0);
    info[root] = {root, std::numeric_limits<edge_T>::max(), 0};
    visited_nodes[root] = true;
    // parent, child
    std::vector<std::pair<node_T, node_T>> cycle_edge_vec{};
    cycle_edge_vec.reserve(edges.size());
    q.push({root, 0});
    while (!q.empty())
    {
        auto curr = q.front();
        node_T current_node = curr.first;
        auto depth = curr.second;
        q.pop();
        for (edge_T i{vertices[current_node]}; i < vertices[current_node + 1]; i++)
        {
            node_T neighbor = edges[i];
            if (!visited_nodes[neighbor])
            {
                visited_nodes[neighbor] = true;
                info[neighbor] = {current_node, i, depth + 1};
                q.push({neighbor, depth + 1});
            }
            else
            {
                if (info[neighbor].distance < depth)
                    continue; // skip back edges to ancestors
                cycle_edge_vec.emplace_back(current_node, neighbor);
            }
        }
    }
    return std::tuple<
        std::vector<std::pair<node_T, node_T>>,
        std::vector<details::Info<node_T, edge_T>>>{cycle_edge_vec, info};
}

// node_T/edge_T is a generic type that indexes nodes/edges e.g unsigned int
template <typename node_T, typename weight_T>
struct Edge
{
    node_T u;
    node_T v;
    weight_T weight;
};

// template <typename G>
// concept graph = requires(G &graph) {
//     graph.num_nodes()->std::integral_type
// };

struct pair_hash
{
    inline std::size_t operator()(const std::pair<int, int> &v) const
    {
        return v.first * 31 + v.second;
    }
};

template <class T>
    requires std::unsigned_integral<T>
std::vector<T> transpose(const std::vector<T> &input, size_t n_rows, size_t n_cols)
{
    assert(input.size() == n_rows * n_cols);
    constexpr size_t B = 8 * sizeof(T);
    // number of output columns (words) needed to store n_rows bits per column
    const size_t N_COLS = (n_rows + B - 1) / B;
    // number of output rows (one per input bit-column)
    const size_t N_ROWS = n_cols * B;
    std::vector<T> output = std::vector<T>(N_ROWS * N_COLS, static_cast<T>(0));
    for (size_t r = 0; r < n_rows; r++)
    {
        for (size_t c = 0; c < n_cols; c++)
        {
            for (size_t bit = 0; bit < B; bit++)
            {
                size_t in_bit_idx = (r * n_cols + c) * B + bit;
                if (get_bit(input, in_bit_idx))
                {
                    // output row corresponds to the input bit-column
                    size_t out_row = c * B + bit; // 0..N_ROWS-1
                    size_t out_word = out_row * N_COLS + (r / B);
                    size_t out_bit = r % B;
                    size_t out_bit_idx = out_word * B + out_bit;
                    set_bit(output, out_bit_idx);
                }
            }
        }
    }
    return output;
}

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
    for (auto cyc : cycles)
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
SetCover construct_set_cover(
    const std::vector<edge_T> &vertices,
    const std::vector<node_T> &edges,
    const std::vector<double> &weights,
    const std::vector<size_t> &link_vertices,
    const std::vector<size_t> &link_edges,
    const std::vector<double> &link_weights)
{
    // generate min cut matrix
    auto [min_cuts, n_min_cuts] = generate_min_cut_matrix(
        vertices,
        edges,
        weights,
        link_vertices,
        link_edges,
        link_weights);

    // for each link, determine which cuts it crosses
    double start = omp_get_wtime();
    auto set_cover = construct_set_cover_csr_ull(
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
    auto [min_cuts, n_min_cuts] = generate_min_cut_matrix(
        vertices,
        edges,
        weights,
        link_vertices,
        link_edges,
        link_weights);

    // for each link, determine which cuts it crosses
    double start = omp_get_wtime();
    auto set_cover = construct_set_cover_bit_matrix_ull(
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
template <typename node_T, typename edge_T>
    requires std::integral<node_T> && std::integral<edge_T>
SetCoverPseudo construct_set_cover_pseudo(
    const std::vector<edge_T> &vertices,
    const std::vector<node_T> &edges,
    const std::vector<double> &weights,
    const std::vector<size_t> &link_vertices,
    const std::vector<size_t> &link_edges,
    const std::vector<double> &link_weights)
{
    auto [min_cuts, n_min_cuts] = generate_min_cut_matrix(
        vertices,
        edges,
        weights,
        link_vertices,
        link_edges,
        link_weights);
    return SetCoverPseudo{min_cuts, n_min_cuts, link_vertices, link_edges, link_weights};
}
// template <typename node_T, typename edge_T>
//     requires std::integral<node_T> && std::integral<edge_T>
// SetCover construct_set_cover1(
//     const std::vector<edge_T> &vertices,
//     const std::vector<node_T> &edges,
//     const std::vector<double> &weights,
//     const std::vector<size_t> &link_vertices,
//     const std::vector<size_t> &link_edges,
//     const std::vector<double> &link_weights)
// {
//     // auto cycle_edges = std::unordered_set<std::pair<node_T, node_T>>{0, pair_hash};
//     // for checking during the second BFS
//     auto cycle_edges = std::vector<bool>(edges.size(), false);
//     // BFS to find cycles and rooted tree
//     std::vector<std::pair<node_T, node_T>> cycle_edge_vec{};
//     std::vector<find_cycles_and_root_tree::Info> info{};
//     auto [cycle_edge_vec, info] = find_cycles_and_root_tree(vertices, edges, weights, link_vertices, link_edges, link_weights);
//     // for each cycle edge, reconstruct the cycle
//     // could rethink data struct
//     auto cycles = std::vector<std::vector<std::pair<node_T, node_T>>>(cycle_edge_vec.size());
//     for (size_t i{}; i < cycle_edge_vec.size(); i++)
//     {
//         auto &cycle = cycles[i];
//         cycle.emplace_back(cycle_edge_vec[i]);
//         node_T u = cycle_edge_vec[i].first;
//         node_T v = cycle_edge_vec[i].second;
//         // mark only one direction, as we only traverse edges in this direction during min cut
//         cycle_edges[info[v].parent_edge_index] = true;
//         if (info[u].distance < info[v].distance)
//         {
//             cycle_edges[info[v].parent_edge_index] = true;
//             v = info[v].parent_node;
//         }
//         while (u != v)
//         {
//             cycle.emplace_back(info[u].parent_node, u);
//             cycle.emplace_back(info[v].parent_node, v);
//             cycle_edges[info[u].parent_edge_index] = true;
//             cycle_edges[info[v].parent_edge_index] = true;
//             u = info[u].parent_node;
//             v = info[v].parent_node;
//         }
//     }
//
//     // iterate over all edges to find min cuts
//     // auto min_cut = global_mincut_simple({{vertices, edges},
//     //                                     weights}); // should substitute with cactus min cut
//     // this part is awful
//     auto min_cut = std::numeric_limits<double>::max();
//     bfs_single_threaded(
//         vertices,
//         edges,
//         static_cast<node_T>(0),
//         [&cycle_edges, &info, &weights, &min_cut](node_T from, node_T to) noexcept
//         {
//             if (!cycle_edges[info[to].parent_edge_index] && weights[info[to].parent_edge_index] < min_cut)
//             {
//                 min_cut = weights[info[to].parent_edge_index];
//             }
//         },
//         [&cycle_edges, &info, &weights, &min_cut](node_T from, node_T to) noexcept
//         {
//             if (!cycle_edges[info[to].parent_edge_index] && weights[info[to].parent_edge_index] < min_cut)
//             {
//                 min_cut = weights[info[to].parent_edge_index];
//             }
//         },
//         [](node_T from, node_T to) noexcept
//         {
//             return false;
//         });
//     for (const auto &cycle : cycles)
//     {
//         for (auto i = 0; i < cycle.size() - 1; i++)
//         {
//             auto edge1 = cycle[i];
//             auto edge1_idx = info[edge1.second].parent_edge_index;
//             for (auto j = i + 1; j < cycle.size(); j++)
//             {
//                 auto edge2 = cycle[j];
//                 auto edge_idx = info[edge2.second].parent_edge_index;
//                 double cycle_cut = weights[edge1_idx] + weights[edge_idx];
//                 if (cycle_cut < min_cut)
//                 {
//                     min_cut = cycle_cut;
//                 }
//             }
//         }
//     }
//
//     // wrong
//     std::vector<uint8_t> min_cuts{};
//     size_t n_min_cuts{};
//     size_t n_vertices = vertices.size() - 1;
//     min_cuts.resize(n_vertices * n_vertices * (n_vertices - 1) / (2 * 8), 0);
//     for (node_T u{}; u < vertices.size() - 1; u++)
//     {
//         for (edge_T i{vertices[u]}; i < vertices[u + 1]; i++)
//         {
//             node_T v = edges[i];
//             if (u < v && !cycle_edges[i])
//             {
//                 if (weights[i] == min_cut)
//                 {
//                     auto i = n_min_cuts * n_vertices + u;
//                     min_cuts[i / 8] &= 1 << (i % 8);
//                     bfs_single_threaded(
//                         vertices,
//                         edges,
//                         u,
//                         [](node_T from, node_T to) noexcept {},
//                         [&min_cuts, n_min_cuts, n_vertices](node_T from, node_T to) noexcept
//                         {
//                             auto i = n_min_cuts * n_vertices + to;
//                             min_cuts[i / 8] &= 1 << (i % 8);
//                         },
//                         [u, v](node_T from, node_T to) noexcept
//                         {
//                             if (from == u && to == v)
//                                 return true;
//                             return false;
//                         });
//                     n_min_cuts++;
//                 }
//             }
//         }
//     }
//     for (const auto &cycle : cycles)
//     {
//         for (size_t i{}; i < cycle.size() - 1; i++)
//         {
//             for (size_t j = i + 1; j < cycle.size(); j++)
//             {
//                 auto edge1 = cycle[i];
//                 auto edge2 = cycle[j];
//                 auto edge1_idx = get_edge_index(vertices, edges, edge1.first, edge1.second);
//                 auto edge2_idx = get_edge_index(vertices, edges, edge2.first, edge2.second);
//                 double cycle_cut = weights[edge1_idx] + weights[edge2_idx];
//                 if (cycle_cut == min_cut)
//                 {
//                     auto i = n_min_cuts * n_vertices + edge1.first;
//                     min_cuts[i / 8] &= 1 << (i % 8);
//                     bfs_single_threaded(
//                         vertices,
//                         edges,
//                         edge1.first,
//                         [](node_T from, node_T to) noexcept {},
//                         [&min_cuts, n_min_cuts, n_vertices](node_T from, node_T to) noexcept
//                         {
//                             auto i = n_min_cuts * n_vertices + to;
//                             min_cuts[i / 8] &= 1 << (i % 8);
//                         },
//                         [&edge1, &edge2](node_T from, node_T to) noexcept
//                         {
//                             if ((from == edge1.first && to == edge1.second) ||
//                                 (from == edge2.first && to == edge2.second) ||
//                                 (from == edge2.second && to == edge2.first))
//                                 return true;
//                             return false;
//                         });
//                     n_min_cuts++;
//                 }
//             }
//         }
//     }
//     // for each link, determine which cuts it crosses
//     std::vector<size_t> a = std::vector<size_t>(link_edges.size() + 1, static_cast<size_t>(0));
//     std::vector<size_t> b{};
//     b.reserve(link_edges.size() * 3);
//     for (size_t i{}; i < link_vertices.size() - 1; i++)
//     {
//         for (size_t j{link_vertices[i]}; j < link_vertices[i + 1]; j++)
//         {
//             auto link = Edge{i, link_edges[j], link_weights[j]};
//             for (size_t k{}; k < n_min_cuts; k++)
//             {
//                 auto u_idx = k * n_vertices + link.u;
//                 auto v_idx = k * n_vertices + link.v;
//                 auto cond = (min_cuts[u_idx / 8] & (1 << (u_idx % 8))) ^ (min_cuts[v_idx / 8] & (1 << (v_idx % 8)));
//                 if (cond)
//                 {
//                     b.emplace_back(k);
//                 }
//             }
//             a[j + 1] = b.size();
//         }
//     }
//     return SetCover{a, b, link_weights};
// }