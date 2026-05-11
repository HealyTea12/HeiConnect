#pragma once

#include <bit>
#include <cassert>
#include <concepts>
#include <cstddef>
#include <limits>
#include <vector>

using ull = unsigned long long;

// require T to be an unsigned integral type
template <class T>
    requires std::unsigned_integral<T>
void inline set_bit(std::vector<T> &bit_vector, size_t index) noexcept
{
    constexpr unsigned B = 8 * sizeof(T);
    bit_vector[index / (B)] |= (static_cast<T>(1) << (index % (B)));
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
    assert(false && "Edge not found: get_edge_index called with non-existent edge");
    return std::numeric_limits<EdgeID>::max(); // max index is reserved for invalid
}

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