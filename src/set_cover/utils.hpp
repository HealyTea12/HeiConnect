#pragma once

#include <vector>
#include <limits>
#include <cstddef>
#include <type_traits>
#include <algorithm>
#include <omp.h>

template <typename T>
concept MinMax = requires(T a, T b) {
    { a < b } -> std::convertible_to<bool>;
    { a > b } -> std::convertible_to<bool>;
    { std::numeric_limits<T>::lowest() } -> std::same_as<T>;
    { std::numeric_limits<T>::max() } -> std::same_as<T>;
};

template <typename T>
    requires MinMax<T>
size_t argmin(const std::vector<T> &v)
{
    if (v.empty())
        throw std::invalid_argument("argmin: input vector must not be empty");

    const size_t n = v.size();
    const int num_threads = omp_get_max_threads();
    std::vector<T> local_min(num_threads, std::numeric_limits<T>::max());
    std::vector<size_t> local_min_indices(num_threads, 0);

#pragma omp parallel
    {
        const int tid = omp_get_thread_num();
#pragma omp for
        for (size_t i = 0; i < n; ++i)
        {
            if (v[i] < local_min[tid])
            {
                local_min[tid] = v[i];
                local_min_indices[tid] = i;
            }
        }
    }

    // find global min among thread-local minima
    T best = local_min[0];
    size_t best_idx = local_min_indices[0];
    for (int t = 1; t < num_threads; ++t)
    {
        if (local_min[t] < best)
        {
            best = local_min[t];
            best_idx = local_min_indices[t];
        }
    }
    return best_idx;
}

template <typename T>
    requires MinMax<T>
T max(const std::vector<T> &v)
{
    if (v.empty())
        throw std::invalid_argument("max: input vector must not be empty");

    // Use a thread-local maximum and combine them safely
    T global_max = std::numeric_limits<T>::lowest();

#pragma omp parallel
    {
        T thread_max = std::numeric_limits<T>::lowest();
#pragma omp for nowait
        for (size_t i = 0; i < v.size(); ++i)
        {
            thread_max = std::max(thread_max, v[i]);
        }
#pragma omp critical
        {
            global_max = std::max(global_max, thread_max);
        }
    }

    return global_max;
}

// argmax: return index of maximum element (parallelized)
template <typename T>
    requires MinMax<T>
size_t argmax(const std::vector<T> &v)
{
    if (v.empty())
        throw std::invalid_argument("argmax: input vector must not be empty");

    const size_t n = v.size();
    const int num_threads = omp_get_max_threads();
    std::vector<T> local_max(num_threads, std::numeric_limits<T>::lowest());
    std::vector<size_t> local_max_indices(num_threads, 0);

#pragma omp parallel
    {
        const int tid = omp_get_thread_num();
#pragma omp for
        for (size_t i = 0; i < n; ++i)
        {
            if (v[i] > local_max[tid])
            {
                local_max[tid] = v[i];
                local_max_indices[tid] = i;
            }
        }
    }

    // find global max among thread-local maxima
    T best = local_max[0];
    size_t best_idx = local_max_indices[0];
    for (int t = 1; t < num_threads; ++t)
    {
        if (local_max[t] > best)
        {
            best = local_max[t];
            best_idx = local_max_indices[t];
        }
    }
    return best_idx;
}