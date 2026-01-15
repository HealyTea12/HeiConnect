#pragma once

#include <string>
#include <array>
#include <stdexcept>

enum class Algorithms
{
    SetCoverGreedySingleThreadedPQ,
    SetCoverGreedySingleThreadedPQBit,
    SetCoverGreedySingleThreadedPQPseudo,
    SetCoverSharpGreedy,
    SetCoverGreedyCheapest,
    SetCoverILP,
    DirectGreedy,
    GWC,
    MSTConnect,
    DirectILP
};

constexpr std::array<std::string_view, 10> ALGORITHM_NAMES = {
    "SetCoverGreedySingleThreadedPQ",
    "SetCoverGreedySingleThreadedPQBit",
    "SetCoverGreedySingleThreadedPQPseudo",
    "SetCoverSharpGreedy",
    "SetCoverGreedyCheapest",
    "SetCoverILP",
    "DirectGreedy",
    "GWC",
    "MSTConnect",
    "DirectILP"};

inline Algorithms algorithm_from_string(const std::string &algo)
{
    for (size_t i = 0; i < ALGORITHM_NAMES.size(); ++i)
    {
        if (ALGORITHM_NAMES[i] == algo)
        {
            return static_cast<Algorithms>(i);
        }
    }
    throw std::invalid_argument("Unknown algorithm: " + algo);
}

inline std::string algorithm_to_string(Algorithms algo)
{
    return std::string(ALGORITHM_NAMES[static_cast<size_t>(algo)]);
}
