#pragma once

#include <cstddef>
#include <filesystem>
#include <string>

namespace dataset_generator
{

enum class OutputType
{
    graph,
    links
};

enum class DistributionType
{
    constant,
    float_uniform,
    integer_uniform
};

struct Config
{
    OutputType output_type;
    std::string generator;
    std::filesystem::path output;
    std::filesystem::path input_graph;

    std::size_t nodes = 0;
    std::size_t cycles = 0;
    std::size_t cycle_length = 0;
    std::size_t min_cycle_size = 0;
    std::size_t max_cycle_size = 0;

    DistributionType distribution = DistributionType::constant;
    double constant_weight = 1.0;
    double float_uniform_lower = 0.0;
    double float_uniform_upper = 1.0;
    long long integer_uniform_lower = 1;
    long long integer_uniform_upper = 9;
    unsigned int seed = 42;
};

} // namespace dataset_generator
