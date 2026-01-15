#pragma once

#include <filesystem>
#include <memory>
#include "algorithm_registry.hpp"

class AlgorithmRunner
{
public:
    virtual ~AlgorithmRunner() = default;
    virtual void run(const std::filesystem::path &graph_file) = 0;
    virtual void print_results(std::ostream &os) = 0;
};

std::unique_ptr<AlgorithmRunner> create_algorithm_runner(Algorithms algorithm);
