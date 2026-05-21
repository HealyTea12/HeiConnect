#pragma once

#include <memory>
#include <filesystem>
#include "algorithm_runner.hpp"
#include "experiment_utils.hpp"

#include "HeiConnect/graph.hpp"
#include "HeiConnect/set_cover/common.hpp"
#include "HeiConnect/sc_reduction/transform_single_builders.hpp"
#include "HeiConnect/greedy.hpp"
#include "HeiConnect/set_cover/set_cover.hpp"
#include "HeiConnect/set_cover/trimmer.hpp"
#include "HeiConnect/ilp.hpp"

template <typename SetCoverType>
using GreedySetCoverPipeline = SolveTrim<
    SetCoverType,
    USSolution,
    GreedySetCoverSolver<SetCoverType, USSolution>,
    SetCoverTrimmer<SetCoverType, USSolution>>;

template <typename SetCoverType>
using CheapestSetCoverPipeline = SolveTrim<
    SetCoverType,
    USSolution,
    CheapestSetCoverSolver<SetCoverType, USSolution>,
    SetCoverTrimmer<SetCoverType, USSolution>>;

template <typename SetCoverType>
using ILPSetCoverPipeline = SolveTrim<
    SetCoverType,
    USSolution,
    SetCoverSolverILP<SetCoverType, USSolution>,
    SetCoverTrimmer<SetCoverType, USSolution>>;

// ==================== Set Cover Algorithms ====================

class SetCoverGreedySingleThreadedPQRunner : public AlgorithmRunner
{
public:
    void run(const std::filesystem::path &graph_file) override;
};

class SetCoverGreedySingleThreadedPQBitRunner : public AlgorithmRunner
{
public:
    void run(const std::filesystem::path &graph_file) override;
};

class SetCoverGreedySingleThreadedPQPseudoRunner : public AlgorithmRunner
{

public:
    void run(const std::filesystem::path &graph_file) override;
};

class SCGWCPseudoAncestryRunner : public AlgorithmRunner
{

public:
    void run(const std::filesystem::path &graph_file) override;
};

class SetCoverGreedyCheapestRunner : public AlgorithmRunner
{

public:
    void run(const std::filesystem::path &graph_file) override;
};

class SetCoverGreedyCheapestBitRunner : public AlgorithmRunner
{
public:
    void run(const std::filesystem::path &graph_file) override;
};

class SetCoverPseudoGreedyCheapestRunner : public AlgorithmRunner
{
    void run(const std::filesystem::path &graph_file) override;
};

class SetCoverILPRunner : public AlgorithmRunner
{
public:
    void run(const std::filesystem::path &graph_file) override;
};

class SetCoverPseudoILPRunner : public AlgorithmRunner
{
public:
    void run(const std::filesystem::path &graph_file) override;
};

class OracleGreedySingleThreadedPQRunner : public AlgorithmRunner
{
public:
    void run(const std::filesystem::path &graph_file) override;
};

class CycGreedySingleThreadedPQRunner : public AlgorithmRunner
{
public:
    void run(const std::filesystem::path &graph_file) override;
};

// ==================== Direct Graph Algorithms ====================

class DirectGreedyRunner : public AlgorithmRunner
{

public:
    void run(const std::filesystem::path &graph_file) override;
};

class GWCRunner : public AlgorithmRunner
{

public:
    void run(const std::filesystem::path &graph_file) override;
};

class MSTConnectRunner : public AlgorithmRunner
{

public:
    void run(const std::filesystem::path &graph_file) override;
};

class DirectILPRunner : public AlgorithmRunner
{

public:
    void run(const std::filesystem::path &graph_file) override;
};

// ==================== Helper functions ====================

inline WeightedCRFGraph<> load_weighted_crf_graph(const std::filesystem::path &graph_file)
{
    return WeightedCRFGraph<>::read_from_file_graphML(graph_file);
}

inline auto generate_unit_weight_links(const WeightedCRFGraph<> &graph)
{
    return graph.generate_links([](size_t u, size_t v)
                                { return 1.0; });
}

inline graph::GraphPair load_graph_pair(const std::filesystem::path &xml_file)
{
    graph::GraphPair g;
    g.read_graph(xml_file.parent_path() / (xml_file.stem().string() + ".graph"), xml_file);
    g.add_links(xml_file.parent_path() / (xml_file.stem().string() + ".links"), 1.f, 0);
    return g;
}

inline double calculate_edge_list_cost(const auto &edge_list)
{
    double cost = 0.0;
    for (const auto &edge : edge_list)
    {
        cost += edge.weight;
    }
    return cost;
}
