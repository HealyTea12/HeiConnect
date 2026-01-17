#pragma once

#include <memory>
#include <filesystem>
#include "algorithm_runner.hpp"
#include "experiment_utils.hpp"

#include "HeiConnect/graph.hpp"
#include "HeiConnect/set_cover/transform_single.hpp"
#include "HeiConnect/greedy.hpp"
#include "HeiConnect/set_cover/set_cover.hpp"
#include "HeiConnect/ilp.hpp"

// ==================== Set Cover Algorithms ====================

class SetCoverGreedySingleThreadedPQRunner : public AlgorithmRunner
{
public:
    struct Result
    {
        double solution_cost;
        size_t solution_size;
        double time_reduction;
        double time_solving;
        double time_total;
    } result;
    void run(const std::filesystem::path &graph_file) override;
    void print_results(std::ostream &os) override
    {
        os << "Solution cost: " << result.solution_cost << "\n";
        os << "Solution size: " << result.solution_size << "\n";
        os << "Time reduction (s): " << result.time_reduction << "\n";
        os << "Time solving (s): " << result.time_solving << "\n";
        os << "Time total (s): " << result.time_total << "\n";
    }
};

class SetCoverGreedySingleThreadedPQBitRunner : public AlgorithmRunner
{
    struct Result
    {
        double solution_cost;
        size_t solution_size;
        double time_reduction;
        double time_solving;
        double time_total;
    } result;

public:
    void run(const std::filesystem::path &graph_file) override;
    void print_results(std::ostream &os) override
    {
        os << "Solution cost: " << result.solution_cost << "\n";
        os << "Solution size: " << result.solution_size << "\n";
        os << "Time reduction (s): " << result.time_reduction << "\n";
        os << "Time solving (s): " << result.time_solving << "\n";
        os << "Time total (s): " << result.time_total << "\n";
    }
};

class SetCoverGreedySingleThreadedPQPseudoRunner : public AlgorithmRunner
{
    struct Result
    {
        double solution_cost;
        size_t solution_size;
        double time_reduction;
        double time_solving;
        double time_total;
    } result;

public:
    void run(const std::filesystem::path &graph_file) override;
    void print_results(std::ostream &os) override
    {
        os << "Solution cost: " << result.solution_cost << "\n";
        os << "Solution size: " << result.solution_size << "\n";
        os << "Time reduction (s): " << result.time_reduction << "\n";
        os << "Time solving (s): " << result.time_solving << "\n";
        os << "Time total (s): " << result.time_total << "\n";
    }
};

class SetCoverGreedyCheapestRunner : public AlgorithmRunner
{
    struct Result
    {
        double solution_cost;
        size_t solution_size;
        double solution_cost_trimmed;
        size_t solution_size_trimmed;
        double time_reduction;
        double time_solving;
        double time_trimming;
        double time_total;
    } result;

public:
    void run(const std::filesystem::path &graph_file) override;
    void print_results(std::ostream &os) override
    {
        os << "Solution cost: " << result.solution_cost << "\n";
        os << "Solution size: " << result.solution_size << "\n";
        os << "Solution cost (trimmed): " << result.solution_cost_trimmed << "\n";
        os << "Solution size (trimmed): " << result.solution_size_trimmed << "\n";
        os << "Time reduction (s): " << result.time_reduction << "\n";
        os << "Time solving (s): " << result.time_solving << "\n";
        os << "Time trimming (s): " << result.time_trimming << "\n";
        os << "Time total (s): " << result.time_total << "\n";
    }
};

class SetCoverPseudoGreedyCheapestRunner : public AlgorithmRunner
{
    struct Result
    {
        double solution_cost;
        size_t solution_size;
        double solution_cost_trimmed;
        size_t solution_size_trimmed;
        double time_reduction;
        double time_solving;
        double time_trimming;
        double time_total;
    } result;
    void run(const std::filesystem::path &graph_file) override;
    void print_results(std::ostream &os) override
    {
        os << "Solution cost: " << result.solution_cost << "\n";
        os << "Solution size: " << result.solution_size << "\n";
        os << "Solution cost (trimmed): " << result.solution_cost_trimmed << "\n";
        os << "Solution size (trimmed): " << result.solution_size_trimmed << "\n";
        os << "Time reduction (s): " << result.time_reduction << "\n";
        os << "Time solving (s): " << result.time_solving << "\n";
        os << "Time trimming (s): " << result.time_trimming << "\n";
        os << "Time total (s): " << result.time_total << "\n";
    }
};

class SetCoverILPRunner : public AlgorithmRunner
{
public:
    struct Result
    {
        double solution_cost;
        size_t solution_size;
        double time_reduction;
        double time_solving;
        double time_total;
    } result;
    void run(const std::filesystem::path &graph_file) override;
    void print_results(std::ostream &os) override
    {
        os << "Solution cost: " << result.solution_cost << "\n";
        os << "Solution size: " << result.solution_size << "\n";
        os << "Time reduction (s): " << result.time_reduction << "\n";
        os << "Time solving (s): " << result.time_solving << "\n";
        os << "Time total (s): " << result.time_total << "\n";
    }
};

// ==================== Direct Graph Algorithms ====================

class DirectGreedyRunner : public AlgorithmRunner
{
    struct Result
    {
        double solution_cost;
        size_t solution_size;
        double time_total;
    } result;

public:
    void run(const std::filesystem::path &graph_file) override;
    void print_results(std::ostream &os) override
    {
        os << "Solution cost: " << result.solution_cost << "\n";
        os << "Solution size: " << result.solution_size << "\n";
        os << "Time total (s): " << result.time_total << "\n";
    }
};

class GWCRunner : public AlgorithmRunner
{
    struct Result
    {
        double solution_cost;
        size_t solution_size;
        double time_total;
    } result;

public:
    void run(const std::filesystem::path &graph_file) override;
    void print_results(std::ostream &os) override
    {
        os << "Solution cost: " << result.solution_cost << "\n";
        os << "Solution size: " << result.solution_size << "\n";
        os << "Time total (s): " << result.time_total << "\n";
    }
};

class MSTConnectRunner : public AlgorithmRunner
{
    struct Result
    {
        double solution_cost;
        size_t solution_size;
        double time_total;
    } result;

public:
    void run(const std::filesystem::path &graph_file) override;
    void print_results(std::ostream &os) override
    {
        os << "Solution cost: " << result.solution_cost << "\n";
        os << "Solution size: " << result.solution_size << "\n";
        os << "Time total (s): " << result.time_total << "\n";
    }
};

class DirectILPRunner : public AlgorithmRunner
{
    struct Result
    {
        double solution_cost;
        size_t solution_size;
        double time_total;
    } result;

public:
    void run(const std::filesystem::path &graph_file) override;
    void print_results(std::ostream &os) override
    {
        os << "Solution cost: " << result.solution_cost << "\n";
        os << "Solution size: " << result.solution_size << "\n";
        os << "Time total (s): " << result.time_total << "\n";
    }
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
