#pragma once

#include <exception>
#include <memory>
#include <optional>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include <gurobi_c++.h>

#include "HeiConnect/conn_aug/reducers/common.hpp"
#include "HeiConnect/pipeline/common.hpp"
#include "HeiConnect/sc_reduction/transform_single_core.hpp"
#include "HeiConnect/sc_reduction/transform_single_csr.hpp"
#include "HeiConnect/sc_reduction/transform_single_oracle_ancestry.hpp"
#include "HeiConnect/sc_reduction/transform_single_partition_matrix.hpp"
#include "HeiConnect/set_cover/set_cover_csr.hpp"
#include "HeiConnect/set_cover/solver_ilp.hpp"

namespace HeiConnect::detail
{

template<typename NodeID, typename LinkEdgeID, typename LinkEdgeWeight>
WeightedCRFGraph<NodeID, LinkEdgeID, LinkEdgeWeight> extend_link_graph_vertices(
    const WeightedCRFGraph<NodeID, LinkEdgeID, LinkEdgeWeight>& link_graph,
    size_t num_vertices)
{
    std::vector<NodeID> vertices = link_graph.graph.vertices;
    vertices.resize(num_vertices + 1, static_cast<NodeID>(link_graph.num_edges()));
    return {{std::move(vertices), link_graph.graph.edges}, link_graph.weights};
}

template<typename GraphType, typename LinkGraphType>
auto construct_block_tree_set_cover(const GraphType& graph, const LinkGraphType& link_graph)
{
    auto [min_cuts, num_min_cuts] = HeiConnect_details::generate_min_cut_matrix(
        graph.graph.vertices,
        graph.graph.edges,
        graph.weights,
        link_graph.graph.vertices,
        link_graph.graph.edges,
        link_graph.weights);
    return HeiConnect_details::construct_set_cover_csr_ull(
        graph.graph.vertices,
        graph.graph.edges,
        graph.weights,
        link_graph.graph.vertices,
        link_graph.graph.edges,
        link_graph.weights,
        min_cuts,
        num_min_cuts);
}

template<typename GraphType, typename LinkGraphType>
class LazyBlockTreeCallback : public GRBCallback
{
public:
    LazyBlockTreeCallback(
        const GraphType& graph,
        const LinkGraphType& link_graph,
        const std::vector<GRBVar>& variables) :
        m_graph(graph),
        m_link_graph(link_graph),
        m_variables(variables)
    {}

    void rethrow_if_failed() const
    {
        if (m_exception)
        {
            std::rethrow_exception(m_exception);
        }
    }

    size_t num_lazy_constraints() const noexcept
    {
        return m_num_lazy_constraints;
    }

protected:
    void callback() override
    {
        if (where != GRB_CB_MIPSOL)
        {
            return;
        }

        try
        {
            separate_solution();
        }
        catch (...)
        {
            m_exception = std::current_exception();
            abort();
        }
    }

private:
    void separate_solution()
    {
        std::vector<size_t> selected_links;
        for (size_t link_id = 0; link_id < m_variables.size(); ++link_id)
        {
            if (getSolution(m_variables[link_id]) > 0.5)
            {
                selected_links.push_back(link_id);
            }
        }

        UnionFind uf(m_graph.num_vertices());
        add_links_to_union_find_frozen_stack(m_graph, m_link_graph, uf, selected_links);

        auto [contracted_graph, contracted_link_graph, original_link_ids] =
            materialize_contractions_preserving_links(m_graph, m_link_graph, uf);
        if (contracted_graph.num_vertices() <= 1)
        {
            return;
        }

        auto [block_tree, cycle_positions] = contracted_graph.cactus_generate_block_tree(0);
        (void)cycle_positions;
        auto block_tree_link_graph = extend_link_graph_vertices(contracted_link_graph, block_tree.num_vertices());
        auto set_cover = construct_block_tree_set_cover(block_tree, block_tree_link_graph);

        std::vector<GRBLinExpr> expressions(set_cover.get_num_elements());
        for (size_t set_id = 0; set_id < set_cover.get_num_sets(); ++set_id)
        {
            const GRBVar& variable = m_variables[original_link_ids[set_id]];
            set_cover.forEachElement(set_id, [&](size_t element_id) { expressions[element_id] += variable; });
        }

        for (auto& expression : expressions)
        {
            addLazy(expression >= 1);
            ++m_num_lazy_constraints;
        }
    }

    const GraphType& m_graph;
    const LinkGraphType& m_link_graph;
    const std::vector<GRBVar>& m_variables;
    std::exception_ptr m_exception;
    size_t m_num_lazy_constraints = 0;
};

}

namespace HeiConnect
{

template<size_t RecordMetricsLevel = 0>
class LazyBlockTreeSolverILP
{
public:
    static constexpr std::string_view name = "Lazy block-tree ILP solve";

    template<typename GraphType, typename LinkGraphType, typename SolutionType>
    bool solve(const GraphType& graph, const LinkGraphType& link_graph, SolutionType& solution)
    {
        if (graph.num_vertices() <= 1)
        {
            m_feasible = true;
            if constexpr (RecordMetricsLevel > 0)
            {
                m_metrics = StageMetrics{
                    {"cost", "0"},
                    {"size", std::to_string(solution.get_solution_size())},
                    {"status", "OPTIMAL"},
                    {"lazy_constraints", "0"},
                };
            }
            return true;
        }

        auto [block_tree, cycle_positions] = graph.cactus_generate_block_tree(0);
        (void)cycle_positions;
        auto block_tree_link_graph = detail::extend_link_graph_vertices(link_graph, block_tree.num_vertices());
        auto set_cover = detail::construct_block_tree_set_cover(block_tree, block_tree_link_graph);
        auto ilp = build_set_cover_ilp_model(set_cover);

        detail::LazyBlockTreeCallback callback(graph, link_graph, ilp.variables);
        ilp.model->set(GRB_IntParam_LazyConstraints, 1);
        ilp.model->setCallback(&callback);
        m_feasible = solve_set_cover_ilp_model(ilp, solution);
        callback.rethrow_if_failed();

        if constexpr (RecordMetricsLevel > 0)
        {
            const int status = ilp.model->get(GRB_IntAttr_Status);
            m_metrics = StageMetrics{
                {"cost", m_feasible ? std::to_string(ilp.model->get(GRB_DoubleAttr_ObjVal)) : ""},
                {"size", std::to_string(solution.get_solution_size())},
                {"status", SetCoverSolverILP<>::grb_get_status_string(status)},
                {"lazy_constraints", std::to_string(callback.num_lazy_constraints())},
            };
        }
        return m_feasible;
    }

    std::optional<StageMetrics> emit_metrics() const
    {
        if constexpr (RecordMetricsLevel > 0)
        {
            return m_metrics;
        }
        else
        {
            return std::nullopt;
        }
    }

    bool is_feasible() const noexcept
    {
        return m_feasible;
    }

private:
    bool m_feasible = false;
    std::optional<StageMetrics> m_metrics;
};

}
