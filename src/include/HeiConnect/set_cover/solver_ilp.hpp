#pragma once

#include <filesystem>
#include <memory>
#include <optional>
#include <string>
#include <string_view>
#include <tuple>
#include <utility>
#include <vector>

#include <gurobi_c++.h>

#include "HeiConnect/pipeline/common.hpp"
#include "HeiConnect/set_cover/solver_base.hpp"

struct SetCoverILPModel
{
    std::shared_ptr<GRBEnv> environment;
    std::shared_ptr<GRBModel> model;
    std::vector<GRBVar> variables;
};

template<typename SetCoverType>
SetCoverILPModel build_set_cover_ilp_model(const SetCoverType& set_cover)
{
    auto environment = std::make_shared<GRBEnv>(true);
    environment->set("LogFile", "set_cover_ilp.log");
    environment->start();

    auto model = std::make_shared<GRBModel>(*environment);
    std::vector<GRBVar> variables(set_cover.get_num_sets());
    for (size_t set_index = 0; set_index < set_cover.get_num_sets(); ++set_index)
    {
        variables[set_index] = model->addVar(
            0.0,
            1.0,
            set_cover.get_set_cost(set_index),
            GRB_BINARY,
            "s" + std::to_string(set_index));
    }

    std::vector<GRBLinExpr> cover_expr(set_cover.get_num_elements());
    for (size_t set_index = 0; set_index < set_cover.get_num_sets(); ++set_index)
    {
        set_cover.forEachElement(set_index, [&](size_t element) { cover_expr[element] += variables[set_index]; });
    }

    for (size_t element = 0; element < set_cover.get_num_elements(); ++element)
    {
        model->addConstr(cover_expr[element] >= 1, "cover_e" + std::to_string(element));
    }

    model->set(GRB_IntAttr_ModelSense, GRB_MINIMIZE);
    return {std::move(environment), std::move(model), std::move(variables)};
}

inline void write_set_cover_ilp_model(SetCoverILPModel& ilp, const std::filesystem::path& output_file)
{
    ilp.model->update();
    ilp.model->write(output_file.string());
}

template<typename SolutionType>
bool solve_set_cover_ilp_model(SetCoverILPModel& ilp, SolutionType& solution)
{
    ilp.model->optimize();

    const int status = ilp.model->get(GRB_IntAttr_Status);
    if (status != GRB_OPTIMAL && status != GRB_SUBOPTIMAL)
        return false;

    for (size_t set_index = 0; set_index < ilp.variables.size(); ++set_index)
    {
        if (ilp.variables[set_index].get(GRB_DoubleAttr_X) > 0.5)
            solution.add_set(set_index);
    }
    return true;
}

struct BuildSetCoverILPStage
{
    static constexpr std::string_view name = "Build set cover ILP";

    template<typename SetCoverType, typename Context>
    auto operator()(std::shared_ptr<const SetCoverType> set_cover, Context context) const
    {
        auto ilp = build_set_cover_ilp_model(*set_cover);
        return std::tuple{std::move(set_cover), std::move(context), std::move(ilp)};
    }
};

template<size_t RecordMetricsLevel = 0>
class SetCoverSolverILP
{
public:
    static constexpr std::string_view name = "ILP solve";

    static std::string grb_get_status_string(int status)
    {
        // TODO: add all
        switch (status)
        {
            case GRB_OPTIMAL: return "OPTIMAL";
            case GRB_INFEASIBLE: return "INFEASIBLE";
            case GRB_UNBOUNDED: return "UNBOUNDED";
            case GRB_INF_OR_UNBD: return "INF_OR_UNBD";
            case GRB_TIME_LIMIT: return "TIME_LIMIT";
            case GRB_ITERATION_LIMIT: return "ITERATION_LIMIT";
            case GRB_NODE_LIMIT: return "NODE_LIMIT";
            case GRB_USER_OBJ_LIMIT: return "USER_OBJ_LIMIT";
            default: return "UNKNOWN";
        }
    }

    template<typename SetCoverT, typename Context>
    auto operator()(std::shared_ptr<const SetCoverT> set_cover, Context context)
    {
        solve(*set_cover, context);
        return std::tuple{std::move(set_cover), std::move(context)};
    }

    template<typename SetCoverT, typename Context>
    auto operator()(std::shared_ptr<const SetCoverT> set_cover, Context context, SetCoverILPModel ilp)
    {
        solve(ilp, context);
        return std::tuple{std::move(set_cover), std::move(context)};
    }

    template<typename SetCoverType, typename SolutionType>
    bool solve(const SetCoverType& set_cover, SolutionType& solution)
    {
        if (set_cover.get_num_elements() == 0)
        {
            m_feasible = true;
            if constexpr (RecordMetricsLevel > 0)
            {
                m_metrics = StageMetrics{
                    {"cost", "0"},
                    {"size", std::to_string(solution.get_solution_size())},
                    {"status", "OPTIMAL"},
                };
            }
            return true;
        }
        auto ilp = build_set_cover_ilp_model(set_cover);
        return solve(ilp, solution);
    }

    template<typename SolutionType>
    bool solve(SetCoverILPModel& ilp, SolutionType& solution)
    {
        m_feasible = solve_set_cover_ilp_model(ilp, solution);
        const int status = ilp.model->get(GRB_IntAttr_Status);
        if constexpr (RecordMetricsLevel > 0)
        {
            m_metrics = StageMetrics{
                {"cost", m_feasible ? std::to_string(ilp.model->get(GRB_DoubleAttr_ObjVal)) : ""},
                {"size", std::to_string(solution.get_solution_size())},
                {"status", grb_get_status_string(status)},
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
