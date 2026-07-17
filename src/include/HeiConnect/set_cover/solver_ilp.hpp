#pragma once

#include <gurobi_c++.h>

#include "HeiConnect/set_cover/solver_base.hpp"

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

    template<typename SetCoverType, typename SolutionType>
    bool solve(const SetCoverType& set_cover, SolutionType& solution)
    {
        GRBEnv env = GRBEnv(true);
        env.set("LogFile", "set_cover_ilp.log");
        env.start();
        GRBModel model = GRBModel(env);

        std::vector<GRBVar> vars(set_cover.get_num_sets());
        for (size_t set_index = 0; set_index < set_cover.get_num_sets(); ++set_index)
        {
            vars[set_index] =
                model.addVar(0.0, 1.0, set_cover.get_set_cost(set_index), GRB_BINARY, "s" + std::to_string(set_index));
        }

        std::vector<GRBLinExpr> cover_expr(set_cover.get_num_elements());
        for (size_t set_index = 0; set_index < set_cover.get_num_sets(); ++set_index)
        {
            set_cover.forEachElement(set_index, [&](size_t element) { cover_expr[element] += vars[set_index]; });
        }

        for (size_t element = 0; element < set_cover.get_num_elements(); ++element)
        {
            model.addConstr(cover_expr[element] >= 1, "cover_e" + std::to_string(element));
        }

        GRBLinExpr objective = 0;
        for (size_t set_index = 0; set_index < set_cover.get_num_sets(); ++set_index)
        {
            objective += set_cover.get_set_cost(set_index) * vars[set_index];
        }
        model.setObjective(objective, GRB_MINIMIZE);
        model.optimize();

        const int status = model.get(GRB_IntAttr_Status);
        m_feasible = (status == GRB_OPTIMAL || status == GRB_SUBOPTIMAL);
        if (!m_feasible)
        {
            return false;
        }

        for (size_t set_index = 0; set_index < set_cover.get_num_sets(); ++set_index)
        {
            if (vars[set_index].get(GRB_DoubleAttr_X) > 0.5)
            {
                solution.add_set(set_index);
            }
        }
        if constexpr (RecordMetricsLevel > 0)
        {
            m_metrics = StageMetrics{
                {"cost", std::to_string(model.get(GRB_DoubleAttr_ObjVal))},
                {"size", std::to_string(solution.get_solution_size())},
                {"status", grb_get_status_string(status)},
            };
        }
        return true;
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
