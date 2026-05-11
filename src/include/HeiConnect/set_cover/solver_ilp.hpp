#pragma once

#include <gurobi_c++.h>

#include "HeiConnect/set_cover/solver_base.hpp"

template <typename SetCoverType, typename SolutionType>
    requires SetCoverCon<SetCoverType> && ForEachElementCon<SetCoverType>
class SetCoverSolverILP
{
public:
    bool solve(const SetCoverType &set_cover, SolutionType &solution)
    {
        GRBEnv env = GRBEnv(true);
        env.set("LogFile", "set_cover_ilp.log");
        env.start();
        GRBModel model = GRBModel(env);

        std::vector<GRBVar> vars(set_cover.get_num_sets());
        for (size_t set_index = 0; set_index < set_cover.get_num_sets(); ++set_index)
        {
            vars[set_index] = model.addVar(0.0, 1.0, set_cover.get_set_cost(set_index), GRB_BINARY, "s" + std::to_string(set_index));
        }

        std::vector<GRBLinExpr> cover_expr(set_cover.get_num_elements());
        for (size_t set_index = 0; set_index < set_cover.get_num_sets(); ++set_index)
        {
            set_cover.forEachElement(set_index, [&](size_t element)
                                     { cover_expr[element] += vars[set_index]; });
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
        return true;
    }

    bool is_feasible() const noexcept
    {
        return m_feasible;
    }

private:
    bool m_feasible = false;
};
