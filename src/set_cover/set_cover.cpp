#include <vector>
#include <unordered_set>
#include <omp.h>
#include <bits/stdc++.h>
#include <queue>
#include <algorithm>

#include "minmax.hpp"
#include "set_cover/set_cover.hpp"

#include <gurobi_c++.h>

void SetCoverSolverGreedyParallel::solve()
{
    std::vector<double> cost_benefit_ratios = std::vector<double>(set_cover.a.size() - 1, 0.);
    // stop when all elements are covered or we've chosen all available sets
    while (m_total_covered_elements < NUM_ELEMENTS && m_solution.size() != set_cover.a.size() - 1)
    {
        // calculate ratios
#pragma omp parallel for
        for (size_t i = 0; i < set_cover.a.size() - 1; i++)
        {
            // skip already selected sets
            if (m_solution.find(i) != m_solution.end())
            {
                cost_benefit_ratios[i] = -std::numeric_limits<double>::infinity();
                continue;
            }
            int covered = 0;
            for (size_t j = set_cover.a[i]; j < set_cover.a[i + 1]; j++)
            {
                if (m_covered_elements[set_cover.b[j]] == 0)
                    covered += 1;
            }
            if (covered > 0)
                cost_benefit_ratios[i] = covered / set_cover.costs[i];
            else
                cost_benefit_ratios[i] = 0;
        }
        auto best_set = argmax(cost_benefit_ratios);
        this->add_set(best_set);
    }
};

// The idea of this method is to maintain a priority queue of sets based on their cost-benefit ration.
// Then, to select the next best set, we pop the queue. But we need to check if the ration has changed.
// The point is to recalculate the ratios as little as possible.
void SetCoverSolverGreedySingleThreadedPQ::solve()
{
    std::priority_queue<std::pair<double, size_t>> pq;
    // create a heap for the indices of all sets based on their cost-benefit ratio
    for (size_t i = 0; i < set_cover.a.size() - 1; i++)
    {
        int covered = 0;
        for (size_t j = set_cover.a[i]; j < set_cover.a[i + 1]; j++)
        {
            if (m_covered_elements[set_cover.b[j]] == 0)
                covered += 1;
        }
        double ratio = (covered > 0) ? (covered / set_cover.costs[i]) : 0.;
        pq.push({ratio, i});
    }
    // compare based on the ratios
    auto compare = [](std::pair<double, size_t> a, std::pair<double, size_t> b)
    {
        return a.first < b.first;
    };
    // stop when all elements are covered or we've chosen all available sets while (m_total_covered_elements < NUM_ELEMENTS && m_solution.size() != set_cover.a.size() - 1)
    while (m_total_covered_elements < NUM_ELEMENTS && m_solution.size() != set_cover.a.size() - 1)
    {
        auto best_set = pq.top();
        pq.pop();
        // check if ratio has changed
        size_t best_set_idx = best_set.second;
        size_t covered = 0;
        for (auto i = set_cover.a[best_set_idx]; i < set_cover.a[best_set_idx + 1]; i++)
        {
            if (m_covered_elements[set_cover.b[i]] == 0)
                covered += 1;
        }
        double ratio = covered / set_cover.costs[best_set_idx];
        // if changed that put back in and try again
        if (ratio < best_set.first)
        {
            pq.push({ratio, best_set_idx});
            continue;
        }
        // otherwise add the set to the solution
        this->add_set(best_set_idx);
    }
};

// assumes that a represents the elements, and b the sets covering them
void SetCoverSolverSharpGreedy::solve()
{
    while (m_total_covered_elements < NUM_ELEMENTS && m_solution.size() != NUM_ELEMENTS)
    {
        for (size_t e = 0; e < set_cover.a.size() - 1; e++)
        {
            if (m_covered_elements[e] > 0)
                continue;
            auto cheapest_set_cost = std::numeric_limits<double>::max();
            size_t cheapest_set = std::numeric_limits<size_t>::max();
            for (auto i = set_cover.a[e]; i < set_cover.a[e + 1]; i++)
            {
                auto set = set_cover.b[i];
                if (set_cover.costs[set] < cheapest_set_cost)
                {
                    cheapest_set_cost = set_cover.costs[set];
                    cheapest_set = set;
                }
            }
            if (cheapest_set == std::numeric_limits<size_t>::max())
            {
                throw std::runtime_error(
                    "SetCoverSolverSharpGreedy: No set covers element " + std::to_string(e) + ". \n" +
                    "Problem instance is unsolvable.");
            }
            this->add_set(cheapest_set);
        }
    }
}

void SetCoverSolverILP::solve()
{
    GRBEnv env = GRBEnv(true);
    env.set("LogFile", "set_cover_ilp.log");
    env.start();
    GRBModel model = GRBModel(env);
    // create variables for each set
    std::vector<GRBVar> vars = std::vector<GRBVar>(set_cover.a.size() - 1);
    for (size_t i = 0; i < set_cover.a.size() - 1; i++)
    {
        vars[i] = model.addVar(0.0, 1.0, set_cover.costs[i], GRB_BINARY, "s" + std::to_string(i));
    }
    // create coverage constraints: for each element e, sum_{sets i covering e} vars[i] >= 1
    // build linear expressions for each element by iterating sets
    std::vector<GRBLinExpr> cover_expr = std::vector<GRBLinExpr>(NUM_ELEMENTS);
    for (size_t i = 0; i < set_cover.a.size() - 1; i++)
    {
        for (size_t j = set_cover.a[i]; j < set_cover.a[i + 1]; j++)
        {
            size_t elem = set_cover.b[j];
            cover_expr[elem] += vars[i];
        }
    }
    for (size_t e = 0; e < NUM_ELEMENTS; e++)
    {
        model.addConstr(cover_expr[e] >= 1, "cover_e" + std::to_string(e));
    }
    // objective
    GRBLinExpr obj = 0;
    for (size_t i = 0; i < set_cover.a.size() - 1; i++)
    {
        obj += set_cover.costs[i] * vars[i];
    }
    model.setObjective(obj, GRB_MINIMIZE);
    model.optimize();
    // check optimization status
    int status = model.get(GRB_IntAttr_Status);
    std::vector<bool> chosen_sets = std::vector<bool>(set_cover.a.size() - 1, false);
    if (status == GRB_OPTIMAL || status == GRB_SUBOPTIMAL)
    {
        for (size_t i = 0; i < set_cover.a.size() - 1; i++)
        {
            double val = vars[i].get(GRB_DoubleAttr_X);
            if (val > 0.5)
                chosen_sets[i] = true;
        }
        m_feasible = true;
        m_chosen_sets = chosen_sets;
    }
    else
    {
        m_feasible = false;
    }
}
