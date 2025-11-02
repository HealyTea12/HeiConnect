#include <vector>
#include <unordered_set>
#include <omp.h>
#include <bits/stdc++.h>

#include "set_cover/utils.hpp"
#include "set_cover/set_cover.hpp"

std::unordered_set<size_t> SetCoverSolver::get_solution() const noexcept
{
    return m_solution;
}

void SetCoverSolver::add_set(size_t set_index) noexcept
{
    for (size_t j{set_cover.a[set_index]}; j < set_cover.a[set_index + 1]; j++)
    {
        if (m_covered_elements[set_cover.b[j]] == 0)
        {
            m_total_covered_elements += 1;
        }
        m_covered_elements[set_cover.b[j]] += 1;
    }
    m_solution.insert(set_index);
}

void SetCoverSolver::remove_set(size_t set_index) noexcept
{
    for (size_t j{set_cover.a[set_index]}; j < set_cover.a[set_index + 1]; j++)
    {
        m_covered_elements[set_cover.b[j]] -= 1;
        if (m_covered_elements[set_cover.b[j]] == 0)
        {
            m_total_covered_elements -= 1;
        }
    }
    m_solution.erase(set_index);
}

void SetCoverSolverGreedyParallel::solve()
{
    std::vector<double> cost_benefit_ratios(set_cover.a.size() - 1, 0.);
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
