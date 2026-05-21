#pragma once

#include <algorithm>
#include <numeric>
#include <random>
#include <vector>

class CostFractionPerturbator
{
public:
    explicit CostFractionPerturbator(double p, std::mt19937_64& rng) : m_rng(rng), m_p(p)
    {}

    template<typename SetCoverType, typename SolutionType>
    void generateMove(const SetCoverType& set_cover, const SolutionType& solution, std::vector<size_t>& move)
    {
        using CostType = typename SetCoverType::SetCost;
        move.clear();
        if (solution.get_solution_size() == 0 || m_p <= 0.0)
        {
            return;
        }
        std::vector<size_t> solution_vec(solution.get_solution().begin(), solution.get_solution().end());
        const CostType solution_cost =
            std::accumulate(solution_vec.begin(), solution_vec.end(), CostType{}, [&](CostType acc, size_t set) {
                return acc + set_cover.get_set_cost(set);
            });
        std::shuffle(solution_vec.begin(), solution_vec.end(), m_rng);
        CostType removed_cost{};
        for (size_t set : solution_vec)
        {
            if (removed_cost / static_cast<double>(solution_cost) >= m_p)
            {
                break;
            }
            move.emplace_back(set);
            removed_cost += set_cover.get_set_cost(set);
        }
    }

    void set_p(double p)
    {
        m_p = p;
    }

    double get_p() const
    {
        return m_p;
    }

private:
    double m_p;
    std::mt19937_64& m_rng;
};