#pragma once
#include <unordered_set>
#include <vector>
#include <stdexcept>

#include "HeiConnect/minmax.hpp"

struct SetCover
{
    const std::vector<size_t> a;
    const std::vector<size_t> b;
    const std::vector<double> costs;

    SetCover(std::vector<size_t> a,
             std::vector<size_t> b,
             std::vector<double> costs)
        : a(std::move(a)), b(std::move(b)), costs(std::move(costs))
    {
        if (this->a.size() - 1 != this->costs.size())
        {
            throw std::invalid_argument("SetCover: size of a and costs must be equal. Every subset should have a cost.");
        }
    };
};

template <typename Derived>
class SetCoverSolver
{
public:
    SetCoverSolver(const SetCover &set_cover)
        : set_cover(set_cover)
    {
        // `b` stores element indices; number of elements is max index + 1
        NUM_ELEMENTS = max(set_cover.b) + 1;
        // initialize covered counters (one per element) to zero
        m_covered_elements = std::vector<int>(NUM_ELEMENTS, 0);
    }
    ~SetCoverSolver() = default;
    void solve()
    {
        static_cast<Derived *>(this)->solve();
    };
    std::unordered_set<size_t> get_solution() const noexcept;
    double get_solution_cost() const
    {
        double total_cost = 0.0;
        for (const auto &set_index : m_solution)
        {
            total_cost += set_cover.costs[set_index];
        }
        return total_cost;
    }

protected:
    void add_set(size_t set_index) noexcept;
    void remove_set(size_t set_index) noexcept;
    void trim_solution();

protected:
    size_t NUM_ELEMENTS;
    SetCover set_cover;
    size_t m_total_covered_elements = 0;
    std::unordered_set<size_t> m_solution{};
    std::vector<bool> m_chosen_sets;
    bool m_feasible = false;
    std::vector<int> m_covered_elements;
};

template <typename Derived>
std::unordered_set<size_t> SetCoverSolver<Derived>::get_solution() const noexcept
{
    return m_solution;
}

template <typename Derived>
void SetCoverSolver<Derived>::add_set(size_t set_index) noexcept
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

template <typename Derived>
void SetCoverSolver<Derived>::remove_set(size_t set_index) noexcept
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

template <typename Derived>
void SetCoverSolver<Derived>::trim_solution()
{
    return;
}
class SetCoverSolverGreedySingleThreadedPQ : public SetCoverSolver<SetCoverSolverGreedySingleThreadedPQ>
{
public:
    void solve();
};
class SetCoverSolverGreedyParallel : public SetCoverSolver<SetCoverSolverGreedyParallel>
{
public:
    void solve();
};

class SetCoverSolverSharpGreedy : public SetCoverSolver<SetCoverSolverSharpGreedy>
{
public:
    void solve();
};

class SetCoverSolverILP : public SetCoverSolver<SetCoverSolverILP>
{
private:
    std::vector<bool> m_chosen_sets;

public:
    void solve();
};

class SetCoverSolverGreedyWideSingleThreaded : public SetCoverSolver<SetCoverSolverGreedyWideSingleThreaded>
{
public:
    void solve();
};