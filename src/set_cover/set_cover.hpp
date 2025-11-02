#include <unordered_set>
#include <vector>
#include <stdexcept>

#include "set_cover/utils.hpp"

struct SetCover
{
    const std::vector<size_t> &a;
    const std::vector<size_t> &b;
    const std::vector<double> &costs;
    SetCover(const std::vector<size_t> &a,
             const std::vector<size_t> &b,
             const std::vector<double> &costs)
        : a(a), b(b), costs(costs)
    {
        if (a.size() - 1 != costs.size())
        {
            throw std::invalid_argument("SetCover: size of a and costs must be equal. Every subset should have a cost.");
        }
    };
};

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
    virtual void solve() = 0;
    std::unordered_set<size_t> get_solution() const noexcept;

protected:
    void add_set(size_t set_index) noexcept;
    void remove_set(size_t set_index) noexcept;

protected:
    size_t NUM_ELEMENTS;
    SetCover set_cover;
    size_t m_total_covered_elements = 0;
    std::unordered_set<size_t> m_solution{};
    std::vector<int> m_covered_elements;
};

class SetCoverSolverGreedyParallel : public SetCoverSolver
{
public:
    SetCoverSolverGreedyParallel(const SetCover &set_cover)
        : SetCoverSolver(set_cover) {};
    void solve() override;
};

// maybe should be in utils