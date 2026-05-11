#pragma once

#include "HeiConnect/set_cover/common.hpp"

#include <ranges>

struct SetCover
{
public:
    using SetID = size_t;
    using ElementID = size_t;
    using SetCost = double;
    SetCover(std::vector<size_t> a,
             std::vector<ElementID> b,
             std::vector<SetCost> costs,
             size_t num_elements)
        : m_a(std::move(a)), m_b(std::move(b)), m_costs(std::move(costs)), m_num_elements(num_elements)
    {
        if (this->m_a.size() - 1 != this->m_costs.size())
        {
            throw std::invalid_argument("SetCover: size of a and costs must be equal. Every subset should have a cost.");
        }
    };
    size_t get_num_sets() const noexcept
    {
        return m_a.size() - 1;
    }
    size_t get_num_elements() const noexcept
    {
        return m_num_elements;
    }
    SetCost get_set_cost(size_t set_index) const
    {
        return m_costs[set_index];
    }
    auto set_elements(size_t set_index) const
    {
        return std::ranges::subrange(m_b.cbegin() + m_a[set_index], m_b.cbegin() + m_a[set_index + 1]);
    }
    std::vector<ElementID>::const_iterator set_begin(size_t set_index) const
    {
        return set_elements(set_index).begin();
    };
    std::vector<ElementID>::const_iterator set_end(size_t set_index) const
    {
        return set_elements(set_index).end();
    }
    void forEachElement(size_t set_index, const std::function<void(ElementID)> &func) const
    {
        for (ElementID idx = m_a[set_index]; idx < m_a[set_index + 1]; ++idx)
        {
            func(m_b[idx]);
        }
    }
    const std::vector<size_t> &get_a() const noexcept
    {
        return m_a;
    }
    const std::vector<ElementID> &get_b() const noexcept
    {
        return m_b;
    }
    const std::vector<SetCost> &get_costs() const noexcept
    {
        return m_costs;
    }

private:
    const std::vector<size_t> m_a;
    const std::vector<ElementID> m_b;
    const std::vector<SetCost> m_costs;
    const size_t m_num_elements;
};

static_assert(SetCoverCon<SetCover>);
