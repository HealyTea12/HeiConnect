#pragma once

#include "HeiConnect/set_cover/common.hpp"

#include <ranges>

template<typename SetIDT = size_t, typename ElementIDT = size_t, typename SetCostT = double>
struct SetCover
{
public:
    using SetID = SetIDT;
    using ElementID = ElementIDT;
    using SetCost = SetCostT;
    SetCover(std::vector<size_t> a, std::vector<ElementID> b, std::vector<SetCost> costs, size_t num_elements) :
        m_a(std::move(a)),
        m_b(std::move(b)),
        m_costs(std::move(costs)),
        m_numElements(num_elements)
    {
        if (this->m_a.size() - 1 != this->m_costs.size())
        {
            throw std::invalid_argument(
                "SetCover: size of a and costs must be equal. Every subset should have a cost.");
        }
    };
    size_t get_num_sets() const noexcept
    {
        return m_a.size() - 1;
    }
    size_t get_num_elements() const noexcept
    {
        return m_numElements;
    }
    SetCost get_set_cost(SetID set_index) const
    {
        return m_costs[set_index];
    }
    auto set_elements(SetID set_index) const
    {
        return std::ranges::subrange(m_b.cbegin() + m_a[set_index], m_b.cbegin() + m_a[set_index + 1]);
    }
    std::vector<ElementID>::const_iterator set_begin(SetID set_index) const
    {
        return set_elements(set_index).begin();
    };
    std::vector<ElementID>::const_iterator set_end(SetID set_index) const
    {
        return set_elements(set_index).end();
    }
    void forEachElement(SetID set_index, const std::function<void(ElementID)>& func) const
    {
        for (ElementID idx = m_a[set_index]; idx < m_a[set_index + 1]; ++idx)
        {
            func(m_b[idx]);
        }
    }
    const std::vector<size_t>& get_a() const noexcept
    {
        return m_a;
    }
    const std::vector<ElementID>& get_b() const noexcept
    {
        return m_b;
    }
    const std::vector<SetCost>& get_costs() const noexcept
    {
        return m_costs;
    }
    template<typename NewSetCost>
    SetCover<SetID, ElementID, NewSetCost> discretize_costs(size_t num_bins) const
    {
        auto dc = discretize_weights<NewSetCost>(m_costs, num_bins);
        std::vector<NewSetCost> new_costs(dc.begin(), dc.end());
        return SetCover<SetID, ElementID, NewSetCost>(m_a, m_b, new_costs, m_numElements);
    }


private:
    const std::vector<size_t> m_a;
    const std::vector<ElementID> m_b;
    const std::vector<SetCost> m_costs;
    const size_t m_numElements;
};

static_assert(SetCoverCon<SetCover<>>);
