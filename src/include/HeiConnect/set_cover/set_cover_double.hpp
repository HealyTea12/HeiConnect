#pragma once

#include "HeiConnect/set_cover/common.hpp"

#include <cmath>
#include <stdexcept>
#include <utility>
#include <concepts>

template <typename SetCoverT1, typename SetCoverT2>
    requires SetCoverCon<SetCoverT1> && SetCoverCon<SetCoverT2> &&
             ForEachElementCon<SetCoverT1> && ForEachElementCon<SetCoverT2>
class SetCoverDouble
{
public:
    using ElementID = size_t;
    using SetID = size_t;
    using SetCost = double;

    SetCoverDouble(SetCoverT1 first, SetCoverT2 second)
        : m_first(std::move(first)),
          m_second(std::move(second)),
          m_numSets(m_first.get_num_sets()),
          m_firstElements(m_first.get_num_elements()),
          m_secondElements(m_second.get_num_elements())
    {
        if (m_first.get_num_sets() != m_second.get_num_sets())
        {
            throw std::invalid_argument("SetCoverDouble: both representations must have the same number of sets.");
        }

        // Ensure that both views describe the same set costs.
        constexpr double eps = 1e-12;
        for (SetID s = 0; s < m_numSets; ++s)
        {
            const double c1 = static_cast<double>(m_first.get_set_cost(s));
            const double c2 = static_cast<double>(m_second.get_set_cost(s));
            if (std::abs(c1 - c2) > eps)
            {
                throw std::invalid_argument("SetCoverDouble: set costs differ between representations.");
            }
        }
    }

    size_t get_num_sets() const noexcept
    {
        return m_numSets;
    }

    size_t get_num_elements() const noexcept
    {
        return m_firstElements + m_secondElements;
    }

    SetCost get_set_cost(SetID set_index) const
    {
        return m_first.get_set_cost(set_index);
    }

    void forEachElement(SetID set_index, const std::function<void(ElementID)> &func) const
    {
        m_first.forEachElement(set_index, func);
        m_second.forEachElement(set_index, [&](ElementID element)
                                { func(m_firstElements + element); });
    }

    const SetCoverT1 &first() const noexcept
    {
        return m_first;
    }

    const SetCoverT2 &second() const noexcept
    {
        return m_second;
    }

    size_t first_num_elements() const noexcept
    {
        return m_firstElements;
    }

private:
    SetCoverT1 m_first;
    SetCoverT2 m_second;
    size_t m_numSets;
    size_t m_firstElements;
    size_t m_secondElements;
};
