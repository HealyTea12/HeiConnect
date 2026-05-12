#pragma once

#include <array>
#include <cstddef>
#include <span>

template<std::size_t MaxMoveSize>
class SmallMove
{
public:
    using SetID = std::size_t;

    SmallMove() = default;

    void push_back(SetID set)
    {
        m_sets[m_size++] = set;
    }

    SetID operator[](std::size_t i) const
    {
        return m_sets[i];
    }

    std::size_t size() const
    {
        return m_size;
    }

    auto begin() const
    {
        return m_sets.begin();
    }

    auto end() const
    {
        return m_sets.begin() + static_cast<std::ptrdiff_t>(m_size);
    }

private:
    std::array<SetID, MaxMoveSize> m_sets{};
    std::size_t m_size = 0;
};