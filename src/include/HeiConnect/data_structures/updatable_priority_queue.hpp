#pragma once
#include <vector>

template <typename T, typename PT>
class UpdatablePriorityQueue
{
public:
    void push(const T &value, const PT &priority)
    {
        m_heap.push_back(std::make_pair(priority, value));
        sift_up(m_heap.size() - 1);
    }

    void pop()
    {
        std::swap(m_heap.front(), m_heap.back());
        m_heap.pop_back();
        sift_down(0);
    }

    const std::pair<PT, T> &top() const
    {
        return m_heap.front();
    }

    bool empty() const
    {
        return m_heap.empty();
    }

private:
    std::vector<std::pair<PT, T>> m_heap;
};