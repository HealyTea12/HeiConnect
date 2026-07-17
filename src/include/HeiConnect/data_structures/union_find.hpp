#pragma once

#include <vector>
#include <numeric>


class UnionFind
{
public:
    explicit UnionFind(size_t n) : m_parent(n)
    {
        std::iota(m_parent.begin(), m_parent.end(), 0);
    }

    size_t find(size_t x)
    {
        while (x != m_parent[x])
        {
            x = m_parent[x];
        }
        return x;
    }

    void unite(size_t a, size_t b)
    {
        a = find(a);
        b = find(b);

        if (a != b)
        {
            m_parent[b] = a;
        }
    }

private:
    std::vector<size_t> m_parent;
};