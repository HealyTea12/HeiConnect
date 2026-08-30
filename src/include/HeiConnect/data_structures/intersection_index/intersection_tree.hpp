#pragma once

#include "HeiConnect/data_structures/intersection_index/baseline.hpp"

#include <algorithm>
#include <limits>
#include <numeric>
#include <type_traits>

template<int RecordStatsLevel = 0, bool PruneByLevel = false>
class IntersectionTreeIdx : public BaseIntersectionIdx<RecordStatsLevel>
{
public:
    using typename BaseIntersectionIdx<RecordStatsLevel>::Index;
    using typename BaseIntersectionIdx<RecordStatsLevel>::Interval;

    IntersectionTreeIdx() = default;

    explicit IntersectionTreeIdx(const std::vector<IntersectionRecord>& records) : m_records(records)
    {
        m_active.resize(records.size(), true);
        m_node_by_interval.resize(records.size());
        std::vector<Index> indices(records.size());
        std::iota(indices.begin(), indices.end(), 0);
        std::sort(indices.begin(), indices.end(), [&](Index left, Index right) {
            return m_records[left].interval < m_records[right].interval;
        });
        m_root = build(indices, 0, indices.size(), -1);
        m_initial_max_end.reserve(m_nodes.size());
        if constexpr (PruneByLevel)
        {
            m_initial_min_level.reserve(m_nodes.size());
        }
        m_initial_levels.reserve(records.size());
        for (const IntersectionRecord& record : records)
        {
            m_initial_levels.push_back(record.level);
        }
        for (const Node& node : m_nodes)
        {
            m_initial_max_end.push_back(node.max_end);
            if constexpr (PruneByLevel)
            {
                m_initial_min_level.push_back(node.min_level);
            }
        }
    }

    std::unique_ptr<BaseIntersectionIdx<RecordStatsLevel>>
    make(const std::vector<IntersectionRecord>& records) const override
    {
        return std::make_unique<IntersectionTreeIdx<RecordStatsLevel, PruneByLevel>>(records);
    }

    void forEachIntersection(std::function<void(Index, Interval)> callback, Interval query, size_t exclusive_level)
        const override
    {
        if constexpr (RecordStatsLevel > 0)
        {
            m_metrics.queries++;
        }
        queryTree(m_root, query, exclusive_level, callback);
    }

    void popInterval(Interval interval) override
    {
        for (Index index = 0; index < m_records.size(); ++index)
        {
            if (m_active[index] && m_records[index].interval == interval)
            {
                popInterval(index);
                return;
            }
        }
    }

    void popInterval(Index index) override
    {
        if (index < m_active.size() && m_active[index])
        {
            if constexpr (RecordStatsLevel > 0)
            {
                m_metrics.intervals_popped++;
            }
            m_active[index] = false;
            int node = m_node_by_interval[index];
            while (node >= 0)
            {
                update(node);
                node = m_nodes[node].parent;
            }
        }
    }

    void activateInterval(Index index, size_t level) override
    {
        if (index >= m_active.size())
        {
            return;
        }

        m_records[index].level = level;
        m_active[index] = true;
        int node = m_node_by_interval[index];
        while (node >= 0)
        {
            update(node);
            node = m_nodes[node].parent;
        }
    }

    void clear() override
    {
        std::fill(m_active.begin(), m_active.end(), false);
        for (size_t node = m_nodes.size(); node > 0; --node)
        {
            update(static_cast<int>(node - 1));
        }
    }

    void reset() override
    {
        std::fill(m_active.begin(), m_active.end(), true);
        for (Index index = 0; index < m_records.size(); ++index)
        {
            m_records[index].level = m_initial_levels[index];
        }
        for (size_t node = 0; node < m_nodes.size(); ++node)
        {
            m_nodes[node].max_end = m_initial_max_end[node];
            if constexpr (PruneByLevel)
            {
                m_nodes[node].min_level = m_initial_min_level[node];
            }
        }
    }

    IntersectionIndexMetrics emit_metrics() const override
    {
        return m_metrics;
    }

private:
    struct Empty
    {};

    struct Node
    {
        Index interval{};
        size_t max_end{};
        [[no_unique_address]] std::conditional_t<PruneByLevel, size_t, Empty> min_level;
        int left{-1};
        int right{-1};
        int parent{-1};
    };

    static constexpr size_t inactive_end = std::numeric_limits<size_t>::min();
    static constexpr size_t inactive_level = std::numeric_limits<size_t>::max();

    int build(const std::vector<Index>& indices, size_t begin, size_t end, int parent)
    {
        if (begin == end)
        {
            return -1;
        }

        const size_t middle = begin + (end - begin) / 2;
        const int node = static_cast<int>(m_nodes.size());
        const Index interval = indices[middle];
        m_nodes.push_back({});
        m_nodes[node].interval = interval;
        m_nodes[node].max_end = std::get<1>(m_records[interval].interval);
        if constexpr (PruneByLevel)
        {
            m_nodes[node].min_level = m_records[interval].level;
        }
        m_nodes[node].parent = parent;
        m_node_by_interval[indices[middle]] = node;
        m_nodes[node].left = build(indices, begin, middle, node);
        m_nodes[node].right = build(indices, middle + 1, end, node);
        update(node);
        return node;
    }

    size_t maxEnd(int node) const
    {
        return node < 0 ? inactive_end : m_nodes[node].max_end;
    }

    size_t minLevel(int node) const
    {
        if constexpr (PruneByLevel)
        {
            return node < 0 ? inactive_level : m_nodes[node].min_level;
        }
        return inactive_level;
    }

    void update(int node)
    {
        const auto index = m_nodes[node].interval;
        const size_t own_end = m_active[index] ? std::get<1>(m_records[index].interval) : inactive_end;
        m_nodes[node].max_end = std::max({own_end, maxEnd(m_nodes[node].left), maxEnd(m_nodes[node].right)});
        if constexpr (PruneByLevel)
        {
            const size_t own_level = m_active[index] ? m_records[index].level : inactive_level;
            m_nodes[node].min_level =
                std::min({own_level, minLevel(m_nodes[node].left), minLevel(m_nodes[node].right)});
        }
    }

    void
    queryTree(int node, Interval query, size_t exclusive_level, const std::function<void(Index, Interval)>& callback)
        const
    {
        if (node < 0 || m_nodes[node].max_end < std::get<0>(query))
        {
            return;
        }
        if constexpr (PruneByLevel)
        {
            if (m_nodes[node].min_level >= exclusive_level)
            {
                if constexpr (RecordStatsLevel > 1)
                {
                    m_metrics.subtrees_pruned_by_level++;
                }
                return;
            }
        }

        if constexpr (RecordStatsLevel > 1)
        {
            m_metrics.candidates_inspected++;
        }

        queryTree(m_nodes[node].left, query, exclusive_level, callback);
        const Index index = m_nodes[node].interval;
        const Interval interval = m_records[index].interval;
        if (m_active[index] && (!PruneByLevel || m_records[index].level < exclusive_level) &&
            std::get<0>(interval) <= std::get<1>(query) && intervals_touch_or_cross(interval, query))
        {
            if constexpr (RecordStatsLevel > 0)
            {
                m_metrics.callbacks++;
            }
            callback(index, interval);
        }
        if (std::get<0>(interval) <= std::get<1>(query))
        {
            queryTree(m_nodes[node].right, query, exclusive_level, callback);
        }
    }

    std::vector<IntersectionRecord> m_records;
    std::vector<Node> m_nodes;
    std::vector<char> m_active;
    std::vector<int> m_node_by_interval;
    std::vector<size_t> m_initial_max_end;
    std::vector<size_t> m_initial_min_level;
    std::vector<size_t> m_initial_levels;
    int m_root{-1};
    mutable IntersectionIndexMetrics m_metrics;
};

template<int RecordStatsLevel = 0>
using WeightedIntersectionTreeIdx = IntersectionTreeIdx<RecordStatsLevel, true>;
