#pragma once

#include "HeiConnect/data_structures/intersection_index/baseline.hpp"

#include <stdexcept>
#include <algorithm>
#include <tuple>
#include <cstdint>
#include <functional>
#include <limits>
#include <numeric>
#include <queue>
#include <utility>
#include <vector>


/**
 * Cycle Domination Baseline.
 * Inputs:
 *  - links: a vector of tuples of [u,v,w] representing unique links on a cycle graph.
 *           They are normalized so that 0 <= u < v < n_nodes.
 * - n_nodes: the number of nodes in the cycle graph.
 *
 * This algorithm will performa a "bounded" Dijkstra from each node in the cycle graph.
 * The goal of the Dijkstra is to calculate the shortest distance to each other node in the cycle, following the link
 * graph. The link intersection graph is a graph where conceptually each node represents a link where you pay its weight
 * to enter it. There are also nodes representing the cycle nodes, which have 0 cost.
 *
 * Hence, this algorithm performs a Dijkstra on the link intersection graph starting at each cycle node.
 * It calculates the maximum weight of a link incident to the source node and uses it as a cutoff.
 * Once it relaxes a link, it may pop all neighbouring links off the intersectionIndex, and their final distance is
 * calculated. Once it relaxes a link, it has finalized the distance to that link's endpoints, so it can mark them as
 * complete and their distance is calculated. Once all endpoints are complete, the Dijkstra can terminate early.
 *
 * After each run of Dijkstra, we check if any of the links incident to the source node have a shorter distance to their
 * other endpoint. if so, we mark that link as removable.
 */
struct CycleReductionMetrics
{
    size_t sources{};
    size_t priority_queue_pops{};
    size_t links_enqueued{};
    size_t intersection_candidates_enqueued{};
    size_t intersection_candidates_already_explored{};
    size_t intersection_candidates_rejected_by_cutoff{};
    size_t termination_by_cutoff{};
    size_t termination_by_completion{};
    size_t termination_by_empty_queue{};
    size_t maximum_queue_size{};
    size_t maximum_pops_per_link{};
    size_t completed_vertices_at_stop{};
    size_t possible_priority_queue_pops{};
    IntersectionIndexMetrics intersection_index;

    void add(const CycleReductionMetrics& other)
    {
        sources += other.sources;
        priority_queue_pops += other.priority_queue_pops;
        links_enqueued += other.links_enqueued;
        intersection_candidates_enqueued += other.intersection_candidates_enqueued;
        intersection_candidates_already_explored += other.intersection_candidates_already_explored;
        intersection_candidates_rejected_by_cutoff += other.intersection_candidates_rejected_by_cutoff;
        termination_by_cutoff += other.termination_by_cutoff;
        termination_by_completion += other.termination_by_completion;
        termination_by_empty_queue += other.termination_by_empty_queue;
        maximum_queue_size = std::max(maximum_queue_size, other.maximum_queue_size);
        maximum_pops_per_link = std::max(maximum_pops_per_link, other.maximum_pops_per_link);
        completed_vertices_at_stop += other.completed_vertices_at_stop;
        possible_priority_queue_pops += other.possible_priority_queue_pops;
        intersection_index.add(other.intersection_index);
    }
};

template<int RecordStatsLevel = 0, typename NodeID, typename Weight>
auto cycle_domination_baseline(
    const std::vector<std::tuple<NodeID, NodeID, Weight>>& links,
    uint64_t n_nodes,
    const BaseIntersectionIdx<RecordStatsLevel>& intersection_index_type,
    CycleReductionMetrics* output_metrics = nullptr,
    bool reuse_intersection_index = true)
{
    using LinkID = size_t;
    using QueueEntry = std::pair<Weight, LinkID>;

    const LinkID invalid_link = std::numeric_limits<LinkID>::max();

    std::vector<std::vector<LinkID>> incident(n_nodes);
    std::vector<std::tuple<size_t, size_t>> intervals;
    std::vector<Weight> sorted_weights;
    intervals.reserve(links.size());
    sorted_weights.reserve(links.size());

    for (LinkID id = 0; id < links.size(); ++id)
    {
        const auto& [u, v, weight] = links[id];
        incident[static_cast<size_t>(u)].push_back(id);
        incident[static_cast<size_t>(v)].push_back(id);
        intervals.emplace_back(static_cast<size_t>(u), static_cast<size_t>(v));
        sorted_weights.push_back(weight);
    }
    std::sort(sorted_weights.begin(), sorted_weights.end());

    std::vector<IntersectionRecord> intersection_records;
    intersection_records.reserve(links.size());
    for (LinkID id = 0; id < links.size(); ++id)
    {
        const Weight weight = std::get<2>(links[id]);
        const size_t level = static_cast<size_t>(
            std::lower_bound(sorted_weights.begin(), sorted_weights.end(), weight) - sorted_weights.begin());
        intersection_records.push_back({intervals[id], level});
    }

    std::vector<char> removable(links.size(), false);
    std::vector<char> explored(links.size());
    std::vector<char> complete(n_nodes);
    std::vector<Weight> vertex_distance(n_nodes);
    std::vector<Weight> value_to_beat(n_nodes);
    std::vector<LinkID> candidate_link(n_nodes);
    auto intersection_index = intersection_index_type.make(intersection_records);
    CycleReductionMetrics metrics;
    std::vector<size_t> pops_per_link;
    if constexpr (RecordStatsLevel > 1)
    {
        pops_per_link.resize(links.size());
    }

    for (size_t source = 0; source < n_nodes; ++source)
    {
        if (incident[source].empty())
        {
            continue;
        }

        if constexpr (RecordStatsLevel > 0)
        {
            metrics.sources++;
            metrics.possible_priority_queue_pops += links.size();
        }
        if (reuse_intersection_index)
        {
            intersection_index->reset();
        }
        else
        {
            intersection_index = intersection_index_type.make(intersection_records);
        }

        std::fill(explored.begin(), explored.end(), false);
        std::fill(complete.begin(), complete.end(), false);
        std::fill(candidate_link.begin(), candidate_link.end(), invalid_link);
        size_t n_complete_vertices = 0;

        Weight cutoff{};
        for (LinkID id : incident[source])
        {
            const auto& [u, v, weight] = links[id];
            const size_t target = static_cast<size_t>(u) == source ? static_cast<size_t>(v) : static_cast<size_t>(u);
            value_to_beat[target] = weight;
            candidate_link[target] = id;
            cutoff = std::max(cutoff, weight);
        }

        std::priority_queue<QueueEntry, std::vector<QueueEntry>, std::greater<QueueEntry>> queue;

        auto complete_endpoint = [&](NodeID endpoint, Weight value) {
            const size_t vertex = static_cast<size_t>(endpoint);
            if (complete[vertex])
            {
                return;
            }

            n_complete_vertices++;
            complete[vertex] = true;
            vertex_distance[vertex] = value;

            if (candidate_link[vertex] != invalid_link && vertex_distance[vertex] < value_to_beat[vertex])
            {
                removable[candidate_link[vertex]] = true;
            }
        };

        for (LinkID id : incident[source])
        {
            explored[id] = true;
            intersection_index->popInterval(id);
            queue.emplace(std::get<2>(links[id]), id);
            if constexpr (RecordStatsLevel > 0)
            {
                metrics.links_enqueued++;
            }
        }

        if constexpr (RecordStatsLevel > 1)
        {
            metrics.maximum_queue_size = std::max(metrics.maximum_queue_size, queue.size());
        }

        bool stopped_early = false;
        while (!queue.empty())
        {
            const auto [current_distance, current] = queue.top();
            queue.pop();
            if constexpr (RecordStatsLevel > 0)
            {
                metrics.priority_queue_pops++;
            }
            if constexpr (RecordStatsLevel > 1)
            {
                pops_per_link[current]++;
            }

            const auto& [u, v, weight] = links[current];
            complete_endpoint(u, current_distance);
            complete_endpoint(v, current_distance);

            if (current_distance >= cutoff || n_complete_vertices == n_nodes)
            {
                if constexpr (RecordStatsLevel > 0)
                {
                    stopped_early = true;
                    if (n_complete_vertices == n_nodes)
                    {
                        metrics.termination_by_completion++;
                    }
                    else
                    {
                        metrics.termination_by_cutoff++;
                    }
                }
                break;
            }

            std::vector<LinkID> neighbours;
            std::vector<LinkID> rejected_by_cutoff;
            const Weight weight_limit = cutoff - current_distance;
            const size_t exclusive_level = static_cast<size_t>(
                std::lower_bound(sorted_weights.begin(), sorted_weights.end(), weight_limit) - sorted_weights.begin());
            intersection_index->forEachIntersection(
                [&](LinkID next, typename BaseIntersectionIdx<RecordStatsLevel>::Interval) {
                    const Weight next_weight = std::get<2>(links[next]);
                    if (explored[next])
                    {
                        if constexpr (RecordStatsLevel > 0)
                        {
                            metrics.intersection_candidates_already_explored++;
                        }
                    }
                    else if (next_weight >= cutoff - current_distance)
                    {
                        if constexpr (RecordStatsLevel > 0)
                        {
                            metrics.intersection_candidates_rejected_by_cutoff++;
                        }
                        rejected_by_cutoff.push_back(next);
                    }
                    else
                    {
                        explored[next] = true;
                        neighbours.push_back(next);
                        if constexpr (RecordStatsLevel > 0)
                        {
                            metrics.intersection_candidates_enqueued++;
                        }
                    }
                },
                intervals[current],
                exclusive_level);

            for (LinkID rejected : rejected_by_cutoff)
            {
                intersection_index->popInterval(rejected);
            }

            for (LinkID next : neighbours)
            {
                intersection_index->popInterval(next);
                queue.emplace(current_distance + std::get<2>(links[next]), next);
                if constexpr (RecordStatsLevel > 0)
                {
                    metrics.links_enqueued++;
                }
            }
            if constexpr (RecordStatsLevel > 1)
            {
                metrics.maximum_queue_size = std::max(metrics.maximum_queue_size, queue.size());
            }
        }
        if constexpr (RecordStatsLevel > 0)
        {
            if (!stopped_early)
            {
                metrics.termination_by_empty_queue++;
            }
        }
        if constexpr (RecordStatsLevel > 1)
        {
            metrics.completed_vertices_at_stop += n_complete_vertices;
        }
    }

    if constexpr (RecordStatsLevel > 1)
    {
        if (!pops_per_link.empty())
        {
            metrics.maximum_pops_per_link = *std::max_element(pops_per_link.begin(), pops_per_link.end());
        }
    }
    if constexpr (RecordStatsLevel > 0)
    {
        metrics.intersection_index = intersection_index->emit_metrics();
    }

    std::vector<int> result;
    for (LinkID id = 0; id < links.size(); ++id)
    {
        if (removable[id])
        {
            result.push_back(static_cast<int>(id));
        }
    }

    if constexpr (RecordStatsLevel > 0)
    {
        if (output_metrics != nullptr)
        {
            *output_metrics = metrics;
        }
    }
    return result;
}

template<int RecordStatsLevel = 0, typename NodeID, typename Weight>
auto cycle_domination_baseline(
    const std::vector<std::tuple<NodeID, NodeID, Weight>>& links,
    uint64_t n_nodes,
    CycleReductionMetrics* output_metrics = nullptr)
{
    const BaselineIntersectionIdx<RecordStatsLevel> intersection_index_type;
    return cycle_domination_baseline<RecordStatsLevel>(links, n_nodes, intersection_index_type, output_metrics);
}

/**
 * Single-pass cycle domination.
 *
 * The Dijkstra tree is rooted at the globally cheapest link. For every input link {u,v}, the cheapest tree path between
 * a link incident to u and a link incident to v is a valid connected replacement. The input link is removable when
 * that path is strictly cheaper than the input link.
 *
 * Unlike the baseline, this algorithm deliberately keeps exploring after every cycle vertex has been reached because
 * later links add paths to the tree used to certify removals. The search only stops at the global weight cutoff or when
 * the queue is empty.
 */
template<int RecordStatsLevel = 0, typename NodeID, typename Weight>
auto cycle_domination_single_pass(
    const std::vector<std::tuple<NodeID, NodeID, Weight>>& links,
    uint64_t n_nodes,
    const BaseIntersectionIdx<RecordStatsLevel>& intersection_index_type,
    CycleReductionMetrics* output_metrics = nullptr)
{
    using LinkID = size_t;
    using QueueEntry = std::pair<Weight, LinkID>;

    const LinkID invalid_link = std::numeric_limits<LinkID>::max();
    CycleReductionMetrics metrics;
    if (links.empty())
    {
        if constexpr (RecordStatsLevel > 0)
        {
            if (output_metrics != nullptr)
            {
                *output_metrics = metrics;
            }
        }
        return std::vector<int>{};
    }

    std::vector<std::tuple<size_t, size_t>> intervals;
    std::vector<Weight> sorted_weights;
    std::vector<std::vector<LinkID>> incident(n_nodes);
    intervals.reserve(links.size());
    sorted_weights.reserve(links.size());
    for (LinkID id = 0; id < links.size(); ++id)
    {
        const auto& [u, v, weight] = links[id];
        intervals.emplace_back(static_cast<size_t>(u), static_cast<size_t>(v));
        sorted_weights.push_back(weight);
        incident[static_cast<size_t>(u)].push_back(id);
        incident[static_cast<size_t>(v)].push_back(id);
    }
    std::sort(sorted_weights.begin(), sorted_weights.end());

    std::vector<IntersectionRecord> intersection_records;
    intersection_records.reserve(links.size());
    for (LinkID id = 0; id < links.size(); ++id)
    {
        const Weight weight = std::get<2>(links[id]);
        const size_t level = static_cast<size_t>(
            std::lower_bound(sorted_weights.begin(), sorted_weights.end(), weight) - sorted_weights.begin());
        intersection_records.push_back({intervals[id], level});
    }

    const LinkID cheapest_link = static_cast<LinkID>(std::min_element(
        links.begin(),
        links.end(),
        [](const auto& left, const auto& right) { return std::get<2>(left) < std::get<2>(right); }) - links.begin());
    const Weight cutoff = sorted_weights.back();

    std::vector<char> explored(links.size(), false);
    std::vector<char> complete(n_nodes, false);
    std::vector<LinkID> predecessor(links.size(), invalid_link);
    std::vector<Weight> distance(links.size());
    std::vector<LinkID> tree_order;
    tree_order.reserve(links.size());
    std::vector<size_t> pops_per_link;
    if constexpr (RecordStatsLevel > 1)
    {
        pops_per_link.resize(links.size());
    }

    auto intersection_index = intersection_index_type.make(intersection_records);
    std::priority_queue<QueueEntry, std::vector<QueueEntry>, std::greater<QueueEntry>> queue;
    explored[cheapest_link] = true;
    predecessor[cheapest_link] = cheapest_link;
    distance[cheapest_link] = std::get<2>(links[cheapest_link]);
    intersection_index->popInterval(cheapest_link);
    queue.emplace(distance[cheapest_link], cheapest_link);
    tree_order.push_back(cheapest_link);

    if constexpr (RecordStatsLevel > 0)
    {
        metrics.sources = 1;
        metrics.links_enqueued = 1;
        metrics.possible_priority_queue_pops = links.size();
    }
    if constexpr (RecordStatsLevel > 1)
    {
        metrics.maximum_queue_size = queue.size();
    }

    size_t n_complete_vertices = 0;
    bool stopped_by_cutoff = false;
    while (!queue.empty())
    {
        const auto [current_distance, current] = queue.top();
        queue.pop();
        if constexpr (RecordStatsLevel > 0)
        {
            metrics.priority_queue_pops++;
        }
        if constexpr (RecordStatsLevel > 1)
        {
            pops_per_link[current]++;
        }

        const auto& [u, v, weight] = links[current];
        for (const NodeID endpoint : {u, v})
        {
            const size_t vertex = static_cast<size_t>(endpoint);
            if (!complete[vertex])
            {
                complete[vertex] = true;
                n_complete_vertices++;
            }
        }

        if (current_distance >= cutoff)
        {
            stopped_by_cutoff = true;
            if constexpr (RecordStatsLevel > 0)
            {
                metrics.termination_by_cutoff++;
            }
            break;
        }

        std::vector<LinkID> neighbours;
        std::vector<LinkID> rejected_by_cutoff;
        const Weight weight_limit = cutoff - current_distance;
        const size_t exclusive_level = static_cast<size_t>(
            std::lower_bound(sorted_weights.begin(), sorted_weights.end(), weight_limit) - sorted_weights.begin());
        intersection_index->forEachIntersection(
            [&](LinkID next, typename BaseIntersectionIdx<RecordStatsLevel>::Interval) {
                const Weight next_weight = std::get<2>(links[next]);
                if (explored[next])
                {
                    if constexpr (RecordStatsLevel > 0)
                    {
                        metrics.intersection_candidates_already_explored++;
                    }
                }
                else if (next_weight >= cutoff - current_distance)
                {
                    if constexpr (RecordStatsLevel > 0)
                    {
                        metrics.intersection_candidates_rejected_by_cutoff++;
                    }
                    rejected_by_cutoff.push_back(next);
                }
                else
                {
                    explored[next] = true;
                    predecessor[next] = current;
                    distance[next] = current_distance + next_weight;
                    neighbours.push_back(next);
                    tree_order.push_back(next);
                    if constexpr (RecordStatsLevel > 0)
                    {
                        metrics.intersection_candidates_enqueued++;
                    }
                }
            },
            intervals[current],
            exclusive_level);

        for (LinkID rejected : rejected_by_cutoff)
        {
            intersection_index->popInterval(rejected);
        }
        for (LinkID next : neighbours)
        {
            intersection_index->popInterval(next);
            queue.emplace(distance[next], next);
            if constexpr (RecordStatsLevel > 0)
            {
                metrics.links_enqueued++;
            }
        }
        if constexpr (RecordStatsLevel > 1)
        {
            metrics.maximum_queue_size = std::max(metrics.maximum_queue_size, queue.size());
        }
    }

    if constexpr (RecordStatsLevel > 0)
    {
        if (!stopped_by_cutoff)
        {
            metrics.termination_by_empty_queue++;
        }
        metrics.intersection_index = intersection_index->emit_metrics();
    }
    if constexpr (RecordStatsLevel > 1)
    {
        metrics.completed_vertices_at_stop = n_complete_vertices;
        metrics.maximum_pops_per_link = *std::max_element(pops_per_link.begin(), pops_per_link.end());
    }

    const Weight infinity = std::numeric_limits<Weight>::max();
    std::vector<char> removable(links.size(), false);
    std::vector<Weight> tree_distance(links.size());
    std::vector<Weight> vertex_distance(n_nodes);
    for (size_t source = 0; source < n_nodes; ++source)
    {
        if (incident[source].empty())
        {
            continue;
        }

        std::fill(tree_distance.begin(), tree_distance.end(), infinity);
        for (LinkID id : incident[source])
        {
            if (explored[id])
            {
                tree_distance[id] = std::get<2>(links[id]);
            }
        }

        for (auto iterator = tree_order.rbegin(); iterator != tree_order.rend(); ++iterator)
        {
            const LinkID current = *iterator;
            const LinkID parent = predecessor[current];
            if (current == parent || tree_distance[current] >= cutoff)
            {
                continue;
            }

            const Weight parent_weight = std::get<2>(links[parent]);
            if (parent_weight < cutoff - tree_distance[current])
            {
                tree_distance[parent] =
                    std::min(tree_distance[parent], tree_distance[current] + parent_weight);
            }
        }

        for (LinkID current : tree_order)
        {
            const LinkID parent = predecessor[current];
            if (current == parent || tree_distance[parent] >= cutoff)
            {
                continue;
            }

            const Weight current_weight = std::get<2>(links[current]);
            if (current_weight < cutoff - tree_distance[parent])
            {
                tree_distance[current] =
                    std::min(tree_distance[current], tree_distance[parent] + current_weight);
            }
        }

        std::fill(vertex_distance.begin(), vertex_distance.end(), infinity);
        for (LinkID id : tree_order)
        {
            const auto& [u, v, weight] = links[id];
            vertex_distance[static_cast<size_t>(u)] =
                std::min(vertex_distance[static_cast<size_t>(u)], tree_distance[id]);
            vertex_distance[static_cast<size_t>(v)] =
                std::min(vertex_distance[static_cast<size_t>(v)], tree_distance[id]);
        }

        for (LinkID id : incident[source])
        {
            const auto& [u, v, weight] = links[id];
            const size_t target = static_cast<size_t>(u) == source ? static_cast<size_t>(v) : static_cast<size_t>(u);
            if (vertex_distance[target] < weight)
            {
                removable[id] = true;
            }
        }
    }

    std::vector<int> result;
    for (LinkID id = 0; id < links.size(); ++id)
    {
        if (removable[id])
        {
            result.push_back(static_cast<int>(id));
        }
    }

    if constexpr (RecordStatsLevel > 0)
    {
        if (output_metrics != nullptr)
        {
            *output_metrics = metrics;
        }
    }
    return result;
}

template<int RecordStatsLevel = 0, typename NodeID, typename Weight>
auto cycle_domination_single_pass(
    const std::vector<std::tuple<NodeID, NodeID, Weight>>& links,
    uint64_t n_nodes,
    CycleReductionMetrics* output_metrics = nullptr)
{
    const BaselineIntersectionIdx<RecordStatsLevel> intersection_index_type;
    return cycle_domination_single_pass<RecordStatsLevel>(links, n_nodes, intersection_index_type, output_metrics);
}

// /*
// template<typename NodeID, typename Weight>
// auto full_cycle_link_domination(const std::vector<std::tuple<NodeID, NodeID, Weight>>& links, uint64_t n_nodes)
// {
//     auto comp = [](const auto& a, const auto& b) {
//         return std::get<0>(a) < std::get<0>(b);
//     };
//     std::priority_queue<std::tuple<Weight, NodeID, NodeID>> link_queue;
//     std::sort(links.begin(), links.end(), comp);
//     for (uint64_t i{}; i < n_nodes; i++)
//     {
//         Weight upper_limit = 0;
//         auto start = std::lower_bound(links.begin(), links.end(), std::make_tuple(i, NodeID{}, Weight{}), comp);
//         auto end = std::upper_bound(links.begin(), links.end(), std::make_tuple(i, NodeID{}, Weight{}), comp);
//         for (auto it = start; it != end; it++)
//         {
//             const auto& [u, v, w] = *it;
//             if (w > upper_limit)
//             {
//                 upper_limit = w;
//             }
//         }
//     }
//     return;
// }
//
// template<typename NodeID, typename Weight>
// struct Link
// {
//     NodeID u;
//     NodeID v;
//     Weight weight;
// };
//
// /**
//  * Returns the input indices of all links l = {a,b} for which there is
//  * an a-to-b chain of touching links with total weight *strictly* less
//  * than weight(l).
//  *
//  * Two links touch if they share an endpoint or cross on the cycle.
//  * i.e l1 = {a,b} and l2 = {c,d} touch if a=c, a=d, b=c, b=d, or a<c<b<d or c<a<d<b.
//  *
//  * Assumptions:
//  *   - 0 <= link.u < link.v < n;
//  *   - links are unique;
//  *   - weights are nonnegative.
//  */
// template<typename NodeID, typename Weight>
// std::vector<int> cycle_link_domination_simple_baseline(int n, const std::vector<Link<NodeID, Weight>>& links)
// {
//     constexpr Weight INF = std::numeric_limits<Weight>::max();
//
//     constexpr Weight NO_TARGET = -1;
//     const size_t INVALID_LINK_ID = std::numeric_limits<size_t>::max();
//
//     const size_t m = links.size();
//
//     /*
//      * Link ids of links incident to each cycle vertex.
//      * TODO: Flatten into csr
//      */
//     std::vector<std::vector<size_t>> incident(n);
//
//     for (size_t id = 0; id < m; ++id)
//     {
//         incident[links[id].u].push_back(id);
//         incident[links[id].v].push_back(id);
//     }
//
//     /*
//      * Link indices sorted by their smaller endpoint.
//      */
//     std::vector<size_t> sortedLinks(m);
//     std::iota(sortedLinks.begin(), sortedLinks.end(), 0ull);
//
//     std::sort(sortedLinks.begin(), sortedLinks.end(), [&](int first, int second) {
//         if (links[first].u != links[second].u)
//         {
//             return links[first].u < links[second].u;
//         }
//
//         return links[first].v < links[second].v;
//     });
//
//     auto touch = [&](size_t first, size_t second) {
//         const CycleLink& a = links[first];
//         const CycleLink& b = links[second];
//
//         const bool shareEndpoint = a.u == b.u || a.u == b.v || a.v == b.u || a.v == b.v;
//
//         const bool cross = (a.u < b.u && b.u < a.v && a.v < b.v) || (b.u < a.u && a.u < b.v && b.v < a.v);
//
//         return shareEndpoint || cross;
//     };
//
//     std::vector<char> removable(m, false);
//
//     std::vector<Weight> linkDistance(m);
//
//     /*
//      * For a fixed source s:
//      *
//      * valueToBeat[v] = weight of the unique link {s,v},
//      *                  or NO_TARGET if there is no such link.
//      *
//      * candidateLink[v] = input index of {s,v}.
//      */
//     std::vector<Weight> valueToBeat(n);
//     std::vector<int> candidateLink(n);
//
//     using QueueEntry = std::pair<Weight, size_t>; // (distance, link id)
//
//     /*
//      * Run the link-graph Dijkstra once from every cycle vertex that is
//      * incident to at least one link.
//      */
//     for (size_t source = 0; source < n; ++source)
//     {
//         if (incident[source].empty())
//         {
//             continue;
//         }
//         // std::fill(linkDistance.begin(), linkDistance.end(), INF);
//         std::fill(valueToBeat.begin(), valueToBeat.end(), NO_TARGET);
//         std::fill(candidateLink.begin(), candidateLink.end(), INVALID_LINK_ID);
//         /*
//          * No distance greater than or equal to cutoff can eliminate
//          * any link incident to this source.
//          */
//         Weight cutoff = 0;
//
//         for (size_t id : incident[source])
//         {
//             // which one is target/source setting both might be faster, but who cares
//             const size_t target;
//             if (links[id].u == source)
//             {
//                 target = links[id].v;
//             }
//             else
//             {
//                 target = links[id].u;
//             }
//
//             valueToBeat[target] = links[id].weight;
//             candidateLink[target] = id;
//
//             cutoff = std::max(cutoff, links[id].weight);
//         }
//
//         std::priority_queue<QueueEntry, std::vector<QueueEntry>, std::greater<QueueEntry>> pq;
//
//         /*
//          * Reaching a link also reaches both of its endpoints.
//          *
//          * The distance need not be finalized before this check: every
//          * tentative label corresponds to an actual link chain, so any
//          * value strictly below valueToBeat[vertex] is already a valid
//          * certificate of removability.
//          * After relaxing an endpoint once, it is actually never necessary to relax it again.
//          * Once all endpoints have been relaxed, the algorithm can terminate early.
//          */
//         auto relaxEndpoint = [&](NodeID vertex, Weight value) {
//             if (valueToBeat[vertex] != NO_TARGET && value < valueToBeat[vertex])
//             {
//                 removable[candidateLink[vertex]] = true;
//             }
//         };
//
//         auto relaxLink = [&](size_t id, Weight value) {
//             /*
//              * Values at least cutoff cannot strictly beat any candidate
//              * link belonging to this source run.
//              */
//             if (value >= cutoff || value >= linkDistance[id])
//             {
//                 return;
//             }
//
//             linkDistance[id] = value;
//             pq.emplace(value, id);
//
//             relaxEndpoint(links[id].u, value);
//             relaxEndpoint(links[id].v, value);
//         };
//
//         /*
//          * The virtual source vertex reaches an incident link by paying
//          * that link's weight.
//          */
//         for (int id : incident[source])
//         {
//             relaxLink(id, links[id].weight);
//         }
//
//         /*
//          * Continue only while the minimum queue value can still beat at
//          * least one direct link incident to the source.
//          */
//         while (!pq.empty() && pq.top().first < cutoff)
//         {
//             const auto [distance, current] = pq.top();
//             pq.pop();
//
//             if (distance != linkDistance[current])
//             {
//                 continue; // Stale queue entry.
//             }
//
//             const CycleLink& currentLink = links[current];
//
//             /*
//              * Examine potential neighbors in increasing order of their
//              * smaller endpoint.
//              *
//              * When next.u > currentLink.v, next cannot cross currentLink
//              * and cannot share an endpoint with it, so the scan can stop.
//              */
//             for (int next : sortedLinks)
//             {
//                 const CycleLink& nextLink = links[next];
//
//                 if (nextLink.u > currentLink.v)
//                 {
//                     break;
//                 }
//
//                 if (next == current || !touch(current, next))
//                 {
//                     continue;
//                 }
//
//                 /*
//                  * Require
//                  *
//                  *     distance + nextLink.weight < cutoff.
//                  *
//                  * This form avoids integer overflow.
//                  */
//                 if (nextLink.weight >= cutoff - distance)
//                 {
//                     continue;
//                 }
//
//                 relaxLink(next, distance + nextLink.weight);
//             }
//         }
//     }
//
//     std::vector<int> result;
//
//     for (int id = 0; id < m; ++id)
//     {
//         if (removable[id])
//         {
//             result.push_back(id);
//         }
//     }
//
//     return result;
// }
//
//
// struct CycleLink
// {
//     int u;
//     int v;
//     std::int64_t weight;
// };
//
// /**
//  * Cycle link reduction algorithm based on bounded dijkstra + segment trees for calculating neighbours on
//  * the cycle intersection graph.
//  */
// namespace cycle_link_reduction_detail
// {
//
//     class MaxReporter
//     {
//     public:
//         explicit MaxReporter(int size) : size_(size), base_(1)
//         {
//             while (base_ < std::max(1, size_))
//             {
//                 base_ <<= 1;
//             }
//             tree_.assign(2 * base_, kNegativeInfinity);
//         }
//
//         // Leaves are ordered by the left endpoint of each link.
//         // A leaf stores the right endpoint of the link if it is active.
//         void build(const std::vector<int>& order, const std::vector<char>& active, const std::vector<int>& right)
//         {
//             std::fill(tree_.begin(), tree_.end(), kNegativeInfinity);
//
//             for (int pos = 0; pos < size_; ++pos)
//             {
//                 const int id = order[pos];
//                 if (active[id])
//                 {
//                     tree_[base_ + pos] = right[id];
//                 }
//             }
//
//             for (int node = base_ - 1; node >= 1; --node)
//             {
//                 tree_[node] = std::max(tree_[2 * node], tree_[2 * node + 1]);
//             }
//         }
//
//         void eraseAt(int position)
//         {
//             int node = base_ + position;
//
//             if (tree_[node] == kNegativeInfinity)
//             {
//                 return;
//             }
//
//             tree_[node] = kNegativeInfinity;
//
//             for (node >>= 1; node >= 1; node >>= 1)
//             {
//                 const int newValue = std::max(tree_[2 * node], tree_[2 * node + 1]);
//
//                 if (tree_[node] == newValue)
//                 {
//                     break;
//                 }
//
//                 tree_[node] = newValue;
//             }
//         }
//
//         // Returns a position in [queryLeft, queryRight) whose stored value is
//         // strictly larger than threshold, or -1 if none exists.
//         int findPosition(int queryLeft, int queryRight, int threshold) const
//         {
//             if (queryLeft >= queryRight || tree_[1] <= threshold)
//             {
//                 return -1;
//             }
//
//             return findPositionRecursive(1, 0, base_, queryLeft, queryRight, threshold);
//         }
//
//     private:
//         static constexpr int kNegativeInfinity = -1;
//
//         int size_;
//         int base_;
//         std::vector<int> tree_;
//
//         int
//         findPositionRecursive(int node, int nodeLeft, int nodeRight, int queryLeft, int queryRight, int threshold)
//         const
//         {
//             if (nodeRight <= queryLeft || queryRight <= nodeLeft || tree_[node] <= threshold)
//             {
//                 return -1;
//             }
//
//             if (nodeRight - nodeLeft == 1)
//             {
//                 return nodeLeft < size_ ? nodeLeft : -1;
//             }
//
//             const int middle = nodeLeft + (nodeRight - nodeLeft) / 2;
//
//             int answer = findPositionRecursive(2 * node, nodeLeft, middle, queryLeft, queryRight, threshold);
//
//             if (answer != -1)
//             {
//                 return answer;
//             }
//
//             return findPositionRecursive(2 * node + 1, middle, nodeRight, queryLeft, queryRight, threshold);
//         }
//     };
//
//     class MinReporter
//     {
//     public:
//         explicit MinReporter(int size) : size_(size), base_(1)
//         {
//             while (base_ < std::max(1, size_))
//             {
//                 base_ <<= 1;
//             }
//             tree_.assign(2 * base_, kInfinity);
//         }
//
//         // Leaves are ordered by the right endpoint of each link.
//         // A leaf stores the left endpoint of the link if it is active.
//         void build(const std::vector<int>& order, const std::vector<char>& active, const std::vector<int>& left)
//         {
//             std::fill(tree_.begin(), tree_.end(), kInfinity);
//
//             for (int pos = 0; pos < size_; ++pos)
//             {
//                 const int id = order[pos];
//                 if (active[id])
//                 {
//                     tree_[base_ + pos] = left[id];
//                 }
//             }
//
//             for (int node = base_ - 1; node >= 1; --node)
//             {
//                 tree_[node] = std::min(tree_[2 * node], tree_[2 * node + 1]);
//             }
//         }
//
//         void eraseAt(int position)
//         {
//             int node = base_ + position;
//
//             if (tree_[node] == kInfinity)
//             {
//                 return;
//             }
//
//             tree_[node] = kInfinity;
//
//             for (node >>= 1; node >= 1; node >>= 1)
//             {
//                 const int newValue = std::min(tree_[2 * node], tree_[2 * node + 1]);
//
//                 if (tree_[node] == newValue)
//                 {
//                     break;
//                 }
//
//                 tree_[node] = newValue;
//             }
//         }
//
//         // Returns a position in [queryLeft, queryRight) whose stored value is
//         // strictly smaller than threshold, or -1 if none exists.
//         int findPosition(int queryLeft, int queryRight, int threshold) const
//         {
//             if (queryLeft >= queryRight || tree_[1] >= threshold)
//             {
//                 return -1;
//             }
//
//             return findPositionRecursive(1, 0, base_, queryLeft, queryRight, threshold);
//         }
//
//     private:
//         static constexpr int kInfinity = std::numeric_limits<int>::max();
//
//         int size_;
//         int base_;
//         std::vector<int> tree_;
//
//         int
//         findPositionRecursive(int node, int nodeLeft, int nodeRight, int queryLeft, int queryRight, int threshold)
//         const
//         {
//             if (nodeRight <= queryLeft || queryRight <= nodeLeft || tree_[node] >= threshold)
//             {
//                 return -1;
//             }
//
//             if (nodeRight - nodeLeft == 1)
//             {
//                 return nodeLeft < size_ ? nodeLeft : -1;
//             }
//
//             const int middle = nodeLeft + (nodeRight - nodeLeft) / 2;
//
//             int answer = findPositionRecursive(2 * node, nodeLeft, middle, queryLeft, queryRight, threshold);
//
//             if (answer != -1)
//             {
//                 return answer;
//             }
//
//             return findPositionRecursive(2 * node + 1, middle, nodeRight, queryLeft, queryRight, threshold);
//         }
//     };
//
// } // namespace cycle_link_reduction_detail
//
// /**
//  * Returns the indices of all strictly removable links.
//  *
//  * A link l = {a,b} is returned exactly when there is an a-to-b chain
//  * of links having total weight strictly smaller than l.weight.
//  *
//  * Consecutive links in a chain must either:
//  *   - cross as chords of the cycle, or
//  *   - share an endpoint.
//  *
//  * Requirements:
//  *   - vertices are numbered 0,...,n-1 in cyclic order;
//  *   - link weights are nonnegative;
//  *   - links have distinct endpoints.
//  *
//  * The returned indices are in increasing input order.
//  */
// std::vector<int> removableCycleLinks(int n, const std::vector<CycleLink>& links)
// {
//     using namespace cycle_link_reduction_detail;
//     using Distance = std::int64_t;
//
//     constexpr Distance kInfinity = std::numeric_limits<Distance>::max();
//
//     const int numberOfLinks = static_cast<int>(links.size());
//
//     if (n < 3)
//     {
//         throw std::invalid_argument("A cycle must contain at least three vertices.");
//     }
//
//     if (numberOfLinks == 0)
//     {
//         return {};
//     }
//
//     /*
//      * Coordinate-compress only vertices appearing as link endpoints.
//      * Thus the storage is O(numberOfLinks), even when n is very large.
//      */
//     std::vector<int> coordinates;
//     coordinates.reserve(2 * numberOfLinks);
//
//     for (const CycleLink& link : links)
//     {
//         if (link.u < 0 || link.u >= n || link.v < 0 || link.v >= n)
//         {
//             throw std::invalid_argument("A link endpoint is outside [0,n).");
//         }
//
//         if (link.u == link.v)
//         {
//             throw std::invalid_argument("Self-links are not supported.");
//         }
//
//         if (link.weight < 0)
//         {
//             throw std::invalid_argument("Link weights must be nonnegative.");
//         }
//
//         coordinates.push_back(link.u);
//         coordinates.push_back(link.v);
//     }
//
//     std::sort(coordinates.begin(), coordinates.end());
//     coordinates.erase(std::unique(coordinates.begin(), coordinates.end()), coordinates.end());
//
//     const int numberOfEndpointVertices = static_cast<int>(coordinates.size());
//
//     std::vector<int> left(numberOfLinks);
//     std::vector<int> right(numberOfLinks);
//     std::vector<Distance> weight(numberOfLinks);
//     std::vector<int> degree(numberOfEndpointVertices, 0);
//
//     for (int id = 0; id < numberOfLinks; ++id)
//     {
//         int a = static_cast<int>(
//             std::lower_bound(coordinates.begin(), coordinates.end(), links[id].u) - coordinates.begin());
//
//         int b = static_cast<int>(
//             std::lower_bound(coordinates.begin(), coordinates.end(), links[id].v) - coordinates.begin());
//
//         if (a > b)
//         {
//             std::swap(a, b);
//         }
//
//         left[id] = a;
//         right[id] = b;
//         weight[id] = links[id].weight;
//
//         ++degree[a];
//         ++degree[b];
//     }
//
//     /*
//      * Store endpoint-to-link incidence lists in CSR form.
//      * This represents shared-endpoint adjacencies without materializing
//      * the corresponding cliques.
//      */
//     std::vector<int> incidenceStart(numberOfEndpointVertices + 1, 0);
//
//     for (int vertex = 0; vertex < numberOfEndpointVertices; ++vertex)
//     {
//         incidenceStart[vertex + 1] = incidenceStart[vertex] + degree[vertex];
//     }
//
//     std::vector<int> incidence(2 * numberOfLinks);
//     std::vector<int> incidenceCursor = incidenceStart;
//
//     for (int id = 0; id < numberOfLinks; ++id)
//     {
//         incidence[incidenceCursor[left[id]]++] = id;
//         incidence[incidenceCursor[right[id]]++] = id;
//     }
//
//     /*
//      * The two orderings support implicit crossing-neighbor reporting.
//      */
//     std::vector<int> byLeft(numberOfLinks);
//     std::vector<int> byRight(numberOfLinks);
//     std::vector<int> positionByLeft(numberOfLinks);
//     std::vector<int> positionByRight(numberOfLinks);
//
//     for (int id = 0; id < numberOfLinks; ++id)
//     {
//         byLeft[id] = id;
//         byRight[id] = id;
//     }
//
//     std::sort(byLeft.begin(), byLeft.end(), [&](int first, int second) {
//         return std::tie(left[first], right[first], first) < std::tie(left[second], right[second], second);
//     });
//
//     std::sort(byRight.begin(), byRight.end(), [&](int first, int second) {
//         return std::tie(right[first], left[first], first) < std::tie(right[second], left[second], second);
//     });
//
//     for (int position = 0; position < numberOfLinks; ++position)
//     {
//         positionByLeft[byLeft[position]] = position;
//         positionByRight[byRight[position]] = position;
//     }
//
//     // First position in byLeft whose left endpoint is > value.
//     auto upperLeft = [&](int value) {
//         int low = 0;
//         int high = numberOfLinks;
//
//         while (low < high)
//         {
//             const int middle = low + (high - low) / 2;
//
//             if (left[byLeft[middle]] <= value)
//             {
//                 low = middle + 1;
//             }
//             else
//             {
//                 high = middle;
//             }
//         }
//
//         return low;
//     };
//
//     // First position in byLeft whose left endpoint is >= value.
//     auto lowerLeft = [&](int value) {
//         int low = 0;
//         int high = numberOfLinks;
//
//         while (low < high)
//         {
//             const int middle = low + (high - low) / 2;
//
//             if (left[byLeft[middle]] < value)
//             {
//                 low = middle + 1;
//             }
//             else
//             {
//                 high = middle;
//             }
//         }
//
//         return low;
//     };
//
//     // First position in byRight whose right endpoint is > value.
//     auto upperRight = [&](int value) {
//         int low = 0;
//         int high = numberOfLinks;
//
//         while (low < high)
//         {
//             const int middle = low + (high - low) / 2;
//
//             if (right[byRight[middle]] <= value)
//             {
//                 low = middle + 1;
//             }
//             else
//             {
//                 high = middle;
//             }
//         }
//
//         return low;
//     };
//
//     // First position in byRight whose right endpoint is >= value.
//     auto lowerRight = [&](int value) {
//         int low = 0;
//         int high = numberOfLinks;
//
//         while (low < high)
//         {
//             const int middle = low + (high - low) / 2;
//
//             if (right[byRight[middle]] < value)
//             {
//                 low = middle + 1;
//             }
//             else
//             {
//                 high = middle;
//             }
//         }
//
//         return low;
//     };
//
//     /*
//      * One shortest-path run from a vertex tests every link assigned to
//      * that vertex. Assigning a link to its higher-degree endpoint is a
//      * heuristic that reduces the number of source runs on instances
//      * such as stars.
//      */
//     auto sourceOf = [&](int id) {
//         if (degree[left[id]] != degree[right[id]])
//         {
//             return degree[left[id]] > degree[right[id]] ? left[id] : right[id];
//         }
//
//         return left[id];
//     };
//
//     auto targetOf = [&](int id, int source) {
//         return left[id] == source ? right[id] : left[id];
//     };
//
//     std::vector<int> testOrder(numberOfLinks);
//
//     for (int id = 0; id < numberOfLinks; ++id)
//     {
//         testOrder[id] = id;
//     }
//
//     std::sort(testOrder.begin(), testOrder.end(), [&](int first, int second) {
//         const int firstSource = sourceOf(first);
//         const int secondSource = sourceOf(second);
//
//         if (firstSource != secondSource)
//         {
//             return firstSource < secondSource;
//         }
//
//         return first < second;
//     });
//
//     MaxReporter maxRightReporter(numberOfLinks);
//     MinReporter minLeftReporter(numberOfLinks);
//
//     std::vector<char> active(numberOfLinks, 0);
//     std::vector<Distance> endpointDistance(numberOfEndpointVertices, kInfinity);
//     std::vector<char> removable(numberOfLinks, 0);
//
//     int groupBegin = 0;
//
//     while (groupBegin < numberOfLinks)
//     {
//         const int source = sourceOf(testOrder[groupBegin]);
//
//         int groupEnd = groupBegin;
//         Distance cutoff = 0;
//
//         while (groupEnd < numberOfLinks && sourceOf(testOrder[groupEnd]) == source)
//         {
//             cutoff = std::max(cutoff, weight[testOrder[groupEnd]]);
//             ++groupEnd;
//         }
//
//         /*
//          * We only need paths with total weight < cutoff.
//          * Therefore links of weight >= cutoff cannot be part of a
//          * relevant path.
//          */
//         if (cutoff > 0)
//         {
//             for (int id = 0; id < numberOfLinks; ++id)
//             {
//                 active[id] = static_cast<char>(weight[id] < cutoff);
//             }
//
//             maxRightReporter.build(byLeft, active, right);
//
//             minLeftReporter.build(byRight, active, left);
//
//             std::fill(endpointDistance.begin(), endpointDistance.end(), kInfinity);
//
//             using QueueItem = std::pair<Distance, int>;
//
//             std::priority_queue<QueueItem, std::vector<QueueItem>, std::greater<QueueItem>> queue;
//
//             /*
//              * The first settled neighbor that discovers a link gives
//              * its final shortest-path label: every transition entering
//              * link id adds exactly weight[id].
//              */
//             auto discover = [&](int id, Distance predecessorDistance) {
//                 if (!active[id])
//                 {
//                     return;
//                 }
//
//                 active[id] = 0;
//
//                 maxRightReporter.eraseAt(positionByLeft[id]);
//
//                 minLeftReporter.eraseAt(positionByRight[id]);
//
//                 /*
//                  * Only distances strictly below cutoff matter.
//                  * Writing the test this way also avoids overflow.
//                  */
//                 if (weight[id] < cutoff - predecessorDistance)
//                 {
//                     queue.emplace(predecessorDistance + weight[id], id);
//                 }
//             };
//
//             /*
//              * Expanding an endpoint implicitly traverses all
//              * shared-endpoint adjacencies.
//              */
//             auto expandEndpoint = [&](int vertex, Distance distance) {
//                 if (endpointDistance[vertex] != kInfinity)
//                 {
//                     return;
//                 }
//
//                 endpointDistance[vertex] = distance;
//
//                 for (int position = incidenceStart[vertex]; position < incidenceStart[vertex + 1]; ++position)
//                 {
//                     discover(incidence[position], distance);
//                 }
//             };
//
//             expandEndpoint(source, 0);
//
//             while (!queue.empty())
//             {
//                 const auto [distance, current] = queue.top();
//                 queue.pop();
//
//                 const int p = left[current];
//                 const int q = right[current];
//
//                 expandEndpoint(p, distance);
//                 expandEndpoint(q, distance);
//
//                 /*
//                  * Report every undiscovered link f satisfying
//                  *
//                  *     p < left[f] < q < right[f].
//                  */
//                 const int firstLeftPosition = upperLeft(p);
//                 const int lastLeftPosition = lowerLeft(q);
//
//                 while (true)
//                 {
//                     const int position = maxRightReporter.findPosition(firstLeftPosition, lastLeftPosition, q);
//
//                     if (position == -1)
//                     {
//                         break;
//                     }
//
//                     discover(byLeft[position], distance);
//                 }
//
//                 /*
//                  * Report every undiscovered link f satisfying
//                  *
//                  *     left[f] < p < right[f] < q.
//                  */
//                 const int firstRightPosition = upperRight(p);
//                 const int lastRightPosition = lowerRight(q);
//
//                 while (true)
//                 {
//                     const int position = minLeftReporter.findPosition(firstRightPosition, lastRightPosition, p);
//
//                     if (position == -1)
//                     {
//                         break;
//                     }
//
//                     discover(byRight[position], distance);
//                 }
//             }
//
//             for (int position = groupBegin; position < groupEnd; ++position)
//             {
//                 const int id = testOrder[position];
//                 const int target = targetOf(id, source);
//
//                 if (endpointDistance[target] < weight[id])
//                 {
//                     removable[id] = 1;
//                 }
//             }
//         }
//
//         groupBegin = groupEnd;
//     }
//
//     std::vector<int> answer;
//
//     for (int id = 0; id < numberOfLinks; ++id)
//     {
//         if (removable[id])
//         {
//             answer.push_back(id);
//         }
//     }
//
//     return answer;
// }
// //
