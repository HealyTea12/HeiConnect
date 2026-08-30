#include <gtest/gtest.h>

#include "HeiConnect/conn_aug/reducers/common.hpp"
#include "HeiConnect/conn_aug/reducers/cycle_reducer.hpp"
#include "HeiConnect/data_structures/intersection_index/intersection_tree.hpp"
#include "HeiConnect/data_structures/graph_utils.hpp"

#include <limits>
#include <queue>
#include <random>

static bool cycle_links_touch(
    const std::tuple<int, int, int>& left,
    const std::tuple<int, int, int>& right)
{
    const auto [a, b, left_weight] = left;
    const auto [c, d, right_weight] = right;
    return a == c || a == d || b == c || b == d || (a < c && c < b && b < d) ||
        (c < a && a < d && d < b);
}

static std::vector<int> brute_force_cycle_dominations(const std::vector<std::tuple<int, int, int>>& links)
{
    using QueueEntry = std::pair<int, size_t>;
    const int infinity = std::numeric_limits<int>::max();
    std::vector<int> result;

    for (size_t target = 0; target < links.size(); ++target)
    {
        const auto [source_vertex, target_vertex, target_weight] = links[target];
        std::vector<int> distance(links.size(), infinity);
        std::priority_queue<QueueEntry, std::vector<QueueEntry>, std::greater<QueueEntry>> queue;
        for (size_t id = 0; id < links.size(); ++id)
        {
            if (id == target)
            {
                continue;
            }
            const auto [u, v, weight] = links[id];
            if (u == source_vertex || v == source_vertex)
            {
                distance[id] = weight;
                queue.emplace(weight, id);
            }
        }

        int replacement_weight = infinity;
        while (!queue.empty())
        {
            const auto [current_weight, current] = queue.top();
            queue.pop();
            if (current_weight != distance[current])
            {
                continue;
            }

            const auto [u, v, weight] = links[current];
            if (u == target_vertex || v == target_vertex)
            {
                replacement_weight = current_weight;
                break;
            }

            for (size_t next = 0; next < links.size(); ++next)
            {
                if (next == target || next == current || !cycle_links_touch(links[current], links[next]))
                {
                    continue;
                }
                const int next_weight = std::get<2>(links[next]);
                if (current_weight <= target_weight && next_weight <= target_weight - current_weight &&
                    current_weight + next_weight < distance[next])
                {
                    distance[next] = current_weight + next_weight;
                    queue.emplace(distance[next], next);
                }
            }
        }

        if (replacement_weight <= target_weight)
        {
            result.push_back(static_cast<int>(target));
        }
    }
    return result;
}

TEST(ConnAugReducers, MergesCrossingContractionsOnCycle)
{
    auto graph = create_cycle_graph_undirected(4);
    const auto link_graph = WeightedCRFGraph<>::vec_links_to_csr({{0, 2, 1.0}, {1, 3, 1.0}}, 4);
    const std::vector<size_t> selected_links{0, 1};
    UnionFind uf(graph.num_vertices());

    const auto stats = add_links_to_union_find(graph, link_graph, uf, selected_links);

    EXPECT_EQ(uf.find(0), uf.find(1));
    EXPECT_EQ(uf.find(0), uf.find(2));
    EXPECT_EQ(uf.find(0), uf.find(3));
    EXPECT_EQ(stats.total_merged_nodes, 3);
}

TEST(ConnAugReducers, KeepsNonCrossingContractionsSeparateOnCycle)
{
    auto graph = create_cycle_graph_undirected(4);
    const auto link_graph = WeightedCRFGraph<>::vec_links_to_csr({{0, 2, 1.0}}, 4);
    const std::vector<size_t> selected_links{0};
    UnionFind uf(graph.num_vertices());

    const auto stats = add_links_to_union_find(graph, link_graph, uf, selected_links);

    EXPECT_EQ(uf.find(0), uf.find(2));
    EXPECT_NE(uf.find(0), uf.find(1));
    EXPECT_NE(uf.find(0), uf.find(3));
    EXPECT_EQ(stats.total_merged_nodes, 1);
}

TEST(ConnAugReducers, KeepsNestedContractionsSeparateOnCycle)
{
    auto graph = create_cycle_graph_undirected(4);
    const auto link_graph = WeightedCRFGraph<>::vec_links_to_csr({{0, 3, 1.0}, {1, 2, 1.0}}, 4);
    const std::vector<size_t> selected_links{0, 1};
    UnionFind uf(graph.num_vertices());

    const auto stats = add_links_to_union_find(graph, link_graph, uf, selected_links);

    EXPECT_EQ(uf.find(0), uf.find(3));
    EXPECT_EQ(uf.find(1), uf.find(2));
    EXPECT_NE(uf.find(0), uf.find(1));
    EXPECT_EQ(stats.total_merged_nodes, 2);
}

TEST(ConnAugReducers, FrozenStackMergesCrossingContractionsOnCycle)
{
    auto graph = create_cycle_graph_undirected(4);
    const auto link_graph = WeightedCRFGraph<>::vec_links_to_csr({{0, 2, 1.0}, {1, 3, 1.0}}, 4);
    const std::vector<size_t> selected_links{0, 1};
    UnionFind uf(graph.num_vertices());

    const auto stats = add_links_to_union_find_frozen_stack(graph, link_graph, uf, selected_links);

    EXPECT_EQ(uf.find(0), uf.find(1));
    EXPECT_EQ(uf.find(0), uf.find(2));
    EXPECT_EQ(uf.find(0), uf.find(3));
    EXPECT_EQ(stats.total_merged_nodes, 3);
}

TEST(ConnAugReducers, FrozenStackKeepsNonCrossingContractionsSeparateOnCycle)
{
    auto graph = create_cycle_graph_undirected(4);
    const auto link_graph = WeightedCRFGraph<>::vec_links_to_csr({{0, 2, 1.0}}, 4);
    const std::vector<size_t> selected_links{0};
    UnionFind uf(graph.num_vertices());

    const auto stats = add_links_to_union_find_frozen_stack(graph, link_graph, uf, selected_links);

    EXPECT_EQ(uf.find(0), uf.find(2));
    EXPECT_NE(uf.find(0), uf.find(1));
    EXPECT_NE(uf.find(0), uf.find(3));
    EXPECT_EQ(stats.total_merged_nodes, 1);
}

TEST(ConnAugReducers, FrozenStackKeepsNestedContractionsSeparateOnCycle)
{
    auto graph = create_cycle_graph_undirected(4);
    const auto link_graph = WeightedCRFGraph<>::vec_links_to_csr({{0, 3, 1.0}, {1, 2, 1.0}}, 4);
    const std::vector<size_t> selected_links{0, 1};
    UnionFind uf(graph.num_vertices());

    const auto stats = add_links_to_union_find_frozen_stack(graph, link_graph, uf, selected_links);

    EXPECT_EQ(uf.find(0), uf.find(3));
    EXPECT_EQ(uf.find(1), uf.find(2));
    EXPECT_NE(uf.find(0), uf.find(1));
    EXPECT_EQ(stats.total_merged_nodes, 2);
}

TEST(ConnAugReducers, MaterializedContractionsPreserveParallelLinks)
{
    const auto graph = WeightedCRFGraph<>::vec_links_to_csr(
        {{0, 1, 1.0}, {1, 0, 1.0}, {1, 2, 1.0}, {2, 1, 1.0}},
        3);
    const auto link_graph = WeightedCRFGraph<>::vec_links_to_csr(
        {{0, 2, 2.0}, {0, 1, 1.0}, {1, 2, 3.0}},
        3);
    UnionFind uf(graph.num_vertices());
    uf.unite(0, 1);

    auto [contracted_graph, contracted_link_graph, original_link_ids] =
        materialize_contractions_preserving_links(graph, link_graph, uf);

    EXPECT_EQ(contracted_graph.num_vertices(), 2);
    ASSERT_EQ(contracted_link_graph.num_edges(), 2);
    EXPECT_EQ(contracted_link_graph.graph.edges[0], contracted_link_graph.graph.edges[1]);
    EXPECT_EQ(original_link_ids, (std::vector<size_t>{0, 2}));
}

TEST(ConnAugReducers, IntersectionTreeMatchesBaselineCycleReduction)
{
    std::vector<std::tuple<int, int, int>> links;
    constexpr int cycle_size = 12;
    for (int u = 0; u < cycle_size; ++u)
    {
        for (int v = u + 1; v < cycle_size; ++v)
        {
            links.emplace_back(u, v, (u + v) % cycle_size + 1);
        }
    }

    const BaselineIntersectionIdx<0> baseline;
    const IntersectionTreeIdx<0> intersection_tree;
    const WeightedIntersectionTreeIdx<0> weighted_intersection_tree;
    const auto baseline_result = cycle_domination_baseline(links, cycle_size, baseline);
    EXPECT_EQ(
        baseline_result,
        cycle_domination_baseline(links, cycle_size, intersection_tree));
    EXPECT_EQ(
        baseline_result,
        cycle_domination_baseline(links, cycle_size, weighted_intersection_tree));
}

TEST(ConnAugReducers, GlobalSweepMatchesExactOracleOnFullLinkGraph)
{
    std::vector<std::tuple<int, int, int>> links;
    constexpr int cycle_size = 12;
    for (int u = 0; u < cycle_size; ++u)
    {
        for (int v = u + 1; v < cycle_size; ++v)
        {
            links.emplace_back(u, v, (u + v) % cycle_size + 1);
        }
    }

    const auto exact = brute_force_cycle_dominations(links);
    const IntersectionTreeIdx<0> intersection_tree;
    const WeightedIntersectionTreeIdx<0> weighted_intersection_tree;

    EXPECT_EQ(cycle_domination_global_sweep(links, cycle_size), exact);
    EXPECT_EQ(cycle_domination_global_sweep(links, cycle_size, intersection_tree), exact);
    EXPECT_EQ(cycle_domination_global_sweep(links, cycle_size, weighted_intersection_tree), exact);
}

TEST(ConnAugReducers, GlobalSweepMatchesStrictBaselineWhenEqualityIsImpossible)
{
    std::vector<std::tuple<int, int, int>> links;
    const WeightedIntersectionTreeIdx<0> weighted_intersection_tree;
    constexpr int cycle_size = 6;
    int weight = 1;
    for (int u = 0; u < cycle_size; ++u)
    {
        for (int v = u + 1; v < cycle_size; ++v)
        {
            links.emplace_back(u, v, weight);
            weight *= 2;
        }
    }

    const auto exact = brute_force_cycle_dominations(links);
    EXPECT_EQ(cycle_domination_baseline(links, cycle_size), exact);
    EXPECT_EQ(cycle_domination_global_sweep(links, cycle_size, weighted_intersection_tree), exact);
}

TEST(ConnAugReducers, GlobalSweepMatchesStrictBaselineOnSparseCrossingLinks)
{
    constexpr int cycle_size = 8;
    constexpr int half_cycle = cycle_size / 2;
    std::vector<std::tuple<int, int, int>> links;
    for (int u = 0; u < half_cycle; ++u)
    {
        links.emplace_back(u, u + half_cycle, u + 1);
        links.emplace_back(u, half_cycle + (u + 1) % half_cycle, cycle_size * cycle_size + u + 1);
    }
    const auto exact = brute_force_cycle_dominations(links);
    const WeightedIntersectionTreeIdx<0> weighted_intersection_tree;

    EXPECT_EQ(cycle_domination_baseline(links, cycle_size, weighted_intersection_tree), exact);
    EXPECT_EQ(cycle_domination_global_sweep(links, cycle_size, weighted_intersection_tree), exact);
}

TEST(ConnAugReducers, GlobalSweepMatchesExactOracleOnRandomInstances)
{
    std::mt19937 random_engine(123456);
    const WeightedIntersectionTreeIdx<0> weighted_intersection_tree;
    for (int cycle_size = 3; cycle_size <= 10; ++cycle_size)
    {
        std::bernoulli_distribution include_link(0.65);
        std::uniform_int_distribution<int> link_weight(1, 20);
        for (int instance = 0; instance < 100; ++instance)
        {
            std::vector<std::tuple<int, int, int>> links;
            for (int u = 0; u < cycle_size; ++u)
            {
                for (int v = u + 1; v < cycle_size; ++v)
                {
                    if (include_link(random_engine))
                    {
                        links.emplace_back(u, v, link_weight(random_engine));
                    }
                }
            }

            const auto exact = brute_force_cycle_dominations(links);
            EXPECT_EQ(cycle_domination_global_sweep(links, cycle_size), exact);
            EXPECT_EQ(cycle_domination_global_sweep(links, cycle_size, weighted_intersection_tree), exact);
        }
    }
}

TEST(ConnAugReducers, GlobalSweepMatchesExactOracleExhaustivelyOnFourVertices)
{
    const std::array<std::pair<int, int>, 6> endpoints{
        std::pair{0, 1},
        std::pair{0, 2},
        std::pair{0, 3},
        std::pair{1, 2},
        std::pair{1, 3},
        std::pair{2, 3}};

    constexpr int assignments = 4 * 4 * 4 * 4 * 4 * 4;
    for (int assignment = 0; assignment < assignments; ++assignment)
    {
        int encoded = assignment;
        std::vector<std::tuple<int, int, int>> links;
        for (const auto [u, v] : endpoints)
        {
            const int weight = encoded % 4;
            encoded /= 4;
            if (weight != 0)
            {
                links.emplace_back(u, v, weight);
            }
        }

        EXPECT_EQ(cycle_domination_global_sweep(links, 4), brute_force_cycle_dominations(links));
    }
}

TEST(ConnAugReducers, GlobalSweepFindsDisconnectedEqualCostWitness)
{
    const std::vector<std::tuple<int, int, int>> links{{4, 5, 1}, {0, 1, 4}, {0, 2, 2}, {1, 3, 2}};
    const auto exact = brute_force_cycle_dominations(links);
    const IntersectionTreeIdx<0> intersection_tree;
    const WeightedIntersectionTreeIdx<0> weighted_intersection_tree;

    EXPECT_EQ(exact, (std::vector<int>{1}));
    EXPECT_EQ(cycle_domination_global_sweep(links, 6), exact);
    EXPECT_EQ(cycle_domination_global_sweep(links, 6, intersection_tree), exact);
    EXPECT_EQ(cycle_domination_global_sweep(links, 6, weighted_intersection_tree), exact);
}

TEST(ConnAugReducers, GlobalSweepRecordsRecursiveEqualCostWitness)
{
    std::vector<std::tuple<int, int, int>> links{{0, 1, 3}, {0, 4, 1}, {2, 5, 1}, {1, 3, 1}};
    std::sort(links.begin(), links.end());

    do
    {
        EXPECT_EQ(cycle_domination_global_sweep(links, 6), brute_force_cycle_dominations(links));
    }
    while (std::next_permutation(links.begin(), links.end()));
}

TEST(ConnAugReducers, GlobalSweepDoesNotTreatSingletonAsEqualCostReplacement)
{
    const std::vector<std::tuple<int, int, int>> links{{0, 1, 4}};

    EXPECT_TRUE(cycle_domination_global_sweep(links, 3).empty());
}

TEST(ConnAugReducers, GlobalSweepRecordsSharedEndpointEqualCostWitness)
{
    const std::vector<std::tuple<int, int, int>> links{{0, 4, 4}, {0, 2, 2}, {2, 4, 2}};
    const WeightedIntersectionTreeIdx<0> weighted_intersection_tree;

    EXPECT_EQ(cycle_domination_global_sweep(links, 5), (std::vector<int>{0}));
    EXPECT_EQ(cycle_domination_global_sweep(links, 5, weighted_intersection_tree), (std::vector<int>{0}));
}

TEST(ConnAugReducers, GlobalSweepRejectsNestedNoncrossingWitness)
{
    const std::vector<std::tuple<int, int, int>> links{{0, 3, 4}, {0, 4, 2}, {1, 3, 2}};

    EXPECT_TRUE(cycle_domination_global_sweep(links, 5).empty());
}

TEST(ConnAugReducers, GlobalSweepHandlesWeightCutoffWithoutOverflow)
{
    const auto maximum_weight = std::numeric_limits<uint64_t>::max();
    const std::vector<std::tuple<int, int, uint64_t>> equal_links{
        {0, 1, maximum_weight},
        {0, 2, maximum_weight - 2},
        {1, 3, 2}};
    const std::vector<std::tuple<int, int, uint64_t>> overflowing_links{
        {0, 1, maximum_weight},
        {0, 2, maximum_weight - 1},
        {1, 3, 2}};

    EXPECT_EQ(cycle_domination_global_sweep(equal_links, 4), (std::vector<int>{0}));
    EXPECT_TRUE(cycle_domination_global_sweep(overflowing_links, 4).empty());
}

TEST(ConnAugReducers, GlobalSweepDoesNotStopWhenAllCycleVerticesAreReached)
{
    const std::vector<std::tuple<int, int, int>> links{{0, 1, 1}, {1, 2, 2}, {0, 2, 10}};
    CycleReductionMetrics metrics;

    cycle_domination_global_sweep<2>(links, 3, &metrics);

    EXPECT_EQ(metrics.sources, 1);
    EXPECT_EQ(metrics.termination_by_completion, 0);
    EXPECT_EQ(metrics.termination_by_cutoff + metrics.termination_by_empty_queue, 1);
    EXPECT_EQ(metrics.completed_vertices_at_stop, 3);
}

TEST(ConnAugReducers, RecordsDetailedCycleReductionMetrics)
{
    const std::vector<std::tuple<int, int, int>> links{{0, 2, 1}, {1, 3, 2}, {0, 3, 3}};
    const WeightedIntersectionTreeIdx<2> intersection_tree;
    CycleReductionMetrics metrics;

    cycle_domination_baseline<2>(links, 4, intersection_tree, &metrics);

    EXPECT_EQ(metrics.sources, 4);
    EXPECT_GT(metrics.priority_queue_pops, 0);
    EXPECT_LE(metrics.priority_queue_pops, metrics.possible_priority_queue_pops);
    EXPECT_GT(metrics.intersection_index.queries, 0);
    EXPECT_GT(
        metrics.intersection_index.candidates_inspected + metrics.intersection_index.subtrees_pruned_by_level,
        0);
    EXPECT_GT(metrics.maximum_pops_per_link, 0);
    EXPECT_EQ(
        metrics.intersection_index.callbacks,
        metrics.intersection_candidates_enqueued + metrics.intersection_candidates_already_explored +
            metrics.intersection_candidates_rejected_by_cutoff);
    EXPECT_EQ(metrics.intersection_candidates_already_explored, 0);
    EXPECT_EQ(metrics.intersection_candidates_rejected_by_cutoff, 0);
}
