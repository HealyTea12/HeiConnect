#include <gtest/gtest.h>

#include "HeiConnect/conn_aug/reducers/common.hpp"
#include "HeiConnect/conn_aug/reducers/cycle_reducer.hpp"
#include "HeiConnect/data_structures/intersection_index/intersection_tree.hpp"
#include "HeiConnect/data_structures/graph_utils.hpp"

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
