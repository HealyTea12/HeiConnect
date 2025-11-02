#pragma once

#include <vector>

#include "set_cover/graph.hpp"

// Stoer-Wagner global min-cut algorithm
// Assumptions: undirected graph, no self-loops, positive edge weights
// These assumptions are not checked, user must ensure they are fulfilled to ensure correctness
double global_mincut_simple(const WeightedCRFGraph &graph);