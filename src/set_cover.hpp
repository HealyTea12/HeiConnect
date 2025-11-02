#pragma once
#include "graph.hpp"

namespace solver
{
    std::list<graph::Edge> set_cover(const graph::GraphPair &g);
}