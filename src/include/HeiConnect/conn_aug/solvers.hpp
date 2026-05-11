#pragma once

#include <functional>
#include <utility>

#include "HeiConnect/conn_aug/connectivity_augmentation.hpp"
#include "HeiConnect/set_cover/common.hpp"

namespace HeiConnect::conn_aug
{
    template <typename ConnAugSCReducer, typename InstanceT>
    concept ConnectivityReductionFunction = requires(ConnAugSCReducer reduction_fn, const InstanceT &instance) {
        SetCoverCon<decltype(reduction_fn.reduce(instance))>;
    };

    template <typename SolverFn, typename ReducedT>
    concept SetCoverSolveFunction = requires(SolverFn solver_fn, ReducedT &&reduced, const SolverConfig &config) {
        std::invoke(solver_fn, std::forward<ReducedT>(reduced), config);
    };

    template <typename InstanceT, typename ConnAugSCReducer, typename SetCoverSolver>
        requires ConnectivityReductionFunction<ConnAugSCReducer, InstanceT> &&
                 SetCoverSolverCon<SetCoverSolver, decltype(std::declval<ConnAugSCReducer>().reduce(std::declval<InstanceT>()))>
    auto solve(const InstanceT &instance,
               ConnAugSCReducer &&reducer,
               SetCoverSolver &&solver,
               const SolverConfig &config = {})
    {
        auto reduced_problem = reducer.reduce(instance);
        return solver.solve(reduced_problem, config);
    }
}