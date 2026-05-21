#pragma once
#include <tuple>
#include <utility>

template<class... Stages>
class Pipeline
{
public:
    explicit Pipeline(Stages... stages) : m_stages(std::move(stages)...)
    {}

    template<class InitialState>
    auto run(InitialState state)
    {
        return run_impl<0>(std::move(state));
    }

private:
    std::tuple<Stages...> m_stages;

    template<std::size_t I, class State>
    auto run_impl(State state)
    {
        if constexpr (I == sizeof...(Stages))
        {
            return state;
        }
        else
        {
            auto next_state = std::get<I>(m_stages).run(std::move(state));
            return run_impl<I + 1>(std::move(next_state));
        }
    }
};

template<class... Stages>
Pipeline(Stages...) -> Pipeline<Stages...>;