#pragma once

#include <vector>
#include <functional>

template<typename Context>
concept WeightedContextRequirements = requires(Context context) {
    typename Context::SetCoverType;
    typename Context::SetCoverType::ElementID;
    { context.forEachUncoveredElement(std::function<void(typename Context::SetCoverType::ElementID)>{}) };
};

template<typename Context>
    requires WeightedContextRequirements<Context>
class WeightedElementsContext : public Context
{
public:
    using SetCoverType = typename Context::SetCoverType;
    using ElementID = typename SetCoverType::ElementID;

    explicit WeightedElementsContext(const SetCoverType& set_cover) :
        Context(set_cover),
        m_elementWeights(set_cover.get_num_elements(), 1.0)
    {}

    void add_set(size_t set_index)
    {
        Context::add_set(set_index);
        Context::forEachUncoveredElement([this](ElementID element) {
            m_elementWeights[element] *= 0.5; // Example: halve the weight of uncovered elements when covered
        });
    }

    void remove_set(size_t set_index)
    {
        Context::remove_set(set_index);
        Context::forEachUncoveredElement([this](ElementID element) {
            m_elementWeights[element] *= 2.0; // Example: double the weight of uncovered elements when uncovered again
        });
    }

private:
    std::vector<double> m_elementWeights;
};