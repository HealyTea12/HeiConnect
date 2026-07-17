#include "algorithm_registry.hpp"

#include <algorithm>
#include <utility>

void AlgorithmRegistry::add(AlgorithmEntry entry)
{
    auto name = entry.name;
    auto [it, inserted] = m_entries.emplace(std::move(name), std::move(entry));
    if (!inserted)
    {
        throw std::invalid_argument("Algorithm already registered: " + it->first);
    }
}

void AlgorithmRegistry::add(std::string name, std::string description, AlgorithmFactory factory)
{
    add(AlgorithmEntry{std::move(name), std::move(description), std::move(factory)});
}

std::unique_ptr<AlgorithmRunner> AlgorithmRegistry::create(const std::string_view name, const ParamMap& params) const
{
    auto it = m_entries.find(std::string(name));
    if (it == m_entries.end())
    {
        throw std::invalid_argument("Unknown algorithm: " + std::string(name));
    }
    return it->second.factory(params);
}

std::vector<std::string> AlgorithmRegistry::names() const
{
    std::vector<std::string> result;
    result.reserve(m_entries.size());
    for (const auto& [name, _] : m_entries)
    {
        result.push_back(name);
    }
    std::sort(result.begin(), result.end());
    return result;
}

AlgorithmRegistry& get_global_registry()
{
    static AlgorithmRegistry registry{};
    return registry;
}
