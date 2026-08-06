#pragma once

#include <concepts>
#include <functional>
#include <ostream>

#include "HeiConnect/set_cover/common.hpp"

template<typename SetCoverType>
concept SetCoverWritable = requires(const SetCoverType& sc, typename SetCoverType::SetID set_id) {
    { sc.get_num_sets() } -> std::convertible_to<size_t>;
    { sc.get_num_elements() } -> std::convertible_to<size_t>;
    { sc.forEachElement(set_id, std::function<void(typename SetCoverType::ElementID)>{}) };
    { sc.get_set_cost(set_id) };
};
class SetCoverWriter
{
public:
    enum class Format
    {
        DEFAULT,
        HITTING_SET_PACE_CHALLENGE,
        MINIZINC,
    };

    template<typename SetCoverType>
        requires SetCoverWritable<SetCoverType>
    static void write(const SetCoverType& set_cover, std::ostream& out, Format format = Format::DEFAULT)
    {
        switch (format)
        {
            case Format::DEFAULT: write_set_cover_default(set_cover, out); break;
            case Format::HITTING_SET_PACE_CHALLENGE: write_set_cover_pace_challenge(set_cover, out); break;
            case Format::MINIZINC: write_set_cover_minizinc(set_cover, out); break;
        }
    }

private:
    template<typename SetCoverType>
    static void write_set_cover_minizinc(const SetCoverType& set_cover, std::ostream& out)
    {
        using SetCost = std::remove_cvref_t<decltype(set_cover.get_set_cost(typename SetCoverType::SetID{}))>;
        constexpr bool has_float_costs = std::floating_point<SetCost>;

        /*
        num_sets =  10;
        num_elements = 8;
        costs = [ 19, 16, 18, 13, 15, 19, 15, 17, 16, 15];

        sets = [
          {1,6},
          {2,6,8},
          {1,4,7},
          {2,3,5},
          {2,5},
          {2,3},
          {2,3,4},
          {4,5,8},
          {3,6,8},
          {1,6,7}
        ];
        */
        out << "int: num_sets = " << set_cover.get_num_sets() << ";\n";
        out << "int: num_elements = " << set_cover.get_num_elements() << ";\n\n";
        out << "array[1..num_sets] of " << (has_float_costs ? "float" : "int") << ": costs = [";
        for (size_t set_index{}; set_index < set_cover.get_num_sets(); set_index++)
        {
            if (set_index > 0)
                out << ", ";
            out << set_cover.get_set_cost(set_index);
        }
        out << "];\n\n";
        // Print out the sets
        out << "array[1..num_sets] of set of 1..num_elements: sets = [\n";
        for (size_t set_index{}; set_index < set_cover.get_num_sets(); set_index++)
        {
            if (set_index > 0)
                out << ",\n";
            out << "  {";
            bool first_element = true;
            set_cover.forEachElement(set_index, [&out, &first_element](auto element) {
                if (!first_element)
                    out << ", ";
                out << element + 1;
                first_element = false;
            });
            out << "}";
        }
        out << "\n];\n\n";
        out << "array[1..num_sets] of var 0..1: x;\n\n";
        out << "var " << (has_float_costs ? "float" : "int")
            << ": z = sum(i in 1..num_sets) (x[i] * costs[i]);\n\n";
        out << "solve minimize z;\n\n";
        out << "constraint\n";
        out << "  forall(j in 1..num_elements) (\n";
        out << "    sum(i in 1..num_sets) (x[i] * bool2int(j in sets[i])) >= 1\n";
        out << "  );\n\n";
        out << "output\n";
        out << "[\n";
        out << "  \"cost: \" ++ show(z) ++ \"\\n\" ++\n";
        out << "  \"x: \" ++ show(x) ++ \"\\n\" ++\n";
        out << "  \"sets: \" ++ show(sets) ++ \"\\n\"\n";
        out << "];\n";
    }

    template<typename SetCoverType>
    static void write_set_cover_default(const SetCoverType& set_cover, std::ostream& out)
        requires SetCoverWritable<SetCoverType>
    {
        // Header: "p sc <num_sets> <num_elements>"
        out << "p sc " << set_cover.get_num_sets() << " " << set_cover.get_num_elements() << "\n";
        // Each line represents a set: "<element1> <element2> ... <elementN> <set_weight>"
        for (size_t set_index{}; set_index < set_cover.get_num_sets(); ++set_index)
        {
            std::cout << "s" << set_index << " of " << set_cover.get_num_sets() << "\n";
            set_cover.forEachElement(set_index, [&out](auto element) { out << element << " "; });
            out << set_cover.get_set_cost(set_index) << "\n";
        }
    }

    // Pace challenge wants a hitting set problem equivalent to our set cover
    // so we have to write the inverse of our set cover with 1 based indexing
    template<typename SetCoverType>
    static void write_set_cover_pace_challenge(const SetCoverType& set_cover, std::ostream& out)
    {
        out << "p hs " << set_cover.get_num_sets() << " " << set_cover.get_num_elements() << "\n";
        if constexpr (requires { set_cover.forEachSet(size_t{}, std::function<void(size_t)>{}); })
        {
            for (size_t e{}; e < set_cover.get_num_elements(); e++)
            {
                set_cover.forEachSet(e, [&out](size_t set_index) { out << set_index + 1 << " "; });
                out << "\n";
            }
        }
        else
        {
            for (size_t e{}; e < set_cover.get_num_elements(); e++)
            {
                for (size_t set_index{}; set_index < set_cover.get_num_sets(); set_index++)
                {
                    set_cover.forEachElement(set_index, [&out, e, set_index](auto element) {
                        if (static_cast<size_t>(element) == e)
                        {
                            out << set_index + 1 << " ";
                        }
                    });
                }
                out << "\n";
            }
        }
    }
};
