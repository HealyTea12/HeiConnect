#pragma once

#include <functional>
#include <ostream>

class SetCoverWriter
{
public:
    template <typename SetCoverType>
    static void write(const SetCoverType &set_cover, std::ostream &out)
    {
        write_set_cover_pace_challenge(set_cover, out);
    }

private:
    // Pace challenge wants a hitting set problem equivalent to our set cover
    // so we have to write the inverse of our set cover with 1 based indexing
    template <typename SetCoverType>
    static void write_set_cover_pace_challenge(const SetCoverType &set_cover,
                                               std::ostream &out)
    {
        out << "p " << "hs " << set_cover.get_num_sets() << " " << set_cover.get_num_elements() << "\n";
        if constexpr (requires { set_cover.forEachSet(size_t{}, std::function<void(size_t)>{}); })
        {
            for (size_t e{}; e < set_cover.get_num_elements(); e++)
            {
                set_cover.forEachSet(e, [&out](size_t set_index)
                                     { out << set_index + 1 << " "; });
                out << "\n";
            }
        }
        else
        {
            for (size_t e{}; e < set_cover.get_num_elements(); e++)
            {
                for (size_t set_index{}; set_index < set_cover.get_num_sets(); set_index++)
                {
                    set_cover.forEachElement(set_index, [&out, e, set_index](size_t element)
                                             {
                        if (element == e)
                        {
                            out << set_index + 1 << " ";
                        } });
                }
                out << "\n";
            }
        }
    }
};