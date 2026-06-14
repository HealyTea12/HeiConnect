#pragma once

#include <chrono>

namespace HeiConnect::tools
{

    enum class TimeUnit
    {
        Seconds,
        Milliseconds,
        Microseconds,
        Nanoseconds
    };
    template<typename TimePoint>
    std::string format_duration(TimePoint start, TimePoint end, TimeUnit unit)
    {
        switch (unit)
        {
            case TimeUnit::Seconds:
                return std::to_string(std::chrono::duration_cast<std::chrono::seconds>(end - start).count()) + "s";
            case TimeUnit::Milliseconds:
                return std::to_string(std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count()) +
                    "ms";
            case TimeUnit::Microseconds:
                return std::to_string(std::chrono::duration_cast<std::chrono::microseconds>(end - start).count()) +
                    "us";
            case TimeUnit::Nanoseconds:
                return std::to_string(std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count()) + "ns";
        }
    }
}