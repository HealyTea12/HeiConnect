#pragma once

class RuntimePipeline
{
public:
    virtual ~RuntimePipeline() = default;
    template<typename... Args>
    void run(Args&&... args)
    {}

private:
};