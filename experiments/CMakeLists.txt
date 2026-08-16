find_package(Boost REQUIRED COMPONENTS program_options)
find_package(benchmark REQUIRED)

add_executable(generate_datasets)
target_sources(generate_datasets
PRIVATE
    generate_datasets.cpp
    dataset_generator/cli.cpp
    dataset_generator/generators.cpp)
target_link_libraries(
    generate_datasets
PRIVATE
    heiconnect
    Boost::program_options
)

add_executable(reduce_instance reduce_instance.cpp)
target_link_libraries(
    reduce_instance
PRIVATE
    heiconnect
    Boost::program_options
)

add_executable(experiments)
target_sources(experiments
PRIVATE
    exp.cpp
    algorithm_registry.cpp
    algos/sccsr.cpp
    algos/gwc.cpp
    algos/eilp.cpp
    algos/lazy_block_tree_ilp.cpp
)
target_include_directories(experiments PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(experiments
PRIVATE
    heiconnect
    Boost::program_options
)
add_executable(sc_red_bm)
target_sources(sc_red_bm
PRIVATE
    benchmarks/set_cover_reduction_BM.cpp
)
target_link_libraries(sc_red_bm
PRIVATE
    heiconnect
    Boost::program_options
)

add_executable(vec_BM)
target_sources(vec_BM
PRIVATE
    benchmarks/vec_BM.cpp
)
target_link_libraries(vec_BM
PRIVATE    
    benchmark::benchmark
    benchmark::benchmark_main
) 

add_executable(sc_red_BM)
target_sources(sc_red_BM
PRIVATE
    benchmarks/sc_red_BM.cpp
    benchmarks/cap_BM.cpp
    benchmarks/cyc_covercount_gwc_BM.cpp
)
target_link_libraries(sc_red_BM
PRIVATE    
    heiconnect
    benchmark::benchmark
    benchmark::benchmark_main
) 

add_executable(bulk_link_contraction_BM)
target_sources(bulk_link_contraction_BM
PRIVATE
    benchmarks/bulk_link_contraction_BM.cpp
)
target_link_libraries(bulk_link_contraction_BM
PRIVATE
    heiconnect
    benchmark::benchmark
    benchmark::benchmark_main
)

add_executable(cycle_intersection_index_BM)
target_sources(cycle_intersection_index_BM
PRIVATE
    benchmarks/cycle_intersection_index_BM.cpp
)
target_link_libraries(cycle_intersection_index_BM
PRIVATE
    heiconnect
    benchmark::benchmark
    benchmark::benchmark_main
)
