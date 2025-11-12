include(FetchContent)
FetchContent_Declare(
  argtable
  GIT_REPOSITORY https://github.com/argtable/argtable3
  GIT_TAG master
)
FetchContent_MakeAvailable(argtable)

set(ARGTABLE3_ENABLE_EXAMPLES OFF)
set(ARGTABLE3_ENABLE_TESTS OFF)

add_executable(solver app/solver.cpp config.cpp)
target_link_libraries(solver PRIVATE heiconnect argtable)

add_executable(graphchecker app/graphchecker.cpp app/config.cpp)
target_compile_definitions(graphchecker PRIVATE "-DMODE_GRAPHCHECKER")
target_link_libraries(graphchecker PRIVATE OpenMP::OpenMP_CXX)
install(TARGETS graphchecker DESTINATION bin)