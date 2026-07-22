include(FetchContent)
FetchContent_Declare(
  googletest
  URL https://github.com/google/googletest/archive/03597a01ee50ed33e9dfd640b249b4be3799d395.zip
)
# For Windows: Prevent overriding the parent project's compiler/linker settings
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)

enable_testing()
set(test-source-pattern "test_*.cpp")
file(GLOB test_sources "${test-source-pattern}")
message(STATUS "Found test sources: ${test_sources}")

add_executable(
    heiconnect_tests
    test_set_cover.cpp
    test_graph_augmentation.cpp
    ${PROJECT_SOURCE_DIR}/extern/VieCut/extlib/tlx/tlx/logger.cpp
    ${test_sources}
)
target_link_libraries(
    heiconnect_tests
PRIVATE
    GTest::gtest_main
    heiconnect
)
target_include_directories(
    heiconnect_tests
PRIVATE
    ${PROJECT_SOURCE_DIR}/extern/VieCut/lib
    ${PROJECT_SOURCE_DIR}/extern/VieCut/extlib/tlx
)
enable_testing()

include(GoogleTest)
gtest_discover_tests(heiconnect_tests)
