include(FetchContent)
FetchContent_Declare(
        googletest
        URL https://github.com/google/googletest/archive/03597a01ee50ed33e9dfd640b249b4be3799d395.zip
        URL_HASH SHA256=edd885a1ab32b6999515a880f669efadb80b3f880215f315985fa3f6eca7c4d3
        DOWNLOAD_EXTRACT_TIMESTAMP TRUE

)
FetchContent_MakeAvailable(googletest)

enable_testing()

add_executable(
        MiniCFD_tests

        lingalg_basics_tests.cpp
        solver_tests.cpp
        io_basic_test.cpp
)

target_include_directories(
        MiniCFD_tests

        PUBLIC ${PROJECT_SOURCE_DIR}/src
)

target_compile_definitions(
        MiniCFD_tests
        PUBLIC
        TEST_RESOURCE_DIR="${PROJECT_SOURCE_DIR}/tests/resources"
)

target_link_libraries(
        MiniCFD_tests
        GTest::gtest_main
)

include(GoogleTest)
gtest_discover_tests(MiniCFD_tests)
