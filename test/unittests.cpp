#include <gtest/gtest.h>

#include <pybind11/embed.h>

// #include "unittests/function_test.cpp"
// #include "unittests/sdfg_test.cpp"
// #include "unittests/model_test.cpp"
#include "unittests/auto_schedule_test.cpp"

namespace py = pybind11;

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);

    py::scoped_interpreter guard{};

    return RUN_ALL_TESTS();
}