#pragma once

#include <utility>
#include <memory>
#include <vector>

#include <sdfg_wrapper.h>
#include <auto_scheduler.h>


using namespace daisy;
using namespace daisy::tiramisu;

TEST(TestSuiteAutoScheduler, TestCopyKernel)
{
    // Path to python (please give absolute path)
    const std::string py_cmd_path = "/home/lukas/anaconda3/envs/daisy/bin/python";
    const std::string py_interface_path = "/home/lukas/repos/daisy/backends/tiramisu/model/main.py";

    const size_t beam_size = 10;
    const size_t max_depth = 10;

    auto sdfg = sdfg_wrapper::create("/home/lukas/repos/daisy/backends/tiramisu/test/data/copy_kernel.sdfg");

    auto scheduler = auto_scheduler::create(
        sdfg,
        py_cmd_path,
        py_interface_path, 
        "beam",
        max_depth,
        beam_size
    );

    auto result = scheduler->tune();
    std::cout << result << std::endl;
};

// TEST(TestSuiteAutoScheduler, TestMatmulBEAM)
// {
//     // Path to python (please give absolute path)
//     const std::string py_cmd_path = "/home/lukas/anaconda3/envs/daisy/bin/python";
//     const std::string py_interface_path = "/home/lukas/repos/daisy/backends/tiramisu/model/main.py";

//     const size_t max_depth = 15;
//     const size_t beam_size = 10;

//     auto sdfg = sdfg_wrapper::create("/home/lukas/repos/daisy/backends/tiramisu/test/data/matmul.sdfg");

//     auto scheduler = auto_scheduler::create(
//         sdfg,
//         py_cmd_path,
//         py_interface_path, 
//         "beam",
//         max_depth,
//         beam_size
//     );

//     auto result = scheduler->tune();
// };

TEST(TestSuiteAutoScheduler, TestMatmulMCTS)
{
    // Path to python (please give absolute path)
    const std::string py_cmd_path = "/home/lukas/anaconda3/envs/daisy/bin/python";
    const std::string py_interface_path = "/home/lukas/repos/daisy/backends/tiramisu/model/main.py";

    const size_t max_depth = 20;
    const size_t beam_size = 20;


    auto sdfg = sdfg_wrapper::create("/home/lukas/repos/daisy/backends/tiramisu/test/data/matmul.sdfg");

    auto scheduler = auto_scheduler::create(
        sdfg,
        py_cmd_path,
        py_interface_path, 
        "MCTS",
        max_depth,
        beam_size
    );

    auto result = scheduler->tune();
};

TEST(TestSuiteAutoScheduler, TestHdiff)
{
    // Path to python (please give absolute path)
    const std::string py_cmd_path = "/home/lukas/anaconda3/envs/daisy/bin/python";
    const std::string py_interface_path = "/home/lukas/repos/daisy/backends/tiramisu/model/main.py";

    const size_t beam_size = 15;
    const size_t max_depth = 10;

    auto sdfg = sdfg_wrapper::create("/home/lukas/repos/daisy/backends/tiramisu/test/data/hdiff.sdfg");

    auto scheduler = auto_scheduler::create(
        sdfg,
        py_cmd_path,
        py_interface_path, 
        "beam",
        max_depth,
        beam_size
    );

    auto result = scheduler->tune();
};
