#pragma once

#include <sdfg_wrapper.h>

#include <iostream>

using namespace daisy;

using namespace 
daisy::tiramisu;

TEST(TestSuiteSDFG, TestCreate)
{
    std::shared_ptr<sdfg_wrapper> sdfg = sdfg_wrapper::create(
        "/home/lukas/repos/daisy/backends/tiramisu/test/data/hdiff.sdfg"
    );

    ASSERT_EQ(sdfg->get_name(), "benchmark_hdiff_cutout");
    ASSERT_EQ(sdfg->get_base_ast().previous_optims.size(), 0);
    ASSERT_EQ(sdfg->get_base_ast().new_optims.size(), 0);
    ASSERT_EQ(sdfg->get_function().get_computations().size(), 1);
    ASSERT_EQ(sdfg->get_function().get_buffers().size(), 2);

};
