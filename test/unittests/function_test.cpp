#pragma once

#include <utility>
#include <memory>
#include <vector>

#include <function.h>
#include <computation.h>
#include <buffer.h>

#include <iostream>

using namespace daisy;

using namespace 
daisy::tiramisu;

TEST(TestSuiteFunction, TestStream)
{
    std::shared_ptr<function> fct(
        new function("sdfg")
    );

    /**
     * Define the computation
     * \code
     * for (i=0; i<10; i++)
     *   for (j=0; j<20; j++)
     *      S0[i, j] = 2 * B[i, j] + 1
     * \endcode
     */
    std::string iteration_domain = "{ S0[i,j]: 0<=i<10 and 0<=j<20 }";
    std::string access = "{S0[i, j]->C[i, j]}";
    
    std::shared_ptr<buffer> B(
        new buffer(
            "B",
            std::vector<int>{10, 20},
            primitive_t::p_float64,
            argument_t::a_input
        )
    );
    fct->add_buffer({"B", B});

    std::shared_ptr<buffer> C(
        new buffer(
            "C",
            std::vector<int>{10, 20},
            primitive_t::p_float64,
            argument_t::a_output
        )
    );
    fct->add_buffer({"C", C});

    // Metadata for DL model
    dace_desc desc;
    desc.nb_additions = 1;
    desc.nb_substractions = 0;
    desc.nb_multiplications = 1;
    desc.nb_divisions = 0;

    std::vector<std::vector<int>> b_am = {
        {1, 0, 0},
        {0, 1, 0}
    };
    buffer_accesses am = {
        {"B", b_am}
    };
    desc.accesses = am;

    std::shared_ptr<computation> comp(
        new computation(
            fct.get(),
            iteration_domain,
            access,
            primitive_t::p_float64,
            desc
        )
    );

    // Testing Tiramisu functions
    const std::vector<computation*> comps = fct->get_computations();
    ASSERT_EQ(comps.size(), 1);
    ASSERT_EQ(comp.get(), comps[0]);

    auto buffers = fct->get_buffers();
    ASSERT_EQ(buffers.size(), 2);

    auto buffer = buffers.at("B");
    ASSERT_EQ(buffer->get_name(), "B");
    ASSERT_EQ(buffer->get_n_dims(), 2);
    ASSERT_EQ(buffer->get_argument_type(), argument_t::a_input);

    buffer = buffers.at("C");
    ASSERT_EQ(buffer->get_name(), "C");
    ASSERT_EQ(buffer->get_n_dims(), 2);
    ASSERT_EQ(buffer->get_argument_type(), argument_t::a_output);

    ASSERT_EQ(comp->get_name(), "S0");
    ASSERT_EQ(comp->get_data_type(), primitive_t::p_float64);
    ASSERT_EQ(comp->get_loop_levels_number(), 2);

    std::vector<std::string> level_names = comp->get_loop_level_names();
    ASSERT_EQ(level_names[0], "i");
    ASSERT_EQ(level_names[1], "j");

    ASSERT_TRUE(comp->should_schedule_this_computation());
};
