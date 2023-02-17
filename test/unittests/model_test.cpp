#pragma once

#include <utility>
#include <memory>
#include <vector>
#include <iostream>

#include <nlohmann/json.hpp>

#include <function.h>
#include <computation.h>
#include <ast.h>
#include <evaluator.h>

using namespace daisy;

using namespace 
daisy::tiramisu;

using json = nlohmann::json;

TEST(TestSuiteAST, TestStream)
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

    // Run logic
    fct->perform_full_dependency_analysis();

    // Construct AST
    std::shared_ptr<syntax_tree> ast(new syntax_tree(fct));
    ASSERT_TRUE(ast->ast_is_legal());

    // Generate json
    const std::string py_cmd_path = "/home/lukas/anaconda3/envs/daisy/bin/python";
    const std::string py_interface_path = "/home/lukas/repos/daisy/backends/tiramisu/model/main.py";
    std::shared_ptr<evaluate_by_learning_model> model(
        new evaluate_by_learning_model(py_cmd_path, {py_interface_path})
    );

    json program_repr = json::parse(model->get_program_json(*ast));

    // Check computations
    ASSERT_EQ(program_repr["computations"].size(), 1);
    ASSERT_EQ(program_repr["computations"]["S0"]["real_dimensions"].size(), 2);
    ASSERT_EQ(program_repr["computations"]["S0"]["write_access_relation"], "{ S0[i, j] -> C[i, j] }");
    ASSERT_FALSE(program_repr["computations"]["S0"]["comp_is_reduction"]);
    
    ASSERT_EQ(program_repr["computations"]["S0"]["number_of_additions"], 1);
    ASSERT_EQ(program_repr["computations"]["S0"]["number_of_subtraction"], 0);
    ASSERT_EQ(program_repr["computations"]["S0"]["number_of_multiplication"], 1);
    ASSERT_EQ(program_repr["computations"]["S0"]["number_of_division"], 0);

    ASSERT_EQ(program_repr["computations"]["S0"]["data_type"], "float64");
    ASSERT_EQ(program_repr["computations"]["S0"]["data_type_size"], 8);

    ASSERT_EQ(program_repr["computations"]["S0"]["accesses"].size(), 1);
    ASSERT_NE(program_repr["computations"]["S0"]["accesses"][0]["buffer_id"], -1);
    ASSERT_FALSE(program_repr["computations"]["S0"]["accesses"][0]["access_is_reduction"]);

    json schedule_repr = json::parse(model->get_schedule_json(*ast));
};
