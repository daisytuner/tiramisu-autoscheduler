#include "../include/sdfg_wrapper.h"

#include <pybind11/embed.h>
#include <nlohmann/json.hpp>

#include "../include/function.h"

#include <iostream>

using json = nlohmann::json;

namespace daisy {

namespace tiramisu {

namespace py = pybind11;

sdfg_wrapper::sdfg_wrapper(
    std::string sdfg_path,
    std::string build_folder,
    std::shared_ptr<function> fct
)
: name(fct->get_name()), sdfg_path(sdfg_path), build_folder(build_folder), fct(fct), base_ast(new syntax_tree(this->fct))
{

}

std::string sdfg_wrapper::get_name() const
{
    return this->name;
}

function& sdfg_wrapper::get_function()
{
    return *this->fct;
}

syntax_tree& sdfg_wrapper::get_base_ast()
{
    return *this->base_ast;
}

std::vector<float> sdfg_wrapper::benchmark(syntax_tree& ast)
{
    this->fct->reset_schedules();
    apply_optimizations(ast);

    std::vector<float> measurements = {1.0, 2.0, 3.0, 4.0, 5.0};

    this->fct->reset_schedules();
    return measurements;
}


std::shared_ptr<sdfg_wrapper> sdfg_wrapper::create(std::string sdfg_path)
{
    // py::scoped_interpreter guard{};

    py::module_ daisy_backend = py::module_::import("daisytuner.optimizer.tiramisu.tiramisu_tuner");
    py::object py_result = daisy_backend.attr("encode")(sdfg_path);

    std::string repr_str = py_result.cast<std::string>();
    json repr = json::parse(repr_str);

    // Create function from json
    std::shared_ptr<function> fct(
        new function(repr["name"])
    );

    // Add buffers
    for (auto& item : repr["buffers"]) {
        std::shared_ptr<buffer> buf(
            new buffer(
                item["name"],
                item["dim_sizes"],
                str_to_tiramisu_primitive_type(item["type"]),
                str_to_tiramisu_argument_type(item["argt"])
            )
        );
        fct->add_buffer({buf->get_name(), buf});
    }

    // Add computations
    for (auto& item : repr["computations"]) {
        dace_desc desc;
        desc.nb_additions = item["flops"]["number_of_additions"];
        desc.nb_substractions = item["flops"]["number_of_subtraction"];
        desc.nb_multiplications = item["flops"]["number_of_multiplication"];
        desc.nb_divisions = item["flops"]["number_of_division"];
        
        buffer_accesses accesses;
        for (auto& access : item["buffer_accesses"]) {
            std::vector<std::vector<int>> access_matrix;
            for (auto& row : access["access_matrix"]) {
                std::vector<int> r;
                for (auto& val : row) {
                    r.push_back(val);
                }
                access_matrix.push_back(r);
            }
            accesses.insert({access["buffer"], access_matrix});
        }
        desc.accesses = accesses;

        computation* comp = new computation(
            fct.get(),
            item["iteration_domain"],
            item["access"],
            str_to_tiramisu_primitive_type(item["type"]),
            desc
        );
    }
    fct->prepare_schedules_for_legality_checks(true);
    fct->perform_full_dependency_analysis();

    return std::shared_ptr<sdfg_wrapper>(
        new sdfg_wrapper(
            repr["sdfg_path"],
            repr["build_folder"],
            fct
        )
    );
}

}

}