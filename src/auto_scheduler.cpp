#include "../include/auto_scheduler.h"

#include <chrono>
#include <cmath>
#include <fstream>

#include <nlohmann/json.hpp>

using json = nlohmann::json;

namespace daisy {
    
namespace tiramisu
{

auto_scheduler::auto_scheduler(
    std::shared_ptr<sdfg_wrapper> sdfg,
    std::shared_ptr<search_method> searcher
)
: sdfg(sdfg), searcher(searcher)
{

}

std::string auto_scheduler::tune()
{            
    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

    std::vector<std::pair<std::string, float>> best_scheds = searcher->search(*this->sdfg);
    
    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count();

    // Print some info about the search
    std::cout << "Explored schedules: " << searcher->nb_explored_schedules << std::endl;
    std::cout << "Search time: " << duration << " ms " << std::endl;  

    json res(best_scheds);

    std::string tmp = std::to_string(searcher->nb_explored_schedules) + ";" + res.dump();
    return tmp;
}


std::shared_ptr<auto_scheduler> auto_scheduler::create(
    std::shared_ptr<sdfg_wrapper> sdfg,
    std::string py_cmd_path,
    std::string py_interface_path,
    std::string method,
    size_t max_depth,
    size_t beam_size
)
{
    std::shared_ptr<schedules_generator> scheds_gen(
        new daisy_generator()
    );

    std::shared_ptr<evaluate_by_learning_model> model_eval(
        new evaluate_by_learning_model(py_cmd_path, {py_interface_path})
    );

    std::shared_ptr<search_method> searcher;
    if (method == "MCTS") {
        searcher = std::shared_ptr<search_method>(
            new mcts(scheds_gen, model_eval, max_depth, beam_size)
        );
    } else {
        searcher = std::shared_ptr<search_method>(
            new beam_search(scheds_gen, model_eval, max_depth, beam_size)
        );
    }

    return std::shared_ptr<auto_scheduler>(
        new auto_scheduler(sdfg, searcher)
    );
}

}
}
