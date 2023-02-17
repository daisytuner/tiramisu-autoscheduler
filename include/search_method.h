#pragma once

#include <climits>
#include <cfloat>
#include <vector>
#include <iostream>
#include <memory>

#include "schedules_generator.h"
#include "evaluator.h"
#include "sdfg_wrapper.h"
#include "utils.h"

namespace daisy { 
    
namespace tiramisu
{

class auto_scheduler;

const std::vector<optimization_type> DEFAULT_OPTIMIZATIONS_ORDER = {
    INTERCHANGE,
    TILING,
    PARALLELIZE
};

class search_method
{
friend auto_scheduler;
    
protected:
    std::shared_ptr<schedules_generator> scheds_gen;

    std::shared_ptr<evaluate_by_learning_model> model_eval;

    size_t max_depth;
    
    size_t beam_size;
        
public:
    search_method(
        std::shared_ptr<schedules_generator> scheds_gen,
        std::shared_ptr<evaluate_by_learning_model> model_eval,
        size_t max_depth,
        size_t beam_size
    )
    : scheds_gen(scheds_gen), model_eval(model_eval), max_depth(max_depth), beam_size(beam_size)
    {}

    size_t nb_explored_schedules;

    virtual ~search_method() {}
        
    virtual std::vector<std::pair<std::string, float>> search(sdfg_wrapper& sdfg) = 0;

};

class mcts : public search_method
{
public:

    mcts(
        std::shared_ptr<schedules_generator> scheds_gen,
        std::shared_ptr<evaluate_by_learning_model> model_eval,
        size_t max_depth,
        size_t beam_size
    )
    : search_method(scheds_gen, model_eval, max_depth, beam_size)
    {}
        
    virtual ~mcts() {}
    
    virtual std::vector<std::pair<std::string, float>> search(sdfg_wrapper& sdfg);

};

class beam_search : public search_method
{
private:
    std::vector<std::shared_ptr<syntax_tree>> beam_search_subroutine(sdfg_wrapper& sdfg);

public:

    beam_search(
        std::shared_ptr<schedules_generator> scheds_gen,
        std::shared_ptr<evaluate_by_learning_model> model_eval,
        size_t max_depth,
        size_t beam_size
    )
    : search_method(scheds_gen, model_eval, max_depth, beam_size)

    {} 
        
    virtual ~beam_search() {}

    virtual std::vector<std::pair<std::string, float>> search(sdfg_wrapper& sdfg);
};

}

}
