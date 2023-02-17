#include "../include/search_method.h"

#include <random>


namespace daisy {

namespace tiramisu
{

std::vector<std::pair<std::string, float>> mcts::search(sdfg_wrapper& sdfg)
{
    this->nb_explored_schedules = 0;

    std::default_random_engine rand_generator;
    std::vector<std::shared_ptr<syntax_tree>> Q;
    
    std::shared_ptr<syntax_tree> base_ast = sdfg.get_base_ast().copy_ast();
    base_ast->evaluation = model_eval->evaluate(sdfg, *base_ast);
    base_ast->search_depth = 0;
    this->nb_explored_schedules += 1;
    Q.push_back(base_ast);

    base_ast->print_ast();
    std::cout << base_ast->evaluation << std::endl;

    for (size_t epoch = 0; epoch < this->beam_size; ++epoch)
    {
        std::shared_ptr<syntax_tree> sample = base_ast->copy_ast();
        sample->search_depth = 0;
        for (size_t depth = 0; depth < max_depth; ++depth)
        {
            optimization_type optim_type = DEFAULT_OPTIMIZATIONS_ORDER[depth % DEFAULT_OPTIMIZATIONS_ORDER.size()];
            auto children_ = scheds_gen->generate_schedules(*sample, optim_type);

            std::vector<float> weights;
            std::vector<std::shared_ptr<syntax_tree>> children;
            for (auto& child : children_)
            {
                child->transform_ast();
                if (!child->ast_is_legal()) {
                    continue;
                }

                child->print_ast();

                child->search_depth = depth + 1;
                child->evaluation = model_eval->evaluate(sdfg, *child);

                std::cout << child->evaluation << std::endl;

                weights.push_back(child->evaluation * -1);
                children.push_back(child);

                this->nb_explored_schedules++;
            }

            children.push_back(sample->copy_ast());
            weights.push_back(sample->evaluation * -1.0);

            // Sample an AST
            std::discrete_distribution<int> dist(weights.begin(), weights.end());
            sample = children[dist(rand_generator)];
            Q.push_back(sample);
        }
    }

    std::sort(Q.begin(), Q.end(), [](auto a, auto b) {
        return a->evaluation < b->evaluation;
    });
    Q.resize(std::min(beam_size, Q.size()));

    std::vector<std::pair<std::string, float>> res;
    for (auto& candidate : Q)
    {
        res.push_back({candidate->get_schedule_str(), candidate->evaluation});
    }
    return res;
}

std::vector<std::pair<std::string, float>> beam_search::search(sdfg_wrapper& sdfg)
{
    this->nb_explored_schedules = 0;

    std::vector<std::shared_ptr<syntax_tree>> Q;
    std::vector<std::shared_ptr<syntax_tree>> Q_;

    // Evaluate base AST
    std::shared_ptr<syntax_tree> base_ast = sdfg.get_base_ast().copy_ast();
    base_ast->evaluation = model_eval->evaluate(sdfg, *base_ast);
    base_ast->search_depth = 0;
    Q.push_back(base_ast);

    base_ast->print_ast();
    std::cout << base_ast->evaluation << std::endl;

    this->nb_explored_schedules++;

    // Optimize
    for (size_t depth = 0; depth < max_depth; depth++)
    {
        Q_.clear();
        for (auto& ast : Q)
        {
            if (ast->search_depth < (int) depth) {
                continue;
            }

            for (size_t optim = 0; optim < DEFAULT_OPTIMIZATIONS_ORDER.size(); optim++)
            {
                optimization_type optim_type = DEFAULT_OPTIMIZATIONS_ORDER[optim];
                auto candidates = scheds_gen->generate_schedules(*ast, optim_type);
                for (auto candidate : candidates)
                {
                    candidate->transform_ast();
                    if (!candidate->ast_is_legal()) {
                        continue;
                    }

                    candidate->print_ast();

                    candidate->evaluation = model_eval->evaluate(sdfg, *candidate);
                    candidate->search_depth = depth + 1;        

                    std::cout << candidate->evaluation << std::endl;

                    this->nb_explored_schedules++;
                    Q_.push_back(candidate);
                }
            }
        }
        Q.insert(Q.end(), Q_.begin(), Q_.end());

        // Resize to beam
        std::sort(Q.begin(), Q.end(), [](auto a, auto b) {
            return a->evaluation < b->evaluation;
        });
        Q.resize(std::min(beam_size, Q.size()));

        if (Q_.empty()) {
            break;
        }
    }

    std::sort(Q.begin(), Q.end(), [](auto a, auto b) {
        return a->evaluation < b->evaluation;
    });
    Q.resize(std::min(beam_size, Q.size()));

    std::vector<std::pair<std::string, float>> res;
    for (auto& candidate : Q)
    {
        res.push_back({candidate->get_schedule_str(), candidate->evaluation});
    }
    return res;
}

}

}
