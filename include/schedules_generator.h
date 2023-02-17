#pragma once

#include <vector>
#include <tuple>

#include "ast.h"
#include "evaluator.h"

namespace daisy {

namespace tiramisu
{

const std::vector<int> TILING_FACTORS_DEFAULT_LIST = {2, 4, 6, 8, 13, 16, 26, 32, 64, 128, 256};
/**
 * Generate a set of AST's from a given AST.
 * Inherit this class to implement a new way to generate schedules.
 */
class schedules_generator
{
private:

protected:
    /**
     * A list of tiling factors to apply when tiling is applied.
     */
    std::vector<int> tiling_factors_list;

public:
    schedules_generator(
        std::vector<int> const& tiling_factors_list = TILING_FACTORS_DEFAULT_LIST
    )
    : tiling_factors_list(tiling_factors_list)
    {}

    virtual ~schedules_generator() {}

    /**
     * Given an AST, and an optimization to apply, 
     * generate new ASTs by applying the given optimization.
     */
    virtual std::vector<std::shared_ptr<syntax_tree>> generate_schedules(syntax_tree const& ast, optimization_type optim) = 0;
};

/**
 * Generate all combinations of the following optimizations :
 * Fusion, tiling, interchange.
 */
class daisy_generator : public schedules_generator
{
private:

protected:
    /**
     * Try to apply tiling such as the given node is the first loop to tile, 
     * and then call this method recursively on children of the given node.
     */
    void generate_tilings(ast_node *node, std::vector<std::shared_ptr<syntax_tree>>& states, syntax_tree const& ast);
    
    /**
     * Try to apply interchange by swapping the given node with one of its descendents, 
     * and then call this method recursively on children of the given node.
     */
    void generate_interchanges(ast_node *node, std::vector<std::shared_ptr<syntax_tree>>& states, syntax_tree const& ast);
    
    void generate_parallelizations(ast_node *node, std::vector<std::shared_ptr<syntax_tree>>& states, syntax_tree const& ast);

public:
    daisy_generator(
        std::vector<int> const& tiling_factors_list = TILING_FACTORS_DEFAULT_LIST
    )
    : schedules_generator(tiling_factors_list) {}

    virtual std::vector<std::shared_ptr<syntax_tree>> generate_schedules(syntax_tree const& ast, optimization_type optim);
};


} }
