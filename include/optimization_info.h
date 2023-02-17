#pragma once

#include <vector>

namespace daisy {

namespace tiramisu
{

class syntax_tree;
class ast_node;
class computation;

enum optimization_type
{
    TILING,
    INTERCHANGE,
    PARALLELIZE
};

/**
 * Stores information about an optimization.
 * Check the function apply_optimizations() to see how this structure is used.
 */
struct optimization_info
{
    /**
     * The type of this optimization.
     */
    optimization_type type;
    
    /**
     * The list of computations that this optimization will be applied to.
     */
    std::vector<computation*> comps;
    
    /**
     * This attribute is used when transforming the AST.
     * It indicates the node at which to start the transformation.
     */
    ast_node *node;
    
    /**
     * The number of loop levels that this optimization affects.
     * For example, a 2 level tiling affects 2 loop levels, an interchange
     * affects 2 loop levels.
     */
    int nb_l;
    
    /**
     * The loop levels this optimization affects.
     * nb_l indicates the number of loop levels to consider.
     */
    int l0, l1, l2;
    
    /**
     * Contains the factors of each loop level.
     * For example, if the optimization is a 2 level tiling,
     * l0_fact and l1_fact will contain the tiling factors for each loop level.
     */
    int l0_fact = 0, l1_fact = 0, l2_fact = 0, l3_fact = 0;
};

/**
 * Tag the outermost level of each computation to be parallelized.
 */
void parallelize_outermost_levels(std::vector<computation*> const& comps_list);

/**
 * Apply the optimizations specified by the syntax tree using the Tiramisu API.
 */
void apply_optimizations(syntax_tree const& ast);
    
/**
 * Apply the given optimization using the Tiramisu API.
 */
void apply_optimizations(optimization_info const& optim_info);
    
/**
 * Schedule the computations so as to be in the order specified by the AST.
 */
// void apply_fusions(syntax_tree const& ast);
    
/**
 * A recursive subroutine used by apply_fusions(syntax_tree const& ast).
 */
// computation* apply_fusions(ast_node *node, computation *last_comp, int dimension);

/**
 * Apply parallelization through tiramisu API to the loop levels that correspond to the ast_nodes that are tagged for
 * parallelization in the AST
 */
void apply_parallelization(syntax_tree const& ast);

/**
 * A recursive subroutine used by apply_parallelization(syntax_tree const& ast).
 */
void apply_parallelization(ast_node *node);

/**
 * Prints the optimization information
 */
void print_optim(optimization_info optim);

}
}


