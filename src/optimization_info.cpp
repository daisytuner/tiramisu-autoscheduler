#include "../include/optimization_info.h"

#include <iostream>

#include "../include/ast.h"

namespace daisy {

namespace tiramisu
{

void apply_optimizations(syntax_tree const& ast)
{
    // Check ast.h for the difference between ast.previous_optims and ast.new_optims
    for (optimization_info const& optim_info : ast.previous_optims)
        apply_optimizations(optim_info);
        
    for (optimization_info const& optim_info : ast.new_optims)
        apply_optimizations(optim_info);

    // Parallelization needs to be applied after the other transformations in order to have the accurate loop depth of
    // the tagged ast_nodes
    apply_parallelization(ast);
}

void apply_optimizations(optimization_info const& optim_info)
{
    dace_block block(optim_info.comps);
        
    switch (optim_info.type)
    {
        case optimization_type::TILING:
            if (optim_info.nb_l == 2)
                block.tile(optim_info.l0, optim_info.l1, 
                           optim_info.l0_fact, optim_info.l1_fact);
                
            else if (optim_info.nb_l == 3)
                block.tile(optim_info.l0, optim_info.l1, optim_info.l2,
                           optim_info.l0_fact, optim_info.l1_fact, optim_info.l2_fact);
            break;
                
        case optimization_type::INTERCHANGE:
            block.interchange(optim_info.l0, optim_info.l1);
            break;

        default:
            break;
    }
}

void apply_parallelization(syntax_tree const& ast)
{
    for (ast_node *root : ast.roots)
        apply_parallelization(root);
}

void apply_parallelization(ast_node* node)
{
    // if the ast_node is tagged for parallelization, get the child computations and tag them using tag_parallel_level()
    if (node->parallelized)
    {
        std::vector<computation*> involved_computations;
        node->get_all_computations(involved_computations);
        for (computation* comp: involved_computations)
            comp->tag_parallel_level(node->depth);
    }
    for (ast_node *child : node->children)
        apply_parallelization(child);

}

void parallelize_outermost_levels(std::vector<computation*> const& comps_list)
{
    for (computation *comp : comps_list)
        comp->tag_parallel_level(0);
}

void print_optim(optimization_info optim)
{
    switch(optim.type) {
        case optimization_type::INTERCHANGE:
            std::cout << "Interchange" << " L" << optim.l0 << " " << " L" << optim.l1  << std::endl;
            break;

        case optimization_type::TILING:
            std::cout << "Tiling" << " L" << optim.l0 << " " << optim.l0_fact << " L" << optim.l1 << " " << optim.l1_fact;
            if (optim.nb_l == 3)
                std::cout << " L" << optim.l2 << " " << optim.l2_fact;
            std::cout << std::endl;
            break;

        case optimization_type::PARALLELIZE:
            std::cout << "Parallelize" << " L" << optim.l0 << std::endl;
            break;

        default:
            break;
    }
}
}

}
