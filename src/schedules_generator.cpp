#include "../include/schedules_generator.h"

#include <tuple>
#include <iostream>

#include "../include/evaluator.h"
#include "../include/function.h"

namespace daisy
{

namespace tiramisu
{

std::vector<std::shared_ptr<syntax_tree>> daisy_generator::generate_schedules(syntax_tree const& ast, optimization_type optim)
{
    std::vector<std::shared_ptr<syntax_tree>> states;
    
    switch(optim)
    {

        case optimization_type::TILING:
            for (ast_node *root : ast.roots)
                generate_tilings(root, states, ast);
            
            break;

        case optimization_type::INTERCHANGE:
            for (ast_node *root : ast.roots)
                generate_interchanges(root, states, ast);
                    
            break;

        case optimization_type::PARALLELIZE:
            for (ast_node *root : ast.roots)
                generate_parallelizations(root, states, ast);
                    
            break;

        default:
            break;
    }
    
    return states;
}

void daisy_generator::generate_tilings(ast_node *node, std::vector<std::shared_ptr<syntax_tree>>& states, syntax_tree const& ast)
{
    int branch_depth = node->get_loop_levels_chain_depth();

    // Generate tiling with dimension 2, but don't tile the tiles
    if (
        !endsWith(node->name, "_inner") && !endsWith(node->name, "_outer") && node->depth + 1 < branch_depth && branch_depth <= 5
    )
    {
        for (int tiling_size1 : tiling_factors_list)
        {
            if (!can_split_iterator(node->get_extent(), tiling_size1)) {
                continue;
            }                

            ast_node *node2 = node->children[0];
            if (endsWith(node2->name, "_inner") || endsWith(node2->name, "_outer")) {
                break;
            }
            for (int tiling_size2 : tiling_factors_list)
            {
                if (!can_split_iterator(node2->get_extent(), tiling_size2)) {
                    continue;
                }
                    
                // Copy the AST, and add tiling to the list of optimizations
                std::shared_ptr<syntax_tree> new_ast(
                    new syntax_tree()
                );
                ast_node *new_node = ast.copy_and_return_node(*new_ast, node);
                    
                optimization_info optim_info;
                optim_info.type = optimization_type::TILING;
                optim_info.node = new_node;
                
                optim_info.nb_l = 2;
                optim_info.l0 = node->depth;
                optim_info.l1 = node->depth + 1;
                
                optim_info.l0_fact = tiling_size1;
                optim_info.l1_fact = tiling_size2;
                
                new_node->get_all_computations(optim_info.comps);
                
                new_ast->new_optims.push_back(optim_info);
                states.push_back(new_ast);
                
                // Generate tiling with dimension 3
                if (node->depth + 2 < branch_depth && branch_depth <= 4)
                {
                    ast_node *node3 = node2->children[0];
                    if (endsWith(node3->name, "_inner") || endsWith(node3->name, "_outer")) {
                        break;
                    }

                    for (int tiling_size3 : tiling_factors_list)
                    {
                        if (!can_split_iterator(node3->get_extent(), tiling_size3))
                            continue;
                            
                        // Copy the AST, and add tiling to the list of optimizations
                        std::shared_ptr<syntax_tree> new_ast(
                            new syntax_tree()
                        );
                        ast_node *new_node = ast.copy_and_return_node(*new_ast, node);
                            
                        optimization_info optim_info;
                        optim_info.type = optimization_type::TILING;
                        optim_info.node = new_node;
                        
                        optim_info.nb_l = 3;
                        optim_info.l0 = node->depth;
                        optim_info.l1 = node->depth + 1;
                        optim_info.l2 = node->depth + 2;
                        
                        optim_info.l0_fact = tiling_size1;
                        optim_info.l1_fact = tiling_size2;
                        optim_info.l2_fact = tiling_size3;
                        
                        new_node->get_all_computations(optim_info.comps);
                        
                        new_ast->new_optims.push_back(optim_info);
                        states.push_back(new_ast);
                    }
                }
            }
        }
    }
    
    for (ast_node *child : node->children)
        generate_tilings(child, states, ast);
}

void daisy_generator::generate_interchanges(ast_node *node, std::vector<std::shared_ptr<syntax_tree>>& states, syntax_tree const& ast)
{
    if (node->get_extent() > 1)
    {
        int branch_depth = node->get_loop_levels_chain_depth();
        
        std::vector<ast_node*> lower_levels;
        node->get_shared_nodes_from_outermost(lower_levels);

        for (int i = node->depth + 1; i < branch_depth; ++i)
        {
            auto lower_node = lower_levels[i - node->depth - 1];
            if (endsWith(lower_node->name, "_inner")) {
                break;
            }

            // Copy the AST, and add interchange to the list of optimizations
            std::shared_ptr<syntax_tree> new_ast(
                new syntax_tree()
            );
            ast_node *new_node = ast.copy_and_return_node(*new_ast, node);
            
            optimization_info optim_info;
            optim_info.type = optimization_type::INTERCHANGE;
            optim_info.node = new_node;
                
            optim_info.nb_l = 2;
            optim_info.l0 = node->depth;
            optim_info.l1 = i;
            new_node->get_all_computations(optim_info.comps);
                
            new_ast->new_optims.push_back(optim_info);
            states.push_back(new_ast);
        }
    }
    
    for (ast_node *child : node->children)
        generate_interchanges(child, states, ast);
}

void daisy_generator::generate_parallelizations(ast_node *node, std::vector<std::shared_ptr<syntax_tree>>& states, syntax_tree const& ast)
{
    if (node->get_extent() > 1 && !node->parallelized && !endsWith(node->name, "_inner"))
    {
        std::shared_ptr<syntax_tree> new_ast(
            new syntax_tree()
        );
        ast_node *new_node = ast.copy_and_return_node(*new_ast, node);

        optimization_info optim_info;
        optim_info.type = optimization_type::PARALLELIZE;
        optim_info.node = new_node;

        optim_info.nb_l = 1;        
        optim_info.l0 = new_node->depth;
        optim_info.l0_fact = 0;
        new_node->get_all_computations(optim_info.comps);

        new_ast->new_optims.push_back(optim_info);
        states.push_back(new_ast);
    }
    
    for (ast_node *child : node->children)
        generate_parallelizations(child, states, ast);
}

}
}