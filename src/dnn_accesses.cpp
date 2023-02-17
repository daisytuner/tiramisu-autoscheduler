#include "../include/dnn_accesses.h"

#include <iostream>

namespace daisy {
    
namespace tiramisu
{

std::vector<dnn_iterator> 
dnn_iterator::get_iterators_from_computation(computation const& comp)
{
    std::vector<dnn_iterator> iters_list;
    
    isl_set *iter_domain = comp.get_iteration_domain();
    int nb_iterators = isl_set_dim(iter_domain, isl_dim_set);    
    for (int i = 0; i < nb_iterators; ++i)
    {
        std::string name = isl_set_get_dim_name(iter_domain, isl_dim_set, i);
        int low_bound = utility::get_bound(iter_domain, i, false);
        int up_bound = utility::get_bound(iter_domain, i, true);

        iters_list.push_back(dnn_iterator(name, low_bound, up_bound));
    }
    
    return iters_list;
}

void dnn_access_matrix::transform_matrix_by_skewing(int first_node_depth,int alpha,int beta,int gamma,int sigma)
{
    /**
     * We change bases as it is done in linear algebra
    */
    //we search for 

    int first_target = -1;
    int second_target = -1;

    for(int i=0;i<matrix.size();i++)
    { // for all the buffer dimensions
        if(matrix[i][first_node_depth] == 1)
        {
            if(first_target == -1)
            {
                first_target = i;
            }
            else
            {
                first_target = -2; //disable change
            }
        }
        if(matrix[i][first_node_depth+1] == 1)
        {
            if(second_target == -1)
            {
                second_target = i;
            }
            else
            {
                second_target = -2; //disable change
            }
        }   
    }

    if((first_target >=0) && (second_target >= 0))
    {// both have lines general case access change
        int cst_1 = matrix[first_target][matrix[first_target].size()-1];

        int cst_2 = matrix[second_target][matrix[second_target].size()-1];

        //transform
        int cst_new_1 = cst_1*alpha+cst_2*beta;
        int cst_new_2 = cst_1*gamma+cst_2*sigma;

        matrix[first_target][matrix[first_target].size()-1] = cst_new_1;
        matrix[second_target][matrix[second_target].size()-1] = cst_new_2;

    }
    else
    {
        if(first_target >=0)
        {// special case where only one is concrete buffer mapping (first)
            matrix[first_target][matrix[first_target].size()-1] *=alpha ;
        }

        if(second_target >=0)
        {// special case where only one is concrete buffer mapping (second)
            matrix[second_target][matrix[second_target].size()-1] *=sigma ;
        }

    }


}

void dnn_access_matrix::print_access_matrix() const
{
    std::cout<<buffer_name<<":";
    for(auto& line:this->matrix)
    {

        for(int val:line)
        {
            std::cout<<val<<" ";
        }
        std::cout<<",";
        
    }
    std::cout<<"\n";
}

void dnn_accesses::modify_accesses_by_skewing(int first_node_depth,int alpha,int beta,int gamma,int sigma)
{

    for(auto& access:this->accesses_list)
    {
        access.transform_matrix_by_skewing(first_node_depth,alpha,beta,gamma,sigma);
    }
}

void dnn_accesses::print_all_access() const
{
    for(auto& matrix:this->accesses_list)
    {
        matrix.print_access_matrix();
    }
}

}

}
