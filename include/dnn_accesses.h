#pragma once

#include <string>
#include <vector>

#include "core.h"
#include "computation.h"

namespace daisy {

namespace tiramisu
{

/**
 * Contains information about an iterator.
 * Just a convenient class to simplify working with the ML model.
 */
class dnn_iterator
{
public:
    std::string name;
    int low_bound;
    int up_bound;
    
    dnn_iterator(std::string name, int low_bound, int up_bound)
        : name(name), low_bound(low_bound), up_bound(up_bound) {}
        
    /**
     * Return a list of dnn_iterators from the iterators of the given computation.
     */
    static std::vector<dnn_iterator> get_iterators_from_computation(computation const& comp);
};

/**
 * Contains the access matrix for a given access.
 */
class dnn_access_matrix
{
public:
    int nb_iterators;
    int nb_dims;

    std::vector<std::vector<int>> matrix;
    
    /**
     * The buffer that this matrix accesses.
     */
    std::string buffer_name;

    int buffer_id;
    
    /**
     * Create an empty access matrix (filled with zeros),
     * with the given number of iterators and the given number of dimensions.
     */
    dnn_access_matrix(
        int nb_iterators,
        int nb_dims
    )
    : nb_iterators(nb_iterators), nb_dims(nb_dims), matrix()
    {
        for (int i = 0; i < nb_dims; ++i)
            matrix.push_back(std::vector<int>(nb_iterators + 1, 0));
    };

    /**
     * transforms the matrix by skewing
    */
    void transform_matrix_by_skewing(int first_node_depth,int alpha,int beta,int gamma,int sigma);

    void print_access_matrix() const;

};

/**
 * Contains the list of access matrices for a given computation.
 */
class dnn_accesses
{
public:
    /**
     * The computation from which the accesses have been retrieved.
     */
    computation* comp;
    
    int nb_iterators;
    
    /**
     * A list of matrices, such as each matrix represents an access of "comp".
     */
    std::vector<dnn_access_matrix> accesses_list;
        
    /**
     * Create the list of accesses of the given computation.
     */
    dnn_accesses(
        computation *comp,
        int nb_iterators
    )
    : comp(comp), nb_iterators(nb_iterators)
    {
        for (const auto& access : comp->desc.accesses)
        {
            int nb_d = access.second.size();

            dnn_access_matrix matrix(nb_iterators, nb_d);
            matrix.buffer_name = access.first;
            matrix.matrix = std::vector<std::vector<int>>(access.second);

            accesses_list.push_back(matrix);
        }
    };

    void modify_accesses_by_skewing(int first_node_depth,int alpha,int beta,int gamma,int sigma);

    void print_all_access() const;

};

}
}


