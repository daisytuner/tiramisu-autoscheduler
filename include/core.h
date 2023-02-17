#pragma once

#include <string>
#include <vector>
#include <map>
#include <unordered_map>

#include <isl/ctx.h>
#include <isl/aff.h>
#include <isl/set.h>
#include <isl/map.h>
#include <isl/space.h>
#include <isl/flow.h>
#include <isl/id.h>
#include <isl/constraint.h>
#include <isl/union_map.h>
#include <isl/union_set.h>
#include <isl/ast_build.h>
#include <isl/ilp.h>

namespace daisy {


bool endsWith(const std::string& str, const std::string& suffix);

std::string generate_new_variable_name();
std::string generate_new_computation_name();

struct param_pack_1
{
    int in_dim;
    int out_constant;
};

/**
  * tiramisu data types.
  * "p_" stands for primitive.
  */
enum primitive_t
{
    p_uint8,
    p_uint16,
    p_uint32,
    p_uint64,
    p_int8,
    p_int16,
    p_int32,
    p_int64,
    p_float32,
    p_float64,
    p_boolean,
    p_async,
    p_wait_ptr,
    p_void_ptr,  // Used for raw buffers in cuda_ast
    p_none
};

std::string str_from_tiramisu_type_primitive(primitive_t type);
primitive_t str_to_tiramisu_primitive_type(std::string type);

/**
  * Types of function arguments.
  * "a_" stands for argument.
  */
enum argument_t
{
    a_input,
    a_output,
    a_temporary
};

std::string str_from_tiramisu_type_argument(daisy::argument_t type);
argument_t str_to_tiramisu_argument_type(std::string argt);

/**
  * Add a dimension to the range of a map in the specified position.
  * Assume that the name of the new dimension is equal to the name of the corresponding
  * dimension in the domain of the map.
  * Add a constraint that indicates that the added dim is equal to a constant.
  */
isl_map* isl_map_add_dim_and_eq_constraint(isl_map *map, int dim_pos, int constant);

isl_map* isl_map_align_range_dims(isl_map *map, int max_dim);

int isl_map_get_static_dim(isl_map *map, int dim_pos);

int dynamic_dimension_into_loop_level(int dim);

isl_stat extract_static_dim_value_from_bmap(__isl_take isl_basic_map *bmap, void *user);

/**
 * A class containing utility functions.
 */
class utility
{
public:
    /**
     * Traverse recursively the ISL AST tree
     *
     * \p node represents the root of the tree to be traversed.
     *
     * \p dim is the dimension of the loop from which the bounds have to be
     * extracted.
     *
     * \p upper is a boolean that should be set to true to extract
     * the upper bound and false to extract the lower bound.
     */
     static isl_ast_expr* extract_bound_expression(isl_ast_node *ast, int dim, bool upper);

    /**
     * Return a daisy::expr representing the bound of
     * the dimension \p dim in \p set.  If \p upper is true
     * then this function returns the upper bound otherwise
     * it returns the lower bound.
     *
     * For example, assuming that
     *
     * S = {S[i,j]: 0<=i<N and 0<=j<N and i<M}
     *
     * then
     *
     * get_upper_bound(S, 1)
     *
     * would return N-1, while
     *
     * get_upper_bound(S, 0)
     *
     * would return min(N-1,M-1)
     */
    static int get_bound(isl_set *set, int dim, int upper);

    /**
     * Return the extent of the loop.
     *
     * For example:
     *
     * [N]->{C[i,j]: 0 <= i < N and N = 10}
     *
     * then get_extent(C,0) would return 10.
     *
     */
    static int get_extent(isl_set *set, int dim);

    /**
     * Create a comma separated string that represents the list
     * of the parameters of \p set.
     *
     * For example, if the set is
     *
     * [N,M,K]->{S[i]}
     *
     * this function returns the string "N,M,K".
     */
    static std::string get_parameters_list(isl_set *set);
};

}


