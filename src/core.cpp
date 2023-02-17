#include "../include/core.h"

#include <cassert>
#include <isl/val.h>
#include <isl/ast_build.h>

namespace daisy
{

bool endsWith(const std::string& str, const std::string& suffix)
{
    return str.size() >= suffix.size() && 0 == str.compare(str.size()-suffix.size(), suffix.size(), suffix);
}

int id_counter = 0;

std::string generate_new_variable_name()
{
    return "t" + std::to_string(id_counter++);
}

std::string generate_new_computation_name()
{
    return "C" + std::to_string(id_counter++);
}

std::string str_from_tiramisu_type_primitive(primitive_t type)
{
    switch (type)
    {
    case daisy::p_uint8:
        return "uint8";
    case daisy::p_int8:
        return "int8";
    case daisy::p_uint16:
        return "uint16";
    case daisy::p_int16:
        return "int16";
    case daisy::p_uint32:
        return "uint32";
    case daisy::p_int32:
        return "int32";
    case daisy::p_uint64:
        return "uint64";
    case daisy::p_int64:
        return "int64";
    case daisy::p_float32:
        return "float32";
    case daisy::p_float64:
        return "float64";
    case daisy::p_boolean:
        return "bool";
    case daisy::p_wait_ptr:
        return "wait";
    case daisy::p_void_ptr:
        return "void *";
    default:
        return "";
    }
}

primitive_t str_to_tiramisu_primitive_type(std::string type)
{
    if (type == "uint8") {
        return daisy::p_uint8;
    }
    if (type == "int8") {
        return daisy::p_int8;
    }
    if (type == "uint16") {
        return daisy::p_uint16;
    }
    if (type == "int16") {
        return daisy::p_int16;
    }
    if (type == "uint32") {
        return daisy::p_uint32;
    }
    if (type == "int32") {
        return daisy::p_int32;
    }
    if (type == "uint64") {
        return daisy::p_uint64;
    }
    if (type == "int64") {
        return daisy::p_int64;
    }
    if (type == "float32") {
        return daisy::p_float32;
    }
    if (type == "float64") {
        return daisy::p_float64;
    }
    if (type == "bool") {
        return daisy::p_boolean;
    }
    if (type == "wait") {
        return daisy::p_wait_ptr;
    }
    if (type == "void *") {
        return daisy::p_void_ptr;
    }
    return daisy::p_none;
}

std::string str_from_tiramisu_type_argument(daisy::argument_t type)
{
    switch (type)
    {
    case daisy::a_input:
        return "input";
    case daisy::a_output:
        return "output";
    case daisy::a_temporary:
        return "temporary";
    default:
        // ERROR("Tiramisu type not supported.", true);
        return "";
    }
}

argument_t str_to_tiramisu_argument_type(std::string argt)
{
    if (argt == "input") {
        return daisy::argument_t::a_input;
    } else if (argt == "output") {
        return daisy::argument_t::a_output;
    } else {
        return daisy::argument_t::a_temporary;
    }
}

isl_map* isl_map_add_dim_and_eq_constraint(isl_map *map, int dim_pos, int constant)
{
    map = isl_map_insert_dims(map, isl_dim_out, dim_pos, 1);
    map = isl_map_set_tuple_name(map, isl_dim_out, isl_map_get_tuple_name(map, isl_dim_in));

    isl_space *sp = isl_map_get_space(map);
    isl_local_space *lsp =
        isl_local_space_from_space(isl_space_copy(sp));
    isl_constraint *cst = isl_constraint_alloc_equality(lsp);
    cst = isl_constraint_set_coefficient_si(cst, isl_dim_out, dim_pos, 1);
    cst = isl_constraint_set_constant_si(cst, (-1) * constant);
    map = isl_map_add_constraint(map, cst);

    return map;
}

isl_map *isl_map_align_range_dims(isl_map *map, int max_dim)
{
    assert(map != NULL);
    int mdim = isl_map_dim(map, isl_dim_out);
    //assert(max_dim >= mdim);

    // in case where the max_dim is bigger than this map dimension, we add zeros to the schedule.
    if(max_dim >= mdim)
    {
        //DEBUG(10, daisy::str_dump("Input map:", isl_map_to_str(map)));

        const char *original_range_name = isl_map_get_tuple_name(map, isl_dim_out);

        map = isl_map_add_dims(map, isl_dim_out, max_dim - mdim);

        for (int i = mdim; i < max_dim; i++)
        {
            isl_space *sp = isl_map_get_space(map);
            isl_local_space *lsp = isl_local_space_from_space(sp);
            isl_constraint *cst = isl_constraint_alloc_equality(lsp);
            cst = isl_constraint_set_coefficient_si(cst, isl_dim_out, i, 1);
            map = isl_map_add_constraint(map, cst);
        }

        map = isl_map_set_tuple_name(map, isl_dim_out, original_range_name);

        //DEBUG(10, daisy::str_dump("After alignment, map = ",
        //                            isl_map_to_str(map)));
    }
    else
    {
      // in case where the max_dim is smaller than this map dimension, we project_out (delete) additional dimensions
       
        //DEBUG(10, daisy::str_dump("Input map:", isl_map_to_str(map)));
         map = isl_map_project_out(map,isl_dim_out,max_dim,mdim-max_dim);
        //DEBUG(10, daisy::str_dump("After alignment, map = ",
        //                            isl_map_to_str(map)));
    }

    //DEBUG_INDENT(-4);
    return map;
}

/**
 * Return the value of the static dimension.
 *
 * For example, if we have a map M = {S0[i,j]->[0,0,i,1,j,2]; S0[i,j]->[1,0,i,1,j,3]}
 * and call isl_map_get_static_dim(M, 5, 1), it will return 3.
 */
int isl_map_get_static_dim(isl_map *map, int dim_pos)
{
    // DEBUG_FCT_NAME(3);
    // DEBUG_INDENT(4);

    assert(map != NULL);
    assert(dim_pos >= 0);
    assert(dim_pos <= (signed int) isl_map_dim(map, isl_dim_out));

    // DEBUG(3, daisy::str_dump("Getting the constant coefficient of ",
    //                             isl_map_to_str(map));
    //       daisy::str_dump(" at dimension ");
    //       daisy::str_dump(std::to_string(dim_pos)));

    struct param_pack_1 *data = (struct param_pack_1 *) malloc(sizeof(struct param_pack_1));
    data->out_constant = 0;
    data->in_dim = dim_pos;

    isl_map_foreach_basic_map(isl_map_copy(map),
                              &extract_static_dim_value_from_bmap,
                              data);

    // DEBUG(3, daisy::str_dump("The constant is: ");
    //       daisy::str_dump(std::to_string(data->out_constant)));

    // DEBUG_INDENT(-4);

    return data->out_constant;
}

/**
 * Take a basic map as input, go through all of its constraints,
 * identifies the constraint of the static dimension param_pack_1.in_dim
 * (passed in user) and replace the value of param_pack_1.out_constant if
 * the static dimension is bigger than that value.
 */
isl_stat extract_static_dim_value_from_bmap(__isl_take isl_basic_map *bmap, void *user)
{
    struct param_pack_1 *data = (struct param_pack_1 *) user;

    isl_constraint_list *list = isl_basic_map_get_constraint_list(bmap);
    int n_constraints = isl_constraint_list_n_constraint(list);

    for (int i = 0; i < n_constraints; i++)
    {
        isl_constraint *cst = isl_constraint_list_get_constraint(list, i);
        isl_val *val = isl_constraint_get_coefficient_val(cst, isl_dim_out, data->in_dim);
        if (isl_val_is_one(val)) // i.e., the coefficient of the dimension data->in_dim is 1
        {
            isl_val *val2 = isl_constraint_get_constant_val(cst);
            int const_val = (-1) * isl_val_get_num_si(val2);
            data->out_constant = const_val;
            // DEBUG(3, daisy::str_dump("Dimensions found.  Constant = " +
            //                             std::to_string(const_val)));
        }
    }

    return isl_stat_ok;
}

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
isl_ast_expr* utility::extract_bound_expression(isl_ast_node *node, int dim, bool upper)
{
    assert(node != NULL);
    assert(dim >= 0);

    isl_ast_expr* result;
    if (isl_ast_node_get_type(node) == isl_ast_node_for)
    {
        isl_ast_expr *init_bound = isl_ast_node_for_get_init(node);
        isl_ast_expr *upper_bound = isl_ast_node_for_get_cond(node);
        if (dim == 0)
        {
            if (upper)
            {
                isl_ast_expr *cond = isl_ast_node_for_get_cond(node);

                /**
                  * If we have an expression
                  *  i < N
                  * or an expression
                  *  i <= N - 1
                  *
                  * In both cases, the returned bound should be (N-1).
                  */
                if (isl_ast_expr_get_op_type(cond) == isl_ast_op_lt)
                {
                    // Create an expression of "1".
                    isl_val *one = isl_val_one(isl_ast_node_get_ctx(node));
                    // Add 1 to the ISL ast upper bound to transform it into a strict bound.
                    result = isl_ast_expr_sub(isl_ast_expr_get_op_arg(cond, 1), isl_ast_expr_from_val(one));
                }
                else if (isl_ast_expr_get_op_type(cond) == isl_ast_op_le)
                {
                    result = isl_ast_expr_get_op_arg(cond, 1);
                }
            }
            else
            {
                result = isl_ast_node_for_get_init(node);
            }
        }
        else
        {
            isl_ast_node *body = isl_ast_node_for_get_body(node);
            result = utility::extract_bound_expression(body, dim-1, upper);
            isl_ast_node_free(body);
        }
    }
    else if (isl_ast_node_get_type(node) == isl_ast_node_if)
    {
        isl_ast_expr* cond_bound = isl_ast_node_if_get_cond(node);
        isl_ast_expr* then_bound = utility::extract_bound_expression(isl_ast_node_if_get_then(node), dim, upper);
        if (isl_ast_node_if_has_else(node))
        {
            isl_ast_expr* else_bound = utility::extract_bound_expression(isl_ast_node_if_get_else(node), dim, upper);
            result = else_bound;
        }
        else
            result = then_bound;
    } else {
        // Problem
    }

    return result;
}

/**
 * - Generate code:
 * - Generate time-processor domain.
 * - Generate an ISL AST.
 * - Traverse the tree until the level \p dim.
 * - Extract the bounds of that level.
 * - During the traversal, assert that the loop is fully nested.
 *
 */
int utility::get_bound(isl_set *set, int dim, int upper)
{
    assert(set != NULL);
    assert(dim >= 0);
    assert(dim < isl_space_dim(isl_set_get_space(set), isl_dim_set));
    assert(isl_set_is_empty(set) == isl_bool_false);

    isl_ast_build *ast_build;
    isl_ctx *ctx = isl_set_get_ctx(set);
    ast_build = isl_ast_build_alloc(ctx);

    // Create identity map for set.
    isl_space *sp = isl_set_get_space(set);
    isl_map *sched = isl_map_identity(isl_space_copy(isl_space_map_from_set(sp)));
    sched = isl_map_set_tuple_name(sched, isl_dim_out, "");

    // Generate the AST.
    isl_options_set_ast_build_atomic_upper_bound(ctx, 1);
    isl_options_get_ast_build_exploit_nested_bounds(ctx);
    isl_options_set_ast_build_group_coscheduled(ctx, 1);
    isl_options_set_ast_build_allow_else(ctx, 1);
    isl_options_set_ast_build_detect_min_max(ctx, 1);

    // Computing the polyhedral hull of the input set.
    //DEBUG(3, daisy::str_dump("Computing the polyhedral hull of the input set."));
    //set = isl_set_from_basic_set(isl_set_affine_hull(isl_set_copy(set)));
    //DEBUG(3, daisy::str_dump("The polyhedral hull is: ", isl_set_to_str(set)));

    // Intersect the iteration domain with the domain of the schedule.
    isl_map *map = isl_map_intersect_domain(isl_map_copy(sched), isl_set_copy(set));

    // Set iterator names
    int length = isl_map_dim(map, isl_dim_out);
    isl_id_list *iterators = isl_id_list_alloc(ctx, length);

    for (int i = 0; i < length; i++)
    {
        std::string name;
        if (isl_set_has_dim_name(set, isl_dim_set, i) == true)
            name = isl_set_get_dim_name(set, isl_dim_set, i);
        else
            name = generate_new_variable_name();
        isl_id *id = isl_id_alloc(ctx, name.c_str(), NULL);
        iterators = isl_id_list_add(iterators, id);
    }

    ast_build = isl_ast_build_set_iterators(ast_build, iterators);

    isl_ast_node *node = isl_ast_build_node_from_schedule_map(ast_build, isl_union_map_from_map(map));
    isl_ast_expr* bound = utility::extract_bound_expression(node, dim, upper);
    isl_ast_build_free(ast_build);

    // Convert to int
    isl_val *init_val = isl_ast_expr_get_val(bound);
    int b = (int) isl_val_get_num_si(init_val);
    isl_val_free(init_val);
    
    return b;
}

/**
 * Transform a dynamic schedule dimension into the corresponding loop level.
 *
 * In the example below, the loop level that corresponds
 * to the dynamic dimension 2 is 0, and to the dynamic dimension 4 is 1, ...
 *
 * The first dimension is the duplication dimension, the following
 * dimensions are static, dynamic, static, dynamic, ...
 *
 * Loop level               :    -1         0      1      2
 * Schedule dimension number:        0, 1   2  3   4  5   6  7
 * Schedule:                        [0, 0, i1, 0, i2, 0, i3, 0]
 */
int dynamic_dimension_into_loop_level(int dim)
{
    assert(dim % 2 == 0);
    int level = (dim - 2)/2;
    return level;
}


} // namespace daisy

