#include "../include/computation.h"

#include <cassert>
#include <queue>
#include <iostream>

#include "../include/core.h"
#include "../include/function.h"

namespace daisy {

namespace tiramisu {

/**
 * Transform the loop level into its corresponding dynamic schedule
 * dimension.
 *
 * In the example below, the dynamic dimension that corresponds
 * to the loop level 0 is 2, and to 1 it is 4, ...
 *
 * The first dimension is the duplication dimension, the following
 * dimensions are static, dynamic, static, dynamic, ...
 *
 * Loop level               :    -1         0      1      2
 * Schedule dimension number:        0, 1   2  3   4  5   6  7
 * Schedule:                        [0, 0, i1, 0, i2, 0, i3, 0]
 */
int loop_level_into_dynamic_dimension(int level)
{
    return 1 + (level * 2 + 1);
}

/**
 * Transform the loop level into the first static schedule
 * dimension after its corresponding dynamic dimension.
 *
 * In the example below, the first static dimension that comes
 * after the corresponding dynamic dimension for
 * the loop level 0 is 3, and to 1 it is 5, ...
 *
 * Loop level               :    -1         0      1      2
 * Schedule dimension number:        0, 1   2  3   4  5   6  7
 * Schedule:                        [0, 0, i1, 0, i2, 0, i3, 0]
 */
int loop_level_into_static_dimension(int level)
{
    return loop_level_into_dynamic_dimension(level) + 1;
}


// TODO: fix this function
isl_map *add_eq_to_schedule_map(int dim0, int in_dim_coefficient, int out_dim_coefficient,
                                int const_conefficient, isl_map *sched)
{
    isl_map *identity = isl_set_identity(isl_map_range(isl_map_copy(sched)));
    identity = isl_map_universe(isl_map_get_space(identity));
    isl_space *sp = isl_map_get_space(identity);
    isl_local_space *lsp = isl_local_space_from_space(isl_space_copy(sp));

    // Create a transformation map that transforms the schedule.
    for (int i = 0; i < isl_map_dim (identity, isl_dim_out); i++)
        if (i == dim0)
        {
            isl_constraint *cst = isl_constraint_alloc_equality(isl_local_space_copy(lsp));
            cst = isl_constraint_set_coefficient_si(cst, isl_dim_in, dim0, in_dim_coefficient);
            cst = isl_constraint_set_coefficient_si(cst, isl_dim_out, dim0, -out_dim_coefficient);
            // TODO: this should be inverted into const_conefficient.
            cst = isl_constraint_set_constant_si(cst, -const_conefficient);
            identity = isl_map_add_constraint(identity, cst);
        }
        else
        {
            // Set equality constraints for dimensions
            isl_constraint *cst2 = isl_constraint_alloc_equality(isl_local_space_copy(lsp));
            cst2 = isl_constraint_set_coefficient_si(cst2, isl_dim_in, i, 1);
            cst2 = isl_constraint_set_coefficient_si(cst2, isl_dim_out, i, -1);
            identity = isl_map_add_constraint(identity, cst2);
        }

    isl_map *final_identity = identity;
    // DEBUG_INDENT(-4);

    return sched;
}

int compute_recursively_max_AST_depth(isl_ast_node *node)
{
    assert(node != NULL);

    // DEBUG_FCT_NAME(10);
    // DEBUG_INDENT(4);

    int result = -1;

    // DEBUG(10, daisy::str_dump("Computing maximal AST depth from the following ISL AST node "););
    // DEBUG(10, daisy::str_dump("\n"); daisy::str_dump(std::string(isl_ast_node_to_C_str(node))));

    if (isl_ast_node_get_type(node) == isl_ast_node_block)
    {
        // DEBUG(10, daisy::str_dump("Computing maximal depth from a block."));

        isl_ast_node_list *list = isl_ast_node_block_get_children(node);
        isl_ast_node *child = isl_ast_node_list_get_ast_node(list, 0);
        result = compute_recursively_max_AST_depth(child);

        for (int i = 1; i < isl_ast_node_list_n_ast_node(list); i++)
        {
            child = isl_ast_node_list_get_ast_node(list, i);
            result = std::max(result, compute_recursively_max_AST_depth(child));
        }
    }
    else if (isl_ast_node_get_type(node) == isl_ast_node_for)
    {
        // DEBUG(10, daisy::str_dump("Computing maximal depth from a for loop."));
        isl_ast_node *body = isl_ast_node_for_get_body(node);
        result = compute_recursively_max_AST_depth(body) + 1;
        isl_ast_node_free(body);
    }
    else if (isl_ast_node_get_type(node) == isl_ast_node_user)
    {
        // DEBUG(10, daisy::str_dump("Reached a user node."));
        result = 1;
    }
    else if (isl_ast_node_get_type(node) == isl_ast_node_if)
    {
        // DEBUG(10, daisy::str_dump("Computing maximal depth from an if conditional."));

        result = compute_recursively_max_AST_depth(isl_ast_node_if_get_then(node));

        if (isl_ast_node_if_has_else(node))
            result = std::max(result, compute_recursively_max_AST_depth(isl_ast_node_if_get_else(node)));
    }
    // else
    // {
    //     ERROR("Found an unsupported ISL AST node while computing the maximal AST depth.", true);
    // }

    // DEBUG(3, daisy::str_dump("Current depth = " + std::to_string(result)));
    // DEBUG_INDENT(-4);

    return result;
}


/**
 * Dummy constructor for derived classes.
 */
computation::computation()
{
    this->access = NULL;
    this->schedule = NULL;
    // this->stmt = Halide::Internal::Stmt();
    this->time_processor_domain = NULL;
    this->duplicate_number = 0;

    this->schedule_this_computation = false;
    this->data_type = p_none;
    // this->expression = daisy::expr();

    this->ctx = NULL;

    // this->lhs_access_type = daisy::o_access;
    // this->lhs_argument_idx = -1;
    // this->rhs_argument_idx = -1;
    // this->wait_argument_idx = -1;
    // this->_is_library_call = false;
    // this->wait_access_map = nullptr;
    // this->wait_index_expr = nullptr;

    this->iteration_domain = NULL;
    this->name = "";
    this->fct = NULL;
    this->is_let = false;
}

computation::computation(
      function *fct,
      std::string iteration_domain,
      std::string access,
      primitive_t data_type,
      dace_desc desc
) : fct(fct), data_type(data_type), desc(desc)
{
    init_computation(iteration_domain, access);
}

void computation::init_computation(std::string iteration_space_str, std::string access_str)
{


    // DEBUG(3, daisy::str_dump("Constructing the computation: " + iteration_space_str));

    assert(this->fct != NULL);
    assert(iteration_space_str.length() > 0 && ("Empty iteration space"));

    // Initialize all the fields to NULL (useful for later asserts)
    this->access = NULL;
    this->time_processor_domain = NULL;
    this->duplicate_number = 0;
    this->automatically_allocated_buffer = NULL;
    // In the constructor of computations, we assume that every created
    // computation is the first computation, then, if this computation
    // was created by add_definitions(), we change is_first_definition
    // to false, otherwise we keep it.
    // We do the same for first_definition.
    is_first = true;
    first_definition = NULL;
    this->definitions_number = 1;
    this->definition_ID = 0;

    this->schedule_this_computation = true;

    this->ctx = this->fct->get_isl_ctx();

    iteration_domain = isl_set_read_from_str(ctx, iteration_space_str.c_str());
    this->name_unnamed_iteration_domain_dimensions();
    this->name = std::string(isl_space_get_tuple_name(isl_set_get_space(iteration_domain), isl_dim_type::isl_dim_set));

    number_of_dims = isl_set_dim(iteration_domain, isl_dim_type::isl_dim_set);
    for (unsigned i = 0; i < number_of_dims; i++) {
        if (isl_set_has_dim_name(iteration_domain, isl_dim_type::isl_dim_set, i)) {
            std::string dim_name(isl_set_get_dim_name(iteration_domain, isl_dim_type::isl_dim_set, i));
            this->access_variables.push_back(make_pair(i, dim_name));
        }
    }

    this->fct->add_computation(this);
    this->set_identity_schedule_based_on_iteration_domain();

    this->is_inline = false;
    this->is_let = false;

    // Set the names of output dimensions to be equal to the names of iteration domain schedules.
    std::vector<std::string> nms = this->get_iteration_domain_dimension_names();
    // Rename the dimensions of the schedule domain so that when we set the names of
    // the schedule range dimension to be equal to the names of the domain, we do not
    // get a conflict.
    for (int i = 0; i< this->get_iteration_domain_dimensions_number(); i++)
        this->set_schedule_domain_dim_names({i}, {generate_new_variable_name()});
    for (int i = 0; i< nms.size(); i++)
        this->set_loop_level_names({i}, {nms[i]});

    // If there are computations that have already been defined and that
    // have the same name, check that they have constraints over their iteration
    // domains.
    std::vector<computation *> same_name_computations =
        this->fct->get_computation_by_name(name);
    if (same_name_computations.size() > 1)
    {
        if (isl_set_plain_is_universe(this->get_iteration_domain()))
        {
            // ERROR("Computations defined multiple times should"
            //                 " have bounds on their iteration domain", true);
        }

        for (auto c : same_name_computations)
        {
            if (isl_set_plain_is_universe(c->get_iteration_domain()))
            {
                // ERROR("Computations defined multiple times should"
                //                 " have bounds on their iteration domain", true);
            }
        }
    }

    this->updates.push_back(this);


    this->set_access(access_str);

    // DEBUG_INDENT(-4);
}

const std::string& computation::get_name() const
{
    return this->name;
}

function* computation::get_function() const
{
    return this->fct;
}

primitive_t computation::get_data_type() const
{
    return this->data_type;
}

isl_set* computation::get_iteration_domain() const
{
    return this->iteration_domain;
}

void computation::set_iteration_domain(isl_set *domain)
{
    this->iteration_domain = domain;
}

isl_map* computation::get_schedule()
{
    return this->schedule;
}

isl_ctx* computation::get_ctx() const
{
    return this->ctx;
}

isl_map *computation::get_access_relation() const
{
    return this->access;
}

void computation::set_access(isl_map *access)
{
    assert(access != NULL);

    this->set_access(isl_map_to_str(access));
}

/**
 * Set the access function of the computation.
 *
 * The access function is a relation from computations to buffer locations.
 * \p access_str is a string that represents the relation (in ISL format,
 * http://isl.gforge.inria.fr/user.html#Sets-and-Relations).
 */
void computation::set_access(std::string access_str)
{


    // DEBUG(3, daisy::str_dump("Setting access " + access_str + " for computation " + this->get_name()));

    this->access = isl_map_read_from_str(this->ctx, access_str.c_str());

    /**
     * Set the access relations of all the computations that have the same name
     * (duplicates and updates).
     */
    std::vector<computation *> same_name_computations =
        this->get_function()->get_computation_by_name(this->get_name());

    if (same_name_computations.size() > 1)
        for (auto c : same_name_computations)
        {
            c->access = isl_map_read_from_str(this->ctx, access_str.c_str());
        }

    /**
     * Check that if there are other computations that have the same name
     * as this computation, then the access of all of these computations
     * should be the same.
     */
    std::vector<computation *> computations =
        this->get_function()->get_computation_by_name(this->get_name());
    for (auto c : computations)
        if (isl_map_is_equal(this->get_access_relation(), c->get_access_relation()) == isl_bool_false)
        {
            // ERROR("Computations that have the same name should also have the same access relation.",
            //                 true);
        }

    assert(this->access != nullptr && "Set access failed");

    // DEBUG_INDENT(-4);
}

std::shared_ptr<buffer> computation::get_buffer() const
{
    if (this->access == NULL)
    {
        return nullptr;
    }

    std::string buffer_name = isl_map_get_tuple_name(this->access, isl_dim_out);
    assert((this->get_function()->get_buffers().count(buffer_name) > 0) && ("Buffer does not exist"));
    return this->get_function()->get_buffers().find(buffer_name)->second;
}

bool computation::should_schedule_this_computation() const
{
    return schedule_this_computation;
}

int computation::compute_maximal_AST_depth()
{
    // DEBUG_FCT_NAME(10);
    // DEBUG_INDENT(4);

    this->name_unnamed_time_space_dimensions();
    this->gen_time_space_domain();
    isl_set *set = this->get_trimmed_time_processor_domain();
    assert(set != NULL);

    // DEBUG(10, daisy::str_dump(std::string("Getting the ") +
    //                              " maximal AST depth of the set ",
    //                              isl_set_to_str(set)));

    isl_ast_build *ast_build;
    isl_ctx *ctx = isl_set_get_ctx(set);
    ast_build = isl_ast_build_alloc(ctx);

    // Create identity map for set.
    isl_space *sp = isl_set_get_space(set);
    isl_map *sched = isl_map_identity(isl_space_copy(isl_space_map_from_set(sp)));
    sched = isl_map_set_tuple_name(sched, isl_dim_out, "");

    // Generate the AST.
    // DEBUG(10, daisy::str_dump("Setting ISL AST generator options."));
    isl_options_set_ast_build_atomic_upper_bound(ctx, 1);
    isl_options_get_ast_build_exploit_nested_bounds(ctx);
    isl_options_set_ast_build_group_coscheduled(ctx, 1);
    isl_options_set_ast_build_allow_else(ctx, 1);
    isl_options_set_ast_build_detect_min_max(ctx, 1);

    // Intersect the iteration domain with the domain of the schedule.
    // DEBUG(10, daisy::str_dump("Generating time-space domain."));
    isl_map *map = isl_map_intersect_domain(isl_map_copy(sched), isl_set_copy(set));

    // Set iterator names
    // DEBUG(10, daisy::str_dump("Setting the iterator names."));
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

    int depth = compute_recursively_max_AST_depth(node);
    
    isl_ast_build_free(ast_build);

    // DEBUG(10, daisy::str_dump("The maximal AST depth is : " + std::to_string(depth)));
    // DEBUG_INDENT(-4);

    return depth;
}

int computation::get_iteration_domain_dimensions_number()
{
    assert(iteration_domain != NULL);

    return isl_set_n_dim(this->iteration_domain);
}

std::vector<std::string> computation::get_iteration_domain_dimension_names()
{


    isl_set *iter = this->get_iteration_domain();

    assert(iter != NULL);

    std::vector<std::string> result;

    for (int i = 0; i < this->get_iteration_domain_dimensions_number(); i++)
    {
        if (isl_set_has_dim_name(iter, isl_dim_set, i))
            result.push_back(std::string(isl_set_get_dim_name(iter,
                                                              isl_dim_set, i)));
        else
        {
            // ERROR("All iteration domain dimensions must have "
            //                 "a name.", true);
        }
    }

    assert(result.size() == this->get_iteration_domain_dimensions_number());

    // DEBUG_INDENT(-4);

    return result;
}

void computation::name_unnamed_iteration_domain_dimensions()
{


    isl_set *iter = this->get_iteration_domain();

    assert(iter != NULL);

    for (int i = 0; i < this->get_iteration_domain_dimensions_number(); i++)
    {
        if (isl_set_has_dim_name(iter, isl_dim_set, i) == isl_bool_false)
            iter = isl_set_set_dim_name(iter, isl_dim_set, i,
                                        generate_new_variable_name().c_str());
    }

    this->set_iteration_domain(iter);

    // DEBUG_INDENT(-4);
}


void computation::set_schedule_domain_dim_names(std::vector<int> loop_levels, std::vector<std::string> names)
{


    this->check_dimensions_validity(loop_levels);
    assert(names.size() > 0);
    assert(names.size() == loop_levels.size());

    for (int i = 0; i < loop_levels.size(); i++)
    {
        assert(loop_levels[i] <= isl_map_dim(this->get_schedule(), isl_dim_in));
        this->schedule = isl_map_set_dim_name(this->get_schedule(),
                                              isl_dim_in, loop_levels[i], names[i].c_str());
        // DEBUG(3, daisy::str_dump("Setting the name of the domain of the schedule dimension " + std::to_string(loop_levels[i]) + " into " + names[i].c_str()));
    }

    // DEBUG(3, daisy::str_dump("The schedule after renaming: ", isl_map_to_str(this->get_schedule())));

    // DEBUG_INDENT(-4);
}

int computation::get_loop_levels_number()
{
    int loop_levels_number = ((isl_map_dim(this->get_schedule(), isl_dim_out)) - 2)/2;

    return loop_levels_number;
}

std::vector<std::string> computation::get_loop_level_names()
{
    std::vector<std::string> names;
    for (int i = 0; i < this->get_loop_levels_number(); i++)
    {
        std::string dim_name = isl_map_get_dim_name(this->get_schedule(), isl_dim_out, loop_level_into_dynamic_dimension(i));
        names.push_back(dim_name);
    }

    return names;
}

std::string computation::get_dimension_name_for_loop_level(int loop_level)
{
    int dim = loop_level_into_dynamic_dimension(loop_level);
    std::string name = isl_map_get_dim_name(this->get_schedule(), isl_dim_out, dim);
    assert(name.size() > 0);
    return name;
}

void computation::set_schedule(isl_map *map)
{
    this->schedule = map;
}

void computation::set_identity_schedule_based_on_iteration_domain()
{
    isl_map* sched = this->gen_identity_schedule_for_iteration_domain();
    this->set_schedule(sched);
}


isl_map* computation::gen_identity_schedule_for_iteration_domain()
{
    isl_space *sp = isl_set_get_space(this->get_iteration_domain());
    isl_map *sched = isl_map_identity(isl_space_map_from_set(sp));
    sched = isl_map_intersect_domain(sched, isl_set_copy(this->get_iteration_domain()));
    sched = isl_map_coalesce(sched);

    // Add Beta dimensions.
    for (int i = 0; i < isl_space_dim(sp, isl_dim_out) + 1; i++)
    {
        sched = isl_map_add_dim_and_eq_constraint(sched, 2 * i, 0);
    }

    // Add the duplication dimension.
    sched = isl_map_add_dim_and_eq_constraint(sched, 0, 0);
    return sched;
}

void computation::add_schedule_constraint(std::string domain_constraints, std::string range_constraints)
{



    assert(this->ctx != NULL);
    isl_map *sched = this->get_schedule();

    if (!domain_constraints.empty())
    {
        isl_set *domain_cst = isl_set_read_from_str(this->ctx, domain_constraints.c_str());
        assert(domain_cst != NULL);

        // DEBUG(3, daisy::str_dump("Adding the following constraints to the domain of the schedule : "));
        // DEBUG(3, daisy::str_dump(isl_set_to_str(domain_cst)));
        // DEBUG(3, daisy::str_dump("The schedule is : "));
        // DEBUG(3, daisy::str_dump(isl_map_to_str(sched)));

        sched = isl_map_intersect_domain(isl_map_copy(sched), isl_set_copy(domain_cst));

    }

    if (!range_constraints.empty())
    {
        isl_set *range_cst = isl_set_read_from_str(this->ctx, range_constraints.c_str());

        // DEBUG(3, daisy::str_dump("Adding the following constraints to the range of the schedule : "));
        // DEBUG(3, daisy::str_dump(isl_set_to_str(range_cst)));
        // DEBUG(3, daisy::str_dump("The schedule : ", isl_map_to_str(sched)));

        sched = isl_map_intersect_range(isl_map_copy(sched), isl_set_copy(range_cst));
    }

    this->set_schedule(sched);

    // DEBUG(3, daisy::str_dump("Schedule after transformation : "));
    // DEBUG(3, daisy::str_dump(isl_map_to_str(this->get_schedule())));

    // DEBUG_INDENT(-4);
}

void computation::gen_time_space_domain()
{


    assert(this->get_iteration_domain() != NULL);
    assert(this->get_schedule() != NULL);

    // DEBUG(3, daisy::str_dump("Iteration domain:", isl_set_to_str(this->get_iteration_domain())));

    isl_set *iter = isl_set_copy(this->get_iteration_domain());
    iter = this->intersect_set_with_context(iter);

    // DEBUG(3, daisy::str_dump("Iteration domain Intersect context:", isl_set_to_str(iter)));

    time_processor_domain = isl_set_apply(
                                iter,
                                isl_map_copy(this->get_schedule()));

    // DEBUG(3, daisy::str_dump("Schedule:", isl_map_to_str(this->get_schedule())));
    // DEBUG(3, daisy::str_dump("Generated time-space domain:", isl_set_to_str(time_processor_domain)));

    // DEBUG_INDENT(-4);
}

isl_set *computation::get_time_processor_domain() const
{
    return time_processor_domain;
}

isl_set* computation::get_trimmed_time_processor_domain()
{
    isl_set *tp_domain = isl_set_copy(this->get_time_processor_domain());
    const char *name = isl_set_get_tuple_name(isl_set_copy(tp_domain));
    isl_set *tp_domain_without_duplicate_dim = isl_set_project_out(isl_set_copy(tp_domain), isl_dim_set, 0, 1);
    tp_domain_without_duplicate_dim = isl_set_set_tuple_name(tp_domain_without_duplicate_dim, name);
    return tp_domain_without_duplicate_dim ;
}

isl_set* computation::intersect_set_with_context(isl_set *set)
{


    // Unify the space of the context and the "missing" set so that we can intersect them.
    isl_set *context = isl_set_copy(this->get_function()->get_program_context());
    if (context != NULL)
    {
        isl_space *model = isl_set_get_space(isl_set_copy(context));
        set = isl_set_align_params(set, isl_space_copy(model));
        // DEBUG(10, daisy::str_dump("Context: ", isl_set_to_str(context)));
        // DEBUG(10, daisy::str_dump("Set after aligning its parameters with the context parameters: ",
        //                              isl_set_to_str (set)));

        isl_id *missing_id1 = NULL;
        if (isl_set_has_tuple_id(set) == isl_bool_true)
        {
            missing_id1 = isl_set_get_tuple_id(set);
        }
        else
        {
            std::string name = isl_set_get_tuple_name(set);
            assert(name.size() > 0);
            missing_id1 = isl_id_alloc(this->get_ctx(), name.c_str(), NULL);
        }

        int nb_dims = isl_set_dim(set, isl_dim_set);
        context = isl_set_add_dims(context, isl_dim_set, nb_dims);
        // DEBUG(10, daisy::str_dump("Context after adding dimensions to make it have the same number of dimensions as missing: ",
        //                              isl_set_to_str (context)));
        context = isl_set_set_tuple_id(context, isl_id_copy(missing_id1));
        // DEBUG(10, daisy::str_dump("Context after setting its tuple ID to be equal to the tuple ID of missing: ",
        //                              isl_set_to_str (context)));
        set = isl_set_intersect(set, isl_set_copy(context));
        // DEBUG(10, daisy::str_dump("Set after intersecting with the program context: ",
        //                              isl_set_to_str (set)));
    }

    // DEBUG_INDENT(-4);

    return set;
}

std::vector<int> computation::get_loop_level_numbers_from_dimension_names(
        std::vector<std::string> dim_names)
{


    assert(dim_names.size() > 0);

    std::vector<int> dim_numbers;

    for (auto const& dim : dim_names)
    {
        assert(dim.size()>0);
        
        // DEBUG(10, daisy::str_dump("Searching for the dimension " + dim));

        if (dim == "root")
        {
            int d = computation::root_dimension;
            dim_numbers.push_back(d);
        }
        else
        {
            int d = isl_map_find_dim_by_name(this->get_schedule(), isl_dim_out, dim.c_str());
            // DEBUG(10, daisy::str_dump("Searching in the range of ",
            //                              isl_map_to_str(this->get_schedule())));

            if (d < 0)
            {
                return std::vector<int>{-2};
                // ERROR("Dimension " + dim + " not found.", true);
            }

            // DEBUG(10, daisy::str_dump("Corresponding loop level is " +
            //                              std::to_string(dynamic_dimension_into_loop_level(d))));

            dim_numbers.push_back(dynamic_dimension_into_loop_level(d));
        }
    }

    this->check_dimensions_validity(dim_numbers);
    return dim_numbers;
}

void computation::check_dimensions_validity(std::vector<int> dimensions)
{
    assert(dimensions.size() > 0);

    for (auto const& dim : dimensions)
    {
        // DEBUG(10, daisy::str_dump("Checking the validity of loop level " +
        //                              std::to_string(dim)));

        assert(dim >= computation::root_dimension);

        if (loop_level_into_dynamic_dimension(dim) >=
            isl_space_dim(isl_map_get_space(this->get_schedule()),
                          isl_dim_out))
        {
            // ERROR("The dynamic dimension " +
            //                 std::to_string(loop_level_into_dynamic_dimension(dim)) +
            //                 " is not less than the number of dimensions of the "
            //                 "time-space domain " +
            //                 std::to_string(isl_space_dim(isl_map_get_space(
            //                         this->get_schedule()), isl_dim_out)), true);
        }
    }
}

void computation::name_unnamed_time_space_dimensions()
{
    isl_map *sched = this->get_schedule();

    assert(sched != NULL);

    for (int i = 0; i < this->get_loop_levels_number(); i++)
    {
        if (isl_map_has_dim_name(sched, isl_dim_out, loop_level_into_dynamic_dimension(i)) == isl_bool_false)
            sched = isl_map_set_dim_name(sched, isl_dim_out, loop_level_into_dynamic_dimension(i), generate_new_variable_name().c_str());
    }
    this->set_schedule(sched);
}

void computation::set_loop_level_names(std::vector<std::string> names)
{


    assert(names.size() > 0);

    // DEBUG(3, daisy::str_dump("Number of loop levels: " + std::to_string(this->get_loop_levels_number())));
    // DEBUG(3, daisy::str_dump("Number of names to be set: " + std::to_string(names.size())));

    for (int i = 0; i < names.size(); i++)
    {
        if (isl_map_has_dim_name(this->get_schedule(), isl_dim_out, loop_level_into_dynamic_dimension(i)) == isl_bool_true)
        {
            this->schedule = isl_map_set_dim_name(this->get_schedule(),
                                                  isl_dim_out,
                                                  loop_level_into_dynamic_dimension(i),
                                                  names[i].c_str());
            // DEBUG(3, daisy::str_dump("Setting the name of loop level " + std::to_string(i) + " into " + names[i].c_str()));
        }
    }

    // DEBUG(3, daisy::str_dump("The schedule after renaming: ", isl_map_to_str(this->get_schedule())));

    // DEBUG_INDENT(-4);
}

void computation::set_loop_level_names(std::vector<int> loop_levels, std::vector<std::string> names)
{


    this->check_dimensions_validity(loop_levels);
    assert(names.size() > 0);
    assert(names.size() == loop_levels.size());

    for (int i = 0; i < loop_levels.size(); i++)
    {
        if (loop_level_into_static_dimension(loop_levels[i]) <= isl_map_dim(this->get_schedule(), isl_dim_out))
        {
            this->schedule = isl_map_set_dim_name(this->get_schedule(),
                                                  isl_dim_out,
                                                  loop_level_into_dynamic_dimension(loop_levels[i]),
                                                  names[i].c_str());
            // DEBUG(3, daisy::str_dump("Setting the name of loop level " + std::to_string(loop_levels[i]) + " into " + names[i].c_str()));
        }
    }

    // DEBUG(3, daisy::str_dump("The schedule after renaming: ", isl_map_to_str(this->get_schedule())));

    // DEBUG_INDENT(-4);
}


void computation::assert_names_not_assigned(std::vector<std::string> dimensions)
{
    for (auto const dim: dimensions)
    {
        int d = isl_map_find_dim_by_name(this->get_schedule(), isl_dim_out,
                                         dim.c_str());
        if (d >= 0)
        {
            // ERROR("Dimension " + dim + " is already in use.", true);
        }

        d = isl_map_find_dim_by_name(this->get_schedule(), isl_dim_in,
                                     dim.c_str());
        if (d >= 0)
        {
            // ERROR("Dimension " + dim + " is already in use.", true);
        }
    }
}

void computation::update_names(std::vector<std::string> original_loop_level_names, std::vector<std::string> new_names, int erase_from, int nb_loop_levels_to_erase)
{


    // DEBUG_NO_NEWLINE(3, daisy::str_dump("Original loop level names: "));
    // for (auto n: original_loop_level_names)
    // {
    //     DEBUG_NO_NEWLINE_NO_INDENT(3, daisy::str_dump(n + " "));
    // }
    // DEBUG_NEWLINE(3);

    // DEBUG_NO_NEWLINE(3, daisy::str_dump("New names: "));
    // for (auto n: new_names)
    // {
    //     DEBUG_NO_NEWLINE_NO_INDENT(3, daisy::str_dump(n + " "));
    // }
    // DEBUG_NEWLINE(3);

    // DEBUG(3, daisy::str_dump("Start erasing from: " + std::to_string(erase_from)));
    // DEBUG(3, daisy::str_dump("Number of loop levels to erase: " + std::to_string(nb_loop_levels_to_erase)));

    original_loop_level_names.erase(original_loop_level_names.begin() + erase_from, original_loop_level_names.begin() + erase_from + nb_loop_levels_to_erase);

    // DEBUG_NO_NEWLINE(3, daisy::str_dump("Original loop level names after erasing loop levels: "));
    // for (auto n: original_loop_level_names)
    // {
    //     DEBUG_NO_NEWLINE_NO_INDENT(3, daisy::str_dump(n + " "));
    // }
    // DEBUG_NEWLINE(3);

    original_loop_level_names.insert(original_loop_level_names.begin() + erase_from, new_names.begin(), new_names.end());

    // DEBUG_NO_NEWLINE(3, daisy::str_dump("Original loop level names after inserting the new loop levels: "));
    // for (auto n: original_loop_level_names)
    // {
    //     DEBUG_NO_NEWLINE_NO_INDENT(3, daisy::str_dump(n + " "));
    // }
    // DEBUG_NEWLINE(3);

    this->set_loop_level_names(original_loop_level_names);
    this->name_unnamed_time_space_dimensions();

    // DEBUG(3, daisy::str_dump("Names updated. New names are: "));
    // for (auto n: this->get_loop_level_names())
    // {
    //     DEBUG_NO_NEWLINE_NO_INDENT(3, daisy::str_dump(n + " "));
    // }
    // DEBUG(3, daisy::str_dump(""));

    // DEBUG_INDENT(-4);
}

computation& computation::get_last_update()
{
    return this->get_update(this->get_updates().size()-1);
}

std::vector<computation*>& computation::get_updates() {
    return this->updates;
}

computation& computation::get_update(int i)
{
    return *(this->updates[i]);
}

void computation::tag_parallel_level(int par_dim)
{
    this->get_function()->add_parallel_dimension(this->get_name(), par_dim);
}

void computation::interchange(std::string L0_var, std::string L1_var)
{
    std::vector<int> dimensions = this->get_loop_level_numbers_from_dimension_names({L0_var, L1_var});
    this->check_dimensions_validity(dimensions);
    int L0 = dimensions[0];
    int L1 = dimensions[1];

    this->interchange(L0, L1);
}

void computation::interchange(int L0, int L1)
{


    int inDim0 = loop_level_into_dynamic_dimension(L0);
    int inDim1 = loop_level_into_dynamic_dimension(L1);

    assert(inDim0 >= 0);
    assert(inDim0 < isl_space_dim(isl_map_get_space(this->get_schedule()),
                                  isl_dim_out));
    assert(inDim1 >= 0);
    assert(inDim1 < isl_space_dim(isl_map_get_space(this->get_schedule()),
                                  isl_dim_out));

    isl_map *schedule = this->get_schedule();

    // DEBUG(3, daisy::str_dump("Original schedule: ", isl_map_to_str(schedule)));
    // DEBUG(3, daisy::str_dump("Interchanging the dimensions " + std::to_string(
                                    // L0) + " and " + std::to_string(L1)));

    int n_dims = isl_map_dim(schedule, isl_dim_out);

    std::string inDim0_str = isl_map_get_dim_name(schedule, isl_dim_out, inDim0);
    std::string inDim1_str = isl_map_get_dim_name(schedule, isl_dim_out, inDim1);

    std::vector<isl_id *> dimensions;

    // ------------------------------------------------------------
    // Create a map for the duplicate schedule.
    // ------------------------------------------------------------

    std::string map = "{ " + this->get_name() + "[";

    for (int i = 0; i < n_dims; i++)
    {
        if (i == 0)
        {
            int duplicate_ID = isl_map_get_static_dim(schedule, 0);
            map = map + std::to_string(duplicate_ID);
        }
        else
        {
            if (isl_map_get_dim_name(schedule, isl_dim_out, i) == NULL)
            {
                isl_id *new_id = isl_id_alloc(this->get_ctx(), generate_new_variable_name().c_str(), NULL);
                schedule = isl_map_set_dim_id(schedule, isl_dim_out, i, new_id);
            }

            map = map + isl_map_get_dim_name(schedule, isl_dim_out, i);
        }

        if (i != n_dims - 1)
        {
            map = map + ",";
        }
    }

    map = map + "] ->" + this->get_name() + "[";

    for (int i = 0; i < n_dims; i++)
    {
        if (i == 0)
        {
            int duplicate_ID = isl_map_get_static_dim(schedule, 0);
            map = map + std::to_string(duplicate_ID);
        }
        else
        {
            if ((i != inDim0) && (i != inDim1))
            {
                map = map + isl_map_get_dim_name(schedule, isl_dim_out, i);
                dimensions.push_back(isl_map_get_dim_id(schedule, isl_dim_out, i));
            }
            else if (i == inDim0)
            {
                map = map + inDim1_str;
                isl_id *id1 = isl_id_alloc(this->get_ctx(), inDim1_str.c_str(), NULL);
                dimensions.push_back(id1);
            }
            else if (i == inDim1)
            {
                map = map + inDim0_str;
                isl_id *id1 = isl_id_alloc(this->get_ctx(), inDim0_str.c_str(), NULL);
                dimensions.push_back(id1);
            }
        }

        if (i != n_dims - 1)
        {
            map = map + ",";
        }
    }

    map = map + "]}";

    // DEBUG(3, daisy::str_dump("A map that transforms the duplicate"));
    // DEBUG(3, daisy::str_dump(map.c_str()));

    isl_map *transformation_map = isl_map_read_from_str(this->get_ctx(), map.c_str());


    transformation_map = isl_map_set_tuple_id(
                             transformation_map, isl_dim_in, isl_map_get_tuple_id(isl_map_copy(schedule), isl_dim_out));
    isl_id *id_range = isl_id_alloc(this->get_ctx(), this->get_name().c_str(), NULL);
    transformation_map = isl_map_set_tuple_id(
                             transformation_map, isl_dim_out, id_range);


    // Check that the names of each dimension is well set
    for (int i = 1; i < isl_map_dim(transformation_map, isl_dim_in); i++)
    {
        isl_id *dim_id = isl_id_copy(dimensions[i - 1]);
        transformation_map = isl_map_set_dim_id(transformation_map, isl_dim_out, i, dim_id);
        assert(isl_map_has_dim_name(transformation_map, isl_dim_in, i));
        assert(isl_map_has_dim_name(transformation_map, isl_dim_out, i));
    }

    // DEBUG(3, daisy::str_dump("Final transformation map : ", isl_map_to_str(transformation_map)));

    schedule = isl_map_apply_range(isl_map_copy(schedule), isl_map_copy(transformation_map));

    // DEBUG(3, daisy::str_dump("Schedule after interchange: ", isl_map_to_str(schedule)));

    this->set_schedule(schedule);

    // DEBUG_INDENT(-4);
}

void computation::split(int L0, int sizeX)
{


    int inDim0 = loop_level_into_dynamic_dimension(L0);

    assert(this->get_schedule() != NULL);
    assert(inDim0 >= 0);
    assert(inDim0 < isl_space_dim(isl_map_get_space(this->get_schedule()), isl_dim_out));
    assert(sizeX >= 1);

    isl_map *schedule = this->get_schedule();
    
    int duplicate_ID = isl_map_get_static_dim(schedule, 0);

    schedule = isl_map_copy(schedule);
    schedule = isl_map_set_tuple_id(schedule, isl_dim_out,
                                    isl_id_alloc(this->get_ctx(), this->get_name().c_str(), NULL));


    // DEBUG(3, daisy::str_dump("Original schedule: ", isl_map_to_str(schedule)));
    // DEBUG(3, daisy::str_dump("Splitting dimension " + std::to_string(inDim0)
    //                             + " with split size " + std::to_string(sizeX)));

    std::string inDim0_str;

    std::string outDim0_str = generate_new_variable_name();
    std::string static_dim_str = generate_new_variable_name();
    std::string outDim1_str = generate_new_variable_name();

    int n_dims = isl_map_dim(this->get_schedule(), isl_dim_out);
    std::vector<isl_id *> dimensions;
    std::vector<std::string> dimensions_str;
    std::string map = "{";

    // -----------------------------------------------------------------
    // Preparing a map to split the duplicate computation.
    // -----------------------------------------------------------------

    map = map + this->get_name() + "[";

    for (int i = 0; i < n_dims; i++)
    {
        if (i == 0)
        {
            std::string dim_str = generate_new_variable_name();
            dimensions_str.push_back(dim_str);
            map = map + dim_str;
        }
        else
        {
            std::string dim_str = generate_new_variable_name();
            dimensions_str.push_back(dim_str);
            map = map + dim_str;

            if (i == inDim0)
            {
                inDim0_str = dim_str;
            }
        }

        if (i != n_dims - 1)
        {
            map = map + ",";
        }
    }

    map = map + "] -> " + this->get_name() + "[";

    for (int i = 0; i < n_dims; i++)
    {
        if (i == 0)
        {
            map = map + dimensions_str[i];
            dimensions.push_back(isl_id_alloc(
                                     this->get_ctx(),
                                     dimensions_str[i].c_str(),
                                     NULL));
        }
        else if (i != inDim0)
        {
            map = map + dimensions_str[i];
            dimensions.push_back(isl_id_alloc(
                                     this->get_ctx(),
                                     dimensions_str[i].c_str(),
                                     NULL));
        }
        else
        {
            map = map + outDim0_str + ", " + static_dim_str + ", " + outDim1_str;
            isl_id *id0 = isl_id_alloc(this->get_ctx(),
                                       outDim0_str.c_str(), NULL);
            isl_id *id2 = isl_id_alloc(this->get_ctx(),
                                       static_dim_str.c_str(), NULL);
            isl_id *id1 = isl_id_alloc(this->get_ctx(),
                                       outDim1_str.c_str(), NULL);
            dimensions.push_back(id0);
            dimensions.push_back(id2);
            dimensions.push_back(id1);
        }

        if (i != n_dims - 1)
        {
            map = map + ",";
        }
    }

    map = map + "] : " + dimensions_str[0] + " = " + std::to_string(duplicate_ID) + " and " +
          outDim0_str + " = floor(" + inDim0_str + "/" +
          std::to_string(sizeX) + ") and " + outDim1_str + " = (" +
          inDim0_str + "%" + std::to_string(sizeX) + ") and " + static_dim_str + " = 0}";

    isl_map *transformation_map = isl_map_read_from_str(this->get_ctx(), map.c_str());

    for (int i = 0; i < dimensions.size(); i++)
        transformation_map = isl_map_set_dim_id(
                                 transformation_map, isl_dim_out, i, isl_id_copy(dimensions[i]));

    transformation_map = isl_map_set_tuple_id(
                             transformation_map, isl_dim_in,
                             isl_map_get_tuple_id(isl_map_copy(schedule), isl_dim_out));
    isl_id *id_range = isl_id_alloc(this->get_ctx(), this->get_name().c_str(), NULL);
    transformation_map = isl_map_set_tuple_id(transformation_map, isl_dim_out, id_range);

    // DEBUG(3, daisy::str_dump("Transformation map : ",
    //                             isl_map_to_str(transformation_map)));

    schedule = isl_map_apply_range(isl_map_copy(schedule), isl_map_copy(transformation_map));

    // DEBUG(3, daisy::str_dump("Schedule after splitting: ", isl_map_to_str(schedule)));

    this->set_schedule(schedule);

    // DEBUG_INDENT(-4);
}

bool computation::separateAndSplit(int L0, int v)
{


    // DEBUG(3, daisy::str_dump("Applying separateAndSplit on loop level " + std::to_string(L0) + " with a split factor of " + std::to_string(v)));

    this->gen_time_space_domain();

    // Compute the depth before any scheduling.
    int original_depth = this->compute_maximal_AST_depth();

    // DEBUG(3, daisy::str_dump("Computing upper bound at loop level " + std::to_string(L0)));

    int loop_upper_bound = daisy::utility::get_bound(this->get_trimmed_time_processor_domain(), L0, true);

    // DEBUG(3, daisy::str_dump("Computing lower bound at loop level " + std::to_string(L0)));

    int loop_lower_bound = daisy::utility::get_bound(this->get_trimmed_time_processor_domain(), L0, false);

    // std::string lower_without_cast = std::to_string(loop_lower_bound);
    // while (lower_without_cast.find("cast") != std::string::npos) // while there is a "cast" in the expression
    // {
    //     // Remove "cast" from the string, we do not need it.
    //     // An alternative to this would be to actually mutate the expression N and remove the cast
    //     // operator, but that is more time consuming to implement than replacing the string directly.
    //     int pos = lower_without_cast.find("cast");
    //     lower_without_cast = lower_without_cast.erase(pos, 4);
    // }

    int loop_bound = loop_upper_bound - loop_lower_bound + 1;

    // DEBUG(3, daisy::str_dump("Loop bound for the loop to be separated and split: "); loop_bound.dump(false));

    /*
     * Separate this computation. That is, create two identical computations
     * where we have the constraint
     *     i < v * floor(loop_bound/v)
     * in the first and
     *     i >= v * floor(loop_bound/v)
     * in the second.
     *
     * The first is called the full computation while the second is called
     * the separated computation.
     * The two computations are identical in every thing except that they have
     * two different schedules.  Their schedule restricts them to a smaller domain
     * (the full or the separated domains) and schedule one after the other.
     */
    this->separate(L0, std::to_string(loop_bound), v, std::to_string(loop_lower_bound));

    // Make a copy of the schedule before splitting so that we revert the
    // schedule if splitting did not have any effect (i.e., did not happen).
    isl_map *sc = isl_map_copy(this->get_schedule());

    /**
     * Split the full computation since the full computation will be vectorized.
     */
    //this->get_update(0).split(L0, v);
    this->get_update(0).split_with_lower_bound(L0, v, std::to_string(loop_lower_bound));

    // Compute the depth after scheduling.
    int depth = this->compute_maximal_AST_depth();

    bool split_happened = false;
    if (depth == original_depth)
    {
        // DEBUG(3, daisy::str_dump("Split did not happen."));
        split_happened = false;

	// DEBUG(3, daisy::str_dump("Cancel splitting."));
	this->set_schedule(sc);
    }
    else
    {
         split_happened = true;
        //  DEBUG(3, daisy::str_dump("Split happenned."));
    }

    this->get_function()->align_schedules();

    // DEBUG_INDENT(-4);

    return split_happened;
}

void computation::split_with_lower_bound(int L0, int sizeX, std::string lower_bound)
{


    int inDim0 = loop_level_into_dynamic_dimension(L0);

    assert(this->get_schedule() != NULL);
    assert(inDim0 >= 0);
    assert(inDim0 < isl_space_dim(isl_map_get_space(this->get_schedule()), isl_dim_out));
    assert(sizeX >= 1);

    isl_map *schedule = this->get_schedule();
    int duplicate_ID = isl_map_get_static_dim(schedule, 0);

    schedule = isl_map_copy(schedule);
    schedule = isl_map_set_tuple_id(schedule, isl_dim_out,
                                    isl_id_alloc(this->get_ctx(), this->get_name().c_str(), NULL));


    // DEBUG(3, daisy::str_dump("Original schedule: ", isl_map_to_str(schedule)));
    // DEBUG(3, daisy::str_dump("Splitting dimension " + std::to_string(inDim0)
    //                             + " with split size " + std::to_string(sizeX)));

    std::string inDim0_str;

    std::string outDim0_str = generate_new_variable_name();
    std::string static_dim_str = generate_new_variable_name();
    std::string outDim1_str = generate_new_variable_name();

    int n_dims = isl_map_dim(this->get_schedule(), isl_dim_out);
    std::vector<isl_id *> dimensions;
    std::vector<std::string> dimensions_str;
    std::string map = "{";

    // -----------------------------------------------------------------
    // Preparing a map to split the duplicate computation.
    // -----------------------------------------------------------------

    map = map + this->get_name() + "[";

    for (int i = 0; i < n_dims; i++)
    {
        if (i == 0)
        {
            std::string dim_str = generate_new_variable_name();
            dimensions_str.push_back(dim_str);
            map = map + dim_str;
        }
        else
        {
            std::string dim_str = generate_new_variable_name();
            dimensions_str.push_back(dim_str);
            map = map + dim_str;

            if (i == inDim0)
            {
                inDim0_str = dim_str;
            }
        }

        if (i != n_dims - 1)
        {
            map = map + ",";
        }
    }

    map = map + "] -> " + this->get_name() + "[";

    for (int i = 0; i < n_dims; i++)
    {
        if (i == 0)
        {
            map = map + dimensions_str[i];
            dimensions.push_back(isl_id_alloc(
                                     this->get_ctx(),
                                     dimensions_str[i].c_str(),
                                     NULL));
        }
        else if (i != inDim0)
        {
            map = map + dimensions_str[i];
            dimensions.push_back(isl_id_alloc(
                                     this->get_ctx(),
                                     dimensions_str[i].c_str(),
                                     NULL));
        }
        else
        {
            map = map + outDim0_str + ", " + static_dim_str + ", " + outDim1_str;
            isl_id *id0 = isl_id_alloc(this->get_ctx(),
                                       outDim0_str.c_str(), NULL);
            isl_id *id2 = isl_id_alloc(this->get_ctx(),
                                       static_dim_str.c_str(), NULL);
            isl_id *id1 = isl_id_alloc(this->get_ctx(),
                                       outDim1_str.c_str(), NULL);
            dimensions.push_back(id0);
            dimensions.push_back(id2);
            dimensions.push_back(id1);
        }

        if (i != n_dims - 1)
        {
            map = map + ",";
        }
    }

    /**
     * The main idea is to avoid a case when the iterations dont start from 0.
     * let assume it start from 1, in that case after splitting of \p i into \p i1 & \p i2 with factor 4, the result would be:
     *  in first iteration of \p i1 would would compute \p i2 = [1,2,3,4] then [5,6,7,8]
     *  Whereas in the legacy version in first iteration of \p i1 would would compute \p i2 = [1,2,3]  then [4,5,6,7]
     * 
    */
    inDim0_str =  "("+inDim0_str +" - "+lower_bound+")";

    map = map + "] : " + dimensions_str[0] + " = " + std::to_string(duplicate_ID) + " and " +
          outDim0_str + " = floor(" + inDim0_str + "/" +
          std::to_string(sizeX) + ") and " + outDim1_str + " = (" +
          inDim0_str + " %" + std::to_string(sizeX) + ") and " + static_dim_str + " = 0}";


    isl_map *transformation_map = isl_map_read_from_str(this->get_ctx(), map.c_str());

    for (int i = 0; i < dimensions.size(); i++)
        transformation_map = isl_map_set_dim_id(
                                 transformation_map, isl_dim_out, i, isl_id_copy(dimensions[i]));

    transformation_map = isl_map_set_tuple_id(
                             transformation_map, isl_dim_in,
                             isl_map_get_tuple_id(isl_map_copy(schedule), isl_dim_out));
    isl_id *id_range = isl_id_alloc(this->get_ctx(), this->get_name().c_str(), NULL);
    transformation_map = isl_map_set_tuple_id(transformation_map, isl_dim_out, id_range);

    // DEBUG(3, daisy::str_dump("Transformation map : ",
    //                             isl_map_to_str(transformation_map)));

    schedule = isl_map_apply_range(isl_map_copy(schedule), isl_map_copy(transformation_map));

    // DEBUG(3, daisy::str_dump("Schedule after splitting: ", isl_map_to_str(schedule)));

    this->set_schedule(schedule);

    // DEBUG_INDENT(-4);
}

void computation::separate(int dim, std::string N, int v, std::string lower_bound)
{


    // DEBUG(3, daisy::str_dump("Separating the computation at level " + std::to_string(dim)));

    // DEBUG(3, daisy::str_dump("Generating the time-space domain."));
    // this->gen_time_space_domain();


    //////////////////////////////////////////////////////////////////////////////

    // We create the constraint (i < v*floor(N/v))
    // DEBUG(3, daisy::str_dump("Constructing the constraint (i<v*floor(N/v))"));
    // DEBUG(3, daisy::str_dump("Removing any cast operator in N."));
    std::string N_without_cast = N;
    while (N_without_cast.find("cast") != std::string::npos) // while there is a "cast" in the expression
    {
        // Remove "cast" from the string, we do not need it.
        // An alternative to this would be to actually mutate the expression N and remove the cast
        // operator, but that is more time consuming to implement than replacing the string directly.
        int pos = N_without_cast.find("cast");
        N_without_cast = N_without_cast.erase(pos, 4);
    }

    std::string lower_without_cast = lower_bound;
    while (lower_without_cast.find("cast") != std::string::npos) // while there is a "cast" in the expression
    {
        // Remove "cast" from the string, we do not need it.
        // An alternative to this would be to actually mutate the expression N and remove the cast
        // operator, but that is more time consuming to implement than replacing the string directly.
        int pos = lower_without_cast.find("cast");
        lower_without_cast = lower_without_cast.erase(pos, 4);
    }

    std::string constraint;
    constraint = "";
    for (int i=0; i<isl_map_dim(this->get_schedule(), isl_dim_param); i++)
    {
        if (i==0)
            constraint += "[";
        constraint += isl_map_get_dim_name(this->get_schedule(), isl_dim_param, i);
        if (i!=isl_map_dim(this->get_schedule(), isl_dim_param)-1)
            constraint += ",";
        else
            constraint += "]->";
    }
    constraint += "{" + this->get_name() + "[0,";
    for (int i=1; i<isl_map_dim(this->get_schedule(), isl_dim_out); i++)
    {
        if ((i%2==0) && (isl_map_has_dim_name(this->get_schedule(), isl_dim_out, i)==true))
            constraint += isl_map_get_dim_name(this->get_schedule(), isl_dim_out, i);
        else
            constraint += "o" + std::to_string(i);
        if (i != isl_map_dim(this->get_schedule(), isl_dim_out)-1)
            constraint += ",";
    }
    constraint += "]: ";

    std::string constraint1 = constraint +
                        this->get_dimension_name_for_loop_level(dim) + " < (" + std::to_string(v) + "*( floor((" + N_without_cast + ")/" + std::to_string(v) + ")) +"+lower_without_cast+" )}";
    // DEBUG(3, daisy::str_dump("The constraint is:" + constraint1));

    // We create the constraint (i >= v*floor(N/v))
    // DEBUG(3, daisy::str_dump("Constructing the constraint (i>=v*(floor(N/v)))"));
    std::string constraint2 = constraint +
                    this->get_dimension_name_for_loop_level(dim) + " >= (" + std::to_string(v) + "*(floor((" + N_without_cast + ")/" + std::to_string(v) + ")) +"+lower_without_cast+" )}";
    // DEBUG(3, daisy::str_dump("The constraint is:" + constraint2));

    //////////////////////////////////////////////////////////////////////////////

    isl_set *constraint2_isl = isl_set_read_from_str(this->get_ctx(), constraint2.c_str());

    if (isl_set_is_empty(isl_map_range(isl_map_intersect_range(isl_map_copy(this->get_schedule()), constraint2_isl))) == false)
    {
        // DEBUG(3, daisy::str_dump("The separate computation is not empty."));

        // Create the separated computation.
        // First, create the domain of the separated computation (which is identical to
        // the domain of the original computation). Both also have the same name.
        // TODO: create copy functions for all the classes so that we can copy the objects
        // we need to have this->get_expr().copy()

        std::string domain_str = std::string(isl_set_to_str(this->get_iteration_domain()));
        this->add_definitions(
            domain_str,
            this->should_schedule_this_computation(),
            this->get_data_type(),
            this->get_function()
        );

        // Set the schedule of the newly created computation (separated
        // computation) to be equal to the schedule of the original computation.
        isl_map *new_schedule = isl_map_copy(this->get_schedule());
        this->get_last_update().set_schedule(new_schedule);

        // Create the access relation of the separated computation.
        if (this->get_access_relation() != NULL)
        {
            // DEBUG(3, daisy::str_dump("Creating the access function of the separated computation.\n"));
            this->get_last_update().set_access(isl_map_copy(this->get_access_relation()));

            // DEBUG(3, daisy::str_dump("Access of the separated computation:",
            //                             isl_map_to_str(this->get_last_update().get_access_relation())));
        }

        this->get_last_update().add_schedule_constraint("", constraint2.c_str());

        // Mark the separated computation to be executed after the original (full)
        // computation.
        this->get_last_update().after(*this, dim);

        // DEBUG(3, daisy::str_dump("The separate computation:"); this->get_last_update().dump());
    }
    // else
    // {
    //     DEBUG(3, daisy::str_dump("The separate computation is empty. Thus not added."));
    // }

    this->add_schedule_constraint("", constraint1.c_str());

    // DEBUG(3, daisy::str_dump("The original computation:"); this->dump());

    // DEBUG_INDENT(-4);
}

void computation::add_definitions(std::string iteration_domain_str, bool schedule_this_computation, daisy::primitive_t t, function *fct)
{
    computation *new_c = new computation(
      this->fct,
      iteration_domain_str,
      isl_map_to_str(this->access),
      t,
      this->desc
    );
    new_c->is_first = false;
    new_c->first_definition = this;
    new_c->is_let = this->is_let;
    new_c->definition_ID = this->definitions_number;
    this->definitions_number++;

    this->updates.push_back(new_c);
}

void computation::after(computation &comp, int level)
{


    // DEBUG(3, daisy::str_dump("Scheduling " + this->get_name() + " to be executed after " +
    //                             comp.get_name() + " at level " + std::to_string(level)));

    auto &graph = this->get_function()->sched_graph;

    auto &edges = graph[&comp];

    auto level_it = edges.find(this);

    if (level_it != edges.end())
    {
        if (level_it->second > level)
        {
            level = level_it->second;
        }
    }

    edges[this] = level;

    this->get_function()->starting_computations.erase(this);

    this->get_function()->sched_graph_reversed[this][&comp] = level;

    assert(this->get_function()->sched_graph_reversed[this].size() < 2 &&
            "Node has more than one predecessor.");

    // DEBUG(10, daisy::str_dump("sched_graph[" + comp.get_name() + ", " +
    //                              this->get_name() + "] = " + std::to_string(level)));

    // DEBUG_INDENT(-4);
}

void computation::tile(int L0, int L1, int sizeX, int sizeY)
{


    // Check that the two dimensions are consecutive.
    // Tiling only applies on a consecutive band of loop dimensions.
    assert(L1 == L0 + 1);
    assert((sizeX > 0) && (sizeY > 0));
    assert(this->get_iteration_domain() != NULL);
    this->check_dimensions_validity({L0, L1});

    this->split(L0, sizeX);
    this->split(L1 + 1, sizeY);

    this->interchange(L0 + 1, L1 + 1);

    // DEBUG_INDENT(-4);
}

void computation::tile(int L0, int L1, int L2, int sizeX, int sizeY, int sizeZ)
{


    // Check that the two dimensions are consecutive.
    // Tiling only applies on a consecutive band of loop dimensions.
    assert(L1 == L0 + 1);
    assert(L2 == L1 + 1);
    assert((sizeX > 0) && (sizeY > 0) && (sizeZ > 0));
    assert(this->get_iteration_domain() != NULL);

    this->check_dimensions_validity({L0, L1, L2});

    //  Original loops
    //  L0
    //    L1
    //      L2

    this->split(L0, sizeX); // Split L0 into L0 and L0_prime
    // Compute the new L1 and the new L2 and the newly created L0 (called L0 prime)
    int L0_prime = L0 + 1;
    L1 = L1 + 1;
    L2 = L2 + 1;

    //  Loop after transformation
    //  L0
    //    L0_prime
    //      L1
    //        L2

    this->split(L1, sizeY);
    int L1_prime = L1 + 1;
    L2 = L2 + 1;

    //  Loop after transformation
    //  L0
    //    L0_prime
    //      L1
    //        L1_prime
    //          L2

    this->split(L2, sizeZ);

    //  Loop after transformation
    //  L0
    //    L0_prime
    //      L1
    //        L1_prime
    //          L2
    //            L2_prime

    this->interchange(L0_prime, L1);
    // Change the position of L0_prime to the new position
    int temp = L0_prime;
    L0_prime = L1;
    L1 = temp;

    //  Loop after transformation
    //  L0
    //    L1
    //      L0_prime
    //        L1_prime
    //          L2
    //            L2_prime

    this->interchange(L0_prime, L2);
    // Change the position of L0_prime to the new position
    temp = L0_prime;
    L0_prime = L2;
    L2 = temp;

    //  Loop after transformation
    //  L0
    //    L1
    //      L2
    //        L1_prime
    //          L0_prime
    //            L2_prime

    this->interchange(L1_prime, L0_prime);

    //  Loop after transformation
    //  L0
    //    L1
    //      L2
    //        L0_prime
    //          L1_prime
    //            L2_prime

    // DEBUG_INDENT(-4);
}

void computation::tile(std::string L0, std::string L1, int sizeX, int sizeY, std::string L0_outer, std::string L1_outer, std::string L0_inner, std::string L1_inner)
{


    assert(L0.length() > 0);
    assert(L1.length() > 0);
    assert(L0_outer.length() > 0);
    assert(L1_outer.length() > 0);
    assert(L0_inner.length() > 0);
    assert(L1_inner.length() > 0);

    std::vector<std::string> original_loop_level_names = this->get_loop_level_names();

    this->assert_names_not_assigned({L0_outer, L1_outer,
                                     L0_inner, L1_inner});

    std::vector<int> dimensions =
        this->get_loop_level_numbers_from_dimension_names({L0,
                                                           L1});
    assert(dimensions.size() == 2);

    // DEBUG(3, daisy::str_dump("The loop level that corresponds to " +
    //                             L0 + " is " + std::to_string(dimensions[0])));
    // DEBUG(3, daisy::str_dump("The loop level that corresponds to " +
    //                             L1 + " is " + std::to_string(dimensions[1])));

    this->tile(dimensions[0], dimensions[1], sizeX, sizeY);

    // Replace the original dimension name with new dimension names
    this->update_names(original_loop_level_names, {L0_outer, L1_outer, L0_inner, L1_inner}, dimensions[0], 2);

    // DEBUG_INDENT(-4);
}

bool computation::involved_subset_of_dependencies_is_legal(computation * second)
{



    assert(!this->get_name().empty());
    assert(this->get_function() != NULL);
    assert(!second->get_name().empty());
    
    assert(this->get_function()->dep_read_after_write != NULL);
    
    //extract this => Second dependencies from function

    isl_union_map * read_after_write_dep = isl_union_map_range_factor_domain(
        isl_union_map_copy(this->get_function()->dep_read_after_write));

    isl_union_map * write_after_read_dep = isl_union_map_range_factor_domain(
        isl_union_map_copy(this->get_function()->dep_write_after_read));

    isl_union_map * write_after_write_dep = isl_union_map_range_factor_domain(
        isl_union_map_copy(this->get_function()->dep_write_after_write));


    //construct the maps of the space needed
    std::string space_map = "{" + this->get_name() + "[";

    for(int i = 0 ; i< this->number_of_dims ; i++){

        space_map += "i" + std::to_string(i);

        if(i != (this->number_of_dims - 1))
        {
            space_map += ",";
        }
    }
    space_map += "]->" + second->get_name() + "[";

    for(int i = 0 ; i < second->number_of_dims ; i++){

        space_map+= "i" + std::to_string(i) ;

        if(i != (second->number_of_dims -1 ))
        {
            space_map += "," ;
        }
    }
    
    space_map += "]}" ;

    // DEBUG(3, daisy::str_dump(" the space map is to extract deps : "+space_map));

    isl_space * space = isl_map_get_space(
        isl_map_read_from_str(this->get_ctx(),space_map.c_str())
        );

    isl_map * my_map1 = isl_union_map_extract_map(read_after_write_dep,isl_space_copy(space));

    isl_map * my_map2 = isl_union_map_extract_map(write_after_read_dep,isl_space_copy(space));

    isl_map * my_map3 = isl_union_map_extract_map(write_after_write_dep,isl_space_copy(space));

    // DEBUG(10, daisy::str_dump(" the extracted deps are : "+std::string(isl_map_to_str(my_map1))));

    // DEBUG(10, daisy::str_dump(" the extracted deps are : "+std::string(isl_map_to_str(my_map2))));

    // DEBUG(10, daisy::str_dump(" the extracted deps are : "+std::string(isl_map_to_str(my_map3))));

    std::vector<isl_basic_map *> all_basic_maps ;
    

    auto f = [](isl_basic_map * bmap,void * user) { 

        std::vector<isl_basic_map *>& myName = *reinterpret_cast<std::vector<isl_basic_map*>*>(user);
     
        myName.push_back(bmap);
        return isl_stat_ok;
    };
    
    isl_stat (*fun_ptr)(isl_basic_map * p,void * m) = (f);

    isl_map_foreach_basic_map(my_map1,fun_ptr,(void * ) &all_basic_maps);

    isl_map_foreach_basic_map(my_map2,fun_ptr,(void * ) &all_basic_maps);

    isl_map_foreach_basic_map(my_map3,fun_ptr,(void * ) &all_basic_maps);


    /* ==========================================
        extract schedules of 2 computation
    */
   assert(this->get_function()->get_schedule()!= NULL) ;

   int m1 = isl_map_dim(this->get_schedule(), isl_dim_out);
   int m2 = isl_map_dim(second->get_schedule(), isl_dim_out);

   assert(m1 == m2) ;

    // DEBUG(3, daisy::str_dump(" the current schedule of computation "+this->get_name()+" : "+std::string(isl_map_to_str(this->get_schedule()))));
    // DEBUG(3, daisy::str_dump(" the current schedule of computation "+second->get_name()+" : "+std::string(isl_map_to_str(second->get_schedule()))));

    /* ==========================================
       making schedules comparable by mapping to the same time space
    */
  
    std::string empty = "" ;

    isl_map * this_schedule_unify = isl_map_set_tuple_name(
        isl_map_copy(this->schedule),
        isl_dim_out,
        empty.c_str()
        );


    isl_map * second_schedule_unify = isl_map_set_tuple_name(
        isl_map_copy(second->schedule),
        isl_dim_out,
        empty.c_str()
        );
    
    // DEBUG(3, daisy::str_dump(" first schedule adjusted into timestamp "+std::string(isl_map_to_str(this_schedule_unify))));
    // DEBUG(3, daisy::str_dump(" second schedule adjusted into timestamp "+std::string(isl_map_to_str(second_schedule_unify))));

    

    std::string s0_set = "[";

    for(int i=0 ;i<this->number_of_dims;i++)
    {
        s0_set += "n" + std::to_string(i) ;
        if(i != (this->number_of_dims -1 ))
        {
            s0_set += ",";
        }
    }
    s0_set += "]->{" + this->get_name() + "[";

    for(int i=0 ;i < this->number_of_dims; i++)
    {
        s0_set += "n" + std::to_string(i);
        if(i != (this->number_of_dims -1 ))
        {
            s0_set += "," ;
        }
    }
    s0_set += "]}" ;


    // DEBUG(3, daisy::str_dump(" initial set of first computation is : "+s0_set));

    isl_set * first_set = isl_set_read_from_str(this->get_ctx(),s0_set.c_str()) ;// in S0

    
   /* ==========================================
        always check that S0 is lex inferior than S1 , for that : S1>=S0 need always to be empty
    */
 
    bool overall_corectness = true ;

    // DEBUG(10, daisy::str_dump(" check the respect of previous deps nature start : "));
    
    for (auto& dependency:all_basic_maps)
    {

        // DEBUG(10, daisy::str_dump(" the dependency is : "+std::string(isl_basic_map_to_str(dependency))));

        isl_set * second_set = isl_set_apply(
            isl_set_copy(first_set),
            isl_map_from_basic_map(isl_basic_map_copy(dependency))
            );
        // in S1

        isl_set * time_first = isl_set_apply(isl_set_copy(first_set),isl_map_copy(this_schedule_unify));
        isl_set * time_second = isl_set_apply(isl_set_copy(second_set),isl_map_copy(second_schedule_unify));

        isl_map * result_sup = isl_set_lex_ge_set(
            time_first,
            time_second
        );

        if(isl_map_is_empty(result_sup) == false)
        {
            overall_corectness = false;
            // DEBUG(10, daisy::str_dump(" dependency is wrong by current schedule "));
            break;
        }
        else
        {
                // DEBUG(10, daisy::str_dump(" this dependency is respected by the current schedule  "));          
        }

        isl_set_free(second_set);
        isl_map_free(result_sup);
    }

    // DEBUG_INDENT(-4);

    isl_map_free(this_schedule_unify);
    isl_map_free(second_schedule_unify);

    isl_map_free(my_map1);
    isl_map_free(my_map2);
    isl_map_free(my_map3);

    isl_union_map_free(write_after_read_dep);
    isl_union_map_free(read_after_write_dep);
    isl_union_map_free(write_after_write_dep);

    return overall_corectness;
}

/**
  * Implementation internals.
  *
  * This function gets as input a loop level and translates it
  * automatically to the appropriate schedule dimension by:
  * (1) getting the dynamic schedule dimension that corresponds to
  * that loop level, then adding +1 which corresponds to the first
  * static dimension that comes after the dynamic dimension.
  *
  * Explanation of what static and dynamic dimensions are:
  * In the time-processor domain, dimensions can be either static
  * or dynamic.  Static dimensions are used to order statements
  * within a given loop level while dynamic dimensions represent
  * the actual loop levels.  For example, the computations c0 and
  * c1 in the following loop nest
  *
  * for (i=0; i<N: i++)
  *   for (j=0; j<N; j++)
  *   {
  *     c0;
  *     c1;
  *   }
  *
  * have the following representations in the iteration domain
  *
  * {c0(i,j): 0<=i<N and 0<=j<N}
  * {c1(i,j): 0<=i<N and 0<=j<N}
  *
  * and the following representation in the time-processor domain
  *
  * {c0[0,i,0,j,0]: 0<=i<N and 0<=j<N}
  * {c1[0,i,0,j,1]: 0<=i<N and 0<=j<N}
  *
  * The first dimension (dimension 0) in the time-processor
  * representation (the leftmost dimension) is a static dimension,
  * the second dimension (dimension 1) is a dynamic dimension that
  * represents the loop level i, ..., the forth dimension is a dynamic
  * dimension that represents the loop level j and the last dimension
  * (dimension 4) is a static dimension and allows the ordering of
  * c1 after c0 in the loop nest.
  *
  * \p dim has to be a static dimension, i.e. 0, 2, 4, 6, ...
  */
void computation::after_low_level(computation &comp, int level)
{


    // for loop level i return 2*i+1 which represents the
    // the static dimension just after the dynamic dimension that
    // represents the loop level i.
    int dim = loop_level_into_static_dimension(level);

    // DEBUG(3, daisy::str_dump("Setting the schedule of ");
    //       daisy::str_dump(this->get_name());
    //       daisy::str_dump(" after ");
    //       daisy::str_dump(comp.get_name());
    //       daisy::str_dump(" at dimension ");
    //       daisy::str_dump(std::to_string(dim)));
    // DEBUG(3, daisy::str_dump("Setting the schedule of ");
    //       daisy::str_dump(this->get_name());
    //       daisy::str_dump(" to be equal to the schedule of ");
    //       daisy::str_dump(comp.get_name());
    //       daisy::str_dump(" at all the dimensions before dimension ");
    //       daisy::str_dump(std::to_string(dim)));

    comp.get_function()->align_schedules();

    // DEBUG(3, daisy::str_dump("Preparing to adjust the schedule of the computation ");
    //       daisy::str_dump(this->get_name()));
    // DEBUG(3, daisy::str_dump("Original schedule: ", isl_map_to_str(this->get_schedule())));

    assert(this->get_schedule() != NULL);
    // DEBUG(3, daisy::str_dump("Dimension level in which ordering dimensions will be inserted : ");
    //       daisy::str_dump(std::to_string(dim)));
    assert(dim < (signed int) isl_map_dim(this->get_schedule(), isl_dim_out));
    assert(dim >= computation::root_dimension);

    isl_map *new_sched = NULL;
    for (int i = 1; i<=dim; i=i+2)
    {
        if (i < dim)
        {
            // Get the constant in comp, add +1 to it and set it to sched1
            int order = isl_map_get_static_dim(comp.get_schedule(), i);
            new_sched = isl_map_copy(this->get_schedule());
            new_sched = add_eq_to_schedule_map(i, 0, -1, order, new_sched);
        }
        else // (i == dim)
        {
            // Get the constant in comp, add +1 to it and set it to sched1
            int order = isl_map_get_static_dim(comp.get_schedule(), i);
            new_sched = isl_map_copy(this->get_schedule());
            new_sched = add_eq_to_schedule_map(i, 0, -1, order + 10, new_sched);
        }
        this->set_schedule(new_sched);
    }

    // DEBUG(3, daisy::str_dump("Schedule adjusted: ",
    //                             isl_map_to_str(this->get_schedule())));

    // DEBUG_INDENT(-4);
}

void computation::tile(std::string L0, std::string L1, std::string L2, int sizeX, int sizeY, int sizeZ, std::string L0_outer, std::string L1_outer, std::string L2_outer, std::string L0_inner, std::string L1_inner, std::string L2_inner)
{
    assert(L0.length() > 0);
    assert(L1.length() > 0);
    assert(L2.length() > 0);
    assert(L0_outer.length() > 0);
    assert(L1_outer.length() > 0);
    assert(L2_outer.length() > 0);
    assert(L0_inner.length() > 0);
    assert(L1_inner.length() > 0);
    assert(L2_inner.length() > 0);

    this->assert_names_not_assigned({L0_outer, L1_outer,
                                     L2_outer, L0_inner,
                                     L1_inner, L2_inner});

    std::vector<std::string> original_loop_level_names = this->get_loop_level_names();

    std::vector<int> dimensions =
        this->get_loop_level_numbers_from_dimension_names({L0,
                                                           L1,
                                                           L2});
    assert(dimensions.size() == 3);

    // DEBUG(3, daisy::str_dump("The loop level that corresponds to " +
    //                             L0 + " is " + std::to_string(dimensions[0])));
    // DEBUG(3, daisy::str_dump("The loop level that corresponds to " +
    //                             L1 + " is " + std::to_string(dimensions[1])));
    // DEBUG(3, daisy::str_dump("The loop level that corresponds to " +
    //                             L2 + " is " + std::to_string(dimensions[2])));

    this->tile(dimensions[0], dimensions[1], dimensions[2],
               sizeX, sizeY, sizeZ);

    this->update_names(original_loop_level_names, {L0_outer, L1_outer, L2_outer,
                                                   L0_inner, L1_inner, L2_inner}, dimensions[0], 3);

    // DEBUG_INDENT(-4);
}

dace_block::dace_block(const std::vector<computation*> children) : children(children)
{
    // dace_block is a special child of computation. Don't call parent constructor.
}

void dace_block::interchange(int L0, int L1) {
    for (auto &child : this->children) {
        child->interchange(L0, L1);
    }
}

void dace_block::tile(int L0, int L1, int sizeX, int sizeY) {
    for (auto &child : this->children) {
        child->tile(L0, L1, sizeX, sizeY);
    }
}

void dace_block::tile(int L0, int L1, int L2, int sizeX, int sizeY, int sizeZ) {
    for (auto &child : this->children) {
        child->tile(L0, L1, L2, sizeX, sizeY, sizeZ);
    }
}

}
}
