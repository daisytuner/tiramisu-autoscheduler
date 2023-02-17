#pragma once

#include <string>
#include <vector>
#include <map>
#include <unordered_map>
#include <unordered_set>
#include <list>
#include <memory>

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

#include "core.h"
#include "buffer.h"

namespace daisy {

namespace tiramisu {

class computation_info;
class syntax_tree;
class dnn_accesses;

class function;

typedef std::unordered_map<std::string, std::vector<std::vector<int>>> buffer_accesses;

struct dace_desc
{
  /**
   * Some metrics about the computation.
   */
  int nb_additions;
  int nb_substractions;
  int nb_multiplications;
  int nb_divisions;
  buffer_accesses accesses;
};


class computation
{

  friend dnn_accesses;
  friend computation_info;
  friend syntax_tree;
  friend function;

private:

    /**
      * An ISL context associate with the function.
      */
    isl_ctx *ctx;

    /**
      * The name of this computation.
      * Computation names should not start with _ (an underscore).
      * Names starting with _ are reserved names.
      */
    std::string name;

    /* The number of dimensions in the original definition of the computation. */
    int number_of_dims;

    /**
      * Data type: the type of the value returned by the computation.
      */
    primitive_t data_type;

    /** 
     * Information about the computation
     */
    dace_desc desc;

    /**
      * Access function.  A map indicating how each computation should be stored
      * in memory.  It indicates in which buffer the computation should be stored
      * and which element of the buffer exactly it should be stored.
      */
    isl_map *access;

    /**
      * A vector that contains the list of let statements associated
      * with this computation.
      *
      * A let statement that is associated with the computation is a let statement
      * that will be added just before the computation.  The scope of the variable
      * defined by the let statement is this computation alone. i.e., it is not
      * defined in other computations.
      */
    //std::vector<std::pair<std::string, daisy::expr>> associated_let_stmts;

    /**
     * The buffer attached "automatically" to this computation.
     * If the buffer is not created automatically, this variable will be empty.
     */
    buffer* automatically_allocated_buffer;

    /**
      * Number of definitions added to this computation. Each time the function
      * add_definitions is called, definitions_number is incremented.
      */
    int definitions_number;

    /**
      * The ID of this definition. Each new computation when created has a
      * definition_ID equal to 0.  When a new definition is added, its ID
      * is 1, when a new definition is added, its definition ID is set to 2,
      * ...
      */
    int definition_ID;

    /**
     * An integer indicating the number of duplicates of this computation.
     * We use this to set the duplicate ID of any newly created computation.
     * Whenever a new duplicate is create, this number is incremented.
     */
    int duplicate_number;

    /**
      * If has_multiple_definitions() is true, then this variable contains the
      * computation that was defined first among all the multiple definitions.
      */
    computation *first_definition;

    /**
      * The function where this computation is declared.
      */
    function *fct;

    /**
      * An isl_ast_expr representing the index of the array where the computation
      * will be stored.  This index is computed after the scheduling is done.
      */
    std::vector<isl_ast_expr *> index_expr;

    /**
     * A map between the original names of the iterators of a computation
     * and their transformed form after schedule (also after renaming).
     *
     * If in the original computation, we had
     *
     * {C[i0, i1]: ...}
     *
     * And if in the generated code, the iterators are called c0, c1, c2 and c3 and
     * the loops are tiled, then the map will be
     *
     * {<i0, c0*10+c2>, <i1, c1*10+c3>}.
     */
    std::map<std::string, isl_ast_expr *> iterators_map;

    /**
      * Does this computation represent a let statement ?
      *
      * Let statements should be treated differently:
      * - During Halide code generation a Halide let statement should be
      * created instead of an assignment statement.
      * - A let statement does not have/need an access function because
      * it writes directly to a scalar.
      * - When targeting Halide, let statements should be created after
      * their body is created, because the body is an argument needed
      * for the creation of the let statement.
      */
    bool is_let;

    /**
      * This is variable is set to true if this computation is the first definition.
      * It is set to false if has_multiple_definitions() is true but this computation
      * is not the first one defined.
      * This is useful because in tiramisu, we assume that all the computations that
      * have the same name must have the same memory access. Thus any new definition
      * of a computation must copy the same memory access as the first definition, thus
      * we need to know which computation is the first definition.
      */
    bool is_first;

    /* Is true if the the computation is inline. */
    bool is_inline;

    /**
      * Iteration domain of the computation.
      * In this representation, the order of execution of computations
      * is not specified, the computations are also not mapped to memory.
     */
    isl_set *iteration_domain;

    /**
      * A vector describing the access variables in the original definition of  a computation.
      * For every named dimension, a pair representing the index of the named dimension
      * and the name of the dimension is added to access_variables.
      * E.g. if a computation is defined as S[i, 0, j], access_variables will contain
      * {(0, "i"), (2, "j")}.
      */
    std::vector<std::pair<int, std::string>> access_variables;

    /**
     * A predicate around the computation. The computation is executed
     * only if this predicate is true. This is useful to insert a non-affine
     * condition around the computation.
     */
    //daisy::expr predicate;

    /**
      * The schedules of the computation.
      */
    isl_map * schedule;

    /**
      * TODO: use buffers directly from computations, no need to have
      * bindings.
      *
      * \p schedule_this_computation should be set to true when the computation
      * should be scheduled and when code for the computation should be generated
      * during code generation.
      * It should be set to false when the computation is used to represent a
      * buffer (i.e., the computation is used only as a binding to a buffer).
      * In this case, the computation is not scheduled and no code for the
      * computation is generated.
      */
    bool schedule_this_computation;

    // Private class members that are computed during code generation.

    /**
      * Time-processor domain of the computation.
      * In this representation, the logical time of execution and the
      * processor where the computation will be executed are both
      * specified.
      */
    isl_set *time_processor_domain;


    /**
      * Contains a list of all definitions added to this computation. The 0th definition is
      * always this very computation.
      */
    std::vector<computation*> updates;

    void init_computation(std::string iteration_space_str, std::string access);

protected:
  computation();

public:

    /**
      * \brief Constructor for computations.
      *
      * \details \p iteration_domain is a string that represents the iteration
      * domain of the computation.  The iteration domain should be written
      * in the ISL format (http://barvinok.gforge.inria.fr/barvinok.pdf Section 1.2.1).
      *
      * The iteration domain of a statement is a set that contains
      * all of the execution instances of the statement (a statement in a
      * loop has an execution instance for each loop iteration in which
      * it executes). Each execution instance of a statement in a loop
      * nest is uniquely represented by an identifier and a tuple of
      * integers  (typically, the values of the outer loop iterators).
      *
      * For example, the iteration space of the statement S0 in the following
      * loop nest
      *
      * \code
      * for (i=0; i<2; i++)
      *   for (j=0; j<3; j++)
      *      S0;
      * \endcode
      *
      * is {S0[0,0], S0[0,1], S0[0,2], S0[1,0], S0[1,1], S0[1,2]}
      *
      * S0[0,0] is the execution instance of S0 in the iteration [0,0].
      *
      * The previous set of integer tuples can be compactly described
      * by affine constraints as follows
      *
      * {S0[i,j]: 0<=i<2 and 0<=j<3}
      *
      * In general, the loop nest
      *
      * \code
      * for (i=0; i<N; i++)
      *   for (j=0; j<M; j++)
      *      S0;
      * \endcode
      *
      * has the following iteration domain
      *
      * {S0[i,j]: 0<=i<N and 0<=j<M}
      *
      * This should be read as: the set of points [i,j] such that
      * 0<=i<N and 0<=j<M.
      *
      * The name of the computation in the iteration domain should not
      * start with _ (an underscore).  Names starting with _ are reserved
      * names.
      *
      * \p data_type is the type of the computation, i.e. the type of the expression
      * computed by the computation. Example of types include (p_uint8,
      * p_uint16, p_uint32, ...).
      *
      * \p fct is a pointer to the Tiramisu function where this computation
      * should be added.
      *
      * Bound Inference:
      * The user can declare computations without providing any constraint
      * about the iteration domain, in this case he can rely on bound inference
      * to infer the constraints about each iteration domain. The user needs only
      * to provide constraints over the domains of the last computations (last
      * consumers), and Tiramisu will propagate these constraints to all the
      * chain of computations that precede those consumers.
      * Note that bound inference is not possible if you have multiple definitions
      * of the same computation. In such a case, you should provide constraints
      * over the iteration domain when you declare the computation.
      *
      * Examples about bound inference are provided in test_22 to test_25.
      */
    computation(
      function *fct,
      std::string iteration_domain,
      std::string access,
      primitive_t data_type,
      dace_desc desc
    );

    /**
     * Return the name of the computation.
     */
    const std::string& get_name() const;

    /**
      * Return the function where the computation is declared.
      */
    function* get_function() const;

    /**
      * Get the data type of the computation.
      */
    primitive_t get_data_type() const;

    /**
     * Return the iteration domain of the computation.
     * In this representation, the order of execution of computations
     * is not specified, the computations are also not mapped to memory.
     */
    isl_set* get_iteration_domain() const;
    void set_iteration_domain(isl_set* domain);

    /**
     * Get the schedule of the computation.
     */
    isl_map *get_schedule();

    /**
      * Return the context of the computations.
      */
    isl_ctx *get_ctx() const;

    /**
      * Return the names of iteration domain dimensions.
      */
    std::vector<std::string> get_iteration_domain_dimension_names();

    /**
      * Return the access function of the computation.
      */
    isl_map *get_access_relation() const;

    /**
     * Set the access function of the computation.
     *
     * The access function is a relation from computations to buffer locations.
     * \p access_str is a string that represents the relation (in ISL format,
     * http://isl.gforge.inria.fr/user.html#Sets-and-Relations).
     */
    void set_access(std::string access_str);
    void set_access(isl_map *access);

    std::shared_ptr<buffer> get_buffer() const;

    /**
      * Get the number of dimensions of the iteration
      * domain of the computation.
      */
    int get_iteration_domain_dimensions_number();

    void name_unnamed_iteration_domain_dimensions();

    /**
      * Get the number of loop levels.
      */
    int get_loop_levels_number();

    /**
      * Search the time-space domain (the range of the schedule) and
      * return the loop level numbers that correspond to the dimensions
      * named \p dim.
      * In other words, translate the vector of dimension names (\p dim_names)
      * into loop level numbers. We need to do this because the internal Tiramisu
      * scheduling functions use dimension numbers instead of dimension
      * names (which are used in the user level scheduling functions).
      */
    std::vector<int> get_loop_level_numbers_from_dimension_names(std::vector<std::string> dim_names);

    /**
      * Return the names of loop levels.
      * (i.e., the names of dynamic dimensions in time-space).
      */
    std::vector<std::string> get_loop_level_names();

    std::string get_dimension_name_for_loop_level(int loop_level);


    void check_dimensions_validity(std::vector<int> dimensions);

    /**
      * Return true if the this computation is supposed to be scheduled
      * by Tiramisu.
      */
    bool should_schedule_this_computation() const;

    int compute_maximal_AST_depth();


    /**
      * Set the schedule indicated by \p map.
      *
      * \p map is a string that represents a mapping from the iteration domain
      * to the time-space domain (the ISL format to represent maps is
      * documented in http://barvinok.gforge.inria.fr/barvinok.pdf in Sec 1.2.2).
      *
      * The schedule is a map from the iteration domain to a time space
      * domain. The same name of space should be used for both the range
      * and the domain of the schedule.
      *
      * In general, users can set the schedule using high level functions such
      * as before(), after(), tile(), compute_at(), vectorize(), split(), ...
      * The use of this function is only reserved for advanced users who want
      * a low level control of the schedule.
      *
      * Vectors in the time-space domain have the following form
      *
      * computation_name[redundancy_ID,static,dynamic,static,dynamic,static,dynamic,static,...]
      *
      * The first dimension of the vector is used to indicate the redundancy ID
      * (the notion of the redundancy ID is explained later).
      *
      * The following dimensions are interleaved dimensions: static, dynamic, static,
      * dynamic, ...
      * Dynamic dimensions represent the loop levels, while static dimensions are
      * used to order statements within a given block of statements in a given loop
      * level.
      * For example, the computations c0 and c1 in the following loop nest
      *
      * \code
      * for (i=0; i<N: i++)
      *   for (j=0; j<N; j++)
      *   {
      *     c0;
      *     c1;
      *   }
      *
      * \endcode
      *
      * have the following representations in the iteration domain
      *
      * \code
      * {c0[i,j]: 0<=i<N and 0<=j<N}
      * {c1[i,j]: 0<=i<N and 0<=j<N}
      * \endcode
      *
      * and the following representation in the time-space domain
      *
      * \code
      * {c0[0,0,i,0,j,0]: 0<=i<N and 0<=j<N}
      * {c1[0,0,i,0,j,1]: 0<=i<N and 0<=j<N}
      * \endcode
      *
      * The first dimension (dimension 0) in the time-space
      * domain (the leftmost dimension) is the redundancy ID
      * (in this case it is 0, the meaning of this ID will be explained later).
      * The second dimension (starting from the left) is a static dimension,
      * the third dimension is a dynamic dimension that
      * represents the loop level i, ..., the fifth dimension is a dynamic
      * dimension that represents the loop level j and the last dimension
      * (dimension 5) is a static dimension and allows the ordering of
      * c1 after c0 in the loop nest.
      *
      * To transform the previous iteration domain to the
      * time-space domain, the following schedule should be used:
      *
      * \code
      * {c0[i,j]->c0[0,0,i,0,j,0]: 0<=i<N and 0<=j<N}
      * {c1[i,j]->c1[0,0,i,0,j,1]: 0<=i<N and 0<=j<N}
      * \endcode
      *
      * The first dimension called "redundancy ID" is only meaningful if the
      * computation is redundant. i.e., some parts of the computation are
      * redundantly computed.  Redundant computations are in general used to
      * maximize parallelism and data locality on the expense of doing some
      * computations redundantly.
      * For example, if the two computations c1(i,j) and c2(i,j) both depend
      * on the computation c0(i,j), instead of waiting for c0(i,j) to be
      * computed and then computing c1(i,j) and c2(i,j) in parallel, the thread
      * executing c1(i,j) can compute c0(i,j) by itself and then run c1(i,j).
      * The thread that computes c2(i,j) can do the same and compute c0(i,j)
      * by itself and then compute c2(i,j). In this case the two threads do not
      * need to wait. This is done at the expense of redundant computation since
      * c0(i,j) is computed by both threads.
      *
      * In general redundant computations are useful when tiling stencil
      * computations.  In the context of stencils such a tiling is called
      * "overlapped tiling".  Tiles that depend on results computed by other
      * tiles that run in parallel can compute the boundaries redundantly which
      * allows them to avoid waiting and thus can run in parallel.
      *
      * In Tiramisu, the user can indicate that a chunk of a computation
      * should be computed redundantly. The original computation always has a redundancy
      * ID equal to 0 (which means this is the original computation).
      * The redundant chunk has an ID that is different from 0 and that is
      * used to uniquely identify it.
      *
      * For example if we want to compute all of c0 three times (that is,
      * compute the original computation and compute two redundant computations),
      * we can use the following schedules:
      *
      * The schedule of the original computation:      {c0[i,j]->c0[0,0,i,0,j,0]: 0<=i<N and 0<=j<N}
      * The schedule of the redundant computation N°1: {c0[i,j]->c0[1,0,i,0,j,0]: 0<=i<N and 0<=j<N}
      * The schedule of the redundant computation N°2: {c0[i,j]->c0[2,0,i,0,j,0]: 0<=i<N and 0<=j<N}
      *
      * The schedule of c0 in this case would be three maps that map c0[i,j] to
      * the three different redundant computations in the time-processor domain:
      *
      * \code
      * {c0[i,j]->c0[0,0,i,0,j,0]: 0<=i<N and 0<=j<N;
      *  c0[i,j]->c0[1,0,i,0,j,0]: 0<=i<N and 0<=j<N;
      *  c0[i,j]->c0[2,0,i,0,j,0]: 0<=i<N and 0<=j<N}
      * \endcode
      *
      * The function set_schedule() overrides any other schedule set by the high level
      * scheduling functions.  Currently the user has to choose between using the high
      * level scheduling functions or using this low level set_schedule function. The user
      * cannot mix the use of the two in the same program because they are not compatible.
      */
    // @{
    void set_schedule(isl_map *map);

    /**
      * Set an identity schedule for the computation.
      *
      * This identity schedule is an identity relation created from the
      * iteration domain.
      * This sets the schedule of the original computation
      * and does not set the schedule of any duplicate
      * computation.
      */
    void set_identity_schedule_based_on_iteration_domain();

    isl_map* gen_identity_schedule_for_iteration_domain();

    void add_schedule_constraint(std::string domain_constraints, std::string range_constraints);

    void gen_time_space_domain();
    
    isl_set* get_time_processor_domain() const;

    isl_set* get_trimmed_time_processor_domain();

    isl_set* intersect_set_with_context(isl_set *set);

    /**
      * Assign a name to iteration domain dimensions that do not have a name.
      */
    void name_unnamed_time_space_dimensions();

    void set_schedule_domain_dim_names(std::vector<int> loop_levels, std::vector<std::string> names);

    /**
     * Checks the correctness of a subset of dependencies after applying changes on the schedules (e.g., tiling, skewing, and shifting).
     * The checked subset of dependencies is the set of dependencies mapping from this computation (this) to second computation (second).
     * This methods returns a boolean: True if this subset of dependencies is respected, otherwise False.
     * It relies fully on the dependence analysis result, so the  method \p perform_full_dependency_analysis() must be invoked before.
     * To correctly invoke this method : schedules must be aligned (same out dimension size) and ordered,
     * so invoking \p prepare_schedules_for_legality_checks() method before is mandatory. 
    */
    bool involved_subset_of_dependencies_is_legal(computation * second) ;

    void set_loop_level_names(std::vector<std::string> names);
    void set_loop_level_names(std::vector<int> loop_levels, std::vector<std::string> names);

    void assert_names_not_assigned(std::vector<std::string> dimensions);

    /**
      * Update loop level names. This function should be called after each scheduling operation
      * because scheduling usually changes the original loop level names.
      * This function erases \p nb_loop_levels_to_erase loop level names starting from the
      * loop level \p start_erasing. It then inserts the loop level names \p new_names in
      * \p start_erasing. In other words, it replaces the names of loop levels from
      * \p start_erasing to \p start_erasing + \p nb_loop_levels_to_erase with the loop levels
      * indicated by \p new_names.  This function sets the non erased loop levels to be equal to the
      * original loop level names.
      *
      * \p original_loop_level_names : a vector containing the original loop level names (loop level
      * names before scheduling).
      *
      * \p new_names : the new loop level names.
      *
      * \p start_erasing : start erasing loop levels from this loop level.
      *
      * \p nb_loop_levels_to_erase : number of loop levels to erase.
      *
      * Example. Assuming the original loop levels are {i0, i1, i2, i3}
      *
      * Calling this->update_names({i0, i1, i2, i3}, {x0, x1}, 1, 2) updates the loop levels to become
      * {i0, x0, x1, i3}.
      */
    void update_names(std::vector<std::string> original_loop_level_names, std::vector<std::string> new_names, int start_erasing, int nb_loop_levels_to_erase);


    computation& get_last_update();

    std::vector<computation*>& get_updates();

    computation& get_update(int i);

    /**
      * Schedule this computation to run after the computation \p comp.
      * The computations are placed after each other in the loop level \p level.
      * The outermost loop level is 0.  The root level is computation::root_dimension.
      *
      * For example assuming we have the two computations
      *
      *     {S0[i,j]: 0<=i<N and 0<=j<N} and {S1[i,j]: 0<=i<N and 0<=j<N}
      *
      * In order to make S1 run after S0 in the i loop, one should use
      *
      *     S1.after_low_level(S0,0)
      *
      * which means: S1 is after S0 at the loop level 0 (which is i).
      *
      * The corresponding code is
      *
      * \code
      *     for (i=0; i<N; i++)
      *     {
      *         for (j=0; j<N; j++)
      *             S0;
      *         for (j=0; j<N; j++)
      *             S1;
      *     }
      * \endcode
      *
      * S1.after_low_level(S0,1)
      *
      * means: S1 is after S0 at the loop level 1 (which is j) and would yield
      * the following code
      *
      * \code
      * for (i=0; i<N; i++)
      *   for (j=0; j<N; j++)
      *   {
      *     S0;
      *     S1;
      *   }
      * \endcode
      *
      * S1.after_low_level(S0, computation::root_dimension)
      * means S1 is after S0 at the main program level and would yield
      * the following code
      *
      * \code
      * for (i=0; i<N; i++)
      *   for (j=0; j<N; j++)
      *     S0;
      * for (i=0; i<N; i++)
      *   for (j=0; j<N; j++)
      *     S1;
      * \endcode
      *
      * To specify that this computation is after \p comp in multiple levels,
      * the user can provide those levels in the \p levels vector.
      *
      * S1.after_low_level(S0, {0,1})
      *
      * means that S1 is after S0 in the loop level 0 and in the loop level 1.
      *
      * Note that
      *
      * S1.after_low_level(S0, L)
      *
      * would mean that S1 and S0 share the same loop nests for all the loop
      * levels that are before L and that S1 is after S0 in L only.  S1 is not
      * after S0 in the loop levels that are before L.
      *
      */
    // @{
    void after_low_level(computation &comp, int level);

    /**
      * Tag the loop level \p L to be parallelized.
      */
    void tag_parallel_level(int par_dim);

    /**
      * Interchange (swap) the two loop levels \p L0 and \p L1.
      */
    virtual void interchange(int L0, int L1);
    void interchange(std::string L0, std::string L1);

    /**
      * Split the loop level \p L0 of the iteration space into two
      * new loop levels.
      *
      * \p sizeX is the extent (size) of the inner loop created after
      * splitting.
      */
    //@{
    void split(int L0, int sizeX);
    //@}

    void split_with_lower_bound(int L0, int sizeX, std::string lower_bound);


    void separate(int dim, std::string N, int v, std::string lower_bound);

    void add_definitions(std::string iteration_domain_str, bool schedule_this_computation, daisy::primitive_t t, function *fct);

    void after(computation &comp, int level);


    /**
      * Separate then split the loop \p L0
      * (see daisy::computation::separate and daisy::computation::split)
      *
      * This function returns true if the split actually happened (in some
      * cases, if the loop cannot be split because it is too small for example,
      * then no split happens even after calling the split function, in such
      * cases the function returns false).
      */
    //@{
    bool separateAndSplit(int L0, int sizeX);
    //@}

    /**
      * Tile the two loop levels \p L0 and \p L1 with rectangular
      * tiling. \p sizeX and \p sizeY represent the tile size.
      * \p L0 and \p L1 should be two consecutive loop levels
      * (i.e., \p L0 = \p L1 + 1) and they should satisfy
      * \p L0 > \p L1.
      */
    // @{
    virtual void tile(int L0, int L1, int sizeX, int sizeY);
    virtual void tile(int L0, int L1, int L2, int sizeX, int sizeY, int sizeZ);
    void tile(std::string L0, std::string L1, int sizeX, int sizeY, std::string L0_outer, std::string L1_outer, std::string L0_inner, std::string L1_inner);
    void tile(std::string L0, std::string L1, std::string L2, int sizeX, int sizeY, int sizeZ, std::string L0_outer, std::string L1_outer, std::string L2_outer, std::string L0_inner, std::string L1_inner, std::string L2_inner);
    // @}

    /**
      * Equivalent of computation::root but to be used with scheduling
      * functions that take loop level (integers) as input instead of
      * daisy::var.
      */
    const static int root_dimension = -1;

};

/**
  * A class that represents a group of computations that share schedules up to
  * a level. This class is used for convenience when a set of scheduling
  * commands should be applied to several computations (i.e GPU kernels).
  */
class dace_block: public computation {
private:

    /**
      * List of computations that this block contains. Actual order of the
      * computations is independent of the vector order and determined by the
      * scheduling commands.
      */
    const std::vector<computation *> children;

public:

    /**
      * \brief Constructor for block.
      *
      * \details This constructor creates a block of computations.
      *
      * In Tiramisu each computation should be scheduled separately. If the
      * user has multiple computations and wants to tile them for example,
      * the user has to apply the tiling command on each one of the computations
      * separately which is not convenient. The block class solves this problem.
      * It allows the user to create a block of computations. Any scheduling
      * command applied to the block is automatically applied to all of its
      * computations.
      *
      * \p children is the list of the computations of the block.
      *
      * The actual order of the computations is not determined by the vector
      * order. It is rather determined by the scheduling commands.
      */
    dace_block(const std::vector<computation*> children);

    /**
      * Overriden scheduling methods from computation. These transformations
      * will be applied to all children of the block.
      */
    // @{
    void interchange(int L0, int L1) override;
    //void skew(int i, int j, int a, int b) override;
    //void skew(int i, int j, int a, int b, int c, int d) override;
    void tile(int L0, int L1, int sizeX, int sizeY) override;
    void tile(int L0, int L1, int L2, int sizeX, int sizeY, int sizeZ) override;
    // @}
};  // class block

}

}

