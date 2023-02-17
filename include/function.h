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

class syntax_tree;

class computation;


class function
{
  friend syntax_tree;
  friend computation;

private:
    /**
      * An ISL context associate with the function.
      */
    isl_ctx *ctx;

    /**
      * The name of the function.
      * Function names should not start with _ (an underscore).
      * Names starting with _ are reserved names.
      */
    std::string name;

    /**
      * Body of the function (a vector of computations).
      * The order of the computations in the vector does not have any
      * effect on the actual order of execution of the computations.
      * The order of execution of computations is specified through the
      * schedule.
      */
    std::vector<computation*> body;

    /**
      * A map representing the buffers of the function. Some of these
      * buffers are passed to the function as arguments and some are
      * declared and allocated within the function itself.
      */
    std::map<std::string, std::shared_ptr<buffer>> buffers_list;

    /**
     * The context set of the function, i.e. a set representing the
     * constraints over the parameters.
     * The parameters of a function are the function invariants (constants).
     */
    isl_set *context_set;

    /**
      * The set of all computations that have no computation scheduled before them.
      * Does not include allocation computations created using
      * allocate_and_map_buffer_automatically().
      */
    std::unordered_set<computation*> starting_computations;

    /**
      * isl_union_map that describes the true dependencies considering buffers access (the read after write dependencies),
      * the relation describes for each read access the last previous write access that writes into the same buffer element
      * all isl_maps in the union_isl_map are represented in the format:
      * last_write -> [read_access -> used_buffer]
      * Isl methods (isl_map_range_factor_domain and isl_map_range_factor_range) could be used to the extract a more simple relation.
      * Useful to check correctness with legality checks.
      */

    isl_union_map * dep_read_after_write ;

    /**
      * isl_union_map that describes the false dependencies considering buffers access (the write after write dependencies),
      * the relation describes for each write access the last previous write access that write into the same buffer element
      * all isl_maps in the union_isl_map are represented in the format:
      * last_write -> [write_access -> used_buffer]
      * Isl methods (isl_map_range_factor_domain and isl_map_range_factor_range) could be used to the extract a more simple relation.
      * Useful to check correctness with legality checks.
      */

    isl_union_map * dep_write_after_write ; 

    
    /**
      * isl_union_map that describes the false dependencies considering buffers access (the write after read dependencies),
      * the relation describes for each write access all the previous read access that used a the previous value [a buffer element] before current write 
      * all isl_maps in the union_isl_map are represented in the format:
      * read_access_with_the_previous_value-> [ write_access -> used_buffer]
      * Isl methods (isl_map_range_factor_domain and isl_map_range_factor_range) could be used to the extract a more simple relation.
      * Useful to check correctness with legality checks.
      */

    isl_union_map * dep_write_after_read ;

    /**
     *  A union map that describes the live-in access (e.g., compulation[i,j]-> buffer1[i,j]).
     *  Read accesses that have no write before them (i.e., the value is external).
    */
    isl_union_map * live_in_access ;

    /**
    *  A union map that describes the live-out access (e.g., compulation[i,j]-> buffer1[i,j]).
    *  Write access that do not have another write after them (i.e., last written value into the buffer).
    */
    isl_union_map * live_out_access ;

    /**
      * Keeps track of allocation computations created using
      * allocate_and_map_buffer_automatically() to schedule them during gen_ordering_schedules.
      */
     std::vector<computation*> automatically_allocated;

    /**
      * A vector representing the parallel dimensions around
      * the computations of the function.
      * A parallel dimension is identified using the pair
      * <computation_name, level>, for example the pair
      * <S0, 0> indicates that the loop with level 0
      * (i.e. the outermost loop) around the computation S0
      * should be parallelized.
      */
    std::vector<std::pair<std::string, int>> parallel_dimensions;

    /**
      * Stores all high level scheduling instructions between computations; i.e. if a user calls
      * for example c2.after(c1, L), sched_graph[&c1] would contain the key &c2, and
      * sched_graph[&c1][&c2] = L.
      * At the end of scheduling, the graph should respect the following rules:
      *     - There should be exactly one computation with no computation scheduled before it.
      *     - Each other computation should have exactly one computation scheduled before it.
      * In other words, the graph should be a valid tree.
      * Does not include allocation computations created using
      * allocate_and_map_buffer_automatically().
      */
    std::unordered_map<computation*, std::unordered_map<computation*, int>> sched_graph;

    /**
      * Same as sched_graph, except in reverse order (from after to before).
      */
    std::unordered_map<computation*, std::unordered_map<computation*, int>> sched_graph_reversed;

    /**
     * A boolean set to true if low level scheduling was used in the program.
     * If it is used, then high level scheduling commands such as .before(),
     * .after(), ...
     */
    bool use_low_level_scheduling_commands;

public:

    function(std::string name);

    const std::string &get_name() const;

    isl_set* get_program_context() const;

    isl_ctx *get_isl_ctx() const;

    const std::map<std::string, std::shared_ptr<buffer>> get_buffers() const;
    
    /**
      * Add a buffer to the function.
      * The buffers of the function are either:
      * - buffers passed to the function as arguments, or
      * - buffers that are declared and allocated within the function
      * itself.
      * The first element of the pair is the name of the buffer (it is
      * used as a key), the second element of the pair is a pointer
      * to the buffer.
      */
    void add_buffer(std::pair<std::string, std::shared_ptr<buffer>> buf);

    const std::vector<computation*> get_computations() const;

    /**
      * Add a computation to the function.  The order in which
      * computations are added to the function is not important.
      * The order of execution is specified using the schedule.
      * This doesn't allow computations with duplicate names.
      */
    void add_computation(computation *cpt);

    std::vector<computation*> get_computation_by_name(std::string str) const;

    /**
      * Return the union of all the schedules
      * of the computations of the function.
      */
    isl_union_map *get_schedule() const;

    /**
      * Return the union of all the iteration domains
      * of the computations of the function.
      */
    isl_union_set *get_iteration_domain() const;

    void add_parallel_dimension(std::string stmt_name, int vec_dim);

    void reset_schedules();

    void remove_dimension_tags();

   /**
     * Return true if the usage of high level scheduling comments is valid; i.e. if
     * the scheduling relations formed using before, after, compute_at, etc.. form a tree.
     *
     * More specifically, it verifies that:
     *     - There should be exactly one computation with no computation scheduled before it.
     *     - Each other computation should have exactly one computation scheduled before it.
     */
    bool is_sched_graph_tree();
    bool is_sched_graph_tree_dfs(computation * comp,std::unordered_set<computation *> &visited);
  
    void clear_sched_graph();


  /**
       * Checks if the given fused computations could legally have their loop level \p i as parallel using dependence analysis and legality check.
       * It relies fully on the dependence analysis result, so the  method \p perform_full_dependency_analysis() must be invoked before.
       * To correctly invoke this method : schedules must be aligned (same out dimension size) and ordered,
       * so invoking \p prepare_schedules_for_legality_checks() method before is mandatory. 
    */
    bool loop_parallelization_is_legal(std::string i, std::vector<computation*> fused_computations);

    /**
     * Full check of schedule legality for this function.
     */
    bool check_legality_for_function();

    /**
     * Prepare the schedules of the computations for legality checks for this implicit function by :
     * Aligning the schedules dimensions and generating the order between them.
     * If reset_static_dimesion is set to true then the computations static dimensions would be resetted to 0.
     * The execution order and static dimensions would be redefined back from \p sched_graph. i.e. the static dimensions manually specified with the low level interface would be lost.
     * Resetting static dimensions would allow to make changes to \p sched_graph using after or then, these changes would be correctly evaluated by the legality checks.
     * 
     */
    void prepare_schedules_for_legality_checks(bool reset_static_dimesion);

    /**
     * Performe a full dependency analysis RAW/WAR/WAW. The result is stored in the function's attributes
     * Before invoking this method, the user must call daisy::prepare_schedules_for_legality_checks() and must define the buffer associated with each computation.
     */
    void perform_full_dependency_analysis();

    /**
     * resets all the static beta dimensions in all the computations to Zero.
     * This would allow the execution of fuction.generate_ordering many times without issues.
     * Although, the static beta dimenesion and ordering that are not specified using the scheduling graph (after or then) would be lost.
    */
    void reset_all_static_dims_to_zero();

    /**
     * Modify the schedules of the computations of this function to reflect
     * the order specified using the high level scheduling commands.
     *
     * Commands like .after() and .before() do not directly modify the schedules
     * but rather modify the sched_graph graph.
     */
    void gen_ordering_schedules();

    /**
     * Calculate all the dependencies in the function RAW/WAW/WAR & store in the function's attributes
     * All schedules must be ordered (after or then), and with same length using:
     * 1-gen_ordering_schedules
     * 2-align_schedules
    */
    void calculate_dep_flow();

    isl_union_map* compute_dep_graph();

    int get_max_schedules_range_dim() const;

    void align_schedules();

};

}

}

