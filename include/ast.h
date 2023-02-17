#pragma once

#include <string>
#include <unordered_map>

#include <isl/map.h>

#include "utils.h"
#include "optimization_info.h"
#include "dnn_accesses.h"

namespace daisy {
    
namespace tiramisu
{

class syntax_tree;
class computation;

/**
 * stores the state of the computation's schedule.
*/
class state_computation
{
    friend computation;

    private:
    /**
     * Isl_map that represent the schedule of the computation
    */
    isl_map * current_schedule = NULL;

   /**
    * Computation reference that points to the real computation object.
    * Useful to stage optimizations with setting current_schedule as schedule then use computations features.
   */
    computation * staging_computation = NULL;
    /**
     * Describes wether this state is staged inside the computation or not.
    */
    bool is_state_staged;

    protected:

    public:

    /**
     * constructors
    */
    //@{
    state_computation(computation * reference);
    state_computation(state_computation const& reference);
    //@}

    /**
     * returns the isl_map
    */
   isl_map * get_inner_isl_map() const;

    /**
     * moves current schedule into the computation object as staging area.
     * mermutes the schedules between the computation and this schedule.
    */
    void move_schedule_to_staging();

    /**
     * Gets the schedule back from the computation and store it in this class instance.
    */
    void recover_schedule_from_staging();

    /**
     * get real computation in which we set the current state as a schedule inside it. 
    */
    computation * get_computation_staged();

    /**
     * get real computation without staging a schedule inside it.
    */

    computation * get_computation_unstated() const;

    /**
     * 
    */
    bool is_this_state_staged() const;

    ~state_computation()
    {
        isl_map_free(current_schedule);
    }

};

/**
 * Stores information about a computation.
 */
class computation_info
{
private:

protected:

public:
    /**
     * Pointer to the corresponding computation.
     */
    computation *comp_ptr;
    
    /**
     * List of iterators of the computation.
     */
    std::vector<dnn_iterator> iters;
    
    /**
     * List of accesses of the computation.
     */
    dnn_accesses accesses;
    
    /**
     * Number of dimensions of the output buffer.
     */
    int buffer_nb_dims;
    
    /**
     * True if this computation is a reduction.
     */
    bool is_reduction;

    /**
     * A string representing the ISL write access relation of the computation.
     */
    std::string write_access_relation;

    /**
     * The ID of the buffer where the computation is stored
     */
    int storage_buffer_id;

    /**
     * A string representing the data type of the computation
     */
    std::string data_type_str;

    /**
     * Size of the data type in Bytes
     */
    int data_type_size;

    /**
     * Some metrics about the computation.
     */
    int nb_additions;
    int nb_substractions;
    int nb_multiplications;
    int nb_divisions;
    
    /**
     * Get info about the given computation. The AST is needed to get some info.
     */
    computation_info(computation *comp, syntax_tree *ast);


    /**
     * Copy constructor
    */
    //computation_info(computation_info const& reference);

    /**
     * modifies the accesses by skewing
    */
    void set_accesses_changes_with_skewing(int first_node_depth,int alpha,int beta,int gamma,int sigma);

    /**
     * Returns the size of the computation's data type in Bytes
     */
    int get_data_type_size();
};

/**
 * A node in the AST represents a loop level.
 */
class ast_node
{
private:

protected:

public:
    /**
     * Depth of this loop level.
     */
    int depth;
    
    /**
     * Name of this loop level iterator.
     */
    std::string name;
    
    /**
     * Lower bound of this loop level iterator.
     */
    int low_bound;
    
    /**
     * Upper bound of this loop level iterator.
     */
    int up_bound;

    /**
     * True if the loop level has been parallelized
     */
    bool parallelized = false;

    /**
    * True if the loop level has been skewed
    */
    bool skewed = false;

    /**
     * List of the computations computed at this level.
     */
    std::vector<computation_info> computations;

    /**
     * Structure that holds the state of each computation inside this node.
    */
    std::vector<state_computation> isl_states;

	/**
	 * Next loop levels.
	 */
    std::vector<ast_node*> children;
    
    /**
     * Parent of this loop level.
     */
    ast_node *parent = nullptr;
    
	/**
	 * Create an empty AST node.
	 */
	ast_node() {}

	/**
	 * Create an AST node from the given computation.
	 */
	ast_node(computation *comp, syntax_tree* ast);
        
    ~ast_node()
    {
        for (ast_node* child : children)
            delete child;
    }
    
    /**
     * Return the extent of this loop level.
     */
    int get_extent() const { return up_bound - low_bound + 1; }
    
    /**
     * Copy this node and return the copy.
     */
    ast_node* copy_node() const;

    /**
     * Copy the tree rooted at this node into new_node and return
     * a pointer to the copied version of node_to_find.
     *
     * This function is used if you want to copy a tree, and need the new location
     * of a node.
     */
    ast_node* copy_and_return_node(ast_node *new_node, ast_node *node_to_find) const;
    
    /**
     * Fill the given array with the extents of the innermost loop levels
     * contained in this subtree.
     */
    void get_innermost_extents(std::vector<int>& extents) const;
    
    /**
     * Get the computations located at the innermost loop levels of this subtree.
     */
    void get_innermost_computations(std::vector<computation*>& comps);
    
    /**
     * Fill the given array with the nodes representing the innermost loop levels
     * contained in this subtree.
     */
    void get_innermost_nodes(std::vector<ast_node*>& nodes);
    
    /**
     * Get the root of the tree to which this node belongs to.
     */
    ast_node* get_root_node();
    
    /**
     * Get the node located at the leftmost side of the tree.
     */
    ast_node* get_leftmost_node();
    
    /**
     * Get the node located at the rightmost side of the tree.
     */
    ast_node* get_rightmost_node();

    /**
     * get all the nodes starting from root that have 1 child, 
     * i.e. the shared nodes between all computations
    */
    void get_shared_nodes_from_outermost(std::vector<ast_node*>& shared) const;

    /**
     * Recompute the depth of each node of the tree rooted at
     * this node, with the given depth being the depth of this node.
     */
    void update_depth(int depth);

    /**
     * Fill the given array with all the computations computed 
     * at this level and the levels below.
     */
    void get_all_computations(std::vector<computation*>& comps);

    /**
     * Starting from this node, get the number of nodes that have no computation,
     * and only one child.
     */
    int get_loop_levels_chain_depth() const;

    /**
     * Print the subtree rooted at this node.
     */
    void print_node() const;

    /**
     * Print the subtree of isl_states
    */
    void print_isl_states() const;

    /**
     * prints the computations's accesses of this AST
    */
    void print_computations_accesses() const;

    /**
     * Changes the access within computation_info after applying the skewing optimization
    */
    void transform_accesses_with_skewing(int a,int b);

    void transform_accesses_with_skewing_positive(int a,int b,int c, int d);

    void set_accesses_changes_with_skewing(int first_node_depth,int alpha,int beta,int gamma,int sigma);

    /**
     * create initial isl_states from current computations
    */
    void create_initial_states();

    /**
     * stage the isl_states to the real computations
    */
    void stage_isl_states();

    /**
     * recover the states from the computations
    */
    void recover_isl_states();

    /**
     * pushs all the computations inside this node recursively 
    */
    void collect_all_computation(std::vector<computation_info*>& vector);

    /**
     * get the extent of this node, i.e:
     * return upper_bound - lower_bound
    */
    int get_node_loop_extent() const;
};

class syntax_tree
{
private:
    /**
     * Transform the AST using the order specified by the user
     * with the command "after".
     */
    void order_computations();

protected:

public:
    /**
     * The function represented by the AST.
     */
    std::shared_ptr<function> fct;
    
    /**
      * AST root nodes.
      */
    std::vector<ast_node*> roots;
    
    /**
     * The list of computations contained in this AST.
     */
    std::vector<computation*> computations_list;
    
    /**
     * A mapping between each computation and the node where it is contained.
     */
    std::unordered_map<computation*, ast_node*> computations_mapping;
    
    /**
     * The list of buffers used by the program.
     */
    std::vector<std::string> buffers_list;
    
    /**
     * A mapping between each computation and the name of the buffer where
     * it is stored.
     */
    std::unordered_map<std::string, std::string> buffers_mapping;

    /**
     * An evaluation given by a class of type evaluation_function.
     */
    float evaluation;
    
    /**
     * The depth of this AST in a search method.
     * Used to keep track of the depth reached by a search method.
     */
    int search_depth = 0;
    
    /**
     *
     */
    std::vector<optimization_info> previous_optims;
    
    /**
     *
     */
    std::vector<optimization_info> new_optims;
    
    /**
     * The iterators in JSON format.
     * Use by the class evaluate_by_learning_model.
     */
    std::string iterators_json;
    
    /**
     * The structure represented by this AST in JSON format.
     * Use by the class evaluate_by_learning_model.
     */
    std::string tree_structure_json;
        
    /**
     * Create an empty AST.
     */
    syntax_tree() {}
    
    /**
     * Create an AST from the given function.
     */
    syntax_tree(std::shared_ptr<function> fct);
    
    ~syntax_tree()
    {
        for (ast_node *node : roots)
            delete node;
    }
    
    std::vector<computation*> const& get_computations() const { return computations_list; }
    
    /**
     * Transform the AST by applying the last optimization found.
     */
    void transform_ast();

    /**
     * Transform the AST by applying the given optimization.
     */
    void transform_ast(optimization_info const& opt);

    /**
     * These methods are used to transform the AST given a specific type of optimization.
     */
    void transform_ast_by_tiling(optimization_info const& opt);
    void transform_ast_by_interchange(optimization_info const& opt);
    void transform_ast_by_parallelism(const optimization_info &info);
    
    /**
     * Copy this AST, and return the copy.
     */
    std::shared_ptr<syntax_tree> copy_ast() const;

    /**
     * Copy this AST to new_ast and return a pointer to the copied version of node_to_find.
     * Can be used if you want to copy this AST, and need to find the new location of a node.
     */
    ast_node* copy_and_return_node(syntax_tree& new_ast, ast_node *node_to_find) const;
    
    /**
     * Clear computations_mapping, and recompute the nodes where computations are stored.
     */
    void recompute_computations_mapping();
    
    /**
     * Recursive subroutine used by void recompute_computations_mapping();
     */
    void recompute_computations_mapping(ast_node *node);
    
    /**
     * Return the schedule of this AST.
     */
    std::vector<optimization_info> get_schedule() const;
    
    /**
     * Add the content of new_optims to previous_optims and
     * clear new_optims.
     */
	void clear_new_optimizations();
    
    /**
     * Return the node corresponding to the given loop level of the given computation.
     */
    ast_node* find_node_by_level(computation *comp, int level);
    
    /**
     * Get the extents of the loop levels shared by all computations.
     */
    std::vector<int> get_shared_levels_extents() const;
    
    /**
     * Get the extents of all the innermost loop levels.
     */
    std::vector<int> get_innermost_extents() const;
    
    /*
     * Get the computations located at the innermost loop levels.
     */
    std::vector<computation*> get_innermost_computations();
    
    /**
     * Return the nodes representing the innermost loop levels.
     */
    std::vector<ast_node*> get_innermost_nodes() const;

    /**
     * get all the nodes starting from root that have 1 child, 
     * i.e. the shared nodes between all computations
    */
    void get_shared_nodes_from_outermost(std::vector<ast_node*>& shared) const;
    
    /**
     * Return the position, in the list of the buffers, of the buffer where
     * the given computation is stored.
     */
    int get_buffer_id_from_computation_name(std::string comp_name);
    
    /**
     * Return the position of the given buffer in the list of buffers.
     */
    int get_buffer_id(std::string const& buf_name) const;

    /**
     * Print the AST to stdout.
     */
    void print_ast() const;

    /**
     * prints the computations's accesses of this AST
    */
    void print_computations_accesses() const;

    /**
     * Print information about the applied optimizations
     */
    void print_new_optims() const;

    void print_previous_optims() const;

    void print_isl_states() const;

    void create_initial_isl_state() const;

    /**
     * push isl_states to the real computations to use tiramisu API
    */
    void stage_isl_states() const;

    /**
     * Recover the isl_states from the computation 
    */
    void recover_isl_states() const;

    /**
     * Check the correctness of the list of applied structurel optimizations
    */
    bool ast_is_legal() const;

    /**
     * Encodes the transformations applied to the ast as a string, this is used for saving the exploration trace while
     * sampling schedules
     */
    std::string get_schedule_str();

    /**
     * Predicts if the schedule applied to the ast is worth evaluating and exploring further.
     */
    bool schedule_is_prunable();

    /**
     * Checks if the AST's evaluation can be predicted using manual engineered rules
     */
    bool can_set_default_evaluation();
};

/**
 * The candidate_trace class is used to recursively record the visited candidates during the search space exploration
 */
class candidate_trace
{
private:
    /**
     * A numerical id is assigned for each explored candidate
     */
    int candidate_id;

public:
    int get_candidate_id() const;

private:
    /**
     * The search depth where this candidate was visited
     */
    int exploration_depth;

    /**
     * candidate's evaluation
     */
    float evaluation;

    /**
     * The applied schedule, encoded as a string
     */
    std::string schedule_str;

    /**
     * The child candidates that are derived from this candidate
     */
    std::vector<candidate_trace*> child_candidates;

public:

    candidate_trace(syntax_tree* ast, int candidate_id);

    ~candidate_trace();

    /**
     * Initializes a new child candidate from an AST and adds it to the children list
     */
    void add_child_path(syntax_tree* ast, int candidate_id);

    /**
     * Recursively exports the exploration trace into a JSON format
     */
    std::string get_exploration_trace_json();

    /**
     * A mapping between ASTs and explored child candidates
     */
    std::unordered_map<syntax_tree*, candidate_trace*> child_mappings;

};


} }


