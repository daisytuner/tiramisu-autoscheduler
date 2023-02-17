#include "../include/function.h"

#include <cassert>
#include <queue>

#include "../include/core.h"
#include "../include/computation.h"

namespace daisy {

namespace tiramisu {

/**
 * Compute the accesses of the RHS of the computation
 * \p comp and store them in the accesses vector.
 * If \p return_buffer_accesses is set to true, this function returns access functions to buffers. Otherwise it returns access functions to computations.
 */
void get_rhs_accesses(const function *func, const computation *comp, std::vector<isl_map *> &accesses, bool return_buffer_accesses)
{
    // DEBUG_FCT_NAME(3);
    // DEBUG_INDENT(4);

    // const daisy::expr &rhs = comp->get_expr();
    // if (comp->is_wait()) {
    //     // Need to swap the access map of the operation we wait on
    //     daisy::computation *waitee = func->get_computation_by_name(rhs.get_name())[0];
    //     isl_map *orig = isl_map_copy(waitee->get_access_relation());
    //     waitee->set_access(waitee->wait_access_map);
    //     generator::traverse_expr_and_extract_accesses(func, comp, rhs, accesses, return_buffer_accesses);
    //     waitee->set_access(orig);
    // } else {
    //     generator::traverse_expr_and_extract_accesses(func, comp, rhs, accesses, return_buffer_accesses);
    // }

    // DEBUG_INDENT(-4);
    // DEBUG_FCT_NAME(3);
}

/**
 * Construct a function with the name \p name.
 */
function::function(std::string name)
{
    assert(!name.empty() && ("Empty function name"));

    this->name = name;
    // this->halide_stmt = Halide::Internal::Stmt();
    // this->ast = NULL;
    this->context_set = NULL;
    this->use_low_level_scheduling_commands = false;
    // this->_needs_rank_call = false;

    // Allocate an ISL context.  This ISL context will be used by
    // the ISL library calls within Tiramisu.
    this->ctx = isl_ctx_alloc();
};

const std::string& function::get_name() const
{
    return name;
}

isl_set* function::get_program_context() const
{
    if (context_set != NULL)
    {
        return isl_set_copy(context_set);
    }
    else
    {
        return NULL;
    }
}

isl_ctx* function::get_isl_ctx() const
{
    return ctx;
}

const std::map<std::string, std::shared_ptr<buffer>> function::get_buffers() const
{
    return buffers_list;
}

void function::add_buffer(std::pair<std::string, std::shared_ptr<buffer>> buf)
{
    assert(!buf.first.empty() && ("Empty buffer name."));
    assert((buf.second != NULL) && ("Empty buffer."));

    this->buffers_list.insert(buf);
}

const std::vector<computation*> function::get_computations() const
{
    return body;
}

void function::add_computation(computation *cpt)
{
    // DEBUG_FCT_NAME(3);
    // DEBUG_INDENT(4);

    assert(cpt != NULL);

    this->body.push_back(cpt);
    if (cpt->should_schedule_this_computation())
        this->starting_computations.insert(cpt);

    // DEBUG_INDENT(-4);
}

std::vector<computation*> function::get_computation_by_name(std::string name) const
{
    assert(!name.empty());

    std::vector<computation *> res_comp;
    for (const auto &comp : this->get_computations())
    {
        if (name == comp->get_name())
        {
            res_comp.push_back(comp);
        }
    }

    return res_comp;
}

isl_union_map* function::get_schedule() const
{
    isl_union_map *result = NULL;
    isl_space *space = NULL;

    if (!this->body.empty())
    {
        space = isl_map_get_space(this->body[0]->get_schedule());
    }
    else
    {
        return NULL;
    }

    assert(space != NULL);
    result = isl_union_map_empty(isl_space_copy(space));

    for (const auto &cpt : this->body)
    {
        isl_map *m = isl_map_copy(cpt->get_schedule());
        result = isl_union_map_union(isl_union_map_from_map(m), result);
    }

    result = isl_union_map_intersect_domain(result, this->get_iteration_domain());

    return result;
}


bool function::is_sched_graph_tree_dfs(computation * comp,std::unordered_set<computation *> &visited)
{
    // Do not visit anything that was already returned
    if (visited.find(comp) != visited.end())
        return false;

    visited.insert(comp);

    for (auto &edge: this->sched_graph[comp])
    {
        if (!is_sched_graph_tree_dfs(edge.first, visited))
            return false;
    }

    return true;
}

bool function::is_sched_graph_tree()
{
    // DEBUG_FCT_NAME(3);
    // DEBUG_INDENT(4);

    if (this->starting_computations.size() != 1)
    {
        // DEBUG_INDENT(-4);
        return false;
    }

    // Contains all nodes that have been visited
    std::unordered_set<computation *> visited;

    for (auto &comp: this->starting_computations)
    {
        if (!is_sched_graph_tree_dfs(comp, visited))
        {
            // DEBUG_INDENT(-4);
            return false;
        }
    }

    // DEBUG_INDENT(-4);
    return true;
}


void function::add_parallel_dimension(std::string stmt_name, int vec_dim)
{
    this->parallel_dimensions.push_back({stmt_name, vec_dim});
}

void function::reset_schedules()
{
    for (computation* comp : this->body)
    {
        comp->set_identity_schedule_based_on_iteration_domain();
    }
    
    remove_dimension_tags();
    clear_sched_graph();
}

void function::remove_dimension_tags()
{
    parallel_dimensions.clear();
}

void function::reset_all_static_dims_to_zero()
{   
    // DEBUG_FCT_NAME(3);
    // DEBUG_INDENT(4);

    for(auto computation:this->get_computations())
    {
        isl_map * schedule = isl_map_copy(computation->get_schedule());
        isl_space * range_space = isl_set_get_space(isl_map_range(isl_map_copy(schedule)));

        isl_map * transformation_map = isl_map_universe(isl_space_map_from_set(range_space));

        int m1  = isl_map_dim(schedule, isl_dim_out);

        // pos 0 is static always 0 by default, the we will have static then dynamic for all the rest.
        transformation_map = isl_map_fix_si(transformation_map,isl_dim_out,0,0);

        for(int i=1; i<m1; i++)
        {
            if(i%2 == 1)
            {// case of static dimension, fix position to 0.
                transformation_map = isl_map_fix_si(transformation_map,isl_dim_out,i,0);
            }
            else
            {// equate input and output in case of dynamic dimensions
                transformation_map = isl_map_equate(transformation_map,isl_dim_out,i,isl_dim_in,i);
            }
        }
        // DEBUG(3, daisy::str_dump(" Initial schedule before initialization of beta dimensions : "+std::string(isl_map_to_str(schedule))));
        // DEBUG(3, daisy::str_dump(" Transformation Map : "+std::string(isl_map_to_str(transformation_map))));

        schedule = isl_map_apply_range(schedule,transformation_map);

        // DEBUG(3, daisy::str_dump(" Initialized Schedule : "+std::string(isl_map_to_str(schedule))));

        computation->set_schedule(schedule);

    }

    // DEBUG_INDENT(-4);
}


void function::clear_sched_graph()
{
    sched_graph.clear();
    sched_graph_reversed.clear();
}

int function::get_max_schedules_range_dim() const
{
    int max_dim = 0;
    for (const auto &comp : this->get_computations())
    {
        if(comp->schedule_this_computation){

            isl_map *sched = comp->get_schedule();
            int m = isl_map_dim(sched, isl_dim_out);
            max_dim = std::max(max_dim, m);
        }
    }

    return max_dim;
}


void function::align_schedules()
{
    int max_dim = this->get_max_schedules_range_dim();

    for (auto &comp : this->get_computations())
    {
        
        isl_map *dup_sched = comp->get_schedule();
        assert(dup_sched != NULL);
        dup_sched = isl_map_align_range_dims(dup_sched, max_dim);
        comp->set_schedule(dup_sched);
        comp->name_unnamed_time_space_dimensions();
        
    }
}

bool function::loop_parallelization_is_legal(std::string i, std::vector<computation*> fused_computations)
{
    assert(!this->get_name().empty());
    assert(this->dep_read_after_write != NULL );
    assert(this->dep_write_after_write != NULL );
    assert(this->dep_write_after_read != NULL );
    assert(fused_computations.size()>0);


    computation * first_computation = fused_computations[0];
    
    std::vector<std::string> original_loop_level_names = first_computation->get_loop_level_names();

    std::vector<int> dimensions =
        first_computation->get_loop_level_numbers_from_dimension_names({i});

    first_computation->check_dimensions_validity(dimensions);

    int dim_parallel = dimensions[0];
    int par_dim = 1 + (dim_parallel * 2 + 1);

    // Extracting deps

     isl_union_map * read_after_write_dep = isl_union_map_range_factor_domain(
        isl_union_map_copy(this->dep_read_after_write));

    isl_union_map * write_after_read_dep = isl_union_map_range_factor_domain(
        isl_union_map_copy(this->dep_write_after_read));

    isl_union_map * write_after_write_dep = isl_union_map_range_factor_domain(
        isl_union_map_copy(this->dep_write_after_write));

    isl_union_map * all_deps = isl_union_map_union(
        read_after_write_dep,
        write_after_read_dep
        );

    // all the deps in 1 union map
    all_deps = isl_union_map_union(all_deps,write_after_write_dep);

    // all current schedules in 1 union map
    std::string empty_union = "{}";
    std::string empty_time  = "";

    isl_union_map * schedules = isl_union_map_read_from_str(this->get_isl_ctx(),empty_union.c_str());

    isl_map * schedule_itr = NULL;

    for( auto& computation: fused_computations)
    {
        schedule_itr = isl_map_copy(computation->get_schedule());

        schedule_itr = isl_map_set_tuple_name(schedule_itr,isl_dim_out,empty_time.c_str());

        schedules = isl_union_map_union(schedules,isl_union_map_from_map(schedule_itr));

    }

    // application to discard unused dep & represent them in their time space

    all_deps = isl_union_map_apply_range(all_deps,isl_union_map_copy(schedules));

    all_deps = isl_union_map_apply_domain(all_deps,isl_union_map_copy(schedules));

    if(isl_union_map_is_empty(all_deps))
    {
        return true;
    }

    isl_map * equation_map = isl_map_from_union_map(all_deps);

    bool overall_legality = false;

    /*
        isl_equate adds restriction that both domain and range positions are equal
        we suppose that legality of the lexicographical order is checked elsewhere, so we only need to check for loop caried dependencies.
        if adding equation of == between input set & output set of map for a dimension strictly before the parallel one is empty means: 
            dep is not a carried one for the parallel loop lvl, so parallelism is legal.

        else 
            if all previous equations added does not make the map empty then the last possibility is:
                dep is within the same loop iteration, parallel is true ( true if equate doesn't make the map empty)
                else it's false
                
    */
    for(int i=0;i<par_dim;i++)
    {
        equation_map = isl_map_equate(equation_map,isl_dim_in,i,isl_dim_out,i);


        if(isl_map_is_empty(equation_map))
        {
            overall_legality = true;
            break;
        }
    
    }


    if(!overall_legality)
    {
        isl_map * equation_map_final = isl_map_equate(isl_map_copy(equation_map),isl_dim_in,par_dim,isl_dim_out,par_dim);

        if(isl_map_is_equal(equation_map,equation_map_final) == isl_bool_false)
        {
            overall_legality = false;
        }
        else{
            overall_legality = true;
        }
        isl_map_free(equation_map_final);
    }
    
    isl_map_free(equation_map);
    isl_union_map_free(schedules);

    return overall_legality;
}

void function::prepare_schedules_for_legality_checks(bool reset_static_dimesion)
{
    this->align_schedules();

    if(reset_static_dimesion == true)
    {
        this->reset_all_static_dims_to_zero();
    }

    this->gen_ordering_schedules();
}

bool function::check_legality_for_function()
{
    assert(this->dep_read_after_write!=NULL);

    isl_union_map * all_deps = isl_union_map_range_factor_domain(
        isl_union_map_copy(this->dep_read_after_write));

    all_deps = isl_union_map_union(all_deps,
        isl_union_map_range_factor_domain(isl_union_map_copy(this->dep_write_after_read)));

    all_deps = isl_union_map_union(all_deps, 
        isl_union_map_range_factor_domain(isl_union_map_copy(this->dep_write_after_write)));

    isl_union_map * universe_of_all_deps = isl_union_map_universe(all_deps);

    std::vector<isl_map *> all_basic_maps;
    
    auto f = [](isl_map * bmap,void * user) { 

        std::vector<isl_map *>& myName = *reinterpret_cast<std::vector<isl_map*>*>(user);
     
        myName.push_back(bmap);
        return isl_stat_ok;
    };
    
    isl_stat (*fun_ptr)(isl_map * p,void * m) = (f);

    isl_union_map_foreach_map(universe_of_all_deps,fun_ptr,(void * ) &all_basic_maps);

    isl_set * left_hs = NULL;
    isl_set * right_hs = NULL; // hand side

    computation * left_comp = NULL;
    computation * right_comp = NULL;

    std::string left_computation_name =  "";
    std::string right_computation_name = "";

    bool over_all_legality = true;
    
    for(auto& space_dep:all_basic_maps)
    {
        //DEBUG(3, daisy::str_dump(" the map of deps is  "+std::string(isl_map_to_str(space_dep))));

        left_hs = isl_map_domain(isl_map_copy(space_dep));
        right_hs = isl_map_range(isl_map_copy(space_dep));

        left_computation_name =  isl_space_get_tuple_name(
            isl_set_get_space(left_hs),isl_dim_set);

        right_computation_name =  isl_space_get_tuple_name(
            isl_set_get_space(right_hs),isl_dim_set);

        //DEBUG(3, daisy::str_dump(" checking legality of dependences "+left_computation_name+" -> "+right_computation_name));
        
        left_comp = this->get_computation_by_name(left_computation_name)[0];
        right_comp = this->get_computation_by_name(right_computation_name)[0];

        if( left_comp->involved_subset_of_dependencies_is_legal(right_comp) == false )
        {
            over_all_legality = false;
            break;
        }
    }

    isl_union_map_free(universe_of_all_deps);

    return over_all_legality;
}

void function::perform_full_dependency_analysis()
{

    // align schedules and order schedules
    this->align_schedules();
    this->gen_ordering_schedules();
    // could save default schedules and order here
    this->calculate_dep_flow();
    
}


void function::gen_ordering_schedules()
{
    // DEBUG_FCT_NAME(3);
    // DEBUG_INDENT(4);

    if (this->use_low_level_scheduling_commands)
    {
        // DEBUG(3, daisy::str_dump("Low level scheduling commands were used."));
        // DEBUG(3, daisy::str_dump("Discarding high level scheduling commands."));
        return;
    }

    // this->dump_sched_graph();

    if(this->is_sched_graph_tree())
    {
        // DEBUG(3, daisy::str_dump("this->is_sched_graph_tree(): true."));

        std::priority_queue<int> level_to_check;
        std::unordered_map<int, std::deque<computation *>> level_queue;

        auto current_comp = *(this->starting_computations.begin());

        auto init_sched = automatically_allocated;
        init_sched.push_back(current_comp);

        for (auto it = init_sched.begin(); it != init_sched.end() && it + 1 != init_sched.end(); it++)
            (*(it+1))->after_low_level(**it, computation::root_dimension);

        bool comps_remain = true;
        while(comps_remain)
        {
            for (auto &edge: this->sched_graph[current_comp])
            {
                if (level_queue[edge.second].size() == 0)
                    level_to_check.push(edge.second);

                level_queue[edge.second].push_back(edge.first);
            }

            comps_remain = level_to_check.size() > 0;
            // If we haven't exhausted all computations
            if (comps_remain)
            {
                int fuse_level = level_to_check.top();
                auto next_comp = level_queue[fuse_level].front();
                level_queue[fuse_level].pop_front();

                // assert(this->get_max_iteration_domains_dim() > fuse_level);

                next_comp->after_low_level((*current_comp), fuse_level);

                current_comp = next_comp;
                if (level_queue[fuse_level].size() == 0)
                    level_to_check.pop();
            }
        }
    }
    // else
    // {
    //     DEBUG(3, daisy::str_dump("this->is_sched_graph_tree(): false."));
    // }

    // DEBUG_INDENT(-4);
}

isl_union_set* function::get_iteration_domain() const
{
    isl_union_set *result = NULL;
    isl_space *space = NULL;

    if (!this->body.empty())
    {
        space = isl_set_get_space(this->body[0]->get_iteration_domain());
    }
    else
    {
        return NULL;
    }

    assert(space != NULL);
    result = isl_union_set_empty(space);

    for (const auto &cpt : this->body)
    {
        if (cpt->should_schedule_this_computation())
        {
            isl_set *cpt_iter_space = isl_set_copy(cpt->get_iteration_domain());
            result = isl_union_set_union(isl_union_set_from_set(cpt_iter_space), result);
        }
    }

    return result;
}

void function::calculate_dep_flow()
{
    // DEBUG_FCT_NAME(3);
    // DEBUG_INDENT(4);

    assert(this->get_computations().size() > 0);
    assert(this->get_computations()[0]->get_schedule() != NULL);

    // DEBUG(3, daisy::str_dump(" generating depandencies graph"));

    isl_union_map * ref_res = this->compute_dep_graph();

    if(ref_res == NULL)
    {
        // no deps fill with empty union maps

        std::string str_map = "{}";

        this->dep_read_after_write = isl_union_map_read_from_str(this->get_isl_ctx(),str_map.c_str());

        this->dep_write_after_read = isl_union_map_read_from_str(this->get_isl_ctx(),str_map.c_str());;

        this->dep_write_after_write = isl_union_map_read_from_str(this->get_isl_ctx(),str_map.c_str());;

        this->live_in_access = isl_union_map_read_from_str(this->get_isl_ctx(),str_map.c_str());;

        this->live_out_access = isl_union_map_read_from_str(this->get_isl_ctx(),str_map.c_str());;

        // DEBUG(3, daisy::str_dump(" No deps detected just empty maps "));

        // DEBUG_INDENT(-4);

        return ;
    }

    isl_union_map * ref_graph = isl_union_map_reverse(ref_res);

    

    // DEBUG(3, daisy::str_dump(" the referencing union map is for dependecy analysis: "+std::string(isl_union_map_to_str(ref_graph))));

    int time_space_dim = isl_map_dim(this->get_computations()[0]->get_schedule(), isl_dim_out);

    std::string to_time_space_map_str = "[";

    std::string to_time_space_map_str_2 = "[";

    for(int i=0; i < time_space_dim; i++)
    {
        to_time_space_map_str+="t"+std::to_string(i);
        to_time_space_map_str_2+="t"+std::to_string(i);

        if(i != (time_space_dim - 1))
        {
            to_time_space_map_str+=",";
            to_time_space_map_str_2+=",";
            
        }

    }
    std::string ready_time_str = to_time_space_map_str+"]->" + to_time_space_map_str_2+"]";// without {} yet

    // DEBUG(3, daisy::str_dump(" using to generate time stamp tmp map "+ready_time_str));


    std::string access_start = "{}";

    // S0[i,j] -> buff[i] the writing stmt
    isl_union_map * write_access = isl_union_map_read_from_str(this->get_isl_ctx(),access_start.c_str());

    isl_union_map * isl_schedule = isl_union_map_read_from_str(this->get_isl_ctx(),access_start.c_str());


    std::string identity = "";

    isl_map * isl_identity = NULL;
    
    for(auto& comput : this->get_computations())
    {
        identity = "{"+comput->get_name() +ready_time_str + "}";

        isl_identity = isl_map_read_from_str(this->get_isl_ctx(),identity.c_str());

        // TODO : use default schedule instead when save/restore states is implemented 
        isl_map * corrected = isl_map_apply_range(isl_map_copy(comput->get_schedule()),isl_identity);

        // DEBUG(10, daisy::str_dump(" - > compuatation's schedule to time stamp op result is : "+std::string(isl_map_to_str(corrected))));
        
        isl_schedule = isl_union_map_union(isl_schedule , isl_union_map_from_map(corrected));

        write_access = isl_union_map_union(write_access,isl_union_map_from_map(isl_map_copy(comput->get_access_relation())));
        
    } 

    isl_union_set * iteration_domains = this->get_iteration_domain();

    isl_union_map * write_acccess_without_domain = isl_union_map_copy(write_access);

    write_access = isl_union_map_intersect_domain(write_access, isl_union_set_copy(iteration_domains));

    isl_schedule = isl_union_map_intersect_domain(isl_schedule, isl_union_set_copy(iteration_domains));
    
    isl_union_map * read_access = isl_union_map_apply_range(
        isl_union_map_copy(ref_graph),
        write_acccess_without_domain
    );

    read_access = isl_union_map_intersect_domain(read_access, isl_union_set_copy(iteration_domains));

    //combine reads previous with their access to establish the read access S0[i,j] -> buf2[j] in read 

    // DEBUG(3, daisy::str_dump("the overall function schedule is : "+std::string(isl_union_map_to_str(isl_schedule))));

    // DEBUG(3, daisy::str_dump("the write access for computations is : "+std::string(isl_union_map_to_str(write_access))));

    // DEBUG(3, daisy::str_dump(" The read access for computations : "+std::string(isl_union_map_to_str(read_access))));

    isl_union_access_info *info = isl_union_access_info_from_sink( isl_union_map_copy(read_access));

    info = isl_union_access_info_set_schedule_map(info,isl_union_map_copy(isl_schedule));

    info = isl_union_access_info_set_must_source(info,isl_union_map_copy(write_access));

    isl_union_flow * flow = isl_union_access_info_compute_flow(info);

    //DEBUG(3, daisy::str_dump(" dependency analysis with must for read after write ( no predicats ) result  : "+std::string(isl_union_flow_to_str(flow))));

    isl_union_map * read_after_write_dep = isl_union_flow_get_full_must_dependence(flow);

    isl_union_map * read_from_outside = isl_union_flow_get_must_no_source(flow);

    // DEBUG(3, daisy::str_dump(" read after write True dependencies are in the form { last_write_access -> the read statement } : "+std::string(isl_union_map_to_str(read_after_write_dep))));
       
    // DEBUG(3, daisy::str_dump(" live-in : the computations / statement with these read access have not been written in this function (outside value)  : "+std::string(isl_union_map_to_str(read_from_outside))));
    

    info = isl_union_access_info_from_sink(isl_union_map_copy(write_access));

    info = isl_union_access_info_set_schedule_map(info,isl_union_map_copy(isl_schedule));

    info = isl_union_access_info_set_must_source(info,isl_union_map_copy(write_access));

    flow = isl_union_access_info_compute_flow(info);

    isl_union_map * write_after_write_dep = isl_union_flow_get_full_must_dependence(flow);

    // DEBUG(3, daisy::str_dump(" write after write dependencies are { last_previous_write -> new write stmt } : "+std::string(isl_union_map_to_str(write_after_write_dep))));


    isl_union_map * not_last_writes = isl_union_map_range_factor_range( isl_union_map_copy(write_after_write_dep));
    
    isl_union_map * live_out = isl_union_map_subtract(
        isl_union_map_copy(write_access),
        isl_union_map_copy(not_last_writes)
    );

    live_out = isl_union_map_intersect_domain(live_out, this->get_iteration_domain());

    // DEBUG(3, daisy::str_dump(" live out last access are : "+std::string(isl_union_map_to_str(live_out))));

    isl_union_map * read_without_write_stmt = isl_union_map_subtract(isl_union_map_copy(read_access), isl_union_map_copy(write_access));

    info = isl_union_access_info_from_sink(isl_union_map_copy(write_access));

    info = isl_union_access_info_set_schedule_map(info,isl_union_map_copy(isl_schedule));

    info = isl_union_access_info_set_may_source(info,isl_union_map_copy(read_without_write_stmt));

    info = isl_union_access_info_set_kill(info,isl_union_map_copy(write_access));

    flow = isl_union_access_info_compute_flow(info);

    //DEBUG(3, daisy::str_dump(" dependency analysis for WAR dep : "+std::string(isl_union_flow_to_str(flow))));

    isl_union_map * anti_dependencies = isl_union_flow_get_full_may_dependence(flow);

    // DEBUG(3, daisy::str_dump(" write after read anti_dependencies are in the form { last_previous_read -> new write stmt } : "+std::string(isl_union_map_to_str(anti_dependencies))));

    //DEBUG(3, daisy::str_dump(" the initialisation stmt writes with no previous read before are : "+std::string(isl_union_map_to_str(initialisation_access))));
      
    this->dep_read_after_write = read_after_write_dep;

    this->dep_write_after_read = anti_dependencies;

    this->dep_write_after_write = write_after_write_dep;

    this->live_in_access = read_from_outside;

    this->live_out_access = live_out;
    
    // DEBUG_INDENT(-4);

}

isl_union_map *function::compute_dep_graph() {
    isl_union_map *result = NULL;

    for (const auto &consumer : this->get_computations()) {

        isl_union_map *accesses_union_map = NULL;
        std::vector < isl_map * > accesses_vector;
        get_rhs_accesses(this, consumer, accesses_vector, false);

        if (!accesses_vector.empty()) {
            // Create a union map of the accesses to the producer.
            if (accesses_union_map == NULL) {
                isl_space *space = isl_map_get_space(accesses_vector[0]);
                assert(space != NULL);
                accesses_union_map = isl_union_map_empty(space);
            }

            for (size_t i = 0; i < accesses_vector.size(); ++i) {
                isl_map *reverse_access = isl_map_reverse(accesses_vector[i]);
                accesses_union_map = isl_union_map_union(isl_union_map_from_map(reverse_access),
                                                         accesses_union_map);
            }

            if (result == NULL) {
                result = isl_union_map_copy(accesses_union_map);
                isl_union_map_free(accesses_union_map);
            } else {
                result = isl_union_map_union(result, accesses_union_map);
            }
        }
    }

    return result;
}

}
}
