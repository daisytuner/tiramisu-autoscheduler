#pragma once

#include <memory>

#include "ast.h"
#include "sdfg_wrapper.h"
#include "search_method.h"
#include "utils.h"

namespace daisy {

namespace tiramisu
{

class auto_scheduler
{
private:
    std::shared_ptr<sdfg_wrapper> sdfg;
        
    std::shared_ptr<search_method> searcher;
        
public:

    auto_scheduler(
      std::shared_ptr<sdfg_wrapper> sdfg,
      std::shared_ptr<search_method> searcher
    );

    std::string tune();

    static std::shared_ptr<auto_scheduler> create(
      std::shared_ptr<sdfg_wrapper> sdfg,
      std::string py_cmd_path,
      std::string py_interface_path,
      std::string method,
      size_t max_depth,
      size_t beam_size
    );

};

}

}


