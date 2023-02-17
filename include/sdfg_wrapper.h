#pragma once

#include <string>
#include <memory>

#include "ast.h"
#include "buffer.h"
#include "function.h"
#include "computation.h"

namespace daisy {

namespace tiramisu {

class sdfg_wrapper
{
private:

    std::string name;

    std::string sdfg_path;

    std::string build_folder;

    std::shared_ptr<function> fct;

    std::shared_ptr<syntax_tree> base_ast;

    sdfg_wrapper(
        std::string sdfg_path,
        std::string build_folder,
        std::shared_ptr<function> fct
    );

public:

    std::string get_name() const;

    function& get_function();

    syntax_tree& get_base_ast();

    std::vector<float> benchmark(syntax_tree& ast);

    static std::shared_ptr<sdfg_wrapper> create(std::string sdfg_path);

};

}

}
