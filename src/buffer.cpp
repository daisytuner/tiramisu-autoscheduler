#include "../include/buffer.h"

#include <cassert>
#include <iostream>
#include <vector>

#include "../include/function.h"

namespace daisy {

namespace tiramisu {

buffer::buffer(
    std::string name,
    std::vector<int> dim_sizes,
    daisy::primitive_t type,
    daisy::argument_t argt
)
: argtype(argt), dim_sizes(dim_sizes), name(name), type(type)
{
    assert(!name.empty() && "Empty buffer name");
};


daisy::argument_t buffer::get_argument_type() const
{
    return argtype;
}

const std::string& buffer::get_name() const
{
    return name;
}

int buffer::get_n_dims() const
{
    return this->get_dim_sizes().size();
}

daisy::primitive_t buffer::get_elements_type() const
{
    return type;
}

const std::vector<int> &buffer::get_dim_sizes() const
{
    return dim_sizes;
}

}
}
