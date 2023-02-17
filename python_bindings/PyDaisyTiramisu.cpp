#include <memory>
#include <pybind11/pybind11.h>

#include "../include/sdfg_wrapper.h"
#include "../include/auto_scheduler.h"

namespace py = pybind11;

using namespace daisy;
using namespace daisy::tiramisu;

PYBIND11_MODULE(daisy_tiramisu, m) {
    py::class_<sdfg_wrapper, std::shared_ptr<sdfg_wrapper>>(m, "SDFGWrapper")
        .def(py::init(&sdfg_wrapper::create));

    py::class_<auto_scheduler, std::shared_ptr<auto_scheduler>>(m, "AutoScheduler")
        .def(py::init(&auto_scheduler::create))
        .def("tune", &auto_scheduler::tune);
    
}
