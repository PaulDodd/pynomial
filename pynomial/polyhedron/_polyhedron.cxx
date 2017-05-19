
#include "Polyhedra.h"
#include "Net.h"
#include "module.h"
#include <pynomial/extern/pybind/include/pybind11/pybind11.h>
#include <pynomial/extern/pybind/include/pybind11/stl_bind.h>
#include <pynomial/extern/pybind/include/pybind11/eigen.h>

using namespace pynomial::polyhedron;

//! Define the _hpmc python module exports
PYBIND11_PLUGIN(_polyhedron)
{
    pybind11::module m("_polyhedron");
    export_polyhedron_module(m);
    m.def("intersection", &intersection);
    // export_net_module(m);
    return m.ptr();
}
