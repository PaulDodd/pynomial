
#include "Network.h"
#include "NN.h"
#include "module.h"
#include <pynomial/extern/pybind/include/pybind11/pybind11.h>
#include <pynomial/extern/pybind/include/pybind11/stl_bind.h>
#include <pynomial/extern/pybind/include/pybind11/eigen.h>

using namespace pynomial::graph;

//! Define the _hpmc python module exports
PYBIND11_PLUGIN(_graph)
{
    pybind11::module m("_graph");
    export_graph_module(m);
    export_NN_module(m);
    return m.ptr();
}
