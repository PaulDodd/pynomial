
#include "sphere.h"
#include "module.h"
#include <pynomial/extern/pybind/include/pybind11/pybind11.h>
#include <pynomial/extern/pybind/include/pybind11/stl_bind.h>
#include <pynomial/extern/pybind/include/pybind11/eigen.h>

using namespace pynomial::geometry;


//! Define the _hpmc python module exports
PYBIND11_PLUGIN(_geometry)
{
    pybind11::module m("_geometry");
    export_geometry_module(m);
    return m.ptr();
}
