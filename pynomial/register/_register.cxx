
#include "module.h"
#include <pynomial/extern/pybind/include/pybind11/pybind11.h>
#include <pynomial/extern/pybind/include/pybind11/stl_bind.h>
#include <pynomial/extern/pybind/include/pybind11/eigen.h>

using namespace pynomial::_register;
PYBIND11_PLUGIN(_register)
{
    pybind11::module m("_register");
    export_register_module(m);
    export_transform_module(m);
    return m.ptr();
}
