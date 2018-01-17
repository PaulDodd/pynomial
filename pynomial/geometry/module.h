

#pragma once
#include "sphere.h"
#include <pynomial/extern/pybind/include/pybind11/pybind11.h>
#include <pynomial/extern/pybind/include/pybind11/stl_bind.h>

namespace pynomial{
    namespace geometry{
        void export_geometry_module(pybind11::module& m)
        {
            typedef Sphere<Eigen::Vector3d> sphere_t;
            pybind11::class_<sphere_t, std::shared_ptr< sphere_t > >(m,"sphere")
                .def(pybind11::init< double >())
                .def("project", &sphere_t::project)
                .def("get_radius", &sphere_t::GetRadius)
                .def("_SetRadius", &sphere_t::SetRadius)
            ;
        }
    }
}
