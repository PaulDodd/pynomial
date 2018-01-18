

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
            typedef SphereArc<Eigen::Vector3d, copy_type::deep> sphere_arc_t;
            pybind11::class_<sphere_arc_t, std::shared_ptr< sphere_arc_t > >(m,"sphere_arc")
                .def(pybind11::init< const sphere_t&,
                                    const Eigen::Vector3d&,
                                    const Eigen::Vector3d& >())
                .def("_PlaneNormal", &sphere_arc_t::PlaneNormal)
                .def("_Length", &sphere_arc_t::Length)
                .def("_MinPhi", &sphere_arc_t::MinPhi)
                .def("_MaxPhi", &sphere_arc_t::MaxPhi)
                .def("_MinTheta", &sphere_arc_t::MinTheta)
                .def("_MaxTheta", &sphere_arc_t::MaxTheta)
            ;
        }
    }
}
