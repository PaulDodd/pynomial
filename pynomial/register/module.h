#pragma once
#include "brute_force.h"
#include "transform.h"
#include <pynomial/extern/pybind/include/pybind11/pybind11.h>
namespace pynomial{
    namespace _register{ // NOTE: register is a keyword so we have to use the underscore.
        void export_register_module(pybind11::module& m)
        {
            pybind11::class_<RegisterBruteForce, std::shared_ptr<RegisterBruteForce> >(m,"BruteForce")
                .def(pybind11::init< const matrix&, double >())
                .def("Fit", &RegisterBruteForce::Fit)
                .def("GetRotation", &RegisterBruteForce::GetRotation)
                .def("GetCost", &RegisterBruteForce::GetCost)
                .def("SetNumShuffles", &RegisterBruteForce::SetNumShuffles)
            ;
        }
        void export_transform_module(pybind11::module& m)
        {
            m.def("quaternion_from_rotation_matrix", &make_quaternion_from_rotation_matrix);
        }
    }
}
