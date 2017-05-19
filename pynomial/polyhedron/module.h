#pragma once
#include "Polyhedra.h"
#include "Net.h"
#include <pynomial/extern/pybind/include/pybind11/pybind11.h>
namespace pynomial{
    namespace polyhedron{
        void export_polyhedron_module(pybind11::module& m)
        {
            pybind11::class_<Polyhedron, std::shared_ptr<Polyhedron> >(m,"polyhedron")
                .def(pybind11::init< >())
                .def("Set", (void (Polyhedron::*)(Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>)) &Polyhedron::Set) //(void (Polyhedron::*)(const vector< Eigen::Vector3d >&))
                // .def("Set", (void (Polyhedron::*)(const vector< Eigen::Vector3d >&,const vector< vector<int> >&)) &Polyhedron::Set)
                .def("Dual", &Polyhedron::Dual)
                .def("writePOS", &Polyhedron::writePOS)
                .def("Volume", &Polyhedron::Volume)
                .def("Centroid", &Polyhedron::Centroid)
                .def("Vertices", &Polyhedron::Vertices)
            ;
        }
    }
}
