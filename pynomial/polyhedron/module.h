#pragma once
#include "Polyhedra.h"
#include "Net.h"
#include "minkowski.h"
#include <pynomial/extern/pybind/include/pybind11/pybind11.h>
#include <pynomial/extern/pybind/include/pybind11/stl_bind.h>

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
                .def("GetDeterminant", &Polyhedron::GetDeterminant)
                .def("GetTrace", &Polyhedron::GetTrace)
                .def("FacetAreas", &Polyhedron::FacetAreas)
                .def("GetMergedFaces", &Polyhedron::GetMergedFaces)
                .def("GetMergedFacesMap", &Polyhedron::GetMergedFacesMap)
                .def("GetFacetNormals", &Polyhedron::GetFacetNormals)
                .def("Vertices", &Polyhedron::Vertices)
                .def("Curvature", &Polyhedron::Curvature)
                .def("MergeFacets", &Polyhedron::MergeFacets)
                .def("GetDihedrals", &Polyhedron::GetDihedrals)
                .def("Transform", &Polyhedron::Transform)
            ;
            pybind11::class_<CDihedralAngle, std::shared_ptr<CDihedralAngle> >(m,"_dihedral_angles")
                .def(pybind11::init<size_t, size_t, double>())
                .def("Face0", &CDihedralAngle::GetFace0)
                .def("Face1", &CDihedralAngle::GetFace1)
                .def("Angle", &CDihedralAngle::GetAngle)
                ;
            pybind11::bind_vector<CDihedralAngle>(m,"dihedral_angle_vector");

            pybind11::class_<SlopeDiagram, std::shared_ptr<SlopeDiagram> >(m,"_SDR")
                .def(pybind11::init<const Polyhedron& >())
                .def(pybind11::init<const SlopeDiagram& >())
            ;
            pybind11::class_<SimilarityResult, std::shared_ptr<SimilarityResult> >(m,"_SimilarityResult")
                .def(pybind11::init< >())
                .def_readonly("sigma", &SimilarityResult::sigma)
                .def_readonly("rotationP", &SimilarityResult::rotationP)
                .def_readonly("rotationQ", &SimilarityResult::rotationQ)
            ;

        }
    }
}
