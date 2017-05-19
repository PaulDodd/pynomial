#pragma once
#include "Network.h"
#include "NN.h"
#include <pynomial/extern/pybind/include/pybind11/pybind11.h>
namespace pynomial{
    namespace graph{
        void export_graph_module(pybind11::module& m)
        {
            pybind11::class_<CNetwork, std::shared_ptr<CNetwork> >(m,"graph")
                .def(pybind11::init<int, size_t>())
                .def("add_edge", &CNetwork::AddEdge)
                .def("remove_edge", &CNetwork::RemoveEdge)
                .def("is_edge", &CNetwork::DoesEdgeExist)
                .def("neighbors", &CNetwork::GetNeighbors)
                .def("path_connected", &CNetwork::AreNodesConnected)
                .def("spanning_tree", &CNetwork::KruskalAlgorithm)
            ;
        }

        void export_NN_module(pybind11::module& m)
        {
            pybind11::class_<NN, std::shared_ptr<NN> >(m,"NN")
                .def(pybind11::init<const NN::matrix_t&, unsigned int, double>())
                .def("match", &NN::match)
                .def("num_points", &NN::num_points)
                .def("pair_distance", &NN::pair_distance)
            ;
        }
    }
}
