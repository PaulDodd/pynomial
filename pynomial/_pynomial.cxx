
#include "pynomial.h"
#include <pynomial/extern/pybind/include/pybind11/pybind11.h>
#include <pynomial/extern/pybind/include/pybind11/stl_bind.h>
#include <pynomial/extern/pybind/include/pybind11/eigen.h>
#include "extern/num_util.h"


using namespace std;
using namespace pynomial;

void* py3_import_array()
{
    import_array();
    return NULL;
}

template<class T>
void safe_bind_vector(pybind11::module& m, const char * name)
{
    try{
        pybind11::bind_vector<T>(m,name);
    }catch(std::runtime_error err)
    {
        std::cerr << "warning: generic_type named \'" << name << "\' is already registered" << std::endl;
    }
}
PYBIND11_PLUGIN(_pynomial)
{
    pybind11::module m("_pynomial");

    // setup needed for numpy
    py3_import_array();

    m.attr("__version__") = pybind11::make_tuple(0, 0, 0);


    safe_bind_vector<float>(m,"std_vector_float");
    safe_bind_vector<string>(m,"std_vector_string");
    safe_bind_vector<unsigned int>(m,"std_vector_uint");
    safe_bind_vector<unsigned long>(m,"std_vector_ul");
    safe_bind_vector<int>(m,"std_vector_int");
    // safe_bind_vector<Eigen::VectorXd>(m, 'std_vector_eigen_vectorXd');
    // safe_bind_vector<Eigen::Vector3d>(m, 'std_vector_eigen_vector3d');
    // pybind11::bind_vector<Scalar3>(m,"std_vector_scalar3");
    // pybind11::bind_vector<Scalar4>(m,"std_vector_scalar4");

    // InstallSIGINTHandler();

    // data structures
    return m.ptr();
}
