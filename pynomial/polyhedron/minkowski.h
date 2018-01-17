/*

Here we define a few functions to perform some basics of
convexity theory (support functions) and Minkowski theory (addition, subtraction,
and mixed volumes).

*/

#ifndef pynomial_MINKOWSKI
#define pynomial_MINKOWSKI

#include "Polyhedra.h"
namespace pynomial{
namespace polyhedron {
using std::vector;
using std::string;
using std::cout;
using std::endl;

struct support_function
{

    support_function(const Polyhedron& P) : verts(P.Vertices()){}
    // TODO: needs verification: the support function is maximized on a vertex.
    double operator () (const Eigen::Vector3d& v)
    {
        // h(P, v) = sup { a . v : a in P }
        Eigen::VectorXd dots = verts*v;
        return dots.maxCoeff();
    }
    Eigen::MatrixXd verts;
};

struct minkowski_sum
{
};

struct minkowski_sub
{
};

struct voulume_ABB
{
};


class SlopeDiagram
{
public:
    SlopeDiagram(const Polyhedron& P) : m_P(P)
    {
    }
    struct align
    {
    };
private:
    const Polyhedron& m_P;
    Eigen::MatrixXd m_SphericalPoints;              // correspond to the face normals
    vector< vector<int> > m_SphericalEdges;         // join faces of the polyhedra
    vector< vector<int> > m_SphericalFacets;        // the spherical polygons that make the faces
    std::vector<double> m_FacetArea;                // surface area of the facets
    std::vector<double> m_EdgeWeights;              // edge lengths
};




}} // namespaces
#endif // pynomial_MINKOWSKI
