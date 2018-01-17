
#ifndef pynomial_Polyhedra_h
#define pynomial_Polyhedra_h

#include <vector>
#include <string>
#include <functional>
#include <algorithm>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Sparse>

// Sepcial note for now most of these problems are solved only for the 3D case
// I am leaving the higher dimensional counter parts for later.

namespace pynomial{
namespace geometry{

typedef Eigen::VectorXd DefaultPointType;

template<class PointType>
struct norm_fn : std::function<double(const PointType&)>
{
    double operator()(const PointType& p) { return p.norm(); }
};

template<class PointType>
struct norm2_fn : std::function<double(const PointType&)>
{
    double operator()(const PointType& p) { return p.squaredNorm(); }
};


template<class PointType>
struct dot_fn : std::function<double(const PointType&,const PointType&)>
{
    double operator()(const PointType& a,const PointType& b) { return a.dot(b); }
};

template<class PointType>
struct cross_fn : std::function<PointType(const PointType&, const PointType&)>
{
    PointType operator()(const PointType& a, const PointType& b) { return a.cross(b); }
};


template<class PointType=DefaultPointType>
Eigen::Vector3d Cartesian2Spherical(const PointType& v)
{
    norm_fn<PointType> norm;
    Eigen::Vector3d sph; // r, theta, phi
    sph[0] = norm(v);
    sph[1] = atan2(v[1], v[0]);
    sph[2] = acos(v[2]/sph[0]);
    return sph;
}

template<class PointType=DefaultPointType>
class Sphere
{
    norm_fn<PointType> norm;
public:
// Constrctors
    Sphere(double radius = 1.0) : m_Radius(radius), m_dim(3) {} // TODO: make m_dim a parameter
    Sphere( const Sphere& s) : m_Radius(s.GetRadius()){}
// Accessors
    unsigned int dim() {return m_dim;}
    double GetRadius() { return m_Radius; }
    void SetRadius(double r) { m_Radius = r; }
// Helpers
    PointType project(const PointType& v) { return (m_Radius / norm(v)) * v; }
    // TODO: other useful methods
    // -- uniform points on a sphere
    // -- random uniform points on a sphere
    // -- polygons and geometry?
private:
    double m_Radius;
    unsigned int m_dim;
};

enum copy_type{shallow, deep};

template<class PointType=DefaultPointType, copy_type copy = deep>
class SphereArc
{
    norm_fn<PointType> norm;
    dot_fn<PointType> dot;
    cross_fn<PointType> cross;
public:
    SphereArc(const Sphere<PointType>& s, const PointType& p1, const PointType& p2)
    {
        assert(p1.rows() == dim &&  p2.rows() == dim);
        assert(s.dim() == 3);
        m_Sphere = s;
        #if copy == deep
        m_p1 = m_Sphere.project(p1);
        m_p2 = m_Sphere.project(p2);
        #else
        m_p1 = p1;
        m_p2 = p2;
        #endif
    }

// parameterizations
    // s(t) = c1(t)*p1 + c2(t)*p2
    // PointType Slerp(double t);

    // theta(phi)
    // double Theta(double phi);

    // phi(theta)
    // double Phi(double theta);

// properties
    bool PlaneNormal(Eigen::VectorXd& n)
    {
        // equation of the plane that intersects the sphere is of the form
        // < n, x > = 0
        // where n = cross(p1, p2) / r ^ 2
        double r = m_Sphere.GetRadius();
        n = cross(m_p1/r, m_p2/r);
        double mag = norm(n);
        // TODO: make tolerance a parameter
        if(mag < 1e-4) // the antipodal case retun
        {
            n = Eigen::Vector3d::Zero();
            return false;
        }
        n.normalize();
        return true;
    }

    double Length()
    {
        double r = m_Sphere.GetRadius();
        return r * acos(dot(m_p1,m_p2)/r/r);
    }

    double MinPhi()
    {
        return min(Cartesian2Spherical(m_p1)[2], Cartesian2Spherical(m_p2)[2]);
    }

    double MaxPhi()
    {
        return max(Cartesian2Spherical(m_p1)[2], Cartesian2Spherical(m_p2)[2]);
    }

    double MinTheta()
    {
        return min(Cartesian2Spherical(m_p1)[1], Cartesian2Spherical(m_p2)[1]);
    }

    double MaxTheta()
    {
        return max(Cartesian2Spherical(m_p1)[1], Cartesian2Spherical(m_p2)[1]);
    }

private:
    const Sphere<PointType>& m_Sphere;
    #if copy == deep
    PointType m_p1;
    PointType m_p2;
    #else
    const PointType& m_p1;
    const PointType& m_p2;
    #endif
};

// // This class is the composition of spherical points, arcs, and polygons.
// // the purpose of this class is to help calculate the
// template<class PointType=DefaultPointType>
// class SphereComplex
// {
// public:
//     // public members and
//     sphere(double radius = 1.0) : m_Radius(radius) {}
//
// private:
//     Sphere m_sphere;
//     std::vector<>
//
// };

}} // namespaces


#endif
