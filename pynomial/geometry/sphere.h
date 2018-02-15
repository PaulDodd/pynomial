
#ifndef pynomial_Sphere_h
#define pynomial_Sphere_h

#include <vector>
#include <string>
#include <functional>
#include <algorithm>
#include <iostream>
#include <type_traits>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Sparse>

// Sepcial note for now most of these problems are solved only for the 3D case
// I am leaving the higher dimensional counter parts for later.

namespace pynomial{
namespace geometry{

enum copy_type{shallow, deep};

typedef Eigen::Vector3d DefaultPointType;

template<class PointType>
struct norm_fn : std::function<double(const PointType&)>
{
    double operator()(const PointType& p) const { return p.norm(); }
};

template<class PointType>
struct norm2_fn : std::function<double(const PointType&)>
{
    double operator()(const PointType& p) const { return p.squaredNorm(); }
};


template<class PointType>
struct dot_fn : std::function<double(const PointType&,const PointType&)>
{
    double operator()(const PointType& a,const PointType& b) const { return a.dot(b); }
};

template<class PointType>
struct cross_fn : std::function<PointType(const PointType&, const PointType&)>
{
    PointType operator()(const PointType& a, const PointType& b) const { return a.cross(b); }
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
    double GetRadius() const { return m_Radius; }
    void SetRadius(double r) { m_Radius = r; }
// Helpers
    PointType project(const PointType& v) const { return (m_Radius / norm(v)) * v; }
    // TODO: other useful methods
    // -- uniform points on a sphere
    // -- random uniform points on a sphere
    // -- polygons and geometry?

    static Sphere<PointType> Unit() { return Sphere<PointType>(1.0); }
private:
    double m_Radius;
    unsigned int m_dim;
};

template<class PointType=DefaultPointType, copy_type copy = deep>
class SphereArc
{
    typedef typename std::conditional<copy == deep, PointType, const PointType&>::type point_t;
    norm_fn<PointType> norm;
    dot_fn<PointType> dot;
    cross_fn<PointType> cross;
public:
    SphereArc(const Sphere<PointType>& s, const PointType& p1, const PointType& p2) : m_Sphere(s), m_p1(p1), m_p2(p2)
    {
        // assert(p1.rows() == dim &&  p2.rows() == dim);
        // assert(s.dim() == 3);
        // std::is_same<point_t, const PointType&>::value;
        // typeid(point_t).name();
        // std::cout << "copy type: " << typeid(point_t).name() << " " << std::boolalpha << std::is_same<point_t, const PointType&>::value << std::endl;
    }

// parameterizations
    // s(t) = c1(t)*p1 + c2(t)*p2
    PointType Slerp(double t)
    {
        return DefaultPointType::Zero();
    }

    // theta(phi)
    // double Theta(double phi);

    // phi(theta)
    // double Phi(double theta);

// properties
    bool PlaneNormal(Eigen::Vector3d& n) const
    {
        // equation of the plane that intersects the sphere is of the form
        // < n, x > = 0
        // where n = cross(p1, p2) / r ^ 2
        double r = m_Sphere.GetRadius();
        // std::cout << std::hex << "p1:\n" << &m_p1 << "\np2:\n" << m_p2 << std::endl;
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

    double Length() const
    {
        double r = m_Sphere.GetRadius();
        return r * acos(dot(m_p1,m_p2)/r/r);
    }

    double MinPhi() const
    {
        return std::min(Cartesian2Spherical(m_p1)[2], Cartesian2Spherical(m_p2)[2]);
    }

    double MaxPhi() const
    {
        return std::max(Cartesian2Spherical(m_p1)[2], Cartesian2Spherical(m_p2)[2]);
    }

    double MinTheta() const
    {
        return std::min(Cartesian2Spherical(m_p1)[1], Cartesian2Spherical(m_p2)[1]);
    }

    double MaxTheta() const
    {
        return std::max(Cartesian2Spherical(m_p1)[1], Cartesian2Spherical(m_p2)[1]);
    }

    const PointType& GetP1() const { return m_p1; }

    const PointType& GetP2() const { return m_p2; }

    bool On(const PointType& p, bool verbose=false) const
    {
        double r = m_Sphere.GetRadius();
        PointType ps = m_Sphere.project(p)/r;
        double d1 = r*acos(dot(ps, m_p1)/r);
        double d2 = r*acos(dot(ps, m_p2)/r);
        if(verbose)
        {
            std::cout << r << ", "<< d1 << ", "<< d2 << ", "<< Length() << "\np1:\n" << m_p1 << "\np2:\n" << m_p2 << "\n ps:\n" << ps << std::endl;
            std::cout<< "trangle check: " << fabs(Length() - d1 - d2) << std::endl;
        }
        return fabs(Length() - d1 - d2) < 1e-4; // triangle inequality.
    }

private:
    const Sphere<PointType>& m_Sphere;
    point_t m_p1;
    point_t m_p2;
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
