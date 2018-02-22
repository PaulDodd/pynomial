/*

Here we define a few functions to perform some basics of
convexity theory (support functions) and Minkowski theory (addition, subtraction,
and mixed volumes).

*/

#ifndef pynomial_MINKOWSKI
#define pynomial_MINKOWSKI

#include <array>
#include "Polyhedra.h"
#include "../geometry/sphere.h"
namespace pynomial{
namespace polyhedron {
// using std::vector;
// using std::string;
// using std::cout;
// using std::endl;

template<class T>
inline std::pair<T,T> quad(T a, T b, T c)
{
    double s = b*b - 4.0*a*c;
    if( s < 0 )
    {
        throw std::runtime_error("complex is not supported.");
    }
    s = sqrt(s);
    return std::make_pair<T,T>((s-b)/2.0/a, -(b+s)/2.0/a);
}

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
    PointMatrix verts;
};

struct minkowski_sum
{
};

struct minkowski_sub
{
};

struct voulume_ABB
{
    double operator()(){
        return 0.0;
    }
};


class SlopeDiagram
{
public:
    typedef geometry::DefaultPointType point_t;
    typedef geometry::Sphere<geometry::DefaultPointType> sphere_t;
    typedef geometry::SphereArc<geometry::DefaultPointType, geometry::shallow> arc_t;
    typedef Eigen::Transform<double, 3, Eigen::TransformTraits::Affine> transform_t;
    typedef Eigen::Quaternion<double> quaternion_t;
    typedef enum{type1=0x01, type2=0x02} alignment_type;
public:

    SlopeDiagram(const Polyhedron& P) : m_P(P), m_Sphere(sphere_t::Unit())
    {
        // std::cout << " N - SDR " << std::endl;
        PointMatrix temp;
        m_P.FacetNormals(temp, true);
        // std::cout << " fn - SDR -  "<< temp.rows() << std::endl;
        m_SphericalPoints.reserve(temp.rows());
        for(int i = 0; i < temp.rows(); i++)
        {
            // std::cout << " pn - i -  "<< i << std::endl;
            point_t newp = temp.row(i);
            // newp[0] = temp(i, 0);
            // newp[1] = temp(i, 1);
            // newp[2] = temp(i, 2);
            // std::cout << " pn - SDR -  \n"<< newp << std::endl;
            m_SphericalPoints.push_back(newp);
        }
        // std::cout << " i - SDR -  "<< m_SphericalPoints.size() << std::endl;
        graph::CNetwork face_graph;
        m_P.GetMergedFaceGraph(face_graph);
        // std::cout << " i - SDR -  face off"<< std::endl;
        for(auto it = face_graph.begine(); it.IsValid(); it++)
        {
            // if(it->source >= m_SphericalPoints.size() || it->target >= m_SphericalPoints.size())
            //     std::cout << " i - SDR - 0 " << std::endl;
            arc_t arc(m_Sphere, m_SphericalPoints[it->source], m_SphericalPoints[it->target]);
            // std::cout << "d("<< it->source << ", " << it->target << ") = " << arc.Length() << std::endl;
            m_SphericalEdges.push_back(arc);
        }
        // std::cout << " X - SDR " << std::endl;
    }
    SlopeDiagram(const SlopeDiagram& src) : m_P(src.GetP()), m_Sphere(sphere_t::Unit()), m_SphericalPoints(src.GetPoints())
    {
        // std::cout << "Copy SDR!" << std::endl;
        graph::CNetwork face_graph;
        m_P.GetMergedFaceGraph(face_graph);
        for(auto it = face_graph.begine(); it.IsValid(); it++)
        {
            arc_t arc(m_Sphere, m_SphericalPoints[it->source], m_SphericalPoints[it->target]);
            m_SphericalEdges.push_back(arc);
        }
        // std::cout << "Copy SDR! - out" << std::endl;
    }

    void Transform(const transform_t& t)
    {
        // PointMatrix temp = m_SphericalPoints.transpose();
        // std::cout << t.rotation() << std::endl;
        for(int i = 0; i < m_SphericalPoints.size(); i++)
            {
            // std::cout << "pre: \n" << m_SphericalPoints[i].norm() << std::endl;
            m_SphericalPoints[i] = (t.rotation() * m_SphericalPoints[i]).eval();
            // std::cout << "post: \n" << m_SphericalPoints[i].norm() << std::endl;
            }
    }

    const point_t& GetPoint(const int& i) const { return m_SphericalPoints[i]; }

    const Polyhedron& GetP() const { return m_P; }

    void writePOS(const string& filename, std::ios_base::openmode mode=std::ios_base::out ) const
    {
        std::ofstream file(filename, mode);
        std::stringstream ss, connections;
        ss << "def "<< "sdr" << " \"poly3d " <<  m_SphericalPoints.size() << " ";
        std::string center_sphere  = "def origin \"sphere 0.1 005F5F5F\"";
        for(int i = 0; i < m_SphericalPoints.size(); i++)
            ss << m_SphericalPoints[i][0] << " " << m_SphericalPoints[i][1] << " " << m_SphericalPoints[i][2] << " ";
        ss << "005984FF\"";
        std::string hull  = ss.str();
        file << hull << std::endl;
        file << center_sphere << std::endl;
        file << "sdr" << " 0 0 0 1 0 0 0" << std::endl;
        file << "origin 0 0 0 "<< std::endl;
        file << "eof" << std::endl;
        file.close();
    }

    unsigned int NumPoints() const { return m_SphericalPoints.size(); }
    unsigned int NumArcs() const { return m_SphericalEdges.size(); }

public:
    // Here I make an iterator instead of a full enumerative method to give more
    // flexibility to the user. e.g., they can decide when to exit.
    // do you want to check both types of minimizations?
    class alignment_iterator
    {
        struct _type_1_state{
            std::pair<int,int> p;  // points in P
            std::pair<int,int> q;  // points in Q
            std::pair<int,int> q2; // points in Q
            std::pair<int,int> a;  // arcs in P
            transform_t rot_p;
            bool b_rot_p;
            transform_t rot_q;
            bool b_rot_q;
            bool end() const
            {
                return  p.first == p.second;
            }
            _type_1_state& operator++()
            {
                bool inc = true;
                while(p.first < p.second)
                {
                    while(q.first < q.second)
                    {
                        while(q2.first < q2.second)
                        {
                            if(q.first == q2.first)
                            {
                                q2.first++;
                                continue;
                            }

                            while(a.first+1 < a.second)
                            {
                                if(inc)
                                {
                                    b_rot_q = false;
                                    b_rot_p = false;
                                    ++a.first;
                                }
                                return *this;
                            }
                            a.first = 0;
                            inc = false;
                            ++q2.first;
                        }
                        b_rot_q = true;
                        q2.first = 0;
                        ++q.first;
                    }
                    b_rot_p = true;
                    q.first = 0;
                    ++p.first;
                }
                return *this;
            }
            void print() const { std::cout << "p: " << p.first << ", q: " << q.first << ", q2: " << q2.first << ", a: " << a.first << std::endl; }

            // _type_1_state operator++(int)
            // {
            //     // copy *this. then advance the state then return the copy.
            //     _type_1_state it(*this);
            //     advance_type_1();
            //     return it;
            // }

        };
        struct _type_2_state{
            /* No implementation yet */
        };

        struct _state_info{
            _state_info(SlopeDiagram& _P, SlopeDiagram& _Q, alignment_type typ) : atype(typ), P(_P), Q(_Q)
            {
                rotP.setIdentity();
                rotQ.setIdentity();
                // Type 1 state
                state1.p = std::make_pair(0, P.NumPoints());
                state1.q = std::make_pair(0, Q.NumPoints());
                state1.q2 = std::make_pair(1, Q.NumPoints());
                state1.a = std::make_pair(0, P.NumArcs());
                state1.rot_p.setIdentity();
                state1.b_rot_p = true;
                state1.rot_q.setIdentity();
                state1.b_rot_q = true;
            }

            _state_info(const _state_info& src) : atype(src.atype), P(src.P), Q(src.Q)
            {
            }

            alignment_type atype;
            SlopeDiagram& P;
            SlopeDiagram& Q;
            transform_t rotQ; // continuously keep track of the rotation from the original Q to Q'
            transform_t rotP; // continuously keep track of the rotation from the original P to P'
            _type_1_state state1;
        } state_info;
    public:
        alignment_iterator(SlopeDiagram& P, SlopeDiagram&Q, alignment_type typ=type1) : state_info(P, Q, typ), bfirst(true)
        {
        }
    private:
        /*
        Mixed volume is minimzed by the rotation that is produced by either
        1.  It is fulfilled simultaneously that one spherical point of SDR(Q')
            coincides with a spherical point of SDR(P) and another spherical
            point of SDR(Q') intersects a spherical arc of SDR(P);

        2.  Three spherical points of SDR(Q') intersect three spherical arcs
            of SDR(P) simultaneously.
        */
        void advance_type_1()
        {
            // pseudo-code:
            //  for p in P(points)
            //      rot_p = {p -> (0,0,1)}*rot_p;
            //      for q in Q(points)
            //          rot2 = q -> p (the rotation that rotates point q to p)
            //          for q2 in Q(points)
            //              for a in P(arcs)
            //                  if( phi_min < phi(q2) < phi_max )
            //                      // found a critical rotation!
            //                      theta = q2 |-> a (angle about the z axis to make q2 intersect a)
            //                      rot3 = AngleAxis(theta, (0,0,1))
            //                      rot_q = {rot3*rot2}*rot_q;

            _type_1_state& state = state_info.state1;
            static const Eigen::Vector3d zhat = {0.0, 0.0, 1.0};
            if(!bfirst) // on the first sweep we need to check the initial state.
                ++state;
            bfirst = false;
            brotp = false;
            brotq = false;
            while(!state.end())
            {
                if(state.b_rot_p)
                {
                    // The transform from the current state to the next state.
                    state.rot_p = transform_t(quaternion_t().FromTwoVectors(state_info.P.GetPoint(state.p.first), zhat));
                    // std::cout << "transform : \n" << state.rot_p.linear() << std::endl;
                    // The transformation from the original state to the current state.
                    state_info.rotP = state.rot_p*state_info.rotP;
                    const arc_t& arc = state_info.P.GetEdges()[0];

                    state_info.P.writePOS("testP.pos", std::ios_base::app);
                    point_t normal = point_t::Zero();
                    arc.PlaneNormal(normal);
                    // std::cout << "n - Normal: \n" << normal << std::endl;
                    state_info.P.Transform(state.rot_p);
                    brotp = true;
                    // std::cout << "p - point: \n" << &state_info.P.GetPoint(state.p.first) << std::endl;
                    arc.PlaneNormal(normal);
                    // std::cout << "x - Normal: \n" << normal << std::endl;
                    state_info.P.writePOS("testP.pos", std::ios_base::app);
                }
                if(state.b_rot_q)
                {
                    // The transform from the current state to the next state.
                    state.rot_q = transform_t(quaternion_t().FromTwoVectors(state_info.Q.GetPoint(state.q.first), zhat));
                    // The transformation from the original state to the current state.
                    state_info.rotQ = state.rot_q*state_info.rotQ; // TODO: is there aliasing?
                    state_info.Q.writePOS("testQ.pos", std::ios_base::app);
                    // std::cout << "n - Q point\n" << state_info.Q.GetPoint(state.q.first) << std::endl;
                    state_info.Q.Transform(state.rot_q);
                    brotq = true;
                    // std::cout << "x - Q point\n" << state_info.Q.GetPoint(state.q.first) << std::endl;
                    state_info.Q.writePOS("testQ.pos", std::ios_base::app);
                }
                // state.print();
                const arc_t& arc = state_info.P.GetEdges()[state.a.first];
                const point_t& end1 = arc.GetP1();
                const point_t& end2 = arc.GetP2();
                double zmin  = std::min(end1[2], end2[2]);
                double zmax  = std::max(end1[2], end2[2]);
                point_t point  = state_info.Q.GetPoint(state.q2.first);

                if(zmin < point[2] && zmax > point[2])
                {
                    // std::cout << "arc 1  = \n" << end1 << std::endl;
                    // std::cout << "arc 2  = \n" << end2 << std::endl;
                    // std::cout << "found critical rotation: "<< zmin << " "<< point[2] << " "<< zmax << std::endl;
                    point_t normal = point_t::Zero();
                    point_t point_arc = {0.0,0.0,point[2]};
                    point_t point_arc2 = {0.0,0.0,point[2]};
                    arc.PlaneNormal(normal);
                    double a, b, c, x1, x2, y1, y2;
                    if(fabs(normal[0]) < 1e-4 && fabs(normal[1]) < 1e-4)
                    {
                        // std::cout << "found critical rotation: "<< zmin << " "<< point[2] << " "<< zmax << std::endl;
                        // std::cout << "Normal: \n" << normal << std::endl;
                        // throw std::runtime_error("infinite solutions found");
                        x1 = end1[0];
                        y1 = end1[1];
                        x2 = end2[0];
                        y2 = end2[1];
                    }
                    else if( fabs(normal[0]) < 1e-4 )
                    {
                        x1 = sqrt(1 - (point[2]*point[2]*(1-(normal[2]*normal[2]/normal[1]/normal[1]))));
                        x2 = -x1;
                        y1 = -normal[2]*point[2]/normal[1];
                        y2 = y1;
                    }
                    else if(fabs(normal[1]) < 1e-4)
                    {
                        x1 = -normal[2]*point[2]/normal[0];
                        x2 = x1;
                        y1 = sqrt(1 - (point[2]*point[2]*(1-(normal[2]*normal[2]/normal[0]/normal[0]))));
                        y2 = -y1;
                    }
                    else
                    {
                        a = (normal[1]*normal[1]/normal[0]/normal[0]) + 1;
                        b = 2.0*normal[1]*normal[2]*point[2]/normal[0]/normal[0];
                        c = (((normal[2]*normal[2]/normal[0]/normal[0])+1)*point[2])+1.0;
                        if( b*b < 4.0*a*c) // not sure what to do here.
                        {
                            ++state;
                            continue;
                        }
                        auto result = quad(a,b,c);
                        x1 = -(result.first*normal[1] + point[2]*normal[2])/normal[0];
                        x2 = -(result.second*normal[1] + point[2]*normal[2])/normal[0];
                        y1 = result.first;
                        y2 = result.second;
                    }

                    point_arc[0] = x1;
                    point_arc[1] = y1;

                    point_arc2[0] = x2;
                    point_arc2[1] = y2;

                    if(arc.On(point_arc))
                    {
                        point_t sp = geometry::Cartesian2Spherical(point);
                        point_t spa = geometry::Cartesian2Spherical(point_arc);
                        double delta_theta = spa[1] - sp[1];
                        // std::cout << "delta_theta: " << delta_theta << std::endl;
                        state.rot_q = transform_t(Eigen::AngleAxis<double>(delta_theta, zhat));
                    }
                    else if(arc.On(point_arc2))
                    {
                        point_t sp = geometry::Cartesian2Spherical(point);
                        point_t spa = geometry::Cartesian2Spherical(point_arc);
                        double delta_theta = spa[1] - sp[1];
                        // std::cout << "delta_theta: " << delta_theta << std::endl;
                        state.rot_q = transform_t(Eigen::AngleAxis<double>(delta_theta, zhat));
                    }
                    else{

                        // std::cout << "Neither result is on the arc? logical error? \n ******************** \n"  << std::endl;
                        // arc.On(point_arc, true);
                        // arc.On(point_arc2, true);
                        // std::cout << "\n ******************** \n"  << std::endl;
                        ++state; // not sure what to do here.
                        continue;
                    }


                    // The transform from the current state to the next state.
                    // quaternion_t().FromTwoVectors(state_info.Q.GetPoint(state.q.first), zhat)
                    // state.rot_q = transform_t();
                    // The transformation from the original state to the current state.
                    state_info.rotQ = state.rot_q*state_info.rotQ;
                    state_info.Q.writePOS("testQ.pos", std::ios_base::app);
                    // std::cout   << "point 1: " << point_arc[0] << ", "<< point_arc[1] << ", " << point_arc[2]
                    //             << "\npoint 2: "<< point_arc2[0] << ", "<< point_arc2[1] << ", " << point_arc2[2]
                    //             << std::endl;
                    // std::cout << "******* n - Q point\n" << state_info.Q.GetPoint(state.q2.first) << std::endl;
                    state_info.Q.Transform(state.rot_q);
                    brotq = true;
                    // std::cout << "******* x - Q point\n" << state_info.Q.GetPoint(state.q2.first) << std::endl;
                    // std::cout << "******* x - Q ref point\n" << state_info.Q.GetPoint(state.q.first) << std::endl;
                    // state_info.Q.writePOS("testQ.pos", std::ios_base::app);
                    return;
                }
                ++state;
            }
        }

        void advance_type_2(SlopeDiagram& P, SlopeDiagram& Q, std::vector<transform_t>& Transforms )
        {
            // This one is a bit more tricky
            //

        }
    public:
        alignment_iterator& operator++()
        {
            // advance the state then return a refernce to *this.
            advance_type_1();
            return *this;
        }

        alignment_iterator operator++(int)
        {
            // copy *this. then advance the state then return the copy.
            alignment_iterator it(*this);
            advance_type_1();
            return it;
        }

        bool end() const { return state_info.state1.end(); }

        bool is_rotationP() const { return brotp; }

        bool is_rotationQ() const { return brotq; }

        const transform_t& rotationP() const { return state_info.rotP; }

        const transform_t& rotationQ() const { return state_info.rotQ; }

    private:
        bool bfirst;
        bool brotp;
        bool brotq;
    };

    const std::vector< arc_t >& GetEdges() const { return m_SphericalEdges; }
    const std::vector< point_t >& GetPoints() const { return m_SphericalPoints; }
private:
    const Polyhedron& m_P;
    sphere_t m_Sphere;
    std::vector< point_t >  m_SphericalPoints;                      // correspond to the face normals
    std::vector< arc_t >    m_SphericalEdges;              // join faces of the polyhedra
    // vector< vector<int> > m_SphericalFacets;            // the spherical polygons that make the faces
    // std::vector<double> m_FacetArea;                    // surface area of the facets
    // std::vector<double> m_EdgeWeights;                  // edge lengths
};

struct SimilarityResult
{
    double sigma;
    Eigen::Matrix3d rotationP;
    Eigen::Matrix3d rotationQ;
};

SimilarityResult SimilarityMeasure(/*SimilarityResult& result,*/
                        Polyhedron A,
                        Polyhedron B,
                        float tolerance=0.95,
                        int iteration_max=-1)
{
    int count = 0;
    double sigma = 0.0, sigma_max = 0.0;
    SlopeDiagram P(A);
    SlopeDiagram Q(B);
    double VA = pow(A.Volume(), (1.0/3.0));
    double VB = pow(B.Volume(), (2.0/3.0));
    double VABB = 1.0;
    support_function hA(A);
    Eigen::VectorXd SA = B.FacetAreas(true);
    Eigen::Matrix3d rotA;
    Eigen::Matrix3d rotB;
    PointMatrix vertsA = A.Vertices().transpose();
    for(SlopeDiagram::alignment_iterator it(P,Q); !it.end(); ++it )
    {
        if(it.is_rotationP())
        {
            hA.verts = (it.rotationP().rotation()*vertsA).eval().transpose();
        }
        Eigen::VectorXd sup(Q.NumPoints());
        for(int i = 0; i < Q.NumPoints(); i++)
        {
            sup[i] = hA(Q.GetPoint(i));
        }
        VABB = sup.dot(SA)/3.0;
        sigma = VA*VB/VABB;
        if(sigma > sigma_max)
        {
            sigma_max = sigma;
            rotA = it.rotationP().rotation();
            rotB = it.rotationQ().rotation();
        }

        ++count;
        if((iteration_max > 0 && count > iteration_max) || sigma_max > tolerance)
            break;
    }
    // std::cout << "tried "<< count << " rotations!" << std::endl;
    // std::cout << "sigma max = " << sigma_max <<std::endl;
    // result.reset(new SimilarityResult);
    // result.sigma = static_cast<float>(sigma_max);
    // result.rotationP = rotA;
    // result.rotationQ = rotB;
    return {sigma_max, rotA, rotB};
}


}} // namespaces
#endif // pynomial_MINKOWSKI
