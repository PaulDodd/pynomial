//
//  Polyhedra.h
//  Created by Paul M Dodd on 2/22/13.
//

#ifndef pynomial_Polyhedra_h
#define pynomial_Polyhedra_h

#include <vector>
#include <string>
#include <algorithm>

// #include "SharedInclude.h"
// #include "json_wrapper.h"
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Sparse>
#include <chull/chull.h>
#include "pynomial/graph/Network.h"
#include "pynomial/graph/utils.h"
#include "pynomial/register/brute_force.h"
// #include "../geometry/geometry.h"
//#include "StatDist.h"

namespace pynomial{
namespace polyhedron {
typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> PointMatrix;
using std::vector;
using std::string;
using std::cout;
using std::endl;

class Polyhedron;

class mass_properties_base
{
public:
    mass_properties_base() : m_volume(0.0), m_center_of_mass(0.0, 0.0, 0.0)
        {
        for(unsigned int i = 0; i < 6; i++) m_inertia[i] = 0.0;
        }

    double getVolume() { return m_volume; }

    const Eigen::Vector3d& getCenterOfMass() { return m_center_of_mass; }

    double getCenterOfMassElement(unsigned int i)
        {
        if(i == 0 )
            return m_center_of_mass[0];
        else if(i == 1 )
            return m_center_of_mass[1];
        else if (i == 2)
            return m_center_of_mass[2];
        else
            throw std::runtime_error("index out of range");
        }

    double getInertiaTensor(unsigned int i)
        {
        if(i >= 6 )
            throw std::runtime_error("index out of range");
        return m_inertia[i];
        }

    double getDeterminant()
        {
        Eigen::Vector3d a(m_inertia[0], m_inertia[3], m_inertia[5]), b(m_inertia[3], m_inertia[1], m_inertia[4]), c(m_inertia[5], m_inertia[4], m_inertia[2]);
        return a.dot(b.cross(c));
        }

protected:
    virtual void compute() {throw std::runtime_error("mass_properties::compute() is not implemented for this shape.");}
    double m_volume;
    Eigen::Vector3d m_center_of_mass;
    double m_inertia[6]; // xx, yy, zz, xy, yz, xz
};

template<class Shape>
class mass_properties : public mass_properties_base
{
public:
    mass_properties() : mass_properties_base() {}
};

template <>
class mass_properties< Polyhedron > : public mass_properties_base
{
using mass_properties_base::m_volume;
using mass_properties_base::m_center_of_mass;
using mass_properties_base::m_inertia;

public:
    mass_properties(const Polyhedron& poly);

    mass_properties(const std::vector< Eigen::Vector3d >& p, const std::vector<std::vector<int> >& f) :  points(p), faces(f)
        {
        compute();
        }

    unsigned int getFaceIndex(unsigned int i, unsigned int j) { return faces[i][j]; }

    unsigned int getNumFaces() { return faces.size(); }

protected:
/*
    algorithm taken from
    http://www.geometrictools.com/Documentation/PolyhedralMassProperties.pdf
    NOTE: The moment of inertia tensor is calculated relative to the center of mass of the polyhedron.
*/
    virtual void compute()
        {
        // const std::vector<std::vector<unsigned int> >& faces = convex_hull.getFaces();
        // const std::vector< Eigen::Vector3d >& points = convex_hull.getPoints();
        const double mult[10] = {1.0/6.0 ,1.0/24.0 ,1.0/24.0 ,1.0/24.0 ,1.0/60.0 ,1.0/60.0 ,1.0/60.0 ,1.0/120.0 ,1.0/120.0 ,1.0/120.0};
        double intg[10] = {0,0,0,0,0,0,0,0,0,0}; // order: 1, x, y, z, xˆ2, yˆ2, zˆ2, xy, yz, zx
        for (unsigned int t=0; t<faces.size(); t++)
            {
            //get vertices of triangle
            Eigen::Vector3d v0, v1, v2;
            Eigen::Vector3d a1, a2, d;
            v0 = points[faces[t][0]];
            v1 = points[faces[t][1]];
            v2 = points[faces[t][2]];
            // get edges and cross product of edges
            a1 = v1 - v0;
            a2 = v2 - v0;
            d = a1.cross(a2);

            Eigen::Vector3d temp0, temp1, temp2, f1, f2, f3, g0, g1, g2;
            temp0 = v0 + v1;
            f1 = temp0 + v2;
            temp1 = v0.cwiseProduct(v0);
            temp2 = temp1 + v1.cwiseProduct(temp0);
            f2 = temp2 + v2.cwiseProduct(f1);
            f3 = v0.cwiseProduct(temp1) + v1.cwiseProduct(temp2) + v2.cwiseProduct(f2);
            g0 = f2 + v0.cwiseProduct(f1 + v0);
            g1 = f2 + v1.cwiseProduct(f1 + v1);
            g2 = f2 + v2.cwiseProduct(f1 + v2);

            intg[0] += d[0]*f1[0];
            intg[1] += d[0]*f2[0]; intg[2] += d[1]*f2[1]; intg[3] += d[2]*f2[2];
            intg[4] += d[0]*f3[0]; intg[5] += d[1]*f3[1]; intg[6] += d[2]*f3[2];
            intg[7] += d[0]*(v0[1]*g0[0] + v1[1]*g1[0] + v2[1]*g2[0]);
            intg[8] += d[1]*(v0[2]*g0[1] + v1[2]*g1[1] + v2[2]*g2[1]);
            intg[9] += d[2]*(v0[0]*g0[2] + v1[0]*g1[2] + v2[0]*g2[2]);
            }
        for(unsigned int i = 0; i < 10; i++ )
            {
            intg[i] *= mult[i];
            }

        m_volume = intg[0];

        m_center_of_mass[0] = intg[1];
        m_center_of_mass[1] = intg[2];
        m_center_of_mass[2] = intg[3];
        m_center_of_mass /= m_volume;

        double cx2 = m_center_of_mass[0]*m_center_of_mass[0], cy2 = m_center_of_mass[1]*m_center_of_mass[1], cz2 = m_center_of_mass[2]*m_center_of_mass[2];
        double cxy = m_center_of_mass[0]*m_center_of_mass[1], cyz = m_center_of_mass[1]*m_center_of_mass[2], cxz = m_center_of_mass[0]*m_center_of_mass[2];
        m_inertia[0] = intg[5] + intg[6] - m_volume*(cy2 + cz2);
        m_inertia[1] = intg[4] + intg[6] - m_volume*(cz2 + cx2);
        m_inertia[2] = intg[4] + intg[5] - m_volume*(cx2 + cy2);
        m_inertia[3] = -(intg[7] - m_volume*cxy);
        m_inertia[4] = -(intg[8] - m_volume*cyz);
        m_inertia[5] = -(intg[9] - m_volume*cxz);
        }
private:
    const std::vector< Eigen::Vector3d >& points;
    const std::vector<std::vector<int> >& faces;
};

// class cwless
// {
//     const Eigen::Vector3d& reference;
//     const std::vector<Eigen::Vector3d>& pointset;
//     int refi;
// public:
//     cwless(const Eigen::Vector3d& ref, const std::vector<Eigen::Vector3d>&  points, int r) : reference(ref), pointset(points), refi(r){}
//     bool operator()(int a, int b)
//     {
//         // A few cases to consider first
//         if(a == b)
//             return false;
//         if(a == refi)
//             return true;
//         if(b == refi)
//             return false;
//         Eigen::Vector3d va, vb;
//         va = pointset[a] - pointset[refi];
//         vb = pointset[b] - pointset[refi];
//         return reference.dot(va.cross(vb)) > 0;
//     }
// };


class CDihedralAngle //: public json::CJSONValueObject<CDihedralAngle>
{
    size_t m_Face0;
    size_t m_Face1;
    double m_Angle;
    public:
        CDihedralAngle(size_t f0, size_t f1, double theta) : /*CJSONValueObject<CDihedralAngle>("", this),*/ m_Face0(f0), m_Face1(f1), m_Angle(theta) { /*SetupJSONObject();*/ }
        CDihedralAngle(const CDihedralAngle& src) : /*CJSONValueObject<CDihedralAngle>("", this),*/ m_Face0(0), m_Face1(0), m_Angle(0.0) { CopyFrom(src); }
        ~CDihedralAngle() {}
        CDihedralAngle& CopyFrom(const CDihedralAngle& src) { m_Face0 = src.GetFace0(); m_Face1 = src.GetFace1(); m_Angle = src.GetAngle(); return *this; }
        CDihedralAngle& operator = (const CDihedralAngle& src) { return CopyFrom(src); }
        const size_t& GetFace0() const { return m_Face0; }
        const size_t& GetFace1() const { return m_Face1; }
        const double& GetAngle() const { return m_Angle; }
        void SetFace0(const size_t& x) { m_Face0 = x; }
        void SetFace1(const size_t& x) { m_Face1 = x; }
        void SetAngle(const double& x) { m_Angle = x; }
        // void SetupJSONObject()
        // {
        //     AddUIntegerValue("Face-0", &m_Face0);
        //     AddUIntegerValue("Face-1", &m_Face1);
        //     AddFloatingPointValue("Angle", &m_Angle);
        // }
};

//TODO: compute convex hull.
class Polyhedron //: public json::CJSONValueObject<CPolyhedron>
{
public:
    Polyhedron(string name="Polyhedron") : m_Name(name) //: CJSONValueObject<Polyhedron>("Polyhedron", this)
    {
        m_Graph.AddNetworkType(graph::Undirected);
        m_FaceGraph.AddNetworkType(graph::Undirected);
        // #ifdef c_plus_plus_11
        // utils::rng_util::seed_generator(rng);
        // #endif
        // SetupJSONObject();
    }

    Polyhedron( const std::vector<Eigen::Vector3d>& points,
                string name="Polyhedron") : m_Name(name)
    {
        // compute convex hull and set the other data structures.
        m_Graph.AddNetworkType(graph::Undirected);
        m_FaceGraph.AddNetworkType(graph::Undirected);
        Set(points);
    }

    // Polyhedron( const std::vector<Eigen::Vector3d>& points,
    //             std::vector< std::vector<int> > faces,
    //             string name="Polyhedron") : m_Name(name)
    // {
    //     // precompute other things from this data.
    //     Set(points, faces);
    // }

    ~Polyhedron() { /*cout << "Calling ~Polyhedron" << endl; */}

// JSON Object methods
    // void SetupJSONObject()
    // {
    //     AddNameValuePair<string, json::CJSONValueString>("Name", &m_Name);
    //     AddNameValuePair<vector< vector<int> >, json::CJSONValueArray< vector<int>,
    //                                                              json::CJSONValueArray< int, json::CJSONValueInt> > >("Faces", &m_Faces);
    //     AddNameValuePair<vector< vector<int> >, json::CJSONValueArray< vector<int>,
    //                                                              json::CJSONValueArray< int, json::CJSONValueInt> > >("Edges", &m_Edges);
    //     AddNameValuePair<vector< vector<double> >, json::CJSONValueArray< vector<double>,
    //                                                              json::CJSONValueArray< double, json::CJSONValueDouble> > >("Vertices", &m_Vertices);
    //
    //     AddNameValuePair<vector< vector< vector<double> > >, json::CJSONValueArray< vector< vector< double > >, json::CJSONValueArray< vector<double>,
    //                                                              json::CJSONValueArray< double, json::CJSONValueDouble> > > >("FaceShape", &m_FaceShape);
    //
    //     AddNameValuePair< vector<int>, json::CJSONValueArray< int, json::CJSONValueInt> >("FaceShapeIndices", &m_FaceShapeIndices);
    //
    //     AddNameValuePair<vector< vector<int> >, json::CJSONValueArray< vector<int>,
    //                                                              json::CJSONValueArray< int, json::CJSONValueInt> > >("AdjacentFaceIndices", &m_AdjFaceIndices);
    //
    //     AddNameValuePair<   vector<CDihedralAngle>,
    //                         json::CJSONValueArray<  CDihedralAngle,
    //                                                 json::CJSONValueObject<CDihedralAngle> > >("DihedralAngles", &m_DihedralAngles);
    //     m_Graph.SetupJSONObject();
    //     AddObjectValue("Graph", &m_Graph);
    // }

// Class methodss
    // #ifdef c_plus_plus_11
/*
    template<class RNG>
    void GenerateRandomWeights(graph::CNetwork& rnd, RNG& rng)
    {
//            size_t ct = 0;
        std::uniform_real_distribution<double> u01(0,1);
        auto uniform01 = std::bind(u01, rng);

        if(!rnd.IsInitialized()) rnd.Initialize(int(m_Graph.NumNodes()));
        for(graph::CNetwork::edge_iterator iter(&m_Graph); iter.IsValid(); iter++)
        {
            rnd.AddEdge(size_t(iter->source), size_t(iter->target), uniform01());
        }
    }
    // #endif

    size_t FindEdge(int v1, int v2)
    {
        for(size_t e = 0; e < m_Edges.size(); e++)
        {
            if(utils::IsInVec(m_Edges[e], v1) && utils::IsInVec(m_Edges[e], v2))
            {
                return e;
            }
        }
        return m_Edges.size();
    }
*/
    vector<int> FindFaces(int v1, int v2) const
    {
        vector<int> retval;
        for(int f = 0; size_t(f) < m_Faces.size(); f++)
        {
            if ( IsVertexInFace(f, v1) && IsVertexInFace(f, v2))
            {
                retval.push_back(f);
            }
        }

        return retval;
    }

    bool IsVertexInFace(const int& fid, const int& vertex) const
    {
        size_t fndx = size_t(fid);
        for (size_t n = 0; n < m_Faces[fndx].size(); n++) {
            if(m_Faces[fndx][n] == vertex)
                return true;
        }
        return false;
    }

    void PrintInfo() const
    {
        cout << "Polyhedron : " << m_Name << endl;
        cout << "Faces : " << endl;
        for(size_t i = 0; i < m_Faces.size(); i++)
        {
            for( size_t j = 0; j< m_Faces[i].size(); j++)
            {
                cout << m_Faces[i][j] <<" ";
            }
            cout << endl;
        }
        /*
        cout << "Edges : " << endl;
        for(size_t i = 0; i < m_Edges.size(); i++)
        {
            for( size_t j = 0; j< m_Edges[i].size(); j++)
            {
                cout << m_Edges[i][j] <<" ";
            }
            cout << endl;
        }
        */
        cout << "Vertices : " << endl;
        for(size_t i = 0; i < m_Vertices.size(); i++)
        {
            for( size_t j = 0; j< m_Vertices[i].size(); j++)
            {
                cout << m_Vertices[i][j] <<" ";
            }
            cout << endl;
        }
    }

// Accessors
    const string& GetName() const {return m_Name; }

    void SetName(const string& name) { m_Name = name; }

    size_t NumVertices() const { return m_Vertices.size(); }
/*
    size_t NumEdges() const { return m_Edges.size(); }
*/
    size_t NumFaces() const { return m_Faces.size(); }
/*
    const vector< vector<int> >& GetEdges() { return m_Edges; }

    int GetFaceShapeIndex(size_t face) const { return m_FaceShapeIndices[face]; }

    const vector< vector<double> >& GetFaceShape(size_t face) const { return m_FaceShape[size_t(m_FaceShapeIndices[face])]; } // copy into a matrix?
*/
    const vector< vector<int> >& GetFaces() const { return m_Faces; }

    const vector< Eigen::Vector3d >& GetVertices() const { return m_Vertices; }

    PointMatrix Vertices(bool shift=false, bool scale=false) const
    {
        double s = scale ? pow((1.0/m_volume), (1.0/3.0)) : 1.0;
        PointMatrix verts(m_Vertices.size(), 3);
        for(int i = 0; i < m_Vertices.size(); i++)
            verts.row(i) = shift ? m_Vertices[i] - m_Centroid : m_Vertices[i];
        return s*verts;
    }

    void FacetNormals(PointMatrix& normals, bool bMerged=true) const
    {
        if(bMerged && !m_MergedFaces.size())
            throw std::runtime_error("must merge facets before finding normals");
        // std::cout << "merged faces size = " << m_MergedFaces.size() << std::endl;
        const std::vector< std::vector<int> >& faces = (bMerged) ? m_MergedFaces : m_Faces;
        normals.resize(faces.size(), 3);
        for(int f = 0; f < faces.size(); f++)
        {
            normals.row(f) = chull::getOutwardNormal(m_Vertices, m_Centroid, faces[f]);
        }
        // std::cout << "normals: \n" << normals << std::endl;
        // return normals;
    }

    void MergeFacets(double threshold=1e-6, bool bForce=true)
    {
        if(m_MergedFaces.size() && !bForce)
            return;
        m_MergedFaces.clear(); // clear the data.
        m_MergedFaces.reserve(m_Faces.size());
        double thresh = 1.0 - threshold;
        PointMatrix normals;
        FacetNormals(normals, false);
        PointMatrix sim = normals * normals.transpose();
        m_MergeMap.resize(m_Faces.size(), -1);
        std::vector<bool> found(m_Faces.size(), false);
        std::queue<int> worklist;
        chull::GrahamScan<double> ghull(m_Vertices);
        for(size_t f = 0; f < m_Faces.size(); f++)
        {
            if(found[f])
                continue;
            found[f] = true;
            std::set<int> face;
            graph::CNetwork::vertex_iterator bfs(&m_FaceGraph, f, graph::CNetwork::vertex_iterator::BreadthFirstSearch, [sim, thresh](int i, int j){ return sim(i,j) > thresh; });
            for (/* bfs intialized*/; bfs.IsValid(); ++bfs)
            {
                // std::cout << "face " << f  << " is the same as " << *bfs << std::endl;
                face.insert(m_Faces[*bfs].begin(), m_Faces[*bfs].end());
                found[*bfs] = true;
                m_MergeMap[*bfs] = m_MergedFaces.size();
            }
            if(face.size() > 3)
            {
                std::vector<int> ret;
                ret.reserve(face.size());
                std::back_insert_iterator< std::vector<int> > it(ret);
                ghull.compute(it, face.begin(), face.end(), normals.row(f));
                m_MergedFaces.push_back(ret);
            }
            else if(face.size() == 3)
            {
                std::vector<int> ret(face.begin(), face.end());
                m_MergedFaces.push_back(ret);
            }
            else{
                throw std::runtime_error("face with fewer than 3 vertices found");
            }
        }
        // int i = 0;
        // for(auto face : m_MergedFaces)
        // {
        //     std::cout << "face " << i++ << ": " << std::endl;
        //     for(auto f : face)
        //     {
        //         std:: cout << f << ", ";
        //     }
        //     std::cout << std::endl;
        // }
        // std::cout << "done!" << std::endl;
    }

    const double& Volume() { return m_volume; }

    const Eigen::Vector3d& Centroid() { return m_Centroid; }

    void GetMergedFaceGraph(graph::CNetwork& g) const
    {
        std::map<size_t, size_t> _;
        std::map<size_t, size_t> __;
        g = m_FaceGraph.QuotientGraph(m_MergeMap, _ , __ ); }

    void Set(Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> verts)
    {
        vector< Eigen::Vector3d > vs(verts.rows());
        for(int i = 0; i < verts.rows(); i++)
            vs[i] = verts.row(i);
        Set(vs);
    }

    void Set(const vector< Eigen::Vector3d >& verts)
    {
        chull::ConvexHull<double> hull(verts);
        hull.compute();
        std::vector< int > vmap(verts.size(), -1); // -1 mean not mapped
        m_Faces.reserve(hull.getFaces().size());
        vector< vector<unsigned int> >::const_iterator fbegin = hull.getFaces().cbegin();
        vector< vector<unsigned int> >::const_iterator fend = hull.getFaces().cend();
        int ct = 0;
        for(vector< vector<unsigned int> >::const_iterator fi = fbegin; fi != fend; fi++)
        {
            std::vector<int> face;
            face.reserve(fi->size());
            for(unsigned int vi : *fi) // for each vertex in faces
            {
                if(vmap[vi] < 0)
                {
                    vmap[vi] = ct++;
                }
                face.push_back(vmap[vi]);
                // std::cout << "("<< vi << ", " <<  vmap[vi] << "), ";
            }
            // std::cout << std::endl;
            m_Faces.push_back(face);
        }
        m_Vertices.resize(ct);
        for(int vi = 0; vi < vmap.size(); vi++)
        {
            if(vmap[vi] >= 0)
            {
                m_Vertices[vmap[vi]] = verts[vi];
            }
        }
        sortFaces(); // assumes 3d, triangular faces!

        computeVertFaces();
        computeMassProps();
    }

    void Set(const vector< Eigen::Vector3d >& verts, vector< vector<int> >& faces)
    {
        m_Vertices.clear(); m_Vertices.reserve(verts.size());
        std::map<int, int> vmap;
        int ct = 0;
        for(int i = 0; i < verts.size(); i++)
        {
            if(vmap.find(i) == vmap.end())
            {
                vmap[i] = ct++;
                m_Vertices.push_back(verts[i]);
                for(int j = i+1; j < verts.size(); j++)
                {
                    if((verts[i]-verts[j]).norm() < double(1e-4))
                    {
                        vmap[j] = vmap[i];
                    }
                }
            }
        }
        Eigen::Vector3d vavg = Eigen::Vector3d::Zero();
        for(auto v : m_Vertices)
            vavg += v;
        vavg /= double(m_Vertices.size());
        m_Faces.clear();
        for(size_t i = 0; i < faces.size(); i++)
        {
            std::vector<int> face;
            for(size_t _i = 0; _i < faces[i].size(); _i++) utils::PushUnique(face, vmap[faces[i][_i]]);
            chull::counter_clockwise<double> ccw(m_Vertices[vmap[faces[i][0]]], chull::getOutwardNormal(m_Vertices, vavg, faces[i]));
            auto lessfn = ccw.lessfn<int, decltype(m_Vertices) >(m_Vertices);
            std::sort(face.begin(), face.end(), lessfn);
            std::cout << "sorted - "<< i <<": ";
            for(size_t _i = 0; _i < face.size(); _i++)  std::cout << face[_i] << ", ";
            std::cout << std::endl;

            for(size_t _i = 0; _i < face.size(); _i++)
            {
                std::cout   << "vertex " << face[_i] << std::endl
                            << m_Vertices[face[_i]] << std::endl;
            }

            if(face.size() > 3)
            {
                for(size_t v = 1; v < face.size()-1; v++)
                {
                    // std::cout << "making face: {" << face[0] << ", " << face[v] << ", " <<  face[v+1] << "}"<< std::endl;
                    // std::cout << "making face: {" << m_Vertices[face[0]] << ", " << m_Vertices[face[v]] << ", " <<  m_Vertices[face[v+1]] << "}"<< std::endl;
                    m_Faces.push_back(std::vector<int>({face[0], face[v], face[v+1]}));
                }
            }
            else if(face.size() == 3)
            {
                m_Faces.push_back( face );
            }
            else
            {
                std::cout << "Faces must have at least three vertices." << std::endl;
                std::cout << "face - "<< i <<" ( size "<< faces[i].size() <<"): ";
                for(size_t _i = 0; _i < face.size(); _i++)  std::cout << face[_i] << ", ";
                std::cout << std::endl;
                for(size_t _i = 0; _i < faces[i].size(); _i++)  std::cout << faces[i][_i] << ", ";
                std::cout << std::endl;

                throw(std::runtime_error("Error: Can not set Polyhedron!"));
            }
        }
        sortFaces();
        computeVertFaces();
        computeMassProps();
    }

private:
    void computeMassProps()
    {
        m_Graph.Initialize(m_Vertices.size());
        m_FaceGraph.Initialize(m_Faces.size());
        // requires that the faces are sorted and triangulated.
        mass_properties<Polyhedron> mp(*this);
        m_volume = mp.getVolume();
        m_Centroid = mp.getCenterOfMass();

        // also set the dihedrals and adjacency information.
        m_DihedralAngles.clear();
        for(int i = 0; i < m_Faces.size(); i++)
        {
            int n = m_Faces[i].size();
            for(int v = 0; v < n; v++)
            {
                std::vector<int> faces = FindFaces(m_Faces[i][v], m_Faces[i][(v+1)%n]);
                assert(faces.size() == 2);
                if(faces[0] < i) continue;
                m_Graph.AddEdge(m_Faces[i][v], m_Faces[i][(v+1)%n]);
                Eigen::Vector3d n0 = chull::getOutwardNormal(m_Vertices, m_Centroid, m_Faces[faces[0]]);
                Eigen::Vector3d n1 = -chull::getOutwardNormal(m_Vertices, m_Centroid, m_Faces[faces[1]]);
                n0.normalize();
                n1.normalize();
                double cosd = n0.dot(n1); // cos(-pi) to cos(pi).
                // std::cout << "(" << faces[0] << ", " << faces[1] << "): " << acos(cosd) << std::endl;
                m_DihedralAngles.push_back(CDihedralAngle(faces[0], faces[1], acos(cosd)));
                m_FaceGraph.AddEdge(faces[0], faces[1]);
            }
        }
    }

    void computeVertFaces()
    {
        // Now I need to sort the faces connected to each vertex in counter clockwise order
        m_VertFaces.resize(m_Vertices.size());
        std::vector<Eigen::Vector3d> fcenters;
        for(int f = 0; f < m_Faces.size(); f++)
        {
            Eigen::Vector3d center = Eigen::Vector3d::Zero();
            for(int f_ = 0; f_ < m_Faces[f].size(); f_++)
            {
                center += m_Vertices[m_Faces[f][f_]];
            }
            center /= double(m_Faces[f].size()); // correct because triangles
            fcenters.push_back(center);
        }
        for( int v = 0; v < m_Vertices.size(); v++)
        {
            std::vector<int> vf;
            for(int f = 0; f < m_Faces.size(); f++)
            {
                if(std::find(m_Faces[f].begin(), m_Faces[f].end(), v) != m_Faces[f].end())
                {
                    vf.push_back(f);
                }
            }
            if(vf.size() <= 2)
            {
                std::cout   << "v = " << v << std::endl
                            << "vfs = " << vf.size() << std::endl
                            << "vf: ";
                for(auto ii : vf ) std::cout << ii << ", ";
                std::cout   << std::endl
                            << "faces:" << std::endl;
                for(auto& f : m_Faces ) {
                    for( auto fi : f)
                        std::cout << fi << ", ";
                    std::cout << std::endl;
                }
            }
            assert(vf.size() > 2);
            auto n = m_Vertices[v];
            n.normalize();
            chull::PlaneRn<double, 3> support(n, m_Vertices[v]);
            chull::counter_clockwise<double> ccw(fcenters[0], n);
            for(size_t i = 0; i < fcenters.size(); i++)
                fcenters[i] = support.projection(fcenters[i]);
            auto lessfn = ccw.lessfn<int, decltype(fcenters) >(fcenters);
            std::sort(vf.begin(), vf.end(), lessfn);
            m_VertFaces[v]=vf;
        }
    }

    void sortFace(const std::vector< Eigen::Vector3d >& points, const Eigen::Vector3d& inside_point, std::vector< std::vector<int> >& faces, const unsigned int& faceid, double thresh = 0.0001)
    {
        Eigen::Vector3d a = points[faces[faceid][0]], b = points[faces[faceid][1]], c = points[faces[faceid][2]], n, nout;
        nout = chull::getOutwardNormal(points, inside_point, faces[faceid]);
        n = (b - a).cross(c - a);
        if ( nout.dot(n) < 0 )
            std::reverse(faces[faceid].begin(), faces[faceid].end());
    }

    void sortFaces(double thresh = 1e-6)
    {
        Eigen::Vector3d inside_point(0.0,0.0,0.0);
        for(size_t i = 0; i < m_Vertices.size(); i++)
        {
            inside_point += m_Vertices[i];
        }
        inside_point /= double(m_Vertices.size());

        for( unsigned int f = 0; f < m_Faces.size(); f++ )
            sortFace(m_Vertices, inside_point, m_Faces, f, thresh);
    }

/*
        bool AreAdjacentFaces(const size_t& f, const size_t& g)
        {
            bool bFound = false;
            for(size_t i = 0; i < m_AdjFaceIndices.size() && !bFound; i++)
            {
                if(utils::IsInVec(m_AdjFaceIndices[i], int(f)) && utils::IsInVec(m_AdjFaceIndices[i], int(g)) )
                {
                    bFound = true;
                }
            }
            return bFound;
        }

        int FindConnectingEdge( const int& f1, const int& f2)
        {
            for(size_t e = 0; e < NumEdges(); e++)
            {
                vector<int> faces = FindFaces(m_Edges[e][0], m_Edges[e][1]);
                if(utils::IsInVec(faces, f1) && utils::IsInVec(faces, f2))
                    return int(e);
            }
            return int(NumEdges());
        }

        const vector<CDihedralAngle>& GetDihedrals() const { return m_DihedralAngles; }
*/

public:

    const vector< vector<int> >& AdjFaceIndices() { return m_AdjFaceIndices; }

    const vector<CDihedralAngle>& GetDihedrals() const { return m_DihedralAngles; }

    void Dual(Polyhedron& dual, const Eigen::Vector3d& origin, double radius) const
    {
        // Assume the origin is inside the polyhedron
        // computes the dual polyhedron with respect to the sphere defined by
        // origin and radius.
        // the distances from the origin for the polyhedron (d) and its dual (d') satisfy:
        // d * d' = r^2
        // v = (x1, x2, x3)
        // f = ((n1, n2, n3), (p1,p2,p3)) (unit normals and closest point to the origin)
        // v' = (n1,n2,n3)* r / d((p1,p2,p3))
        // f' = ((x1,x2,x3)/d(v), r/d(v)*((x1,x2,x3)/d(v)))
        // e' = m_AdjFaceIndices
        // m_AdjFaceIndices' = m_Edges
        // faces' = // list of faces that contain the v, then sort faces.

        double r2 = radius*radius;
        std::vector< Eigen::Vector3d > verts(m_Faces.size());
        for(size_t f = 0; f < m_Faces.size(); f++)
        {
            // std::cout << "facet["<< f <<"]: " << m_Faces[f].size() << std::endl;
            assert(m_Faces[f].size() >= 3);
            // std::cout << "facet["<< f <<"]: " << m_Faces[f][0] << ", " << m_Faces[f][1] << ","<< m_Faces[f][2] << std::endl;
            // std::cout << "nVerts: " << m_Vertices.size() << std::endl;
            chull::PlaneRn<double, 3> facet = chull::PlaneRn<double, 3>::Through(m_Vertices[m_Faces[f][0]]-origin, m_Vertices[m_Faces[f][1]]-origin, m_Vertices[m_Faces[f][2]]-origin);
            double dist = -facet.signedDistance(origin);
            // std::cout << "offset: " << facet.offset() << ",  dist: " << dist << std::endl;
            // std::cout << "normal: " << facet.normal() << std::endl;
            if(fabs(facet.offset() - facet.normal().dot(origin)) < 1e-6)
                throw std::runtime_error("Facet runs throught the origin of sphere -- not well posed to take dual");
                // std::cout << m_Vertices[m_Faces[f][0]]  << " \n\n" << m_Vertices[m_Faces[f][1]] << " \n\n" << m_Vertices[m_Faces[f][2]] << " \n\n" << std::endl;

            verts[f] = (facet.normal()*(r2/dist)) + origin;
            // std::cout << "Vert norm: " << verts[f].norm() << std::endl;
        }

        // std::vector< std::vector<int> > temp(m_VertFaces); //TODO: remove this expensive copy
        dual.Set(verts);
        dual.SetName(m_Name + "-Dual");

    }

    void writePOS(const string& filename, std::ios_base::openmode mode=std::ios_base::out) const
    {
        std::ofstream file(filename, mode);
        std::stringstream ss, connections;
        ss << "def "<< m_Name << " \"poly3d " << m_Vertices.size() << " ";
        std::string center_sphere  = "def origin \"sphere 0.1 005F5F5F\"";
        for(vector< Eigen::Vector3d >::const_iterator iter = m_Vertices.cbegin(); iter != m_Vertices.cend(); iter++)
            ss << (*iter)[0] << " " << (*iter)[1] << " " << (*iter)[2] << " ";
        ss << "005984FF\"";
        std::string hull  = ss.str();
        file << hull << std::endl;
        file << center_sphere << std::endl;
        file << m_Name << " 0 0 0 1 0 0 0" << std::endl;
        file << "origin 0 0 0 "<< std::endl;
        file << "eof" << std::endl;
        file.close();
    }

private:
    graph::CNetwork             m_Graph;
    graph::CNetwork             m_FaceGraph;
    string                      m_Name;

    vector< vector<int> >       m_Faces;        // these faces are always triangular
    vector< vector<int> >       m_MergedFaces;  // these faces need not bee triangular
    vector< vector<int> >       m_VertFaces;    // these are faces adjacent to vertex i
    vector<size_t>       m_MergeMap;    // id of the merged face.
    vector< Eigen::Vector3d >   m_Vertices;     // the vertex representaion of the polyhedron
    double                      m_volume;       // the volume of the polyhedron
    Eigen::Vector3d             m_Centroid;     // the Centroid of the polyhedron

    // vector< vector<int> >       m_Edges;
    vector< vector<int> >       m_AdjFaceIndices;
    // vector<int>                 m_FaceShapeIndices;
    vector<CDihedralAngle>      m_DihedralAngles;
    // vector< vector< vector<double> > >     m_FaceShape;
};

// from above...
mass_properties< Polyhedron >::mass_properties(const Polyhedron& poly) :  points(poly.GetVertices()), faces(poly.GetFaces())
{
    compute();
}

inline void intersection(const Polyhedron& A, const Polyhedron& B, const Eigen::Vector3d& inside, Polyhedron& ret)
{
    Polyhedron dA, dB, temp;
    A.Dual(dA, inside, 1.0);
    B.Dual(dB, inside, 1.0);
    int Na, Nb;
    Na = dA.GetVertices().size();
    Nb = dB.GetVertices().size();

    vector< Eigen::Vector3d > verts(Na + Nb);
    // std::cout << "A* verts("<< Na << "): "<< std::endl;
    for(int i = 0; i < Na; i++){
        // std::cout << dA.GetVertices()[i] << "\n" << std::endl;
        verts[i] = dA.GetVertices()[i];
    }
    // std::cout << "B* verts("<< Nb << "): "<< std::endl;
    for(int i = 0; i < Nb; i++){
        // std::cout << dB.GetVertices()[i] << "\n" << std::endl;
        verts[i+Na] = dB.GetVertices()[i];
    }
    try{
        temp.Set(verts);
    } catch(std::runtime_error error)
    {
        std::cerr << "could not compute ConvexHull for dual union" << std::endl;
        throw(error);
    }

    try{
        temp.Dual(ret, inside, 1.0);
    } catch(std::runtime_error error)
    {
        std::cerr << "could not compute dual for dual union" << std::endl;
        throw(error);
    }
}


}

}



#endif
