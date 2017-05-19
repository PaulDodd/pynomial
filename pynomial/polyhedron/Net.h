//
//  Net.h
//  folding
//
//  Created by Paul M Dodd on 5/9/14.
//
//
#if 0
#ifndef folding_Net_h
#define folding_Net_h

#include "SharedInclude.h"
#include "Network.h"
#include "Polyhedra.h"
// #include "Shape.h"



#define ALPHABET "01ABCDEFGHIJKLMNOPQRSTUVWXYZ"
// JSON Parser String constants for CNetInfo class.
#define JSON_NET_INFO_NUM_VERTEX_CONNECTIONS        "NumVertexConnections"
#define JSON_NET_INFO_VERTEX_CONNECTIONS_LIST       "ListVertexConnections"


// JSON Parser String constants for CNet class.
#define JSON_NET_ID                     "NetId"
#define JSON_NET_POLYHEDRON             "Polyhedron"
#define JSON_NET_VERTEX_GRAPH           "VertexGraph"
#define JSON_NET_FACE_GRAPH             "FaceGraph"
#define JSON_NET_CUT_GRAPH              "CutGraph"
#define JSON_NET_FACES                  "Faces"
#define JSON_NET_EDGES                  "Edges"
#define JSON_NET_VERTICES               "Vertices"
#define JSON_NET_NETINFO                "NetInfo"
#define JSON_NET_GLUING                 "Gluing"
#define JSON_NET_HINGE                  "Hinge"
#define JSON_NET_ISOMORPHISMS           "Isomorphisms"
#define JSON_NET_CMP_STRING             "CmpString"
#define JSON_NET_DIHEDRALS              "DihedralAngles"

namespace NetAnalysis{

using namespace std; // TODO: remove this
using json::CJSONValueObject;

inline size_t invalid_index() { return size_t(-1); }

class CNetInfo : public json::CJSONValueObject<CNetInfo>
{

    public:

    // Constructors
        CNetInfo() : CJSONValueObject<CNetInfo>("", this), m_NumVertexConnections(0) { SetupJSONObject(); }

    // Destructors
        ~CNetInfo() {}

    // JSON File I/O
        void SetupJSONObject()
        {
            AddUIntegerValue(JSON_NET_INFO_NUM_VERTEX_CONNECTIONS, &m_NumVertexConnections);
            AddNameValuePair<   vector<size_t>,
                                json::CJSONValueArray<  size_t,
                                                        json::CJSONValueUInt> >(JSON_NET_INFO_VERTEX_CONNECTIONS_LIST, &m_VertexConnections);
        }

    // Accessor Methods
        void            SetVertexConnections() { m_NumVertexConnections = m_VertexConnections.size(); }
        size_t          GetNumVertexConnections() const { return m_NumVertexConnections; }
        vector<size_t>  GetVertexConnections() const { return m_VertexConnections; }

    // overloaded operators
        CNetInfo& operator = (const CNetInfo& src) { return CopyFrom(src);}
    // Helper functions
        CNetInfo& CopyFrom(const CNetInfo& src)
        {
            m_NumVertexConnections = src.GetNumVertexConnections();
            m_VertexConnections = src.GetVertexConnections();
            return *this;
        }

    private:
        size_t          m_NumVertexConnections;
        vector<size_t>  m_VertexConnections;
};

class CNet : public json::CJSONValueObject<CNet>
{
#ifdef c_plus_plus_11
    typedef typename utils::float_vec_is_equal<double, -4> float_vec_compare_func;
#endif
    public:
    // Constrctors
        CNet(bool ignore = false) : CJSONValueObject<CNet>("", this), m_bInitialized(false), m_bIgnoreSymmetry(ignore), m_Id(size_t(-1))/* m_pShape(NULL),*/  { SetupJSONObject(); }

        CNet(const CNet& src) : CJSONValueObject<CNet>("", this), m_bInitialized(false), m_bIgnoreSymmetry(src.GetIgnoreSymmetryFlag()), m_Id(size_t(-1))  /* m_pShape(NULL),*/
        {
            CopyFrom(src);
            SetupJSONObject();
        }

    // Destrctors
        ~CNet() { Destroy(); }

        void Destroy()
        {
//            cout << "Destroying net " << m_Id << "!"<< endl;
//            m_VertexGraph.Destroy();
//            m_FaceGraph.Destroy();
//            m_CuttingGraph.Destroy();
//            cout << "Graphs destroyed " << endl;
        }

    // JSON File I/O
        void SetupJSONObject()
        {
            AddUIntegerValue(JSON_NET_ID, &m_Id);

            AddStringValue(JSON_NET_POLYHEDRON, &m_ParentPolyhedron);

            AddStringValue(JSON_NET_CMP_STRING, &m_PolygonCompareString);

            // m_VertexGraph.SetupJSONObject();
            AddObjectValue(JSON_NET_VERTEX_GRAPH, &m_VertexGraph);

            // m_FaceGraph.SetupJSONObject();
            AddObjectValue(JSON_NET_FACE_GRAPH, &m_FaceGraph);

            AddObjectValue(JSON_NET_CUT_GRAPH, &m_CuttingGraph);

            AddNameValuePair<   vector< vector<size_t> >,
                                json::CJSONValueArray<  vector< size_t >,
                                                        json::CJSONValueArray<  size_t,
                                                                                json::CJSONValueUInt > > >(  JSON_NET_FACES, &m_Faces);
            AddNameValuePair<   vector<vector<size_t> >,
                                json::CJSONValueArray<  vector< size_t >,
                                                        json::CJSONValueArray<  size_t,
                                                                                json::CJSONValueUInt > > >(  JSON_NET_EDGES, &m_Edges);

            AddNameValuePair<   vector<vector<double> >,
                                json::CJSONValueArray<  vector< double >,
                                                        json::CJSONValueArray<  double,
                                                                                json::CJSONValueDouble > > >( JSON_NET_VERTICES, &m_Vertices);

            AddNameValuePair<   vector< bool >, json::CJSONValueArray< bool, json::CJSONValueBool > >( JSON_NET_HINGE, &m_bHinge);

            AddNameValuePair<   vector< vector<size_t> >,
                                json::CJSONValueArray<  vector< size_t >,
                                                        json::CJSONValueArray<  size_t,
                                                                                json::CJSONValueUInt > > >(  JSON_NET_ISOMORPHISMS, &m_Isomorphisms);
            AddNameValuePair<   vector< vector<size_t> >,
                                json::CJSONValueArray<  vector< size_t >,
                                                        json::CJSONValueArray<  size_t,
                                                                                json::CJSONValueUInt > > >(  JSON_NET_GLUING, &m_Gluing);

            AddNameValuePair<   vector<unfold_polyhedra::CDihedralAngle>,
                                json::CJSONValueArray<  unfold_polyhedra::CDihedralAngle,
                                                        json::CJSONValueObject<unfold_polyhedra::CDihedralAngle> > >(JSON_NET_DIHEDRALS, &m_DihedralAngles);
            m_NetInfo.SetupJSONObject();
            AddObjectValue(JSON_NET_NETINFO, &m_NetInfo);
        }

        bool LoadFromFile( const string& Path )
        {
            m_bInitialized = false;
            if(CJSONValueObject<CNet>::LoadFromFile(Path))
            {
            // put any other initialization methods here
                ComputeFreeEdgeLabels();
                m_bInitialized = true;
            }
            return m_bInitialized;
        }

        bool SaveToFile( const string& Path )
        {
            return CJSONValueObject<CNet>::SaveToFile(Path); // call trhough to the default.
        }
#ifdef c_plus_plus_11
        bool LoadFromTXTFile(const string& netpath)
        {
            cout << "Warning! This is depricated. Should use json file format if possible." << endl;
            ifstream netfile;
            printf("Net path: %s\n", netpath.c_str());
            netfile.open(netpath.c_str(), ios_base::in);
            bool bParseSuccess = false;
            if(netfile.is_open())
            {
                char line[3200];
                bool bContinue = true;
                vector < vector< pair<double, double> > > coords;

                while(bContinue)
                {
                    memset(line, 0x00, 3200*sizeof(char));
                    netfile.getline(line, 3200, '\n');
                    string  sline = line;
                    cout << sline<< endl;
                    if(sline == "VERTEX_COORDS")
                    {
                        continue;
                    }
                    else if(sline == "BUILD_ORDER")
                    {
                        bContinue = false;
                    }
                    else
                    {
                        vector<string> svec = utils::SplitString(sline, ";");
                        vector< pair<double, double> > face_coords;
                        vector<size_t> face;
                        for(size_t i = 0; i < svec.size(); i++)
                        {
                            vector<string> stemp = utils::SplitString(svec[i].substr(svec[i].find_first_of('{')+1, svec[i].find_first_of('}') - svec[i].find_first_of('{')-1), ",");
                            cout << stemp[0] << " "<< stemp[1]<< endl;
                            pair <double, double> point;
                            vector<double> vpoint;

                            point.first = atof(stemp[0].c_str());
                            point.second = atof(stemp[1].c_str());
                            face_coords.push_back(point);
                            vpoint.push_back(point.first); vpoint.push_back(point.second);
                            utils::PushUnique<vector<double>, float_vec_compare_func > (m_Vertices, vpoint);

                            size_t vin = utils::FindInVec<vector<double>, float_vec_compare_func>(m_Vertices, vpoint);
                            if(vin >= m_Vertices.size())
                            {
                                cout << "Could not retrieve the vertex."<< vin << " is out of range" << m_Vertices.size() << endl;
                            }
                            face.push_back(vin);
                            // cout << "Number of vertices = " << m_Vertices.size() << endl;
                        }
                        //cout << m_Faces.size() << ": ";
                        for(size_t f = 0; f < face.size(); f++ )
                        {
                            cout << face[f] << " ";
                            vector<size_t> edge;
                            if(face[f] < face[(f+1)%face.size()])
                            {
                                edge.push_back(face[f]); edge.push_back(face[(f+1)%face.size()]);
                            }
                            else
                            {
                                edge.push_back(face[(f+1)%face.size()]); edge.push_back(face[f]);
                            }
                            utils::PushUnique(m_Edges, edge);
                        }
                        cout << endl;
                        coords.push_back(face_coords);
                        m_Faces.push_back(face);
                    }
                }

                ComputeGraphData();
                m_bInitialized = true;

                PrintNet();
                bParseSuccess = true;
            }
            return bParseSuccess;
        }
#endif

    // Class Methods
        size_t FindEdge(const size_t& v1, const size_t& v2) const
        {
            for(size_t e = 0; e < m_Edges.size(); e++)
            {
                if(utils::IsInVec(m_Edges[e], v1) && utils::IsInVec(m_Edges[e], v2))
                {
                    return e;
                }
            }
            return invalid_index();
        }


        template<class TVal>
        size_t FindVertexByCoord(const TVal& point) const
        {
            for (size_t v = 0; v < m_Vertices.size(); v++)
            {
                if( point[0] == m_Vertices[v][0] &&
                    point[1] == m_Vertices[v][1])
                {
                    return v;
                }
            }
            return invalid_index();
        }

        vector<size_t> FindEdgesWithVertex(const size_t& vertex) const
        {
            vector<size_t> edges;
            for(size_t e = 0; e < m_Edges.size(); e++)
            {
                if( m_Edges[e][0] == vertex ||
                    m_Edges[e][1] == vertex )
                {
                    edges.push_back(e);
                }
            }
            return edges;
        }

        vector<size_t> FindFacesWithVertex(const size_t& vertex) const
        {
            assert(vertex < m_Vertices.size());

            vector<size_t> faces;
            for( size_t f = 0; f < m_Faces.size(); f++)
            {
                for( size_t i = 0; i < m_Faces[f].size(); i++)
                {
                    if( m_Faces[f][i] == vertex)
                        faces.push_back(f);
                }
            }
            return faces;
        }

        vector<size_t> FindCommonEdgesInFaces(size_t fc1, size_t fc2) const
        {
            vector<size_t> edges;
            for(size_t f = 0; f < m_Faces[fc1].size(); f++)
            {
                size_t ff = (f+1) % m_Faces[fc1].size();
                size_t edge1 = FindEdge(m_Faces[fc1][f], m_Faces[fc1][ff]);
                if( edge1 == invalid_index() )
                {
                    PrintNet();
                    std::stringstream ss;
                    ss << "Error!! Could not find edge in face "<< fc1 <<" edge("<<f << ", "<< ff <<"): (" << m_Faces[fc1][f] << ", " << m_Faces[fc1][ff] << ")" << endl;
                    throw std::runtime_error(ss.str());
                }
                for(size_t g = 0; g < m_Faces[fc2].size(); g++)
                {
                    size_t gg = (g+1) % m_Faces[fc2].size();
                    size_t edge2 = FindEdge(m_Faces[fc2][g], m_Faces[fc2][gg]);
                    if( edge2 == invalid_index() )
                    {
                        PrintNet();
                        std::stringstream ss;
                        ss << "Error!! Could not find edge (" << m_Faces[fc2][g] << ", " << m_Faces[fc2][gg] << ")"<< endl;
                        throw std::runtime_error(ss.str());
                    }

                    if(edge1 == edge2)
                    {
                        edges.push_back(edge1);
                        break;
                    }
                }
            }
            return edges;
        }

        vector<size_t> FindFacesWithEdge(const size_t& edge) const
        {
            vector<size_t> faces1 = FindFacesWithVertex(m_Edges[edge][0]);
            vector<size_t> faces2 = FindFacesWithVertex(m_Edges[edge][1]);
            vector<size_t> faces, intersect = utils::math_util::Intersection(faces1, faces2);
            // need to make sure that the faces are truly have that edge.
            for(size_t i = 0; i < intersect.size(); i++){
                size_t f = intersect[i];
                vector<size_t> ndx = GetVertexIndexInFaceFromEdge(f, edge);
                if(utils::mod_dist(int(ndx[0]), int(ndx[1]), int(m_Faces[f].size())) == 1)
                    utils::PushUnique(faces, f);
            }

            return faces;
        }

        size_t GetVertexIndexInFace(const size_t& face, const size_t& vert) const
        {
            for(size_t f = 0; f < m_Faces[face].size(); f++)
            {
                if( m_Faces[face][f] == vert )
                {
                    return f;
                }
            }
            cout << "Vertex not found in faces" << endl;
            return invalid_index();
        }

        vector<size_t> GetVertexIndexInFaceFromEdge(const size_t& face, const size_t& edge) const
        {
            vector<size_t> ret;
            ret.push_back(GetVertexIndexInFace(face, m_Edges[edge][0]));
            ret.push_back(GetVertexIndexInFace(face, m_Edges[edge][1]));
            return ret;
        }

        void ComputeGraphData()
        {
//            cout << "Computing Graph Data" << endl;
            m_VertexGraph.Initialize(int(m_Vertices.size()));
            m_VertexGraph.AddNetworkType((Network::Undirected | Network::Simple));

            m_FaceGraph.Initialize(int(m_Faces.size()));
            m_FaceGraph.AddNetworkType((Network::Undirected | Network::Simple));
            m_bHinge.clear();
            for(size_t e = 0; e < m_Edges.size(); e++)
            {
                vector<size_t> faces = FindFacesWithEdge(e);
                if( faces.size() == 1)
                {
                    m_bHinge.push_back(false);
                }
                else if ( faces.size() == 2)
                {
                    m_bHinge.push_back(true);
//                    cout << "Adding edge (" <<faces[0] << ", "<< faces[1]<<") "  << endl;
                    m_FaceGraph.AddEdge(faces[0], faces[1]);
                }
                else if ( faces.size() > 2 )
                {
//                    cout << "Warning! more than 2 faces attached to an edge" << endl;
                    for(size_t i = 0; i < faces.size(); i++)
                    {
                        for( size_t j = i+1; j < faces.size(); j++)
                        {
                            m_FaceGraph.AddEdge(faces[i], faces[j]);
                        }
                    }
                }
                else
                {
                    cout << "There was an error finding the edge in the face. (no face contains this edge) edge[" << e <<"] = {"<< m_Edges[e][0]<<", " << m_Edges[e][1] << "}"<< endl;
                    // throw(std::runtime_error("Edge does not exist"));
                }
//
//                for(size_t f = 0; f < faces.size(); f++)
//                {
//                    cout << faces[f] << " ";
//                }
//                cout << " " << m_bHinge[e] << endl;

                m_VertexGraph.AddEdge(m_Edges[e][0], m_Edges[e][1]);
            }
        }

        void ComputeOuterVertexCycle()
        {
            bool bContinue = true;
            bool bFoundNextEdge = false;
            size_t vc, vn, vm1;
            vc = 0;
            vm1 = 0;
            vn = 0;
            m_OuterVertexCycle.clear();
            m_OuterVertexCycle.push_back(vc);
            while ( bContinue )
            {
                vector<size_t> edges = FindEdgesWithVertex(vc);
                bFoundNextEdge = false;
                for(size_t ei = 0; ei < edges.size() && !bFoundNextEdge; ei++)
                {
                    if( !m_bHinge[edges[ei]] )
                    {

                        vn = (m_Edges[edges[ei]][0] == vc) ? m_Edges[edges[ei]][1] : m_Edges[edges[ei]][0];
                        // cout << "edge " << edges[ei] << " = " << m_Edges[edges[ei]][0] << ", " << m_Edges[edges[ei]][1] << "-> vc = " << vc << ", vn = " << vn << endl;
                        bFoundNextEdge = (vn != vm1);
                    }
                }
                if(bFoundNextEdge)
                {
                    vm1 = vc;
                    vc = vn;
                    bContinue = (vc != 0);
                    if ( bContinue )
                        m_OuterVertexCycle.push_back(vn);

                }
                else
                {
                    cout << "Error could not find next edge!" << endl;
                    break;
                }
            }

//            cout<< "Vertex Cycle: ";
//            for(size_t i = 0; i < m_OuterVertexCycle.size(); i++)
//            {
//                cout<<m_OuterVertexCycle[i] << (i == m_OuterVertexCycle.size()-1 ? "\n" : ", ");
//            }
        }

        void ComputeCompareString(bool bForce = false)
        {
            if(m_PolygonCompareString.length() != 0 && !bForce)
                return;
            m_PolygonCompareString.clear();

            if ( m_OuterVertexCycle.size() == 0)
                ComputeOuterVertexCycle();

            // cout << "WARNING! This assumes that all faces are identical and regular." << endl;

            for( size_t i = 0; i < m_OuterVertexCycle.size(); i++ )
            {
                size_t degree = FindEdgesWithVertex(m_OuterVertexCycle[i]).size();
                m_PolygonCompareString.push_back(ALPHABET[degree]);
            }

            // cout << "Compare String : " << m_PolygonCompareString << endl;
        }

        template<class MatrixType>
        void ComputeCenterOfFace(const size_t& face, MatrixType& center )
        {
            double sumx = 0, sumy = 0;
            for(size_t f = 0; f < m_Faces[face].size(); f++)
            {
                sumx += m_Vertices[m_Faces[face][f]][0];
                sumy += m_Vertices[m_Faces[face][f]][1];
            }
            center[0] = sumx/double(m_Faces[face].size());
            center[1] = sumy/double(m_Faces[face].size());
        }

        void ComputeCentroid()
        {
            // The centroid can be found by triagulation this means for regular polygons we can
            // compute the weighted average of the center of each face. wieght is the area fraction.
        }

        void ComputeConvexhull()
        {

        }

        double FindMaxDistance() const
        {
            double maxdist = 0.0;
            return maxdist;
        }

        double FindMaxDistFromOrigin() const  // from centroid?
        {
            double maxdist = 0.0;
            return maxdist;
        }

    // Accessor Methods
        bool        IsInitialized() const { return m_bInitialized; }

        size_t      NumFaces() const { return m_Faces.size(); }

        size_t      NumEdges() const { return m_Edges.size(); }

        size_t      NumHinges() const { size_t ct = 0; for(size_t i = 0; i < m_bHinge.size(); i++) {ct += size_t(m_bHinge[i]);} return ct;}

        size_t      NumFreeEdges() const {return NumEdges() - NumHinges();}

        size_t      NumVertices() const { return m_Vertices.size(); }

        const vector<vector<size_t> >& GetFaces() const { return m_Faces; }

        const vector<vector<size_t> >& GetEdges() const { return m_Edges; }

        const vector<vector<double> >& GetVerts() const { return m_Vertices; }

        const vector<bool>& GetHingeFlags() const { return m_bHinge; }

        void ComputeFreeEdgeLabels()
        {
            if(!m_FreeEdgeLabels.size())
            {
                for(size_t i = 0; i < m_bHinge.size(); i++)
                    if(!m_bHinge[i])
                        m_FreeEdgeLabels.push_back(i);
            }
        }

        size_t GetFreeEdgeLabel(const size_t& n) const
        {
            // returns the id of the n-th free edge.
            // requires ComputeFreeEdgeLabels() to be called once prior.
            // ComputeFreeEdgeLabels();
            size_t label = invalid_index();
            if(n < m_FreeEdgeLabels.size())
                label = m_FreeEdgeLabels[n];
            return label;
        }

        size_t GetFreeEdgeIndex(const size_t& label) const
        {
            // returns the index of the free edge with label
            // requires ComputeFreeEdgeLabels() to be called once prior.
            // ComputeFreeEdgeLabels();
            for(size_t n = 0; n < m_FreeEdgeLabels.size(); n++)
                if(label == m_FreeEdgeLabels[n])
                    return n;
            return invalid_index();
        }

        size_t GetHingeIndex(const size_t& label) const
        {
            if(m_bHinge[label])
            {
                size_t ct = 0;
                for(size_t i = 0; i < label; i++)
                    if(m_bHinge[i])
                        ct++;
                return ct;
            }
            return invalid_index();
        }

        std::vector< std::vector<size_t> > HingeMapping() const
        {
            std::vector< std::pair<size_t, size_t> > labels;
            std::vector< std::vector<size_t> > map, edge_map;
            edge_map = EdgeMapping();
            for(size_t i = 0; i < m_bHinge.size(); i++)
            {
                if(m_bHinge[i])
                {
                    labels.push_back(pair<size_t, size_t>(i, labels.size()));
                }
            }

            for (size_t i = 0; i < edge_map.size(); i++)
            {
                vector<size_t> v;
                for( size_t j = 0; j < labels.size(); j++)
                {
                    size_t ndx = 0;
                    for (size_t k = 0; k < labels.size(); k++)
                    {
                        if ( labels[k].first == edge_map[i][labels[j].first])
                        {
                            ndx = k;
                        }
                    }
                    v.push_back(labels[ndx].second);
                }
                map.push_back(v);
            }
//############################################################
//            for (size_t i = 0; i < map.size(); i++)
//            {
//                for( size_t j = 0; j < map[i].size(); j++)
//                {
//                    std::cout << map[i][j] << ", ";
//                }
//                std::cout << std::endl;
//            }
//############################################################
            return map;
        }

        std::vector< std::vector<size_t> > FreeEdgeMapping() const
        {
            std::vector< std::pair<size_t, size_t> > labels; // should have used std::map.
            std::vector< std::vector<size_t> > map, edge_map;
            edge_map = EdgeMapping();
            for(size_t i = 0; i < m_bHinge.size(); i++)
            {
                if(!m_bHinge[i])
                {
                    labels.push_back(pair<size_t, size_t>(i, labels.size()));
                }
            }

            for (size_t i = 0; i < edge_map.size(); i++)
            {
                vector<size_t> v;
                for( size_t j = 0; j < labels.size(); j++)
                {
                    size_t ndx = 0;
                    for (size_t k = 0; k < labels.size(); k++)
                    {
                        if ( labels[k].first == edge_map[i][labels[j].first])
                        {
                            ndx = k;
                        }
                    }
                    v.push_back(labels[ndx].second);
                }
                map.push_back(v);
            }
//############################################################
//            for (size_t i = 0; i < map.size(); i++)
//            {
//                for( size_t j = 0; j < map[i].size(); j++)
//                {
//                    std::cout << map[i][j] << ", ";
//                }
//                std::cout << std::endl;
//            }
//############################################################
            return map;
        }
        void ComputeIsomorphisms()
        {
            if(m_Isomorphisms.size() == 0)
            {
                Network::CIsomorphism<> matcher(&m_VertexGraph, &m_VertexGraph);
                matcher.ComputeMatching();
                if(matcher.NumMappings())
                {
                    m_Isomorphisms.resize(matcher.NumMappings());
                    m_Isomorphisms.insert(m_Isomorphisms.begin(), matcher.begin(), matcher.end());
                }
            }
        }

        std::vector< std::vector<size_t> > EdgeMapping() const
        {
            // std::cout << "N -- EdgeMapping() -- " << std::endl;
            std::vector< std::vector<size_t> > edge_map;

            size_t niso = (m_bIgnoreSymmetry && m_Isomorphisms.size()) ? 1 : m_Isomorphisms.size();
            // std::cout << " niso =  "<< niso << std::endl;
            for(size_t iso = 0; iso < niso; iso++)
            {
                vector<size_t> map;
                for(size_t e = 0; e < m_Edges.size(); e++)
                {

                    vector<size_t> v1 = FindEdgesWithVertex(m_Isomorphisms[iso][m_Edges[e][0]]);
                    vector<size_t> v2 = FindEdgesWithVertex(m_Isomorphisms[iso][m_Edges[e][1]]);

                    vector<size_t> intersection = utils::math_util::Intersection(v1, v2);

                    if ( intersection.size() == 1)
                    {
                        map.push_back(intersection[0]);
                    }
                    else
                    {
                        std::cout << "error! inconsistent json!" << std::endl;
                        throw 0x000002; // need to throw a real error.
                    }
                }
                edge_map.push_back(map);
            }
//############################################################
//            for (size_t i = 0; i < edge_map.size(); i++)
//            {
//                for( size_t j = 0; j < edge_map[i].size(); j++)
//                {
//                    std::cout << edge_map[i][j] << ", ";
//                }
//                std::cout << std::endl;
//            }
//            std::cout << "X -- EdgeMapping() -- " << std::endl;
//############################################################

            return edge_map;
        }

        // const BASE_POLYHEDRA*           GetShape() const { return m_pShape; }

        const Network::CNetwork* GetVertexGraph() const  { return &m_VertexGraph; }

        Network::CNetwork* GetVertexGraph()  { return &m_VertexGraph; }

        Network::CNetwork* GetFaceGraph() { return &m_FaceGraph; }

        const Network::CNetwork* GetFaceGraph() const { return &m_FaceGraph; }

        const Network::CNetwork* GetCuttingGraph() const { return &m_CuttingGraph; }

        size_t GetId() const { return m_Id; }

        const string& GetParentPolyhedron() const { return m_ParentPolyhedron; }

        const vector<size_t>& GetOuterVertexCyle() const  { return m_OuterVertexCycle; }

        const string& GetCompareString()const { return m_PolygonCompareString; }

        const vector< vector< size_t > >& GetGluing() const { return m_Gluing; }

        const vector< vector< size_t > >& GetIsomorphisms() const { return m_Isomorphisms; }

        const CNetInfo& GetNetInfo() const { return m_NetInfo; }

        bool GetIgnoreSymmetryFlag() const { return m_bIgnoreSymmetry; }

        const vector<unfold_polyhedra::CDihedralAngle>& GetDihedrals() const { return m_DihedralAngles; }

        string GetNameString() const
        {
            string name;
            if (m_ParentPolyhedron == "Tetrahedron") {
                name = m_ParentPolyhedron+"Net"+utils::NumberToString(m_Id, 1, '0');
            }
            else if (m_ParentPolyhedron == "Cube") {
                name = m_ParentPolyhedron+"Net"+utils::NumberToString(m_Id, 2, '0');
            }
            else if (m_ParentPolyhedron == "Octahedron") {
                name = m_ParentPolyhedron+"Net"+utils::NumberToString(m_Id, 2, '0');
            }
            else if (m_ParentPolyhedron == "Dodecahedron") {
                name = m_ParentPolyhedron+"Net"+utils::NumberToString(m_Id, 5, '0');
            }
            else if (m_ParentPolyhedron == "Icosahedron") {
                name = m_ParentPolyhedron+"Net"+utils::NumberToString(m_Id, 5, '0');
            }
            else
            {
                name = m_ParentPolyhedron+"Net"+utils::NumberToString(m_Id);
            }

            return name;
        }

        vector<double> GetHingeDihedrals()
        {
            vector<double> ret;
            for(size_t h = 0; h < m_bHinge.size(); h++)
            {
                if(m_bHinge[h])
                {
                    vector<size_t> fcs = FindFacesWithEdge(h);
                    for(size_t d = 0; d < m_DihedralAngles.size(); d++)
                    {
                        if(utils::IsInVec(fcs, m_DihedralAngles[d].GetFace0()) && utils::IsInVec(fcs, m_DihedralAngles[d].GetFace1()))
                        {
                            ret.push_back(m_DihedralAngles[d].GetAngle());
                        }
                    }
                }
            }
            return ret;
        }

        //void SetShape(BASE_POLYHEDRA* pShape) { m_pShape = pShape; }

        void SetFaces(const vector<vector<size_t> >& faces, bool bSetEdges = true )
        {
            m_Faces = faces;
            if(bSetEdges)
            {
                for(size_t f = 0; f < m_Faces.size(); f++){
                    for(size_t g = 0; g < m_Faces[f].size(); g++){
                        vector<size_t> edge;
                        edge.push_back(min(m_Faces[f][g], m_Faces[f][(g+1) % m_Faces[f].size()]));
                        edge.push_back(max(m_Faces[f][g], m_Faces[f][(g+1) % m_Faces[f].size()]));
                        utils::PushUnique(m_Edges, edge);
                    }
                }
            }
        }

        void SetVertices(const vector<vector<double> >& verts) { m_Vertices = verts; }

        void SetDihedrals(const vector<unfold_polyhedra::CDihedralAngle>& adata) { m_DihedralAngles = adata; }

        void SetVertexGraph(Network::CNetwork& g) { m_VertexGraph = g; }

        void SetFaceGraph(Network::CNetwork& g) { m_FaceGraph = g; }

        void SetCuttingGraph(Network::CNetwork& g) { m_CuttingGraph = g; }

        void AddGluedEdgePair(vector<size_t>& pair) { m_Gluing.push_back(pair); }

        void SetNetId(size_t ndx) { m_Id = ndx;}

        void SetPolyhedron(const string& poly) { m_ParentPolyhedron = poly;}
#ifdef c_plus_plus_11
        void SetDihedralsFromGluingEdges(size_t edge1, size_t edge2, std::vector<double>& dihedrals, std::vector<int>& cmp) const
        {
            size_t first, second;
            if(edge1 < edge2)
            {
                first = edge1;
                second = edge2;
            }
            else if( edge2 < edge1 )
            {
                first = edge2;
                second = edge1;
            }
            else
            {
                return; // in the case that edge1 == edge2 nothing needs to be done.
            }
            std::vector<size_t> edge({first, second});
            std::vector<size_t> redge({second, first});
            if(!utils::IsInVec(m_Gluing, edge) && !utils::IsInVec(m_Gluing, redge))
            {
                std::cout << "edge1 = " << edge1 << ", edge2 = " << edge2 << endl;
                throw std::runtime_error("I don't know how to glue those edges");
            }
            vector<size_t> fc1, fc2;
            fc1 = FindFacesWithEdge(first);
            fc2 = FindFacesWithEdge(second);
            if(fc1.size() != 1 || fc2.size() != 1)
                throw std::runtime_error("trying to glue a hinge or an edge that does not belong to a face");


            std::vector< std::vector< size_t > > paths = m_FaceGraph.ShortestPath(fc1[0], fc2[0]);
            if(!paths.size())
                throw std::runtime_error("could not find a path between faces.");
            std::vector<size_t> edgesFound;
            std::set<size_t> hinges;
            for(size_t i = 0; i < paths[0].size()-1; i++) // there should only be one.
            {
                // std::cout << "checking path "<< paths[0][i] << " --> " << paths[0][i+1] << endl;
                edgesFound = FindCommonEdgesInFaces(paths[0][i], paths[0][i+1]);
                hinges.insert(edgesFound.begin(), edgesFound.end());
            }
            for(std::set<size_t>::iterator it = hinges.begin(); it != hinges.end(); it++ )
            {
                std::vector<size_t> faces = FindFacesWithEdge(*it);
                // std::cout << "checkpoint - 1/2 -- " << faces.size() << endl;
                for(size_t i = 0; i < m_DihedralAngles.size() && faces.size()==2; i++)
                {
                    if  (   (m_DihedralAngles[i].GetFace0() == faces[0] && m_DihedralAngles[i].GetFace1() == faces[1]) ||
                            (m_DihedralAngles[i].GetFace0() == faces[1] && m_DihedralAngles[i].GetFace1() == faces[0])
                        )
                    {
                        size_t ndx = GetHingeIndex(*it);
                        // std::cout << "ndx = " << ndx << ", edge = " << *it << std::endl;
                        if(ndx != invalid_index() && ndx < dihedrals.size())
                        {
                            dihedrals[ndx] = m_DihedralAngles[i].GetAngle();
                            cmp[ndx] = 1;
                        }
                        else
                            throw std::runtime_error("could not get the hinge index!");
                    }
                }
            }
        }
#endif
        void PrintNet() const
        {
            printf("-- %s %lu --\n", m_ParentPolyhedron.c_str(), m_Id);
            printf("Vertices(%lu): \n", m_Vertices.size());
            for(size_t v = 0; v < m_Vertices.size(); v++)
                printf("%7.3f, %7.3f \n", m_Vertices[v][0], m_Vertices[v][1]);

            printf("Faces(%lu): \n", m_Faces.size());
            for(size_t f = 0; f < m_Faces.size(); f++)
                for(size_t i = 0; i < m_Faces[f].size(); i++)
                    printf("%3lu%s", m_Faces[f][i], (i ==m_Faces[f].size()-1 ?  "\n" : ", "));

            printf("Edges(%lu): \n", m_Edges.size());
            for(size_t e = 0; e < m_Edges.size(); e++)
                for(size_t i = 0; i < m_Edges[e].size(); i++)
                    printf("%3lu%s", m_Edges[e][i], (i == m_Edges[e].size()-1 ?  (m_bHinge.size() ? (m_bHinge[e] ? "  True\n" : "  False\n") : "  False\n") : ", ")); //
        }
        /*
        *   Should actually verify this but basically what happens is that mathematica return the vertices of the faces in counter clockwise order (from outward
        *   pointing normal) we will then use this to map the edges and vertices. (all we have to do is reverse the order.)
        */
        void CreateNewFace( const size_t& fcndx,
                            const vector< vector<double> >& FaceShape,
                            const size_t& EdgeId,
                            unfold_polyhedra::CPolyhedron* pPoly,
                            map<size_t, vector<size_t> >& vertex_map,
                            map<size_t, vector<size_t> >& edge_map )
        {
            vector<size_t> face, edge, ndxMap;
            vector<int> polyface;
            polyface = pPoly->GetFaces()[fcndx];
            // Create a new face @ the edge corresponding to EdgeId
            if( EdgeId < NumEdges())
            {
                vector<size_t> nfc = FindFacesWithEdge(EdgeId);
                bool bLeft = false;

                EigenMatrixHelper::Vector edge_vec(2), newEdge(2), translate(2), center(2), cross(2), cmx(2);
                EigenMatrixHelper::matrix newCoords(2, FaceShape.size()), rot(2,2);

                ComputeCenterOfFace(nfc[0], center);
                cmx[0] = center[0] - m_Vertices[m_Edges[EdgeId][0]][0]; // assume the reference is the first ndx
                cmx[1] = center[1] - m_Vertices[m_Edges[EdgeId][0]][1];
                edge_vec(0) = m_Vertices[m_Edges[EdgeId][1]][0] - m_Vertices[m_Edges[EdgeId][0]][0];
                edge_vec(1) = m_Vertices[m_Edges[EdgeId][1]][1] - m_Vertices[m_Edges[EdgeId][0]][1];

                double z = cmx[0]*edge_vec[1] - cmx[1]*edge_vec[0];
                bLeft  = z < 0;
                // cout << " face = " << m_Faces.size() << " edge = " << EdgeId << " Flip =" <<boolalpha  << bLeft<< endl;
                newEdge(0) = FaceShape[1][0] - FaceShape[0][0];
                newEdge(1) = FaceShape[1][1] - FaceShape[0][1];

                size_t ref = 0;
                if(bLeft){ // the center is to the left
                    ref = 1;
                    edge_vec(0) = m_Vertices[m_Edges[EdgeId][0]][0] - m_Vertices[m_Edges[EdgeId][1]][0];
                    edge_vec(1) = m_Vertices[m_Edges[EdgeId][0]][1] - m_Vertices[m_Edges[EdgeId][1]][1];
                    face.push_back(m_Edges[EdgeId][1]);
                    face.push_back(m_Edges[EdgeId][0]);

                }
                else
                {
                    face.push_back(m_Edges[EdgeId][0]);
                    face.push_back(m_Edges[EdgeId][1]);
                }

                translate(0) = m_Vertices[m_Edges[EdgeId][ref]][0];
                translate(1) = m_Vertices[m_Edges[EdgeId][ref]][1];

                for(size_t c = 0; c < FaceShape.size(); c++)
                {
                    newCoords(0, int(c)) = FaceShape[c][0];
                    newCoords(1, int(c)) = FaceShape[c][1];
                }

                double cosine = newEdge.dot(edge_vec);
                double z2 = newEdge[0]*edge_vec[1] - newEdge[1]*edge_vec[0];
                bool bLeftAgain  = z2 < 0;
//                cout << "cos = " << cosine << "acos = " << acos(cosine) << endl;
                double theta = bLeftAgain ? M_PI*2.0 - utils::arccos(cosine) : utils::arccos(cosine);
//                cout << " face = " << m_Faces.size() << " edge = " << EdgeId << " Flip =" <<boolalpha  << bLeft << " Theta = " << theta << " Double Flip = "<< bLeftAgain << endl;

                cosine = cos(theta);
                double sine = sin(theta);

                rot(0,0) = cosine;  rot(0,1) = -sine;
                rot(1,0) = sine;    rot(1,1) = cosine;
                newCoords = rot*newCoords;
//                cout << "Rot Coords : " << endl << newCoords << endl;
                for(int i = 0; i < newCoords.cols(); i++) { newCoords(0, i)+=translate(0); newCoords(1, i)+=translate(1); }

//                cout << "New Coords : " << endl << newCoords << endl;
//                cout << "Edge : " << endl << "(" << m_Vertices[face[0]][0] <<", "<< m_Vertices[face[0]][1] << ")" << endl
//                                          << "(" << m_Vertices[face[1]][0] <<", "<< m_Vertices[face[1]][1] << ")" << endl;

                utils::float_is_equal<double, -2> compare_func;
                assert(compare_func(m_Vertices[face[0]][0], newCoords(0,0)) && compare_func(m_Vertices[face[0]][1], newCoords(1,0)));
                assert(compare_func(m_Vertices[face[1]][0], newCoords(0,1)) && compare_func(m_Vertices[face[1]][1], newCoords(1,1)));

                for(int i = 2; i < newCoords.cols(); i++)
                {
                    face.push_back(NumVertices());
                    vector<double> temp;
                    temp.push_back(newCoords(0,i));
                    temp.push_back(newCoords(1,i));
                    m_Vertices.push_back(temp);
                }

                map<size_t, vector<size_t> >::iterator iter = vertex_map.begin();
                size_t vndx1 = invalid_index(), vndx2 = invalid_index();
                for( ; iter != vertex_map.end(); iter++)
                {
                    if(utils::IsInVec(iter->second, face[0])) vndx1 = iter->first;
                    if(utils::IsInVec(iter->second, face[1])) vndx2 = iter->first;
                }

//                cout << "ndx1 = "<< vndx1 << " ndx2 = "<< vndx2 << endl;
                size_t n1, n2;
                n1 = utils::FindInVec(polyface, int(vndx1));
                n2 = utils::FindInVec(polyface, int(vndx2));


                assert(vndx1 != vndx2 && vndx1 != invalid_index() && vndx2 != invalid_index() && n1 != n2);
                int n = abs(int(n2)-int(n1)) > 1 ? n2 > n1 ? -1 : 1 : int(n2)-int(n1);
                for( size_t i = 0; i < face.size(); i++)
                {
                    int m = (int(n1) + int(i)*n) % int(polyface.size());
                    size_t imap =  size_t(m < 0 ? m + int(polyface.size()) : m);
                    assert(imap < polyface.size());
                    ndxMap.push_back(imap);
//                    cout << face[i] << " -> " << polyface[imap] << " i = "<< i << " imap = "<< imap << endl;
                    utils::PushUnique(vertex_map[size_t(polyface[imap])], face[i]);
                }
            }
            else // create new default face.
            {
                for(size_t i = 0; i < FaceShape.size(); i++)
                {
                    size_t fndx1 = size_t((int(polyface.size()) - int(i)-1 < 0) ?   int(i)-1 : int(polyface.size()) - int(i)-1);
                    ndxMap.push_back(fndx1);
                    vertex_map[size_t(polyface[fndx1])].push_back(m_Vertices.size());
                    face.push_back(m_Vertices.size());
                    m_Vertices.push_back(FaceShape[i]);
                }
            }

            for(size_t i = 0; i < face.size(); i++)
            {
                if(face[i] < face[(i+1) % face.size()])
                {
                    edge.push_back(face[i]);
                    edge.push_back(face[(i+1) % face.size()]);
                }
                else
                {
                    edge.push_back(face[(i+1) % face.size()]);
                    edge.push_back(face[i]);
                }
                utils::PushUnique(m_Edges, edge);
                edge.clear();

//                cout << " verts("<< ndxMap[i]<< ", " << ndxMap[(i+1)%ndxMap.size()]<< ") = " << polyface[ndxMap[i]] << ", "<<  polyface[ndxMap[(i+1)%ndxMap.size()]]<< endl;
                size_t eid = pPoly->FindEdge(polyface[ndxMap[i]], polyface[ndxMap[(i+1)%ndxMap.size()]]);
//                cout << "Edge "<< eid << " : " << polyface[ndxMap[i]] << ", "<<  polyface[ndxMap[(i+1)%ndxMap.size()]]<< endl;
//                cout << "mapping edge " << FindEdge(face[i], face[(i+1) % face.size()]) << " -> "<< eid << endl;

                utils::PushUnique(edge_map[eid], FindEdge(face[i], face[(i+1) % face.size()]));

            }
            m_Faces.push_back(face);
        }

    private:
        bool IsEqual(const CNet& cmp) const
        {
            // must convert the polygon to a string and comapare a string in a circular fasion.
            bool bEqual = false;
            utils::circular_string_equal_to eq;
            if(!cmp.GetCompareString().length() || !m_PolygonCompareString.length() )
            {
                cout << "Must compute string before comparing nets!"  << endl;
                bEqual = false;
            }
            else
            {
                bEqual =        eq(m_PolygonCompareString, cmp.GetCompareString(), m_PolygonCompareString.length())
                            ||  eq(utils::InvertString(m_PolygonCompareString), cmp.GetCompareString(), m_PolygonCompareString.length());

            }
            return bEqual;
        }

        CNet& CopyFrom(const CNet& src)
        {
            m_Id = src.GetId();

            m_ParentPolyhedron = src.GetParentPolyhedron();
            m_OuterVertexCycle = src.GetOuterVertexCyle();
            m_PolygonCompareString = src.GetCompareString();

            m_VertexGraph = *src.GetVertexGraph();
            m_FaceGraph = *src.GetFaceGraph();
            m_CuttingGraph = *src.GetCuttingGraph();

            m_Faces = src.GetFaces();
            m_Edges = src.GetEdges();
            m_Vertices = src.GetVerts();
            m_bHinge = src.GetHingeFlags();

            m_Gluing = src.GetGluing();
            m_Isomorphisms = src.GetIsomorphisms();

            m_NetInfo = src.GetNetInfo();
            m_DihedralAngles = src.GetDihedrals();

            m_bInitialized = true;
            return *this;
        }

    public:
        bool operator == (const CNet& cmp) const { return IsEqual(cmp); }

        CNet& operator = (const CNet& src) { return CopyFrom(src); }

    private:
        bool                        m_bInitialized;                 // true if the net has been initialized
        bool                        m_bIgnoreSymmetry;
        size_t                      m_Id;                           // Id of the net
        string                      m_ParentPolyhedron;             // Name of the polyhedron

        vector<size_t>              m_OuterVertexCycle;             // path around the polygon used to build the compare string
        string                      m_PolygonCompareString;         // string representing the interior angle at each vertex in the cycle

        // BASE_POLYHEDRA*         m_pShape;
        Network::CNetwork           m_VertexGraph;                  // Connectivity of the vertices
        Network::CNetwork           m_FaceGraph;                    // Connectivity of the faces
        Network::CNetwork           m_CuttingGraph;                 // Cutting graph from which the net was derived

        vector<vector<size_t> >     m_Faces;                        // Holds the id's of the vertices that belong to each face
        vector<vector<size_t> >     m_Edges;                        // Holds the id's of the vertices that belong to each edge
        vector<vector<double> >     m_Vertices;                     // Coordinates of the vertices
        vector< bool >              m_bHinge;                       // i-th element is true if the i-th edge is shared by 2 faces

        vector< vector<size_t> >    m_Gluing;                      // Holds the edge-id pairs for edges that are "glued" together when folded.
        vector< vector<size_t> >    m_Isomorphisms;                // All the automorphisms for the net (rotational symmetry) computed from python.

        CNetInfo                    m_NetInfo;                      // Other information calculated from the net

        vector<unfold_polyhedra::CDihedralAngle> m_DihedralAngles; // dihedral angle between faces of the folded polyhedron

        std::vector<size_t>         m_FreeEdgeLabels;

    // Values That I Still need
    // 1. Gluing -- what edges need to be bonded to fold the polyhedron (check!)
    // 2. Ideal states/network of states.
    //      a. update a centralized library with the states used.
    //          i. to compare those states I will need to write an Graph Isomorphism finder??
    //      b. build the network
    // 3. Hinges            (check!)
    // 4. Cutting Graph     (check!)
};

};



































#endif
#endif
