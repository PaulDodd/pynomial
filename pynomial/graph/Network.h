//
//  Network.h
//  SharedFolding
//
//  Created by Paul M Dodd on 11/9/13.
//  Copyright (c) 2013 Paul Dodd. All rights reserved.
//

// TODOs:
//  2. Random network.
//      --> mark newman book for details.
//
//

#ifndef SharedFolding_Network_h
#define SharedFolding_Network_h

// #include "SharedInclude.h"
#include <iostream>
#include <vector>
#include <map>
#include <set>
#include <stack>
#include <stdexcept>
#include <eigen3/Eigen/Sparse>
#include "Matrix.h"

#define INVALID_DISTANCE    (-1)
namespace pynomial{
namespace graph{
// using namespace std; // TODO: remove this
// using json::CJSONValueObject;
//using namespace Eigen;

enum NetworkError{
    ErrorSuccess    = 0x00000000,
    ErrorGeneral    = 0x00000001
};

enum NetworkType{
    General     = 0x00000000,
    Undirected  = 0x00000001,
    Directed    = 0x00000002,
    Weighted    = 0x00000004,
    Simple      = 0x00000008,
    Bipartite   = 0x00000010,
    Multigraph  =(0x00000020 | Weighted)
};


#define JSON_NAME_ADJMAT        "AdjMat"
#define JSON_NAME_NETWORK_TYPE  "NetworkType"

using std::vector;
using std::string;
using std::queue;
using std::pair;

class AdjacencyList
{
    private:
        // make so these lists are not copyable.
        AdjacencyList( const AdjacencyList& );
        AdjacencyList& operator = (const AdjacencyList&);

    public:
        AdjacencyList(const EigenMatrixHelper::spMatrix& adjMatrix, bool bDirected) : m_bUpdate(true), m_bPredecessors(bDirected), m_AdjMatrix(adjMatrix) { }

        ~AdjacencyList() {}

        void UpdateCheck(size_t src, size_t tgt)
        {
            m_bUpdate = (m_bUpdate || !m_Successors.size() ? true : !std::binary_search(m_Successors[src].begin(), m_Successors[src].end(), tgt) && (m_bPredecessors ? !std::binary_search(m_Predecessors[tgt].begin(), m_Predecessors[tgt].end(), src) : true));
        }

        void MarkToUpdate() { m_bUpdate = true; }

        bool ShouldUpdate() const { return m_bUpdate; }

        bool Update()
        {
            if(!ShouldUpdate())
                return false;

            BuildLists();

            return true;
        }

        const std::vector< size_t >& Successors(size_t i) const
        {
            AssertValidState();
            AssertSuccessorBounds(i);
            return m_Successors[i];
        }

        const std::vector< double >& Weights(size_t i) const
        {
            AssertValidState();
            AssertSuccessorBounds(i);
            return m_Weights[i];
        }

        const std::vector< size_t >& Predecessors(size_t i) const
        {
            AssertValidState();
            AssertPredecessorBounds(i);
            return m_Predecessors[i];
        }

        bool DoesEdgeExist(size_t src, size_t tgt) const
        {
            AssertValidState();
            AssertSuccessorBounds(src);
            return std::binary_search(m_Successors[src].begin(), m_Successors[src].end(), tgt);
        }

        void Reset() { m_Successors.clear(); m_Predecessors.clear(); m_Weights.clear(); m_bUpdate = true; }

        void PrintListStats() const
        {
            std::cout << "num successors.... " << m_Successors.size() << std::endl;
            std::cout << "num predecessors.. " << m_Predecessors.size() << std::endl;
            std::cout << "is directed....... " << (m_bPredecessors ? "true" : "false") << std::endl;
            std::cout << "needs update...... " << (m_bUpdate ? "true" : "false") << std::endl;
            std::cout << "adj mat addr...... " << std::hex << &m_AdjMatrix << std::dec << std::endl;
            std::cout << "adj mat rows...... " << m_AdjMatrix.rows() << std::endl;
            std::cout << "adj mat cols...... " << m_AdjMatrix.cols() << std::endl;
        }

    private:
        void AssertValidState() const
        {
            if( ShouldUpdate() )
                throw std::runtime_error("error: adjacency list is not in a valid state!");
        }

        void AssertPredecessorBounds(const size_t& i) const
        {
            AssertValidState();
            if( i < m_Predecessors.size() )
            {
                std::stringstream ss;
                ss << "error: trying to access node out of range of adjacency list! " << i << " >= " << m_Predecessors.size();
                throw std::runtime_error(ss.str());
            }

        }

        void AssertSuccessorBounds(const size_t& i) const
        {
            AssertValidState();
            if( i >= m_Successors.size() )
            {
                std::stringstream ss;
                ss << "error: trying to access node out of range of adjacency list! " << i << " >= " << m_Successors.size();
                throw std::runtime_error(ss.str());
            }
        }

        void BuildLists()
        {
            // overall complexity is O(m + n*k*log(k))
            // clear the lists.
            Reset();
            m_Successors.resize(size_t(m_AdjMatrix.rows()), std::vector<size_t>());
            m_Weights.resize(size_t(m_AdjMatrix.rows()), std::vector<double>());

            if(m_bPredecessors)
                m_Predecessors.resize(size_t(m_AdjMatrix.rows()), std::vector<size_t>());


            // these loops will be O(m) number of edges.
            for (int k = 0; k < m_AdjMatrix.outerSize(); ++k)
            {
                for (EigenMatrixHelper::spMatrix::InnerIterator it(m_AdjMatrix,k); it; ++it)
                {
                    size_t col = size_t(it.col()), row = size_t(it.row());
                    m_Successors[col].push_back(row);
                    // m_Weights[col].push_back(it.value());

                    if(m_bPredecessors)
                        m_Predecessors[row].push_back(col);
                }
            }


            // sort the list since this is not guaranteed. O(n*k*log(k))
            for(size_t i = 0; i < m_Successors.size(); i++)
            {
                std::sort(m_Successors[i].begin(), m_Successors[i].end());
                if(m_bPredecessors)
                    std::sort(m_Predecessors[i].begin(), m_Predecessors[i].end());
                m_Weights[i].reserve(m_Successors[i].size());
                for(size_t j = 0; j < m_Successors[i].size(); j++)
                {
                    m_Weights[i].push_back(m_AdjMatrix.coeff(m_Successors[i][j],i)); // this is a bit expensive.
                }
            }


            m_bUpdate = false;
        }

    private:
        bool                m_bUpdate;              //! true if we need to recompute the lists
        bool                m_bPredecessors;        //! true if we need to keep track of the predecessors list (directed graph)

        std::vector< std::vector<size_t> > m_Successors;
        std::vector< std::vector<double> > m_Weights;
        std::vector< std::vector<size_t> > m_Predecessors;

        const EigenMatrixHelper::spMatrix& m_AdjMatrix;
};


//Implementation right now is for general networks only and not for bipartite.
class CNetwork //: public json::CJSONValueObject<CNetwork>
{
    public:
        class edge_info
        {
            public:
                edge_info() : source(-1), target(-1), weight(0) { }

                edge_info(Eigen::SparseMatrix<double>::InnerIterator& it) : source(-1), target(-1), weight(0) { CopyFrom(it); }

                edge_info(const edge_info& src) : source(-1), target(-1), weight(0) { CopyFrom(src); }

                edge_info& CopyFrom(Eigen::SparseMatrix<double>::InnerIterator& it)
                {
                    if(it)
                    {
                        source = int(it.col());
                        target = int(it.row());
                        weight = it.value();
                    }
                    return *this;
                }

                edge_info& CopyFrom(const edge_info& src)
                {
                    source = src.source;
                    target = src.target;
                    weight = src.weight;
                    return *this;
                }

                edge_info& operator=(const edge_info& src) { return CopyFrom(src); }

                bool operator<(const edge_info& cmp) const { LessThan lt; return lt(*this, cmp); }

                bool operator>(const edge_info& cmp) const { GreaterThan gt; return gt(*this, cmp); }
            // public members.
                int source;
                int target;
                double weight;
            public:
                struct LessThan : std::binary_function<edge_info, edge_info, bool>{ bool operator()(const edge_info& lhs, const edge_info& rhs){ return lhs.weight < rhs.weight; }};
                struct GreaterThan : std::binary_function<edge_info, edge_info, bool>{ bool operator()(const edge_info& lhs, const edge_info& rhs){ return lhs.weight > rhs.weight; }};
        };

        template<class TVal>
        class iterator_base : public std::iterator<std::forward_iterator_tag, TVal>
        {
            public:
                iterator_base(const CNetwork* pNetwork) : m_pNetwork(pNetwork) { }
                ~iterator_base() { }
            protected:
                const CNetwork* const m_pNetwork;
        };

        class edge_iterator : public iterator_base<edge_info>
        {
            public:
                edge_iterator(const CNetwork* pNetwork) : iterator_base<edge_info>(pNetwork), m_OuterIndex(0), m_MatrixIter(NULL)
                {

                    NextIter();
                    if(m_MatrixIter)
                    {
                        NextEdge(false);
                    }

                }
                ~edge_iterator()
                {
                    DestroyIter();
                }

                void DestroyIter()
                {
                    if(m_MatrixIter)
                    {
                        delete m_MatrixIter;
                        m_MatrixIter = NULL;
                    }
                }

                void NextEdge(bool bIncrement = true)
                {
                    // std::cout << "N NextEdge" << std::endl;
                    AdvanceIter(bIncrement); // set it up for the next one.
                    if(m_MatrixIter)
                    {
                        m_edge.CopyFrom((*m_MatrixIter));
                    }
                    else
                    {
                        m_edge = edge_info(); // sets to the invalid values.
                    }

                }

                void NextIter()
                {
                    DestroyIter();
                    if(m_OuterIndex < size_t(m_pNetwork->GetAdjMatrix().outerSize()))
                    {
        //                std::cout << "processing outer index " << m_OuterIndex << " of "<< m_pNetwork->GetAdjMatrix().outerSize() << std::endl;
        //                std::cout << m_pNetwork->GetAdjMatrix().toDense()<< std::endl;
                        m_MatrixIter = new Eigen::SparseMatrix<double>::InnerIterator(m_pNetwork->GetAdjMatrix(), int(m_OuterIndex++));
        //                std::cout << "done " << std::endl;
                    }
                }

                void AdvanceIter(bool bIncrement = true)
                {
                    bool bAdvanceIterator = m_MatrixIter;
                    while ( bAdvanceIterator && m_MatrixIter )
                    {
                        bAdvanceIterator = bIncrement ? !(++(*m_MatrixIter)) : !(*m_MatrixIter);
                        // std::cout << "valid iter = " << boolalpha << bAdvanceIterator << " incremented = "<< bIncrement << std::endl;
                        bIncrement = true;
                        if(!bAdvanceIterator) // valid iterator
                        {
                            if(m_pNetwork->IsUndirected())
                            {
                                bAdvanceIterator = m_MatrixIter->row() <= m_MatrixIter->col(); // just need the lower triangle.
                            }

                        }
                        else
                        {
                            NextIter();
                            bIncrement = false;
                        }
                    }
                }

                edge_iterator& Next()
                {
                    NextEdge();
                    return *this;
                }

                edge_iterator& operator ++ ( int ) // unused variable doesn't get a name.
                {
                    // std::cout << "N ++" << std::endl;
                    return Next();
                }

                edge_info* operator -> ()
                {
                    return &m_edge;
                }

                const edge_info& operator*()
                {
                    return m_edge;
                }

                bool IsValid()
                {
                    // std::cout << " is valid: " << boolalpha <<(m_edge.source != -1 && m_edge.target != -1) << std::endl;
                    return (m_edge.source != -1 && m_edge.target != -1);
                }

            private:
                size_t      m_OuterIndex;
                Eigen::SparseMatrix<double>::InnerIterator* m_MatrixIter;
                edge_info   m_edge;
        };

        class vertex_iterator
        {
            public:
                enum AlgorithmType {BreadthFirstSearch, DepthFirstSearch};
            private:
                std::vector<size_t> m_nodes;
                const CNetwork* m_pG;
                std::vector<size_t>::iterator m_Iter;
            public:
                vertex_iterator(const CNetwork* pG, size_t src_node, AlgorithmType type = BreadthFirstSearch) : m_pG(pG)
                {
                    if(type == BreadthFirstSearch)
                    {
                        // set nodes to BFS
                        std::queue<size_t> processList;     // stack of nodes to be processed.
                        size_t currentNode;
                        std::vector<bool> found(m_pG->NumNodes(), false);
                        processList.push(src_node);
                        found[src_node] = true;
                        while(!processList.empty())
                        {
                            currentNode = processList.front();
                            m_nodes.push_back(currentNode);
                            processList.pop();

                            const std::vector<size_t>& neighbors = m_pG->GetNeighbors(currentNode);

                            for(std::vector<size_t>::const_iterator n = neighbors.cbegin(); n != neighbors.cend(); n++)
                            {
                                if( !found[*n] ) // distance is not set. means this is a discovered node.
                                {
                                    found[*n] = true;
                                    processList.push(*n);
                                }
                            }
                        }

                    }
                    else if (type == DepthFirstSearch)
                    {
                        // set nodes to DFS
                        std::stack<size_t> processList;     // stack of nodes to be processed.
                        size_t currentNode;
                        std::vector<bool> found(m_pG->NumNodes(), false);
                        processList.push(src_node);
                        found[src_node] = true;
                        while(!processList.empty())
                        {
                            currentNode = processList.top();
                            m_nodes.push_back(currentNode);
                            processList.pop();
                            const std::vector<size_t>& neighbors = m_pG->GetNeighbors(currentNode);
                            for(std::vector<size_t>::const_iterator n = neighbors.cbegin(); n != neighbors.cend(); n++)
                            {
                                if( !found[*n] ) // distance is not set. means this is a discovered node.
                                {
                                    found[*n] = true;
                                    processList.push(*n);

                                }
                            }
                        }
                    }
                    else
                    {
                        // Error condition
                    }
                    m_Iter = m_nodes.begin();
                }

                vertex_iterator& operator ++() { m_Iter++; return *this; }

                const size_t& operator *() { return *m_Iter; }

                bool IsValid() { return m_Iter != m_nodes.end(); }

                const std::vector<size_t>& SortedVertices() { return m_nodes; }
        };

        static CNetwork Cycle(unsigned int n, size_t netType = Undirected)
        {
            CNetwork cycle(n, netType);
            for(size_t i = 0; i < n; i++) cycle.AddEdge( (i+1)%n, i ); // i -> i+1
            return cycle;
        }

    public:
        CNetwork() : /*CJSONValueObject<CNetwork>("", this),*/ m_bInitialized(false), m_NetworkMask(General), m_AdjList(m_AdjMat, false)//, m_json_adj_mat(JSON_NAME_ADJMAT, &m_AdjMat)
        { /*SetupJSONObject();*/ }

        CNetwork(int nodes, size_t netType =  General) : /*CJSONValueObject<CNetwork>("", this),*/ m_bInitialized(false), m_NetworkMask(netType), m_AdjList(m_AdjMat, IsDirected())//, m_json_adj_mat(JSON_NAME_ADJMAT, &m_AdjMat)
        {
            Initialize(nodes);
            /*SetupJSONObject();*/
        }

        CNetwork(const CNetwork& src) : /*CJSONValueObject<CNetwork>("", this),*/ m_bInitialized(false), m_AdjList(m_AdjMat, src.IsDirected())//, m_json_adj_mat(JSON_NAME_ADJMAT, &m_AdjMat)
        {
            CopyFrom(src);
            /*SetupJSONObject();*/
            MakeAdjacencyList();
        }

        CNetwork(const EigenMatrixHelper::spMatrix& adj, size_t mask) : /*CJSONValueObject<CNetwork>("", this),*/ m_bInitialized(true), m_NetworkMask(mask), m_AdjMat(adj), m_AdjList(m_AdjMat, IsDirected())//, m_json_adj_mat(JSON_NAME_ADJMAT, &m_AdjMat)
        {
            /*SetupJSONObject();*/
            MakeAdjacencyList();
        }

    #ifdef c_plus_plus_11
        template<template<typename ... > class Container >
        CNetwork(const Container<Network::CNetwork::edge_info>& edges, int nodes, size_t netType) : /*CJSONValueObject<CNetwork>("", this),*/ m_bInitialized(false), m_NetworkMask(netType), m_AdjList(m_AdjMat, IsDirected())//, m_json_adj_mat(JSON_NAME_ADJMAT, &m_AdjMat)
        {
            SetFromEdgeInfo(edges, nodes, netType);
            /*SetupJSONObject();*/
        }

        template<template<typename ... > class Container >
        void SetFromEdgeInfo(const Container<Network::CNetwork::edge_info>& edges, int nodes, size_t netType)
        {
            m_NetworkMask = netType;
            Initialize(nodes);
            for(auto edge : edges) AddEdge(edge.target, edge.source, edge.weight);
            MakeAdjacencyList();
        }

    #endif

        ~CNetwork() { Destroy(); }

        void Destroy()
        {
//            std::cout << "Calling CNetwork::Destroy() "<< this << std::endl;
            // CJSONValueObject<CNetwork>::Destroy();
            // m_json_adj_mat.Destroy();
        }

        void Initialize(int nodes)
        {
            m_AdjMat.resize(nodes, nodes); // this will zero out the matrix.
            m_AdjList.MarkToUpdate();
            m_bInitialized = true;
        }

    // Network IO using JSON reader/writer.
    /*
        virtual void SetupJSONObject()
        {
            // std::cout << "Calling CNetwork::SetupJSONObject" << std::endl;
            m_json_adj_mat.SetupJSONObject();
            AddObjectValue(JSON_NAME_ADJMAT, &m_json_adj_mat);
            AddUIntegerValue(JSON_NAME_NETWORK_TYPE, &m_NetworkMask);
        }

        bool LoadFromFile(const string& Path)
        {
            json::CJSONParser json;
            json.LoadFromFile(Path);    // opens the JSON file and loads the data into buffer.
            m_bInitialized = json.ParseObject(this); // Parses all the data for this object.

            if( m_bInitialized )
                MakeAdjacencyList();

            return m_bInitialized;
        }

        bool DumpToFile(const string& Path)
        {
            bool bDumpSuccess = false;
            json::CJSONParser json;
            bDumpSuccess = json.DumpObjectToFile(Path, this);

            if(!bDumpSuccess)
                std::cout << "Error!! Could not dump file: " << Path << std::endl;

            return bDumpSuccess;

        }

        bool Parse(const json_t *pVal)
        {
            bool bParseSuccess = CJSONValueObject<CNetwork>::Parse(pVal);
            m_bInitialized = bParseSuccess && m_AdjMat.rows();
            // std::cout << "Parsed ? " << boolalpha << m_bInitialized << ", " << bParseSuccess << std::endl;
            if( m_bInitialized )
                MakeAdjacencyList();


            return bParseSuccess;
        }
*/
    // General Network Methods

        bool MakeAdjacencyList() { return m_AdjList.Update(); }

        void AddEdge(size_t nodei, size_t nodej, double weight = 1.0)
        {
            // AddEdge will add another edge between nodes i and j.
            // cases:
            //      -- Undirected:  A(i,j) and A(j,i) is incremented by weight.
            //                      i == j case and weight != 1.0 then A(i,i) is incremented only once.
            //      -- Directed:    A(i,j) is incremented by weight.
            //                      note that i <- j
            // std::cout << nodei << "<-" << nodej << " num_nodes = " << NumNodes() << std::endl;

            assert(m_bInitialized && nodei < NumNodes() && nodej < NumNodes());

            AssertMaskIsValid();
            int i = int(nodei);
            int j = int(nodej);

            m_AdjList.UpdateCheck(nodej, nodei);

            if(!IsWeighted())
                weight = 1.0;

            if(IsUndirected())
            {
                if(IsSimple())
                {
                    if(i != j) // no self-loops in simple graph
                    {
                        m_AdjMat.coeffRef(i,j) = weight;
                        m_AdjMat.coeffRef(j,i) = weight;
                    }
                }
                else
                {
                    // Need to look at Newman's book to make sure this is correct.
                    if(i == j && weight != 1.0)
                    {
                        m_AdjMat.coeffRef(i,j) += weight;
                    }
                    else
                    {
                        m_AdjMat.coeffRef(i,j) += weight;
                        m_AdjMat.coeffRef(j,i) += weight;
                    }
                }
            }
            else if(IsDirected())
            {
                if(IsSimple() && i != j){
                    m_AdjMat.coeffRef(i,j) = weight;
                }
                else{
                    m_AdjMat.coeffRef(i,j) += weight;
                }
            }
            else
            {
                throw ErrorGeneral;
            }

            MakeAdjacencyList();
        }

        void InsertEdge(const size_t& src, const size_t& tgt, double weight = 1.0)
        {
            assert(m_bInitialized && src < NumNodes() && tgt < NumNodes());

            AssertMaskIsValid();
            int t = int(tgt);
            int s = int(src);
            m_AdjList.MarkToUpdate();   // mark the adj for update.
            if(!IsWeighted())
                weight = 1.0;

            if(IsUndirected())
            {
                if(IsSimple())
                {
                    if(s != t) // no self-loops in simple graph
                    {
                        m_AdjMat.insert(t, s) = weight;
                        // m_AdjMat.coeffRef(i,s) = weight;
                        // m_AdjMat.coeffRef(j,i) = weight;
                    }
                }
                else
                {
                    // Need to look at Newman's book to make sure this is correct.
                    if(s == t && weight != 1.0)
                    {
                        // m_AdjMat.coeffRef(i,j) += weight;
                        if(!IsWeighted())
                            weight = 2.0;
                        m_AdjMat.insert(t, s) = weight;
                    }
                    else
                    {
                        m_AdjMat.insert(t, s) = weight;
                        // m_AdjMat.coeffRef(i,j) += weight;
                        // m_AdjMat.coeffRef(j,i) += weight;
                    }
                }
            }
            else if(IsDirected())
            {
                if(IsSimple() && s != t){
                    m_AdjMat.insert(t, s) = weight;
                    // m_AdjMat.coeffRef(i,j) = weight;
                }
                else{
                    m_AdjMat.insert(t, s) = weight;
                    // m_AdjMat.coeffRef(i,j) += weight;
                }
            }
            else
            {
                throw ErrorGeneral;
            }
        }

        template< class InputIterator, class UnaryOperator >
        void AddPath(InputIterator first, InputIterator last, UnaryOperator func)
        {
            while(first != last)
            {
                size_t src = func(*first++);
                if(first != last)
                    AddEdge(func(*first), src);
            }

            MakeAdjacencyList();
        }


        void RemoveEdge(size_t nodei, size_t nodej, double weight = 1.0)
        {
            AssertMaskIsValid();
            int i = int(nodei);
            int j = int(nodej);

            weight = (IsMultiGraph() && (weight < m_AdjMat.coeffRef(i,j))) ? m_AdjMat.coeffRef(i,j) - weight : 0.0; // Allow weighted graphs to do this too?

            if(IsUndirected())
            {
                m_AdjMat.coeffRef(i,j) = weight;
                m_AdjMat.coeffRef(j,i) = weight;

            }
            else if(IsDirected())
            {
                m_AdjMat.coeffRef(i,j) = weight;

            }
            else
            {
                throw ErrorGeneral;
            }

            if(fabs(weight)  < SMALL){
                m_AdjList.MarkToUpdate();   // mark the adj for update.
                PruneNetwork(SMALL);        // remove it from the sparse matrix lists.
                m_AdjList.Update();
            }
        }

        bool DoesEdgeExist(size_t src, size_t tgt) const
        {
//            // checks if there is a link between nodei -> nodej
//            bool bExists = false;
//
//            for (int k=0; k<m_AdjMat.outerSize() && !bExists; ++k)
//            {
//                for (Eigen::SparseMatrix<double>::InnerIterator it(m_AdjMat,k); it && !bExists; ++it)
//                {
//                    if((size_t(it.col()) == src && size_t(it.row()) == tgt) && fabs(it.value()) > 0)
//                    {
//                        bExists = true;
//                    }
//                }
//            }
//
            return m_AdjList.DoesEdgeExist(src, tgt);
        }

        std::vector< std::vector<size_t> > GetEdgePairList() const
        {
            std::vector< std::vector<size_t> > edges;

            for (int k=0; k<m_AdjMat.outerSize(); ++k)
            {
                for (Eigen::SparseMatrix<double>::InnerIterator it(m_AdjMat,k); it; ++it)
                {
                    if(IsUndirected() && it.row() <= it.col())
                    {
                        std::vector<size_t> edge;
                        edge.push_back(size_t(it.row())); edge.push_back(size_t(it.col()));
                        edges.push_back(edge);
                    }
                    else if(IsDirected())
                    {
                        std::vector<size_t> edge;
                        edge.push_back(size_t(it.row())); edge.push_back(size_t(it.col()));
                        edges.push_back(edge);
                    }
                }
            }
            return edges;
        }

        // Algorithms to implement:
        // 1. Breadth first search
        // 2. Dijkstra's Algorithm
        //  (shortest / longest path algorithms)
        // 3. Vietbo's Algorithm -- probably not what we want.
        // 4. Min Cut, Max Flow
        // 5. Isomorphism finding Algorithms
        // 6. Degree distributions algorithms
        // 6.5 Centrality measures
        // 7. Similarity measures
        // 8. proper cycle.

        void FindCycles()
        {

        /*
        *   pseudo code:
        *
        *
            for each node, n
                for each neighbor of n, ni
                    for each neighbor of n, nj
                        compute shortest path from ni to nj in G - {n}
                        put in library of paths based using minimal rotation to make this n*log(n)


        *   Complexity: O(n*<k^2>*(c*log(c))*O(shortest path algorithm))
        *
        */

        }

        template<class ArrayType>
        void BreadthFirstSearch(const size_t& source, const size_t& target, CNetwork* pPathTree, ArrayType& Distances) const
        {
            // This algorithm is implemented using the discription on p.320 from Mark Newman's book: Networks, An Indtroduction.
            // to find all the distances then set target = source.
            // TODO: figure out how to deal with the target.

            queue<size_t> processList;      // queue of nodes to be processed.
            size_t dist = 0;                // holds the current distance from source
            bool bTargetFound = false;
            size_t currentNode;

            // Step 1:  Initialize data structures.
            processList.push(source);
            for(size_t n = 0; n < NumNodes(); n++)
                Distances[n] = INVALID_DISTANCE;
            Distances[source] = 0;


            while(!processList.empty())
            {
                // Step 2:  Read next node from queue.
                currentNode = processList.front();
                processList.pop();

                // Step 3:  Find the distance.
                dist = Distances[currentNode];

                // Step 4:  Go to each neighbor, update the distance if it is not set
                //          Add discovered node to queue
                //          Update PathTree with a directed edge currentNode <- Discovered Node
                const std::vector<size_t>& neighbors = GetNeighbors(currentNode);

                for(std::vector<size_t>::const_iterator n = neighbors.cbegin(); n != neighbors.cend(); n++)
                {
                    if(Distances[*n] == INVALID_DISTANCE) // distance is not set. means this is a discovered node.
                    {
                        processList.push(*n);
                        Distances[*n] = dist + 1;
                        if(pPathTree)
                            pPathTree->AddEdge(currentNode, *n);
                    }
                    else if(Distances[*n] == (dist + 1)) // then there exists another path from n->s of length d+1 through currentNode.
                    {
                        if(pPathTree)
                            pPathTree->AddEdge(currentNode, *n);
                    }

                    if(*n == target && source != target) // found target node, we can stop early.
                    {
                        bTargetFound = true; // TODO: use this to stop early.
                    }
                }
                // step 5:  Repeat from Step 2.
            }
        }

        void DepthFirstSearch( const size_t& source ) const
        {

            std::stack<size_t> processList;     // stack of nodes to be processed.
            size_t currentNode;
            std::vector<bool> found(NumNodes(), false);

            // Step 1:  Initialize data structures.
            processList.push(source);
            found[source] = true;
            while(!processList.empty())
            {
                // Step 2:  Read next node from queue.
                currentNode = processList.top();
                processList.pop();

                const vector<size_t>& neighbors = GetNeighbors(currentNode);

                for(vector<size_t>::const_iterator n = neighbors.cbegin(); n != neighbors.cend(); n++)
                {
                    if( !found[*n] ) // distance is not set. means this is a discovered node.
                    {
                        found[*n] = true;
                        processList.push(*n);
                    }
                }
            }
        }

        vector< vector<size_t> > ShortestPath(const size_t& source, const size_t& target) const
        {
            CNetwork shortestPaths(int(NumNodes()), (Directed));
            // std::unique_ptr<size_t []> pDistances(new size_t[NumNodes()]);
            size_t* pDistances = new size_t[NumNodes()];
            try{
            BreadthFirstSearch(source, target, &shortestPaths, pDistances); // will just find one of the shortest paths.
            }
            catch(...)
            {
                std::cout << "Error occured in Shortest path. " << std::endl;
            }
            delete [] pDistances;

            shortestPaths.MakeAdjacencyList();

            return shortestPaths.FindPaths(source, target);
        }

        CNetwork QuotientGraph(const std::vector<size_t>& Partitions, std::map<size_t, size_t>& PartitionMap, std::map<size_t, size_t>& VertexMap) const
        {
            if(Partitions.size() != NumNodes()) // all nodes must be specified
                return CNetwork(0,0);

            int numPartitions = 0;
            for(size_t i = 0; i < Partitions.size(); i++)
            {
                std::map<size_t, size_t>::iterator f = PartitionMap.find(Partitions[i]);
                if(f == PartitionMap.end())
                    PartitionMap.insert(pair<size_t, size_t>(Partitions[i], numPartitions++));
            }
            for(size_t i = 0; i < Partitions.size(); i++) VertexMap.insert(pair<size_t, size_t>(i, PartitionMap[Partitions[i]]));


            CNetwork qgraph(numPartitions, m_NetworkMask); // same type of network.

            for(edge_iterator e(this); e.IsValid(); e++) // now set all of the edges.
            {
                size_t src = PartitionMap[Partitions[size_t(e->source)]], tgt = PartitionMap[Partitions[size_t(e->target)]];
                // std::cout << "found edge " << src << " -> " << tgt << std::endl;
                if(src != tgt) // no self-loops by definition
                    qgraph.AddEdge(tgt, src); // this will sum the weights if not a simple graph.
            }
            return qgraph;
        }

        template<class ForwardIterator>
        void SubGraph(ForwardIterator first_node, ForwardIterator last_node, CNetwork& ret)
        {
            // TODO: this seems round about so maybe we can make it better somehow.

            assert(std::is_sorted(first_node, last_node));
            size_t V = std::distance(first_node, last_node); // number of vertices.
            if(V == NumNodes())
            {
                ret = *this;
                return;
            }

            Eigen::SparseMatrix<double, Eigen::RowMajor> subRowAdjMat(m_AdjMat);
            int i = 0, removed = 0;
            ForwardIterator it = first_node;
            while(it != last_node)
            {
                if(i < *it) {
                    EigenMatrixHelper::RemoveRow(subRowAdjMat, (i-removed)); i++; removed++;
                }
                else if(*it < i) {
                    it++;
                }
                else {
                    it++; i++;
                }
            }
            while( size_t(i) < NumNodes() ) { EigenMatrixHelper::RemoveRow(subRowAdjMat, (i-removed)); i++; removed++; }
            assert( subRowAdjMat.rows() == V );

            // Now to remove the columns
            EigenMatrixHelper::spMatrix subAdjMat(subRowAdjMat);
            i = 0, removed = 0;
            it = first_node;
            while(it != last_node)
            {
                if(i < *it) {
                    EigenMatrixHelper::RemoveColumn(subAdjMat, (i-removed)); i++; removed++;
                }
                else if(*it < i) {
                    it++;
                }
                else {
                    it++; i++;
                }
            }
            while( size_t(i) < NumNodes() ) { EigenMatrixHelper::RemoveColumn(subAdjMat, (i-removed)); i++; removed++; }
            assert( subAdjMat.cols() == V );

            ret.Initialize(subAdjMat.rows());
            ret.SetAdjMatrix(subAdjMat);
            ret.MakeAdjacencyList();
        }
        template<class ForwardIterator>
        void SubGraphEx(ForwardIterator first_node, ForwardIterator last_node, CNetwork& ret)
        {
            /*
            From the Eigen Doc the fastest way to fill the sparse matrix is as follows
            1: SparseMatrix<double> mat(rows,cols);         // default is column major
            2: mat.reserve(VectorXi::Constant(cols,6));
            3: for each i,j such that v_ij != 0
            4:   mat.insert(i,j) = v_ij;                    // alternative: mat.coeffRef(i,j) += v_ij;
            5: mat.makeCompressed();                        // optional
            */
            size_t V = std::distance(first_node, last_node); // number of vertices.
            if(V == NumNodes())
            {
                ret = *this;
                return;
            }
            EigenMatrixHelper::spMatrix mat(V,V);
            mat.reserve(Eigen::VectorXi::Constant(V, V < 100 ? V : 100));
            ForwardIterator it = first_node;
            int s = 0;
            while(it != last_node)
            {
                const std::vector<size_t>& n =  m_AdjList.Successors(*it);
                const std::vector<double>& w =  m_AdjList.Weights(*it);
                for(size_t i = 0; i < n.size(); i++){
                    ForwardIterator find = std::lower_bound(first_node,last_node,n[i]);
                    if(find!=last_node && !(n[i]<*find)){
                        int t = std::distance(first_node, find);
                        mat.insert(t, s) = w[i];
                    }
                }
                it++;
                s++;
            }
            mat.makeCompressed();
            ret.Initialize(mat.rows());
            ret.SetAdjMatrix(mat);
            ret.MakeAdjacencyList();
        }

        const vector<size_t>& GetNeighbors(const size_t& n, bool bPredecessor = false) const
        {
            // If bPredecessor is true then get the in-neighbors (predecessors), otherwise get the out-neighbors (successors).
            // for undirected graphs the Predecessors == Successors and bPredecessor is ignored.
            if(bPredecessor && IsDirected())
            {
                return m_AdjList.Predecessors(n);
            }
            else
            {
                return m_AdjList.Successors(n);
            }
        }

        bool KruskalAlgorithm(CNetwork& SpanningForest); // inline impl below.

        bool IsGraphConnected()
        {
            vector<size_t> distances(NumNodes(), size_t(INVALID_DISTANCE));
            BreadthFirstSearch(0, 0, NULL, distances);
            return !utils::IsInVec(distances, size_t(INVALID_DISTANCE));
        }

        template<class Container>
        bool IsConnected(const Container& src, const Container& tgt )
        {
            MakeAdjacencyList();
            bool bConnected = false;
            for(size_t s = 0; s < src.size() && !bConnected; s++)
            {
                for(size_t t = 0; t < tgt.size() && !bConnected; t++)
                {
                    bConnected = AreNodesConnected(src[s], tgt[t]);
                }
            }
            return bConnected;
        }

        bool AreNodesConnected(const size_t& source, const size_t& target) const
        {
            vector<size_t> distances;
            distances.resize(NumNodes(), size_t(INVALID_DISTANCE));
            BreadthFirstSearch(source, target, NULL, distances);
            return distances[target] != size_t(INVALID_DISTANCE);
        }


        template < class OutputIterator >
        void ComputeComponents(OutputIterator ret)
        {
            vector<bool> found(NumNodes(), false);
            size_t components = 0;
            for(size_t i = 0; i < NumNodes(); i++)
            {
                if(found[i])
                    continue;

                found[i] = true;
                vector<size_t> distances(NumNodes(), size_t(INVALID_DISTANCE));
                BreadthFirstSearch(i, i, NULL, distances);

                for(size_t j = 0; j < distances.size(); j++)
                {
                    if(distances[j] != size_t(INVALID_DISTANCE))
                    {
                        found[j] = true;
                        OutputIterator iter = ret;
                        std::advance(iter, j);
                        *iter = components;
                    }
                }
                components++;
            }
        }

    public:
        void FindPathsR(const size_t& node, const size_t& target, vector< vector<size_t> >& paths, vector<size_t> path = vector<size_t>()) const
        {
            path.push_back(node);
            if(node == target) // base case.
            {
                paths.push_back(path);
            }
            else
            {
            // depth first search
            const std::vector<size_t>& n = GetNeighbors(node);
            for(size_t i = 0; i < n.size(); i++)
                {
                FindPathsR(n[i], target, paths, path); // recursive call
                }
            }
        }
        // probably want to change this or change the name of this function.
        // good that it is private.
        vector< vector<size_t> > FindPaths(const size_t& source, const size_t& target) const
        {
            if (!IsDirected()) // expecting the DAG from BFS.
            {
                std::cout << "FindPaths only works for directed acyclic graphs." << std::endl;
                return vector< vector<size_t> >();
            }
            else if(source == target) // need to verify that source is the "leaf" node.
            {
                std::cout << "**WARNING** source node is equal to target node. " << std::endl;
                return vector< vector<size_t> >();
            }

            vector< vector<size_t> > paths;
            std::vector<bool> pFound(NumNodes(), false);

            queue<size_t> processList;
            size_t node = 0;
            size_t ct =0;

//            for (size_t b = 0; b < NumNodes(); b++)
//                pFound[b] = false;
            pFound[target] = true;

            vector<size_t> neighbors = GetNeighbors(target);  // Gets the upstream neighbors.
            for (vector<size_t>::const_iterator n = neighbors.cbegin(); n != neighbors.cend(); n++)
            {
                vector<size_t> temp;
                if(AreNodesConnected(*n, source)) // They must be connected to be added otherwise incomplete paths will be added.
                {
                    temp.push_back(target);
                    temp.push_back(*n);
                    paths.push_back(temp);
                    processList.push(*n);
                }
                pFound[*n] = true;
            }

            while (!processList.empty())
            {
                node = processList.front();
                // std::cout << "Items left to process "<< processList.size() << " num paths = "<< paths.size() << std::endl;
                processList.pop();

                neighbors = GetNeighbors(node);  // Gets the down stream neighbors.
                ct = 0;
                for (vector<size_t>::iterator n = neighbors.begin(); n != neighbors.end(); n++)
                {
                    if(!AreNodesConnected(*n, source)) // They must be connected to be added otherwise incomplete paths will be added.
                        continue;
                    // std::cout << "ct =  "<< ct << " of " <<  neighbors.size() << " num paths = "<< paths.size() << std::endl;
                    for (size_t p = 0; p < paths.size(); p++)
                    {
                        // std::cout << "paths["<< p << "].size() = " << paths[p].size() << std::endl;
                        size_t len = paths[p].size()-1;
                        if (paths[p][len] == node)
                        {
                            if (ct == (neighbors.size()-1))
                            {
                                // just append to the path already in the list.
                                paths[p].push_back(*n);
                            }
                            else if (paths.size() < 1000)
                            {
                                // copy the path and append to the new path.
                                vector<size_t> copyPath = paths[p];
                                copyPath.push_back(*n);
                                paths.push_back(copyPath);
                            }
                        }
                    }
                    processList.push(*n); // Always add so there will be an infinite loop if the Graph is not DAG
                    ct++;
                }
            }
            return paths;
        }

    public:
    // Centralities
        double ClosenessCentrality(size_t node)
        {


#ifdef c_plus_plus_11
            std::unique_ptr<size_t []> pDistances(new size_t[NumNodes()]);
            double sum = 0.0;
            size_t ct = 0;
            BreadthFirstSearch(node, node, NULL, pDistances); // will just find one of the shortest paths.

            for(size_t n = 0; n < NumNodes(); n++)
            {
                if(pDistances[n] != size_t(INVALID_DISTANCE))
                {
                    ct++;
                    sum += pDistances[n];
                }
            }

            std::cout << "Found " << ct << " vertices in the component of vertex " << node << "." << std::endl;
            return double(ct)/sum;
#else
            std::cout << "Warning*** using a function that requires c++11 in c98!!" << std::endl;
            return 0;
#endif

        }


        // very inefficient way of calculating this.
        size_t GetDegree(const size_t& ndx, bool bOut = true)
        {
            EigenMatrixHelper::spMatrix degreeDist;
            EigenMatrixHelper::spMatrix ones = (Eigen::VectorXd::Ones(int(NumNodes()))).sparseView();

            if(IsDirected() && bOut)
            {
                degreeDist = m_AdjMat.transpose()*ones;
            }
            else
            {
                degreeDist = m_AdjMat*ones;
            }
            return size_t(degreeDist.coeffRef(int(ndx),0));
        }

    // Mask Methods.
        void AssertMaskIsValid() const
        {
            bool expr = IsTypeSet();
            expr = expr && (IsDirected() != IsUndirected());
            assert(expr);
        }

        bool IsTypeSet() const    { return m_NetworkMask !=  General; }

        bool IsDirected() const   { return (m_NetworkMask &  Directed) ==  Directed; }

        bool IsUndirected() const { return (m_NetworkMask &  Undirected) ==  Undirected; }

        bool IsSimple() const     { return (m_NetworkMask &  Simple) == Simple; }

        bool IsWeighted() const   { return (m_NetworkMask & Weighted) == Weighted; }

        bool IsMultiGraph() const { return (m_NetworkMask & Multigraph) == Multigraph; }

        void AddNetworkType(int type) { m_NetworkMask |= size_t(type); }

    // Accessor methods
        size_t      GetNetworkMask()    const   { return m_NetworkMask; }

        bool        IsInitialized()     const   { return m_bInitialized; }

        size_t      NumNodes() const { return size_t(m_AdjMat.rows()); }
        size_t      NumEdges() const
        {
            AssertMaskIsValid();
            if(IsDirected())
            {
                return size_t(m_AdjMat.nonZeros());
            }
            else if (IsUndirected())
            {
                if(IsSimple())
                {
                    if(m_AdjMat.nonZeros()%2 != 0)
                        std::cout << m_AdjMat << std::endl;
                    assert( m_AdjMat.nonZeros()%2 == 0);
                    return size_t(m_AdjMat.nonZeros() / 2); // none on the diag and the matrix is symmetric.
                }
                else{
                    int self_loops = m_AdjMat.diagonal().nonZeros();
                    int total_edges = int((m_AdjMat.nonZeros()-self_loops)/2) + self_loops;
                    return size_t(total_edges);
                }
            }
            assert(false); // should never make it here.
            return 0;
        }
        const AdjacencyList& GetAdjList() const { return m_AdjList; }
        const EigenMatrixHelper::spMatrix& GetAdjMatrix() const { return m_AdjMat; }

        EigenMatrixHelper::spMatrix& GetAdjMatrix() { return m_AdjMat; }

        void SetAdjMatrix(const EigenMatrixHelper::spMatrix& adj) { m_AdjMat = adj; }

    // Sparse Matrix methods
        void PruneNetwork(const double& PruneVal = 0.0 ) { m_AdjMat.prune(PruneVal); }

        template< class ContainerType >
        void Reserve(const ContainerType& sizes) { m_AdjMat.reserve(sizes); }

    // Printing Functions.
        void PrintNetwork() const
        {
            //PruneNetwork();
            printf("Network Mask = %zu \nNumber of Nodes = %zu \n", m_NetworkMask, NumNodes());
            for (int k=0; k<m_AdjMat.outerSize(); ++k)
            {
                for (Eigen::SparseMatrix<double>::InnerIterator it(m_AdjMat,k); it; ++it)
                {
                    printf("(%li <- %li) weight = %f \n", it.row(), it.col(), it.value());
                }
            }
            if(m_AdjMat.nonZeros() == 0)
            {
                printf("There are no links to report. \n");
            }
            m_AdjList.PrintListStats();
            std::cout << "network adj matrix addr : " << std::hex << &m_AdjMat << std::dec << std::endl;
            std::cout << "num rows : " << m_AdjMat.rows() << std::endl;
            std::cout << "num cols : " << m_AdjMat.cols() << std::endl;
            // std::cout << m_AdjMat<<endl;
        }

        bool PrintToFile(FILE* pFile)
        {
            bool bSuccess = pFile!=NULL;
            if(pFile){
                for (int k=0; k<m_AdjMat.outerSize(); ++k)
                {
                    for (Eigen::SparseMatrix<double>::InnerIterator it(m_AdjMat,k); it; ++it)
                    {
                        fprintf(pFile, "%li -> %li, %f \n", it.col(), it.row(), it.value());
                    }
                }
            }
            return bSuccess;
        }

    private: // should be protected?
     // Helper Functions
        bool IsEqual( const CNetwork& cmp)
        {
            bool bEqual = (m_AdjMat.rows() == cmp.GetAdjMatrix().rows()) && (m_AdjMat.cols() == cmp.GetAdjMatrix().cols()); // dim is equal.
            if(bEqual)
            {
                EigenMatrixHelper::spMatrix result = m_AdjMat - cmp.GetAdjMatrix();
                result.prune(SMALL);
                bEqual = result.nonZeros()==0;
                //cout << "There were "<< result.nonZeros() << " elements that did not match." << std::endl;
            }
            return bEqual;
        }

        CNetwork& CopyFrom(const CNetwork& src)
        {
            m_NetworkMask = src.GetNetworkMask();
            m_AdjMat = src.GetAdjMatrix();
            m_AdjList.MarkToUpdate();
            MakeAdjacencyList();
            m_bInitialized = src.IsInitialized();
            return *this;
        }

    public:
    // Overloaded Operators
        bool operator == ( const CNetwork& cmp)
        {
            return IsEqual(cmp);
        }

        CNetwork& operator = (const CNetwork& src) { return CopyFrom(src); }

        double operator()(const size_t& src, const size_t& tgt) const { return m_AdjMat.coeff(int(tgt), int(src)); }

    protected:
        bool                            m_bInitialized;
        size_t                          m_NetworkMask;
        EigenMatrixHelper::spMatrix     m_AdjMat;
        AdjacencyList                   m_AdjList;

    // JSON parser helper object
        // EigenMatrixHelper::CJSONValueSparseMatrix m_json_adj_mat;


};

class DefaultSemanticFeasibility : public std::binary_function<size_t, size_t, bool>
{
    public:
        DefaultSemanticFeasibility() {}
        virtual bool operator()(const size_t&, const size_t&) { return true; }
};

template< class Compare >
class EdgeFeasibilityFunction : public std::binary_function<size_t, size_t, bool>
{
    Compare m_Compare;
    public:
        EdgeFeasibilityFunction(const Compare& c = Compare()) : m_Compare(c) {}
        virtual ~EdgeFeasibilityFunction() {}
        virtual bool operator()(const CNetwork::edge_info& e1, const CNetwork::edge_info& e2) { return m_Compare(e1.weight, e2.weight); }
};

// Some special cases.
typedef EdgeFeasibilityFunction< std::equal_to<double> > DefaultEdgeFeasibilityFunction;

#ifdef c_plus_plus_11
template<int power = -6, size_t constant = 1>
using ApproxEdgeFeasibilityFunction =  EdgeFeasibilityFunction< utils::float_is_equal<double, power, constant> >;
#else
template<int power = -6, size_t constant = 1>
struct ApproxEdgeFeasibilityFunction{ typedef EdgeFeasibilityFunction< utils::float_is_equal<double, power, constant> > type; };
#endif
typedef EdgeFeasibilityFunction< utils::TrueFunction > IgnoreEdgeFeasibilityFunction;

template< typename SemanticFeasibility = DefaultSemanticFeasibility, typename EdgeFeasibility = IgnoreEdgeFeasibilityFunction >
class CIsomorphism
{
    class vf2_state // should I take this out of this class?
    {

        public:
            vf2_state(const CNetwork* pg1, const CNetwork* pg2)
            {
                n1 = pg1->NumNodes();
                n2 = pg2->NumNodes();
                depth = 0;

                core1.resize(pg1->NumNodes(), vf2_state::NULL_NODE());
                core2.resize(pg2->NumNodes(), vf2_state::NULL_NODE());

                in1.resize(pg1->NumNodes(), 0);
                out1.resize(pg1->NumNodes(), 0);

                in2.resize(pg2->NumNodes(), 0);
                out2.resize(pg2->NumNodes(), 0);
            }
            vf2_state(const vf2_state& src)
            {
                core1 = src.core1;
                core2= src.core2;
                in1 = src.in1;
                in2 = src.in2;
                out1 = src.out1;
                out2 = src.out2;

                n1 = src.n1;
                n2 = src.n2;

                depth = src.depth;

            }
            void restore() {} // not sure if we need to do this yet...
            bool ShouldTerminate()
            {
                size_t n = 0;
                for(size_t i = 0; i < core1.size(); i++)
                    if(core1[i] != vf2_state::NULL_NODE()) // assumes that each entry is unique.
                        n++;
                return n == n2;
            }

            void inclusion(const std::pair<size_t, size_t>& nodes, const CNetwork* pg1, const CNetwork* pg2)
            {
                // add in an assert so that we can make sure all of the node mappings will be unique (at least in debug mode)
                assert(core1[nodes.first] == vf2_state::NULL_NODE() && core2[nodes.second] == vf2_state::NULL_NODE());

                ++depth; // keep track of the current depth.

                // now add the next pair to the state.
                core1[nodes.first] = nodes.second;
                core2[nodes.second] = nodes.first;

                if(in1[nodes.first] == 0)
                    in1[nodes.first] = depth;
                if(out1[nodes.first] == 0)
                    out1[nodes.first] = depth;

                if(in2[nodes.second] == 0)
                    in2[nodes.second] = depth;
                if(out2[nodes.second] == 0)
                    out2[nodes.second] = depth;

                // Also update the in/out vectors for this state.
                CalculateT(pg1, n1, core1, in1, out1);
                CalculateT(pg2, n2, core2, in2, out2);

            }

            bool IsInTin1(const size_t& node1) { return in1[node1] == depth && core1[node1] == NULL_NODE(); }
            bool IsInTin2(const size_t& node2) { return in2[node2] == depth && core2[node2] == NULL_NODE(); }

            bool IsInTout1(const size_t& node1) { return out1[node1] == depth && core1[node1] == NULL_NODE(); }
            bool IsInTout2(const size_t& node2) { return out2[node2] == depth && core2[node2] == NULL_NODE(); }

            bool IsInNc1(const size_t& node1) { return out1[node1] == 0; }
            bool IsInNc2(const size_t& node2) { return out2[node2] == 0; }

            static const size_t NULL_NODE() { return size_t(-1); }


            std::vector<size_t> core1;        // holds the current mapping. g1 -> g2
            std::vector<size_t> core2;        // holds the current mapping. g2 -> g1
            std::vector<size_t> in1;          // the ith element is nonzero if m_core1[i] != NULL_NODE or i is in if T_1^{in}.
            std::vector<size_t> in2;          //     -- holds the value of the depth that the node was added to the mapping.
            std::vector<size_t> out1;
            std::vector<size_t> out2;

        private:
            size_t n1;
            size_t n2;
            size_t depth;

            void CalculateT(const CNetwork* pg, const size_t& n, std::vector<size_t>& core0, std::vector<size_t>& in0, std::vector<size_t>& out0)
            {
                // update the other vectors here too.
                for(size_t v = 0; v < n; v++)
                {
                    if(core0[v] == vf2_state::NULL_NODE())      // v is not in M(s)
                    {
                        continue;
                    }
                    else                                        // v is in M(s)
                    {
                        const vector<size_t>& predecessors = pg->GetNeighbors(v, true); // could optimize based on graph type.
                        const vector<size_t>& successors = pg->GetNeighbors(v, false);
                        for(size_t p = 0; p < predecessors.size(); p++)
                        {
                            if(in0[predecessors[p]] == 0 && core0[predecessors[p]] == vf2_state::NULL_NODE())
                            {
                                in0[predecessors[p]] = depth;
                            }
                        }
                        for(size_t s = 0; s < successors.size(); s++)
                        {
                            if(out0[successors[s]] == 0 && core0[predecessors[s]] == vf2_state::NULL_NODE())
                            {
                                out0[successors[s]] = depth;
                            }
                        }
                    }
                }
            }
    };

    public:

        typedef std::vector< vector< size_t > >::iterator iso_iterator;
        static const unsigned int max_count = 1000;
        CIsomorphism(const CNetwork* pg1, const CNetwork* pg2, SemanticFeasibility semantics = SemanticFeasibility(), EdgeFeasibility edge_check = EdgeFeasibility()) : m_pG1(pg1), m_pG2(pg2), m_SemanticCheck(semantics), m_EdgeFeasibility(edge_check), m_CurrentState(m_pG1, m_pG2), m_bswapped(false)
        {
            if(m_pG1->NumNodes() < m_pG2->NumNodes())
            {
                std::swap(m_pG1, m_pG2);
                m_bswapped = true;
                m_CurrentState = vf2_state(m_pG1, m_pG2);
            }
        }

        ~CIsomorphism() {}
        static bool is_null(size_t i) { return i == vf2_state::NULL_NODE(); }
        bool IsIsomorphic() { return ((m_pG1->NumNodes() == m_pG2->NumNodes()) && m_Mappings.size()); }

        // bool FindMappingSVD(EigenMatrixHelper::spMatrix& P)
        // {
        //     bool bIsomorphic = false;
        //     if(m_pG1->NumNodes() > 0 && (m_pG1->NumNodes() == m_pG2->NumNodes()))
        //     {
        //         Eigen::JacobiSVD<EigenMatrixHelper::matrix> svd_g1(m_pG1->GetAdjMatrix().toDense(), Eigen::ComputeFullU);
        //         Eigen::JacobiSVD<EigenMatrixHelper::matrix> svd_g2(m_pG2->GetAdjMatrix().toDense(), Eigen::ComputeFullU);
        //         // A = USV**T
        //         EigenMatrixHelper::spMatrix U1 = EigenMatrixHelper::Dense2Sparse(svd_g1.matrixU());
        //         EigenMatrixHelper::spMatrix U2 = EigenMatrixHelper::Dense2Sparse(svd_g2.matrixU());
        //         P = U2*(U1.transpose());
        //
        //         std::cout << "Permutation Matrix = " << std::endl;
        //         std::cout << P.toDense() << std::endl;
        //         bIsomorphic = EigenMatrixHelper::IsPermuationMatrix(P);
        //     }
        //     return bIsomorphic;
        // }
        void ComputeMatching()
        {
            Match();
        }

        size_t NumMappings() const { return m_Mappings.size(); }

        iso_iterator begin() { return m_Mappings.begin(); }
        iso_iterator end() { return m_Mappings.end(); }

    private:
    // Match()
    // algorithm taken from:
    // Cordella, Foggia, Sansone, and Vento. "A (Sub)Graph Isomorphism Algorithm for Matching Large Graphs". IEEE TPAMI 26.10.1367-72 (2004).
        void  Match()
        {
            if(m_CurrentState.ShouldTerminate())
            {
                if(m_bswapped)
                    m_Mappings.push_back(m_CurrentState.core2); // corresponds to the input g1->g2
                else
                    m_Mappings.push_back(m_CurrentState.core1);
                // std::cout << "mappings size: " << m_Mappings.size() << std::endl;
            }
            else
            {
                vector< pair<size_t, size_t> > candidates = ComputeCandidates();
                for(size_t i = 0; i < candidates.size(); i++)
                {
                    if(syntactic_feasibility(candidates[i].first, candidates[i].second)
                       && semantic_feasibility(candidates[i].first, candidates[i].second))
                    {
                        if(m_Mappings.size()>=max_count)
                            break;
                        vf2_state old_state(m_CurrentState); // this memory could add up in the recursion. should we update/restore the same memory?
                        m_CurrentState.inclusion(candidates[i], m_pG1, m_pG2);
                        Match();
                        m_CurrentState = old_state;
                    }
                }
            }
        }

        vector< pair<size_t, size_t> > ComputeCandidates()
        {
            // computes the candidates of the current state.
            vector< pair<size_t, size_t> > candidates;
            bool bFoundT2Pair = false;
            size_t min_t2 = m_pG2->NumNodes();
            // We want all the nodes in T1 = out1 - core1

            for(size_t t2 = 0; t2 < m_CurrentState.out2.size() && !bFoundT2Pair; t2++)
            {
                if(m_CurrentState.out2[t2] != 0 && m_CurrentState.core2[t2] == vf2_state::NULL_NODE()) // not in M2 but in Out means it is in T2
                {
                    if(t2 < min_t2)
                    {
                        bFoundT2Pair = true;
                        min_t2 = t2;
                    }
                }
            }

            for(size_t t1 = 0; t1 < m_CurrentState.out1.size() && bFoundT2Pair; t1++)
            {
                if(m_CurrentState.out1[t1] != 0 && m_CurrentState.core1[t1] == vf2_state::NULL_NODE()) // not in M1 but in Out means it is in T1
                {
                    candidates.push_back(pair<size_t, size_t>(t1, min_t2));
                }
            }

            if(candidates.empty()) // either T1out or T2out was empty.
            {
                bFoundT2Pair = false; // try with the T^{in} sets.
                min_t2 = m_pG2->NumNodes();
                for(size_t t2 = 0; t2 < m_CurrentState.in2.size() && !bFoundT2Pair; t2++)
                {
                    if(m_CurrentState.in2[t2] != 0 && m_CurrentState.core2[t2] == vf2_state::NULL_NODE()) // not in M2 but in Out means it is in T2
                    {
                        if(t2 < min_t2)
                        {
                            bFoundT2Pair = true;
                            min_t2 = t2;
                        }
                    }
                }

                for(size_t t1 = 0; t1 < m_CurrentState.in1.size() && bFoundT2Pair; t1++)
                {
                    if(m_CurrentState.in1[t1] != 0 && m_CurrentState.core1[t1] == vf2_state::NULL_NODE()) // not in M1 but in Out means it is in T1
                    {
                        candidates.push_back(pair<size_t, size_t>(t1, min_t2));
                    }
                }
            }

            if(candidates.empty()) // The graphs are disconnected!
            {
                bFoundT2Pair = false; // try with the T^{in} sets.
                min_t2 = m_pG2->NumNodes();
                for(size_t t2 = 0; t2 < m_CurrentState.out2.size() && !bFoundT2Pair; t2++)
                {
                    if(m_CurrentState.core2[t2] == vf2_state::NULL_NODE()) // not in M2
                    {
                        if(t2 < min_t2)
                        {
                            bFoundT2Pair = true;
                            min_t2 = t2;
                        }
                    }
                }
                for(size_t t1 = 0; t1 < m_CurrentState.out1.size() && bFoundT2Pair; t1++)
                {
                    if(m_CurrentState.core1[t1] == vf2_state::NULL_NODE()) // not in M1
                    {
                        candidates.push_back(pair<size_t, size_t>(t1, min_t2));
                    }
                }
            }

            return candidates; // will this be too large?
        }

    // inclusion rules:
        bool inclusion_check(const vector<size_t>& neighbors, const size_t& m)
        {
            assert(std::is_sorted(neighbors.begin(), neighbors.end()));
            return std::binary_search(neighbors.begin(), neighbors.end(), m);
        }

        bool pred_rule(const size_t& node_g1, const size_t& node_g2, const vector<size_t>& pred1, const vector<size_t>& pred2)
        {
            // For every n' \in M1 \cap Pred(G1, node_g1) check if there exists
            //           m' \in M2 \cap Pred(G2, node_g2) such that
            //           n' -> node_g1 == m' -> node_g2
            bool bFeasible = true;
            for(size_t n = 0; n < pred1.size() && bFeasible; n++)
            {
                size_t m = m_CurrentState.core1[pred1[n]];
                if( m != vf2_state::NULL_NODE()) // pred1[n] is in M1
                {
                    // predi -> node_gi
                    bFeasible = inclusion_check(pred2, m) && edge_feasibility(pred1[n], node_g1, m, node_g2);
                }
            }

            for(size_t m = 0; m < pred2.size() && bFeasible; m++)
            {
                size_t n = m_CurrentState.core2[pred2[m]];
                if( n != vf2_state::NULL_NODE()) // pred2[m] is in M2
                {
                    // predi -> node_gi
                    bFeasible = inclusion_check(pred1, n) && edge_feasibility(n, node_g1, pred2[m], node_g2);
                }
            }
            return bFeasible;
        }

        bool succ_rule(const size_t& node_g1, const size_t& node_g2, const vector<size_t>& succ1, const vector<size_t>& succ2 )
        {
            // For every n' \in M1 \cap Succ(G1, node_g1) check if there exists
            //           m' \in M2 \cap Succ(G2, node_g2) such that
            //           node_g1 -> n' ==  node_g2 -> m'

            bool bFeasible = true;
            for(size_t n = 0; n < succ1.size() && bFeasible; n++)
            {
                size_t m = m_CurrentState.core1[succ1[n]];
                if( m != vf2_state::NULL_NODE()) // succ1[n] is in M1
                {
                    // node_gi -> predi
                    bFeasible = inclusion_check(succ2, m) && edge_feasibility(node_g1, succ1[n], node_g2, m);
                }
            }

            for(size_t m = 0; m < succ2.size() && bFeasible; m++)
            {
                size_t n = m_CurrentState.core2[succ2[m]];
                if( n != vf2_state::NULL_NODE()) // pred2[m] is in M2
                {
                    // node_gi -> predi
                    bFeasible = inclusion_check(succ1, n) && edge_feasibility(node_g1, n, node_g2, succ2[m]);
                }
            }
            return bFeasible;
        }

        bool in_rule(const vector<size_t>& pred1, const vector<size_t>& pred2, const vector<size_t>& succ1, const vector<size_t>& succ2)
        {
            size_t succCard1 = 0, succCard2 = 0;
            size_t predCard1 = 0, predCard2 = 0;
            for (size_t i = 0; i < succ1.size(); i++){
                if(m_CurrentState.IsInTin1(succ1[i])){
                    succCard1++;
                }
            }
            for (size_t i = 0; i < succ2.size(); i++){
                if(m_CurrentState.IsInTin2(succ2[i])){
                    succCard2++;
                }
            }
            for (size_t i = 0; i < pred1.size(); i++){
                if(m_CurrentState.IsInTin1(pred1[i])){
                    predCard1++;
                }
            }
            for (size_t i = 0; i < pred2.size(); i++){
                if(m_CurrentState.IsInTin2(pred2[i])){
                    predCard2++;
                }
            }

            return (succCard1 >= succCard2) && (predCard1 >= predCard2);
        }

        bool out_rule(const vector<size_t>& pred1, const vector<size_t>& pred2, const vector<size_t>& succ1, const vector<size_t>& succ2)
        {
          size_t succCard1 = 0, succCard2 = 0;
          size_t predCard1 = 0, predCard2 = 0;
          for (size_t i = 0; i < succ1.size(); i++){
              if(m_CurrentState.IsInTout1(succ1[i])){
                  succCard1++;
              }
          }
          for (size_t i = 0; i < succ2.size(); i++){
              if(m_CurrentState.IsInTout2(succ2[i])){
                  succCard2++;
              }
          }
          for (size_t i = 0; i < pred1.size(); i++){
              if(m_CurrentState.IsInTout1(pred1[i])){
                  predCard1++;
              }
          }
          for (size_t i = 0; i < pred2.size(); i++){
              if(m_CurrentState.IsInTout2(pred2[i])){
                  predCard2++;
              }
          }
          return (succCard1 >= succCard2) && (predCard1 >= predCard2);
        }

        bool new_rule(const vector<size_t>& pred1, const vector<size_t>& pred2, const vector<size_t>& succ1, const vector<size_t>& succ2)
        {
            size_t succCard1 = 0, succCard2 = 0;
            size_t predCard1 = 0, predCard2 = 0;
            for (size_t i = 0; i < succ1.size(); i++){
                if(m_CurrentState.IsInNc1(succ1[i])){
                    succCard1++;
                }
            }
            for (size_t i = 0; i < succ2.size(); i++){
                if(m_CurrentState.IsInNc2(succ2[i])){
                    succCard2++;
                }
            }
            for (size_t i = 0; i < pred1.size(); i++){
                if(m_CurrentState.IsInNc1(pred1[i])){
                    predCard1++;
                }
            }
            for (size_t i = 0; i < pred2.size(); i++){
                if(m_CurrentState.IsInNc2(pred2[i])){
                    predCard2++;
                }
            }
            return (succCard1 >= succCard2) && (predCard1 >= predCard2);
        }
    // feasibility rules:
        bool edge_feasibility(const size_t& src_g1, const size_t& tgt_g1, const size_t& src_g2, const size_t& tgt_g2)
        {
            CNetwork::edge_info edge_g1, edge_g2;
            edge_g1.source = src_g1; edge_g1.target = tgt_g1; edge_g1.weight = (*m_pG1)(src_g1, tgt_g1);
            edge_g2.source = src_g2; edge_g2.target = tgt_g2; edge_g2.weight = (*m_pG2)(src_g2, tgt_g2);
            return m_EdgeFeasibility(edge_g1, edge_g2);
        }

        bool semantic_feasibility(const size_t& node_g1, const size_t& node_g2)
        {
            return m_SemanticCheck(node_g1, node_g2);
        }

        bool syntactic_feasibility(const size_t& node_g1, const size_t& node_g2)
        {
            const vector<size_t>& pred1 = m_pG1->GetNeighbors(node_g1, true);
            const vector<size_t>& pred2 = m_pG2->GetNeighbors(node_g2, true);

            const vector<size_t>& succ1 = m_pG1->GetNeighbors(node_g1, false);
            const vector<size_t>& succ2 = m_pG2->GetNeighbors(node_g2, false);

            return (pred_rule(node_g1, node_g2, pred1, pred2) &&
                    succ_rule(node_g1, node_g2, succ1, succ2) &&
                    in_rule(pred1, pred2, succ1, succ2)   &&
                    out_rule(pred1, pred2, succ1, succ2)  &&
                    new_rule(pred1, pred2, succ1, succ2));
        }

    private:
        const CNetwork* m_pG1;
        const CNetwork* m_pG2;
        SemanticFeasibility m_SemanticCheck;
        EdgeFeasibility     m_EdgeFeasibility;
        vector< vector< size_t > > m_Mappings; // store all the mappings? or do we need to make this a generator class.

        // members specific to the VF2 algorithm
        // these vectors define the current state of the algorithm.
        vf2_state m_CurrentState;
        bool m_bswapped;
};

// Implementation to aviod circular references.
// This may be pulled into the class now that the iterator definition is in the class.
inline bool CNetwork::KruskalAlgorithm(CNetwork& SpanningForest)
{
//    std::cout << "N - KruskalAlgorithm" << std::endl;
    bool bSuccess = true;
    SpanningForest.Initialize(int(NumNodes()));
    std::multiset<CNetwork::edge_info> edges;
    vector<size_t> nodes_tree_id;
//    size_t treeCt = 0;

    for( CNetwork::edge_iterator iter(this); iter.IsValid(); iter++) edges.insert(*iter);      // sort edges by weight
    for( size_t i = 0; i < NumNodes(); i++) nodes_tree_id.push_back(i);             // Each vertex is its own tree.

    for( std::multiset<CNetwork::edge_info>::iterator iter = edges.begin(); iter != edges.end(); iter++)
    {
//        std::cout << "src: " << iter->source << " tgt: "<< iter->target << " weight: "<< iter->weight << std::endl;
        size_t src = size_t(iter->source);
        size_t tgt = size_t(iter->target);
        if(nodes_tree_id[src] != nodes_tree_id[tgt])
        {
            SpanningForest.AddEdge(src, tgt, iter->weight);

            // "join" trees together.
            size_t set_id, replace_id;
            set_id = nodes_tree_id[src] < nodes_tree_id[tgt] ? nodes_tree_id[src] : nodes_tree_id[tgt];
            replace_id = nodes_tree_id[src] > nodes_tree_id[tgt] ? nodes_tree_id[src] : nodes_tree_id[tgt];
            for(size_t i = 0; i < nodes_tree_id.size(); i++ )
            {
                if(nodes_tree_id[i] == replace_id) nodes_tree_id[i] = set_id;
            }
        }
    }
    return bSuccess;
}


inline bool IsCycle(const CNetwork& net, unsigned int n)
{
    CNetwork cycle = CNetwork::Cycle(n, net.GetNetworkMask());
    CIsomorphism<> iso(&net, &cycle);
    iso.ComputeMatching();
    return iso.IsIsomorphic();
}


}
}
#endif
