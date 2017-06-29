
#pragma once

#include "Network.h"
#include <pynomial/extern/nanoflann/include/nanoflann.hpp>
#include "../register/brute_force.h"
namespace pynomial{
namespace graph{

template<typename T>
struct edit_path{
    // Some constanst to store some state information.
    static const size_t null;//     = -1; see below.
    static const size_t notset;//   = -2;
    // member varaibles
    std::vector<size_t> forward; // will have g1.NumNodes()
    std::vector<size_t> reverse; // will have g2.NumNodes()
    T cost;
    // some helper methods
    edit_path(size_t n1, size_t n2, const T& c) : forward(n1, notset), reverse(n2, notset), cost(c) {}

    edit_path(const edit_path<T>& other) : forward(other.forward), reverse(other.reverse), cost(other.cost) {}

    edit_path(edit_path<T>&& other) : forward(std::move(other.forward)), reverse(std::move(other.reverse)), cost(other.cost) {}

    bool operator < (const edit_path<T>& other){ return cost < other.cost; }

    bool operator > (const edit_path<T>& other){ return cost > other.cost; }

    size_t& operator[](const size_t& i){ return forward[i]; }

    size_t& operator()(const size_t& i){ return reverse[i]; }

    bool complete() const
    {
        for(size_t val : forward)
            if(val == notset)
                return false;
        for(size_t val : reverse)
            if(val == notset)
                return false;
        return true;
    }

    size_t count_mapped() const
    {
        return std::count(forward.begin(), forward.end(), notset);
    }

    edit_path<T>& flip()
    {
        std::vector<size_t> f = forward;
        forward = reverse;
        reverse = f;
        return *this;
    }

};
template< typename T>
const size_t edit_path<T>::null = -1;
template<typename T>
const size_t edit_path<T>::notset = -2;


template<typename T>
class GraphEditDistance
{
public:
    typedef edit_path<T> edit_path_t;
    typedef std::function<T(const edit_path_t&, const CNetwork&, const CNetwork&)> CostFunction;
private:
        CostFunction    m_Cost;
        CostFunction    m_Heuristic;
public:
    GraphEditDistance() {}

    edit_path_t operator()(const CNetwork& g1, const CNetwork& g2)
    {
        if(g1.NumNodes() > g2.NumNodes())
        {
            return (*this)(g2, g1).flip();
        }
        std::set< edit_path_t > open;
        for(size_t i = 0; g2.NumNodes(); i++)
        {
            edit_path_t path(g1.NumNodes(), g2.NumNodes(), 0.0);
            path[0] = i; // set the forward and reverse maps.
            path(i) = 0;
            path.cost = m_Cost(path, g1, g2);
            open.insert(path);
        }
        edit_path_t path(g1.NumNodes(), g2.NumNodes(), 0.0);
        path[0] = edit_path_t::null; // deletion of the node.
        path.cost = m_Cost(path, g1, g2);
        open.insert(path);
        edit_path_t min_path;
        size_t n = g1.NumNodes();
        while(true)
        {
            if(open.begin() == open.end())
            {
                std::cout << "ERROR!! Could not find a minimum edit path!" << std::endl;
                throw(std::runtime_error("ERROR!! Could not find a minimum edit path!"));
            }
            min_path = *open.begin();
            open.erase(open.begin());
            if( min_path.complete() ) // we found it!
            {
                break;
            }
            size_t k = min_path.count_mapped();
            if(k < n)
            {
                for(size_t v2 = 0; v2 < min_path.reverse.size(); v2++)
                {
                    if(min_path.reverse[v2] == edit_path_t::notset) // we have not mapped this one yet.
                    {
                        edit_path_t new_path(min_path);
                        new_path[k] = v2;
                        new_path(v2) = k;
                        new_path.cost = m_Cost(new_path, g1, g2);
                        open.insert(std::move(new_path)); // use move semantics to avoid another copy.
                    }
                }
                edit_path_t new_path(min_path);
                new_path[k] = edit_path_t::null; // delete the node.
                new_path.cost = m_Cost(new_path, g1, g2);
                open.insert(std::move(new_path)); // use move semantics to avoid another copy.
            }
            else
            {
                // Here we have mapped all nodes in g1 -> g2
                // but there are some nodes still in g2 that must be deleted.
                edit_path_t new_path(min_path);
                for(size_t v2 = 0; v2 < min_path.reverse.size(); v2++)
                {
                    if(min_path.reverse[v2] == edit_path_t::notset)
                    {
                        new_path(v2) = edit_path_t::null; // delete the nodes.
                    }
                }
                assert(new_path.complete()); // now all nodes will be mapped.
                new_path.cost = m_Cost(new_path, g1, g2);
                open.insert(std::move(new_path));
            }
        }
        return min_path;
    }

    void SetCost( CostFunction func ) { m_Cost = func; }

    void SetHeuristic( CostFunction func ) { m_Heuristic = func; }
};

class CostFunction {
/*
Some notes and properties that should be held by the cost functions.


1.  c(e) >= 0, for all node and edge edit operations e.
2.  c(e) > 0, for all node and edge deletions and insertions e.
3.  c(u -> w) <= c(u -> v) + c(v -> w)
    c(u -> ε) <= c(u -> v) + c(v -> 0)
    c(0 -> v) <= c(0 -> u) + c(u -> v)

2 and 3 make sure we dont make superfluous substitutions.

*/
};

template< typename Label>
class MinkowskiCost : CostFunction
{
    double m_tau;
public:
    MinkowskiCost(double tau) : m_tau(tau) {}
    double operator ()(const GraphEditDistance<double>::edit_path_t& path, const CNetwork& g1, const CNetwork& g2)
    {
        double cost = 0.0;

        return cost;
    }
};

class NN
{
public:
    typedef Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic>  matrix_t;
    typedef nanoflann::KDTreeEigenMatrixAdaptor< matrix_t >  kd_tree_t;

    NN(const matrix_t& points, unsigned int npoints, double r_cut) : m_dim(0), m_npoints(npoints), m_network(points.rows(), NetworkType::Undirected | NetworkType::Weighted), m_tree(nullptr), m_rcut(r_cut)
    {
        std::cout << "NN with " << npoints << std::endl;
        m_bInitialized = false;
        if(points.rows() > 0)
        {
            if(points.cols() > 0)
            {
                m_dim = points.cols();
                m_points = matrix_t::Zero(points.rows(), points.cols());
                m_bInitialized = true;
            }
        }
        if(m_bInitialized)
        {
            m_points = points;
            // for(size_t i = 0; i < points.size(); i++)
            // {
            //     if(points[i].size() != m_dim) throw std::runtime_error("each point must be the same dimension");
            //     for(size_t j = 0; j < points[i].size(); j++)
            //     {
            //         m_points(i,j) = points[i][j];
            //     }
            // }
            std::unique_ptr<double [] > query_pt(new double[m_dim]);
            m_tree = new kd_tree_t(m_dim, m_points, 20);
            m_tree->index->buildIndex();
            nanoflann::SearchParams params;
            std::vector<std::pair<long,double> >   ret_matches;
            for(size_t i = 0; i < npoints; i++)
            {
                ret_matches.clear();
                for(size_t j = 0; j < m_dim; j++) query_pt[j] = points(i,j);
                params.sorted = false;
                size_t nMatches = m_tree->index->radiusSearch(&query_pt[0], m_rcut*m_rcut, ret_matches, params);
                std::cout << "n =" << nMatches << " for particle "<< i << std::endl;
                for(size_t j = 0; j < nMatches; j++) {
                    if(ret_matches[j].first > i)
                    {
                        // std::cout << i << ": idx["<< j << "]=" << ret_matches[j].first << " dist["<< j << "]=" << ret_matches[j].second << std::endl;
                        // double dist =
                        m_network.AddEdge(i, ret_matches[j].first, ret_matches[j].second); //TODO: make this more efficient
                    }
                }
            }
        }
    }
    ~NN()
    {
        if(m_tree)
            delete m_tree;
    }
    void compute();
    void match()
    {

        //

        // if(m_npoints <=1) throw std::runtime_error("matching requires more than one point");
        // m_dist.clear();
        // m_dist.reserve(m_npoints*(m_npoints-1)/2);
        // GraphDistance<double> d;
        // for(size_t i = 0; i < m_npoints-1; i++)
        // {
        //     CNetwork gi;
        //     gi.AddNetworkType(m_network.GetNetworkMask());
        //     vector<size_t> ni = m_network.GetNeighbors(i);
        //
        //     ni.push_back(i);
        //     std::sort(ni.begin(), ni.end());
        //
        //     // for(size_t k = 0; k < ni.size(); k++)
        //     // {
        //     //     std::cout << "w("<< i << ", "<< ni[k]<< ")="<<m_network.GetAdjMatrix().coeff(ni[k],i) << std::endl;
        //     // }
        //     m_network.SubGraphEx(ni.begin(), ni.end(), gi);
        //     // m_network.SubGraph(ni.begin(), ni.end(), gc);
        //     std::cout << "**************** SubGraph "<< i <<" ************* " << std::endl;
        //     gi.PrintNetwork();
        //     // std::cout << "**************** check ************* " << std::endl;
        //     // gc.PrintNetwork();
        //     std::cout << "processing particle " << i << " with "<< ni.size() << " neighbors" <<'\r' << std::flush;
        //     // std::cout << "processing particle " << i << '\r' << std::flush;
        //     for(size_t j = i+1; j < m_npoints; j++)
        //     {
        //         // break;
        //         CNetwork gj;
        //         gj.AddNetworkType(m_network.GetNetworkMask());
        //         vector<size_t> nj = m_network.GetNeighbors(j);
        //         nj.push_back(j);
        //         std::sort(nj.begin(), nj.end());
        //         m_network.SubGraphEx(nj.begin(), nj.end(), gj);
        //         std::cout << "**************** SubGraph "<< j <<" ************* " << std::endl;
        //         gj.PrintNetwork();
        //
        //         double dist = d(gi,gj);
        //         std::cout << "processing particle " << i << " with "<< j << " dist = "<< dist <<'\r' << std::endl;
        //         m_dist.push_back(dist);
        //         if(j>i+3) return;
        //         // return;
        //     }
        //     // if(i > 1)
        //     //     break;
        // }
        // std::cout << std::endl;
    }
    size_t num_points() {return m_npoints; }
    const std::vector<double>& pair_distance() {return m_dist; }
protected:
    size_t      m_dim;
    size_t      m_npoints;
    CNetwork    m_network;
    matrix_t    m_points;
    kd_tree_t*  m_tree;
    bool m_bInitialized;
    double      m_rcut;
    std::vector<double> m_dist;
    std::vector<double> m_rmsd;
};

}}



/*

Eigen::Matrix<num_t,Dynamic,Dynamic>  mat(nSamples,dim);

	const num_t max_range = 20;

	// Generate points:
	generateRandomPointCloud(mat, nSamples,dim, max_range);

//	cout << mat << endl;

	// Query point:
	std::vector<num_t> query_pt(dim);
	for (size_t d=0;d<dim;d++)
		query_pt[d] = max_range * (rand() % 1000) / num_t(1000);


	// ------------------------------------------------------------
	// construct a kd-tree index:
	//    Some of the different possibilities (uncomment just one)
	// ------------------------------------------------------------
	// Dimensionality set at run-time (default: L2)
	typedef KDTreeEigenMatrixAdaptor< Eigen::Matrix<num_t,Dynamic,Dynamic> >  my_kd_tree_t;

	// Dimensionality set at compile-time
//	typedef KDTreeEigenMatrixAdaptor< Eigen::Matrix<num_t,Dynamic,Dynamic>, SAMPLES_DIM>  my_kd_tree_t;

	// Dimensionality set at compile-time: Explicit selection of the distance metric: L2
//	typedef KDTreeEigenMatrixAdaptor< Eigen::Matrix<num_t,Dynamic,Dynamic>, SAMPLES_DIM,nanoflann::metric_L2>  my_kd_tree_t;

	// Dimensionality set at compile-time: Explicit selection of the distance metric: L2_simple
//	typedef KDTreeEigenMatrixAdaptor< Eigen::Matrix<num_t,Dynamic,Dynamic>, SAMPLES_DIM,nanoflann::metric_L2_Simple>  my_kd_tree_t;

	// Dimensionality set at compile-time: Explicit selection of the distance metric: L1
//	typedef KDTreeEigenMatrixAdaptor< Eigen::Matrix<num_t,Dynamic,Dynamic>, SAMPLES_DIM,nanoflann::metric_L1>  my_kd_tree_t;

	my_kd_tree_t   mat_index(dim , mat, 10  );
	mat_index.index->buildIndex();

	// do a knn search
	const size_t num_results = 3;
	vector<size_t>   ret_indexes(num_results);
	vector<num_t> out_dists_sqr(num_results);

	nanoflann::KNNResultSet<num_t> resultSet(num_results);

	resultSet.init(&ret_indexes[0], &out_dists_sqr[0] );
	mat_index.index->findNeighbors(resultSet, &query_pt[0], nanoflann::SearchParams(10));

	std::cout << "knnSearch(nn="<<num_results<<"): \n";
	for (size_t i=0;i<num_results;i++)
		std::cout << "ret_index["<<i<<"]=" << ret_indexes[i] << " out_dist_sqr=" << out_dists_sqr[i] << endl;


        */
