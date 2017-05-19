
#pragma once

#include "Network.h"
#include <pynomial/extern/nanoflann/include/nanoflann.hpp>
namespace pynomial{
namespace graph{

template<typename T>
class GraphDistance
{
public:
    GraphDistance() {}
    T operator()(const CNetwork& g1, const CNetwork& g2)
    {
        T dist = 10.0;
        const CNetwork* pg1 = &g1, *pg2 = &g2;
        if(pg1->NumNodes() < pg2->NumNodes())
        {
            std::swap(pg1, pg2);
        }


        CIsomorphism<> iso(pg1,pg2); // g1 -> g2
        iso.ComputeMatching();
        std::cout << "g1 = " << g1.NumNodes() << ", g2 = " << g2.NumNodes() << ", mappings = " << iso.NumMappings() << std::endl;

        for(CIsomorphism<>::iso_iterator it = iso.begin(); it != iso.end(); it++)
        {
            T d = 0.0;
            for(CNetwork::edge_iterator edge(pg1); edge.IsValid(); edge++) // TODO: think about this more and make sure it is logically correct.
            {
                size_t src = it->at(edge->source), tgt = it->at(edge->target);
                if(CIsomorphism<>::is_null(src) || CIsomorphism<>::is_null(tgt))
                {
                    std::cout << "w = " << edge->weight << std::endl;
                    d += edge->weight;
                }
                else{
                    double w = pg2->GetAdjMatrix().coeff(tgt, src);
                    std::cout << "diff = " << fabs(edge->weight - w) << std::endl;
                    d += fabs(edge->weight - w);
                }

            }
            // for(size_t i = 0; i < it->size(); i++)
            // {
            //     if(CIsomorphism<>::is_null(it->at(i)))
            //     {
            //         std::cout << "null" << ", ";
            //     }
            //     else
            //     {
            //         std::cout << it->at(i) << ", ";
            //     }
            // }
            // std::cout <<" d = "<< d << std::endl;
            if(it == iso.begin())
            {
                dist = d;
                // mapping = *it;
            }
            else if (dist > d)
            {
                dist = d;
                // mapping = *it; // potentially doing many copies.
            }
        }
        return dist;
    }
    std::vector<size_t> mapping;
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
        if(m_npoints <=1) throw std::runtime_error("matching requires more than one point");
        m_dist.clear();
        m_dist.reserve(m_npoints*(m_npoints-1)/2);
        GraphDistance<double> d;
        for(size_t i = 0; i < m_npoints-1; i++)
        {
            CNetwork gi;
            gi.AddNetworkType(m_network.GetNetworkMask());
            vector<size_t> ni = m_network.GetNeighbors(i);

            ni.push_back(i);
            std::sort(ni.begin(), ni.end());

            // for(size_t k = 0; k < ni.size(); k++)
            // {
            //     std::cout << "w("<< i << ", "<< ni[k]<< ")="<<m_network.GetAdjMatrix().coeff(ni[k],i) << std::endl;
            // }
            m_network.SubGraphEx(ni.begin(), ni.end(), gi);
            // m_network.SubGraph(ni.begin(), ni.end(), gc);
            std::cout << "**************** SubGraph "<< i <<" ************* " << std::endl;
            gi.PrintNetwork();
            // std::cout << "**************** check ************* " << std::endl;
            // gc.PrintNetwork();
            std::cout << "processing particle " << i << " with "<< ni.size() << " neighbors" <<'\r' << std::flush;
            // std::cout << "processing particle " << i << '\r' << std::flush;
            for(size_t j = i+1; j < m_npoints; j++)
            {
                // break;
                CNetwork gj;
                gj.AddNetworkType(m_network.GetNetworkMask());
                vector<size_t> nj = m_network.GetNeighbors(j);
                nj.push_back(j);
                std::sort(nj.begin(), nj.end());
                m_network.SubGraphEx(nj.begin(), nj.end(), gj);
                std::cout << "**************** SubGraph "<< j <<" ************* " << std::endl;
                gj.PrintNetwork();

                double dist = d(gi,gj);
                std::cout << "processing particle " << i << " with "<< j << " dist = "<< dist <<'\r' << std::endl;
                m_dist.push_back(dist);
                if(j>i+3) return;
                // return;
            }
            // if(i > 1)
            //     break;
        }
        std::cout << std::endl;
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
