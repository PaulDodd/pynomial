#pragma once
// stdlib include
#include <iostream>
#include <vector>
#include <random>
#include <algorithm>
#include <chrono>
// eigen include
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Sparse>
// external includes
#include <pynomial/extern/nanoflann/include/nanoflann.hpp>

namespace pynomial{
//
namespace _register{

using std::vector;

typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> matrix;

// TODO: make these functions templates and move to a different location.
inline matrix CenterOfMass(const matrix& P)
{
    // Assumes that P = (v**T) if v is a column vector.  or in other notation  P = [x1, y1, z1; ...]
    matrix cm(1, P.cols());
    for(int i =0; i < P.cols(); i++)
        cm(0,i) = P.col(i).sum()/double(P.rows());
    //cout << "cm = \n"<< cm << endl;

    return cm;
}

inline matrix Translate(const matrix& vec, const matrix& P)
{
    matrix trans = matrix::Zero(P.rows(), P.cols());
    for(int i = 0; i < P.rows(); i++)
        trans.row(i) = P.row(i)+vec;
    return trans;
}

inline double RMSD(const matrix& P, const matrix& Q)
{
    matrix pmq = P-Q;
    return pmq.norm()/double(P.rows());
}

inline void KabschAlgorithm(const matrix& P, const matrix& Q, matrix& Rotation)
{
    // Preconditions: P and Q have been translated to have the same center of mass.
    matrix A = P.transpose()*Q;
    Eigen::JacobiSVD<matrix> svd(A, Eigen::ComputeFullU | Eigen::ComputeFullV );
    // A = USV**T
    matrix U = svd.matrixU();
    matrix V = svd.matrixV();

    double det = (V*U.transpose()).determinant();

    if(det < 0)
    {
        //std::cout << "Vb = \n" << V << std::endl<< std::endl;
        V.col(V.cols()-1)*= -1.0;
        //std::cout << "Va = \n" << V << std::endl<< std::endl;
    }
    //Solve for the rotation matrix
    Rotation= V*U.transpose();
}

inline void AlignVectorSets(matrix& P,matrix& Q, matrix* pRotation = NULL)
{
    // Aligns p with q.
    // both p and q will be changed in this operation.

    matrix rotation;
    //Translate both p,q to origin.
    P = Translate(-CenterOfMass(P), P);
    Q = Translate(-CenterOfMass(Q), Q);
    KabschAlgorithm(P, Q, rotation); // Find the rotation.
    P = (rotation*P.transpose()).transpose();  // Apply the transform

    if(pRotation) // optionally copy the roation.
        *pRotation = rotation;
}

class RegisterBruteForce  // : public Register
{
    typedef nanoflann::KDTreeEigenMatrixAdaptor< matrix >  kd_tree_t;
    public:
        RegisterBruteForce(const matrix& pts, double threshold=1e-6) : m_data(pts), m_rmsd(0.0), m_tol(threshold), m_shuffles(1), m_tree(nullptr)
        {
            m_dim = m_data.cols();
            assert(m_dim == 2 || m_dim == 3);
            m_tree = new kd_tree_t(m_dim, m_data, 20);
            m_tree->index->buildIndex();
        }
        ~RegisterBruteForce()
        {
            if(m_tree)
                delete m_tree;
        }
        bool Fit(const matrix& pts)
        {
            // Here we find a matching, that finds the
            const matrix& points = pts;
            matrix p, q, r;
            int Np = points.rows();
            int Nd = m_data.rows();
            assert(points.cols() == m_data.cols());
            // if(N != m_data.rows())
            // {
            //     std::cout << "brute force matching requires the same number of points" << std::endl;
            //     return false;
            // }
            //
            RandomNumber<std::mt19937_64> rng;
            double rmsd_min = -1.0;
            // std::cout << "making "<< m_shuffles << " trial moves" << std::endl;
            for(size_t shuffles = 0; shuffles < m_shuffles; shuffles++)
            {
                int p0 = 0, p1 = 0, p2 = 0;
                if(shuffles == 0) // This is so you can register the current map if you set shuffles=1
                    {
                    p0 = 0, p1 = 1, p2 = 2;
                    }

                while( p0 == p1 || p0 == p2 || p1 == p2)
                {
                    p0 = rng.random_int(0,Nd-1);
                    p1 = rng.random_int(0,Nd-1);
                    p2 = rng.random_int(0,Nd-1);
                }

                size_t comb[3] = {0, 1, 2};
                p.resize(3,m_data.cols());
                p.row(0) = m_data.row(p0);
                p.row(1) = m_data.row(p1);
                p.row(2) = m_data.row(p2);
                q.resize(3,m_data.cols());
                do
                {
                    do {
                        q.row(0) = points.row(comb[0]);
                        q.row(1) = points.row(comb[1]);
                        q.row(2) = points.row(comb[2]);

                        KabschAlgorithm(p, q, r);
                        double dist = AlignedRSMDTree(points, r);
                        if(dist < rmsd_min || rmsd_min < 0.0)
                        {
                            std::cout << "rmsd = " <<  dist << "          \r" << std::flush;
                            rmsd_min = dist;
                            m_rotation = r;
                            m_rmsd = rmsd_min;
                            if(rmsd_min < m_tol)
                            {
                                return true;
                            }
                        }
                    } while ( std::next_permutation(comb,comb+3) );
                }while(NextCombination(comb, Np, 3)); // O(Np^3)
            }
            return false;
        }

        matrix GetRotation() { return m_rotation; }

        double GetCost() { return m_rmsd; }

        void SetNumShuffles(size_t s) { m_shuffles = s; }

        double HausdorffDist(const matrix& points, const matrix& rot)
        {
            // computes the hausdorff distance:
            // h(X, Y) = max {sup_X inf_Y d(x,y), sup_Y inf_X d(x,y)}
            // h(X,Y) is the hausdorff distance
            // d(x,y) is the euclidian distance
            // x,y are points in the repective sets X, Y.
            assert(rot.cols() && rot.rows());
            double dist = 0.0;
            nanoflann::SearchParams params;
            matrix transformed = rot*points.transpose();
            transformed.transposeInPlace();
            // sup_X inf_Y d(x,y)
            // std::cout << "data: \n" << m_data << "\npoints: \n" << points << "\n " << "\nrot: \n " << rot << "\ntransformed: \n"<< transformed << std::endl;
            for(int r = 0; r < points.rows(); r++)
            {
                double distsq;
                size_t ret_index;
                nanoflann::KNNResultSet<double> resultSet(1);
                resultSet.init(&ret_index, &distsq);
                Eigen::Ref<const Eigen::VectorXd> ref(transformed.row(r));
                m_tree->index->findNeighbors(resultSet, ref.data(), params);
                dist = fmax(dist, distsq);
            }
            // sup_Y inf_X d(x,y)
            kd_tree_t other(m_dim, transformed, 20);
            other.index->buildIndex();
            for(int r = 0; r < m_data.rows(); r++)
            {
                double distsq;
                size_t ret_index;
                nanoflann::KNNResultSet<double> resultSet(1);
                resultSet.init(&ret_index, &distsq);
                Eigen::Ref<const Eigen::VectorXd> ref(m_data.row(r));
                other.index->findNeighbors(resultSet, ref.data(), params);
                dist = fmax(dist, distsq);
            }
            return sqrt(dist); // could just return the sqrd distance but I want to make sure to compare tol to d rather than d^2
        }

        double AlignedRSMD(const matrix& points, const matrix& rot)
        {
            // As named we will do the brute force algorithm here as well.
            // we can optimize it later.
            assert(points.rows() == m_data.rows());
            double rmsd = 0.0;
            for(int r = 0; r < m_data.rows(); r++)
            {
                double dist = -1.0;
                for(int pr = 0; pr < points.rows(); pr++)
                {
                    Eigen::VectorXd pfit = rot*(points.row(pr).transpose());
                    Eigen::VectorXd delta = pfit - m_data.row(r).transpose();

                    // here we should really resolve the mapping and use that.
                    // this can map 2 points to one.
                    if(delta.norm() < dist || dist < 0.0)
                        dist = delta.norm();
                }
                rmsd += dist*dist;
            }

            return sqrt(rmsd/double(points.rows()));
        }

        double AlignedRSMDTree(const matrix& points, const matrix& rot, size_t num_search = 20)
        {
        // This is not exactly what we want to do. The correct matching
        // should be found through the hungarian algorithm, which is O(N^3).
        // this is a O(N*O(search)), note I am not sure what O(search) is but should be
        // less than O(N^2) so this will be faster but will probably give incorrect answers
        // So be aware!
        // the idea here is based on iterative closest point. we will match each point in m_data to the nearest
        // point in points that has not been matched yet. In the case that m_data is not the same size as points
        // we will neglect the extra points when calculating the RMSD.
        // This is also the case when large errors can occur, since the matching will be dependent on the
        // order of the points.

            if( num_search > m_data.rows() )
                num_search = m_data.rows();

            int num_mapped = 0; // little check sum to see if we found all we could
            double rmsd = 0.0;
            std::vector<bool> found(m_data.rows(), false);
            nanoflann::SearchParams params;
            for(int r = 0; r < points.rows(); r++)
            {
                double distsq = -1.0;
                vector<size_t>   ret_indexes(num_search);
                vector<double> out_dists_sqr(num_search);
                nanoflann::KNNResultSet<double> resultSet(num_search);
                resultSet.init( &ret_indexes[0], &out_dists_sqr[0]);

                Eigen::VectorXd query = rot*(points.row(r).transpose()); // query point.
                Eigen::Ref<const Eigen::VectorXd> ref(query);
                m_tree->index->findNeighbors(resultSet, ref.data(), params);
                for (size_t i = 0; i < num_search; i++)
                {
                    if(!found[ret_indexes[i]])
                    {
                        num_mapped++;
                        distsq = out_dists_sqr[i];
                        found[ret_indexes[i]] = true;
                        break;
                    }
                }

                if(distsq < 0.0)
                {
                    distsq = 0.0; // could not map point means that
                }
                rmsd += distsq;
            }
            // assert(num_mapped == std::min(points.rows(), m_data.rows()));
            return sqrt(rmsd/double(num_mapped));
        }


    private:
        // TODO: move this into its own class
        // TODO: class Combination(N, k), class Permutation(N)
        inline bool NextCombination(size_t* comb, int N, int k)
        {
            //    returns next combination.
            if(k == 0 || N == 0 || !comb)
                return false;

            bool bRetVal = false;

            for(int i = k-1; i >= 0; i--) {
                if(comb[i] + 1 < size_t(N+i-k+1)) {
                    comb[i]++;
                    for (int j = i+1; j < k; j++) {
                        comb[j] = comb[j-1]+1;
                    }
                    bRetVal = true;
                    break;
                }
            }

            return bRetVal;
        }
        template<class RNG>
        class RandomNumber
        {
        public:
            RandomNumber() { seed_generator(); }
            int random_int(int a, int b)
            {
                std::uniform_int_distribution<int> distribution(a,b);
                return distribution(m_generator);
            }
        private:
            inline void seed_generator(const size_t& n = 100)
            {
                std::vector<size_t> seeds;
                try {
                    std::random_device rd;
                    for(size_t i = 0; i < n; i++)
                        seeds.push_back(rd());
                } catch (...) {
                    std::cout << "random_device is not available..." << std::endl;
                    seeds.push_back(size_t(std::chrono::system_clock::now().time_since_epoch().count()));
                    // seeds.push_back(size_t(getpid()));
                }
                std::seed_seq seq(seeds.begin(), seeds.end());
                m_generator.seed(seq);
            }
            RNG m_generator;
        };

    private:
        unsigned int m_dim;
        matrix m_data;
        matrix m_rotation;
        matrix m_translation;
        double m_rmsd;
        double m_tol;
        size_t m_shuffles;
        kd_tree_t* m_tree;
};



}
}
