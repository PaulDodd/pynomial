//
// Matrix.h
//

#ifndef SHARED_FOLDING_MATRIX_h
#define SHARED_FOLDING_MATRIX_h

#include <assert.h>
#include <math.h>
#include <ctime>
#include <complex>
#include <cstdlib>
#include <fstream>
#include <string>
#include <cstring>

#include <vector>
#include <map>
#include <set>
#include <queue>
#include <algorithm>

#include <iostream>
#include <iomanip>
#include <dirent.h>
#include <climits>
#include <memory>
#include <ctype.h>
#include <utility>
#include <chrono>
#include <random>
#include <sys/types.h>
#include <unistd.h>
#include "utils.h"
// #include "eigen3/Eigen/Dense"
// #include "eigen3/Eigen/Sparse"
#include <eigen3/Eigen/Sparse>
#include <eigen3/Eigen/Dense>
// #include "utils.h"
// #include "json_wrapper.h"

namespace EigenMatrixHelper{
using std::vector;
using std::min;
using std::max;
using std::complex;

// Datatypes
typedef unsigned long long  uint64;

// Global #defines
#define SMALL           (1.0e-12)  // Should set this to the computer epsilon
#define SMALL_FLOAT     (1.0e-6)
#define LARGE           (1.0e10)
#define BOUNDED(x, xmin, xmax)          (x >= xmin && x <= xmax)



inline const std::string  RESET() { return  "\033[0m"; }
inline const std::string  BLACK() { return    "\033[30m"; }      /* Black */
inline const std::string  RED() { return      "\033[31m"; }      /* Red */
inline const std::string  GREEN() { return    "\033[32m"; }      /* Green */
inline const std::string  YELLOW() { return   "\033[33m"; }      /* Yellow */
inline const std::string  BLUE() { return     "\033[34m"; }      /* Blue */
inline const std::string  MAGENTA() { return  "\033[35m"; }      /* Magenta */
inline const std::string  CYAN() { return     "\033[36m"; }      /* Cyan */
inline const std::string  WHITE() { return    "\033[37m"; }     /* White */
inline const std::string  BOLDBLACK() { return    "\033[1m\033[30m"; }      /* Bold Black */
inline const std::string  BOLDRED() { return      "\033[1m\033[31m"; }      /* Bold Red */
inline const std::string  BOLDGREEN() { return    "\033[1m\033[32m"; }      /* Bold Green */
inline const std::string  BOLDYELLOW() { return   "\033[1m\033[33m"; }      /* Bold Yellow */
inline const std::string  BOLDBLUE() { return     "\033[1m\033[34m"; }      /* Bold Blue */
inline const std::string  BOLDMAGENTA() { return  "\033[1m\033[35m"; }      /* Bold Magenta */
inline const std::string  BOLDCYAN() { return     "\033[1m\033[36m"; }      /* Bold Cyan */
inline const std::string  BOLDWHITE() { return    "\033[1m\033[37m"; }      /* Bold White */


//#define EIGEN_DEFAULT_DENSE_INDEX_TYPE
typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> matrix;
typedef Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic> ComplexMatrix;
typedef Eigen::Matrix<std::complex<double>, Eigen::Dynamic, 1> ComplexVector;
typedef Eigen::VectorXd Vector;
typedef Eigen::SparseMatrix<double, Eigen::ColMajor> spMatrix;
typedef Eigen::SparseMatrix<std::complex<double>, Eigen::ColMajor> ComplexSparseMatrix;


//inline void Crossf(const Eigen::VectorXd& v1, const Eigen::VectorXd& v2, Eigen::VectorXd& ret)
//{
//    ret(0,0) = v1(1,0)*v2(2,0) - v1(2,0)*v2(1,0);
//    ret(1,0) = v1(2,0)*v2(0,0) - v1(0,0)*v2(2,0);
//    ret(2,0) = v1(0,0)*v2(1,0) - v1(1,0)*v2(0,0);
//}

//inline double VectorDot(CMatrix<double> v1, CMatrix<double> v2)
//{
//    // v1 and v2 are column vectors
//    CMatrix<double> ret(1,1);
//    ret = v1.Transpose()*v2;
//    return ret(0,0);
//}
inline void RemoveRow(Eigen::SparseMatrix<double, Eigen::RowMajor>& mat, int rowToRemove)
{
    int numRows = int(mat.rows())-1;
    int numCols = int(mat.cols());

    if( rowToRemove < numRows )
        mat.middleRows(rowToRemove,numRows-rowToRemove) = mat.middleRows(rowToRemove+1, numRows-rowToRemove).eval();

    mat.conservativeResize(numRows,numCols);
}

inline void RemoveRow(spMatrix& mat, int rowToRemove)
{
    Eigen::SparseMatrix<double, Eigen::RowMajor> rowMatrix(mat);
    RemoveRow(rowMatrix, rowToRemove);
    mat = rowMatrix;
}



inline void RemoveColumn(spMatrix& mat, int colToRemove)
{
    int numRows = int(mat.rows());
    int numCols = int(mat.cols())-1;

    if( colToRemove < numCols )
        mat.middleCols(colToRemove,numCols-colToRemove) = mat.middleCols(colToRemove+1, numCols-colToRemove).eval();

    mat.conservativeResize(numRows,numCols);
}

template <class MVal>
inline MVal DeltaMatrix(const vector<size_t>& ndx, const int& dim)
{
    MVal delta(dim, dim);
    delta.setZero();
    for(size_t i = 0; i < ndx.size(); i++)
            delta.coeffRef((int) ndx[i],(int) ndx[i]) = 1.0;
    return delta;
}

template <class MVal>
inline void DiagonalMatrix(const Vector& diag, MVal& ret)
{
    ret.resize(int(diag.rows()), int(diag.rows()));
    ret.setZero();
    for(int i = 0; i < diag.rows(); i++)
        ret.coeffRef(i, i) = diag(i,0);
}

inline void SparseKeepPositives(spMatrix& mat)
{
    for (int k=0; k<mat.outerSize(); ++k)
        for (Eigen::SparseMatrix<double>::InnerIterator it(mat,k); it; ++it)
            mat.coeffRef(it.row(), it.col()) = max(0.0, it.value());

    mat.prune(SMALL);
}

inline spMatrix Dense2Sparse(const matrix& A)
{
    spMatrix spA(int(A.rows()), int(A.cols()));
    std::cout << "dense 2 sparse" << std::endl;
    for(int n = 0; n < A.rows(); n++)
    {
        for(int m = 0; m < A.cols(); m++)
        {
            if(fabs(A(n,m)) > SMALL)
            {
                spA.coeffRef(n, m) = A(n,m);
            }
        }
    }

    return spA;
}

inline ComplexSparseMatrix Real2Complex(const spMatrix& A)
{
    ComplexSparseMatrix C;
    C.resize(A.rows(), A.cols());
    for(int k = 0; k < A.outerSize(); k++)
    {
        for (Eigen::SparseMatrix<double>::InnerIterator it(A,k); it; ++it)
        {
            C.coeffRef(it.row(), it.col()) = complex<double>(it.value(), 0.0);
        }
    }
    return C;
}

//template<class RealMatrixType, class ComplexMatrixType>
//inline RealMatrixType Complex2Real(const ComplexMatrixType& A)
//{
//    RealMatrixType C;
//    C.resize(A.rows(), A.cols());
//    for(int k = 0; k < A.outerSize(); k++)
//    {
//        for (Eigen::SparseMatrix<double>::InnerIterator it(A,k); it; ++it)
//        {
//            C.coeffRef(it.row(), it.col()) = complex<double>(it.value(), 0.0);
//        }
//    }
//    return C;
//}

///////////////////////////////////////////////////////////////////////
// Optimization Problems
///////////////////////////////////////////////////////////////////////

inline matrix CenterOfMass(const matrix& P)
{
    // Assumes that P = (v**T) if v is a column vector.  or in other notation  P = [x1, y1, z1; ...]
    matrix cm(1, P.cols());
    for(int i =0; i < P.cols(); i++)
        cm(0,i) = P.col(i).sum()/double(P.rows());
    //cout << "cm = \n"<< cm << std::endl;

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
        // std::cout << "Vb = \n" << V << std::endl<< std::endl;
        V.col(V.cols()-1)*= -1.0;
        // std::cout << "Va = \n" << V << std::endl<< std::endl;
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
    P = EigenMatrixHelper::Translate(-EigenMatrixHelper::CenterOfMass(P), P);
    Q = EigenMatrixHelper::Translate(-EigenMatrixHelper::CenterOfMass(Q), Q);
    KabschAlgorithm(P, Q, rotation); // Find the rotation.
    P = (rotation*P.transpose()).transpose();  // Apply the transform

    if(pRotation) // optionally copy the roation.
        *pRotation = rotation;

}

inline double LineSegmentClosestDistance(const matrix& v1, const matrix& v2, const matrix& u1, const matrix& u2)
{
    // Finds the closest point between the line segments V(t) = (v2-v1)t + v1 and U(s) = (u2-u1)*s + u1
    // and returns the distance between those points d(t*, s*).
    // if the line segments are close to paralell (v2-v1)\cdot(u2-u1) then the min{d(0,0), d(1,1), d(0,1), d(1,0), d(0.5, 0.5)} is returned.
    double distance = 0.0, cos_theta = 0.0, s, t;
    double normU, normV, normUV, normABV, normABU; // u.u, v.v, u.v
    Vector u, v, ab, d;
    ab = v1-u1;
    u = u2-u1;
    v = v2-v1;
    normU = u.dot(u);
    normV = v.dot(v);
    normUV = u.dot(v);
    normABV = ab.dot(v);
    normABU = ab.dot(u);
    cos_theta = normUV/sqrt(normU*normV);

    if(fabs(cos_theta-1) < SMALL || fabs(cos_theta+1) < SMALL)  // the lines are parallel
    {
        double dt = LARGE, ds = LARGE, de = LARGE; // I don't really like this. (think of a better way to do this)
        // std::cout << "Parallel lines detected "<< std::endl;
        // check the two end points of the shorter vector against the longer one.
        if(normV < normU) // v must be the longer vector for this to work.
        {
            return LineSegmentClosestDistance(u1, u2, v1, v2);
        }
        t = -v.dot(ab)/normV;
        ab = v1-u2;
        s = -v.dot(ab)/normV;
        // printf("found t = %f and s = %f \n", t, s);
        bool bs = (s >= 0.0 && s <= 1.0);
        bool bt = (t >= 0.0 && t <= 1.0);
        bool bCheckEndPoints = !(bs && bt);

        if(bt){ // t is valid
            d = (v*t + v1)-u1;
            dt = sqrt(d.dot(d));
        }
        if(bs){ // s is valid
            d = (v*s + v1)-u2;
            ds = sqrt(d.dot(d));
        }
        if(bCheckEndPoints){ // check endpoints.
            // printf("Checking end points. \n");
            double d1, d2, d3, d4;
            d = v1-u1; d1 = sqrt(d.dot(d)); // t = 0, s = 0
            d = v1-u2; d2 = sqrt(d.dot(d)); // t = 0, s = 1
            d = v2-u1; d3 = sqrt(d.dot(d)); // t = 1, s = 0
            d = v2-u2; d4 = sqrt(d.dot(d)); // t = 1, s = 1
            de = min(min(d1, d2), min(d3, d4));
        }

        distance = min(de, min(dt, ds));
        assert(distance < LARGE); // if the distance is greater than LARGE then this function works.

        // printf("distance = %f\n", distance);

    }
    else
    {
        s = ((normABV*normUV) - (normABU*normV)) / ((normUV*normUV) - (normU*normV));
        t = (normUV/normV)*s - (normABV/normV);

        // printf("found t = %f and s = %f \n", t, s);
        if(s < 0.0 || s > 1.0 || t < 0.0 || t > 1.0) // check all of the boundary cases.
        {
            // printf("Boundary case detected! \n");
            size_t amin = 0;
            double t0, t1, s0, s1, x;
            vector<double> dvec;
            x = -normABV/normV;
            t0 = (BOUNDED( x, 0.0, 1.0) ? x : (x < 0 ? 0 : 1.0));    d = (v*t0 + v1)-u1;         dvec.push_back(sqrt(d.dot(d)));
            x = (normUV-normABV)/normV;
            t1 = (BOUNDED( x, 0.0, 1.0) ? x : (x < 0 ? 0 : 1.0));    d = (v*t1 + v1)-(u + u1);   dvec.push_back(sqrt(d.dot(d)));
            x = normABU/normU;
            s0 = (BOUNDED( x, 0.0, 1.0) ? x : (x < 0 ? 0 : 1.0));    d = (v1)-(u*s0 + u1);       dvec.push_back(sqrt(d.dot(d)));
            x = (normUV+normABU)/normU;
            s1 = (BOUNDED( x, 0.0, 1.0) ? x : (x < 0 ? 0 : 1.0));    d = (v + v1)-(u*s1 + u1);   dvec.push_back(sqrt(d.dot(d)));
            amin = utils::argmin(dvec);
            if(amin == 0)           {t = t0; s = 0;}
            else if(amin == 1)      {t = t1; s = 1;}
            else if(amin == 2)      {t = 0; s = s0;}
            else if(amin == 3)      {t = 1; s = s1;}
            else                    {std::cout << "ERROR!" << std::endl;}
            assert(!(s < 0.0 || s > 1.0 || t < 0.0 || t > 1.0));

        }
        d = (v*t + v1)-(u*s + u1);
        distance = sqrt(d.dot(d));
        // printf("Shortest Distance d(t = %f, s = %f)= %f \n", t, s, distance);
    }

    return distance;
}

template <class MVal>
inline bool IsPermuationMatrix(const MVal& P)
{
    bool bPerm = true;
    for(int n = 0; n < P.rows() && bPerm; n++)
    {
        if (!(fabs(P.row(n).toDense().maxCoeff()-1.0) < SMALL && fabs(P.row(n).toDense().minCoeff()) < SMALL && fabs(P.row(n).sum() - 1.0) < SMALL))
        {
            bPerm = false;
        }
    }

    for(int m = 0; m < P.cols() && bPerm; m++)
    {
        if (!(fabs(P.col(m).toDense().maxCoeff()-1.0) < SMALL && fabs(P.col(m).toDense().minCoeff()) < SMALL && fabs(P.col(m).sum() - 1.0) < SMALL))
        {
            bPerm = false;
        }
    }

    return bPerm;
}


///////////////////////////////////////////////////////////////////////
// Eigenvalue and eigenvector approximation
///////////////////////////////////////////////////////////////////////

template<class MatrixType, class VectorType>
inline void PowerMethod(const MatrixType& mat, VectorType& EigenVec, double& EigenVal, size_t max_iter = 1000, double tol = SMALL)
{
    // Finds the dominant Eigen value and vector.  However it has been shown that depending on the initial guess you
    // can converge to a different eigen value and vector pair so dont be blind.

    matrix x0 = matrix::Random(mat.rows(), 1); // Could replace with ones vector.
    matrix x1(mat.rows(), 1),temp(mat.rows(),1), err(mat.rows(),1);

    bool bConverged = false;
    size_t k = 0;
    double max0, maxe, max1, lambda = 0;
    double m0, m1, m2, q, eps;

    while ( k++ < max_iter && !bConverged)
    {
        x1 = mat*x0;

        // Calculate the current Rayleigh Quotient
        m0 = (x0.transpose()*x0)(0,0);
        m1 = (x0.transpose()*x1)(0,0);
        m2 = (x1.transpose()*x1)(0,0);

        q = (m1/m0);
        eps = sqrt(((m2/m0) - (q*q)));

        // now the rest is the power method.
        max0 = x0.maxCoeff();
        max1 = x1.maxCoeff();
        err = x0 - x1;
        maxe = err.norm();
        x0 = x1/max1;
        lambda = max1/max0;   // could also so do norm1/norm0
//        if( m0 < SMALL || max0 < SMALL || max1 < SMALL)
//        {
//            std::cout << mat << std::endl;
//            std::cout << x0 << std::endl;
//            std::cout << x1 << std::endl;
//            std::cout << "*** WARNING *** division by a small number detected 1/"<<m0 << std::endl;
//            std::cout << "*** WARNING *** division by a small number detected 1/"<<max0 << std::endl;
//            std::cout << "*** WARNING *** division by a small number detected 1/"<<max1 << std::endl;
//            throw std::runtime_error("division by a small number detected");
//        }
        if (fabs(maxe) < tol) {
            bConverged = true;
        }
    }

    EigenVal = lambda;
    EigenVec = x1/x1.norm();  // Return a unit vector.

    // std::cout << "v("<<EigenVal<<") = \n"<< EigenVec << std::endl;
}

///////////////////////////////////////////////////////////////////////
// JSON Objects for file I/O
///////////////////////////////////////////////////////////////////////
/*
#if __cplusplus >= 201103L
#include <tuple>
typedef tuple<int, int, double> SparseMatrixElement;

#else
typedef vector< double > SparseMatrixElement; // just store the indices as double... [row, col, val]
#endif // C++11

using json::CJSONValueObject;
class CJSONValueSparseMatrix : public json::CJSONValueObject<CJSONValueSparseMatrix>
{
    public:
        CJSONValueSparseMatrix(const string& name, EigenMatrixHelper::spMatrix* pval) : CJSONValueObject<CJSONValueSparseMatrix>(name, this), m_Value(pval) {}
#ifdef c_plus_plus_11
        CJSONValueSparseMatrix(const CJSONValueSparseMatrix& src) = delete;
#else
    private:
        CJSONValueSparseMatrix(const CJSONValueSparseMatrix& src);
    public:
#endif
        ~CJSONValueSparseMatrix(){}

        void Destroy() {json::CJSONValueObject<CJSONValueSparseMatrix>::Destroy();}
    // Overloaded Function
        void SetupJSONObject()
        {
            // std::cout << "Calling CJSONValueSparseMatrix::SetupJSONObject" << std::endl;
            AddIntegerValue("rows", &m_rows);
            AddIntegerValue("cols", &m_cols);

#if __cplusplus >= 201103L
            AddNameValuePair<vector< SparseMatrixElement >, json::CJSONValueArray<  SparseMatrixElement,
                                                                                    json::CJSONValueTuple<  SparseMatrixElement,
                                                                                                            json::CJSONValueInt,
                                                                                                            json::CJSONValueInt,
                                                                                                            json::CJSONValueDouble > > >("matrix", &m_data);
#else
            AddNameValuePair<vector< SparseMatrixElement >, json::CJSONValueArray<  SparseMatrixElement,
                                                                                    json::CJSONValueArray<  double,
                                                                                                            json::CJSONValueDouble > > >("matrix", &m_data);
#endif // C++11

        }

        bool Parse (const json_t* pVal)
        {
            // std::cout << "Calling CJSONValueSparseMatrix::Parse" << std::endl;
            ClearBuffer();
            bool bParseSuccess = CJSONValueObject<CJSONValueSparseMatrix>::Parse(pVal);

            m_Value->resize(m_rows, m_cols);
            CopyFromBuffer();
            ClearBuffer();
            return bParseSuccess;
        }

        bool Dump (json_t*& pRet)
        {
            // std::cout << "Calling CJSONValueSparseMatrix::Dump" << std::endl;
            ClearBuffer();
            m_rows = m_Value->rows();
            m_cols = m_Value->cols();
            CopyToBuffer();
            bool bDumpSuccess = CJSONValueObject<CJSONValueSparseMatrix>::Dump(pRet);
            ClearBuffer();
            return bDumpSuccess;
        }

    // Copy Operations
        void CopyFromBuffer()
        {
            int row, col;
            double val;
            for(size_t i = 0; i < m_data.size(); i++)
            {
#if __cplusplus >= 201103L
                row = std::get<0>(m_data[i]);
                col = std::get<1>(m_data[i]);
                val = std::get<2>(m_data[i]);
#else
                row = int(m_data[i][0]);
                col = int(m_data[i][1]);
                val = m_data[i][2];
#endif // C++11

                (*m_Value).coeffRef(row, col) = val;

            }
        }

        void CopyToBuffer()
        {
            ClearBuffer();
            for (int k=0; k<m_Value->outerSize(); ++k)
            {
                for (Eigen::SparseMatrix<double>::InnerIterator it(*m_Value,k); it; ++it)
                {
#if __cplusplus >= 201103L
                    //cout << "Making tuple: (" << it.row() << ", "<< it.col() << ", " << it.value() << ")" << std::endl;
                    m_data.push_back(std::make_tuple(it.row(), it.col(), it.value()));
#else
                    vector< double > temp;
                    temp.push_back(double(it.row()));
                    temp.push_back(double(it.col()));
                    temp.push_back(double(it.value()));
                    m_data.push_back(temp);
#endif // C++11
                }
            }
        }

        void ClearBuffer()
        {
            m_data.clear();
        }

    private:
        int                                 m_rows;
        int                                 m_cols;     // number of
        vector< SparseMatrixElement >       m_data;     // will read and write in and out of this buffer
        EigenMatrixHelper::spMatrix*        m_Value;
};


*/

///////////////////////////////////////////////////////////////////////
// end EigenMatrixHelper
///////////////////////////////////////////////////////////////////////
}





#endif
