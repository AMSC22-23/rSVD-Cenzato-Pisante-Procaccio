#ifndef SVD_WITH_PM_HPP
#define SVD_WITH_PM_HPP

#include <Eigen/Dense>
#include <iostream>
#include <random>

#include "QR_Decomposition.hpp"

using Vector = Eigen::VectorXd;
using Matrix = Eigen::MatrixXd;

class SVD{
    public:
    /* Constructor : 
            epsilon : precision  */
        SVD(double epsilon) : 
        m_epsilon(epsilon) {}


    /*Given an m × n matrix A, a target number k of singular vectors, and an exponent q 
    (say q = 1 or q = 2), this procedure computes an approximate rank-2k factorization UΣV∗, 
    where U and V are orthonormal, and Σ is nonnegative and diagonal.
        Stage A:
        1. Generate an n × 2k Gaussian test matrix Ω.
        2. Form Y=(AA*)^qAΩ by multiplying alternately with A and A*
        3. Construct a matrix Q whose columns form an orthonormal basis for the range of Y.
        Stage B:
        4. Form B = Q∗A
        5. Compute an SVD of the small matrix: B = U_hat ΣV*
        6. Set U = Q U_hat */
    std::tuple<Matrix, Vector, Matrix> rsvd(Matrix A, int k);


    /* SVD using QR factorization :
        Input:
            A (m x n) : matrix
        Outputs:
            U (m x m) : matrix whose coloumns are left singular vectors of A 
                        [eigenvectors of A*At]
            s (n)     : vector the containing singular values of A
            V (n x n) : matrix whose coloumns are right singular vectors of A 
                        [eigenvectors of At*A] */
    std::tuple<Matrix, Vector, Matrix> svd_with_qr(Matrix A);


    /* Reduced SVD using the Power Method :
        Input:
            A (m x n) : matrix
        Outputs:
            U (m x m) : matrix whose coloumns are left singular vectors of A 
                        [eigenvectors of A*At]
            s (n)     : vector containing the singular values of A
            V (n x n) : matrix whose coloumns are right singular vectors of A 
                        [eigenvectors of At*A] */
    std::tuple<Matrix, Vector, Matrix> svd_with_PM(Matrix A);


    /* Computes the pseudo-inverse of a matrix A (m x n) using SVD */
    Matrix pseudoinverse(const Matrix A){
        int m = A.rows(), n = A.cols();
        auto[U,s,V]=svd_with_PM(A); 
        
        Matrix S_inv = Matrix::Zero(n,n);
        for(size_t i=0; i<n; i++){
            S_inv(i,i) = 1 / s[i];
        }
        return V * S_inv * U.transpose();
    }

    /* Computes the rank of the matrix A */
    int compute_rank(Matrix A);


    /* Destructor */
        ~SVD() = default;


    private:

    /* Power Method :
        Input :
            B (n x n) = At * A
        Output :
            v (n) : x converges to the left singular vector of A
            [x is a vector generated randomically with normal distribution]  */
    Vector PowerMethod(const Matrix B){
        Vector x = genvec(B.cols());        //initial guess
        Vector xold(x.size());
        double err = 1.;
        x /= norm(x);
        while( err > m_epsilon ){
            xold=x;
            x = B * x;
            x = x / norm(x);
            err = norm(xold-x);
        }
        return x;
    }


    /* Generates random vector of n elem with normal distribution */
    Vector genvec(const int n){
        Vector v(n);
        //std::default_random_engine gen;
        std::random_device rd ;
        std::knuth_b reng{rd ()};
        std::normal_distribution<> dice(0.,1.);
        for(int i=0;i<n;i++){
            v[i] = dice(reng);
        }
        return v;
    }


    /*Generates m x n Gaussian matrix M*/
    Matrix genmat(const int m, const int n){
        Matrix M(m,n);
        std::random_device rd ;
        std::knuth_b reng{rd ()};
        std::normal_distribution<> dice(0.,1.);
        for(int i=0;i<m;i++){
            for(int j=0;j<n;j++)
                M(i,j) = dice(reng);
        }
        return M;
    }


    /* Returns the norm of a vector*/
    double norm(const Vector v){
        double norm=0;
        for(int i=0;i<v.size();i++){
            norm += v[i] * v[i];
        }
        return sqrt(norm);
    }


    /* Returns an identity matrix (n x n) */
    Matrix eye(const int n){
        Matrix A=Matrix::Zero(n,n);
        for(int i=0; i<n; i++){
            A(i,i)=1;
        }
        return A;
    }

    const double m_epsilon;
};

#endif // SVD_WITH_PM_HPP