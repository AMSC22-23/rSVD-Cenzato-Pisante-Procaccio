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
            U (m x m) : right singular vectors
            s (n)     : vector containing singular values
            V (n x n) : left singular vectors
        Initialization T0 = A and S0 = AT
        For k = 1, 2, · · ·(repeat until convergence)
            Tk−1 = Uk Rk, Sk−1 = Vk Zk (QR Factorization)
            Tk = Rk Vk and Sk = Zk Uk */
    std::tuple<Matrix, Vector, Matrix> svd_with_qr(Matrix A);


    /* reduced SVD using the Power Method :
        Input:
            A (m x n) : matrix
        Outputs:
            U (m x n) : right singular vectors
            s (n)     : vector containing singular values
            V (n x n) : left singular vectors */
    std::tuple<Matrix, Vector, Matrix> svd_with_PM(Matrix A);


    /* Compute the pseudo-inverse of a matrix A (m x n) using SVD */
    Matrix pseudoinverse(const Matrix A){
        int m = A.rows(), n = A.cols();
        auto[U,s,V]=svd_with_PM(A); //it's the reduced SVD!
        //inverse of sigma
        Matrix S_inv = Matrix::Zero(n,n);
        for(size_t i=0; i<n; i++){
            S_inv(i,i) = 1 / s[i];
        }
        return V * S_inv * U.transpose();
    }

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


    /* Frobenius norm : ||A|| = sqrt(sigma1^2 + ... + sigmar^2),
                    where A is diagonal with singular values 
        In this function the input v contains the singular values of A*/
    double FrobeniusNorm(const Vector v){
        double aux = 0;
        for(size_t i=0;i<v.size();i++){
            aux += v[i] * v[i];
        }
        return sqrt(aux);
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