#ifndef SVD_HPP
#define SVD_HPP

#include "utils.hpp"
#include "QR_Decomposition.hpp"

class SVD{
    public:
    /* Constructor : 
            epsilon : precision  */
        SVD(const double epsilon = 1e-8) : 
        m_epsilon(epsilon) {}

    /* It computes the reduced SVD using the Power Method :
        Input:
            A (m x n) : matrix
        Outputs:
            U (m x m) : matrix whose coloumns are left singular vectors of A 
                        [eigenvectors of A*At]
            s (n)     : vector containing the singular values of A
            V (n x n) : matrix whose coloumns are right singular vectors of A 
                        [eigenvectors of At*A] */
    std::tuple<Matrix, Vector, Matrix> svd_with_PM(Matrix A);


    /* Another algorithm to compute the SVD using the power method */
    std::tuple<Matrix, Vector, Matrix> svd_with_PM2(Matrix A);


    /* Computes the pseudo-inverse of a matrix A (m x n) using SVD */
    Matrix pseudoinverse(const Matrix A){
        auto[U,s,V]=svd_with_PM(A);
        size_t k =s.rows();        
        for(size_t i=0; i<k; i++){ 
            #ifdef EIGEN
            s[i] /= s[i];
            #else
            s(i,0) = 1 / s(i,0);         
            #endif
        }

        return mult_SVD(V,s,U);
    }


    /* This procedure computes an approximate rank-k factorization UΣV∗, 
    where U and V are orthonormal, and Σ is nonnegative and diagonal.
    Inputs:
        A = m x n matrix
        r = target rank
        p = oversampling parameter (default = 5)
        q = exponent of power iteration (default = 1) */
    std::tuple<Matrix, Vector, Matrix> rsvd(const Matrix &A, int r, int p = 5, int q = 1); 


    /* Multiplication in parallel to obtain A (m x n) from the svd,
        Input:
            U = matrix (m x r) 
            s = vector (r)
            V = matrix (n x r) */
    Matrix mult_SVD(Matrix U, Vector s, Matrix V);

    /* Generates m x n Gaussian matrix */
    Matrix genmat(const int m, const int n);
    

    /* Destructor */
        ~SVD() = default;


    private:

    /* Power Method :
        Input :
            B (n x n) = At * A
        Output :
            v (n) : x converges to the left singular vector of A
            [x is a vector generated randomically with normal distribution]  */
    Vector PowerMethod(const Matrix &B){
        Vector x = genmat(B.cols(),1);        //initial guess
        Vector xold(x.rows());
        double err = 1.;
        x = x * (1 / x.norm());
        while( err > m_epsilon ){
            xold = x;
            x = B * x;
            x = x * (1 / x.norm());
            err = (xold-x).norm();
        }
        return x; 
    }

    std::tuple <Vector, Vector> PowerMethod2 (const Matrix &A){
        Vector u(A.rows()), v = genmat(A.cols(),1), v_old(A.cols());
        double err = 1.;
        v = v * (1 / v.norm());

        while ( err > m_epsilon ){
            v_old = v;
            u = A * v;
            u = u * (1/u.norm());
            v = A.transpose() * u;
            v = v * (1 / v.norm());
            err = (v_old-v).norm();
        }
        return std::make_tuple(u,v);
    }

    const double m_epsilon;
};

#endif // SVD_HPP