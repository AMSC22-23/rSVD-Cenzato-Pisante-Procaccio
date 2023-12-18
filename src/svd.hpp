#ifndef SVD_HPP
#define SVD_HPP

#include <random>

#include "utils.hpp"
#include "QR_Decomposition.hpp"

#include <fstream>
#include <sstream> 
#include <iomanip>

class SVD{
    public:
    /* Constructor : 
            epsilon : precision  */
        SVD(double epsilon) : 
        m_epsilon(epsilon) {}
//@note This class is not default constructible. It is better to have defautl constructible classes.
//      Here it is sufficient to give a default value to epsilon in the constructor
//      explicit SVD(double epsilon=1.e-8)
//      It is also better make the constructor explicit to avoid implicit conversions.
    void exportmatrix(const Matrix& A, std::string outputFileName){
        // Write the matrix to the file
        std::ofstream outputFile(outputFileName);
        if (outputFile.is_open()) {
            int rows = A.rows(), cols = A.cols();
            // Write dimensions to the first row
            outputFile << rows << " " << cols << std::endl;

            // Write matrix data
            for (int i = 0; i < rows; ++i) {
                for (int j = 0; j < cols; ++j) {
                    outputFile << std::setw(8) << std::fixed << std::setprecision(4) << A(i,j) << " ";
                }
                outputFile << std::endl;
            }
            std::cout << "Computed matrix has been written to "<< outputFileName << std::endl;

            // Close the file
            outputFile.close();
        } else {
            std::cerr << "Error opening file for writing." << std::endl;
        }

    }


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


    /* Computes the pseudo-inverse of a matrix A (m x n) using SVD */
    Matrix pseudoinverse(const Matrix A){
        int m = A.rows(), n = A.cols();
        int k = (m > n) ? n : m;
        auto[U,s,V]=svd_with_PM(A); 
        
        for(int i=0; i<k; i++)
            s(i,1) = 1 / s(i,1);

        return mult(V,s,U);
    }


    /* This procedure computes an approximate rank-k factorization UΣV∗, 
    where U and V are orthonormal, and Σ is nonnegative and diagonal.
    Inputs:
        A = m x n matrix
        r = target rank
        p = oversampling parameter
        q = exponent of power iteration */
    std::tuple<Matrix, Vector, Matrix> rsvd(Matrix A, int r, int p, int q);  


    /* SVD using QR algorithm :
        Input:
            A (m x n) : matrix
        Outputs:
            U (m x m) : matrix whose coloumns are left singular vectors of A 
                        [eigenvectors of A*At]
            s (n)     : vector the containing singular values of A
            V (n x n) : matrix whose coloumns are right singular vectors of A 
                        [eigenvectors of At*A] */
    std::tuple<Matrix, Vector, Matrix> svd_with_qr(Matrix A);


    /* Multiplication to obtain A (m x n) from the svd,
        Input:
            U = matrix (m x n) -- va bene anche per U (m x m)
            s = vector (n)
            V = matrix (n x n) 
        [Also used to calculate the inverse of A (n x m)]*/
    Matrix mult(Matrix U, Vector s, Matrix V){
        int m = U.rows(), n = V.rows(), k = s.rows();
        //if (m == n) n = V.rows();       // to compute the inverse
        //int k = (m > n) ? n : m;
        Matrix A = Matrix::Zero(m,n);
        //!note with the Eigen you can do A = U * s.asDiagonal() * V.transpose();
        for(int r=0; r<m; r++)
            for(int c=0; c<n; c++)
                for(int i=0; i<k; i++)
                    #ifdef EIGEN
                    A(r,c) += s[i] * U(r,i) * V(c,i);
                    #else
                    A(r,c) += s(i,1) * U(r,i) * V(c,i);
                    #endif
        return A;
    }


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
        Vector x = genmat(B.cols(),1);        //initial guess
        Vector xold(x.rows());
        double err = 1.;
        //@note why not implement the *= operator and do x *= 1 / x.norm() ?
        //      with Eigen it is already possible.
        x = x * (1 / x.norm());
        while( err > m_epsilon ){
            xold = x;
            x = B * x;
            x = x * (1 / x.norm());
            err = (xold-x).norm();
        }
        return x;
    }


    /* Generates m x n Gaussian matrix M */
    Matrix genmat(const int m, const int n){
        Matrix M(m,n);
        //std::default_random_engine gen;
        std::random_device rd ;
        std::knuth_b reng{rd ()};
        std::normal_distribution<> dice(0.,1.);
        for(int i=0;i<m;i++){
            for(int j=0;j<n;j++)
                M(i,j) = dice(reng);
        }
        return M;
    }

    const double m_epsilon;
};

#endif // SVD_HPP