#ifndef SVD_WITH_PM_HPP
#define SVD_WITH_PM_HPP

#include <Eigen/Dense>
#include <iostream>
#include <random>

using Vector=Eigen::VectorXd;
using Matrix=Eigen::MatrixXd;

class SVD{
    public:
    /* Constructor : 
            epsilon : precision of power method  */
        SVD(double epsilon) : 
        m_epsilon(epsilon) {}

    /* Generates random vector of n elem with normal distribution */
    Vector genvec(int n){
        Vector v(n);
        //std::default_random_engine gen;
        std::random_device rd ;
        std::knuth_b reng{rd ()};
        std::normal_distribution<> dice(0);
        for(int i=0;i<n;i++){
            v(i) = dice(reng);
        }
        return v;
    }

    /* SVD with Power Method :
            A (m x n) : the matrix on which the SVD is performed
            s (n)     : vector containing singular values
            U (m x m) : right singular vectors
            V (n x n) : left singular vectors*/
    void svd_with_PM(Matrix A, Matrix &U, Vector &s, Matrix &V){
        int n=A.cols(), m=A.rows(), i=0;
        Matrix B(n,n);
        Vector u(m), v(n), x0(n);
        double sigma=1.;
        //double h = log(4*log(2*n/m_delta)/(m_epsilon*m_delta))/(2*m_landa);
        while((abs(sigma) > m_epsilon) && i<n){                      
            x0 = genvec(n);            
            //std::cout<<x0<<"\n\n";
            B =  A.transpose() * A;
            v = PowerMethod(B, x0);
            sigma = norm(A * v);
            u = A * v / sigma;

            U.col(i) = u;           
            s[i]=sigma;
            V.col(i) = v; 
            
            A = A - sigma * u * v.transpose(); 
            std::cout<<A<<std::endl<<std::endl;
            i++;
        }
    }

    /* Destructor */
        ~SVD() = default;


    private:

    /* Power Method :
            B (n x n) = At * A
            x (n) = initial guess
            Output :
            v (n) : x converges to the left singular vector of A  */
    Vector PowerMethod(Matrix B, Vector x){
        Vector xold=Vector::Ones(x.size());
        while( norm(xold-x) > m_epsilon ){
            xold=x;
            x = B * x;
            x = x / norm(x);
        }
        return x;
    }

    double norm(Vector v){
        double norm=0;
        for(int i=0;i<v.size();i++){
            norm += v[i] * v[i];
        }
        return sqrt(norm);
    }

    const double m_epsilon;
};

#endif // SVD_WITH_PM_HPP