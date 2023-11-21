#ifndef SVD_WITH_PM_HPP
#define SVD_WITH_PM_HPP

#include <Eigen/Dense>
#include <iostream>
#include <random>

using Vector=Eigen::VectorXd;
using Matrix=Eigen::MatrixXd;

class svd{
    public:
    /* Constructor 
            inputs : precision (eps) with probability 1-delta 
                     and landa */
        svd(double eps, double delta, double landa) : 
        m_eps(eps), m_delta(delta), m_landa(landa) {}

    /* Power Method */
    void PM(Matrix A,double h,Vector x, Vector &u, Vector &v, double &sigma){
        for(int i=0;i<int(h);i++){
            x = A.transpose() * A * x;
            x = x / x.norm();
        }
        v = x;
        sigma = (A * v).norm();
        u = A * v / sigma;
    }

    /* Generates random vector of n elem with normal distribution */
    Vector genvec(int n){
        Vector v(n);
        std::default_random_engine gen;
        std::normal_distribution<> dice(0);
        for(int i=0;i<n;i++){
            v(i) = dice(gen);
        }
        return v;
    }

    /* SVD with Power Method :
        input: A, U, s, V
        A is the matrix on which the SVD is performed
        s is as a vector containing singular values
        U and V are matrices containing right and left singular vectors*/
    void svd_with_PM(Matrix A, Matrix &U, Vector &s,Matrix &V){
        int n=A.cols();
        int m=A.rows();
        for(int i=0;i<n;i++){                     // n < m 
            Vector x0 = genvec(n);                //initial guess 
            Vector u(m), v(n);
            double sigma;
            double h = log(4*log(2*n/m_delta)/(m_eps*m_delta))/(2*m_landa);
            //double h=5;
            std::cout<<h<<std::endl;
            PM(A,h,x0,u,v,sigma);

            U.col(i) = u;           
            s[i]=sigma;
            V.col(i) = v; 
            
            A = A - sigma * u * v.transpose(); 
            std::cout<<h<<std::endl;
            std::cout<<U<<std::endl;
        }
    }

    /* Destructor */
        ~svd() = default;

    private:
        const double m_eps;
        const double m_delta;
        const double m_landa;
};

#endif // SVD_WITH_PM_HPP