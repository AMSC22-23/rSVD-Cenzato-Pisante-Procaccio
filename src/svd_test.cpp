#include <Eigen/Dense>
#include <iostream>

using Vector=Eigen::VectorXd;
using Matrix=Eigen::MatrixXd;

#include "svd_with_pm.hpp"

int main(){

    // Test matrix
    Matrix A(4,3);
    A << 1,2,3,
        4,5,6,
        7,8,9,
        10,11,12;

    /*Matrix A(5,5);
    A << 2,-1,0,0,0,
        -1,2,-1,0,0,
        0,-1,2,-1,0,
        0,0,-1,2,-1,
        0,0,0,-1,2;*/

    int m=A.rows(), n=A.cols();

    // Definitions of matrices
    Matrix U(m,m),V(n,n);
    Vector s(n);
    double m_eps=1e-10, m_delta = 0.1, m_landa = 0.1; //da capire come settare parametri

    svd obj(m_eps, m_delta, m_landa);
    obj.svd_with_PM(A,U,s,V);

    std::cout<<"\nSVD :\n";
    std::cout<<U<<std::endl;
    std::cout<<s<<std::endl;
    std::cout<<V<<std::endl;

    return 0;
}