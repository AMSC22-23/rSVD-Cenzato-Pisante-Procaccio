#include "svd.hpp"

//g++ -I${mkEigenInc} svd_test.cpp svd.cpp QR_Decomposition.cpp -o prova

int main(){

    // Test matrices

    /*Matrix A(4,3);
    A << 1,2,3,
        4,5,6,
        7,8,9,
        10,11,12;*/

    /*Matrix A(5,5);
    A << 2,-1,0,0,0,
        -1,2,-1,0,0,
        0,-1,2,-1,0,
        0,0,-1,2,-1,
        0,0,0,-1,2;*/

    Matrix A(5,3);
    A <<0.8147, 0.0975, 0.1576,
        0.9058, 0.2785, 0.9706,
        0.1270, 0.5469, 0.9572,
        0.9134, 0.9575, 0.4854,
        0.6324, 0.9649, 0.8003;

    /*Matrix A(4,4);
    A <<4,3,2,1,
        3,4,3,2,
        2,3,4,3,
        1,2,3,4;*/
 

    int m=A.rows(), n=A.cols();
    Matrix U(m,n), V(n,n);
    Vector s(n);
    double m_eps=1e-13;        
    SVD obj(m_eps);

    std::tie(U,s,V) = obj.svd_with_PM(A);
    std::cout<<"\nReduced SVD with Power Method:\n";
    std::cout<<U<<std::endl;
    std::cout<<s<<std::endl;
    std::cout<<V<<std::endl;
    std::cout<<"\nA = U * S * Vt :"<<std::endl;
    std::cout<<obj.mult(U,s,V)<<std::endl;

    /*std::cout<<"\nPseudo-inverse :\n";
    std::cout<<obj.pseudoinverse(A)<<std::endl;*/

    std::tie(U,s,V) = obj.svd_with_qr(A);
    std::cout<<"\nSVD with QR :\n";
    std::cout<<U<<std::endl;
    std::cout<<s<<std::endl;
    std::cout<<V<<std::endl;
    std::cout<<"\nA = U * S * Vt :"<<std::endl;
    std::cout<<obj.mult(U,s,V)<<std::endl; 

    return 0;
}