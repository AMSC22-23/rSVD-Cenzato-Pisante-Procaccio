#include "svd_with_pm.hpp"

int main(){

    // Test matrix
    /*Matrix A(4,3);
    A << 1,2,3,
        4,5,6,
        7,8,9,
        10,11,12;*/

    Matrix A(5,5);
    A << 2,-1,0,0,0,
        -1,2,-1,0,0,
        0,-1,2,-1,0,
        0,0,-1,2,-1,
        0,0,0,-1,2;

    int m=A.rows(), n=A.cols();

    Matrix U(m,m),V(n,n),S=Matrix::Zero(m,n);
    Vector s(n);
    double m_eps=1e-12;        

    SVD obj(m_eps);
    obj.svd_with_PM(A,U,s,V);

    std::cout<<"\nSVD :\n";
    std::cout<<U<<std::endl;
    std::cout<<s<<std::endl;
    std::cout<<V<<std::endl;

    for(int i=0;i<A.cols();i++){
        S(i,i)=s[i];
    }

    std::cout<<"\n"<<U*S*V.transpose()<<std::endl;

    return 0;
}