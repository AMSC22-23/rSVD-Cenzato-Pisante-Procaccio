
#include "svd.hpp"

//g++ -I${mkEigenInc} svd_test.cpp svd.cpp QR_Decomposition.cpp -o prova

Matrix eye(const int n){
        Matrix A=Matrix::Zero(n,n);
        for(int i=0; i<n; i++){
            A(i,i)=1;
        }
        return A;
}

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

    /*Matrix A(5,3);
    A <<0.8147, 0.0975, 0.1576,
        0.9058, 0.2785, 0.9706,
        0.1270, 0.5469, 0.9572,
        0.9134, 0.9575, 0.4854,
        0.6324, 0.9649, 0.8003;*/

    Matrix A(4,4);
    A<< 4, 5, 6, 7,
        5, 4, 5, 6,
        6, 5, 4, 5,
        7, 6, 5, 4;

    int m=A.rows(), n=A.cols();

    Matrix U = eye(n),V(n,n);
    Vector s(n);
    QR_Decomposition obj_qr;
    int nmax = 20;

    //std::cout<<A<<std::endl<<std::endl;
    //QR algorithm to find eigenvalues of B
    for(size_t i=0;i<nmax;i++){
        auto [Q,R] = obj_qr.Givens_solve(A);
        A = R * Q;
        std::cout<<R<<std::endl<<std::endl;
        std::cout<<Q<<std::endl<<std::endl;
        std::cout<<A<<std::endl<<std::endl;
        U = U * Q;
    }

    //std::cout<<U*A*U.transpose()<<std::endl;


   /*Matrix U(m,n), V(n,n), S=Matrix::Zero(n,n);
    Vector s(n);
    double m_eps=1e-12;        

    SVD obj(m_eps);*/

    /*std::tie(U,s,V) = obj.svd_with_PM(A);
    for(int i=0;i<n;i++){
        S(i,i)=s[i];
    }
    std::cout<<"\nSVD with PM:\n";
    std::cout<<U<<std::endl;
    std::cout<<s<<std::endl;
    std::cout<<V<<std::endl;
    std::cout<<"\nA = U * S * Vt :"<<std::endl;
    std::cout<<U*S*V.transpose()<<std::endl;*/

    /*std::cout<<"\nPseudo-inverse :\n";
    std::cout<<obj.pseudoinverse(A)<<std::endl;*/

    //std::tie(U,s,V) = obj.svd_with_qr(A);
    /*
    for(int i=0;i<n;i++){
        S(i,i)=s[i];
    }
    std::cout<<"\nSVD with QR :\n";
    std::cout<<U<<std::endl;
    std::cout<<s<<std::endl;
    std::cout<<V<<std::endl;
    std::cout<<"\nA = U * S * Vt :"<<std::endl;
    std::cout<<U*S*V.transpose()<<std::endl;*/

    return 0;
}