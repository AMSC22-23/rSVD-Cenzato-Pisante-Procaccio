#include "svd.hpp"

std::tuple<Matrix, Vector, Matrix> SVD::svd_with_PM(Matrix A){
    int n=A.cols(), m=A.rows(), i=0;
    int k = (m > n) ? n : m;
    Matrix B(n,n), U(m,k), V(n,k);
    Vector u(m), v(n), s(k);
    double sigma=1.;
    while((sigma > m_epsilon) && i<k){                                
        B =  A.transpose() * A;
        v = PowerMethod(B);
        sigma = norm(A * v); 
        u = A * v / sigma;

        V.col(i) = v; 
        U.col(i) = u;           
        s[i]=sigma;
            
        A = A - sigma * u * v.transpose(); 
        i++;
    }
    return std::make_tuple(U,s,V);
}


std::tuple<Matrix, Vector, Matrix> SVD::rsvd(Matrix A, int r, int p, int q){  
    int m=A.rows(), n=A.cols(), k = r + p;
    Matrix Z(m,k) ,Omega = genmat(n,k);
    QR_Decomposition QR;

    Z = A * Omega;                     
    for(size_t i=0; i<q; i++){
        Z = A * (A.transpose() * Z);
    } 
    auto [Q,R] = QR.Givens_solve(Z);

    auto Y = Q.transpose() * A;  

    auto [U,s,V] = svd_with_PM(Y);
    U = Q * U;

    return std::make_tuple(U,s,V);
}


std::tuple<Matrix, Vector, Matrix> SVD::svd_with_qr(Matrix A){   
    int m = A.rows(), n = A.cols();
    int k = (m > n) ? n : m;
    Matrix B = A.transpose() * A, R_err;
    Matrix V = eye(n),U(m,m), I=eye(n);
    Vector s(n),s_old(n);
    QR_Decomposition obj_qr;
    int nmax = 5 * n, i = 0;
    double err = 1., epsilon = 1e-3;

    //QR algorithm to find eigenvalues of B (n x n)
    while( err > epsilon && i < 19){
        auto [Q,R] = obj_qr.Givens_solve(B);
        obj_qr.setQR_for_svd(Q,R);
        B = R * Q;
        V = V * Q;

        for(size_t i=0; i<n; i++){        
            s[i] = sqrt(B(i,i));
        }

        err = (Q-I).norm()/(n*n);
        i++;
        if(i==18) exportmatrix(R,"R.txt");
    }
    std::cout<<i<<" , err = "<<err<<std::endl;

    for(size_t i=0; i<k; i++){ 
        U.col(i) = A * V.col(i) / s[i];
    }

    return std::make_tuple(U,s,V);
}