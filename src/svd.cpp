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
        sigma = (A * v).norm();
        u = A * v * (1/sigma);

        #ifdef EIGEN
        V.col(i)=v;
        U.col(i)=u;
        s[i]=sigma;
        #else
        V.col(i,v); 
        U.col(i,u);           
        s(i,1)=sigma;
        #endif

        A = A - (sigma * u * v.transpose()); 
        i++;
    }
    return std::make_tuple(U,s,V);
}


std::tuple<Matrix, Vector, Matrix> SVD::rsvd(Matrix A, int r, int p, int q){  
    int m=A.rows(), n=A.cols(), k = r + p;
    Matrix Z(m,k) , P = genmat(n,k), Y(k,n);
    QR_Decomposition QR;
    exportmatrix(P,"P.txt");

    Z = A * P;                             // m x k                 
    for(int i=0; i<q; i++){
        Z = A * (A.transpose() * Z);
    } 
    exportmatrix(Z,"Z.txt"); 

    auto [Q,R] = QR.Givens_solve(Z);
    QR.setQR_for_svd(Q,R);
    exportmatrix(Q,"Q_before.txt");
    exportmatrix(R,"R.txt");
    #ifdef EIGEN
    Q.topLeftCorner(m,k);
    #else
    Q.trimCols(m-k);                       // m x k
    #endif
    
    exportmatrix(Q,"Q.txt");
    
    Y = Q.transpose() * A;                // k x n

    auto [Uy,s,V] = svd_with_PM(Y);
    auto U = Q * Uy;

    return std::make_tuple(U,s,V);
}


std::tuple<Matrix, Vector, Matrix> SVD::svd_with_qr(Matrix A){   
    int m = A.rows(), n = A.cols();
    int min = (m > n) ? n : m;
    //Matrix B = A.transpose() * A;
    Matrix V(n,n), U(m,m), S = A.transpose(), Q;
    Vector s(n),s_old(n);
    QR_Decomposition obj_qr;
    int nmax = 5 * n, cont = 0;
    double err = 1., epsilon = 1e-10;

    int flag=0;

    V.setIdentity();
    U.setIdentity();
    //QR algorithm to find eigenvalues of A (m x n)
    while( cont < nmax && flag==0 ){
        std::tie(Q,S) = obj_qr.Givens_solve(S.transpose());
        obj_qr.setQR_for_svd(Q,S);
        U = U * Q;
        std::tie(Q,S) = obj_qr.Givens_solve(S.transpose());
        obj_qr.setQR_for_svd(Q,S);
        V = V * Q;

        /*s_old = s;
        for(size_t i=0; i<n; i++){        
            s[i] = sqrt(A(i,i));
        }*/
        //err = (s_old-s).norm()/n;

        flag = 1;
        int k=0;
        while(k<n && flag==1){
            int j=k+1;
            while(flag ==1 && j<n){
                if(k<n-1 && S(k,j) > epsilon) {
                    flag = 0;
                    break;
                }
                j++;
            }
            j=k-1;
            while(flag ==1 && j<n){
                if(k>0 && S(k,j) > epsilon)  {
                    flag = 0;
                    break;
                }
                j++;
            }
            k++;
        }

        cont++;
    }
    std::cout<<cont<<" , err = "<<err<<std::endl;

    for(int i=0; i<min; i++){        
        s(i,1) = sqrt(abs(S(i,i)));//@note again! use std::abs
    }

    /*for(size_t i=0; i<k; i++){ 
        U.col(i) = A * V.col(i) / s[i];
    }*/

    return std::make_tuple(U,s,V);
}