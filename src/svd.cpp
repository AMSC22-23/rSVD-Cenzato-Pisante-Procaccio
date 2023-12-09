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
    Matrix B = A.transpose() * A, B_old = B;
    Matrix V = eye(n),U(m,m);
    Vector s(n), s_old(n);
    QR_Decomposition obj_qr;
    int nmax = 5 * n, i = 0;
    double err = 1.;
    for(size_t i=0; i<n; i++){            
        s[i] = sqrt(B(i,i));
    }

    //QR algorithm to find eigenvalues of B (n x n)
    while( err > m_epsilon && i < nmax){
        auto [Q,R] = obj_qr.Givens_solve(B);
        B = R * Q;
        V = V * Q; 

        err = norm(B_old-B);
        i++;
    }

    for(size_t i=0; i<n; i++){
        U.col(i) = A * V.col(i) / s[i];
    }

    return std::make_tuple(U,s,V);
}


int SVD::compute_rank(Matrix A) {       // costa O(n^3)
    int n=A.cols(), m=A.rows();
    int rank = 0;
    std::vector<bool> row_selected(n, false);
    for (size_t i = 0; i < m; ++i) {
        int j=0;
        while((row_selected[j] || abs(A(j,i)) < m_epsilon) && j<n)
            j++;

        if (j != n) {
            ++rank;
            row_selected[j] = true;
            for (int p = i + 1; p < m; ++p)
                A(j,p) /= A(j,i);
            for (int k = 0; k < n; ++k) {
                if (k != j && abs(A(k,i)) > m_epsilon) {
                    for (int p = i + 1; p < m; ++p)
                        A(k,p) -= A(j,p) * A(k,i);
                }
            }
        }
    }
    return rank;
}