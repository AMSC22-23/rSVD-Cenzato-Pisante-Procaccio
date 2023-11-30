#include "svd.hpp"

    std::tuple<Matrix, Vector, Matrix> SVD::rsvd(Matrix A, int k){      //da vedere bene!
        int m=A.rows(), n=A.cols();
        Matrix U(m,m), V(n,n), Omega = genmat(n, 2*k), Y(m,2*k), B(m,n);
        Vector s(n);
        QR_Decomposition QR;

        Y = (A * A.transpose()) * A * Omega;        
        auto [Q,R]=QR.Givens_solve(Y);
        B = Q * A;
        std::tie(U,s,V) = svd_with_PM(B);
        U = Q * U;

        return std::make_tuple(U,s,V);
    }


    std::tuple<Matrix, Vector, Matrix> SVD::svd_with_qr(Matrix A){      //non funziona ancora!
        int m=A.rows(), n=A.cols();
        Matrix U = eye(m), V = eye(n);
        Matrix Q1(m,m), Q2(n,n), R2 = A.transpose(), R1;
        Vector s=genvec(n), sold(n);
        QR_Decomposition QR;
        double err=1; 

        /*while(err > m_epsilon){          
            std::tie(Q1,A) = QR.Givens_solve(A);
            std::tie(Q2,R2) = QR.Givens_solve(R2);
            A = A * Q2;
            R2 = R2 * Q1; 
            U = U * Q1;
            V = V * Q2;
            sold=s;
            for(size_t i=0; i<n;i++){
            s[i] = A(i,i);
            }
            err = FrobeniusNorm(sold-s);       //c'è modo migliore?
        }*/

        Matrix B = A.transpose() * A;
        while( err > m_epsilon ){
            std::tie(Q1,R1) = QR.Givens_solve(B);
            B = R1 * Q1;
            U = U * Q1;
            sold=s;
            for(size_t i=0; i<n;i++){
            s[i] = B(i,i);
            }
            err = FrobeniusNorm(sold-s);
        }

        Vector v(n),u(m);
        for(size_t i=0; i<n; i++){
            u = U[i];
            v = A * u / s[i];
            V.col(i)=v;
        }

        return std::make_tuple(U,s,V);
    }



    std::tuple<Matrix, Vector, Matrix> SVD::svd_with_PM(Matrix A){
        int n=A.cols(), m=A.rows(), i=0;
        Matrix B(n,n), U(m,n), V(n,n);
        Vector u(m), v(n), s(n);
        double sigma=1.;
        while((sigma > m_epsilon) && i<n){                                
            B =  A.transpose() * A;
            v = PowerMethod(B);
            sigma = norm(A * v);
            u = A * v / sigma;

            U.col(i) = u;           
            s[i]=sigma;
            V.col(i) = v; 
            
            A = A - sigma * u * v.transpose(); 
            i++;
        }
        return std::make_tuple(U,s,V);
    }


int SVD::compute_rank(Matrix A) {       // costa O(n^3) !!!!!
    int n=A.cols(), m=A.rows();
    int rank = 0, j;
    std::vector<bool> row_selected(n, false);
    for (int i = 0; i < m; ++i) {
        for (j = 0; j < n; ++j) {
            if (!row_selected[j] && abs(A(j,i)) > m_epsilon)
                break;
        }

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