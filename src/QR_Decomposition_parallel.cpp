#include "QR_Decomposition.hpp"


using Matrix=Eigen::MatrixXd;
using Vector=Eigen::VectorXd;


std::tuple<Matrix, Matrix> QR_Decomposition::Givens_solve_parallel(Matrix A){
    
    int m=A.rows();
    int n=A.cols();

    /**
     * Set Q back to the identity, set R equal to A
    */
    Q.resize(m, m);
    Q.setIdentity();

    R=A;
    /**
    * Assembling the matrix Q by applying the Givens rotation at each 
    * iteration and applying each component of Q to R in order to make it triangular
    */
    
    for (int j = 0;j<n;j++){
        for(int i=m-1;i>j;i--){

            /**
                * Givens rotation gets applied by calculating the values of:
                * c: cosine of the angle of rotation
                * s: sine of the angle of rotation
                * 
                * The idea is to calculate the rotations on smaller vectors, iterating on the cells 
                * below the diagonal to pull them to zero
            */
            
            double a=R.coeffRef(i-1,j);
            double b=R.coeffRef(i,j);
            double c,s;
        
            if (abs(a)>abs(b) ){
                if(a!=0.0){
                    int segno = std::signbit(a) ? -1 : 1;
                    c = segno / sqrt(1+(b/a)*(b/a));
                    s  = abs(c/a)*b;
                    } else{
                        c=0.0;
                        s=(a >= 0.0 ? 1.0 : -1.0);
                    }
                }else if (b!=0.0)  {
                    int segno = (std::signbit(b) ? -1 : 1);
                    s = segno / sqrt(1+(a/b)*(a/b));
                    c=abs(s/b)*a;
                    } else{
                        s=0.0;
                        c=(b >= 0.0 ? 1.0 : -1.0);
                }

        
            /**
                * Instead of creating the Givens matrix, I do directly the computation on Q and R
                * In addition I use a temporal variable to avoid changing the matrix before the computation
                * Since there is dependencies within the for loop, I use atomic add
            */
            double tmp = 0.0;
            
                
                #pragma omp sections 
                {
                    #pragma omp section
                    {
                        #pragma omp parallel for shared(R,c,s) num_threads(2)
                        for (int k = 0; k < n; k++) {
                            tmp = c * R.coeffRef(i - 1, k) + s * R.coeffRef(i, k);
                            R.coeffRef(i, k) = -s * R.coeffRef(i - 1, k) + c * R.coeffRef(i, k);
                            R.coeffRef(i - 1, k) = tmp;
                        }
                    }
                
                    
                    #pragma omp section
                    {
                        #pragma omp parallel for shared(Q,c,s) num_threads(2)
                        for (int k = 0; k < m; k++) {
                            tmp = Q.coeffRef(k, i - 1) * c + Q.coeffRef(k, i) * s;
                            Q.coeffRef(k, i) = Q.coeffRef(k, i - 1) * -s + Q.coeffRef(k, i) * c;
                            Q.coeffRef(k, i - 1) = tmp;
                        }
                    }
                    
                }
            R.coeffRef(i, j) = 0.;
        }
    }
    
    return std::make_tuple(Q,R);
}

std::tuple<Matrix, Matrix> QR_Decomposition::QR_parallel(Matrix A){
    
    int m=A.rows();
    int n=A.cols();

    /**
     * Set Q back to the identity, set R equal to A
    */
    Q.resize(m, m);
    Q.setIdentity();

    R=A;
    /**
    * Assembling the matrix Q by applying the Givens rotation at each 
    * iteration and applying each component of Q to R in order to make it triangular
    */
    #pragma omp parallel for num_threads(4) ordered
    for (int j = 0;j<n;j++){
        #pragma omp ordered
        for(int i=m-1+2*omp_get_thread_num();i>j;i--){
            if (i<m){
                std::cout<<"thread "<<omp_get_thread_num()<<"inizia"<<std::endl;
            /**
                * Givens rotation gets applied by calculating the values of:
                * c: cosine of the angle of rotation
                * s: sine of the angle of rotation
                * 
                * The idea is to calculate the rotations on smaller vectors, iterating on the cells 
                * below the diagonal to pull them to zero
            */
            
            double a=R.coeffRef(i-1,j);
            double b=R.coeffRef(i,j);
            double c,s;
        
            if (abs(a)>abs(b) ){
                if(a!=0.0){
                    int segno = std::signbit(a) ? -1 : 1;
                    c = segno / sqrt(1+(b/a)*(b/a));
                    s  = abs(c/a)*b;
                    } else{
                        c=0.0;
                        s=(a >= 0.0 ? 1.0 : -1.0);
                    }
                }else if (b!=0.0)  {
                    int segno = (std::signbit(b) ? -1 : 1);
                    s = segno / sqrt(1+(a/b)*(a/b));
                    c=abs(s/b)*a;
                    } else{
                        s=0.0;
                        c=(b >= 0.0 ? 1.0 : -1.0);
                }

        
            /**
                * Instead of creating the Givens matrix, I do directly the computation on Q and R
                * In addition I use a temporal variable to avoid changing the matrix before the computation
                * Since there is dependencies within the for loop, I use atomic add
            */
            double tmp = 0.0;
            
                
                
            #pragma omp parallel for shared(R,c,s) num_threads(2)
            for (int k = j; k < n; k++) {
                tmp = c * R.coeffRef(i - 1, k) + s * R.coeffRef(i, k);
                R.coeffRef(i, k) = -s * R.coeffRef(i - 1, k) + c * R.coeffRef(i, k);
                R.coeffRef(i - 1, k) = tmp;
            }
        
    
        
            #pragma omp parallel for shared(Q,c,s) num_threads(2)
            for (int k = j; k < m; k++) {
                tmp = Q.coeffRef(k, i -1) * c + Q.coeffRef(k, i) * s;
                Q.coeffRef(k, i) = Q.coeffRef(k, i - 1) * -s + Q.coeffRef(k, i) * c;
                Q.coeffRef(k, i-1) = tmp;
            }
                    
            R.coeffRef(i, j) = 0.;
        }
        }
    }
    
    return std::make_tuple(Q,R);
}    



void setQR_for_svd_parallel(Matrix Q, Matrix R){
int m=Q.rows();
int n=R.rows();
    #pragma omp parallel for num_threads(2)
    for (int i = 0; i < n - 1; i++) {
        if (R.coeffRef(i, i) > 0) {
            #pragma omp parallel for num_threads(2) shared(R)
            for (int j = i; j < n; j++) {
                R.coeffRef(i, j) *= -1;
            }
            #pragma omp parallel for num_threads(2) shared(Q)
            for (int j = 0; j < m; j++) {
                Q.coeffRef(j, i) *= -1;
            }
        }
    }

}


