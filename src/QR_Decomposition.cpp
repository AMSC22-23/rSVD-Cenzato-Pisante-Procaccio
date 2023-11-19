#include "QR_Decomposition.hpp"


void QR_Decomposition::Givens_solve(Eigen::MatrixXd A){

    int m=A.rows();
    int n=A.cols();
    
    /**
     * Initialize matrix Q (size m x m), matrix R(m x n) and matrix of rotations G(m x m)
    */
   Q.resize(m,m);
   Eigen::MatrixXd G(m,m);
        for(int i=0;i<m;i++){
            Q.coeffRef(i,i)=1;
            G.coeffRef(i,i)=1;
        }
   R=A;

   /**
    * Assembling the matrix Q by applying the Givens rotation at each 
   */
    
        for (int j = 0;j<n;j++){
            for(int i=m-1;i>j;i--){
                
                /**
                 * G initialized back to identity
                */
                for(int k=0;k<m;k++)    {   G.coeffRef(k,k)=1;  }
                    
                /**
                 * Givens rotation gets applied by calculating the values of:
                 * c: cosine of the angle of rotation
                 * s: sine of the angle of rotation
                 * r: length of the vector in R^2
                 * 
                 * The idea is to calculate the rotations on smaller vectors, iterating on the cells 
                 * below the diagonal to pull them to zero
                */
                
                double a=R.coeffRef(i-1,j);
                double b=R.coeffRef(i,j);
                double c,s;
                       
                                
                c = a / sqrt(pow(a,2) + pow(b,2));
                s = b / sqrt(pow(a,2) + pow(b,2));
                            

                /**
                 * Apply c,s to the Givens matrix G_i
                */

                G.coeffRef(i-1,i-1)=c;
                G.coeffRef(i-1,i)=s;
                G.coeffRef(i,i-1)=-s;
                G.coeffRef(i,i)=c;

                /**
                 * Compute Q, R
                */
                Q = Q*Eigen::MatrixXd(G.transpose());
                R = G*R;

                /**
                * Forcing rotated component to zero to avoid floating point approximations
                */
                R.coeffRef(i,j)=0;

                G=Eigen::MatrixXd::Zero(m, m);;

            }
        }
    }


