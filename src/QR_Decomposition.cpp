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
        }
   R=A;

   /**
    * Assembling the matrix Q by applying the Givens rotation at each 
    * iteration and applying each component of Q to R in order to make it triangular
   */
    
        for (int j = 0;j<n;j++){
            for(int i=m-1;i>j;i--){
                
                /**
                 * G initialized back to identity
                */
                G=Eigen::MatrixXd::Zero(m, m);
                for(int k=0;k<m;k++)    {   
                  G.coeffRef(k,k)=1;  
                  }
                    
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

                

            }
        }
    }

void QR_Decomposition::HouseHolder_solve(Eigen::MatrixXd A){

    int m=A.rows();
    int n=A.cols();

    /**
     * Initialize v,u
    */
    Eigen::VectorXd u(m),v(m);

    /**
     * Initialize matrix Q (size m x m), matrix R(m x n) and rotation matrix P(m x m)
    */
    Eigen::MatrixXd I(m,m);
    for(int i=0;i<m;i++){
        I.coeffRef(i,i)=1;
    }
    Eigen::MatrixXd Q=I;
    Eigen::MatrixXd P=I;

    Eigen::MatrixXd R=A;

    /**
     * Starting the computation of Q,R
    */
    double mag=0.0;
    double alpha=0.0;
    for(int j=0;j<n;j++){
        /**
         * Initialize u,v to zero at each iteration-i
        */
       u=Eigen::VectorXd::Zero(m);
       v=Eigen::VectorXd::Zero(m);

       /**
        * evaluating each component of the matrix R
       */
       mag=0.0;
      for(int i=j;i<m;i++){
        u.coeffRef(i)=R.coeffRef(i,j);
        mag+=u.coeffRef(i)*u.coeffRef(i);
      }
      mag=sqrt(mag);
      alpha = (u.coeffRef(j) < 0) ? mag : -mag ;

      mag=0.0;
      for(int i=0;i<m;i++){
        v.coeffRef(i)= (j == i) ? (u.coeffRef(i) + alpha) : u.coeffRef(i);
        alpha = (u.coeffRef(j) < 0) ? mag : -mag ;
        mag+=v.coeffRef(i)*v.coeffRef(i);
      }
      mag=sqrt(mag);



      for (int i = j; i < m; i++) v.coeffRef(i) /= mag;
      

    /**
     * Computing P at the j-th iterate and applying the rotation to R,Q
    */
      P=I-2.0*v*v.transpose();
      R=P*R;
      std::cout<<R<<std::endl<<std::endl;
      Q=Q*P;
    }
    
}


