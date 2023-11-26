#include "QR_Decomposition.hpp"

using Matrix=Eigen::MatrixXd;
using Vector=Eigen::VectorXd;


std::tuple<Matrix, Matrix> QR_Decomposition::Givens_solve(Matrix A){
    
    int m=A.rows();
    int n=A.cols();

    /**
        * Initialize matrix Q (size m x m), matrix R(m x n) and matrix of rotations G(m x m)
    */
    Q.resize(m,m);
    //Q.coeffRef(m-1,m-1)=1;
    for(int i=0;i<m;++i) Q.coeffRef(i,i)=1.;

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
                * r: length of the vector in R^2
                * 
                * The idea is to calculate the rotations on smaller vectors, iterating on the cells 
                * below the diagonal to pull them to zero
            */
            
            double a=R.coeffRef(i-1,j);
            double b=R.coeffRef(i,j);
            double c,s;
      
            if (abs(a)>abs(b) ){
                if(a!=0.0){
                    c = 1 / sqrt(1+(b/a)*(b/a));
                    s  = -c*b/a;
                } else{
                    c=0.0;
                    s=1.0;
                }  
            } else if (b!=0.0)  {
                s = -1 / sqrt(1+(a/b)*(a/b));
                c=-s*a/b;
            } else{
                s=1.0;
                c=0.0;
            }

        
            /**
                * Instead of creating the Givens matrix, I do directly the computation on Q and R
                * In addition I use a temporal variable to avoid changing the matrix before the computation
            */
            double tmp=0.0;
            for(int k=0;k<n;k++){
                tmp=c*R.coeffRef(i-1,k)-s*R.coeffRef(i,k);
                R.coeffRef(i,k)=s*R.coeffRef(i-1,k)+c*R.coeffRef(i,k); 
                R.coeffRef(i-1,k)=tmp;
            }
                /**
            * Forcing rotated component to zero to avoid floating point approximations
            */
            
            R.coeffRef(i,j)=0;

            /**
             * Computation of Q, at each iterate
            */
            //Q.coeffRef(i-1,i-1)=1;
            for(int k=0;k<m;k++){
                tmp=Q.coeffRef(k,i-1)*c+Q.coeffRef(k,i)*-s;
                Q.coeffRef(k,i)=Q.coeffRef(k,i-1)*s+Q.coeffRef(k,i)*c; 
                Q.coeffRef(k,i-1)=tmp;
            }
            


        }
    }
    return std::make_tuple(Q,R);
}






    void QR_Decomposition::HouseHolder_solve(Matrix A){

    int m=A.rows();
    int n=A.cols();

    /**
        * Initialize v,u
    */
        Vector u(m),v(m);

    /**
        * Initialize matrix Q (size m x m), matrix R(m x n) and rotation matrix P(m x m)
    */
        Matrix I(m,m);
    for(int i=0;i<m;i++){
        I.coeffRef(i,i)=1;
    }

        Matrix Q=I;
    Matrix P=I;

    Matrix R=A;

    /**
        * Starting the computation of Q,R
    */
    double mag=0.0;
    double alpha=0.0;
    for(int j=0;j<n;j++){
        /**
            * Initialize u,v to zero at each iteration-i
        */


        u=Vector::Zero(m);
        v=Vector::Zero(m);

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

