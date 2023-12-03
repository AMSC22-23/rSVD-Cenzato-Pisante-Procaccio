#include <iostream>
#include "fullMatrix.hpp"
#include "QR_Decomposition.hpp"

int main(){



    // Assegnazione di valori agli elementi
 // Test matrix
 /*
 Matrix A(4,4);
    A<< 4,3,2,1,
3,4,3,2,
2,3,4,3,
1,2,3,4;
  */
 Matrix A(4,4);
 A<<121.0991,   15.0583,    1.3624,   -0.3454,
   15.0583,   13.4599,   1.7818,  -0.6114,
    1.3624,    1.7818,   0.8663,   -0.3217,
   -0.3454,   -0.6114,   -0.3217,    0.5747;  
/* 
Matrix A(5,3);
A<<0.8147, 0.0975, 0.1576,
0.9058, 0.2785, 0.9706,
0.1270, 0.5469, 0.9572,
0.9134, 0.9575, 0.4854,
0.6324, 0.9649, 0.8003;
*/ 

/*
Matrix A(3,3);
    A<<12, -51, 4,
6, 167, -68,
-4, 24, -41;
*/
/*
    Matrix A(5,5);
    A << 2,-1,0,0,-1.e-12,
        -1,2,-1,0,0,
        0,-1,2,-1,0,
        1.e-12,-1.e-12,-1,2,-1,
        1.e-12,1.e-12,1.e-12,-1,2;
 */   
    QR_Decomposition QR_A;
    
    auto [Q,R]=QR_A.Givens_solve(A);
    
    std::cout<<R<<std::endl;
    std::cout<<Q<<std::endl;


	return 0;
}
