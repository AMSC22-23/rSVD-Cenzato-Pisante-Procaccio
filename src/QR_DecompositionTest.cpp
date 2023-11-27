#include <iostream>
#include "fullMatrix.hpp"
#include "QR_Decomposition.hpp"

int main(){



    // Assegnazione di valori agli elementi
 // Test matrix
    /*    
    Matrix A(5,3);
    A<<0.8147, 0.0975, 0.1576,
0.9058, 0.2785, 0.9706,
0.1270, 0.5469, 0.9572,
0.9134, 0.9575, 0.4854,
0.6324, 0.9649, 0.8003;
*/
/*
Matrix A(4,3);
    A << 1,2,3,
        4,5,6,
        7,8,9,
        10,11,12;
*/
    Matrix A(5,5);
    A << 2,-1,0,0,-1.e-12,
        -1,2,-1,0,0,
        0,-1,2,-1,0,
        1.e-12,-1.e-12,-1,2,-1,
        1.e-12,1.e-12,1.e-12,-1,2;
    
    QR_Decomposition QR_A;
    
    auto [Q,R]=QR_A.Givens_solve(A);
     
    std::cout<<R<<std::endl;
    std::cout<<Q<<std::endl;
	std::cout<<Q*R<<std::endl;

    auto [Q1,R1]=QR_A.HouseHolder_solve(A);
    std::cout<<"---------------------"<<std::endl;
    std::cout<<R1<<std::endl;
    std::cout<<Q1<<std::endl;
	std::cout<<Q1*R1<<std::endl;

	return 0;
}
