#include <iostream>
#include "fullMatrix.hpp"
#include "QR_Decomposition.hpp"

int main(){


Eigen::MatrixXd A(5,3);
    // Assegnazione di valori agli elementi
     
    A<<0.8147, 0.0975, 0.1576,
0.9058, 0.2785, 0.9706,
0.1270, 0.5469, 0.9572,
0.9134, 0.9575, 0.4854,
0.6324, 0.9649, 0.8003;
    
    QR_Decomposition QR_A;
    
    auto [Q,R]=QR_A.Givens_solve(A);
     
    std::cout<<R<<std::endl;
    std::cout<<Q<<std::endl;


	return 0;
}