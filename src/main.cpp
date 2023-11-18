#include <iostream>
#include <Eigen/Dense>
#include "fullMatrix.hpp"
#include "QR_Decomposition.hpp"


int main(){

Eigen::MatrixXd A(7,5);

    // Assegnazione di valori agli elementi
    A<< 
	1, 2, 3, 7, 19,
	4, 5, 6, 7, 10,
    9, 10, 13, 12, 18,
	29, 35, 42, 15, 2,
	12, 45, 343, 1000, 2,
	12, 68, 89, 100, 45,
	1, 3, 13, 71, 98;

    
    QR_Decomposition QR_A(A);
    QR_A.Givens_solve(QR_A.getA());
    Eigen::MatrixXd R=QR_A.getR();

    std::cout<<R<<std::endl;





	return 0;
}
