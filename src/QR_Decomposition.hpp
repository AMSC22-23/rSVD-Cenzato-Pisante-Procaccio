#pragma once

#include <iostream>
#include <Eigen/Dense>
#include <cmath>
#include <omp.h>


using Matrix=Eigen::MatrixXd;
using Vector=Eigen::VectorXd;

class QR_Decomposition{
public:
    /**
     * Build the constructor
    */
    QR_Decomposition() {};

    /**
     * Function Givens_solve uses Givens to find the QR factorization 
    */
    std::tuple<Matrix, Matrix> Givens_solve(const Matrix A);
    /**
     * Function HouseHolder_solve uses HouseHolder to find the QR factorization 
    */
    std::tuple<Matrix, Matrix> HouseHolder_solve(const Matrix A);

    /**
     * Parallel method for setting R,Q for svd with QR
    */
    void setQR_for_svd(Matrix Q, Matrix R){};

    /**
     * Parallel method for QR Decomposition
    */
    std::tuple<Matrix, Matrix> Givens_solve_parallel(const Matrix A);

    /**
     * Parallel method for setting R,Q for svd with QR
    */
    void setQR_for_svd_parallel(Matrix Q, Matrix R){};
    /**
     *I call the distructor
    */
~QR_Decomposition() = default;

private:
    Matrix Q,R;

};