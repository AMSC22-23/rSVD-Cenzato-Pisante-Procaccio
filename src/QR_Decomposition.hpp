#pragma once

#include <iostream>
#include <Eigen/Dense>
#include <cmath>

using Matrix=Eigen::MatrixXd;

class QR_Decomposition{
public:
    /**
     * Build the constructor
    */
    QR_Decomposition(const Matrix& A_)    :   A(A_)
    {};

    /**
     * Function Givens_solve uses Givens to find the QR factorization 
    */
    void Givens_solve(const Matrix A);
    /**
     * Function HouseHolder_solve uses HouseHolder to find the QR factorization 
    */
    void HouseHolder_solve(Matrix A);
    /**
     * Getter functions to matrix A,R,Q
    */
    const Matrix getA() const{
        return A;
    };

    const Matrix getR() const{
        return R;
    };

    const Matrix getQ() const{
        return Q;
    };

/**
     *I call the distructor
    */
~QR_Decomposition() = default;

private:
    const Matrix A;
    Matrix Q,R;

};