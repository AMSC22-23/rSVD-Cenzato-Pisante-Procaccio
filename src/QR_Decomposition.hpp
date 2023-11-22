#pragma once

#include <iostream>
#include <Eigen/Dense>
#include <cmath>

class QR_Decomposition{
public:
    /**
     * Build the constructor
    */
    QR_Decomposition(const Eigen::MatrixXd& A_)    :   A(A_)
    {};

    /**
     * Function Givens_solve uses Givens to find the QR factorization 
    */
    void Givens_solve(const Eigen::MatrixXd A);
    /**
     * Function HouseHolder_solve uses HouseHolder to find the QR factorization 
    */
    void HouseHolder_solve(Eigen::MatrixXd A);
    /**
     * Getter functions to matrix A,R,Q
    */
    Eigen::MatrixXd getA() const{
        return A;
    };

    Eigen::MatrixXd getR() const{
        return R;
    };

    Eigen::MatrixXd getQ() const{
        return Q;
    };

/**
     *I call the distructor
    */
~QR_Decomposition() = default;

private:
    const Eigen::MatrixXd A;
    Eigen::MatrixXd Q,R;

};