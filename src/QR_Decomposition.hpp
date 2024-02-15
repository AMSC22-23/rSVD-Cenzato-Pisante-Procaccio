#pragma once

#include "utils.hpp"


class QR_Decomposition{
public:
    /**
     * Build the constructor
    */
   //#note The constructor is not needed here. The synthetic constructor provided by the compiler is enough.
    QR_Decomposition() {};

    /**
     * Function Givens_solve uses Givens to find the QR factorization 
    */
    std::tuple<Matrix, Matrix> Givens_solve(const Matrix A);
    /**
     * Function HouseHolder_solve uses HouseHolder to find the QR factorization 
    */
   //@note Since the types are the same you could have also used std::array<Matrix,2> instead of std::tuple<Matrix,Matrix>
   // just a matter of taste.
    std::tuple<Matrix, Matrix> HouseHolder_solve(const Matrix A);
    std::tuple<Matrix, Matrix> HouseHolder_solve_2(Matrix A);

    /**
     * Parallel method for setting R,Q for svd with QR
    */
    void setQR_for_svd(Matrix Q, Matrix R){};

    /**
     * Parallel method for QR Decomposition
    */
    std::tuple<Matrix, Matrix> Givens_solve_parallel(const Matrix A);
    std::tuple<Matrix, Matrix> HouseHolder_solve_parallel(const Matrix &A);
    std::tuple<Matrix, Matrix> HouseHolder_solve_2_parallel(const Matrix &A);

    /**
     * Parallel method for setting R,Q for svd with QR
    */
    void setQR_for_svd_parallel(Matrix Q, Matrix R){};
    /**
     *I call the distructor
    */
   //@note The destructor is not needed here. The synthetic destructor provided by the compiler is enough.
~QR_Decomposition() = default;

private:
    Matrix Q,R;

};