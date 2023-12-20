#include "svd.hpp"

std::tuple<Matrix, Vector, Matrix> SVD::svd_with_PM(Matrix A)
{
    int n = A.cols(), m = A.rows(), i = 0;
    int k = (m > n) ? n : m;
    Matrix B(n, n), U(m, k), V(n, k);
    Vector u(m), v(n), s(k);
    double sigma = 1.;

    while ((sigma > m_epsilon) && i < k)
    {
        B = A.transpose() * A;
        v = PowerMethod(B);
        sigma = (A * v).norm();
        u = A * v * (1 / sigma);

#ifdef EIGEN
        V.col(i) = v;
        U.col(i) = u;
        s[i] = sigma;
#else
        V.col(i, v);
        U.col(i, u);
        s(i, 0) = sigma;
#endif

        A = A - (sigma * u * v.transpose());
        i++;
    }
    return std::make_tuple(U, s, V);
}

std::tuple<Matrix, Vector, Matrix> SVD::svd_with_PM2(Matrix A)
{
    int n = A.cols(), m = A.rows(), i = 0;
    int k = (m > n) ? n : m;
    Matrix U(m, k), V(n, k);
    Vector s = Vector::Zero(k, 1);
    double sigma = 1.;

    while ((sigma > m_epsilon) && i < k)
    {
        auto [u, v] = PowerMethod2(A);

#ifdef EIGEN
        s[i] = (A * v).norm();
        V.col(i) = v;
        U.col(i) = u;
        A = A - (s[i] * u * v.transpose());
#else
        s(i, 0) = (A * v).norm();
        V.col(i, v);
        U.col(i, u);
        A = A - (s(i, 0) * u * v.transpose());
#endif

        i++;
    }
    return std::make_tuple(U, s, V);
}

std::tuple<Matrix, Vector, Matrix> SVD::rsvd(Matrix A, int r, int p, int q)
{
    int m = A.rows(), n = A.cols(), k = r + p;
    Matrix Z(m, k), P = genmat(n, k), Y(k, n), U(m, k);
    QR_Decomposition QR;

    Z = A * P; // m x k
    for (int i = 0; i < q; i++)
    {
        Z = A * (A.transpose() * Z);
    }
    
    auto [Q, R] = QR.Givens_solve_parallel(Z);
    QR.setQR_for_svd_parallel(Q, R);

#ifdef EIGEN
    Q = Q.topLeftCorner(m, k);
#else
    Q.trimCols(m - k); // m x k
#endif

    Y = Q.transpose() * A; // k x n

    auto [Uy, s, V] = svd_with_PM(Y);
    U = Q * Uy;

    return std::make_tuple(U, s, V);
}
/*
std::tuple<Matrix, Vector, Matrix> SVD::svd_with_qr(const Matrix &A)
{
    int m = A.rows(), n = A.cols();
    int min = (m > n) ? n : m;
    Matrix V(n, min), U(m, min), S;
    Vector s(min);
    QR_Decomposition obj_qr;
    int nmax = 40, i = 0;
    double err = 1., epsilon = 1e-6;

    // block SVD algorithm
    V.setIdentity();
    while (err > epsilon && i < nmax)
    {
        auto [Q, R] = obj_qr.Givens_solve_parallel(A * V);
        U = Q.topLeftCorner(m, min);

        auto [Q2, R2] = obj_qr.Givens_solve_parallel(A.transpose() * U);
        V = Q2.topLeftCorner(n, min);

#pragma omp simd
        for (int i = 0; i < min; i++)
        {
#ifdef EIGEN
            s[i] = std::sqrt(std::abs(R2(i, i)));
#else
            s(i, 0) = std::sqrt(std::abs(R2(i, i)));
#endif
        }

        err = (A - mult_parallel(U, s, V)).norm();
        if (i < 4)
            std::cout << i << " - " << err << std::endl;
        i++;
    }

    return std::make_tuple(U, s, V);
}*/

Matrix SVD::mult_parallel(Matrix U, Vector s, Matrix V)
{
    int m = U.rows(), n = V.rows(), k = s.rows();
    Matrix A = Matrix::Zero(m, n);

#pragma omp parallel for collapse(2) shared(A, s, U, V)
    for (int r = 0; r < m; r++)
    {
        for (int c = 0; c < n; c++)
        {
            double temp = 0.0;
#pragma omp parallel for reduction(+ : temp)
            for (int i = 0; i < k; i++)
            {
#ifdef EIGEN
                temp += s[i] * U(r, i) * V(c, i);
#else
                temp += s(i, 0) * U(r, i) * V(c, i);
#endif
            }
            A(r, c) = temp;
        }
    }
    return A;
}