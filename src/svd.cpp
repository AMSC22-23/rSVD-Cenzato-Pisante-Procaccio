#include "svd.hpp"

std::tuple<Matrix, Vector, Matrix> SVD::svd_with_PM(Matrix A)
{
    int n = A.cols(), m = A.rows(), i = 0;
    int k = (m > n) ? n : m;
    Matrix B(n, n), U(m, k), V(n, k);
    Vector s(k), u(m), v(n);

    double sigma = 1.;
    while ((sigma > m_epsilon) && i < k)
    {
        B = A.transpose() * A;
        v = PowerMethod(B);

        sigma = (A * v).norm();

        if (sigma > m_epsilon)
        {
            u = A * v * (1 / sigma);

#ifdef EIGEN
            V.col(i) = v;
            U.col(i) = u;
            s[i] = sigma;
#else
            V.col(i, v);
            s(i, 0) = sigma;
            U.col(i, u);
#endif
            A = A - (sigma * u * v.transpose());
            i++;
        }
    }

    // if A not full rank
    if (i < k)
    {
        if (i == 0)
        {
            return std::make_tuple(Matrix::Zero(m, 1), Vector::Zero(1, 1), Matrix::Zero(n, 1));
        }
#ifdef EIGEN
        return std::make_tuple(U.topLeftCorner(m, i), s.head(i), V.topLeftCorner(n, i));
#else
        U.trimCols(k - i);
        V.trimCols(k - i);
        s.trimRows(k - i);
#endif
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

    // if A not full rank
    if (i < k)
    {
        if (i == 0)
        {
            return std::make_tuple(Matrix::Zero(m, 1), Vector::Zero(1, 1), Matrix::Zero(n, 1));
        }
#ifdef EIGEN
        return std::make_tuple(U.topLeftCorner(m, i), s.head(i), V.topLeftCorner(n, i));
#else
        U.trimCols(k - i);
        V.trimCols(k - i);
        s.trimRows(k - i);
#endif
    }

    return std::make_tuple(U, s, V);
}

std::tuple<Matrix, Vector, Matrix> SVD::rsvd(const Matrix &A, int r, int p, int q)
{
    int m = A.rows(), n = A.cols(), k = r + p;
    Matrix U(m, k), P = genmat(n, k), V(k, n), Q(m, m);
    Vector s(k);
    QR_Decomposition QR;

    U = A * P; // m x k
    for (int i = 0; i < q; i++)
    {
        U = A * (A.transpose() * U);
    }

    std::tie(Q, U) = QR.Givens_solve_parallel(U);
    QR.setQR_for_svd_parallel(Q, U);

#ifdef EIGEN
    P.resize(m,k);
    P = Q.topLeftCorner(m, k);
    Q = P;
#else
    Q.trimCols(m - k); // m x k
#endif

    V = Q.transpose() * A; // k x n

    P.resize(k, k);
    std::tie(P, s, V) = svd_with_PM(V);
    U = Q * P;

    return std::make_tuple(U, s, V);
}

Matrix SVD::mult_SVD(Matrix U, Vector s, Matrix V)
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

Matrix SVD::preprocess(Matrix &X)
{
// Mean over row
#pragma omp parallel for
    for (int i = 0; i < X.rows(); i++)
    {
        double media = 0;
        for (int j = 0; j < X.cols(); j++)
        {
            media += X(i, j);
        }
        media /= X.cols();
        for (int j = 0; j < X.cols(); j++)
        {
            X(i, j) /= media;
        }
    }

    // Covariance matrix
    return (1.0 / (X.cols() - 1)) * X.transpose() * X;
}
