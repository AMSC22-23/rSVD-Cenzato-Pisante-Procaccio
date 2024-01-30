#include "QR_Decomposition.hpp"

std::tuple<Matrix, Matrix> QR_Decomposition::Givens_solve_parallel(const Matrix A)
{

    int m = A.rows();
    int n = A.cols();
    /**
     * Set Q back to the identity, set R equal to A
     */
    Q.resize(m, m);
    Q.setIdentity();
    R = A;
    /**
     * Assembling the matrix Q by applying the Givens rotation at each
     * iteration and applying each component of Q to R in order to make it triangular
     */
    double c, s, a, b, tmp;
#pragma omp parallel num_threads(4)
    {
        for (int j = 0; j < n; j++)
        {
            for (int i = m - 1; i > j; i--)
            {
                if (R(i, j) != 0)
                {
#pragma omp single
                    {
                        /**
                         * Givens rotation gets applied by calculating the values of:
                         * c: cosine of the angle of rotation
                         * s: sine of the angle of rotation
                         *
                         * The idea is to calculate the rotations on smaller vectors, iterating on the cells
                         * below the diagonal to pull them to zero
                         */

                        a = R(i - 1, j);
                        b = R(i, j);

                        if (std::abs(a) > std::abs(b))
                        {
                            if (a != 0.0)
                            {
                                int segno = std::signbit(a) ? -1 : 1;
                                c = segno / std::sqrt(1 + (b / a) * (b / a));
                                s = std::abs(c / a) * b;
                            }
                            else
                            {
                                c = 0.0;
                                s = (a >= 0.0 ? 1.0 : -1.0);
                            }
                        }
                        else if (b != 0.0)
                        {
                            int segno = (std::signbit(b) ? -1 : 1);
                            s = segno / std::sqrt(1 + (a / b) * (a / b));
                            c = std::abs(s / b) * a;
                        }
                        else
                        {
                            s = 0.0;
                            c = (b >= 0.0 ? 1.0 : -1.0);
                        }

                        /**
                         * Instead of creating the Givens matrix, I do directly the computation on Q and R
                         * In addition I use a temporal variable to avoid changing the matrix before the computation
                         * Since there is dependencies within the for loop, I use atomic add
                         */

                        tmp = 0.0;
                    }

#pragma omp for private(tmp)
                    for (int k = 0; k < n; k++)
                    {
                        tmp = c * R.coeffRef(i - 1, k) + s * R.coeffRef(i, k);
                        R.coeffRef(i, k) = -s * R.coeffRef(i - 1, k) + c * R.coeffRef(i, k);
                        R.coeffRef(i - 1, k) = tmp;
                    }

#pragma omp for private(tmp)
                    for (int k = 0; k < m; k++)
                    {
                        tmp = Q.coeffRef(k, i - 1) * c + Q.coeffRef(k, i) * s;
                        Q.coeffRef(k, i) = Q.coeffRef(k, i - 1) * -s + Q.coeffRef(k, i) * c;
                        Q.coeffRef(k, i - 1) = tmp;
                    }

                    R.coeffRef(i, j) = 0.;
                }
            }
        }
    }
    return std::make_tuple(Q, R);
}

void setQR_for_svd_parallel(Matrix Q, Matrix R)
{
    int m = Q.rows();
    int n = R.rows();
#pragma omp parallel for num_threads(2)
    for (int i = 0; i < n - 1; i++)
    {
        if (R(i, i) > 0)
        {
#pragma omp parallel for num_threads(2) shared(R)
            for (int j = i; j < n; j++)
            {
                R(i, j) *= -1;
            }
#pragma omp parallel for num_threads(2) shared(Q)
            for (int j = 0; j < m; j++)
            {
                Q(j, i) *= -1;
            }
        }
    }
}

std::tuple<Matrix, Matrix> QR_Decomposition::HouseHolder_solve_2_parallel(const Matrix &A)
{
    int m = A.rows();
    int n = A.cols();
    Vector u(m);

    // Initialize matrix Q (size m x m), matrix R(m x n) and rotation matrix P(m x m)
    Matrix I(m, m);
    I.setIdentity();

    Matrix Q = I;
    Matrix P = I;
    Matrix R = A;

    double mag, mag2;

    for (int j = 0; j < std::min(m, n); j++)
    {
        mag = 0.;
        u = Vector::Zero(m, 1);
#pragma omp parallel for reduction(+ : mag)
        for (int i = j; i < m; i++)
        {
            u(i, 0) = R(i, j);
            mag += u(i, 0) * u(i, 0);
        }

        mag2 = mag - u(j, 0) * u(j, 0);
        mag = sqrt(mag);

        u(j, 0) = (u(j, 0) < 0) ? u(j, 0) + mag : u(j, 0) - mag;
        mag2 += u(j, 0) * u(j, 0);
        mag2 = sqrt(mag2);

        if (mag2 >= 0.0000000001)
        {
#pragma omp parallel for simd
            for (int i = j; i < m; i++)
                u(i, 0) /= mag2;

            // Computing P at the j-th iterate and applying the rotation to R,Q          
            P = I - 2.0 * u * u.transpose();
            R = P * R;
            Q = Q * P;
        }
    }

#pragma omp parallel for collapse(2)
    for (int j = 0; j < n; j++)
    {
        for (int i = j + 1; i < m; i++)
        {
            R(i, j) = 0.;
        }
    }

    return std::make_tuple(Q, R);
}

std::tuple<Matrix, Matrix> QR_Decomposition::HouseHolder_solve_parallel(const Matrix &A)
{
    int m = A.rows();
    int n = A.cols();

    /**
     * Initialize matrix Q (size m x m), matrix R(m x n) and rotation matrix P(m x m)
     */
    Matrix Q(m, m);
    Q.setIdentity();
    Matrix R = A;

    double normx = 0.;
    double u1;
    Vector w(m);
    double tau;
    Vector tmp_R(n);
    Vector tmp_Q(m);

    for (int j = 0; j < std::min(n, m); j++)
    {
        int s = ((R(j, j) >= 0) ? -1 : 1);
#pragma omp parallel
        {
#pragma omp for reduction(+ : normx)
            for (int i = j; i < m; i++)
            {
                normx += R(i, j) * R(i, j);
            }
#pragma omp single
            {
                normx = std::sqrt(normx);
                u1 = R(j, j) - s * normx;
            }
#pragma omp for
            for (int i = j + 1; i < m; i++)
            {
                w(i, 0) = R(i, j) / u1;
            }

#pragma omp single
            {

                w(j, 0) = 1;
                tau = -s * u1 / normx;
            }
/**
 * Computation of Q and R
 */
#pragma omp for
            for (int l = j; l < n; l++)
            {
                for (int i = j; i < m; i++)
                {

                    tmp_R(l, 0) += w(i, 0) * R(i, l);
                }
            }

#pragma omp for collapse(2)
            for (int i = j; i < m; i++)
            {

                for (int l = j; l < n; ++l)
                {

                    R(i, l) -= (tau * w(i, 0)) * tmp_R(l, 0);
                }
            }

#pragma omp for
            for (int k = 0; k < m; k++)
            {
                for (int l = j; l < m; ++l)
                {

                    tmp_Q(k, 0) += Q(k, l) * w(l, 0);
                }
            }

#pragma omp for collapse(2)
            for (int k = 0; k < m; ++k)
            {
                for (int l = j; l < m; ++l)
                {

                    Q(k, l) -= tmp_Q(k, 0) * w(l, 0) * tau;
                }
            }
        }
        tmp_R = 0 * tmp_R;
        tmp_Q = 0 * tmp_Q;

        normx = 0.0;
    }

    return std::make_tuple(Q, R);
}
