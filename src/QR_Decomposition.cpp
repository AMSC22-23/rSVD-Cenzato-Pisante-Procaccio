#include "QR_Decomposition.hpp"

std::tuple<Matrix, Matrix> QR_Decomposition::Givens_solve(Matrix A)
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
    double c, s,a,b,tmp;
    for (int j = 0; j < n; j++)
    {
        for (int i = m - 1; i > j; i--)
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
                    c = segno / sqrt(1 + (b / a) * (b / a));
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
                s = segno / sqrt(1 + (a / b) * (a / b));
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
             */
            tmp = 0.0;
            for (int k = 0; k < n; k++)
            {
                tmp = c * R(i - 1, k) + s * R(i, k);
                R(i, k) = -s * R(i - 1, k) + c * R(i, k);
                R(i - 1, k) = tmp;
            }

            /**
             * Forcing rotated component to zero to avoid floating point approximations
             */

            R(i, j) = 0;

            /**
             * Computation of Q, at each iterate
             */
            for (int k = 0; k < m; k++)
            {
                tmp = Q(k, i - 1) * c + Q(k, i) * s;
                Q(k, i) = Q(k, i - 1) * -s + Q(k, i) * c;
                Q(k, i - 1) = tmp;
            }
        }
    }

    return std::make_tuple(Q, R);
}

std::tuple<Matrix, Matrix> QR_Decomposition::HouseHolder_solve(Matrix A)
{

    int m = A.rows();
    int n = A.cols();

    /**
     * Initialize v,u
     */
    Vector u(m), v(m);

    /**
     * Initialize matrix Q (size m x m), matrix R(m x n) and rotation matrix P(m x m)
     */
    Matrix I(m, m);
    I.setIdentity();

    Matrix Q = I;
    Matrix P = I;
    Matrix R = A;

    /**
     * Starting the computation of Q,R
     */
    double mag = 0.0;
    double alpha = 0.0;
    /**
     * Initialize u,v to zero at each iteration-i
     */
    u = Vector(m);
    v = Vector(m);

    for (int j = 0; j < std::min(m,n); j++)
    {

#ifdef EIGEN
        u.setZero();
        v.setZero();
#else
        u = Vector(m);
        v = Vector(m);
#endif

        /**
         * evaluating each component of the matrix R
         */
        mag = 0.0;
        for (int i = j; i < m; i++)
        {
/**
 * Computation of u
 */
#ifdef EIGEN
            u(i) = R(i, j);
            mag += u(i) * u(i);
#else
            u(i, 0) = R(i, j);
            mag += u(i, 0) * u(i, 0);
#endif
        }
        mag = sqrt(mag);
#ifdef EIGEN
        alpha = (u(j) < 0) ? mag : -mag;
#else
        alpha = (u(j, 0) < 0) ? mag : -mag;
#endif

        mag = 0.0;

        for (int i = j; i < m; i++)
        {
/**
 * Computation ov v
 */
#ifdef EIGEN
            v(i) = (j == i) ? (u(i) + alpha) : u(i);
            mag += v(i) * v(i);
#else
            v(i, 0) = (j == i) ? (u(i, 1) + alpha) : u(i, 0);
            mag += v(i, 0) * v(i, 0);
#endif
        }

        mag = sqrt(mag);

        if (mag < 0.0000000001)
            continue;

        for (int i = j; i < m; i++)
            v(i, 0) /= mag;

        /**
         * Computing P at the j-th iterate and applying the rotation to R,Q
         */
        P = I - 2.0 * v * v.transpose();
        R = P * R;
        Q = Q * P;

        /**
         * Force j-th col of R to zero
         */
    }
    for (int j = 0; j < n; j++)
    {
        for (int i = j + 1; i < m; i++)
        {
            R(i, j) = 0.;
        }
    }
    return std::make_tuple(Q, R);
}

void setQR_for_svd(Matrix Q, Matrix R)
{
    int m = Q.rows();
    int n = R.rows();
    for (int i = 0; i < n - 1; i++)
    {
        if (R(i, i) > 0)
        {
            for (int j = i; j < n; j++)
            {
                R(i, j) *= -1;
            }
            for (int j = 0; j < m; j++)
            {
                Q(j, i) *= -1;
            }
        }
    }
}