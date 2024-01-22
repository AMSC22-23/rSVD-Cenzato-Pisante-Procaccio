#include "QR_Decomposition.hpp"

std::tuple<Matrix, Matrix> QR_Decomposition::Givens_solve_parallel(Matrix A)
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
#pragma omp parallel num_threads(8)
{
    for (int j = 0; j < n; j++)
    {
        for (int i = m - 1; i > j; i--)
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

std::tuple<Matrix, Matrix>  QR_Decomposition::Givens_solve_mpi(const Matrix A, const unsigned int rank, const unsigned int size)
{

    int m = A.rows();
    int n = A.cols();
    /**
     * Set Q back to the identity, set R equal to A
     */

    Matrix Q(m, m);
    Q.setIdentity();
    Matrix R = A;
    /**
     * Assembling the matrix Q by applying the Givens rotation at each
     * iteration and applying each component of Q to R in order to make it triangular
     */
    double c, s, a, b, tmp;

    int syncflag = 0;
    MPI_Status status;

    for (int tmp_rank = rank; tmp_rank < m; tmp_rank += size)
    {
        MPI_Barrier(MPI_COMM_WORLD);
        if (tmp_rank != 0 && (tmp_rank % 4 != 0))
        {
            if (rank == 0)
                std::cout.flush() << "FATTO" << std::endl;
            MPI_Recv(&syncflag, 1, MPI_INT, rank - 1, 0, MPI_COMM_WORLD, &status);
        }
        else if (tmp_rank != 0 && tmp_rank % 4 == 0)
        {
            MPI_Recv(&syncflag, 1, MPI_INT, rank - 3, 0, MPI_COMM_WORLD, &status);
        }

        if (tmp_rank <= (2 * syncflag))
        {
            std::cout.flush() << "Processore " << tmp_rank << " inizia." << std::endl;
            for (int i = m - 1; i > tmp_rank; i--)
            {
                std::cout.flush() << "Processore " << tmp_rank << "--> cella " << i << " " << tmp_rank << std::endl;

                /**
                 * Givens rotation gets applied by calculating the values of:
                 * c: cosine of the angle of rotation
                 * s: sine of the angle of rotation
                 *
                 * The idea is to calculate the rotations on smaller vectors, iterating on the cells
                 * below the diagonal to pull them to zero
                 */

                a = R(i - 1, tmp_rank);
                b = R(i, tmp_rank);

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
                if (syncflag - 2 * tmp_rank == 0 || syncflag - 2 * tmp_rank == 1)
                {
                    syncflag++;
                    std::cout.flush() << "Processore " << tmp_rank << " syncflag aumentata a " << syncflag << std::endl;
                    MPI_Bcast(&syncflag, 1, MPI_INT, rank, MPI_COMM_WORLD);
                }

                if (syncflag - 2 * tmp_rank == 2)
                {
                    if (tmp_rank == 0 || (tmp_rank + 1) % 4 != 0)
                    {
                        MPI_Send(&syncflag, 1, MPI_INT, rank + 1, 0, MPI_COMM_WORLD);
                    }
                    else if (tmp_rank != 0 && (tmp_rank + 1) % 4 == 0)
                    {
                        MPI_Send(&syncflag, 1, MPI_INT, rank - 3, 0, MPI_COMM_WORLD);
                    }

                    std::cout.flush() << "Processore " << tmp_rank + 1 << " avvisato " << std::endl;
                }

                for (int k = 0; k < n; k++)
                {
                    tmp = c * R.coeffRef(i - 1, k) + s * R.coeffRef(i, k);
                    R.coeffRef(i, k) = -s * R.coeffRef(i - 1, k) + c * R.coeffRef(i, k);
                    R.coeffRef(i - 1, k) = tmp;
                    MPI_Bcast(&R.coeffRef(i - 1, k), 1, MPI_DOUBLE, rank, MPI_COMM_WORLD);
                    MPI_Bcast(&R.coeffRef(i, k), 1, MPI_DOUBLE, rank, MPI_COMM_WORLD);
                }

                for (int k = 0; k < m; k++)
                {
                    tmp = Q.coeffRef(k, i - 1) * c + Q.coeffRef(k, i) * s;
                    Q.coeffRef(k, i) = Q.coeffRef(k, i - 1) * -s + Q.coeffRef(k, i) * c;
                    Q.coeffRef(k, i - 1) = tmp;
                    MPI_Bcast(&Q.coeffRef(k, i - 1), 1, MPI_DOUBLE, rank, MPI_COMM_WORLD);
                    MPI_Bcast(&Q.coeffRef(k, i), 1, MPI_DOUBLE, rank, MPI_COMM_WORLD);
                }

                MPI_Comm local_comm;
                MPI_Comm_split(MPI_COMM_WORLD, 2 * tmp_rank <= syncflag, rank, &local_comm);

                int local_rank;
                MPI_Comm_rank(local_comm, &local_rank);

                if (syncflag >= 2)
                    MPI_Barrier(local_comm);
                std::cout.flush() << "Superato la barriera" << rank << std::endl;
            }
        }
    }

    return std::make_tuple(Q, R);
}

std::tuple<Matrix, Matrix> QR_Decomposition::HouseHolder_solve_parallel(Matrix A)
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
    double mag, alpha;
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
#pragma omp for
        for (int i = 0; i < m; i++)
        {
            u(i, 0) = 0.;
            v(i, 0) = 0.;
        }
#endif

        /**
         * evaluating each component of the matrix R
         */
        mag = 0.0;

        for (int i = j; i < m; i++)
        {
#ifdef EIGEN
            u(i) = R(i, j);
            mag += u(i) * u(i);
#else
            u(i, 0) = R(i, j);
            mag += u(i, 0) * u(i, 0);
#endif
        }
        mag = std::sqrt(mag);
#ifdef EIGEN
        alpha = (u(j) < 0) ? mag : -mag;
#else
        alpha = (u(j, 0) < 0) ? mag : -mag;
#endif

        mag = 0.0;

        for (int i = j; i < m; i++)
        {
#ifdef EIGEN
            v(i) = (j == i) ? (u(i) + alpha) : u(i);
            mag += v(i) * v(i);
#else
            v(i, 0) = (j == i) ? (u(i, 0) + alpha) : u(i, 0);
            mag += v(i, 0) * v(i, 0);
#endif
        }
        mag = std::sqrt(mag);
        if (mag < 0.0000000001)
            continue;

        for (int i = j; i < m; i++)
            v(i, 0) /= mag;

        /**
         * Computing P at the j-th iterate and applying the rotation to R,Q
         */
        P = I - 2.0 * v * v.transpose();

#pragma omp sections
        {
#pragma omp section
            {
#pragma omp critical
                R = P * R;
            }
#pragma omp section
            {
#pragma omp critical
                Q = Q * P;
            }
        }
    }
    
#pragma omp parallel for collapse(2) num_threads(4)
    for (int j = 0; j < n; j++)
    {
        for (int i = j + 1; i < m; i++)
        {
            R(i, j) = 0.;
        }
    }

    return std::make_tuple(Q, R);
}
