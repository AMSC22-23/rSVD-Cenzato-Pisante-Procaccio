#include <mpi.h>
#include <Eigen/Dense>
#include <tuple>
#include <iostream>
#include <cstdlib>
#include <ctime>
#include <chrono>

using Matrix = Eigen::MatrixXd;
using Vector = Eigen::VectorXd;
/*
std::tuple<int *, int *> calculateScattervParams(int rows, int cols, int size)
{
    int *sendcounts = new int[size];
    int *displs = new int[size];

    int colPerProcess = cols / size;
    int remainingCols = cols % size;
    int currentDispl = 0;

    for (int i = 0; i < size; ++i)
    {
        sendcounts[i] = colPerProcess;
        if (i < remainingCols)
        {
            sendcounts[i]++;
        }

        displs[i] = currentDispl;
        currentDispl += sendcounts[i];
    }

    return std::make_tuple(sendcounts, displs);
}
*/
std::tuple<double, double> findSine_Cosine(Matrix R, int i, int tmp_rank)
{
    /**
     * Givens rotation gets applied by calculating the values of:
     * c: cosine of the angle of rotation
     * s: sine of the angle of rotation
     *
     * The idea is to calculate the rotations on smaller vectors, iterating on the cells
     * below the diagonal to pull them to zero
     */
    double c, s, a, b;
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
    return std::make_tuple(s, c);
}

void applyGivensRotation(Matrix &R, Matrix &Q, int row, double c, double s)
{
    double tmp;

    for (int k = 0; k < R.cols(); k++)
    {
        tmp = c * R(row - 1, k) + s * R(row, k);
        R(row, k) = -s * R(row - 1, k) + c * R(row, k);
        R(row - 1, k) = tmp;
    }

    for (int k = 0; k < Q.rows(); k++)
    {
        tmp = Q(k, row - 1) * c + Q(k, row) * s;
        Q(k, row) = Q(k, row - 1) * -s + Q(k, row) * c;
        Q(k, row - 1) = tmp;
    }
}

std::tuple<Matrix, Matrix> Givens_solve_mpi(const Matrix &A, const int rank, const int size)
{

    int m = A.rows();
    int n = A.cols();
    /**
     * Set Q back to the identity, set R equal to A
     */
    // NOTA
    Matrix Q(m, m);
    Matrix R(m, n);
    Matrix localR(2, n);
    Matrix localQ(m, 2);

    R = A;
    Q.setIdentity();
    /**
     * Assembling the matrix Q by applying the Givens rotation at each
     * iteration and applying each component of Q to R in order to make it triangular
     */

    int syncflag = 0;
    MPI_Status status;
    int flag = 0;
    int tmp = 0;

    if (rank - 1 < (std::min(m, n) + 1) / 2)
    {
        /**
         * Start to iterate over A
         */
        while (rank > 1 && syncflag < 2 * (rank - 1))
        {
            MPI_Recv(&syncflag, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
            MPI_Recv(R.data(), m * n, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &status);
            MPI_Recv(Q.data(), m * m, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &status);
        }
        if (rank != 0)
        {
            for (int tmp_rank = rank - 1; tmp_rank < std::min(m, n) - 2; tmp_rank += size - 1)
            {
                /**
                 * Passed the MPI_Recv, the tmp_rank processor starts
                 */
                tmp = 0;
                while (tmp_rank > size - 1 && tmp != -rank + 1)
                {
                    MPI_Recv(&syncflag, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
                    tmp--;
                    MPI_Send(&tmp, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
                    MPI_Recv(R.data(), m * n, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &status);
                    MPI_Recv(Q.data(), m * m, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &status);
                }

                for (int i = m - 1; i > tmp_rank; i--)
                {
                    /**
                     * Synchronize the computation at each row
                     */
                    MPI_Recv(&syncflag, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
                    MPI_Send(&i, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
                    MPI_Send(&tmp_rank, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);

                    if (R(i, tmp_rank) != 0)
                    {
                        auto [s, c] = findSine_Cosine(R, i, tmp_rank);
                        applyGivensRotation(R, Q, i, c, s);

                        localR = R.block(i - 1, 0, 2, n);
                        localQ = Q.block(0, i - 1, m, 2);

                        MPI_Send(localR.data(), 2 * n, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
                        MPI_Send(localQ.data(), 2 * m, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
                    }

                    if (std::min(n, m) - 2 - tmp_rank < size && i == tmp_rank + 1)
                    {
                        flag++;
                    }

                    MPI_Send(&flag, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);

                    if (flag == 1)
                    {
                        MPI_Send(&rank, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
                    }
                    flag = 0;

                    MPI_Recv(R.data(), m * n, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &status);
                    MPI_Recv(Q.data(), m * m, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &status);
                }
            }
        }

        /**
         * Increase syncflag to set the starting of the next tmp_rank process
         */
        else
        {
            int i, j, tmp;
            int count = 0;
            Vector rm_rank(size);

            for (int k = 0; k < size; k++)
            {
                if (k - 1 < (std::min(m, n) + 1) / 2)
                {
                    rm_rank(k) = k;
                }
                else
                {
                    rm_rank(k) = 0;
                    count++;
                    std::cout.flush() << "Stai sprecando delle risorse importanti. Processo " << k << " fermato!" << std::endl;
                }
            }

            while (count != size - 1)
            {

                syncflag++;

                for (int l = 1; l < size; l++)
                {
                    if (rm_rank(l) != 0)
                    {
                        MPI_Send(&syncflag, 1, MPI_INT, l, 0, MPI_COMM_WORLD);
                    }
                }
                for (int l = 1; l < size && syncflag > 2 * (l - 1); l++)
                {
                    if (rm_rank(l) != 0)
                    {
                        MPI_Recv(&i, 1, MPI_INT, l, 0, MPI_COMM_WORLD, &status);
                        flag = 0;
                        if (i >= 0)
                        {
                            MPI_Recv(&j, 1, MPI_INT, l, 0, MPI_COMM_WORLD, &status);

                            if (R(i, j) != 0)
                            {
                                MPI_Recv(localR.data(), 2 * n, MPI_DOUBLE, l, 0, MPI_COMM_WORLD, &status);
                                MPI_Recv(localQ.data(), 2 * m, MPI_DOUBLE, l, 0, MPI_COMM_WORLD, &status);

                                for (int k = j; k < R.cols(); k++)
                                {
                                    R(i - 1, k) = localR(0, k);
                                    R(i, k) = localR(1, k);
                                }
                                for (int k = j; k < Q.rows(); k++)
                                {
                                    Q(k, i - 1) = localQ(k, 0);
                                    Q(k, i) = localQ(k, 1);
                                }
                            }
                            MPI_Recv(&flag, 1, MPI_INT, l, 0, MPI_COMM_WORLD, &status);
                        }

                        if (flag == 1)
                        {
                            count++;
                            MPI_Recv(&tmp, 1, MPI_INT, l, 0, MPI_COMM_WORLD, &status);
                            rm_rank(tmp) = 0;
                            std::cout.flush() << "Processo " << tmp << " ha terminato" << std::endl;

                            MPI_Send(R.data(), m * n, MPI_DOUBLE, tmp, 0, MPI_COMM_WORLD);
                            MPI_Send(Q.data(), m * m, MPI_DOUBLE, tmp, 0, MPI_COMM_WORLD);
                        }
                    }
                }

                for (int l = 1; l < size; l++)
                {
                    if (rm_rank(l) != 0)
                    {
                        MPI_Send(R.data(), m * n, MPI_DOUBLE, l, 0, MPI_COMM_WORLD);
                        MPI_Send(Q.data(), m * m, MPI_DOUBLE, l, 0, MPI_COMM_WORLD);
                    }
                }
            }
        }

        if (rank == 0)
        {
            for (int j = std::min(n, m) - 2; j < n; j++)
            {
                for (int i = m - 1; i > j; i--)
                {
                    auto [s, c] = findSine_Cosine(R, i, j);
                    applyGivensRotation(R, Q, i, c, s);
                }
            }
            std::cout.flush() << "Processo " << rank << " ha terminato" << std::endl;
        }
    }

    MPI_Barrier(MPI_COMM_WORLD);

    return std::make_tuple(Q, R);
}

int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int m = 500;
    int n = 500;

    Matrix A(m, n);
    for (int i = 0; i < m; ++i)
    {
        for (int j = 0; j < n; ++j)
        {
            // Assegna un valore casuale compreso tra 1 e 100 alla cella (i, j)
            A(i, j) = rand() % 100 + 1;
        }
    }

    /*
 if (rank == 0)
 {

     A = Eigen::MatrixXd::Zero(m, n);
     for (int i = 0; i < m; ++i)
     {
         A(i, i) = 2.0; // Elementi diagonali
         if (i < m - 1)
         {
             A(i, i + 1) = -1.0; // Elementi sopra la diagonale
             A(i + 1, i) = -1.0; // Elementi sotto la diagonale
         }
     }
 }
 */

    MPI_Bcast(A.data(), m * n, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    std::cout << "Partiamo" << std::endl;
    auto start_givens = std::chrono::high_resolution_clock::now();

    MPI_Barrier(MPI_COMM_WORLD);
    auto [Qgp, Rgp] = Givens_solve_mpi(A, rank, size);

    auto end_givens = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration_g;
    duration_g = (end_givens - start_givens);

    if (rank == 0)
    {
        std::cout << "durata:" << std::endl;
        std::cout << duration_g.count() << std::endl;
    }

    /*
        if (rank == 0)
        {
            std::cout << "R=" << std::endl;
            std::cout << Rgp << std::endl;
        }
    */

    MPI_Finalize();
    return 0;
}