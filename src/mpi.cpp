#include <mpi.h>
#include <Eigen/Dense>
#include <tuple>
#include <iostream>

using Matrix = Eigen::MatrixXd;
using Vector = Eigen::VectorXd;

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

void applyGivensRotation(Matrix &R, Matrix &Q, int row, double c, double s, int rank)
{
    double tmp;

    for (int k = 0; k < R.cols(); k++)
    {
        tmp = c * R(row - 1, k) + s * R(row, k);
        R(row, k) = -s * R(row - 1, k) + c * R(row, k);
        R(row - 1, k) = tmp;
        MPI_Bcast(&R.coeffRef(row - 1, k), 1, MPI_DOUBLE, rank, MPI_COMM_WORLD);
        MPI_Bcast(&R.coeffRef(row, k), 1, MPI_DOUBLE, rank, MPI_COMM_WORLD);
    }

    for (int k = 0; k < Q.rows(); k++)
    {
        tmp = Q(k, row - 1) * c + Q(k, row) * s;
        Q(k, row) = Q(k, row - 1) * -s + Q(k, row) * c;
        Q(k, row - 1) = tmp;
        MPI_Bcast(&Q.coeffRef(k, row - 1), 1, MPI_DOUBLE, rank, MPI_COMM_WORLD);
        MPI_Bcast(&Q.coeffRef(k, row), 1, MPI_DOUBLE, rank, MPI_COMM_WORLD);
    }
}

std::tuple<int *, int *> calculateGathervParams(int rows, int cols, int size)
{
    int *recvcounts = new int[size];
    int *displs = new int[size];

    int colPerProcess = cols / size;
    int remainingCols = cols % size;
    int currentDispl = 0;

    for (int i = 0; i < size; ++i)
    {
        recvcounts[i] = colPerProcess;
        if (i < remainingCols)
        {
            recvcounts[i]++;
        }

        displs[i] = currentDispl;
        currentDispl += recvcounts[i];
    }

    return std::make_tuple(recvcounts, displs);
}

std::tuple<Matrix, Matrix> Givens_solve_mpi(Matrix R, int rank, int size)
{

    int m = R.rows();
    int n = R.cols();
    /**
     * Set Q back to the identity, set R equal to A
     */

    Matrix Q(m, m);
    Q.setIdentity();
    /**
     * Assembling the matrix Q by applying the Givens rotation at each
     * iteration and applying each component of Q to R in order to make it triangular
     */
    double tmp;
    int syncflag = 0;
    MPI_Status status;

    int *sendcounts;
    int *displs;
    std::tie(sendcounts, displs) = calculateScattervParams(m, n, size);

    Matrix localR(m, sendcounts[rank]);

    MPI_Scatterv(R.data(), sendcounts, displs, MPI_DOUBLE, localR.data(), localR.size(), MPI_DOUBLE, 0, MPI_COMM_WORLD);

    for (int tmp_rank = rank; tmp_rank < n; tmp_rank += size)
    {

        if (tmp_rank != 0 && (tmp_rank % 4 != 0))
        {
            MPI_Recv(&syncflag, 1, MPI_INT, rank - 1, 0, MPI_COMM_WORLD, &status);
        }
        else if (tmp_rank != 0 && tmp_rank % 4 == 0)
        {
            std::cout.flush() << "FATTO" << std::endl;
            MPI_Recv(&syncflag, 1, MPI_INT, 3, 0, MPI_COMM_WORLD, &status);
        }

        if (tmp_rank <= (2 * syncflag))
        {
            std::cout.flush() << "Processore " << rank << " inizia." << std::endl;
            for (int i = m - 1; i > tmp_rank; i--)
            {
                std::cout.flush() << "Processore " << rank << "--> cella " << i << " " << tmp_rank << std::endl;

                auto [s, c] = findSine_Cosine(R, i, tmp_rank);

                /**
                 * Increase syncflag to set the starting of the next tmp_rank process
                 */
                if (syncflag - 2 * tmp_rank == 0 || syncflag - 2 * tmp_rank == 1)
                {

                    syncflag++;
                    std::cout.flush() << "Processore " << rank << " syncflag aumentata a " << syncflag << std::endl;
                    MPI_Bcast(&syncflag, 1, MPI_INT, rank, MPI_COMM_WORLD);

                    if (syncflag - 2 * tmp_rank == 2)
                    {
                        /**
                         * tmp_rank processor tells tmp_rank + 1 processor to start the computation
                         */
                        if (tmp_rank == 0 || (tmp_rank + 1) % 4 != 0)
                        {
                            MPI_Send(&syncflag, 1, MPI_INT, rank + 1, 0, MPI_COMM_WORLD);
                        }
                        else if (tmp_rank != 0 && (tmp_rank + 1) % 4 == 0)
                        {
                            std::cout.flush() << "SONO " << rank << std::endl;
                            MPI_Send(&syncflag, 1, MPI_INT, rank - 3, 0, MPI_COMM_WORLD);
                        }

                        std::cout.flush() << "Processore " << rank + 1 << " avvisato " << std::endl;
                    }
                }
                if (tmp_rank == n - 2)
                {
                    MPI_Send(&syncflag, 1, MPI_INT, rank + 1, 0, MPI_COMM_WORLD);
                }

                applyGivensRotation(R, Q, i, c, s, rank);
            }
        }
    }

    std::cout.flush() << rank << " ha finito l'esecuzione" << std::endl;
    MPI_Barrier(MPI_COMM_WORLD);
    if (rank == 0)
    {
        std::cout << "Siamo pronti a ricostituire R" << std::endl;
    }
    // Calcolare parametri di MPI_Gatherv
    int *recvcounts;
    int *displs_gather;
    std::tie(recvcounts, displs_gather) = calculateGathervParams(m, n, size);

    /**
     * Memory allocation of R
    */
    Matrix reconstructedR(m, n);


    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Gatherv(localR.data(), localR.size(), MPI_DOUBLE, reconstructedR.data(), recvcounts, displs_gather, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    return std::make_tuple(Q, reconstructedR);
}

int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int m = 12;
    int n = 12;

    Matrix A(m, n);
    if (rank == 0)
    {
        for (int j = 0; j < n; j++)
        {
            for (int i = 0; i < m; i++)
            {
                A(i, j) = 0.5 * i + std::exp(j) + j * j - i * i * i;
            }
        };

        MPI_Bcast(&A, 1, MPI_DOUBLE, rank, MPI_COMM_WORLD);
    }
    Matrix Q;
    Matrix R;
    if (rank == 0)
    {
        std::cout << "Inizio esecuzione" << std::endl;
    }
    MPI_Barrier(MPI_COMM_WORLD);

    auto [Qgp, Rgp] = Givens_solve_mpi(A, rank, size);

    MPI_Barrier(MPI_COMM_WORLD);
    std::cout.flush() << "Abbiamo superato la barriera finale" << std::endl;
    MPI_Barrier(MPI_COMM_WORLD);

    if (rank == 0)
    {
        std::cout << "R=" << std::endl;
        std::cout << R << std::endl;

        std::cout << "Q=" << std::endl;
        std::cout << Q << std::endl;
    }

    MPI_Finalize();
    return 0;
}