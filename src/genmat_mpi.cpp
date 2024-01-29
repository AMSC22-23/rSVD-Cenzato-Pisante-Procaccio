#include "utils.hpp"
#include <mpi.h>

// mpic++ -I${mkEigenInc} -std=c++20 -Wall test_mpi.cpp -o test -DEIGEN
// mpirun -n 8 ./test 1000 1000

void genmat(const int m, const int n, Matrix &global_matrix)
{
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // Calculate the number of rows to generate per process
    int rows_per_process = m / size;
    int extra_rows = m % size;

    // Distribute the remaining rows among the first processors
    int rows_local = (rank < extra_rows) ? rows_per_process + 1 : rows_per_process;

    // Determine the starting row for this process
    int offset = rank * rows_per_process + std::min(rank, extra_rows);

    // Generate local matrix
    Matrix local_matrix = global_matrix.block(offset, 0, rows_local, n);

    std::random_device rd;
    std::mt19937 reng(rd() + rank);
    std::normal_distribution<double> dice(0.0, 1.0);

    for (int i = 0; i < rows_local; i++)
    {
        for (int j = 0; j < n; j++)
        {
            local_matrix(i, j) = dice(reng);
        }
    }

    // Gather the local matrices into a single matrix on the root process
    std::vector<int> recv_counts(size);     // local rows for each process
    std::vector<int> displacements(size);
    for (int i = 0; i < size; ++i)
    {
        recv_counts[i] = (i < extra_rows) ? (rows_per_process + 1) * n : rows_per_process * n; 
        displacements[i] = (i < extra_rows) ? (i * n * (rows_per_process + 1)) : (extra_rows * n *
                            (rows_per_process + 1) + (i - extra_rows) * n * rows_per_process);
    }

    // Gather local results
    MPI_Gatherv(local_matrix.data(), rows_local * n, MPI_DOUBLE, global_matrix.data(), recv_counts.data(),
                displacements.data(), MPI_DOUBLE, 0, MPI_COMM_WORLD);
}

int main(int argc, char **argv)
{
    MPI_Init(NULL, NULL);

    int m = 10, n = 10;
    // Check if command-line arguments are provided
    if (argc >= 3)
    {
        m = std::stoi(argv[1]);
        n = std::stoi(argv[2]);
    }
    Matrix A = Matrix::Zero(m,n);

    auto start = std::chrono::high_resolution_clock::now();
    genmat(m, n, A);
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration_pm = end - start;

    std::cout << "Time of execution: " << duration_pm.count() << " s" << std::endl;

    MPI_Finalize();

    return 0;
}