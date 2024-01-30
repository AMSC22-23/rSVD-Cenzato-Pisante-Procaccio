#include "QR_Decomposition.hpp"
void exportmatrix(Matrix A, std::string outputFileName)
{
    // Write the matrix to the file
    std::ofstream outputFile(outputFileName);
    if (outputFile.is_open())
    {
        int rows = A.rows(), cols = A.cols();
        // Write dimensions to the first row
        outputFile << rows << " " << cols << std::endl;

        // Write matrix data
        for (int i = 0; i < rows; ++i)
        {
            for (int j = 0; j < cols; ++j)
            {
                outputFile << std::setw(8) << std::fixed << std::setprecision(4) << A(i, j) << " ";
            }
            outputFile << std::endl;
        }
        std::cout << "Computed matrix has been written to " << outputFileName << std::endl;

        // Close the file
        outputFile.close();
    }
    else
    {
        std::cerr << "Error opening file for writing." << std::endl;
    }
}

int main(int argc, char **argv)
{
    int m = (argc >= 2) ? std::stoul(argv[1]) : 60;
    int n = (argc >= 3) ? std::stoul(argv[2]) : 60;

    Vector sparsity(15), n_full(15), dur_g(100), dur_p(100), dur_hh(100);
    /*
        sparsity(0, 0) = 5;
        sparsity(1, 0) = 10;
        sparsity(2, 0) = 20;
        sparsity(3, 0) = 50;
        sparsity(4, 0) = 75;
        sparsity(5, 0) = 100;
        sparsity(6, 0) = 125;
        sparsity(7, 0) = 150;
        sparsity(8, 0) = 200;
        sparsity(9, 0) = 275;
        sparsity(10, 0) = 350;
        sparsity(11, 0) = 500;
        sparsity(12, 0) = 750;
        sparsity(13, 0) = 1000;
        sparsity(14, 0) = 1500;
        */

    Matrix A(1000, 1000);

    for (int k = 10; k < 500; k += 10)
    {
        A = Eigen::MatrixXd::Zero(1000, 1000);
        for (int i = 0; i < 1000; ++i)
        {
            A(i, i) = 2.0; // Elementi diagonali
            for (int j = 0; j < k; j++)
            {
                if (i + j < 1000)
                {
                    A(i, i + j) = -1.0; // Elementi sopra la diagonale
                    A(i + j, i) = -1.0; // Elementi sotto la diagonale
                }
            }
        }

        QR_Decomposition QR_A;

        /**
         * Serial execution with Givens
         */

        auto start_givens = std::chrono::high_resolution_clock::now();
        auto [Qg, Rg] = QR_A.Givens_solve(A);
        auto end_givens = std::chrono::high_resolution_clock::now();

        std::chrono::duration<double> duration_g;
        duration_g = (end_givens - start_givens);
        dur_g(k/10, 0) = duration_g.count();

        auto start_hh = std::chrono::high_resolution_clock::now();
        auto [Qh, Rh] = QR_A.HouseHolder_solve(A);
        auto end_hh = std::chrono::high_resolution_clock::now();

        std::chrono::duration<double> duration_hh;
        duration_hh = (end_hh - start_hh);
        dur_hh(k/10, 0) = duration_hh.count();

        std::cout << k / 10 << " iterazione fatta" << std::endl;
    }
    exportmatrix(dur_g, "./../python/dur_givens.txt");
    exportmatrix(dur_hh, "./../python/dur_householder.txt");
    exportmatrix(sparsity, "./../python/zeros.txt");

    return 0;
}
