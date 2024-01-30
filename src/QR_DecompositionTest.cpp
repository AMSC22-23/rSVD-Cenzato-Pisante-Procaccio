#include "QR_Decomposition.hpp"
#include <cstdlib>

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

    Vector m_full(60), n_full(60), dur_g(60), dur_p(60), dur_hh1(60), dur_hh1p(60), dur_hh2(60), dur_hh2p(60), speedupg(60), speeduphh1(60), speeduphh2(60);

    for (int i = 0; i < 60; i++)
    {
        m_full(i, 0) = (i + 1) * 50;
    }

    std::srand(std::time(0));
    for (int i = 0; i < 60; i++)
    {
        int m = m_full(i, 0);
        int n = m_full(i, 0);
        Matrix A(m, n);

        for (int j = 0; j < n; j++)
        {
            for (int i = 0; i < m; i++)
            {
                A(i, j) = rand() % 10;
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
        dur_g(i, 0) = duration_g.count();

        /**
         * Serial execution with HouseHolder
         */

        auto start_serial = std::chrono::high_resolution_clock::now();
        auto [Q, R] = QR_A.HouseHolder_solve(A);
        auto end_serial = std::chrono::high_resolution_clock::now();

        std::chrono::duration<double> duration_hh1 = end_serial - start_serial;
        dur_hh1(i, 0) = duration_hh1.count();

        /**
         * Serial execution with Givens

    */
        auto start_givensp = std::chrono::high_resolution_clock::now();
        auto [Qgp, Rgp] = QR_A.Givens_solve_parallel(A);
        auto end_givensp = std::chrono::high_resolution_clock::now();

        std::chrono::duration<double> duration_gp;
        duration_gp = (end_givensp - start_givensp);
        dur_p(i, 0) = duration_gp.count();

        /**
         * Parallel execution of HouseHolder on OpenMP
         */

        std::cout << std::endl;

        auto start_hh1p = std::chrono::high_resolution_clock::now();
        auto [Qhh1p, Rhh1p] = QR_A.HouseHolder_solve_parallel(A);
        auto end_hh1p = std::chrono::high_resolution_clock::now();

        std::chrono::duration<double> duration_hh1p = end_hh1p - start_hh1p;
        dur_hh1p(i, 0) = duration_hh1p.count();


        std::cout << "Time of execution serial givens: " << duration_g.count() << " s" << std::endl;
        std::cout << "Time of execution parallel givens: " << duration_gp.count() << " s" << std::endl;
        std::cout << "Time of execution serial HouseHolder 1: " << duration_hh1.count() << " s" << std::endl;
        std::cout << "Time of execution parallel HouseHolder 1: " << duration_hh1p.count() << " s" << std::endl;

        double SpeedUpG = duration_g.count() / duration_gp.count();
        double SpeedUpHH1 = duration_hh1.count() / duration_hh1p.count();

        speedupg(i, 0) = SpeedUpG;
        speeduphh1(i, 0) = SpeedUpHH1;

        std::cout << "Speed Up Givens: " << SpeedUpG << std::endl;
        std::cout << "Speed Up HouseHolder 1: " << SpeedUpHH1 << std::endl;
    }

    exportmatrix(dur_g, "./../python/dur_g.txt");
    exportmatrix(dur_p, "./../python/dur_p.txt");
    exportmatrix(dur_hh1, "./../python/dur_hh1.txt");
    exportmatrix(dur_hh1p, "./../python/dur_hh1p.txt");
    exportmatrix(m_full, "./../python/sizem.txt");
    exportmatrix(speedupg, "./../python/SPg.txt");
    exportmatrix(speeduphh1, "./../python/SPhh1.txt");


    return 0;
}
