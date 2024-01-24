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

    Vector m_full(15), n_full(15), dur_g(15), dur_p(15),dur_hh(15);
/*
    m_full(0, 0) = 5;
    m_full(1, 0) = 10;
    m_full(2, 0) = 20;
    m_full(3, 0) = 50;
    m_full(4, 0) = 75;
    m_full(5, 0) = 100;
    m_full(6, 0) = 125;
    m_full(7, 0) = 150;
    m_full(8, 0) = 200;
    m_full(9, 0) = 275;
    m_full(10, 0) = 350;
    m_full(11, 0) = 500;
    m_full(12, 0) = 750;
    m_full(13, 0) = 1000;
    m_full(14, 0) = 1500;

    n_full(0, 0) = 5;
    n_full(1, 0) = 10;
    n_full(2, 0) = 20;
    n_full(3, 0) = 50;
    n_full(4, 0) = 75;
    n_full(5, 0) = 100;
    n_full(6, 0) = 125;
    n_full(7, 0) = 150;
    n_full(8, 0) = 200;
    n_full(9, 0) = 275;
    n_full(10, 0) = 350;
    n_full(11, 0) = 500;
    n_full(12, 0) = 750;
    n_full(13, 0) = 1000;
    n_full(14, 0) = 1500;

    for (int i = 0; i < 15; i++)
    {
        int m = m_full(i, 0);
        int n = n_full(i, 0);
        Matrix A(m, n);*/
        /*
          for (int j = 0; j < n; j++)
          {
              for (int i = 0; i < m; i++)
              {
                  A(i, j) = 0.5*i+std::exp(j)+j*j-i*i*i;
              }
          }
      */
        // Inizializzazione degli elementi della matrice tridiagonale
        Matrix A(m,n);
        for (int i = 0; i < m; ++i)
        {
            A(i, i) = 2.0; // Elementi diagonali
            if (i < m - 1)
            {
                A(i, i + 1) = -1.0; // Elementi sopra la diagonale
                A(i + 1, i) = -1.0; // Elementi sotto la diagonale
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
        //dur_g(i, 0) = duration_g.count();

        std::cout << "Serial Givens:" << std::endl;
        std::cout << "R=" << std::endl;
#ifdef EIGEN
        std::cout << Rg << std::endl;
#else
        Rg.print(std::cout);
#endif
        std::cout << "Q=" << std::endl;
#ifdef EIGEN
        std::cout << Qg << std::endl;
#else
        Qg.print(std::cout);
#endif

        auto start_hh = std::chrono::high_resolution_clock::now();
        auto [Qh, Rh] = QR_A.HouseHolder_solve_2(A);
        auto end_hh = std::chrono::high_resolution_clock::now();

        std::chrono::duration<double> duration_hh;
        duration_hh = (end_hh - start_hh);
        //dur_hh(i, 0) = duration_hh.count();

        std::cout << "HouseHolder:" << std::endl;
        std::cout << "R=" << std::endl;
#ifdef EIGEN
        std::cout << Rh << std::endl;
#else
        Rh.print(std::cout);
#endif

        std::cout << "Q=" << std::endl;

#ifdef EIGEN
        std::cout << Qh << std::endl;
#else
        Qh.print(std::cout);
#endif

        /**
         * Serial execution with HouseHolder


        auto start_serial = std::chrono::high_resolution_clock::now();
        auto [Q, R] = QR_A.HouseHolder_solve(A);
        auto end_serial = std::chrono::high_resolution_clock::now();

        std::chrono::duration<double> duration_s = end_serial - start_serial;
     */
        /*
        std::cout<<"Serial HouseHolder:"<<std::endl;
        std::cout<<"R="<<std::endl;
        #ifdef EIGEN
            std::cout<<R<<std::endl;
        #else
            R.print(std::cout);
        #endif
        std::cout<<"Q="<<std::endl;
        #ifdef EIGEN
            std::cout<<Q<<std::endl;
        #else
            Q.print(std::cout);
        #endif
        */

        /**
         * Serial execution with Givens

    */
        auto start_givensp = std::chrono::high_resolution_clock::now();
        auto [Qgp, Rgp] = QR_A.Givens_solve_parallel(A);
        auto end_givensp = std::chrono::high_resolution_clock::now();

        std::chrono::duration<double> duration_gp;
        duration_gp = (end_givensp - start_givensp);
        //dur_p(i, 0) = duration_gp.count();
        /*
            std::cout << std::endl;
            std::cout << "Parallel Givens:" << std::endl;
            std::cout << "R=" << std::endl;
        #ifdef EIGEN
            std::cout << Rgp << std::endl;
        #else
            Rgp.print(std::cout);
        #endif
            std::cout << "Q=" << std::endl;
        #ifdef EIGEN
            std::cout << Qgp << std::endl;
        #else
            Qgp.print(std::cout);
        #endif
        */

        /**
         * Parallel execution of HouseHolder on OpenMP


        std::cout << std::endl;

        auto start_parallel = std::chrono::high_resolution_clock::now();
        auto [Qp, Rp] = QR_A.HouseHolder_solve_parallel(A);
        auto end_parallel = std::chrono::high_resolution_clock::now();

        std::chrono::duration<double> duration_p = end_parallel - start_parallel;
    */
        /*
        std::cout<<"Parallel HouseHolder:"<<std::endl;
        std::cout<<"R="<<std::endl;
        #ifdef EIGEN
            std::cout<<Rp<<std::endl;
        #else
            Rp.print(std::cout);
        #endif
        std::cout<<"Q="<<std::endl;
        #ifdef EIGEN
            std::cout<<Qp<<std::endl;
        #else
            Qp.print(std::cout);
        #endif
        */
    
    /*
    std::cout << "Time of execution serial givens: " << duration_g.count() << " s" << std::endl;
    // std::cout << "Time of execution serial HouseHolder: " << duration_s.count() << " s" << std::endl;
    std::cout << "Time of execution parallel givens: " << duration_gp.count() << " s" << std::endl;
    // std::cout << "Time of execution parallel: " << duration_p.count() << " s" << std::endl;

    double SpeedUpG = duration_g.count() / duration_gp.count();
    // double SpeedUpHH = duration_s.count() / duration_p.count();
    std::cout << "Speed Up Givens: " << SpeedUpG << std::endl;
    // std::cout << "Speed Up HouseHolder: " << SpeedUpHH << std::endl;
*/

    //exportmatrix(dur_g, "./../python/dur_g.txt");
    //exportmatrix(dur_p, "./../python/dur_p.txt");
    //exportmatrix(m_full, "./../python/sizem.txt");
    //exportmatrix(n_full, "./../python/sizen.txt");

    return 0;
}
