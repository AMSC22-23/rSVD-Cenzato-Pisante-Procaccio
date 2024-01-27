#include "svd.hpp"

// g++ -I${mkEigenInc} svd_test.cpp svd.cpp QR_Decomposition_parallel.cpp -o svd -Wall -DEIGEN -fopenmp -DPARALLEL
// or to test with our matrix class
// g++ -std=c++20 svd_test.cpp svd.cpp QR_Decomposition_parallel.cpp -o svd -Wall -fopenmp -DPARALLEL

void exportmatrix(Matrix A, std::string outputFileName);

int main(int argc, char **argv)
{
    int flag = (argc > 1) ? std::stoi(argv[1]) : 0;
    // flag = 0: computes everything, 1: pm1 & pm2, 2: pm1 & rsvd, 3: pseudoinverse.

    /*const char* filePath;
    if (argc > 2) {
        filePath = argv[2];
        std::cout << "Attempting to open file: " << filePath << std::endl;
    } else {
        filePath = "test_matrices/matrix2.txt";  // Default file path
        std::cout << "No file path provided. Using default file: " << filePath << std::endl;
    }

    std::ifstream file(filePath);
    if (!file.is_open())
    {
        std::cerr << "Error opening file." << std::endl;
        return 1;
    }

    // Take the dimensions of the matrix
    int m, n;
    file >> m >> n;
    Matrix A(m, n);
    // Read the matrix data
    for (int i = 0; i < m; ++i)
    {
        for (int j = 0; j < n; ++j)
        {
            file >> A(i, j);
        }
    }
    // Close the file
    file.close();*/

    int m = 10, n = 10;
    // Check if command-line arguments are provided
    if (argc > 2)
    {
        m = std::stoi(argv[2]);
        n = std::stoi(argv[3]);
    }
    std::cout << "\nGaussian matrix " << m << " x " << n << std::endl;

    SVD obj;
    Matrix A = obj.genmat(m, n);

    std::cout << "\nSVD with Power Method:\n";

    // 1-st algorithm using B = At * A
    auto start = std::chrono::high_resolution_clock::now();
    auto [U_pm, s_pm, V_pm] = obj.svd_with_PM(A);
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration_pm = end - start;
    std::cout << "Time of execution power method 1-st algorithm: " << duration_pm.count() << " s" << std::endl;
    std::cout << "pm1 : || A - U * S * Vt || = " << (A - obj.mult_SVD(U_pm, s_pm, V_pm)).norm() << std::endl;
    // exportmatrix(U_pm, "U_pm.txt");
    exportmatrix(s_pm.transpose(), "s_pm.txt");
    // exportmatrix(V_pm.transpose(), "Vt_pm.txt");

    if (!flag || flag == 1)
    {
        // 2-nd algorithm using A
        start = std::chrono::high_resolution_clock::now();
        auto [U_pm2, s_pm2, V_pm2] = obj.svd_with_PM2(A);
        end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> duration_pm2 = end - start;
        std::cout << "Time of execution power method 2-nd algorithm: " << duration_pm2.count() << " s" << std::endl;
        std::cout << "pm2 : || A - U * S * Vt || = " << (A - obj.mult_SVD(U_pm2, s_pm2, V_pm2)).norm() << std::endl;
        std::cout << "Difference eigenvalues = " << (s_pm2 - s_pm).norm() << std::endl;
    }

    if (!flag || flag == 2)
    {
        int r = (argc > 3) ? std::stoi(argv[4]) : 5;
        start = std::chrono::high_resolution_clock::now();
        auto [U_rsvd, s_rsvd, V_rsvd] = obj.rsvd(A, r, 0, 1);
        end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> duration_rsvd = end - start;
        std::cout << "\nrSVD ( r = " << r << " ):\n";
        std::cout << "Time of execution rSVD algorithm: " << duration_rsvd.count() << " s" << std::endl;
        std::cout << "Norm of A - U * S * Vt = " << (A - obj.mult_SVD(U_rsvd, s_rsvd, V_rsvd)).norm() << std::endl;

        double SpeedUp = duration_pm.count() / duration_rsvd.count();
        std::cout << "Speed Up randomized: " << SpeedUp << std::endl;
#ifdef EIGEN
        std::cout << "Norm of difference first r eigenvalues = " << (s_rsvd - s_pm.head(s_rsvd.rows())).norm() / s_rsvd.rows() << std::endl;
#endif

        // exportmatrix(U_rsvd, "U_rsvd.txt");
        exportmatrix(s_rsvd.transpose(), "s_rsvd.txt");
        // exportmatrix(V_rsvd.transpose(), "Vt_rsvd.txt");
    }

    if (flag == 3)
    {
        std::cout << "\nPseudo-inverse :\n";
        start = std::chrono::high_resolution_clock::now();
        auto A_inv = obj.pseudoinverse(A);
        end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> duration = end - start;
        std::cout << "Time of execution to compute the pseudo-inverse: " << duration.count() << " s" << std::endl;
        // exportmatrix(A_inv,"inv.txt");
    }

    return 0;
}

void exportmatrix(Matrix A, std::string outputFileName)
{
    // Write the matrix to the file
    std::ofstream outputFile(outputFileName);
    if (outputFile.is_open())
    {
        outputFile << A.rows() << " " << A.cols() << std::endl;
        outputFile << A;
        std::cout << "Computed matrix has been written to " << outputFileName << std::endl;
        // Close the file
        outputFile.close();
    }
    else
    {
        std::cerr << "Error opening file for writing." << std::endl;
    }
}