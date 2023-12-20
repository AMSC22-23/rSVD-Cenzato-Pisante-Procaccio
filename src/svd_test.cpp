#include "svd.hpp"

// g++ -I${mkEigenInc} svd_test.cpp svd.cpp QR_Decomposition_parallel.cpp -o svd -Wall -DEIGEN -fopenmp
// or to test with our matrix class
// g++ svd_test.cpp svd.cpp QR_Decomposition_parallel.cpp -o prova -Wall

// In parallel:
// g++ svd_test.cpp svd.cpp QR_Decomposition_parallel.cpp -o prova -Wall -fopenmp

void exportmatrix(Matrix A, std::string outputFileName){
    // Write the matrix to the file
    std::ofstream outputFile(outputFileName);
    if (outputFile.is_open()) {
        int rows = A.rows(), cols = A.cols();
        // Write dimensions to the first row
        outputFile << rows << " " << cols << std::endl;

        // Write matrix data
        for (int i = 0; i < rows; ++i) {
            for (int j = 0; j < cols; ++j) {
                outputFile << std::setw(8) << std::fixed << std::setprecision(4) << A(i,j) << " ";
            }
            outputFile << std::endl;
        }
        std::cout << "Computed matrix has been written to "<< outputFileName << std::endl;

        // Close the file
        outputFile.close();
    } else {
        std::cerr << "Error opening file for writing." << std::endl;
    }

}

int main(int argc, char ** argv){
    std::ifstream file("test_matrices/matrix2.txt");           //file with dim and then matrix
    if (!file.is_open()) {
        std::cerr << "Error opening file." << std::endl;
        return 1;
    } 

    // Take the dimensions of the matrix
    int m, n;
    file >> m >> n;
    Matrix A(m,n);
    // Read the matrix data
    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < n; ++j) {
            file >> A(i,j);
        }
    }
    // Close the file
    file.close();
        
    SVD obj;
    
    
    std::cout<<"\nSVD with Power Method:\n";

    // 1-st algorithm using B = At * A
    auto start = std::chrono::high_resolution_clock::now();
    auto [U_pm,s_pm,V_pm] = obj.svd_with_PM(A);
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration_pm = end  - start;

    // 2-nd algorithm using A
    start = std::chrono::high_resolution_clock::now();
    auto [U_pm2,s_pm2,V_pm2] = obj.svd_with_PM2(A);
    end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration_pm2 = end  - start;

    std::cout<<"Time of execution power method 1 algorithm: "<< duration_pm.count() << " s" << std::endl;
    std::cout<<"Time of execution power method 2 algorithm: "<< duration_pm2.count() << " s" << std::endl;
    std::cout<<"pm1 : || A - U * S * Vt || = "<<(A-obj.mult_parallel(U_pm,s_pm,V_pm)).norm()<<std::endl;
    std::cout<<"pm2 : || A - U * S * Vt || = "<<(A-obj.mult_parallel(U_pm2,s_pm2,V_pm2)).norm()<<std::endl;
    std::cout<<"Difference eigenvalues = "<<(s_pm2-s_pm).norm()<<std::endl;
    exportmatrix(U_pm,"U_pm.txt");
    exportmatrix(s_pm.transpose(),"s_pm.txt");
    exportmatrix(V_pm.transpose(),"Vt_pm.txt");


    //Test for svd multiplication in parallel
    /*auto start_pm = std::chrono::high_resolution_clock::now();
    std::cout<<"|| A - U * S * Vt || = "<<(A-obj.mult(U_pm,s_pm,V_pm)).norm()<<std::endl;
    auto end_pm = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration_pm = end_pm  - start_pm;
    std::cout<<"Time of execution serial: "<< duration_pm.count() << " s" << std::endl;

    start_pm = std::chrono::high_resolution_clock::now();
    std::cout<<"|| A - U * S * Vt || = "<<(A-obj.mult_parallel(U_pm,s_pm,V_pm)).norm()<<std::endl;
    end_pm = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration_pmp = end_pm  - start_pm;
    std::cout<<"Time of execution parallel: "<< duration_pmp.count() << " s" << std::endl;
    double SpeedUp = duration_pm.count()/duration_pmp.count();
    std::cout << "Speed Up: "<< SpeedUp << std::endl;*/

/*Test with Eigen
Reduced SVD with Power Method:
|| A - U * S * Vt || = 1.60171e-12
Time of execution serial: 0.0645291 s
|| A - U * S * Vt || = 1.60171e-12
Time of execution parallel: 0.0066921 s
Speed Up: 9.64258*/


    std::cout<<"\nPseudo-inverse :\n";
    start = std::chrono::high_resolution_clock::now();
    auto A_inv = obj.pseudoinverse(A);
    end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = end  - start;
    std::cout<<"Time of execution to compute the pseudo-inverse: "<< duration.count() << " s" << std::endl;
    //exportmatrix(obj.pseudoinverse(A),"inv.txt");

     
    /*auto start_qr = std::chrono::high_resolution_clock::now();
    auto [U_qr,s_qr,V_qr] = obj.svd_with_qr(A);
    auto end_qr = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration_qr = end_qr  - start_qr;
    std::cout<<"\nSVD with QR :\n";
    std::cout<<"Time of execution QR algorithm: "<< duration_qr.count() << " s" << std::endl;
    std::cout<<"Norm of A - U * S * Vt = "<<(A-obj.mult(U_qr,s_qr,V_qr)).norm()<<std::endl;
    exportmatrix(U_qr,"U_qr.txt"); 
    exportmatrix(s_qr.transpose(),"s_qr.txt");
    exportmatrix(V_qr.transpose(),"Vt_qr.txt");*/
    

    int r=15, p=0, q=1;
    start = std::chrono::high_resolution_clock::now();
    auto [U_rsvd,s_rsvd,V_rsvd] = obj.rsvd(A,r,p,q);
    end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration_rsvd = end  - start;
    std::cout<<"\nrSVD ( r = "<< r <<" ):\n";
    std::cout<<"Time of execution rSVD algorithm: "<< duration_rsvd.count() << " s" << std::endl;
    std::cout<<"Norm of A - U * S * Vt = "<<(A-obj.mult_parallel(U_rsvd,s_rsvd,V_rsvd)).norm()<<std::endl;

    double SpeedUp = duration_pm.count()/duration_rsvd.count();
    std::cout << "Speed Up randomized: "<< SpeedUp << std::endl;
    std::cout << "Norm of difference eigenvalues = " << (s_rsvd - s_pm.head(s_rsvd.rows())).norm()/(r+p) << std::endl;

    exportmatrix(U_rsvd,"U_rsvd.txt");
    exportmatrix(s_rsvd.transpose(),"s_rsvd.txt");
    exportmatrix(V_rsvd.transpose(),"Vt_rsvd.txt");

    return 0;
}