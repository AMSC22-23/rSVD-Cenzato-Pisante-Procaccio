#include "svd.hpp"

//g++ -I${mkEigenInc} svd_test.cpp svd.cpp QR_Decomposition.cpp -o prova

#include <fstream>
#include <sstream> 
#include <iomanip>

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

int main(){
    std::ifstream file("test_matrices/matrix4.txt");           //file with dim and then matrix
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

    double m_eps=1e-13;        
    SVD obj(m_eps);
    
    
    auto [U_pm,s_pm,V_pm] = obj.svd_with_PM(A);
    std::cout<<"\nReduced SVD with Power Method:\n";
    std::cout<<"Norm of A - U * S * Vt = "<<(A-obj.mult(U_pm,s_pm,V_pm)).norm()<<std::endl;
    exportmatrix(U_pm,"U_pm.txt");
    exportmatrix(s_pm.transpose(),"s_pm.txt");
    exportmatrix(V_pm.transpose(),"Vt_pm.txt");

    /*std::cout<<"\nPseudo-inverse :\n";
    std::cout<<obj.pseudoinverse(A)<<std::endl;*/

    /*auto [U_qr,s_qr,V_qr] = obj.svd_with_qr(A);
    std::cout<<"\nSVD with QR :\n";
    std::cout<<"Norm of A - U * S * Vt = "<<(A-obj.mult(U_qr,s_qr,V_qr)).norm()<<std::endl;
    exportmatrix(U_qr,"U_qr.txt"); 
    exportmatrix(s_qr.transpose(),"s_qr.txt");
    exportmatrix(V_qr.transpose(),"Vt_qr.txt");*/
    

    /*int r=40, p=5, q=1;
    auto [U_rsvd,s_rsvd,V_rsvd] = obj.rsvd(A,r,p,q);
    std::cout<<"\nrSVD :\n";
    std::cout<<"Norm of A - U * S * Vt = "<<(A-obj.mult(U_rsvd,s_rsvd,V_rsvd)).norm()<<std::endl;
    exportmatrix(U_rsvd,"U_rsvd.txt");
    exportmatrix(s_rsvd.transpose(),"s_rsvd.txt");
    exportmatrix(V_rsvd.transpose(),"Vt_rsvd.txt");*/

    /*std::cout<<"\nComparison between power method and rSVD :\n";
    std::cout<<"Norm of difference U = "<<(U_pm-U_rsvd).norm()<<std::endl;
    std::cout<<"Norm of difference s = "<<(s_pm-s_rsvd).norm()<<std::endl;
    std::cout<<"Norm of difference V = "<<(V_pm-V_rsvd).norm()<<std::endl;*/

    /*std::cout<<"\nComparison between power method and qr algorithm :\n";
    std::cout<<"Norm of difference U = "<<(U_pm-U_qr).norm()<<std::endl;
    std::cout<<"Norm of difference s = "<<(s_pm-s_qr).norm()<<std::endl;
    std::cout<<"Norm of difference V = "<<(V_pm-V_qr).norm()<<std::endl;*/

    return 0;
}