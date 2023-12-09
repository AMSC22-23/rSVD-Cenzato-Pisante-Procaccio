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
                outputFile << std::setw(8) << std::fixed << std::setprecision(2) << A(i,j) << " ";
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
    std::ifstream file("test_matrices/matrix2.txt");           //file with dim and then matrix
    if (!file.is_open()) {
        std::cerr << "Error opening file." << std::endl;
        return 1;
    }

    int m, n;
    file >> m >> n;
    Matrix A(m,n);
    // Read the matrix data
    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < n; ++j) {
            file >> A(i,j);
        }
    }
    file.close();

    double m_eps=1e-13;        
    SVD obj(m_eps);

    //std::tie(U,s,V) = obj.svd_with_PM(A);
    auto [Up,sp,Vp] = obj.svd_with_PM(A);
    std::cout<<"\nReduced SVD with Power Method:\n";
    std::cout<<"Norm of A - U * S * Vt = "<<(A-obj.mult(Up,sp,Vp)).norm()<<std::endl;
    exportmatrix(Up,"U_pm.txt");
    exportmatrix(sp.transpose(),"s_pm.txt");
    exportmatrix(Vp.transpose(),"Vt_pm.txt");

    /*std::cout<<Up<<std::endl;
    std::cout<<sp<<std::endl;
    std::cout<<Vp<<std::endl;*/

    /*std::cout<<"\nPseudo-inverse :\n";
    std::cout<<obj.pseudoinverse(A)<<std::endl;*/

    /*std::tie(U,s,V) = obj.svd_with_qr(A);
    std::cout<<"\nSVD with QR :\n";
    std::cout<<"Norm of A - U * S * Vt = "<<(A-obj.mult(U,s,V)).norm()<<std::endl;*/
    //std::cout<<U<<std::endl;
    //std::cout<<s<<std::endl;
    //std::cout<<V<<std::endl; 

    int r=3, p=5, q=1;
    auto [U,s,V] = obj.rsvd(A,r,p,q);
    std::cout<<"\nrSVD :\n";
    std::cout<<"Norm of A - U * S * Vt = "<<(A-obj.mult(U,s,V)).norm()<<std::endl;
    /*std::cout<<U<<std::endl;
    std::cout<<s<<std::endl;
    std::cout<<V<<std::endl;*/
    exportmatrix(U,"U_rsvd.txt");
    exportmatrix(s.transpose(),"s_rsvd.txt");
    exportmatrix(V.transpose(),"Vt_rsvd.txt");

    std::cout<<"\nComparison between power method and rSVD :\n";
    std::cout<<"Norm of difference U = "<<(Up-U).norm()<<std::endl;
    std::cout<<"Norm of difference s = "<<(sp-s).norm()<<std::endl;
    std::cout<<"Norm of difference V = "<<(Vp-V).norm()<<std::endl;
    //colonne di U e V vengono con stessi numeri in entrambi i metodi ma a volte hanno senso opposto!

    return 0;
}