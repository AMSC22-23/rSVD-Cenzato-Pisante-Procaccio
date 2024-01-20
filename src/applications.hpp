#include "svd.hpp"

class APPLICATIONS
{
    public:
    APPLICATIONS()
    {}

    Matrix pca(const Matrix &A, const int r){
        Matrix X=A;
        Matrix C(X.cols(),X.cols());
        SVD obj;
        C = obj.preprocess(X);
        auto[U,s,V] = obj.rsvd(C,r);

        return U.transpose() * X;
    }

    Matrix image_compression(const Matrix &R, const Matrix &G, const Matrix &B, unsigned int r){
        Matrix X;
        X = grayscale(R,G,B);
        SVD obj;
        auto [U,s,V] = obj.rsvd(X,r,0);
        return obj.mult_SVD(U,s,V);
    }

    Matrix grayscale(const Matrix &R, const Matrix &G, const Matrix &B){
        size_t m = R.rows(), n = R.cols();
        Matrix X(m,n);
        for(size_t i=0;i<m;i++){
            for(size_t j=0;j<n;j++){
                X(i,j) = (R(i,j) + G(i,j) + B(i,j))/3;                
            }
        }
        return X;
    }

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
};