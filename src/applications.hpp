#include "svd.hpp"
#include "include/stb_image.h"
#include "include/stb_image_write.h"

class APPLICATIONS
{
public:
    APPLICATIONS()
    {}

    /* Reduction of dimensionality of the matrix A:
        It returns a matrix with the first r (+ 5 = oversampling parameter)
        principal components axes.*/
    Matrix pca(const Matrix &A, const int r)
    {
        Matrix X = A;
        size_t m = X.rows(), n = X.cols();
        SVD obj;

// Center X
#pragma omp parallel for
        for (size_t i = 0; i < m; i++)
        {
            double media = 0.;
            for (size_t j = 0; j < n; j++)
            {
                media += X(i, j);
            }
            media /= n;
            for (size_t j = 0; j < n; j++)
            {
                X(i, j) -= media;
            }
        }

// Compute rSVD
        auto [U, s, V] = obj.rsvd(X, r);


// Compute Principal Components Matrix T.
/*#ifdef EIGEN
        return U * s.asDiagonal();
#else
        Matrix T(m, s.rows());
        #pragma omp parallel for
        for(size_t i = 0; i<s.rows(); i++)
            for(size_t j = 0; j<X.rows(); j++)
                T(j,i) = s(i,0) * U(j,i);
        return T;
#endif*/
    return U.transpose() * X;
    }

    Matrix image_compression(const stbi_uc *R, int channels, int channel, int Heigth, int Width, int r, int p)
    {
        Matrix X;
        X = ExctractComponentLuminosity(R, channels, channel, Heigth, Width);

        SVD obj;
        auto [U, s, V] = obj.rsvd(X, r, p, 0);
        return obj.mult_SVD(U, s, V);
    }

    Matrix ExctractComponentLuminosity(const stbi_uc *image, int channels, int colour, int Heigth, int Width)
    {
        Matrix RgbComponent(Heigth, Width);
#pragma omp parallel for collapse(2)
        for (int i = 0; i < Heigth; i++)
        {
            for (int j = 0; j < Width; j++)
            {
                int index = (i * Width + j) * channels;

                /**
                 * Exctraction of luminosity for each component
                 */
                unsigned char rgb = image[index + colour];
                RgbComponent(i, j) = static_cast<double>(rgb) / 255.;
            }
        }

        return RgbComponent;
    }

    void backward_conversion(stbi_uc *image, const Matrix CompressedComponent, int channels, int channel, int Height, int Width)
    {
#pragma omp parallel for collapse(2)
        for (int i = 0; i < Height; ++i)
        {
            for (int j = 0; j < Width; ++j)
            {
                image[(i * Width + j) * channels + channel] = static_cast<stbi_uc>(std::min(255.0, std::max(0.0, CompressedComponent(i, j) * 255.0)));
            }
        }
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