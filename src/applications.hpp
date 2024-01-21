#include "svd.hpp"

class APPLICATIONS
{
public:
    APPLICATIONS()
    {
    }

    Matrix pca(const Matrix &A, const int r)
    {
        Matrix X = A;
        Matrix C(X.cols(), X.cols());
        SVD obj;
        C = obj.preprocess(X);
        auto [U, s, V] = obj.rsvd(C, r);

        return U.transpose() * X;
    }

    Matrix image_compression(const stbi_uc *R, int channels, int Heigth, int Width)
    {
        Matrix X;
        X = grayscale(R, channels, Heigth, Width);

        SVD obj;
        int r = 5;
        int p = 10;
        auto [U, s, V] = obj.rsvd(X, r, p);
        return obj.mult_SVD(U, s, V);
    }

    Matrix grayscale(const stbi_uc *image, int channels, int Heigth, int Width)
    {
        Matrix grayscaleImage(Heigth, Width);
        #pragma omp parallel for collapse(2)
        for (int i = 0; i < Heigth; i++)
        {
            for (int j = 0; j < Width; j++)
            {
                int index = (i * Width + j) * channels;

                /**
                 * Extraction of color components from an image in RGB format
                 * We store the three components in an array
                 */
                unsigned char red = image[index];
                unsigned char green = image[index + 1];
                unsigned char blue = image[index + 2];

                /**
                    Set the grayscale value by using the luminosity formula for each channel and normalize by the max value = 255
                */
                // grayscaleImage(i, j) = static_cast<double>(0.21 * red + 0.72 * green + 0.07 * blue) / 255.;

                /**
                 * AVG method
                 */
                grayscaleImage(i, j) = static_cast<double>((red + green + blue) / 3) / 255.;

                /**
                 * LIGHTNESS method
                 */
                //grayscaleImage(i, j) = static_cast<double>((std::max({red, green, blue}) + std::min({red, green, blue})) / 2) / 255.;
            }
        }

        return grayscaleImage;
    }

    void backward_conversion(stbi_uc *image, Matrix Compressed, int channels, int Height, int Width)
    {
        // Convert the Eigen matrix back to image data
        #pragma omp parallel for collapse(2)
        for (int i = 0; i < Height; ++i)
        {
            for (int j = 0; j < Width; ++j)
            {
                // Clip values to ensure they are within the valid range [0, 255]
                image[(i * Width + j) * channels] = static_cast<stbi_uc>(std::min(255.0, std::max(0.0, Compressed(i, j) * 255.0)));
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