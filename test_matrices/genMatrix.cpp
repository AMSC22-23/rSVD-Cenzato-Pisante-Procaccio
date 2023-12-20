#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <random>

int main() {
    // Specify the dimensions of the matrix
    int rows = 6;
    int cols = 5;

    // Generate a random matrix
    std::vector<std::vector<double>> matrix(rows, std::vector<double>(cols));
    
    // Set up a random number generator
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> dist(1.0, 100.0);  // can adjust the range

    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            matrix[i][j] = dist(gen);
        }
    }

    // Write the matrix to a file
    std::ofstream outputFile("matrix.txt");
    if (outputFile.is_open()) {
        // Write dimensions to the first row
        outputFile << rows << " " << cols << std::endl;

        // Write matrix data
        for (int i = 0; i < rows; ++i) {
            for (int j = 0; j < cols; ++j) {
                outputFile << std::setw(8) << std::fixed << std::setprecision(2) << matrix[i][j] << " ";
            }
            outputFile << std::endl;
        }

        std::cout << "Random matrix has been written to matrix.txt" << std::endl;

        // Close the file
        outputFile.close();
    } else {
        std::cerr << "Error opening file for writing." << std::endl;
    }

    return 0;
}