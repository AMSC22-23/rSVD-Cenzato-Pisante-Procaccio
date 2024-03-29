// g++ -I${mkEigenInc} cancer.cpp svd.cpp QR_Decomposition_parallel.cpp -o cancer -Wall -DEIGEN -fopenmp -DPARALLEL
// or to test with our matrix class
// g++ cancer.cpp svd.cpp QR_Decomposition_parallel.cpp -o cancer -Wall -fopenmp -DPARALLEL

#include "applications.hpp"

Matrix readCSV(const std::string& filename);
std::vector<unsigned int> readCSV_grp(const std::string& filename);
void exportmatrix(Matrix A, std::string outputFileName);

int main()
{
    std::string ovariancancer_obs_path  = "../data/ovariancancer_obs.csv";
    Matrix A = readCSV(ovariancancer_obs_path );    // 216 x 4000

    //std::string ovariancancer_grp_path  = "../python/data/ovariancancer_grp.csv";
    //std::vector<unsigned int> labels = readCSV_grp(ovariancancer_grp_path);

    APPLICATIONS obj;
    Vector explained_variance = obj.explained_variance(A);
    std::cout << "Explained variance 1-st eigenvalue = " << explained_variance(0,0) << std::endl;

    auto start = std::chrono::high_resolution_clock::now();
    Matrix T = obj.pca(A, 20);
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = end - start;
    exportmatrix(T,"T.txt");
    std::cout << "Time of execution PCA = " << duration.count() << " s" << std::endl;

    return 0;
}

Matrix readCSV(const std::string& filename) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error opening file: " << filename << std::endl;
        return Matrix();
    }
    else{
        std::cout << "Opening file: "<< filename << std::endl;
    }

    std::vector<std::vector<double>> data;
    std::string line;

    while (std::getline(file, line)) {
        std::istringstream iss(line);
        std::vector<double> row;
        std::string value;

        while (std::getline(iss, value, ',')) {
            row.push_back(std::stod(value));
        }

        data.push_back(row);
    }

    // Determine the dimensions of the matrix
    int rows = data.size();
    int cols = (rows > 0) ? data[0].size() : 0;

    // Create a matrix and fill it with data
    Matrix result(rows, cols);

    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            result(i, j) = data[i][j];
        }
    }

    return result;
}

std::vector<unsigned int> readCSV_grp(const std::string& filename) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error opening file: " << filename << std::endl;
        return std::vector<unsigned int>();
    }
    else{
        std::cout << "Opening file: "<< filename << std::endl;
    }

    std::vector<unsigned int> data;
    std::string line;

    while (std::getline(file, line)) {
        if(!line.empty()){      
            if(line == "Cancer"){
                data.push_back(1);
            }
            else{
                data.push_back(0);
            }
        }
    }

    return data;
}

void exportmatrix(Matrix A, std::string outputFileName)
{
    // Write the matrix to the file
    std::ofstream outputFile(outputFileName);
    if (outputFile.is_open())
    {
		outputFile<<A;
        std::cout << "Computed matrix has been written to " << outputFileName << std::endl;
        // Close the file
        outputFile.close();
    }
    else
    {
        std::cerr << "Error opening file for writing." << std::endl;
    }
}
