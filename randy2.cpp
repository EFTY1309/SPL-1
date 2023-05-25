#include <iostream>
#include <vector>
#include <cstdlib>
#include <ctime>
#include <fstream>

int main() {
    // Set the seed for the random number generator
    std::srand(std::time(0));

    int rows, columns;
    std::cout << "Enter the number of rows: ";
    std::cin >> rows;
    std::cout << "Enter the number of columns: ";
    std::cin >> columns;

    // Create a 2D vector to store the matrix
    std::vector<std::vector<int>> matrix(rows, std::vector<int>(columns));

    // Generate random non-negative elements greater than 100 and less than 1000
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < columns; ++j) {
            matrix[i][j] = std::rand() % 100;  // Generates a random number between 101 and 999
        }
    }

    // Open the output file
    std::ofstream outputFile("H.txt");

    // Check if the file was successfully opened
    if (!outputFile) {
        std::cout << "Failed to open the file." << std::endl;
        return 1;
    }

    // Print the matrix to the file
    outputFile << "Matrix:" << std::endl;
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < columns; ++j) {
            outputFile << matrix[i][j] << " ";
        }
        outputFile << std::endl;
    }

    // Close the file
    outputFile.close();

    std::cout << "Output is stored in the file 'spl.txt'." << std::endl;

    return 0;
}
