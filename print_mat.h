#include <iostream>
#pragma once
#include <fstream>
#include <iomanip>
using namespace std;

void print_two_matrix(double **m1, double **m2, int m, int k, int n)
{
      std::ofstream file("result.txt");
    if (!file) {
        std::cout << "Unable to open the file." << std::endl;
        return;
    }
// Set output format to store two digits after the decimal point
    file << std::fixed << std::setprecision(2);

    // Write matrices to the file
    file<<"First matrix"<<std::endl;
    file << m << " " << k << std::endl;
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < k; j++) {
            file << m1[i][j] << " ";
        }
        file << std::endl;
    }

    file << std::endl; // Separate the two matrices

    file<<"Second matrix"<<std::endl;
    file << k << " " << n << std::endl;
    for (int i = 0; i < k; i++) {
        for (int j = 0; j < n; j++) {
            file << m2[i][j] << " ";
        }
        file << std::endl;
    }
    
    file.close();
    std::cout << "Matrices stored successfully in result.txt." << std::endl;

}
// when we print a declared 2d array
void print_matrix(double **V, int m, int n)
{
    for (int i = 0; i < m; i++)
    {
        for (int j = 0; j < n; j++)
        {
            printf("%.2lf\t", V[i][j]);
        }
        printf("\n");
    }
    cout << endl;
}


