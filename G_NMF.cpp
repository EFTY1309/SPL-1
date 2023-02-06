#include <iostream>
#include <math.h>
#include "print_matrix.h"
#include "matrix_operation.h"
#include "matrix_factorizations.h"
using namespace std;

void G_NMF()
{

double *matrix[N],*AdjM[N],*degree[N];
    int row, col, i, j, k;

    cout<<"Enter the number of row and column for matrix:"<<endl;
    // Taking input
    cin >> row >> col; // m*n matrix

    for (i = 0; i < row; i++)
        matrix[i] = (double *)malloc(col * sizeof(double));

    cout<<"Enter matrix:"<<endl;
    for (i = 0; i < row; i++)
    {
        for (j = 0; j < col; j++)
            cin >> matrix[i][j];
    }

    for (i = 0; i < row; i++)
        AdjM[i] = (double *)malloc(col * sizeof(double));

    AdjacentMatrix(AdjM,matrix,row,col);

    for (i = 0; i < row; i++)
        degree[i] = (double *)malloc(col * sizeof(double));

    DegreeMatrix(degree,AdjM,row,col);

    print_matrix(matrix, row, col);
    cout<<"Adjacent Matrix is:"<<endl;
    print_matrix(AdjM, row, col);

    cout<<"Degree Matrix is:"<<endl;
    print_matrix(degree, row, col);

    // Generating two matrices using user input dimension and random number generator
    //  Dimension of broken matrix
    cout<<"Enter the dimension in which you want to break:";
    cin >> k;
    double *W[N], *H[N]; // broken down in m*k and k*n matrix

    for (i = 0; i < row; i++)
        W[i] = (double *)malloc(k * sizeof(double));
    for (i = 0; i < k; i++)
        H[i] = (double *)malloc(col * sizeof(double));

    cout<<"Enter W matrix:"<<endl;
    for (i = 0; i < row; i++)
    {
        for (j = 0; j < k; j++)
        {
            cin>>W[i][j];
        }
    }

    cout<<"Enter H matrix:"<<endl;
    for (i = 0; i < k; i++)
    {
        for (j = 0; j < col; j++){

            cin>>H[i][j];
        }
    }

    cout<<"Show initial two matrices";
    print_two_matrix(W, H, row, k, col);
    // multiplication
    double *V[N];
    for (i = 0; i < row; i++)
        V[i] = (double *)malloc(col * sizeof(double));
    multiply(V, W, H, row, k, col);

    //printf("V= \n");

    int counter = 1;
    printf("Initial cost: ");
    // cost function
    double cost = cost_function(matrix, V, row, col);



}
