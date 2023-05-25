#include <iostream>
#include <math.h>
#pragma once
#include "print_mat.h"
#include "matrix_operation.h"
using namespace std;


void Update_H(double **W, double **H, double **V, int row, int k, int col)
{
    double *transpose_W[N], *numerator[N],*transpose_WV[N];

    // allocating transpose_W and numerator
    for (int i = 0; i < k; i++)
        transpose_W[i] = (double *)malloc(row * sizeof(double));

    for (int i = 0; i < k; i++)
        numerator[i] = (double *)malloc(col * sizeof(double));

    for(int i=0;i<k;i++)
        transpose_WV[i]=(double *)malloc(col * sizeof(double));


    transpose(W, transpose_W, row, k); // WT

    multiply(transpose_WV, transpose_W, V, k, row, col);// WT*V

    double *Hd[N],*AdjM[N],*degree[N],*lemdaHd[N];

    for(int i=0;i<k;i++)
        Hd[i]=(double *)malloc(row * sizeof(double));

    for (int i = 0; i < row; i++)
        AdjM[i] = (double *)malloc(row * sizeof(double));

    for(int i=0;i<k;i++)
        lemdaHd[i]=(double *)malloc(row * sizeof(double));

    AdjacentMatrix(AdjM,V,row,col);
    //cout<<"adjacent matrix is"<<endl;
    //print_matrix(AdjM,row,col);

    for (int i = 0; i < row; i++)
        degree[i] = (double *)malloc(row * sizeof(double));

    DegreeMatrix(degree,AdjM,row,col);
    //cout<<"degree matrix is"<<endl;
    //print_matrix(degree,row,col);

    multiply(Hd,H,degree,k,col,row);//Hd
    constMultiplication(lemdaHd,Hd,lemda,k,col);
    add_element_wise(numerator,transpose_WV,lemdaHd,k,col);

    //denominator part

     double *den_part1[N], *denominator[N],*transpose_WWH[N],*HA[N],*lemdaHA[N];

    for (int i = 0; i < k; i++)
        den_part1[i] = (double *)malloc(k * sizeof(double));

     for (int i = 0; i < k; i++)
        transpose_WWH[i] = (double *)malloc(col * sizeof(double));

    for (int i = 0; i < k; i++)
        denominator[i] = (double *)malloc(col * sizeof(double));

     for (int i = 0; i < k; i++)
        HA[i] = (double *)malloc(col * sizeof(double));

     for (int i = 0; i < k; i++)
        lemdaHA[i] = (double *)malloc(col * sizeof(double));

    multiply(den_part1, transpose_W, W, k, row, k); // WT*W
    multiply(transpose_WWH, den_part1, H, k, k, col);//(WT*W)*H

    multiply(HA,H,AdjM,k,col,col);
    constMultiplication(lemdaHA,HA,lemda,k,col);
    add_element_wise(denominator,transpose_WWH,lemdaHA,k,col);


    double *updated_H[N];    

    for (int i = 0; i < k; i++)
        updated_H[i] = (double *)malloc(col * sizeof(double));
    divide_element_wise(updated_H, numerator, denominator, k, col);

    double *ans_H[N];

    for (int i = 0; i < k; i++)
        ans_H[i] = (double *)malloc(col * sizeof(double));
    multiply_element_wise(ans_H, H, updated_H, k, col);

    copy_matrix(ans_H, H, k, col);
    free_matrix(ans_H, k);
    free_matrix(transpose_W, k);
    free_matrix(numerator, k);
    free_matrix(transpose_WV, k);
    free_matrix(lemdaHd,k);

    free_matrix(Hd, k);
    free_matrix(AdjM, row);
    free_matrix(degree, row);
    free_matrix(den_part1, k);
    free_matrix(denominator, k);
    free_matrix(transpose_WWH, k);

    free_matrix(HA, k);
    free_matrix(updated_H, k);
    free_matrix(lemdaHA, k);
}

void Update_W(double **W, double **H, double **V, int row, int k, int col)
{
    double *HT[N], *numerator[N];

    for (int i = 0; i < col; i++)
        HT[i] = (double *)malloc(k * sizeof(double));


    for (int i = 0; i < row; i++)
        numerator[i] = (double *)malloc(k * sizeof(double));


    //cout<<"D matrix is:"<<endl;
    //print_matrix(D,row,col);

    transpose(H, HT, k, col); // HT
    multiply(numerator, V, HT, row, col, k); // V*HT

    double *H_transpose_H[N];



     for (int i = 0; i < k; i++)
        H_transpose_H[i] = (double *)malloc( k * sizeof(double));
      // HD
     multiply(H_transpose_H, H, HT, k, col, k); // H_transpose_H

    double *denominator[N];

    for (int i = 0; i < row; i++)
        denominator[i] = (double *)malloc(k * sizeof(double));


    multiply(denominator, W, H_transpose_H, row, k, k);

    double *updated_W[N];

    for (int i = 0; i < row; i++)
        updated_W[i] = (double *)malloc(k * sizeof(double));

    divide_element_wise(updated_W, numerator, denominator, row, k);

    double *ans_W[N];

    for (int i = 0; i < row; i++)
        ans_W[i] = (double *)malloc(k * sizeof(double));

    multiply_element_wise(ans_W, W, updated_W, row, k);

    copy_matrix(ans_W, W, row, k);
    free_matrix(ans_W, row);
    free_matrix(updated_W, row);
    free_matrix(denominator, row);
    free_matrix(HT, col);
    free_matrix(numerator, row);

    free_matrix(H_transpose_H, k);
}

void G_NMF()
{

    double *matrix[N];
    int row, col, i, j, k;

     char choice;
     std::srand(std::time(0));

    cout << "Enter 'M' for manual input or 'F' to read from a file: ";
    cin >> choice;

    if (choice == 'M' || choice == 'm') {

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

    } else if (choice == 'F' || choice == 'f') {
        ifstream inputFile("matrix.txt");

        if (!inputFile.is_open()) {
            cout << "Failed to open the file." << endl;
            return;
        }

        // Read the number of rows and columns
        inputFile >> row >> col;

        // Allocate memory for the matrix
        for (i = 0; i < row; i++) {
            matrix[i] = new double[col];
        }

        // Read the matrix elements
        for (i = 0; i < row; i++) {
            for (j = 0; j < col; j++) {
                inputFile >> matrix[i][j];
            }
        }

        inputFile.close();
    } else {
        cout << "Invalid choice. Exiting." << endl;
        return;
    }

    // Rest of the code remains unchanged
    // ...

    // Deallocate memory for the matrix
    //print_matrix(matrix, row, col);
    //  Dimension of broken matrix
    cout<<"Enter the dimension in which you want to break:";
    cin >> k;
    double *W[N], *H[N]; // broken down in m*k and k*n matrix

    cout << "Enter 'M' to manually input W matrix or 'F' to read from a file: ";
    char ch0;
    cin>>ch0;

    if (ch0 == 'M' || ch0 == 'm') {
        cout << "Enter W matrix:" << endl;

        for (i = 0; i < row; i++)
            W[i] = new double[k];

        for (i = 0; i < row; i++) {
            for (j = 0; j < k; j++) {
                cin >> W[i][j];
            }
        }
    }

    else if (ch0 == 'F' || ch0 == 'f') {
        ifstream wFile("W.txt");

        if (!wFile.is_open()) {
            cout << "Failed to open the W matrix file." << endl;
            return;
        }

        // Read the number of rows and columns for W matrix
        int wRow, wCol;
        wFile >> wRow >> wCol;

        if (wRow != row || wCol != k) {
            cout << "Invalid dimensions for W matrix. Exiting." << endl;
            return;
        }

        for (i = 0; i < row; i++)
            W[i] = new double[k];

        // Read the W matrix elements
        for (i = 0; i < row; i++) {
            for (j = 0; j < k; j++) {
                 wFile >> W[i][j];
            }
        }

        wFile.close();

    }
    else {
        cout << "Invalid choice. Exiting." << endl;
        return;
    }
    cout << "Enter 'M' to manually input W matrix or 'F' to read from a file: ";
    char ch1;
    cin>>ch1;
    
    if (ch1 == 'M' || ch1 == 'm') {
        cout << "Enter H matrix:" << endl;

        for (i = 0; i < k; i++)
            H[i] = new double[col];

        for (i = 0; i < k; i++) {
            for (j = 0; j < col; j++) {
                cin >> H[i][j];
            }
        }

    } else if (ch1== 'F' || ch1 == 'f') {
        ifstream hFile("H.txt");

        if (!hFile.is_open()) {
            cout << "Failed to open the H matrix file." << endl;
            return;
        }

        // Read the number of rows and columns for H matrix
        int hRow, hCol;
        hFile >> hRow >> hCol;

        if (hRow != k || hCol != col) {
            cout << "Invalid dimensions for H matrix. Exiting." << endl;
            return;
        }

        for (i = 0; i < k; i++)
            H[i] = new double[col];

        // Read the H matrix elements
        for (i = 0; i < k; i++) {
            for (j = 0; j < col; j++) {
                hFile >> H[i][j];
            }
        }

        hFile.close();

    } else {
        cout << "Invalid choice. Exiting." << endl;
        return;
    }

    //cout<<"Show initial two matrices";
    //print_two_matrix(W, H, row, k, col);
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
    double initial_cost = cost;
    double prev_cost = 0;
    printf("Updating costs:\n ");
    while (cost > EPSILON)
    {

        if ((counter % 2) == 0)
        {
            Update_H(W, H, matrix, row, k, col);
            //cout << "New H" << endl;
        }
        else
        {
            Update_W(W, H, matrix, row, k, col);
           // cout << "New W" << endl;
            //print_matrix(W,row,k);
        }
        counter++;
        multiply(V, W, H, row, k, col);
        cost = cost_function(matrix, V, row, col);
        if (fabs(prev_cost - cost) <= EPSILON)
        {
            printf("Reached relative minima\n");
            break;
        }
        else
        {
            prev_cost = cost;
        }

        // local minima reached need to stop by calculating difference with previous error
    }
    printf("Factorization done!");
    printf("The beginning cost was: %lf\n", initial_cost);
    printf("The final cost was: %lf\n", cost);
    printf("Total number of iterations before arriving at result: %d\n", counter);
    printf("The broken down matrix:\n ");
    print_two_matrix(W, H, row, k, col);


}
