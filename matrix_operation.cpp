#include <iostream>
#include <math.h>
#include "matrix_operation.h"
using namespace std;

void multiply(double **V, double **m1, double **m2, int m, int k, int n)
{
    for (int i = 0; i < m; i++)
    {
        for (int j = 0; j < n; j++)
        {
            V[i][j] = 0;
            for (int x = 0; x < k; x++)
            {
                V[i][j] += m1[i][x] * m2[x][j];
            }
        }
    }
}


void FindD(double **D,double **X,double **ab_D,int row,int col)
{
    double a1[N];
    for(int i=0;i<row;i++)
        for(int j=0;j<col;j++)
    {
        D[i][j]=0;
    }

    for(int j=0;j<col;j++)
    {
        double power=0;
        for(int i=0;i<row;i++)
        {
            power+=pow((X[i][j]-ab_D[i][j] ),2);
        }
        a1[j]=sqrt(power);
    }


    for(int i=0;i<row;i++)
    {
       for(int j=0;j<col;j++)
        {
            if(i==j) D[i][j]=1/a1[j];
        }
    }

}

void multiply_element_wise(double **V, double **m1, double **m2, int m, int n)
{
    for (int i = 0; i < m; i++)
    {
        for (int j = 0; j < n; j++)
        {
            V[i][j] = m1[i][j] * m2[i][j];
        }
    }
}
void divide_element_wise(double **V, double **m1, double **m2, int m, int n) // divide element wise
{
    for (int i = 0; i < m; i++)
    {
        for (int j = 0; j < n; j++)
        {
            V[i][j] = m1[i][j] / m2[i][j];
        }
    }
}

void copy_matrix(double **from, double **to, int row, int col)
{
    for (int i = 0; i < row; i++)
    {
        for (int j = 0; j < col; j++)
        {
            to[i][j] = from[i][j];
        }
    }
}

void transpose(double **input_matrix, double **transpose_matrix, int row, int col)
{
    for (int i = 0; i < row; i++)
    {
        for (int j = 0; j < col; j++)
        {
            transpose_matrix[j][i] = input_matrix[i][j];
        }
    }

}


double cost_function(double **initial_matrix, double **current, int row, int col)
{
    double cost = 0.0; // epsilon
    double sum = 0.0;

    for (int i = 0; i < row; i++)
    {
        for (int j = 0; j < col; j++)
        {
            // square of the difference between initial and final for each index
            sum += pow((initial_matrix[i][j] - current[i][j]), 2);
        }
    }
    // root of the sum of squares epsilon
    cost = sqrt(sum);

    cout << "Cost is :" << cost << endl;

    return cost;
}

void free_matrix(double **matrix, int row)
{
    for (int i = 0; i < row; i++)

    {
        free(matrix[i]);
    }
}


