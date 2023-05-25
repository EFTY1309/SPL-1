#include <iostream>
#pragma once
#include <math.h>
#include <stdlib.h>
#define STANDARD_EPSILON 2.2204460492503130808472633361816E-16
#define EPSILON 0.1
const int N = 10000;
const int alpha=0.01;
const int beta=0.1;
const int gamma=10;
const int lemda=0.0000;
const int INT_MAX1=1e5;

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

void AdjacentMatrix(double **Adj_M,double **X,int row,int col)
{
    double res;
    for(int i=0;i<row;i++)
        for(int j=0;j<row;j++)
    {
        Adj_M[i][j]=0;
    }

    for(int i=0;i<row;i++)
    {
        for(int j=0;j<col;j++)
        {
            if(i==j)
            {
             continue;
             if(j== row-1) break;

            }
            res=0;

            for(int k=0;k<col;k++)
            {
              res+=pow((X[i][k]-X[j][k]),2);
            }
            Adj_M[i][j]=sqrt(res);
        }
    }
}

void DegreeMatrix(double **degree,double **AdjM,int row,int col)
{
    for(int i=0;i<row;i++)
        for(int j=0;j<row;j++)
    {
        degree[i][j]=0;
    }

    for(int i=0;i<row;i++){
        double min_element = INT_MAX;
        for(int j=0;j<row;j++)
    {
       if (AdjM[i][j] < min_element and AdjM[i][j] !=0) {
        min_element = AdjM[i][j];

      }
    }
    degree[i][i]=min_element;

    }
}


void FindD(double **D,double **X,double **ab_D,int row,int col)
{
    double a1[N];
    for(int i=0;i<col;i++)
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


    for(int i=0;i<col;i++)
    {
       for(int j=0;j<col;j++)
        {
            if(i==j) D[i][j]=1/a1[j];
        }
    }

}

void add_element_wise(double **V, double **m1, double **m2, int m, int n)
{
    for (int i = 0; i < m; i++)
    {
        for (int j = 0; j < n; j++)
        {
            V[i][j] = m1[i][j]+m2[i][j];
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

void FindG(double **G,double **X,double **ab_G,int row,int col)
{
    double a1[N];
    for(int i=0;i<row;i++)
        for(int j=0;j<col;j++)
    {
        G[i][j]=0;
    }

    for(int j=0;j<col;j++)
    {
        double power=0;
        for(int i=0;i<row;i++)
        {
            power+=pow((X[i][j]-ab_G[i][j] ),2);
        }
        a1[j]=sqrt(power);
    }

    for(int i=0;i<row;i++)
    {
       for(int j=0;j<col;j++)
        {
            if(i==j) G[i][j]=0.5*a1[j];
        }
    }

}

void FindP(double **P,double **W,int row,int col)
{
    double a1[N];
    for(int i=0;i<col;i++)
        for(int j=0;j<col;j++)
    {
        P[i][j]=0;
    }

    for(int j=0;j<col;j++)
    {
        double power=0;
        for(int i=0;i<row;i++)
        {
            power+=pow((W[i][j] ),2);
        }
        a1[j]=sqrt(power);
    }


    for(int i=0;i<col;i++)
    {
       for(int j=0;j<col;j++)
        {
            if(i==j) P[i][j]=0.5*a1[j];
        }
    }

}

void constMultiplication(double **X,double **P,double cons,int row,int col)
{
  for(int i=0;i<row;i++)
        for(int j=0;j<col;j++)
  {
      X[i][j]=cons*P[i][j];
  }
}
void Ecalculation(double**E,int row,int col)
{

    for(int i=0;i<row;i++)
        for(int j=0;j<row;j++) E[i][j]=0;

    for(int i=0;i<row;i++)
        for(int j=0;j<row;j++)
        {
            if(i==j) E[i][j]=1;
        }
}




