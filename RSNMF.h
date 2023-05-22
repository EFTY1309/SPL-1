#include <iostream>
#pragma once
#include <math.h>
#include "print_mat.h"
#include "matrix_operation.h"
using namespace std;

void update_H1(double **W, double **H, double **V, int row, int k, int col)
{

   double *transpose_W[N], *numerator[N],*transpose_WV[N],*transpose_WVG[N],*G[N],*WH[N];
    cout<<"jai na";
    // allocating transpose_W and numerator
    for (int i = 0; i < k; i++)
        transpose_W[i] = (double *)malloc(row * sizeof(double));

    for (int i = 0; i < k; i++)
        numerator[i] = (double *)malloc(col * sizeof(double));

    for(int i=0;i<k;i++)
        transpose_WV[i]=(double *)malloc(col * sizeof(double));

     for(int i=0;i<k;i++)
        transpose_WVG[i]=(double *)malloc(col * sizeof(double));

    for(int i=0;i<row;i++)
        G[i]=(double *)malloc(col * sizeof(double));

    for(int i=0;i<row;i++)
        WH[i]=(double *)malloc(col * sizeof(double));


    multiply(WH,W, H, row, k, col);
    FindG(G,V,WH,row,col);
    free_matrix(WH, row);


    //cout<<"D matrix is:"<<endl;
    //print_matrix(D,row,col);

    transpose(W, transpose_W, row, k); // WT
    multiply(transpose_WV, transpose_W, V, k, row, col);// WT*V
    multiply(transpose_WVG, transpose_WV, G, k, col, row);// WT*VD
    free_matrix(transpose_WV, k);

    double *AdjM[N],*degree[N],*HAdj[N],*alphaHAdj[N],*num1[N],*betaH[N];

    for (int i = 0; i < row; i++)
        AdjM[i] = (double *)malloc(col * sizeof(double));

    for (int i = 0; i < row; i++)
        degree[i] = (double *)malloc(col * sizeof(double));

    for (int i = 0; i < k; i++)
        HAdj[i] = (double *)malloc(col * sizeof(double));

    for (int i = 0; i < k; i++)
        alphaHAdj[i] = (double *)malloc(col * sizeof(double));

    for (int i = 0; i < k; i++)
        num1[i] = (double *)malloc(col * sizeof(double));

    for (int i = 0; i < k; i++)
        betaH[i] = (double *)malloc(col * sizeof(double));

    AdjacentMatrix(AdjM,V,row,col);
    DegreeMatrix(degree,AdjM,row,col);

    multiply(HAdj,H,AdjM,k,col,row);
    free_matrix(AdjM, row);

    constMultiplication(alphaHAdj,HAdj,alpha,k,col);
    free_matrix(HAdj, k);

    add_element_wise(num1,transpose_WVG,alphaHAdj,k,col);
    free_matrix(transpose_WVG, k);
    free_matrix(alphaHAdj, k);


    constMultiplication(betaH,H,beta,k,col);
    add_element_wise(numerator,num1,betaH,k,col);
    free_matrix(num1, k);
    free_matrix(betaH, k);


    double *den_part1[N], *denominator[N],*transpose_WW[N],*transpose_WWH[N],*transpose_WWHG[N],*Hdegree[N],*alphaHdegree[N],*En[N],*HEn[N],*betaHE[N],*E[N];

    for (int i = 0; i < k; i++)
        den_part1[i] = (double *)malloc(col * sizeof(double));

    for (int i = 0; i < k; i++)
        transpose_WW[i] = (double *)malloc(k * sizeof(double));

    for (int i = 0; i < k; i++)
        transpose_WWH[i] = (double *)malloc(col * sizeof(double));

    for (int i = 0; i < k; i++)
        transpose_WWHG[i] = (double *)malloc(col * sizeof(double));

    for (int i = 0; i < k; i++)
        Hdegree[i] = (double *)malloc(col * sizeof(double));

    for (int i = 0; i < k; i++)
        alphaHdegree[i] = (double *)malloc(col * sizeof(double));

    for (int i = 0; i < k; i++)
        denominator[i] = (double *)malloc(col * sizeof(double));

     for (int i = 0; i < row; i++)
        En[i] = (double *)malloc(col * sizeof(double));

     for (int i = 0; i < k; i++)
        HEn[i] = (double *)malloc(col * sizeof(double));

     for (int i = 0; i < k; i++)
        betaHE[i] = (double *)malloc(col * sizeof(double));

    for (int i = 0; i < row; i++)
        E[i] = (double *)malloc(col * sizeof(double));


    multiply(transpose_WW, transpose_W, W, k, row, k); // WT*W
    free_matrix(transpose_W, k);

    multiply(transpose_WWH, transpose_WW, H, k, k, col);//(WT*W)*H
    free_matrix(transpose_WW, k);

    multiply(transpose_WWHG, transpose_WWH, G, k, col, row);//(WT*W)*HD
    free_matrix(G, row);
    free_matrix(transpose_WWH, k);

    multiply(Hdegree,H,degree,k,col,col);
    free_matrix(degree, row);

    constMultiplication(alphaHdegree,Hdegree,alpha,k,col);
    free_matrix(Hdegree, k);

    add_element_wise(den_part1,transpose_WWHG,alphaHdegree,k,col);
    free_matrix(transpose_WWHG, k);
    free_matrix(alphaHdegree, k);

    Ecalculation(E,row,col);
    constMultiplication(En,E,1/row,row,col);
    free_matrix(E, row);

    multiply(HEn,H,En,k,col,col);
    free_matrix(En, row);

    constMultiplication(betaHE,HEn,beta,k,col);
    free_matrix(HEn, k);

    add_element_wise(denominator,den_part1,betaHE,k,col);
    free_matrix(den_part1, k);
    free_matrix(betaHE, k);

    double *updated_H[N];                           // the term that is to be multiplied with H

    for (int i = 0; i < k; i++)
        updated_H[i] = (double *)malloc(col * sizeof(double));
    divide_element_wise(updated_H, numerator, denominator, k, col);

    double *ans_H[N];

    for (int i = 0; i < k; i++)
        ans_H[i] = (double *)malloc(col * sizeof(double));
    multiply_element_wise(ans_H, H, updated_H, k, col);

    copy_matrix(ans_H, H, k, col);
    free_matrix(ans_H, k);
    free_matrix(updated_H, k);
   
    free_matrix(numerator, k);
    
    

   
    
    
    

    
    
    
    

    
    free_matrix(denominator, k);
    
    

    
    
   
   

    
    
    

}


void update_W1(double **W, double **H, double **V, int row, int k, int col)
{
    double *HT[N], *numerator[N],*VG[N],*G[N],*WH[N];

    for (int i = 0; i < col; i++)
        HT[i] = (double *)malloc(k * sizeof(double));

    for (int i = 0; i < row; i++)
        VG[i] = (double *)malloc(col * sizeof(double));

    for (int i = 0; i < row; i++)
        numerator[i] = (double *)malloc(k * sizeof(double));

    for(int i=0;i<row;i++)
        G[i]=(double *)malloc(col * sizeof(double));

    for(int i=0;i<row;i++)
        WH[i]=(double *)malloc(col * sizeof(double));


    multiply(WH,W, H, row, k, col);

    //cout<<"WH matrix is:"<<endl;
    //print_matrix(WH,row,col);
    FindG(G,V,WH,row,col);
    free_matrix(WH, row);
    //cout<<"G matrix is:"<<endl;
   // print_matrix(G,row,col);

    //cout<<"D matrix is:"<<endl;
   // print_matrix(D,row,col);

    transpose(H, HT, k, col);
     //print_matrix(HT,col,k); // HT
    multiply(VG, V, G, row, col, col);
    free_matrix(G, row);
    //cout<<"VG matrix is:"<<endl;
    //print_matrix(VG,row,col);
    multiply(numerator, VG, HT, row, col, k); 
    free_matrix(VG, row);// VG*HT

    double *HG[N],*HG_transpose_H[N],*WHG_transpose_H[N],*WP[N],*gammaWP[N],*P[N];

     for (int i = 0; i < k; i++)
        HG[i] = (double *)malloc( col * sizeof(double));

     for (int i = 0; i < k; i++)
        HG_transpose_H[i] = (double *)malloc( k * sizeof(double));

     for(int i=0;i<row;i++)
         WHG_transpose_H[i] = (double *)malloc( k * sizeof(double));

     for(int i=0;i<k;i++)
        P[i] = (double *)malloc( k * sizeof(double));

      for(int i=0;i<row;i++)
        WP[i] = (double *)malloc( k * sizeof(double));

      for(int i=0;i<row;i++)
        gammaWP[i] = (double *)malloc( k * sizeof(double));


    multiply(HG, H, G, k, col, row); // HD
    multiply(HG_transpose_H, HG, HT, k, col, k);// HD_transpose_H
    free_matrix(HT, col);
    free_matrix(HG, k);

    multiply(WHG_transpose_H,W,HG_transpose_H,row,k,k);
    free_matrix(HG_transpose_H, k);

    FindP(P,W,row,k);
    multiply(WP,W,P,row,k,k);
    free_matrix(P, k);
    constMultiplication(gammaWP,WP,gamma,row,k);
    free_matrix(WP, row);



    double *denominator[N];
    for (int i = 0; i < row; i++)
        denominator[i] = (double *)malloc(k * sizeof(double));

    add_element_wise(denominator,WHG_transpose_H ,gammaWP , row, k);
    free_matrix(gammaWP, row);

    free_matrix(WHG_transpose_H, row);

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
    
    free_matrix(numerator, row);
    

    
    
   
    


    
    
    
    
    free_matrix(denominator, row);
    free_matrix(updated_W, row);

}



void RSNMF()
{

    double *matrix[N];
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


    // Generating two matrices using user input dimension and random number generator

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

    double initial_cost = cost;
    double prev_cost = 0;
    printf("Updating costs:\n ");
    while (cost > EPSILON)
    {

        if ((counter % 2) == 0)
        {
           update_W1(W, H, matrix, row, k, col);
            print_matrix(W,row,k);

        }
        else
        {
          // cout<<row<<k<<col;
            update_H1(W, H, matrix, row, k, col);
            cout << "New H" << endl;
            print_matrix(H,k,col);
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








/*
double *transpose_W[N], *numerator[N],*transpose_WV[N],*transpose_WVG[N],*G[N],*WH[N];
    cout<<"jai na";
    // allocating transpose_W and numerator
    for (int i = 0; i < k; i++)
        transpose_W[i] = (double *)malloc(row * sizeof(double));

    for (int i = 0; i < k; i++)
        numerator[i] = (double *)malloc(col * sizeof(double));

    for(int i=0;i<k;i++)
        transpose_WV[i]=(double *)malloc(col * sizeof(double));

     for(int i=0;i<k;i++)
        transpose_WVG[i]=(double *)malloc(col * sizeof(double));

    for(int i=0;i<row;i++)
        G[i]=(double *)malloc(col * sizeof(double));

    for(int i=0;i<row;i++)
        WH[i]=(double *)malloc(col * sizeof(double));


    multiply(WH,W, H, row, k, col);
    FindG(G,V,WH,row,col);


    //cout<<"D matrix is:"<<endl;
    //print_matrix(D,row,col);

    transpose(W, transpose_W, row, k); // WT
    multiply(transpose_WV, transpose_W, V, k, row, col);// WT*V
    multiply(transpose_WVG, transpose_WV, G, k, col, row);// WT*VD


    double *AdjM[N],*degree[N],*HAdj[N],*alphaHAdj[N],*num1[N],*betaH[N];

    for (int i = 0; i < row; i++)
        AdjM[i] = (double *)malloc(col * sizeof(double));

    for (int i = 0; i < row; i++)
        degree[i] = (double *)malloc(col * sizeof(double));

    for (int i = 0; i < k; i++)
        HAdj[i] = (double *)malloc(col * sizeof(double));

    for (int i = 0; i < k; i++)
        alphaHAdj[i] = (double *)malloc(col * sizeof(double));

    for (int i = 0; i < k; i++)
        num1[i] = (double *)malloc(col * sizeof(double));

    for (int i = 0; i < k; i++)
        betaH[i] = (double *)malloc(col * sizeof(double));

    AdjacentMatrix(AdjM,V,row,col);
    DegreeMatrix(degree,AdjM,row,col);

    multiply(HAdj,H,AdjM,k,col,row);
    constMultiplication(alphaHAdj,HAdj,alpha,k,col);
    add_element_wise(num1,transpose_WVG,alphaHAdj,k,col);
    constMultiplication(betaH,H,beta,k,col);
    add_element_wise(numerator,num1,betaH,k,col);


    double *den_part1[N], *denominator[N],*transpose_WW[N],*transpose_WWH[N],*transpose_WWHG[N],*Hdegree[N],*alphaHdegree[N],*En[N],*HEn[N],*betaHE[N],*E[N];

    for (int i = 0; i < k; i++)
        den_part1[i] = (double *)malloc(col * sizeof(double));

    for (int i = 0; i < k; i++)
        transpose_WW[i] = (double *)malloc(k * sizeof(double));

    for (int i = 0; i < k; i++)
        transpose_WWH[i] = (double *)malloc(col * sizeof(double));

    for (int i = 0; i < k; i++)
        transpose_WWHG[i] = (double *)malloc(col * sizeof(double));

    for (int i = 0; i < k; i++)
        Hdegree[i] = (double *)malloc(col * sizeof(double));

    for (int i = 0; i < k; i++)
        alphaHdegree[i] = (double *)malloc(col * sizeof(double));

    for (int i = 0; i < k; i++)
        denominator[i] = (double *)malloc(col * sizeof(double));

     for (int i = 0; i < row; i++)
        En[i] = (double *)malloc(col * sizeof(double));

     for (int i = 0; i < k; i++)
        HEn[i] = (double *)malloc(col * sizeof(double));

     for (int i = 0; i < k; i++)
        betaHE[i] = (double *)malloc(col * sizeof(double));

    for (int i = 0; i < row; i++)
        E[i] = (double *)malloc(col * sizeof(double));


    multiply(transpose_WW, transpose_W, W, k, row, k); // WT*W
    multiply(transpose_WWH, transpose_WW, H, k, k, col);//(WT*W)*H
    multiply(transpose_WWHG, transpose_WWH, G, k, col, row);//(WT*W)*HD
    multiply(Hdegree,H,degree,k,col,col);
    constMultiplication(alphaHdegree,Hdegree,alpha,k,col);
    add_element_wise(den_part1,transpose_WWHG,alphaHdegree,k,col);
    Ecalculation(E,row,col);
    constMultiplication(En,E,1/row,row,col);
    multiply(HEn,H,En,k,col,col);
    constMultiplication(betaHE,HEn,beta,k,col);
    add_element_wise(denominator,den_part1,betaHE,k,col);


    double *updated_H[N];                           // the term that is to be multiplied with H

    for (int i = 0; i < k; i++)
        updated_H[i] = (double *)malloc(col * sizeof(double));
    divide_element_wise(updated_H, numerator, denominator, k, col);

    double *ans_H[N];

    for (int i = 0; i < k; i++)
        ans_H[i] = (double *)malloc(col * sizeof(double));
    multiply_element_wise(ans_H, H, updated_H, k, col);

    copy_matrix(ans_H, H, k, col);
    free_matrix(ans_H, k);
    free_matrix(updated_H, k);
    free_matrix(transpose_W, k);
    free_matrix(numerator, k);
    free_matrix(transpose_WV, k);
    free_matrix(transpose_WVG, k);

    free_matrix(G, row);
    free_matrix(WH, row);
    free_matrix(AdjM, row);
    free_matrix(degree, row);

    free_matrix(HAdj, k);
    free_matrix(alphaHAdj, k);
    free_matrix(num1, k);
    free_matrix(betaH, k);

    free_matrix(den_part1, k);
    free_matrix(denominator, k);
    free_matrix(transpose_WW, k);
    free_matrix(transpose_WWH, k);

    free_matrix(transpose_WWHG, k);
    free_matrix(Hdegree, k);
    free_matrix(alphaHdegree, k);
    free_matrix(En, row);

    free_matrix(HEn, k);
    free_matrix(betaHE, k);
    free_matrix(E, row);

    */
