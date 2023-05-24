
#include <iostream>
#include "l21_norm.h"
#include "G_NMF.h"
#include "RSNMF.h"
#include "LSNMF.h"

#include <windows.h>

int main()
{

    printf("Welcome to matrix factorization!\n");
    printf("1. matrix factorization using Lee and Sungs Multiplicative update\n");
    printf("2. matrix factorization using l21 norm\n");
    printf("3. matrix factorization using Graph Laplacian Matrix\n");
    printf("4. matrix factorization using Robust structured Non negative matrix factorization\n");

    int choice;
    printf("\nEnter your choice: ");
    scanf("%d", &choice);
    if (choice >= 1 && choice <= 5)
    {
        switch (choice)
        {

        case 1:
            printf("You have chosen lee and Seiungs non-negative matrix factorization\n");
            LSNMF();
            break;
        case 2:
            printf("You have chosen l21 norm matrix factorization\n");
            l21_norm();
            break;
        case 3:
            printf("You have chosen Graph regularized Non negative matrix factorization\n");
            G_NMF();
            break;

        case 4:
            printf("You have chosen Robust structured Non negative matrix factorization\n");
            RSNMF();
            break;

        }
    }
}

