#include <iostream>
#include "matrix_factorizations.h"

int main()
{
    printf("Welcome to matrix factorization!\n");
    printf("1. Matrix factorization using l2,1 norm");

    int choice;
    printf("\nEnter your choice: ");
    scanf("%d", &choice);
    if (choice >= 1 && choice <= 5)
    {
        switch (choice)
        {
        case 1:
            printf("You have chosen l21 norm matrix factorization\n");
            l21_norm();
            break;

        }
    }
}

/*
1
5 5
1 3 4 5 6
1 2 4 5 1
1 2 3 5 6
6 7 8 9 10
11 12 13 14 15
3
2 3 4
3 4 5
6 7 8
9 10 11
2 4 5
3 4 5 6 7
1 2 3 4 5
2 7 8 9 10*/

