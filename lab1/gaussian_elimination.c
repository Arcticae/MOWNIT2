
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define num double

#define true 1
#define false 0
const num trunk = 6;
num eps;

void calc_eps() {
    eps = 1 / pow(10, trunk);
}

num getTrunked(num what, num places) {
    return (floor((what * pow(10, places) + 0.5)) / pow(10, places));
}

int gaussian_elimination_1(int n, num **AB, num *X) {
    int i, j, k;
    num m, s;


    for (i = 0; i < n - 1; i++) {
        for (j = i + 1; j < n; j++) {
            if (fabs(AB[i][i]) < eps) return false;
            m = -AB[j][i] / AB[i][i];
            for (k = i + 1; k <= n; k++)
                AB[j][k] += m * AB[i][k];
        }
    }


    for (i = n - 1; i >= 0; i--) {
        s = AB[i][n];
        for (j = n - 1; j >= i + 1; j--)
            s -= AB[i][j] * X[j];
        if (fabs(AB[i][i]) < eps) return false;
        X[i] = s / AB[i][i];
    }
    return true;
}

void show_euclides_distance(num *result, num *exact_result, int size) {
    printf("Euclides distance from the right result is : \n");
    for (int i = 0; i < size; i++) {
        printf("%.16lf ", fabs(result[i] - exact_result[i]));
    }
    printf("\n");
}

void show_max_distance(num *result, num *exact_result, int size) {
    printf("Max metric distance from the right result is : \n");
    num max = 0.0, tmp;
    for (int i = 0; i < size; i++) {
        if ((tmp = fabs(result[i] - exact_result[i])) > max)
            max = tmp;
    }
    printf("\t %.16lf \t", max);
    printf("\n");
}

int main(int argc, char **argv) {
    if (argc < 2) {
        fprintf(stderr, "Not enough arguments supplied to program, give me the size of matrix.\n");
        exit(1);
    }
    int size = atoi(argv[1]);
    calc_eps();
    int i, j;
    num **A = (num **) malloc(size * sizeof(num *));
    num *B = malloc(size * sizeof(num));
    for (i = 0; i < size; i++)
        A[i] = malloc(size * sizeof(num));
    printf("Give A matrix parameters:\n");
    for (i = 0; i < size; i++)
        for (j = 0; j < size; j++) {
            scanf("%lf", &A[i][j]);
        }
    printf("Give B matrix parameters: \n");
    for (i = 0; i < size; i++)
        scanf("%lf", &B[i]);

    num **AB = (num **) malloc(size * sizeof(num *));
    for (i = 0; i < size; i++)
        AB[i] = (num *) malloc((size + 1) * sizeof(num));

    for (i = 0; i < size; i++) {
        for (j = 0; j < size + 1; j++) {
            if (j != size)
                AB[i][j] = A[i][j];
            else
                AB[i][j] = B[i];
        }
    }
    num *result = (num *) malloc(size * sizeof(num));

    for (i = 0; i < size; i++) {
        for (j = 0; j < size + 1; j++) {
            AB[i][j] = getTrunked(AB[i][j], trunk);
            printf("%.lf ", AB[i][j]);
        }
        printf("\n");
    }
    if (gaussian_elimination_1(size, AB, result) == false) {
        printf("Unsolvable matrixes given\n");
        exit(1);
    }
    printf("Results:\n");
    for (i = 0; i < size; i++)
    {
        result[i] = getTrunked(result[i],trunk);
        printf("%.16lf  ", result[i]);
    }

    printf("Do you want to compare the results? \n Y/N");
    char choice;
    scanf(" %c", &choice);
    if (choice == 'Y' || choice == 'y') {
        num right_result[size];
        printf("Give me the right result: \n");
        for (int j = 0; j < size; j++) {
            scanf("%lf", &right_result[j]);
        }
        show_euclides_distance(result, right_result, size);
        show_max_distance(result, right_result, size);
    }

    return 0;

}