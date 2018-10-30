#define _DEFAULT_SOURCE

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <sys/time.h>
#include <unistd.h>
#include <sys/resource.h>
#define num double


#define true 1
#define false 0
num eps;
num precision = 17.0;

void calc_eps() {
    eps = 1 / (num) pow(10, precision);
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
    num sum_of_squares = 0.0;
    for (int i = 0; i < size; i++) {
        sum_of_squares += pow(fabs(result[i] - exact_result[i]), 2.0);
    }
    printf("Euclides distance from the right result is : \n\t%.16lf\n", sqrt(sum_of_squares));
    printf("\n");
}

void show_relative_error(num *result, num *exact_result, int size) {
    num mean = 0.0;

    for (int i = 0; i < size; i++) {
        mean += fabs(result[i] - exact_result[i]);
    }
    mean /= (num) size;
    printf("Mean relative error of the result is : \n\t%.16lf\n", mean);
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

num *zad1_calculate_b(num **A, num *X, int n) {
    num *B = (num *) malloc(n * sizeof(num *));

    for (int i = 0; i < n; i++)
        B[i] = 0.0;

    for (int i = 0; i < n; i++)
        for (int k = 0; k < n; k++)
            B[i] += A[i][k] * X[k];


    return B;
}

num *zad3_solve_thomas(num **AB, int size) {
    num *l = malloc(size * sizeof(num));
    num *u = malloc(size * sizeof(num));
    num *a = malloc(size * sizeof(num));
    num *c = malloc(size * sizeof(num));

    num *x = malloc(size * sizeof(num));
    num *y = malloc(size * sizeof(num));

    u[0] = AB[0][0];
    c[0] = AB[0][1];
    for (int i = 1; i < size - 1; i++) {
        a[i] = AB[i][i - 1];
        c[i] = AB[i][i + 1];

        l[i] = (a[i] / u[i - 1]);
        u[i] = AB[i][i] - (l[i] * c[i - 1]);
    }
    a[size - 1] = AB[size - 1][size - 2];
    l[size - 1] = (a[size - 1] / u[size - 2]);
    u[size - 1] = AB[size - 1][size - 1] - (l[size - 1] * c[size - 2]);


    // L * y = b
    // U * x = y

    y[0] = AB[0][size];
    for (int i = 1; i < size; i++) {
        y[i] = (AB[i][size] - (l[i] * y[i - 1]));
    }


    x[size - 1] = (y[size - 1] / u[size - 1]);
    for (int i = size - 2; i >= 0; i--) {
        x[i] = ((y[i] - (c[i] * x[i + 1])) / u[i]);
    }


    return x;
}

num **zad1_generate_a(int n) {
    num **A = (num **) malloc(n * sizeof(num *));
    for (int i = 0; i < n; i++)
        A[i] = malloc(n * sizeof(num));

    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++) {
            A[i][j] = (i == 0 ? 1.0 : 1.0 / ((num) i + 1.0 + (num) j + 1.0 - 1.0));
        }
    return A;
}

num **zad2_generate_a(int n) {
    num **A = (num **) malloc(n * sizeof(num *));
    for (int i = 0; i < n; i++)
        A[i] = malloc(n * sizeof(num));

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if (j >= i) {
                A[i][j] = (2.0 * ((num) i + 1.0)) / ((num) j + 1.0);
            } else {
                A[i][j] = A[j][i];
            }
        }
    }
    return A;
}

num **zad3_generate(int n, num m, num k) {
    num **A = (num **) malloc(n * sizeof(num *));
    for (int i = 0; i < n; i++)
        A[i] = malloc(n * sizeof(num));

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if (i == j)
                A[i][j] = -m * ((num) i + 1.0) - k;
            if (i == j - 1)
                A[i][j] = (num) i + 1.0;
            if (i == j + 1 && i > 0)
                A[i][j] = m / ((num) i + 1.0);
            if (j < i - 1 || j > i + 1)
                A[i][j] = 0.0;
        }
    }
    return A;

}

num *zad1_get_x(int n) {
    num *X = malloc(n * sizeof(num));
    num set[2] = {1.0, -1.0};
    for (int i = 0; i < n; i++) {
        X[i] = set[rand() % 2];
    }
    return X;
}


num **concat_AB_arrays(num **A, num *B, int size) {
    int i, j;
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
    return AB;
}



int main(int argc, char **argv) {
    srand((unsigned int) time(NULL));

    calc_eps();
    //   for (int size = 100; size < 10000; size += 100) {
    int i, j, precision = 17.0;
    int size;
    printf("What size of the matrix?\n");
    scanf("%d", &size);
    num **A, **AB;
    num *B;
    num *X, *Xprim;
    Xprim = malloc(size * sizeof(num));

    /*
    //Uncomment segment for zad 1 contents
    printf("Zad1 A:\n");
    A = zad1_generate_a(size);
    //truncate_2d(A, precision, size);
    print_2d_matrix(A, size);

    printf("Zad1 X:\n");
    X = zad1_get_x(size);
    //truncate_1d(X, precision, size);
    print_1d_matrix(X, size);

    printf("Zad1 B:\n");
    B = zad1_calculate_b(A, X, size);
    print_1d_matrix(B, size);
    */
    /*
    //Uncomment for zadanie 2
    printf("-------MATRIX SIZE %d ------------\n" ,size);
    //printf("Zad2 A:\n");
    A = zad2_generate_a(size);
    truncate_2d(A, precision, size);
    //print_2d_matrix(A, size);

    //printf("Zad2 X:\n");
    X = zad1_get_x(size);
    truncate_1d(X, precision, size);
    // print_1d_matrix(X,size);

    //printf("Zad2 B:\n");
    B = zad1_calculate_b(A, X, size);
    // print_1d_matrix(B, size);
    */
    num m = 4.0;
    num k = 5.0;

    struct timeval *realtime = malloc(4 * sizeof(struct timeval));
    //realtime[0-1] = thomas's
    //realtime[2-3] = gaussian algorithm

    //Uncomment for zadanie 3
    A = zad3_generate(size, m, k);
    //print_2d_matrix(A, size);

    X = zad1_get_x(size);
    //  print_1d_matrix(X, size);

    B = zad1_calculate_b(A, X, size);
    // print_1d_matrix(B, size);


    AB = concat_AB_arrays(A, B, size);
    num **ABprim = concat_AB_arrays(A, B, size);
    num *Xbis;


    gettimeofday(&realtime[0], NULL);
    Xbis = zad3_solve_thomas(ABprim, size);
    gettimeofday(&realtime[1], NULL);
    gettimeofday(&realtime[2], NULL);
    if (gaussian_elimination_1(size, AB, Xprim) == false) {
        printf("Unsolvable matrixes given\n");
        exit(1);
    }
    gettimeofday(&realtime[3], NULL);

    printf("Gaussian differences:\n\n");
    show_euclides_distance(X, Xprim, size);
    show_max_distance(X, Xprim, size);
    show_relative_error(X, Xprim, size);

    printf("\nThomas algorithm differences: \n\n");
    show_euclides_distance(X, Xbis, size);
    show_max_distance(X, Xbis, size);
    show_relative_error(X, Xbis, size);


    struct timeval *thomas_time = malloc(sizeof(struct timeval));
    timersub(&realtime[1], &realtime[0], thomas_time);
    struct timeval *gaussian_time = malloc(sizeof(struct timeval));
    timersub(&realtime[3], &realtime[2], gaussian_time);

    printf("Time of thomas's algorithm: %ld.%06ld\n", thomas_time->tv_sec, thomas_time->tv_usec);
    printf("Time of gaussian algorithm: %ld.%06ld\n", gaussian_time->tv_sec, gaussian_time->tv_usec);


    return 0;

}