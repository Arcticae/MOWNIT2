#include <iostream>
#include <cmath>
#include <armadillo>

using namespace std;
using namespace arma;

#define num double

void show_euclides_distance(num *result, num *exact_result, int size) {
    num sum_of_squares = 0.0;
    for (int i = 0; i < size; i++) {
        sum_of_squares += pow(fabs(result[i] - exact_result[i]), 2.0);
    }
    printf("Euclides distance from the right result is : \n\t%.16lf\n", sqrt(sum_of_squares));
    printf("\n");
}

double spec_radius(Mat<double> A, int size) {
    vec N(size);
    for (int i = 0; i < size; i++) {
        N(i) = 1 / A(i, i);
    }
    Mat<double> M(size, size);
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            if (j == i) M(i, j) = 0.0;
            else M(i, j) = -(A(i, j) * N(i));
        }
    }
    cx_vec eigens = eig_gen(M);
    double result = max(abs(eigens));

    return result;
}

num mean_relative_error(num *result, num *exact_result, int size) {
    num mean = 0.0;

    for (int i = 0; i < size; i++) {
        mean += fabs(result[i] - exact_result[i]);
    }
    mean = mean / (num) size;

    return mean;
}

num max_distance(num *result, num *exact_result, int size) {
    num max = 0.0, tmp;
    for (int i = 0; i < size; i++) {
        if ((tmp = fabs(result[i] - exact_result[i])) > max)
            max = tmp;
    }
    return max;
}

num *add_1d_to_1d(num *A, num *B, int n) {
    num *Aprim = new num[n];
    for (int i = 0; i < n; i++) {
        Aprim[i] = A[i] + B[i];
    }
    return Aprim;
}

num **mul_2d_with_2d(num **A, num **B, int n) {
    num **RES = new num *[n];
    for (int i = 0; i < n; ++i) {
        RES[i] = new num[n];
        for (int j = 0; j < n; ++j) {
            RES[i][j] = 0.0;
            for (int k = 0; k < n; ++k) {
                RES[i][j] += A[i][k] * B[k][j];
            }
        }
    }
    return RES;
}

num *mul_2d_with_1d(num **A, num *X, int n) {
    num *B = (num *) malloc(n * sizeof(num));

    for (int i = 0; i < n; i++)
        B[i] = 0.0;

    for (int i = 0; i < n; i++)
        for (int k = 0; k < n; k++)
            B[i] += A[i][k] * X[k];


    return B;
}

num *zad1_get_x(int n) {
    num *X = new num[n];
    num set[2] = {1.0, -1.0};
    for (int i = 0; i < n; i++) {
        X[i] = set[rand() % 2];
    }
    return X;
}

void print_1d_matrix(num *B, int n) {
    cout << endl;
    for (int i = 0; i < n; i++)cout << B[i] << " ";
    cout << endl;
}

void print_2d_matrix(num **A, int n) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++)
            cout << A[i][j] << " ";

        cout << endl;
    }
}

num **zad1_get_A(int n) {
    num **A = new num *[n];
    for (int i = 0; i < n; i++) {
        A[i] = new num[n];
    }

    for (int i = 0; i < n; i++) {
        A[i][i] = 6.0;

        for (int j = 0; j < n; j++) {
            if (j > i)
                A[i][j] = pow(-1, (num) j) * (0.5 / (num) j);
            if (i == j + 1)
                A[i][j] = (0.5 / (num) (i));
            if (j < i - 1)
                A[i][j] = 0.0;
        }
    }


    return A;
}
Mat<num> convert_to_mat(num**A,int n){
    Mat<num> res(n,n);
    for(int i = 0;i<n;i++){
        for(int j =0 ;j<n;j++){
            res(i,j)= A[i][j];
        }
    }
    return res;
}

void add_2d_with_2d(num **A, num **B, int n) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            A[i][j] += B[i][j];
        }
    }
}

num **calculate_N(num **D, int n) {
    num **N = new num *[n];

    for (int i = 0; i < n; i++) {
        N[i] = new num[n];
        for (int j = 0; j < n; j++) {
            N[i][j] = 0.0;
            if (i == j)
                N[i][j] = 1 / D[i][j];
        }
    }
    return N;
}

num **calculate_M(num **D, num **L, num **U, int n) {
    num **M = calculate_N(D, n);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            M[i][j] = -M[i][j];
        }
    }

    return M;
}


num *jacobi_iteration(num **A, num *B, num *X, int n) {
    num **L = new num *[n];
    num **D = new num *[n];
    num **U = new num *[n];
    for (int i = 0; i < n; i++) {
        D[i] = new num[n];
        L[i] = new num[n];
        U[i] = new num[n];

        for (int j = 0; j < n; j++) {
            D[i][j] = 0.0;
            L[i][j] = 0.0;
            U[i][j] = 0.0;
            if (i == j)
                D[i][j] = A[i][j];
            if (i < j)
                L[i][j] = A[i][j];
            if (i > j)
                U[i][j] = A[i][j];
        }
    }
    num **M = calculate_M(D, L, U, n);
    num **N = calculate_N(D, n);
    num *ret = add_1d_to_1d(mul_2d_with_1d(M, X, n), mul_2d_with_1d(N, B, n), n);
    for (int i = 0; i < n; i++) {
        delete M[i], N[i], L[i], D[i], U[i];
    }
    delete M, N, L, D, U;
    return ret;
}


int main() {
    srand(static_cast<unsigned int>(time(nullptr)));
    int siz, step;
    num eps;
    cout << "Biggest Size of the matrix?";//has to be changed in a quest every time for measurements
    cin >> siz;
    cout << "Step in size?";
    cin >> step;
    cout << "Epsillion?";
    cin >> eps;


    for (int n = step; n <= siz; n += step) {
        cout << "==============SIZE " << n << "===============" << endl;
        num **A = zad1_get_A(n);
        num *Xreal = zad1_get_x(n);
        num *B = mul_2d_with_1d(A, Xreal, n);
        num *X = new num[n];
        num *Xbis = new num[n];
        for (int i = 0; i < n; i++) {
            X[i] = 0.0;
            Xbis[i] = 0.0;

        }
        cout << endl;

        //zad1 begin
        int iter = 0;
        num rel_error;
        num max_dist;
        num rel_diff;
        do {

            if (iter % 2 == 0) {
                Xbis = jacobi_iteration(A, B, X, n);
                rel_error = mean_relative_error(mul_2d_with_1d(A, Xbis, n), B, n);
                rel_diff = mean_relative_error(X, Xbis, n);
                max_dist = max_distance(mul_2d_with_1d(A, Xbis, n), B, n);
            } else {
                X = jacobi_iteration(A, B, Xbis, n);
                rel_error = mean_relative_error(mul_2d_with_1d(A, X, n), B, n);
                rel_diff = mean_relative_error(X, Xbis, n);
                max_dist = max_distance(mul_2d_with_1d(A, X, n), B, n);
            }
            iter++;
            cout << "Relative error: " << rel_error << " Relative difference: " << rel_diff << " Max distance: "
                 << max_dist << endl;

        } while ((rel_error > eps || iter < 2) && iter < 100);
        cout<<endl<<"Spectral: "<< spec_radius(convert_to_mat(A,n),n);
        cout << "There were " << iter << " iterations" << endl << "------RESULT-------" << endl;
        if ((iter - 1) % 2 == 0) {
            print_1d_matrix(Xbis, n);

        } else {
            print_1d_matrix(X, n);
        }
        cout << endl << "MEAN RELATIVE OF AX-b IS:" << rel_error;
    }
    return 0;
}