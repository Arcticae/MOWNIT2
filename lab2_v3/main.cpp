#include <iostream>
#include <cmath>
#include <vector>
#include <ctime>
#include <armadillo>

#define num double
using namespace std;
using namespace arma;

void show_euclides_distance(vec result, vec exact_result, int size) {
    num sum_of_squares = 0.0;
    for (int i = 0; i < size; i++) {
        sum_of_squares += pow(fabs(result(i) - exact_result(i)), 2.0);
    }
    printf("Euclides distance from the right result is : \n\t%.16lf\n", sqrt(sum_of_squares));
    printf("\n");
}

num mean_relative_error(vec result, vec exact_result, int size) {
    num mean = 0.0;

    for (int i = 0; i < size; i++) {
        mean += fabs(result(i) - exact_result(i));
    }
    mean = mean / (num) size;

    return mean;
}



vec jacobi2(Mat<num> A, vec B, vec start, num ro, int *iter, int size) {


    vec N(size);
    for (int i = 0; i < size; i++) {
        N(i) = 1 / A(i, i);
    }
    Mat<num> M(size, size);
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            if (j == i) M(i, j) = 0.0;
            else M(i, j) = -(A(i, j) * N(i));
        }
    }
    vec X = zeros(size);
    vec tmp(size);
    int it = 1;
    while (1) {
        for (int i = 0; i < size; i++) {
            tmp(i) = N(i) * B(i);
            for (int j = 0; j < size; j++) {
                tmp(i) += M(i, j) * X(j);
            }
        }
        if (max(abs(A * tmp - B)) < ro) {
            for (int i = 0; i < size; i++) {
                X(i) = tmp(i);
            }
            (*iter) = it;
            break;
        }
        for (int i = 0; i < size; i++) {
            X(i) = tmp(i);
        }
        it++;
    }
    return X;
}

vec jacobi1(Mat<num> A, vec B, vec start, num ro, int *iter, int n) {


    vec N(n);
    for (int i = 0; i < n; i++) {
        N(i) = 1 / A(i, i);
    }
    Mat<num> M(n, n);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if (j == i) M(i, j) = 0.0;
            else M(i, j) = -(A(i, j) * N(i));
        }
    }
    vec X = zeros(n);
    vec tmp(n);
    int it = 1;
    while (1) {
        for (int i = 0; i < n; i++) {
            tmp(i) = N(i) * B(i);
            for (int j = 0; j < n; j++) {
                tmp(i) += M(i, j) * X(j);
            }
        }
        if (max(abs(tmp - X)) < ro) {
            for (int i = 0; i < n; i++) {
                X(i) = tmp(i);
            }
            (*iter) = it;
            break;
        }
        for (int i = 0; i < n; i++) {
            X(i) = tmp(i);
        }
        it++;
    }
    return X;
}

num spec_radius(Mat<num> A, int size) {
    vec N(size);
    for (int i = 0; i < size; i++) {
        N(i) = 1 / A(i, i);
    }
    Mat<num> M(size, size);
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            if (j == i) M(i, j) = 0.0;
            else M(i, j) = -(A(i, j) * N(i));
        }
    }
    cx_vec eigens = eig_gen(M);
    num result = max(abs(eigens));

    return result;
}

Mat<num> Matrix_A(int n) {
    Mat<num> A(n, n);
    for (int i = 0; i < n; i++) {
        A(i, i) = 6.0;
        for (int j = 0; j < n; j++) {
            if (j > i)
                A(i, j) = pow(-1, (num) j) * (0.5 / (num) j);
            if (i == j + 1)
                A(i, j) = (0.5 / (num) (i));
            if (j < i - 1)
                A(i, j) = 0.0;
        }

    }
    return A;
}

vec Matrix_X(int n) {
    vec X(n);
    num set[2] = {-1.0, 1.0};
    for (int i = 0; i < n; i++) {
        X(i) = set[rand() % 2];
    }

    return X;
}


int main(int argc, char *argv[]) {
    srand(time(nullptr));

    int maxsize;
    int steplength;
    int beggining;
    num eps;

    cout << "Max Size of matrix?" << endl;
    cin >> maxsize;
    cout<<"Step length?"<<endl;
    cin>> steplength;
    cout<<"Beggining size?"<<endl;
    cin>>beggining;

    eps = 0.0000000000001; //smallest possible*1000, then step by 100times this
    for(int i = 0 ;i<4;i++){
        cout<<"==========================EPS"<<eps<<"================"<<endl;
    for(int size = beggining ; size <= maxsize ; size+=steplength){
        cout<<"==========================SIZE "<<size<<"================"<<endl;
    vec X = Matrix_X(size);
    Mat<num> A = Matrix_A(size);
    vec B(size);
    B = A * X;
    vec X_calc_1;
    vec X_calc_2;
    int iterations1;
    int iterations2;
    cout.precision(10);
    X_calc_1 = jacobi1(A, B, zeros<vec>(size), eps, &iterations1, size);
    X_calc_2 = jacobi2(A, B, zeros<vec>(size), eps, &iterations2, size);
    cout << "First type of ending condition: " << endl;
    cout << endl << "Max error is : " << max(abs(X_calc_1 - X)) << endl;
    cout<<"Mean relative error is :"<< mean_relative_error(X_calc_1,X,size)<<endl;
    show_euclides_distance(X_calc_2,X,size);
    cout << "Number of iterations was " << iterations1 << endl;

    cout << endl << "Second type of ending condtition: " << endl;
    cout << endl << "The error is: " << max(abs(X_calc_2 - X)) << endl;
    cout<<"Mean relative error is :"<< mean_relative_error(X_calc_2,X,size)<<endl;
    show_euclides_distance(X_calc_2,X,size);
    cout << "Number of iterations was " << iterations2 << endl;

    cout << "Spectral Radius of iteration matrix is: " << spec_radius(A, size) << endl;
    cout << "---------------------------------------------------------------" << endl;
    }
    eps*=100;
    }
    return 0;
}