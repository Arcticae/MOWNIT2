//
// Created by timelock on 29.10.18.
//

#include <iostream>

#include <stdlib.h>
#include <armadillo>
#define num float
num m=4.0;
num k=5.0;
using namespace arma;
mat zad1_generate_a(int n) {
    mat* A = new mat(n,n);
    mat X = *A;
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++) {
            X(i,j) = (i == 0 ? 1.0 : 1.0 / ((num) i + 1.0 + (num) j + 1.0 - 1.0));
        }
    return X;
}

mat zad2_generate_a(int n) {
    mat *A = new mat(n,n);
    mat X= *A;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if (j >= i) {
                X(i,j)= (2.0 * ((num) i + 1.0)) / ((num) j + 1.0);
            } else {
                X(i,j)= X(j,i);
            }
        }
    }
    return X;
}

mat zad3_generate(int n) {
    mat* A = new mat(n,n);
    mat X = *A;

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if (i == j)
                X(i,j)= -m * ((num) i + 1.0) - k;
            if (i == j - 1)
                X(i,j) = (num) i + 1.0;
            if (i == j + 1 && i > 0)
                X(i,j) = m / ((num) i + 1.0);
            if (j < i - 1 || j > i + 1)
                X(i,j) = 0.0;
        }
    }
    return X;

}

int main(){
    int a;
    std::cin >> a;


    mat matrix = zad1_generate_a(a);
    mat matrix2 = zad2_generate_a(a);
    mat matrix3 = zad3_generate(a);

    cout << "Zad1 cond" << cond(matrix)<<endl;
    cout << "Zad2 cond" << cond(matrix2)<<endl;
    cout << "Zad3 cond" << cond(matrix3)<<endl;
}