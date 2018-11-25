#include <iostream>
#include <cmath>
#include <fstream>
#include <cstring>

#define num double
using namespace std;
typedef struct point {
    num x;
    num y;
};

num lagrange_interpolate(point *points, int n, num x) {
    num result = 0.0;
    num partial;
    num top, bot;
    for (int i = 0; i < n; i++) {
        top = 1.0, bot = 1.0;
        for (int j = 0; j < n; j++) {
            if (i != j) {
                top *= (x - points[j].x);
                bot *= (points[i].x - points[j].x);
                partial = top / bot;
            }
        }
        result += points[i].y * partial;
    }
    return result;
}
//our function

num inline f_x(num x) {
    return exp(-3 * sin(2 * x));
}

num inline df_x(num x) {
    return -6 * f_x(x) * (num) cos(2 * x);
}

num hermite_interpolate(point *points, int n, num x) {
    point z[n * 2];
    num matrix[2 * n][2 * n];
    memset(matrix, 0, sizeof(matrix));
    for (int i = 0; i < n; i++) {
        z[i] = points[i];
        z[2 * i] = points[i];
    }
    for (int i = 0; i < 2 * n; i++) {
        for (int j = 0; j < n; j++) {
            if (j == 0)
                matrix[i][j] = (z[i].y);
            else if ((j == 1) && i % 2 == 1)
                matrix[i][j] = df_x(z[i].y);
            else {
                matrix[i][j] = matrix[i][j - 1] - matrix[i - 1][j - 1];
                matrix[i][j] = matrix[i][j] / (z[i].x - z[i - j].x);
            }
        }
    }
    num result = 0.0;
    num tmp = 1.0;
    for (int i = 0; i < n; i++) {
        result += matrix[i][i] * tmp;
        tmp *= x - z[i].x;
    }
    return result;
}

num newton_interpolate(point *points, int n, num x) {
    num b[n];
    for (int j = 0; j < n; j++) {
        b[j] = points[j].y;
    }
    for (int j = 1; j < n; j++) {
        for (int k = n - 1; k >= j; k--) {
            b[k] = (b[k] - b[k - 1]) / (points[k].x - points[k - j].x);
        }
    }

    num v[n];
    for (int i = 0; i < n; i++) {
        v[i] = b[i];
    }

    for (int i = n - 2; i >= 0; i--) {
        v[i] = v[i + 1] * (x - points[i].x) + b[i];
    }
    return v[0];
}


int main() {
    //constants
    int n = 10;
    cout << "Give the number of interpolation nodes" << endl;
    cin >> n;
    num a, b;
    //Given
    a = 0.0;
    b = 3 * M_PI;
    num step = b / (num) n;
    point even_points[n];
    point chebyshev_points[n];
    cout << "Given point coordinates" << endl;
    //Points For constant distribution
    for (num i = a, j = 0; i <= b, j < n - 1; i += step, j += 1) {
        even_points[(int) j].x = (num) i;
        cout << "X: " << even_points[(int) j].x;
        even_points[(int) j].y = f_x(i);
        cout << " Y: " << even_points[(int) j].y << endl;
    }

    even_points[n - 1].y = (num) f_x(b);
    even_points[n - 1].x = (num) b;
    //Points for chebyshev nodes
    for (int i = 1; i <= n; i++) {
        chebyshev_points[i - 1].x = (1.0 / 2.0 * (a + b) +
                                     1.0 / 2.0 * (b - a) * (num) cos((2.0 * (num) i - 1.0) * M_PI / (2.0 * (num) n)));
        chebyshev_points[i - 1].y = f_x(chebyshev_points[i - 1].x);
    }

    ofstream lagrange_results, real_results, newton_results, interpolation_points, hermite_results;
    interpolation_points.open("/home/timelock/Pulpit/interpolation_points.csv");
    lagrange_results.open("/home/timelock/Pulpit/lagrange_resutlts.csv");
    newton_results.open("/home/timelock/Pulpit/newton_results.csv");
    real_results.open("/home/timelock/Pulpit/real_results.csv");
    hermite_results.open("/home/timelock/Pulpit/hermite_results.csv");

    //UNCOMMENT FOR PART WITH EVEN POINTS
    /*
    interpolation_points << "X,Y" << endl;
    for (int i = 0; i < n; i++) {
        interpolation_points << even_points[i].x << "," << even_points[i].y << endl;
    }
    lagrange_results << "X,Y" << endl;
    for (num i = a; i <= b; i += 0.01) {
        lagrange_results << i << "," << lagrange_interpolate(even_points, n, i) << endl;
    }

    newton_results << "X,Y" << endl;
    for (num i = a; i <= b; i += 0.01) {
        newton_results << i << "," << newton_interpolate(even_points, n, i) << endl;
    }
    real_results << "X,Y" << endl;
    for (num i = a; i <= b; i += 0.01) {
        real_results << i << "," << f_x(i) << endl;
    }
    hermite_results << "X,Y" << endl;

    for (num i = a; i <= b; i += 0.01) {
        hermite_results << i << "," << hermite_interpolate(even_points, n, i) << endl;
    }
     */
    //END UNCOMMENT FOR PART WITH EVEN POINTS

    //UNCOMMENT FOR PART WITH CHEBYSHEV POINTS

    interpolation_points << "X,Y" << endl;
    for (int i = 0; i < n; i++) {
        interpolation_points << chebyshev_points[i].x << "," << chebyshev_points[i].y << endl;
    }
    lagrange_results << "X,Y" << endl;
    for (num i = a; i <= b; i += 0.01) {
        lagrange_results << i << "," << lagrange_interpolate(chebyshev_points, n, i) << endl;
    }

    newton_results << "X,Y" << endl;
    for (num i = a; i <= b; i += 0.01) {
        newton_results << i << "," << newton_interpolate(chebyshev_points, n, i) << endl;
    }
    real_results << "X,Y" << endl;
    for (num i = a; i <= b; i += 0.01) {
        real_results << i << "," << f_x(i) << endl;
    }
    hermite_results << "X,Y" << endl;

    for (num i = a; i <= b; i += 0.01) {
        hermite_results << i << "," << hermite_interpolate(chebyshev_points, n, i) << endl;
    }

    //END UNCOMMENT FOR PART WITH CHEBYSHEV POINTS


    hermite_results.close();
    interpolation_points.close();
    lagrange_results.close();
    real_results.close();
    return 0;
}