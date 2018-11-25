#include <iostream>
#include <cmath>
#include <fstream>

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
    return exp(-3*sin(2*x));
}

num newton_interpolate(point * points, int n, num x){
    num b[n];
    for(int j  = 0 ; j< n ; j++){
        b[j] = points[j].y;
    }
    for(int j = 1 ; j<n ; j++ ){
        for(int k = n-1; k>=j ; k--){
            b[k] = (b[k] - b[k-1]) / (points[k].x - points[k-j].x);
        }
    }

    num v[n];
    for(int i = 0; i < n ; i++ ){
        v[i]= b[i];
    }

    for(int i = n-2 ;i >= 0 ;i-- ){
        v[i] = v[i+1]*(x-points[i].x) + b[i];
    }
    return v[0];
}


int main() {
    //constants
    int n = 10;
    cout<<"Give the number of interpolation nodes"<<endl;
    cin>>n;
    num a,b;
    //Given
    a = 0.0;
    b = 3*M_PI;
    num step = b/(num) n;
    point given_points[n];
    cout << "Given point coordinates" << endl;
    //Points For constant distribution
    for (num i = a,  j = 0; i <= b,j<n-1; i+=step, j+=1) {
        given_points[(int) j].x =(num) i;
        cout << "X: " << given_points[(int) j].x;
        given_points[(int) j].y = f_x(i);
        cout << " Y: " << given_points[(int) j].y << endl;
    }
    //Points For chebyshev
    for(num i = a, j = 0 ; i <)
    given_points[n-1].x=(num)b;
    given_points[n-1].y=(num)f_x(b);
    ofstream lagrange, real_results,newton,interpolation_points;
    interpolation_points.open("/home/timelock/Pulpit/interpolation_points.csv");
    lagrange.open("/home/timelock/Pulpit/lagrange_resutlts.csv");
    newton.open("/home/timelock/Pulpit/newton_results.csv");
    real_results.open("/home/timelock/Pulpit/real_results.csv");

    interpolation_points<<"X,Y"<<endl;
    for(int i = 0 ;i < n ;i++){
        interpolation_points<<given_points[i].x<<","<<given_points[i].y<<endl;
    }
    lagrange << "X,Y" << endl;
    for (num i = a ; i <= b ; i+=0.01) {
        lagrange << i << "," << lagrange_interpolate(given_points, n, i) << endl;
    }

    newton << "X,Y" << endl;
    for (num i = a ; i <= b ; i+=0.01) {
        newton<< i << "," << newton_interpolate(given_points, n, i) << endl;
    }
    real_results << "X,Y" << endl;
    for (num i = a ; i <= b ; i+=0.01) {
        real_results << i << "," << f_x(i) << endl;
    }

    lagrange.close();
    real_results.close();
    return 0;
}