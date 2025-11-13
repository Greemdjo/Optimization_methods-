#include<iostream>
#include<vector>
#include <cmath>
#include "Point_lab3.h"
using std::cin;
using std::cout;
using std::vector;
using std::endl;

double function(Point X) {
	double x = X[0], y = X[1];
	return exp(2 + x*x + y*y) + x + y;
}
Point gradientF(Point X) {
	double x = X[0], y = X[1];
	Point result({ 1 + 2 * x * exp(2 + x * x + y * y), 1 + 2 * y * exp(2 + x * x + y * y) });	return result;
}
double Armicho_rule(Point x_k_plus_1, Point x_k, double alpha, double gamma, double theta,
    double(*function)(Point), Point(*gradientF)(Point), size_t k) {

    while (function(x_k - alpha * gradientF(x_k)) - function(x_k) >
        -gamma * alpha * pow(gradientF(x_k).norma(), 2)) {
        alpha = theta * alpha;
    }
    return alpha;
}

Point gradient_descent_method(double eps, double alpha, double gamma, double theta,
    Point x0, double(*function)(Point), Point(*gradientF)(Point)) 
    {
        size_t k = 0;
        Point x_k1 = x0, x_k = x0;

        while (gradientF(x_k).norma() > eps) {
            alpha = Armicho_rule(x_k1, x_k, alpha, gamma, theta, function, gradientF, k);
            x_k1 = x_k - alpha * gradientF(x_k);
            x_k = x_k1;
            k = k + 1;
        }
        cout << "K = " << k << endl;
        return x_k1;
    }

int main() {
	double alpha = 0.1;
	double gamma = 0.2;
	double theta = 0.2;
	double eps = 0.01;
	Point x0 = vector<double>{ 0, 0 };

    cout << "Method Gradient Descent:" << endl;
    Point point_min1 = gradient_descent_method(eps, alpha, gamma, theta, x0, function, gradientF);
    cout << "x1_min: " << point_min1[0] << " x2_min: " << point_min1[1] << endl;


}

