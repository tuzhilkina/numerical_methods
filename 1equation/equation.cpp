#include <iostream>
#include <math.h>
#include <iomanip>

#define min(a,b) ((a) <= (b) ? (a) : (b))
#define max(a,b) ((a) > (b) ? (a) : (b))

int sign(const double& x) {
    if (x < 0) return -1;
    if (x == 0) return 0;
    if (x > 0) return 1;
}

double func(const double& x) {
    return 2 * x / (x * x + 1) + exp(x);
}

double bisect(double a, double b, const double& eps) {
    int epoch(0);
    double x(0);
    std::cout << std::setw(5) << "Epoch"
        << std::setw(10) << "a"
        << std::setw(10) << "b"
        << std::setw(10) << "x" << "\n";
    while (b - a > eps) {
        x = (a + b) / 2;
        std::cout << std::setw(5) << epoch
            << std::setw(10) << a
            << std::setw(10) << b
            << std::setw(10) << x << "\n";
        (sign(func(a)) != sign(func(x))) ? b = x : a = x;
        ++epoch;
    }
    x = (a + b) / 2;
    std::cout << "\nx = " << x << "\n\n\n";
    return x;
}

double secant(double x1, double x2, const double& eps) {
    int epoch(0);
    double x3(0);
    std::cout << std::setw(5) << "Epoch"
        << std::setw(10) << "x1"
        << std::setw(10) << "x2"
        << std::setw(10) << "x" << "\n";
    while (abs(x2 - x1) > eps) {
        x3 = x2 - func(x2) * (x2 - x1) / (func(x2) - func(x1));
        std::cout << std::setw(5) << epoch
            << std::setw(10) << min(x1, x2)
            << std::setw(10) << max(x1, x2)
            << std::setw(10) << x3 << "\n";
        x1 = x2;
        x2 = x3;
        ++epoch;
    }
    std::cout << "\nx = " << x3 << "\n\n\n";
    return x3;
}

double false_pos(double x0, double x1, const double& eps) {
    std::cout << std::setw(5) << "Epoch"
        << std::setw(10) << "x0"
        << std::setw(10) << "x1"
        << std::setw(10) << "x" << "\n";
    double delx(x1 - x0), x(0);
    double f0(func(x0)), f1(func(x1));
    int epoch(0);
    while (delx > eps) {
        x = x1 * f0 / (f0 - f1) + x0 * f1 / (f1 - f0);
        double f(func(x));
        std::cout << std::setw(5) << epoch
            << std::setw(10) << x0
            << std::setw(10) << x1
            << std::setw(10) << x << "\n";
        if (f0 * f > 0) {
            delx = x - x0;
            x0 = x;
            f0 = f;
        }
        else {
            delx = x1 - x;
            x1 = x;
            f1 = f;
        }
        ++epoch;
    }
    std::cout << "\nx = " << x << "\n\n\n";
    return x;
}

int main() {
    bisect(-2, 0, 0.01);
    bisect(-2, 0, 0.001);
    bisect(-2, 0, 0.0001);

    secant(-2, 0, 0.01);
    secant(-2, 0, 0.001);
    secant(-2, 0, 0.0001);

    false_pos(-2, 0, 0.01);
    false_pos(-2, 0, 0.001);
    false_pos(-2, 0, 0.0001);
}