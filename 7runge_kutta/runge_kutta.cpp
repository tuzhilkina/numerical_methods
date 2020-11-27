#include <iostream>
#include <cmath>
#include <vector>
using namespace std;

double y_orig(const double& x, const double& y) {
    return x * y * y - log(x * x) + log(y);
}

double secant(const double& x, double x1, double x2, const double& eps) {
    double x3(0);
    while (abs(x2 - x1) > eps) {
        x3 = x2 - y_orig(x, x2) * (x2 - x1) / (y_orig(x, x2) - y_orig(x, x1));
        x1 = x2;
        x2 = x3;
    }
    return x3;
}

double f(const double& x, const double& y) {
	return (2 * y - x * y * y * y) / (x + 2 * x * x * y * y);
}

void RK(const double& a, const double& b, const double& h, const double& eps) {
	vector<double> x, y, k1, k2, k3, k4;
	double n((b - a) / h);

	x.push_back(a);
	y.push_back(secant(a, a, b, 0.1 * eps));
	printf("\n   n = %2d,   x = %.1f:   ", 0, x.back());
	printf("|%.3f - %.3f| = %.3f\n", y.back(), y.back(), abs(y.back() - y.back()));

	for (int i(1); i <= n; ++i) {
		x.push_back(x[0] + i * h);
		k1.push_back(h * f(x.back(), y.back()));
		k2.push_back(h * f(x.back() + h / 2, y.back() + h * k1.back() / 2));
		k3.push_back(h * f(x.back() + h / 2, y.back() + h * k2.back() / 2));
		k4.push_back(h * f(x.back() + h, y.back() + k3.back()));
		y.push_back(y.back() + (k1.back() + 2 * k2.back() + 2 * k3.back() + k4.back()) / 6);

		double y0(secant(x.back(), a, b, 0.1 * eps));

		printf("   n = %2d,   x = %.1f:   ", i, x.back());
		printf("|%.3f - %.3f| = %.3f\n", y.back(), y0, abs(y.back() - y0));
	}
}

int main() {
	double a(1), b(2), h(0.05);
	RK(a, b, h, 0.01);
}