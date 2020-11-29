#include <iostream>
#include <math.h>
#include <vector>
using namespace std;

double trapezoid(const double& a, const double& b, const size_t& n) {
	cout << "\n";
	double h((b - a) / n);
	vector<double> x, y;
	x.push_back(a);
	double sum_y(0);

	for (size_t i(0); i <= n; ++i)
		x.push_back(h + x.back());
	for (size_t i(0); i <= n; ++i) {
		y.push_back(sqrt(x[i] * x[i] + 2 * x[i] + 5));
		printf("   n = %zd:    x = %10.7f   y = %10.7f\n", i, x[i], y[i]);
		if ((i != 0) && (i != n))
			sum_y += y.back();
	}

	double I(h * ((y[0] + y[n]) / 2 + sum_y));
	printf("\n   I = %8.5f\n\n", I);

	return I;
}

double simpson(const double& a, const double& b, const size_t& n) {
	cout << "\n";
	double h((b - a) / n);
	vector<double> x, y;
	x.push_back(a);

	for (size_t i(0); i <= n; ++i)
		x.push_back(h + x.back());
	for (size_t i(0); i <= n; ++i) {
		y.push_back(sqrt(x[i] * x[i] + 2 * x[i] + 5));
		printf("   n = %zd:    x = %10.7f   y = %10.7f\n", i, x[i], y[i]);
	}

	double sigma1(0), sigma2(0);
	for (size_t i(1); i <= n - 1; i += 2)
		sigma1 += y[i];
	for (size_t i(2); i <= n - 2; i += 2)
		sigma2 += y[i];

	double I((h / 3) * (y[0] + y[n] + 4 * sigma1 + 2 * sigma2));
	printf("\n   I = %10.7f\n\n", I);

	return I;
}

double gauss(const double& a, const double& b, const size_t& n) {
	cout << "\n";
	vector<double> t{ -0.906180, -0.538469, 0, 0.538469, 0.906180 };
	vector<double> A { 0.236927, 0.478629, 0.568889, 0.478629, 0.236927 };

	double I(0);
	vector<double> x;
	for (size_t i(1); i <= n; ++i) {
		x.push_back((a + b) / 2 + ((b - a) * t[i - 1]) / 2);
		double f(sqrt(x.back() * x.back() + 2 * x.back() + 5));
		printf("   n = %zd:    x = %10.7f   f = %10.7f\n", i, x.back(), f);
		I += A[i - 1] * f;
	}
	I = (b - a) * I / 2;
	printf("\n   I = %10.7f\n\n", I);

	return I;
}

int main() {
	setlocale(LC_ALL, "Russian");
	double It(trapezoid(-1, 1, 6));
	double Is(simpson(-1, 1, 6));
	double Ig(gauss(-1, 1, 5));

	double I(4.59117429878528);
	printf("\n\n%10.7f   %10.7f   %10.7f", It, Is, Ig);
	printf("\n\n%10.7f   %10.7f   %10.7f", abs(It - I), abs(Is - I), abs(Ig - I));

}