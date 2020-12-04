#define _USE_MATH_DEFINES // for C++
#include <cmath>
#include <iostream>
#include <vector>
#include <string>
using namespace std;

#define max(a,b) ((a) > (b) ? (a) : (b))

class matrix {
public:
	explicit matrix() = default;
	matrix(const matrix&) = default;
	~matrix() = default;
	matrix& operator=(const matrix&) = default;
	matrix(const size_t& n_, const vector<double>& data_, const vector<double>& b_);

	double& A(const size_t& row_i, const size_t& col_i);
	const double& A(const size_t& row_i, const size_t& col_i) const;
	double& B(const size_t& i);
	const double& B(const size_t& i) const;

	void print() const;

	vector<double> TD_method() const;

	size_t n{ 0 };

private:

	vector<double> data;
	vector<double> b;
};

matrix::matrix(const size_t& n_, const vector<double>& data_, const vector<double>& b_)
	: n(n_)
	, data(data_)
	, b(b_) {
	if (data.size() != n * n) {
		printf("Неверные количество элементов матрицы или ее размер!\n");
		throw invalid_argument("Change the input vector");
	}
	if (b.size() != n) {
		printf("Неверное количество элементов вектора!\n");
		throw invalid_argument("Change the input vector");
	}
}

double& matrix::A(const size_t& row_i, const size_t& col_i) {
	if ((row_i > n) || (col_i > n) || (row_i <= 0) || (col_i <= 0))
		throw out_of_range("Invalid index");
	return data[(row_i - 1) * n + (col_i - 1)];
}
const double& matrix::A(const size_t& row_i, const size_t& col_i) const {
	if ((row_i > n) || (col_i > n) || (row_i <= 0) || (col_i <= 0))
		throw out_of_range("Invalid index");
	return data[(row_i - 1) * n + (col_i - 1)];
}

double& matrix::B(const size_t& i) {
	if ((i > n) || (i <= 0))
		throw out_of_range("Invalid index");
	return b[i - 1];
}
const double& matrix::B(const size_t& i) const {
	if ((i > n) || (i <= 0))
		throw out_of_range("Invalid index");
	return b[i - 1];
}

void matrix::print() const {
	printf("\n");
	for (size_t i(1); i <= n; ++i) {
		printf("   |");
		for (size_t j(1); j <= n; ++j)
			printf("%6.3f  ", A(i, j));
		printf("|      |%6.3f |\n", b[i - 1]);
	}
	printf("\n");
}

vector<double> matrix::TD_method() const {
	// δ_1 = - d_1 / c_1 
	vector<double> delta{ -A(1, 2) / A(1, 1) };
	// λ_1 = r_1 / c_1
	vector<double> lambda{ B(1) / A(1, 1) };

	// Прямая прогонка
	for (size_t i(2); i <= n - 1; ++i) {
		// δ_i =   - d_i        / (c_i     + b_i        * δ_(i-1))
		double d_i(-A(i, i + 1) / (A(i, i) + A(i, i - 1) * delta.back()));
		// λ_i =   (r_i  - b_i         * λ_(i-1))       / (c_i     + b_i         * δ_(i-1))
		double l_i((B(i) - A(i, i - 1) * lambda.back()) / (A(i, i) + A(i, i - 1) * delta.back()));

		delta.push_back(d_i);
		lambda.push_back(l_i);
	}

	vector<double> x{ b };
	// x_n   = (r_n      - b_n         * λ_(n-1))       / (c_n     + b_n         * δ_(n-1)0)
	x.back() = (b.back() - A(n, n - 1) * lambda.back()) / (A(n, n) + A(n, n - 1) * delta.back());

	// Обратная прогонка
	for (size_t i(n - 1); i >= 1; --i) 
		x[i - 1] = delta[i - 1] * x[i] + lambda[i - 1];

	return x;
}


double P(const double& x) {
	return 0;
}

double Q(const double& x) {
	return 1;
}

double F(const double& x) {
	return 1;
}


vector<double> finite_dif(const double& a, const double& b, const size_t& n) {
	vector<double> x, p, q, f;
	double h((b - a) / n);
	x.push_back(a);
	for (size_t i(1); i <= n; ++i)
		x.push_back(x[0] + i * h);
	for (size_t i(0); i <= n; ++i) {
		p.push_back(P(x[i]));
		f.push_back(F(x[i]));
		q.push_back(Q(x[i]));
	}

	vector<double> data;
	vector<double> fi;

	double alpha0(1), alpha1(0), A(0), beta0(1), beta1(0), B(0);

	double S((2 - h * h * q[1]) / (1 - h * p[1] / 2));
	double T(-(1 + h * p[1] / 2) / (1 - h * p[1] / 2));
	double U(h * h * f[1] / (1 - h * p[1] / 2));
	double temp(2 * h * alpha0 - 3 * alpha1);
	double L(temp * S + 4 * alpha1);
	double M(temp * T - alpha1);
	double N(2 * A * h - temp * U);

	double V(-(1 - h * p[n - 1] / 2) / (1 + h * p[n - 1] / 2));
	double W((2 - h * h * q[n - 1]) / (1 + h * p[n - 1] / 2));
	double Z(h * h * f[n - 1] / (1 + h * p[n - 1] / 2));
	temp = (2 * h * beta0 + 3 * beta1);
	double P(temp * V + beta1);
	double Q(temp * W - 4 * beta1);
	double R(2 * B * h - temp * Z);

	size_t m(n - 1);
	data.push_back(L);	data.push_back(M); fi.push_back(N);
	for (size_t i(1); i <= m - 2; ++i)
		data.push_back(0);

	for (size_t i(1); i <= m - 2; ++i) {
		for (size_t j(1); j <= i - 1; ++j)
			data.push_back(0);
		data.push_back(1 - h * p[i + 1] / 2);
		data.push_back(-(2 - h * h * q[i]));
		data.push_back(1 + h * p[i - 1] / 2);
		for (size_t j(1); j <= m - i - 2; ++j)
			data.push_back(0);

		fi.push_back(h * h * f[i]);
	}
	for (size_t i(1); i <= m - 2; ++i)
		data.push_back(0);
	data.push_back(P); data.push_back(Q); fi.push_back(R);


	matrix SLAE{ m, data, fi };
	SLAE.print();
	vector<double> y0(SLAE.TD_method());
	vector<double> y;
	y.push_back(S * y0[0] + T * y0[1] + U);
	y.insert(y.end(), y0.begin(), y0.end());
	y.push_back(V * y0[n - 3] + W * y0[n - 2] + Z);

	for (size_t i(0); i < y.size(); ++i) {
		printf("(%.5f, %.5f)", x[i], y[i]);
		if (i != y.size() - 1) cout << ", ";
	}
	printf("\n\n\n");
	return y;
}


int main() {
	vector<double> yn(finite_dif(0, M_PI, 6));
	vector<double> y2n(finite_dif(0, M_PI, 12));

	for (size_t i(0); i < yn.size(); ++i) {
		printf("   n = %2zd:    ", i);
		printf("|%8.5f - %8.5f| = %.5f\n", yn[i], y2n[2*i], abs(yn[i] - y2n[2 * i]));
	}
}

