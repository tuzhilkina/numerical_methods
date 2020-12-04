#include <iostream>
#include <vector>
#include <iomanip>
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
	vector<double> Jacobi(const double& eps);
	vector<double> Seidel(const double& eps);
	bool criteria();

	size_t n{ 0 };

private:
	double delta(const vector<double>& x0, const vector<double>& x1, size_t& epoch);

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
			printf("%10.5f   ", A(i, j));
		printf("|        |%10.5f   |\n", b[i - 1]);
	}
	printf("\n");
}

bool matrix::criteria() {
	for (size_t i(1); i <= n; ++i) {
		if (A(i, i) == 0) {
			printf("   Строка %zd: в матрице присутствует нулевой диагональный элемент!\n\n\n\n\n", i);
			return false;
		}
		double sum_ij(0);
		for (size_t j(1); j <= n; ++j) {
			if (i != j)
				sum_ij += abs(A(i, j));
		}
		if (abs(A(i, i)) <= sum_ij) {
			printf("   Строка %zd: не соблюдено условие диагонального преобладания!\n\n\n\n\n", i);
			return false;
		}
	}
}

double matrix::delta(const vector<double>& x0, const vector<double>& x1, size_t& epoch) {
	printf("\n   x%zd = (", epoch);
	for (const auto& i : x0)
		printf("%10.5f   ", i);
	++epoch;
	printf(")\n   x%zd = (", epoch);
	for (const auto& i : x1)
		printf("%10.5f   ", i);

	double max_delta(abs(x1[0] - x0[0]));
	printf(")\n   d%zd = (%10.5f   ", epoch, max_delta);
	for (size_t i(1); i < n; ++i) {
		double delta(abs(x0[i] - x1[i]));
		printf("%10.5f   ", delta);
		if (delta > max_delta)
			max_delta = delta;
	}
	printf(")\n");
	return max_delta;
}

vector<double> matrix::Jacobi(const double& eps) {
	vector<double> x1(b);
	for (auto& i : x1)
		i = 1.7E+308;
	if (!criteria()) return x1;

	vector<double> x0(b);

	bool stop(false);
	size_t epoch(0);
	while (!stop) {
		for (auto& i : x1)
			i = 0;

		for (size_t i(1); i <= n; ++i) {
			for (size_t j(1); j <= i - 1; ++j) {
				x1[i - 1] -= A(i, j) * x0[j-1] / A(i, i);
			}
			for (size_t j(i + 1); j <= n; ++j) {
				x1[i - 1] -= A(i, j) * x0[j-1] / A(i, i);
			}
			x1[i - 1] += B(i) / A(i, i);
		}

		if (delta(x0, x1, epoch) < eps)
			stop = true;
		if (epoch == 1000)
			stop = true;
		x0 = x1;
	}

	printf("\n   x = |");
	for (const auto& i : x1)
		printf("%10.5f   ", i);
	printf("|\n\n\n\n\n");
	return x1;
}

vector<double> matrix::Seidel(const double& eps) {
	vector<double> x0(b);
	vector<double> x1(b);
	for (auto& i : x1)
		i = 1.7E+308;

	if (!criteria()) return x1;

	bool stop(false);
	size_t epoch(0);
	x1 = b;
	while (!stop) {
		for (size_t i(1); i <= n; ++i) {
			x1[i - 1] = 0;
			for (size_t j(1); j <= i - 1; ++j) {
				x1[i - 1] -= A(i, j) * x1[j - 1] / A(i, i);
			}
			for (size_t j(i + 1); j <= n; ++j) {
				x1[i - 1] -= A(i, j) * x1[j - 1] / A(i, i);
			}
			x1[i - 1] += B(i) / A(i, i);
		}

		if (delta(x0, x1, epoch) < eps)
			stop = true;
		if (epoch == 1000)
			stop = true;
		x0 = x1;
	}

	printf("\n   x = |");
	for (const auto& i : x1)
		printf("%10.5f   ", i);
	printf("|\n\n\n\n\n");
	return x1;
}


int main() {
	setlocale(LC_ALL, "Russian");
	matrix A{ 5, vector<double>{
		22.87454,   6.35792,   5.66149,  -4.21491,   -0.575884,
		-3.81329,  28.3874,   -4.00311,   8.53389,   -9.16745,
		-2.59194,   9.72777, -28.67855,   1.45604,    4.5732,
		-5.41612,  -4.64095,   8.08344, -25.21043,   -0.0540179,
		-2.08411,   5.56688,   1.51708,   8.31599,  -24.0991851
	}, vector<double> { 2.12867, -3.81024, 9.84924, -1.08676, -8.5815 } };
	A.print();

	vector<double> j(A.Jacobi(0.001));
	vector<double> s(A.Seidel(0.001));
	printf("\n   xj = (");
	for (const auto& i : j)
		printf("%10.5f   ", i);
	printf(")\n   xs = (");
	for (const auto& i : s)
		printf("%10.5f   ", i);
	printf(")\n   d  = (");
	for (size_t i(0); i < j.size(); ++i) {
		double delta(abs(j[i] - s[i]));
		printf("%10.5f   ", delta);
	}
	printf(")\n");


	matrix B{ 5, vector<double>{
		22.87454,   6.35792,   5.66149,  -4.21491,   -0.575884,
		-3.81329,  28.3874,   -4.00311,   8.53389,   -9.16745,
		-2.59194,   9.72777, -38.67855,   1.45604,    4.5732,
		-5.41612,  -4.64095,   8.08344,  -5.21043,   -0.0540179,
		-2.08411,   5.56688,   1.51708,   8.31599,   -0.0991851
    }, vector<double> { 2.12867, -3.81024, 9.84924, -1.08676, -8.5815 } };
	B.print();
	B.Jacobi(0.01);

	matrix C{ 5, vector<double>{
		 22.87454,   6.35792,   5.66149,  -4.21491,   -0.575884,
		-3.81329,   0,        -4.00311,   8.53389,   -9.16745,
		-2.59194,   9.72777,  -8.67855,   1.45604,    4.5732,
		-5.41612,  -4.64095,   8.08344,  -5.21043,   -0.0540179,
		-2.08411,   5.56688,   1.51708,   8.31599,   -0.0991851
    }, vector<double> { 2.12867, -3.81024, 9.84924, -1.08676, -8.5815 } };
	C.print();
	C.Jacobi(0.01);
}