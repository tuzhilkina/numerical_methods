#include <iostream>
#include <vector>
#include <iomanip>
#include <string>

class matrix {
public:
	matrix() = default;
	matrix(const matrix&) = default;
	~matrix() = default;
	matrix& operator=(const matrix&) = default;
	matrix(const int& rows, const int& cols, const std::vector<double>& data_);

	double& at(const int& row_i, const int& col_i);
	const double& at(const int& row_i, const int& col_i) const;

	void print() const;
	void print(const std::vector<double>& b) const;
	void print(const int& n_rows) const;


	matrix LU_make() const;
	std::vector<double> LU_method(const std::vector<double>& b) const;
	std::vector<double> check(const std::vector<double>& x, const std::vector<double>& b) const;
	void cut_matrix();
	std::vector<double> JR_method(const std::vector<double>& b) const;
	std::vector<double> TD_method(const std::vector<double>& b) const;

private:
	void JR_step(std::vector<double>& b);
	int n_rows{ 0 };
	int n_cols{ 0 };
	std::vector<double> data;
};

matrix::matrix(const int& rows, const int& cols, const std::vector<double>& data_) : n_rows(rows), n_cols(cols), data(data_) {
	if (data.size() != n_rows * n_cols)
		throw std::invalid_argument("Change the input vector");
}
double& matrix::at(const int& row_i, const int& col_i) {
	if ((row_i > n_rows) || (col_i > n_cols) || (row_i <= 0) || (col_i <= 0)) {
		throw std::out_of_range("Invalid index");
	}
	return data[(row_i - 1) * n_cols + (col_i - 1)];
}
const double& matrix::at(const int& row_i, const int& col_i) const {
	if ((row_i > n_rows) || (col_i > n_cols) || (row_i <= 0) || (col_i <= 0)) {
		throw std::out_of_range("Invalid index");
	}
	return data[(row_i - 1) * n_cols + (col_i - 1)];
}

void matrix::print() const {
	for (int i(1); i <= n_rows; ++i) {
		for (int j(1); j <= n_cols; ++j)
			std::cout << std::setw(9) << std::setprecision(3) << at(i, j);
		std::cout << "\n";
	}
	std::cout << "\n";
}
void matrix::print(const std::vector<double>& b) const {
	for (int j(1); j <= n_cols; ++j)
		std::cout << std::setw(8) << "ai" << j;
	std::cout << std::setw(18) << std::setprecision(3) << "b\n";

	for (int i(1); i <= n_rows; ++i) {
		for (int j(1); j <= n_cols; ++j)
			std::cout << std::setw(9) << std::setprecision(3) << at(i, j);
		std::cout << std::setw(18) << std::setprecision(3) << b[i-1] << "\n";
	}
	std::cout << "\n";
}

matrix matrix::LU_make() const {
	// 1. Операция факторизации исходной матрицы
	matrix LU{ *this };
	for (int i(1); i <= n_rows; ++i)
		for (int j(1); j <= n_rows; ++j) {
			if ((i == 1) && (j == 1)) continue;
			// Формула 4.1
			if (i >= j)
				for (int s(1); s < j; ++s)
					LU.at(i, j) -= LU.at(i, s) * LU.at(s, j);
			// Формула 4.2
			else if (i < j) {
				for (int s(1); s < i; ++s)
					LU.at(i, j) -= LU.at(i, s) * LU.at(s, j);
				LU.at(i, j) /= LU.at(i, i);
			}
		}
	LU.print();
	return LU;
}
std::vector<double> matrix::LU_method(const std::vector<double>& b) const {
	matrix LU{ LU_make() };
	std::vector<double> x{ b };

	// 2. Решение системы Ly=b
	for (int i(1); i <= n_rows; ++i) {
		for (int j(1); j < i; ++j)
			x[i - 1] -= LU.at(i, j) * x[j - 1];
		x[i - 1] /= LU.at(i, i);
	}

	for (const auto& it : x)
		std::cout << it << " ";
	std::cout << "\n";

	// 3. Решение системы Ux=y
	for (int i(n_rows); i >= 1; --i) {
		for (int j(n_rows); j > i; --j) {
			x[i - 1] -= LU.at(i, j) * x[j - 1];
		}
	}

	std::cout << "x = ( ";
	for (const auto& it : x)
		std::cout << std::setw(6) << it << " ";
	std::cout << ")\n\n";

	check(x, b);

	return x;
}
std::vector<double> matrix::check(const std::vector<double>& x, const std::vector<double>& b) const {
	std::vector<double> b0{ x };
	for (auto& it : b0) it = 0;

	// Умножение исходной матрицы на найденный вектор x
	for (int i(1); i <= n_rows; ++i) {
		for (int j(1); j <= n_rows; ++j) {
			b0[i - 1] += at(i, j) * x[j - 1];
		}
	}

	std::cout << "  Ax  ";
	for (const auto& it : b0)
		std::cout << std::setw(12) << it;
	std::cout << "\n  b   ";
	for (const auto& it : b)
		std::cout << std::setw(12) << it;
	std::cout << "\n  Ax-b";
	
	// Сравнение результатов
	for (int i(0); i < b0.size(); ++i) {
		b0[i] -= b[i];
		std::cout << std::setw(12) << b0[i] << " ";
	}
	std::cout << "\n\n\n";
	return b0;
}

void matrix::JR_step(std::vector<double>& b) {
	// Промежуточный шаг
	for (int k(1); k <= n_rows - 1; ++k) {
		double c = at(1,     1) / sqrt(at(1, 1) * at(1, 1) + at(k + 1, 1) * at(k + 1, 1));
		double s = at(k + 1, 1) / sqrt(at(1, 1) * at(1, 1) + at(k + 1, 1) * at(k + 1, 1));

		for (int j(1); j <= n_cols; ++j) {
			matrix temp{ *this };
			at(k + 1, j) = -s * temp.at(1, j) + c * temp.at(k + 1, j);
			at(1, j) = c * temp.at(1, j) + s * temp.at(k + 1, j);
		}
		std::vector<double> temp{ b };
		b[1 - 1] = c * temp[1 - 1] + s * temp[k + 1 - 1];
		b[k + 1 - 1] = -s * temp[1 - 1] + c * temp[k + 1 - 1];
	}
	std::cout << "after JR_step\n";
	print(b);
}
void matrix::cut_matrix() {
	// Удаление 1-й строки и 1-го столбца
	std::vector<double> cut_data;
	for (int i(2); i <= n_rows; ++i)
		for (int j(2); j <= n_cols; ++j)
			cut_data.push_back(at(i,j));
	n_rows -= 1;
	n_cols -= 1;
	data = cut_data;
}
std::vector<double> matrix::JR_method(const std::vector<double>& b) const {
	matrix work{ *this };
	matrix result{ *this };
	for (auto& it : result.data)
		it = 0;

	std::vector<double> work_b{ b };
	std::vector<double> result_b{ b };
	for (auto& it : result_b)
		it = 0;

	work.print(work_b);

	for (int i(1); i <= 5; ++i) {
		// Преобразование i-го уравнения
		work.JR_step(work_b);

		// Сохранение получившегося уравнения в измененную матрицу
		int k(1);
		for (int j(i); j <= n_cols; ++j) {
			result.at(i, j) = work.at(1,k);
			k += 1;
		}
		result_b[i-1] = work_b[0];
		
		// Получение новой системы для обработки следующего уравнения
		work.cut_matrix();
		for (int j(1); j < work_b.size(); ++j)
			work_b[j - 1] = work_b[j];
		work_b.pop_back();

		std::cout << "Result of " << i << " iteration\n";
		result.print(result_b);
		std::cout << "---------------------------------------------------------------------------------------------\n";

	}

	// Решение преобразованной системы
	std::vector<double> x{ result_b };
	for (int i(n_rows); i >= 1; --i) {
		for (int j(n_rows); j > i; --j) {
			x[i - 1] -= result.at(i, j) * x[j - 1];
		}
		x[i - 1] /= result.at(i, i);
	}

	std::cout << "x = ( ";
	for (const auto& it : x)
		std::cout << it << " ";

	std::cout << ")\n\n";

	check(x, b);
	return x;
}

std::vector<double> matrix::TD_method(const std::vector<double>& b) const {
	// δ_1 = - d_1 / c_1 
	std::vector<double> delta{ -at(1, 2) / at(1, 1) };
	// λ_1 = r_1 / c_1
	std::vector<double> lambda{ b[0] / at(1, 1) };
	
	// Прямая прогонка
	for (int i(2); i <= n_rows - 1; ++i) {
		// δ_i =   - d_i          / (c_i      + b_i          * δ_(i-1))
		double d_i(- at(i, i + 1) / (at(i, i) + at(i, i - 1) * delta.back()));
		// λ_i =   (r_i      - b_i          * λ_(i-1))       / (c_i      + b_i          * δ_(i-1))
		double l_i((b[i - 1] - at(i, i - 1) * lambda.back()) / (at(i, i) + at(i, i - 1) * delta.back()));
		
		delta.push_back(d_i);
		lambda.push_back(l_i);
	}

	std::cout << "\n\n  delta  = ( ";
	for (const auto& it : delta)
		std::cout << std::setw(6) << it << " ";
	std::cout << ")\n";
	std::cout << "  lambda = ( ";
	for (const auto& it : lambda)
		std::cout << std::setw(6) << it << " ";
	std::cout << ")\n\n";

	std::vector<double> x{ b };
	// x_n   = (r_n      - b_n                    * λ_(n-1))      / (c_n                + b_n                    * δ_(n-1)0)
	x.back() = (b.back() - at(n_rows, n_cols - 1) * lambda.back()) / (at(n_rows, n_cols) + at(n_rows, n_cols - 1) * delta.back());

	// Обратная прогонка
	for (int i(n_rows - 1); i >= 1; --i) {
		x[i - 1] = delta[i - 1] * x[i] + lambda[i - 1];
	}

	std::cout << "  x = ( ";
	for (const auto& it : x)
		std::cout << std::setw(6) << it << " ";
	std::cout << ")\n\n";

	check(x, b);

	return x;
}


int main() {
	matrix A{ 5, 5, std::vector<double>{
		 2.87454,   6.35792,   5.66149,  -4.21491,   -0.575884,
		-3.81329,   8.3874,   -4.00311,   8.53389,   -9.16745,
		-2.59194,   9.72777,  -8.67855,   1.45604,    4.5732,
		-5.41612,  -4.64095,   8.08344,  -5.21043,   -0.0540179,
		-2.08411,   5.56688,   1.51708,   8.31599,   -0.0991851
	} };
	std::vector<double> b{ 2.12867, -3.81024, 9.84924, -1.08676, -8.5815 };
	matrix A1{ 5, 5, std::vector<double>{
	     2.87454,   6.35792,   0,         0,          0,
		-3.81329,   8.3874,   -4.00311,   0,          0,
		 0,         9.72777,  -8.67855,   1.45604,    0,
		 0,         0,         8.08344,  -5.21043,   -0.0540179,
		 0,         0,         0,         8.31599,   -0.0991851
    } };

	A.LU_method(b);
	A.JR_method(b); 
	A1.TD_method(b);
}