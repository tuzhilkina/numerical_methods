#include <iostream>
#include <vector>
using namespace std;

struct point {
	point() = default;
	point(const point&) = default;
	~point() = default;
	point& operator=(const point&) = default;

	double x{ 0 };
	double y{ 0 };

	ostream& write_to(ostream& ostrm) const {
		ostrm << "(" << x << ", " << y << ")";
		return ostrm;
	}
};

inline ostream& operator<<(ostream& ostrm, const point& rhs) {
	return rhs.write_to(ostrm);
}

void Lagrange(const vector<point>& data) {
	size_t n(data.size() - 1);
	vector<double> denum;

	cout << "Интерполяционный многочлен Лагранжа:\nL_" << n + 1 << "(x) = ";
	for (size_t j(0); j <= n; ++j) {
		denum.push_back(1);
		for (size_t i(0); i <= n; ++i) {
			if (j == i) continue;

			if (data[i].x > 0) cout << "(x - " << data[i].x << ") * ";
			else if (data[i].x == 0) cout << "x * ";
			else cout << "(x + " << abs(data[i].x) << ") * ";

			denum.back() *= data[i].x - data[j].x;
		}
		(data[j].y > 0) ?
			(cout << data[j].y) :
			(cout << "(-" << abs(data[j].y) << ")");
		(denum.back() > 0) ?
			(cout << " / " << denum.back()) :
			(cout << " / (" << denum.back() << ")");
		if (j != n) cout << " + ";
	}
	cout << "\n";
}

double calc_L(const double& x) {
	double y((x + 1) * x * (x - 1) * (x - 2) * (x - 3) * (x - 4) * (-9.44) / 720
		+ (x + 2) * x * (x - 1) * (x - 2) * (x - 3) * (x - 4) * 1.77 / (-120)
		+ (x + 2) * (x + 1) * (x - 1) * (x - 2) * (x - 3) * (x - 4) * (-0.06) / 48
		+ (x + 2) * (x + 1) * x * (x - 2) * (x - 3) * (x - 4) * (-6.17) / (-36)
		+ (x + 2) * (x + 1) * x * (x - 1) * (x - 3) * (x - 4) * 5.48 / 48
		+ (x + 2) * (x + 1) * x * (x - 1) * (x - 2) * (x - 4) * 4.92 / (-120)
		+ (x + 2) * (x + 1) * x * (x - 1) * (x - 2) * (x - 3) * (-9.66) / 720);
	cout << "   (" << x << ", " << y << ")\n";
	return y;
}

void Newton(const vector<point>& data) {
	size_t n(data.size() - 1);

	cout << "Интерполяционный многочлен Ньютона:\nP_" << n + 1 << "(x) = ";
	cout << data[0].y << " + ";

	vector<double> delta;
	// Разделенные разности 0-го порядка совпадают со значениями функции в узлах
	for (const auto& it : data)
		delta.push_back(it.y);

	for (int m(1); m <= n; ++m) {
		// Разделенные разности k-го порядка определяются через разделенные разности порядка k-1
		vector<double> delta_new;
		for (int i(0); i <= n - m; ++i) {
			delta_new.push_back((delta[i + 1] - delta[i]) / (data[i + m].x - data[i].x));
		}
		delta = delta_new;

		(delta[0] > 0) ?
			(cout << delta[0] << " * ") :
			(cout << "(-" << abs(delta[0]) << ") * ");
		
		for (int i(0); i <= m - 1; ++i) {
			if (data[i].x > 0) cout << "(x - " << data[i].x << ")";
			else if (data[i].x == 0) cout << "x";
			else cout << "(x + " << abs(data[i].x) << ")";
			if (i != m - 1) cout << " * ";
		}
		if (m != n) cout << " + ";
	}
	cout << "\n\n\n\n\n";
}

double calc_R(const double& x) {
	double y(-9.44 + 11.21 * (x + 2) + (-6.52) * (x + 2) * (x + 1) + 1.46 * (x + 2) * (x + 1) * x 
		+ 0.553333 * (x + 2) * (x + 1) * x * (x - 1) 
		+ (-0.544083) * (x + 2) * (x + 1) * x * (x - 1) * (x - 2) 
		+ 0.202028 * (x + 2) * (x + 1) * x * (x - 1) * (x - 2) * (x - 3));
	cout << "   (" << x << ", " << y << ")\n";
	return y;
}

int main() {
	setlocale(LC_ALL, "Russian");

	vector<point> data{ {-2, -9.44}, {-1, 1.77}, {0, -0.06}, {1,-6.17},
		{2, 5.48}, {3, 4.92}, {4, -9.66} };
	vector<point> A{ {-1.5, -14.1014}, {-0.75, -0.931596}, {0, 0}, {0.75, 0.931596}, {1.5, 14.1014} };
	
	size_t n(data.size() - 1);
	cout << "\n\nСеточная функция:\n";
	for (size_t i(0); i <= n; ++i) {
		cout << data[i];
		(i != n) ? (cout << ", ") : (cout << "\n\n\n");
	}


	//Lagrange(data);

	//cout << "   Вычисление значения многочлена в точках:\n";
	//calc_L(-1.3);
	//calc_L(1.75);
	//cout << "\n\n\n\n";

	Newton(data);

	cout << "   Вычисление значения многочлена в точках:\n";
	calc_L(-1.3);
	calc_L(1.75);
	cout << "\n\n\n\n";
}
