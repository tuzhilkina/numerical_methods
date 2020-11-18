#include <iostream>
#include <iomanip>
#include <string>

struct point {
	point() = default;
	point(const point&) = default;
	~point() = default;
	point& operator=(const point&) = default;

	double x{ 0 };
	double y{ 0 };

	double f() const {
		return x * x * x - 3 * x * y * y - x * x + y * y + x - 1;
	}

	double g() const {
		return 3 * x * x * y - y * y * y - 2 * x * y + y;
	}

	std::ostream& write_to(std::ostream& ostrm) const {
		std::string str = "(" + std::to_string(x) + ", " + std::to_string(y) + "),";
		ostrm << str;
		return ostrm;
	}
};


inline std::ostream& operator<<(std::ostream& ostrm, const point& rhs) {
	return rhs.write_to(ostrm);
}

int sign(const double& x) {
	if (x < 0) return -1;
	if (x == 0) return 0;
	if (x > 0) return 1;
}


double S(const point& P, const point& Q, const point& R) {
	double a1(Q.x - P.x);
	double a2(R.x - P.x);
	double b1(Q.y - P.y);
	double b2(R.y - P.y);
	return std::abs(0.5 * (a1 * b2 - a2 * b1));
}

bool same_z_signs(const point& L, const point& R) {
	if (sign(L.f()) != sign(R.f()))
		return false;
	if (sign(L.g()) != sign(R.g()))
		return false;
	return true;
}

bool same_quater_pos(const point& L, const point& R) {
	if (sign(L.x) != sign(R.x))
		return false;
	if (sign(L.y) != sign(R.y))
		return false;
	return true;
}

bool same_quadrant(const point& L, const point& R) {
	if (!same_z_signs(L, R))
		return false;
	if (!same_quater_pos(L, R))\
		return false;
	return true;
}

void delete_point(point& P, point& Q, point& R, point& M) {
	if (same_quadrant(P, M)) {
		P = Q;
		Q = R;
		R = M;
	}
	else if (same_quadrant(Q, M)) {
		Q = R;
		R = M;
	}
	else if (same_quadrant(R, M)) {
		R = M;
	}
	else {
		double s0 = S(Q, R, M);
		double s1 = S(P, R, M);
		double s2 = S(P, Q, M);

		if ((s0 >= s1) && (s0 >= s2)) {
			P = M;
		}
		else if ((s1 >= s0) && (s0 >= s2)) {
			Q = M;
		}
		else if ((s2 >= s0) && (s2 >= s1)) {
			R = M;
		}
	}
}

double false_pos(point P, point Q, point R, const double& eps) {
	point M0;
	size_t epoch(0);
	bool flag(true);

	{
		std::cout.width(5);	    std::cout << std::left << "#";
		std::cout.width(30);	std::cout << std::left << "P";
		std::cout.width(30);	std::cout << std::left << "Q";
		std::cout.width(30);	std::cout << std::left << "R";
		std::cout.width(30);	std::cout << std::left << "M" << "\n";
	}


	while (flag) {
		// Определить три точки, лежащие на поверхностях z=f(x,y) и z=g(x,y)
		double zfp = P.f();		double zfq = Q.f();		double zfr = R.f();
		double zgp = P.g();		double zgq = Q.g();		double zgr = R.g();

		// Коэффициенты уравнения плоскости πf
		double Af = (Q.y - P.y) * (zfr - zfp) - (R.y - P.y) * (zfq - zfp);
		double Bf = (R.x - P.x) * (zfq - zfp) - (Q.x - P.x) * (zfr - zfp);
		double Cf = (Q.x - P.x) * (R.y - P.y) - (R.x - P.x) * (Q.y - P.y);
		double Df = -P.x * Af - P.y * Bf - zfp * Cf;

		// Коэффициенты уравнения плоскости πg
		double Ag = (Q.y - P.y) * (zgr - zgp) - (R.y - P.y) * (zgq - zgp);
		double Bg = (R.x - P.x) * (zgq - zgp) - (Q.x - P.x) * (zgr - zgp);
		double Cg = (Q.x - P.x) * (R.y - P.y) - (R.x - P.x) * (Q.y - P.y);
		double Dg = -P.x * Ag - P.y * Bg - zgp * Cg;

		// Найти точку пересечения плоскостей πf, πg и плоскости z=0
		double dt = Af * Bg - Ag * Bf;
		point M{ (-Df * Bg + Dg * Bf) / dt, (-Af * Dg + Ag * Df) / dt };
		double zfm(M.f());
		double zgm(M.g());

		// Проверка выполнения критерия (|f(M_i) – f(M_(i-1))| + |g(M_i) – g(M_(i-1))|) < ε
		double dv = std::abs(zfm) + std::abs(zgm);
		M0 = M;

		flag = (std::abs(zfm) > eps || std::abs(zgm) > eps);
		if (epoch % 10 == 0 || flag == false)
		{
			std::cout.width(5);	    std::cout << std::left << epoch;
			std::cout.width(30);	std::cout << std::left << P;
			std::cout.width(30);	std::cout << std::left << Q;
			std::cout.width(30);	std::cout << std::left << R;
			std::cout.width(30);	std::cout << std::left << M << "\n";
		}

		delete_point(P, Q, R, M);
		++epoch;
	}
	return 0;
}

int main() {
	point P{ -2, 2 };
	point Q{ -2, -1 };
	point R{ 7, 2 };

	false_pos(P, Q, R, 0.01);
}
