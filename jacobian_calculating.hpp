#include <cassert>
#include <vector> 
#include <adept_source.h>
#include "jacobian_coloring.hpp"

using adept::adouble;
using adept::aVector;

std::vector <double> vjp(adept::Stack &stack, aVector (*f)(const aVector &x),
	const aVector &x, aVector &y, const std::vector <double> &vy) {
	for (int i = 0; i < y.size(); i++) y[i].set_gradient(vy[i]);
	stack.compute_adjoint();
	std::vector <double> result(x.size());
	for (int i = 0; i < x.size(); i++) result[i] = x[i].get_gradient();
	return stack.clear_gradients(), result;
}

std::vector <double> jvp(adept::Stack &stack, aVector (*f)(const aVector &x),
	const aVector &x, aVector &y, const std::vector <double> &vx) {
	for (int i = 0; i < x.size(); i++) x[i].set_gradient(vx[i]);
	stack.compute_tangent_linear();
	std::vector <double> result(y.size());
	for (int i = 0; i < y.size(); i++) result[i] = y[i].get_gradient();
	return stack.clear_gradients(), result;
}

sparse_matrix sparse_jacobian_AD(aVector (*f)(const aVector &x),
	std::vector <jacobian_query> V, std::vector <double> _x, int n, int m) {
	adept::Stack stack; aVector x(n);
	sparse_matrix result;
	result.reserve(n);
	for (int i = 0; i < n; i++) x[i] = _x[i];
	stack.new_recording();
	aVector y = f(x);
	for (jacobian_query i : V) {
		if (i.is_row) {
			std::vector <double> vx(m, 0);
			for (int j : i.list) vx[j] = 1;
			std::vector <double> t = vjp(stack, f, x, y, vx);
			for (std::pair <int, int> j : i.grid)
				result.emplace_back(j.first, j.second, t[j.second]);
		}
		else {
			std::vector <double> vx(n, 0);
			for (int j : i.list) vx[j] = 1;
			std::vector <double> t = jvp(stack, f, x, y, vx);
			for (std::pair <int, int> j : i.grid)
				result.emplace_back(j.first, j.second, t[j.first]);
		}
	}
	return result;
}

sparse_matrix sparse_jacobian_FD(std::vector <double> (*f)(const std::vector <double> &x),
	std::vector <jacobian_query> V, std::vector <double> x, int n, int m, bool central = false, double epsilon = 1e-6) {
	sparse_matrix result;
	std::vector <double> y = f(x);
	for (jacobian_query i : V) {
		assert(!i.is_row);
		if (central) {
			std::vector <double> x1 = x, x2 = x;
			for (int j : i.list) x1[j] += epsilon, x2[j] -= epsilon;
			std::vector <double> y1 = f(x1), y2 = f(x2);
			for (std::pair <int, int> j : i.grid)
				result.emplace_back(j.first, j.second, (y1[j.first] - y2[j.first]) / (2 * epsilon));
		}
		else {
			std::vector <double> x1 = x;
			for (int j : i.list) x1[j] += epsilon;
			std::vector <double> y1 = f(x1);
			for (std::pair <int, int> j : i.grid)
				result.emplace_back(j.first, j.second, (y1[j.first] - y[j.first]) / epsilon); 
		}
	}
	return result;
}

