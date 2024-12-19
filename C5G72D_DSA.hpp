#ifndef _C5G72D_DSA_H
#define _C5G72D_DSA_H_

#include "basic.hpp"
#include "C5G72D_sn_frame.hpp"
#include "jacobian_calculating.hpp"
#include "linear_equation_solver.hpp"
#include "order.hpp"

#include <vector>

namespace C5G72D_DSA {
	using namespace C5G72D_sn_frame;

	vector_double_3d Dx, Dy;
	linear_equation_solver_Ck <float> linear_solver;
	Order order;

	int T;

	int current_g, Nx0, Ny0, Ng0;
	std::vector <double> dx0, dy0;
	std::vector <std::vector <int> > block0;

	int id0(int x, int y) {return x * Ny0 + y;}

	sparse_structure get_diffusion_function_sparse() {
		sparse_structure G(Ng0, std::vector <int> ());
		for (int i = 0; i < Nx0; i++)
			for (int j = 0; j < Ny0; j++) {
				G[id0(i, j)].push_back(id0(i, j));
				if (i != 0) G[id0(i, j)].push_back(id0(i - 1, j));
				if (i + 1 != Nx0) G[id0(i, j)].push_back(id0(i + 1, j));
				if (j != 0) G[id0(i, j)].push_back(id0(i, j - 1));
				if (j + 1 != Ny0) G[id0(i, j)].push_back(id0(i, j + 1));
			}
		return G;
	}

	template <typename Type1, typename Type2>
	Type1 diffusion_function_basic(const Type1 &x) {
		assert(x.size() == Ng0);
		Type1 y(Ng0);
		for (int i = 0; i < Nx0; i++)
			for (int j = 0; j < Ny0; j++) {
				y[id0(i, j)] = (sigma_t[block0[i][j]][current_g] - sigma_s[block0[i][j]][current_g][current_g]) * x[id0(i, j)];
				double dL = dx0[i] * dx0[i], dR = dx0[i] * dx0[i], dU = dy0[j] * dy0[j], dD = dy0[j] * dy0[j];
				if (i > 0) dL = ((dx0[i - 1] + dx0[i]) / 2) * dx0[i];
				if (i + 1 < Nx) dR = ((dx0[i + 1] + dx0[i]) / 2) * dx0[i];
				if (j > 0) dU =  ((dy0[j - 1] + dy0[j]) / 2) * dy0[j];
				if (j + 1 < Ny) dD = ((dy0[j + 1] + dy0[j]) / 2) * dy0[j];
				if (i != 0) y[id0(i, j)] += Dx[current_g][i - 1][j] * (x[id0(i, j)] - x[id0(i - 1, j)]) / dL;
				if (j != 0) y[id0(i, j)] += Dy[current_g][i][j - 1] * (x[id0(i, j)] - x[id0(i, j - 1)]) / dU;
				y[id0(i, j)] += (Dx[current_g][i][j] / dR + Dy[current_g][i][j] / dD) * x[id0(i, j)];
				if (i + 1 != Nx0) y[id0(i, j)] -= Dx[current_g][i][j] * x[id0(i + 1, j)] / dR;
				if (j + 1 != Ny0) y[id0(i, j)] -= Dy[current_g][i][j] * x[id0(i, j + 1)] / dD;
			}
		return y;
	}
	auto diffusion_function = diffusion_function_basic <std::vector <double>, double>;
	auto diffusion_function_AD = diffusion_function_basic <aVector, adouble>; 

	void calc_Dx_Dy() {
		for (int g = 0; g < G; g++)
			for (int i = 0; i < Nx0; i++)
				for (int j = 0; j < Ny0; j++)
					Dx[g][i][j] = Dy[g][i][j] = 1 / (3 * sigma_t[block0[i][j]][g]); 
	}

	void init() {
		std::ignore = scanf("%d", &T);
		Nx0 = Nx / T, Ny0 = Ny / T; Ng0 = Nx0 * Ny0;
		dx0 = std::vector <double> (Nx0, 0);
		dy0 = std::vector <double> (Ny0, 0);
		for (int i = 0; i < Nx0; i++)
			for (int j = 0; j < T; j++)
				dx0[i] += dx[T * i + j];
		for (int i = 0; i < Ny0; i++)
			for (int j = 0; j < T; j++)
				dy0[i] += dy[T * i + j];
		block0 = std::vector <std::vector <int> > (Nx0, std::vector <int> (Ny0));
		for (int i = 0; i < Nx0; i++)
			for (int j = 0; j < Ny0; j++)
				block0[i][j] = block[T * i][T * j];
		Dx = vector_double_3d(G, vector_double_2d(Nx0, std::vector <double> (Ny0)));
		Dy = vector_double_3d(G, vector_double_2d(Nx0, std::vector <double> (Ny0)));
		linear_solver.init_K(G);
		order = Order(order_2D_square::get_order(Nx0, Ny0));
		sparse_structure sparse = get_diffusion_function_sparse();
		printf("init structure!\n");
		linear_solver.init_structure(Ng0, order.X(sparse));
		std::vector <jacobian_query> query = jacobian_coloring::solve(Ng0, Ng0, sparse, false);
		calc_Dx_Dy();
		for (current_g = 0; current_g < G; current_g++) {
			printf("init_value : g = %d\n", current_g);
			sparse_matrix matrix = sparse_jacobian_AD(diffusion_function_AD, query, std::vector <double> (Ng0, 0), Ng0, Ng0);
			linear_solver.init_value(order.X(matrix), current_g);
		}
		printf("end.\n");
	}

	std::vector <double> get_err(vector_double_2d &err, vector_double_2d &delta_X, int g) {
		std::vector <double> result = err[g];
		for (int g0 = 0; g0 < G; g0++) {
			if (g0 == g) continue;
			for (int i = 0; i < Nx0; i++)
				for (int j = 0; j < Ny0; j++)
					result[id0(i, j)] += delta_X[g0][id0(i, j)] * sigma_s[block0[i][j]][g0][g];
		}
		std::vector <double> tmp = diffusion_function(delta_X[g]);
		for (int i = 0; i < Ng0; i++) result[i] -= tmp[i];
		return result;
	}

	double check(std::vector <double> current_err) {
		double sum = 0;
		for (int i = 0; i < Ng0; i++) sum += current_err[i] * current_err[i];
		return sqrt(sum / Ng0);
	}

	vector_double_2d compress(vector_double_2d A) {
		vector_double_2d B(G, std::vector <double> (Ng0, 0));
		for (int g = 0; g < G; g++)
			for (int i = 0; i < Nx0; i++)
				for (int j = 0; j < Ny0; j++) {
					for (int ii = 0; ii < T; ii++)
						for (int jj = 0; jj < T; jj++)
							B[g][id0(i, j)] += A[g][id(T * i + ii, T * j + jj)];
					B[g][id0(i, j)] /= (T * T);
				}
		return B;
	}

	vector_double_2d interpolate(vector_double_2d A) {
		vector_double_2d B(G, std::vector <double> (Ng));
		for (int g = 0; g < G; g++)
			for (int i = 0; i < Nx0; i++)
				for (int j = 0; j < Ny0; j++)
					for (int ii = 0; ii < T; ii++)
						for (int jj = 0; jj < T; jj++)
							B[g][id(T * i + ii, T * j + jj)] = A[g][id0(i, j)];
		return B;
	}

	void DSA_iterate(vector_double_2d &integral_X, double k, vector_double_2d &err) {
		vector_double_2d err0 = compress(err);
		vector_double_2d delta_X0(G, std::vector <double> (Ng0, 0));
		double maxn = 0;
		for (int i = 0; i < G; i++)
			for (double j : err0[i])
				if (std::fabs(j) > maxn) maxn = std::fabs(j);
		for (int i = 0; i < G; i++)
			printf("(%+.0le %+.0le)%c",
			 *std::min_element(err0[i].begin(), err0[i].end()),
			 *std::max_element(err0[i].begin(), err0[i].end()),
			 i + 1 == G ? '\n' : ' ');
		int cnt = 0;
		while (true) {
			double max_eps = 1; bool update = false;
			for (current_g = 0; current_g < G; current_g++) {
				std::vector <double> current_err = get_err(err0, delta_X0, current_g);
				if (check(current_err) < maxn * 1e-3) continue;
				printf("%3d: g = %d err = %.5le -> ", ++cnt, current_g, check(current_err));
				std::vector <double> tmp = order.Y(linear_solver.solve(order.X(current_err), current_g));
				for (int i = 0; i < Ng0; i++) delta_X0[current_g][i] += tmp[i];
				printf("%.5le\n", check(get_err(err0, delta_X0, current_g)));
				update = true;
			}
			if (!update) break;
		}
		for (int i = 0; i < G; i++)
			printf("(%+.0le %+.0le)%c",
			 *std::min_element(delta_X0[i].begin(), delta_X0[i].end()),
			 *std::max_element(delta_X0[i].begin(), delta_X0[i].end()),
			 i + 1 == G ? '\n' : ' ');
		vector_double_2d delta_X = interpolate(delta_X0);
		for (int i = 0; i < G; i++)
			for (int j = 0; j < Ng; j++)
				integral_X[i][j] -= delta_X[i][j];
		printf("DSA end!\n");
	}
}

#endif
