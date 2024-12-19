#ifndef _C5G72D_sn_sweeper_H_
#define _C5G72D_sn_sweeper_H_

#include "basic.hpp"
#include "C5G72D_sn_frame.hpp"

/*
	0 - dx : - dy : -
	1 - dx : + dy : -
	2 - dx : - dy : +
	3 - dx : + dy : +
*/

namespace C5G72D_sn_sweeper {
	vector_double_2d Y1, Y2;
	int d; 

	void sweep(int g, std::vector <double> X,
	 std::vector <double> &Y, int index, int m2, int n2,
	 int l1_s, int l1_t, int l1_d, int l2_s, int l2_t, int l2_d) {
		using namespace C5G72D_sn_frame;
		double px = 2 * direction_data[0][d], py = 2 * direction_data[1][d], pz = direction_data[3][d] / 8;
		Y1[index] = (index == m2 ? std::vector <double> (Ny, 0) : Y1[m2]);
		Y2[index] = (index == n2 ? std::vector <double> (Nx, 0) : Y2[n2]);
		for (int i = l1_s; i != l1_t + l1_d; i += l1_d)
			for (int j = l2_s; j != l2_t + l2_d; j += l2_d) {
				double u = px / dx[i], v = py / dy[j];
				double key = 2 * (X[id(i, j)] + u * Y1[index][j] + v * Y2[index][i]) / (sigma_t[block[i][j]][g] + u + v);
				Y[id(i, j)] += pz * key;
				Y1[index][j] = key - Y1[index][j];
				Y2[index][i] = key - Y2[index][i];
			}
	}

	std::vector <double> C5G72D_sn_sweep(std::vector <double> X, int g) {
		using namespace C5G72D_sn_frame;
		using namespace C5G72D_sn_sweeper;
		std::vector <double> integral(Ng, 0);
		Y1 = vector_double_2d(4, std::vector <double> (Ny)); // (i, j + 1/2)
		Y2 = vector_double_2d(4, std::vector <double> (Nx)); // (i + 1/2, j)
		for (d = 0; d < ND / 4; d++) {
			sweep(g, X, integral, 0, 0, 0, Nx - 1,      0, -1, Ny - 1,      0, -1);
			sweep(g, X, integral, 1, 0, 1,      0, Nx - 1,  1, Ny - 1,      0, -1);
			sweep(g, X, integral, 2, 2, 0, Nx - 1,      0, -1,      0, Ny - 1,  1);
			sweep(g, X, integral, 3, 2, 1,      0, Nx - 1,  1,      0, Ny - 1,  1);
		}
		return integral;
	}
}

#endif

