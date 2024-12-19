#include "basic.hpp"
#include "C5G72D_sn_frame.hpp"
#include "C5G72D_sn_sweeper.hpp"
#include "C5G72D_DSA.hpp"
#include "anderson.hpp"

Timer timer;

namespace C5G72D_source_iteration {
	bool converge;
	double epsilon;
	anderson <float> AA;

	double calc_error(std::vector <double> X, std::vector <double> Y) {
		double result = 0;
		assert(X.size() == Y.size());
		for (int i = 0; i < X.size(); i++)
			result += (X[i] - Y[i]) * (X[i] - Y[i]);
		return sqrt(result / X.size());
	}

	double calc_error_2d(vector_double_2d X, vector_double_2d Y) {
		double result = 0; int num = 0;
		for (int i = 0; i < X.size(); i++) {
			num += X[i].size();
			for (int j = 0; j < X[i].size(); j++)
				result += (X[i][j] - Y[i][j]) * (X[i][j] - Y[i][j]);
		}
		return sqrt(result / num);
	}

	std::vector <double> fold(vector_double_2d integral_X, double k) {
		std::vector <double> T;
		for (int g = 0; g < integral_X.size(); g++)
			T.insert(T.end(), integral_X[g].begin(), integral_X[g].end());
		return T.push_back(1 / k), T;
	}

	std::pair <vector_double_2d, double> unfold(std::vector <double> T) {
		using namespace C5G72D_sn_frame; 
		vector_double_2d integral_X; double k;
		for (int g = 0; g < G; g++)
			integral_X.push_back(std::vector <double> (T.begin() + g * Ng, T.begin() + (g + 1) * Ng));
		return std::make_pair(integral_X, 1 / T.back());
	}

	std::vector <double> get_scattering(vector_double_2d integral_X, int g) {
		using namespace C5G72D_sn_frame;
		std::vector <double> scattering(Ng, 0);
		for (int i = 0; i < Nx; i++)
			for (int j = 0; j < Ny; j++)
				for (int g1 = 0; g1 < G; g1++)
					scattering[id(i, j)] += sigma_s[block[i][j]][g1][g] * integral_X[g1][id(i, j)];
		return scattering;
	}

	std::vector <double> get_fission(vector_double_2d integral_X, int g, double k = 1) {
		using namespace C5G72D_sn_frame;
		std::vector <double> fission(Ng, 0);
		for (int i = 0; i < Nx; i++)
			for (int j = 0; j < Ny; j++) {
				for (int g0 = 0; g0 < G; g0++)
					fission[id(i, j)] += nu[block[i][j]][g0] * sigma_f[block[i][j]][g0] * integral_X[g0][id(i, j)];
				fission[id(i, j)] *= chi[g] / k;
			}
		return fission;
	}

	vector_double_2d get_fission(vector_double_2d integral_X, double k = 1) {
		using namespace C5G72D_sn_frame;
		vector_double_2d fission(G);
		for (int g = 0; g < G; g++)
			fission[g] = get_fission(integral_X, g, k);
		return fission;
	}

	double calc_fission_sum(vector_double_2d &integral_X) {
		using namespace C5G72D_sn_frame;
		double sum = 0;
		for (int i = 0; i < Nx; i++)
			for (int j = 0; j < Ny; j++)
				for (int g = 0; g < G; g++)
					sum += nu[block[i][j]][g] * sigma_f[block[i][j]][g] * integral_X[g][id(i, j)];
		return sum;
	}

	std::vector <double> combine(std::vector <double> X, std::vector <double> Y) {
		for (int i = 0; i < X.size(); i++) X[i] += Y[i];
		return X;
	}

	vector_double_2d combine(vector_double_2d X, vector_double_2d Y) {
		for (int i = 0; i < X.size(); i++)
			for (int j = 0; j < X[i].size(); j++)
				X[i][j] += Y[i][j];
		return X; 
	}

	std::vector <double> get_source(vector_double_2d integral_X, int g, double k) {
		return combine(get_scattering(integral_X, g), get_fission(integral_X, g, k));
	}

	void inner_iterate(vector_double_2d &integral_X, double k, vector_double_2d &err) {
		using namespace C5G72D_sn_frame;
		err.resize(G);
		vector_double_2d fission = get_fission(integral_X, k);
		for (int g = 0; g < G; g++) {
			err[g] = get_scattering(integral_X, g);
			integral_X[g] = C5G72D_sn_sweeper::C5G72D_sn_sweep(combine(err[g], fission[g]), g);
			std::vector <double> ().swap(fission[g]);
		}
		for (int g = 0; g < G; g++) {
			std::vector <double> tmp = get_scattering(integral_X, g);
			for (int i = 0; i < Ng; i++) err[g][i] -= tmp[i];
		}
	}

	void iterate(vector_double_2d &integral_X, double &k) {
		double tmp = calc_fission_sum(integral_X);
		vector_double_2d err;
		inner_iterate(integral_X, k, err);
		printf("DSA_iterate! time = %.2lfs\n", timer.time());
		C5G72D_DSA::DSA_iterate(integral_X, k, err);
		k *= calc_fission_sum(integral_X) / tmp;
	}

	std::vector <double> iterate(std::vector <double> T) {
		vector_double_2d integral_X; double k;
		std::tie(integral_X, k) = unfold(T);
		iterate(integral_X, k);
		return fold(integral_X, k);
	}

	void solve(vector_double_2d &integral_X, double &k, double _epsilon, int max_iter) {
		C5G72D_DSA::init();
		converge = false; epsilon = _epsilon;
		std::vector <double> T = fold(integral_X, k);
		AA.init(T.size(), 8, T, iterate(T));
		AA.set_error_width(T.size() - 1);
		printf("iter %3d : k = %.10le error = %.6le time = %.2lfs\n", 0, k, AA.error_list.back(), timer.time());
		for (int iter = 1; iter <= max_iter; iter++) {
			T = AA.calc_next();
			std::vector <double> fT = iterate(T);
			AA.update(T, fT);
			printf("iter %3d : k = %.10le error = %.6le time = %.2lfs\n", iter, 1 / T.back(), AA.error_list.back(), timer.time());
			if (AA.error_list.back() < epsilon) {T = fT; converge = true; break;}
		}
		std::tie(integral_X, k) = unfold(T);
	}

	void output_power_distribution(const char *file_name, vector_double_2d &integral_X, int cx, int cy, double k) {
		using namespace C5G72D_sn_frame; 
		static const int L = 34, total = 1056;
		
		static const double reference_1[L / 2][L / 2] = {
			{2.203, 2.207, 2.217, 2.229, 2.232, 2.236, 2.189, 2.151, 2.125, 2.058, 2.002, 1.953, 1.857, 1.756, 1.634, 1.485, 1.283},
			{0.000, 2.218, 2.244, 2.280, 2.307, 2.378, 2.257, 2.216, 2.255, 2.121, 2.064, 2.077, 1.921, 1.798, 1.655, 1.493, 1.282},
			{0.000, 0.000, 2.317, 2.450, 2.476, 0.000, 2.389, 2.341, 0.000, 2.242, 2.183, 0.000, 2.064, 1.934, 1.714, 1.513, 1.287},
			{0.000, 0.000, 0.000, 0.000, 2.503, 2.464, 2.298, 2.248, 2.298, 2.152, 2.100, 2.150, 2.082, 0.000, 1.816, 1.541, 1.294},
			{0.000, 0.000, 0.000, 0.000, 2.404, 2.434, 2.281, 2.234, 2.284, 2.139, 2.086, 2.126, 2.002, 1.975, 1.834, 1.562, 1.298},
			{0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 2.354, 2.308, 0.000, 2.211, 2.153, 0.000, 2.031, 1.944, 0.000, 1.615, 1.303},
			{0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 2.220, 2.180, 2.231, 2.088, 2.033, 2.061, 1.905, 1.816, 1.771, 1.532, 1.277},
			{0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 2.141, 2.192, 2.052, 1.997, 2.023, 1.868, 1.780, 1.739, 1.507, 1.258},
			{0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 2.102, 2.045, 0.000, 1.914, 1.822, 0.000, 1.541, 1.248},
			{0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 1.969, 1.917, 1.943, 1.795, 1.711, 1.674, 1.452, 1.213},
			{0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 1.868, 1.895, 1.754, 1.675, 1.636, 1.419, 1.186},
			{0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 1.794, 1.721, 0.000, 1.439, 1.166},
			{0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 1.692, 1.675, 1.565, 1.338, 1.117},
			{0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 1.479, 1.263, 1.069},
			{0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 1.319, 1.176, 1.014},
			{0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 1.080, 0.954},
			{0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.878}
		};

		static const double reference_2[L / 2][L / 2] = {
			{1.312, 1.061, 0.936, 0.864, 0.811, 0.768, 0.712, 0.662, 0.620, 0.569, 0.523, 0.484, 0.437, 0.398, 0.377, 0.409, 0.597},
			{1.295, 1.342, 1.170, 1.092, 1.045, 1.044, 0.916, 0.848, 0.838, 0.731, 0.670, 0.656, 0.565, 0.503, 0.469, 0.515, 0.589},
			{1.289, 1.320, 1.177, 1.175, 1.118, 0.000, 0.950, 0.870, 0.000, 0.755, 0.688, 0.000, 0.606, 0.542, 0.474, 0.506, 0.586},
			{1.292, 1.332, 1.263, 0.000, 1.111, 1.128, 0.955, 0.878, 0.876, 0.757, 0.696, 0.703, 0.587, 0.000, 0.513, 0.512, 0.585},
			{1.294, 1.357, 1.288, 1.180, 1.140, 1.083, 0.926, 0.852, 0.852, 0.735, 0.675, 0.677, 0.608, 0.537, 0.517, 0.522, 0.585},
			{1.297, 1.425, 0.000, 1.281, 1.142, 0.000, 0.984, 0.901, 0.000, 0.785, 0.714, 0.000, 0.620, 0.577, 0.000, 0.549, 0.585},
			{1.273, 1.331, 1.224, 1.145, 1.035, 1.037, 0.897, 0.830, 0.833, 0.718, 0.657, 0.653, 0.559, 0.520, 0.488, 0.513, 0.578},
			{1.254, 1.310, 1.200, 1.121, 1.016, 1.017, 0.884, 0.819, 0.821, 0.709, 0.649, 0.643, 0.550, 0.511, 0.480, 0.507, 0.573},
			{1.243, 1.362, 0.000, 1.190, 1.073, 0.000, 0.943, 0.868, 0.000, 0.758, 0.688, 0.000, 0.589, 0.541, 0.000, 0.531, 0.570},
			{1.211, 1.266, 1.164, 1.087, 0.986, 0.991, 0.860, 0.798, 0.802, 0.692, 0.634, 0.630, 0.538, 0.500, 0.470, 0.496, 0.561},
			{1.185, 1.241, 1.142, 1.072, 0.973, 0.975, 0.847, 0.786, 0.789, 0.682, 0.626, 0.622, 0.533, 0.497, 0.466, 0.491, 0.554},
			{1.164, 1.286, 0.000, 1.163, 1.042, 0.000, 0.905, 0.832, 0.000, 0.728, 0.664, 0.000, 0.578, 0.539, 0.000, 0.515, 0.550},
			{1.120, 1.185, 1.134, 1.040, 1.011, 0.968, 0.830, 0.768, 0.772, 0.668, 0.615, 0.619, 0.557, 0.492, 0.476, 0.481, 0.539},
			{1.077, 1.125, 1.080, 0.000, 0.961, 0.982, 0.837, 0.774, 0.776, 0.674, 0.622, 0.631, 0.529, 0.000, 0.465, 0.464, 0.530},
			{1.035, 1.082, 0.985, 0.997, 0.957, 0.000, 0.824, 0.761, 0.000, 0.668, 0.613, 0.000, 0.544, 0.489, 0.429, 0.457, 0.525},
			{1.005, 1.090, 0.984, 0.940, 0.914, 0.922, 0.817, 0.764, 0.759, 0.668, 0.616, 0.606, 0.525, 0.470, 0.439, 0.477, 0.534},
			{1.013, 0.908, 0.854, 0.818, 0.786, 0.756, 0.709, 0.667, 0.631, 0.583, 0.541, 0.503, 0.458, 0.419, 0.397, 0.420, 0.574}
		};

		static const double reference_3[L / 2][L / 2] = {
			{0.794, 0.790, 0.772, 0.751, 0.726, 0.701, 0.659, 0.621, 0.588, 0.545, 0.505, 0.471, 0.429, 0.395, 0.374, 0.390, 0.501},
			{0.000, 0.826, 0.831, 0.826, 0.810, 0.808, 0.739, 0.697, 0.680, 0.612, 0.568, 0.547, 0.486, 0.442, 0.414, 0.424, 0.527},
			{0.000, 0.000, 0.864, 0.900, 0.886, 0.000, 0.801, 0.754, 0.000, 0.665, 0.617, 0.000, 0.537, 0.490, 0.440, 0.439, 0.533},
			{0.000, 0.000, 0.000, 0.000, 0.891, 0.852, 0.768, 0.722, 0.710, 0.637, 0.594, 0.582, 0.540, 0.000, 0.466, 0.444, 0.530},
			{0.000, 0.000, 0.000, 0.000, 0.838, 0.825, 0.748, 0.706, 0.694, 0.623, 0.581, 0.567, 0.513, 0.493, 0.463, 0.442, 0.519},
			{0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.753, 0.711, 0.000, 0.630, 0.586, 0.000, 0.509, 0.474, 0.000, 0.445, 0.506},
			{0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.689, 0.652, 0.643, 0.578, 0.539, 0.523, 0.465, 0.431, 0.426, 0.412, 0.483},
			{0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.619, 0.610, 0.549, 0.512, 0.496, 0.441, 0.409, 0.404, 0.392, 0.460},
			{0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.543, 0.505, 0.000, 0.436, 0.403, 0.000, 0.385, 0.439},
			{0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.489, 0.456, 0.443, 0.393, 0.365, 0.361, 0.350, 0.411},
			{0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.426, 0.413, 0.368, 0.342, 0.337, 0.327, 0.384},
			{0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.360, 0.335, 0.000, 0.315, 0.358},
			{0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.326, 0.313, 0.296, 0.281, 0.328},
			{0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.271, 0.256, 0.300},
			{0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.242, 0.236, 0.277},
			{0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.232, 0.266},
			{0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.286}
		};

		vector_double_2d result(L, std::vector <double> (L));
		vector_double_2d reference(L, std::vector <double> (L));

		for (int i = 0; i < L; i++)
			for (int j = 0; j < L; j++) {
				result[i][j] = 0;
				for (int ii = i * cx; ii < (i + 1) * cx; ii++)
					for (int jj = j * cy; jj < (j + 1) * cy; jj++)
						for (int g = 0; g < G; g++)
							result[i][j] += sigma_f[block[ii][jj]][g] * integral_X[g][id(ii, jj)] * dx[ii] * dy[jj];
				for (int ii = i * cx; ii < (i + 1) * cx; ii++)
					for (int jj = j * cx; jj < (j + 1) * cy; jj++)
						if (block[ii][jj] == 4 || block[ii][jj] == 5) result[i][j] = -1;
			}

		{
			double average = 0;
			for (int i = 0; i < L; i++)
				for (int j = 0; j < L; j++)
					if (result[i][j] >= 0) average += result[i][j];
			average /= total;
			for (int i = 0; i < L; i++)
				for (int j = 0; j < L; j++)
					if (result[i][j] >= 0) result[i][j] /= average;
		}

		for (int i = 0; i < L; i++)
			for (int j = 0; j < L; j++)
				if (result[i][j] < 0) reference[i][j] = -1;
				else {
					if (i < L / 2 && j < L / 2)
						reference[i][j] = (i < j ? reference_1[i][j] : reference_1[j][i]);
					else if (i >= L / 2 && j >= L / 2)
						reference[i][j] = (i < j ? reference_3[i - L / 2][j - L / 2] : reference_3[j - L / 2][i - L / 2]);
					else
						reference[i][j] = (i < j ? reference_2[i][j - L / 2] : reference_2[j][i - L / 2]);
				}

		FILE *file = fopen(file_name, "w");
		for (int i = 0; i < L; i++)
			for (int j = 0; j < L; j++)
				fprintf(file, "%.10lf%c", result[i][j] == -1 ? 0 : result[i][j], j + 1 == L ? '\n' : ' ');
		fclose(file);
		double MAX = 1, MIN = 1, SUM1 = 0, SUM2 = 0, SUM3 = 0, AVG = 0, RMS = 0, MRE = 0;
		for (int i = 0; i < L; i++)
			for (int j = 0; j < L; j++) {
				if (result[i][j] < 0) continue;
				if (result[i][j] > MAX) MAX = result[i][j];
				if (result[i][j] < MIN) MIN = result[i][j];
				if (i < L / 2 && j < L / 2) SUM1 += result[i][j];
				if (i >= L / 2 && j < L / 2) SUM2 += result[i][j];
				if (i >= L / 2 && j >= L / 2) SUM3 += result[i][j];
				double err = std::fabs(result[i][j] - reference[i][j]) / reference[i][j] * 100;
				AVG += err; RMS += err * err; MRE += err * reference[i][j];
			}
		AVG /= total; RMS = sqrt(RMS / total); MRE /= total;

		printf(" ---------- COMPARE WITH STANDARD RESULT (AZTRAN) ---------- \n");
		printf("KEFF      : %-8.6lf (REF = %-9.6lf ERR = %+-8.2lfpcm)\n" , k   , 1.184946, (k    - 1.184946) / 1.184946 * 1e5);
		printf("MAX       : %-8.3lf (REF = %-9.3lf ERR = %+-8.3lf %% )\n", MAX , 2.503   , (MAX  - 2.503   ) / 2.503    * 100);
		printf("MIN       : %-8.3lf (REF = %-9.3lf ERR = %+-8.3lf %% )\n", MIN , 0.232   , (MIN  - 0.232   ) / 0.232    * 100);
		printf("INNER_UO2 : %-8.3lf (REF = %-9.3lf ERR = %+-8.3lf %% )\n", SUM1, 493.703 , (SUM1 - 493.703 ) / 493.703  * 100);
		printf("MOX       : %-8.3lf (REF = %-9.3lf ERR = %+-8.3lf %% )\n", SUM2, 211.339 , (SUM2 - 211.339 ) / 211.339  * 100);
		printf("OUTER_UO2 : %-8.3lf (REF = %-9.3lf ERR = %+-8.3lf %% )\n", SUM3, 139.617 , (SUM3 - 139.617 ) / 139.617  * 100);
		printf("AVG = %.3lf RMS = %.3lf MRE = %.3lf\n", AVG, RMS, MRE);
		printf(" ----------------------------------------------------------- \n");
	}
}

int main() {
	timer.init();
	int cx, cy, n;
	std::ignore = freopen("input.txt", "r", stdin);
	std::ignore = scanf("%d%d%d", &cx, &cy, &n); 
	C5G72D_sn_frame::init(cx, cy, n);
	vector_double_2d integral_X(C5G72D_sn_frame::G, std::vector <double> (C5G72D_sn_frame::Ng, 1));
	double k = 1;
	C5G72D_source_iteration::solve(integral_X, k, 1e-8, 1e4);
	C5G72D_source_iteration::output_power_distribution("C5G72D_power_distribution.txt", integral_X, cx, cy, k);
	return 0;
}

