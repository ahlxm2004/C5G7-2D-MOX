#include <vector>

template <typename T>
struct anderson {
	const int D = 10;

	int n, m, max_m, ne, error_width;
	std::vector <double> F, G;
	std::vector <T> pre;
	std::vector <std::vector <T> > Xs, Gk;
	std::vector <std::vector <double> > R;
	std::vector <double> error_list; 

	void update(std::vector <double> X, std::vector <double> f_X) {
		if (m == max_m) {
			if ((--m) == -1) pre.clear();
			else {
				for (int i = 0; i < m; i++)
					Xs[i].swap(Xs[i + 1]), Gk[i].swap(Gk[i + 1]);
			}
			for (int i = 0; i < m; i++)
				for (int j = 0; j < m; j++) R[i][j] = R[i + 1][j + 1];
		}
		if ((int)Gk.size() <= m)
			Gk.push_back(std::vector <T> (ne)),
			Xs.push_back(std::vector <T> (n));
		if (m != -1) {
			#pragma omp parallel for
			for (int i = 0; i < n; i++) Xs.back()[i] = f_X[i] - F[i];
			#pragma omp parallel for
			for (int i = 0; i < ne; i++) Gk.back()[i] = Xs.back()[i * D] - X[i * D] + pre[i];
			for (int i = 0; i <= m; i++) {
				double sum = 0;
				#pragma omp parallel for reduction(+: sum)
				for (int j = 0; j < ne; j++)
					sum += Gk[i][j] * Gk.back()[j];
				R[i][m] = R[m][i] = sum * n / ne;
			}
		}
		G.resize(ne); pre.resize(ne);
		#pragma omp parallel for
		for (int i = 0; i < n; i += D) G[i / D] = f_X[i] - X[i], pre[i / D] = X[i];
		m++; F = f_X;

		error_list.push_back(0);
		for (int i = 0; i < error_width; i++)
			error_list.back() += (f_X[i] - X[i]) * (f_X[i] - X[i]);
		error_list.back() = sqrt(error_list.back() / error_width);
	}

	void init(int _n, int _max_m, std::vector <double> X, std::vector <double> f_X) {
		n = _n, max_m = _max_m; ne = (n - 1) / D + 1;
		if (max_m > n) max_m = n;
		m = -1, Xs.clear(), Gk.clear(), pre.clear(), error_list.clear();
		error_width = n; update(X, f_X);
		R = std::vector <std::vector <double> > (max_m, std::vector <double> (max_m, 0));
	}

	void set_error_width(int _error_width) {
		error_width = _error_width;
	}

	void update_m(int d) {
		max_m += d; R.resize(max_m);
		for (int i = 0; i < max_m; i++) R[i].resize(max_m);
	}

	std::vector <double> Gauss_elimination(std::vector <double> b) {
		std::vector <std::vector <double> > R0 = R;
		for (int i = 0; i < m; i++) {
			int t = i;
			for (int j = i + 1; j < m; j++)
				if (fabs(R0[j][i]) > fabs(R0[t][i])) t = j;
			for (int j = i; j < m; j++) std::swap(R0[i][j], R0[t][j]);
			std::swap(b[i], b[t]);
			double u = 1.0 / R0[i][i];
			for (int j = i + 1; j < m; j++) R0[i][j] *= u;
			R0[i][i] = 1; b[i] *= u;
			for (int j = i + 1; j < m; j++) {
				u = -R0[j][i]; b[j] += u * b[i];
				for (int k = i; k < m; k++) R0[j][k] += u * R0[i][k];
			}
		}
		for (int i = m - 1; i >= 0; i--)
			for (int j = i + 1; j < m; j++) b[i] -= b[j] * R0[i][j];
		return b;
	}

	std::vector <double> calc_next() {
		if (Gk.empty()) return F;
		std::vector <double> next = F;
		std::vector <double> gamma(m);
		for (int i = 0; i < m; i++) {
			double sum = 0;
			#pragma omp parallel for reduction(+: sum)
			for (int j = 0; j < ne; j++)
				sum += Gk[i][j] * G[j];
			gamma[i] = sum * n / ne;
		}
		gamma = Gauss_elimination(gamma);
		for (int i = 0; i < m; i++)
			#pragma omp parallel for
			for (int j = 0; j < n; j++)
				next[j] -= gamma[i] * Xs[i][j];
		return next;
	}
};
