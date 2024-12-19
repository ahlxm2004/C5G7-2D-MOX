#include <cmath>
#include <cstdio>
#include "basic.hpp"

template <typename Type>
struct linear_equation_solver_C {
	std::vector <std::vector <int> > M, sub_list; 
	std::vector <std::vector <Type> > T1, T2;
	std::vector <Type> T3;
	double error_rms, error_max;
	int n;
	long long nM, nOper;

	void update_error(Type x, Type y) {
		Type tmp = (fabs(y) < 1e-8 ? x : fabs(x / y - 1));
		error_rms += tmp * tmp;
		if (tmp > error_max) error_max = tmp;
	}

	void init_structure() { // init structure from M
		sub_list = std::vector <std::vector <int> > (n, std::vector <int> ());
		std::vector <int> tmp(n);
		std::vector <std::vector <int> > list(n, std::vector <int> ());
		for (int i = 0; i < n; i++)
			for (int j : M[i])
				if (i > j) list[j].push_back(i);
		int sum = 0;
		for (int i = 0; i < n; i++) {
			sum += M[i].size();
			for (int j = 0; j < M[i].size(); j++)
				while (j < M[i].size() && M[i][j] <= i)
					M[i][j] = M[i].back(), M[i].pop_back();
			for (int j : M[i]) tmp[j] = 0;
			for (int j : list[i]) {
				sub_list[j].push_back(i);
				for (int k : M[j]) tmp[k] = j;
				for (int k : M[i])
					if (tmp[k] != j) {
						M[j].push_back(k);
						if (j > k) list[k].push_back(j);
					}
			}
		}
		nM = nOper = 0;
		for (int i = 0; i < n; i++) nM += M[i].size();
		for (int i = 0; i < n; i++) nOper += sub_list[i].size() + M[i].size() + 1;
	}

	void init_value(sparse_matrix G) { // init value from G
		std::vector <std::vector <std::pair <int, Type> > > CSR(n);
		T1.resize(n), T2.resize(n), T3.resize(n);
		std::vector <Type> V(n);
		Type coef;
		for (sparse_matrix_element i : G)
			CSR[i.x].emplace_back(i.y, i.key);
		for (int i = 0; i < n; i++) {
			for (std::pair <int, Type> j : CSR[i]) V[j.first] = j.second;
			T1[i].resize(sub_list[i].size());
			for (int _ = 0; _ < sub_list[i].size(); _++) {
				int j = sub_list[i][_];
				coef = -V[j];
				for (int k = 0; k < M[j].size(); k++)
					V[M[j][k]] += coef * T2[j][k];
				V[j] = 0; T1[i][_] = coef;
			}
			coef = 1.0 / V[i]; T3[i] = coef;
			V[i] = 0; T2[i].resize(M[i].size());
			for (int j = 0; j < M[i].size(); j++)
				T2[i][j] = V[M[i][j]] * coef, V[M[i][j]] = 0;
		}
	}

	void init_structure(int _n, sparse_structure _M) {
		n = _n, M = _M, init_structure();
	}

	void init(int _n, sparse_matrix G) {
		n = _n;
		M.clear(); M.resize(n);
		for (auto i : G) M[i.x].push_back(i.y);
		init_structure();
		init_value();
	}

	std::vector <double> solve(std::vector <double> b) {
		std::vector <double> X = b, _b(n, 0);
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < sub_list[i].size(); j++)
				X[i] += X[sub_list[i][j]] * T1[i][j];
			X[i] *= T3[i];
		}
		for (int i = n - 1; i >= 0; i--)
			for (int j = 0; j < M[i].size(); j++)
				X[i] -= X[M[i][j]] * T2[i][j];
		return X;
	}
};

template <typename Type>
struct linear_equation_solver_Ck {
	int K;
	std::vector <std::vector <int> > M, sub_list; 
	std::vector <int> mL, mU;
	std::vector <int> L_id, U_id;
	std::vector <std::vector <Type> > _T1, _T2;
	std::vector <std::vector <Type> > _T3;
	int n;
	int nA;

	void init_structure() { // init structure from M
		for (int i = 0; i < n; i++) nA += M[i].size();
		std::vector <std::vector <int> > sub_list(n, std::vector <int> ());
		std::vector <int> tmp(n);
		std::vector <std::vector <int> > list(n, std::vector <int> ());
		for (int i = 0; i < n; i++)
			for (int j : M[i])
				if (i > j) list[j].push_back(i);
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < M[i].size(); j++)
				while (j < M[i].size() && M[i][j] <= i)
					M[i][j] = M[i].back(), M[i].pop_back();
			for (int j : M[i]) tmp[j] = 0;
			for (int j : list[i]) {
				sub_list[j].push_back(i);
				for (int k : M[j]) tmp[k] = j;
				for (int k : M[i])
					if (tmp[k] != j) {
						M[j].push_back(k);
						if (j > k) list[k].push_back(j);
					}
			}
		}
		L_id.resize(n + 1); L_id[0] = 0;
		U_id.resize(n + 1); U_id[0] = 0;
		for (int i = 1; i <= n; i++)
			L_id[i] = L_id[i - 1] + sub_list[i - 1].size(),
			U_id[i] = U_id[i - 1] + M[i - 1].size();
		mL.resize(L_id[n]); mU.resize(U_id[n]);
		for (int i = 0; i < n; i++)
			for (int j = 0; j < sub_list[i].size(); j++)
				mL[L_id[i] + j] = sub_list[i][j];
		for (int i = 0; i < n; i++)
			for (int j = 0; j < M[i].size(); j++)
				mU[U_id[i] + j] = M[i][j];
		std::vector <std::vector <int> > ().swap(sub_list);
		std::vector <std::vector <int> > ().swap(M);
		printf("A (%d * %d, nnz = %d) = L (nnz = %d) * U (nnz = %d)\n", n, n, nA, L_id[n] + n, U_id[n] + n);
	}

	void init_value(sparse_matrix G, int k) { // init value from G
		std::vector <std::vector <std::pair <int, Type> > > CSR(n);
		std::vector <Type> &T1 = _T1[k], &T2 = _T2[k], &T3 = _T3[k];
		T1.resize(L_id[n]), T2.resize(U_id[n]), T3.resize(n);
		std::vector <Type> V(n, 0);
		Type coef;
		for (sparse_matrix_element i : G)
			CSR[i.x].emplace_back(i.y, i.key);
		for (int i = 0; i < n; i++) {
			for (std::pair <int, Type> j : CSR[i]) V[j.first] = j.second;
			for (int _ = L_id[i]; _ < L_id[i + 1]; _++) {
				int j = mL[_]; coef = -V[j];
				for (int k = U_id[j]; k < U_id[j + 1]; k++) V[mU[k]] += coef * T2[k];
				V[j] = 0; T1[_] = coef;
			}
			coef = 1.0 / V[i]; T3[i] = coef; V[i] = 0;
			for (int j = U_id[i]; j < U_id[i + 1]; j++)
				T2[j] = V[mU[j]] * coef, V[mU[j]] = 0;
		}
	}

	void init_K(int _K) {K = _K; _T1.resize(K); _T2.resize(K); _T3.resize(K);}

	void init_structure(int _n, sparse_structure _M) {
		n = _n, M = _M, init_structure();
	}

	std::vector <double> solve(std::vector <double> b, int k) {
		std::vector <Type> &T1 = _T1[k], &T2 = _T2[k], &T3 = _T3[k];
		std::vector <double> X = b, _b(n, 0);
		for (int i = 0, l = 0; i < n; i++) {
			while (l < L_id[i + 1]) X[i] += X[mL[l]] * T1[l], l++;
			X[i] *= T3[i];
		}
		for (int i = n - 1, l = U_id[n]; i >= 0; i--)
			while (l > U_id[i]) l--, X[i] -= X[mU[l]] * T2[l];
		return X;
	}
};

struct linear_equation_solver_D {
	sparse_structure M;
	std::vector <std::vector <std::pair <int, double> > > sub_list;
	std::vector <std::vector <std::pair <int, double> > > T1;
	std::vector <std::vector <double> > T2;
	std::vector <double> T3;
	std::vector <int> CUT;
	double error_rms, error_max;
	std::vector <bool> over;
	std::vector <int> bel;
	int n;
	long long nG, nM, nOper;

	void update_error(double x, double y) {
		double tmp = (fabs(y) < 1e-8 ? x : fabs(x / y - 1));
		error_rms += tmp * tmp;
		if (tmp > error_max) error_max = tmp;
	}

	void init(int _n, sparse_matrix G, std::vector <int> _CUT) {
		printf("init!\n");
		double coef;
		n = _n; CUT = _CUT;
		std::vector <int> tmp(n);
		std::vector <double> V(n);
		T1.clear(), T2.clear(), T3.clear(), sub_list.clear(), M.clear();
		T1.resize(n), T2.resize(n), T3.resize(n), sub_list.resize(n), M.resize(n);
		std::vector <std::vector <int> > list(n, std::vector <int> ());
		bel.resize(n), over.resize(n);
		for (int cut = 0; cut + 1 < CUT.size(); cut++)
			for (int i = CUT[cut]; i < CUT[cut + 1]; i++) bel[i] = cut;
		for (int cut = 0; cut + 1 < CUT.size(); cut++) {
			int l = CUT[cut], r = CUT[cut + 1];
			printf("[%d %d) : ", l, r);
			for (sparse_matrix_element i : G)
				if (bel[i.x] == cut && bel[i.y] == cut)
					list[i.y].push_back(i.x), M[i.x].push_back(i.y);
			for (int i = l; i < r; i++) over[i] = false;
			for (int i = l; i < r; i++) {
				if (i == l || i + 1 == r || i % 20000 == 0) printf("%d ", i);
				over[i] = true;
				for (int j = 0; j < M[i].size(); j++) {
					while (j < M[i].size() && over[M[i][j]])
						M[i][j] = M[i].back(), M[i].pop_back();
				}
				
				for (int j : M[i]) tmp[j] = 0;
				for (int j : list[i]) {
					if (over[j]) continue;
					sub_list[j].push_back(std::make_pair(i, 0.));
					for (int k : M[j]) tmp[k] = j;
					for (int k : M[i])
						if (tmp[k] != j) M[j].push_back(k), list[k].push_back(j);
				}
			}
			printf("|| ");
			for (int i = l; i < r; i++) T1[i].clear(), over[i] = false;
			for (sparse_matrix_element i : G)
				if (bel[i.x] == cut && bel[i.y] == cut) T1[i.x].emplace_back(i.y, i.key);
			for (int i = l; i < r; i++) {
				if (i == l || i + 1 == r || i % 20000 == 0) printf("%d ", i);
				for (std::pair <int, double> j : T1[i]) V[j.first] = j.second;
				std::vector <std::pair <int, double> > ().swap(T1[i]);
				for (auto &j : sub_list[i]) {
					coef = -V[j.first];
					for (int k = 0; k < M[j.first].size(); k++)
						V[M[j.first][k]] += coef * T2[j.first][k];
					V[j.first] = 0; j.second = coef;
				}
				coef = 1.0 / V[i]; T3[i] = coef;
				V[i] = 0; T2[i].resize(M[i].size());
				for (int j = 0; j < M[i].size(); j++)
					T2[i][j] = V[M[i][j]] * coef, V[M[i][j]] = 0;
			}
			printf("\n");
		}
		for (sparse_matrix_element i : G)
			if (bel[i.x] > bel[i.y]) sub_list[i.x].emplace_back(i.y, -i.key);
		nG = G.size(), nM = nOper = 0;
		for (int i = 0; i < n; i++) nM += M[i].size();
		for (int i = 0; i < n; i++) nOper += sub_list[i].size() + M[i].size() + 1;
	}

	std::vector <double> solve(std::vector <double> b) {
		std::vector <double> X = b, _b(n, 0);
		for (int cut = 0; cut + 1 < CUT.size(); cut++) {
			int l = CUT[cut], r = CUT[cut + 1];
			for (int i = l; i < r; i++) {
				for (auto j : sub_list[i])
					X[i] += X[j.first] * j.second;
				X[i] *= T3[i];
			}
			for (int i = r - 1; i >= l; i--)
				for (int j = 0; j < M[i].size(); j++)
					X[i] -= X[M[i][j]] * T2[i][j];
		}
		return X;
	}
};

