#include <cstdio>
#include "basic.hpp"

template <typename Type>
struct linear_equation_solver_k {
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
