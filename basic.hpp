#ifndef _BASIC_H_
#define _BASIC_H_

#include <omp.h>
#include <vector>
#include <cmath>
#include <ctime>

template <typename Type>
struct row_element_basic {
	int col; Type key;
	row_element_basic(int _col = 0, Type _key = 0) : col(_col), key(_key) {}
	bool operator < (const row_element_basic &b) const {
		return col < b.col;
	}
};

typedef row_element_basic <double> row_element;
typedef row_element_basic <float> row_element_f;

struct sparse_matrix_element {
	int x, y; double key;
	sparse_matrix_element(int _x = 0, int _y = 0, double _key = 0) : x(_x), y(_y), key(_key) {}
	bool operator < (const sparse_matrix_element &b) const {
		return x < b.x || (x == b.x && y < b.y);
	}
};

typedef std::vector <sparse_matrix_element> sparse_matrix;

void reduce_sparse_matrix(sparse_matrix &A) {
	for (int i = 0; i < A.size(); i++)
		if (fabs(A[i].key) < 1e-16) std::swap(A[i], A.back()), A.pop_back(), i--;
}

typedef std::vector <std::vector <double> > full_matrix;

template <typename Type>
struct compressed_matrix_basic {
	int n;
	std::vector <Type> C, K, inv_K;
	std::vector <unsigned> V, T;

	compressed_matrix_basic() {n = 0; C = K = inv_K = {}; V = T = {};}

	void init(int _n, std::vector <std::vector <row_element_basic <Type> > > &A, int nnz) {
		n = _n;
		C.resize(nnz), V.resize(nnz), T.resize(n + 1), K.resize(n); inv_K = {};
		for (int i = 0, sum = 0; i < n; i++) {
			T[i] = sum;
			for (row_element_basic <Type> j : A[i])
				if (i == j.col) K[i] = j.key;
				else V[sum] = j.col, C[sum] = j.key, sum++;
		}
		T[n] = nnz;
	}

	compressed_matrix_basic(int _n, std::vector <std::vector <row_element_basic <Type> > > &A, int nnz = 0) {
		if (nnz == 0) {
			#pragma omp parallel for reduction(+: nnz)
			for (int i = 0; i < _n; i++) nnz += A[i].size();
		}
		init(_n, A, nnz);
	}

	compressed_matrix_basic(int _n, sparse_matrix &A) {
		std::vector <std::vector <row_element_basic <Type> > > B(_n, std::vector <row_element_basic <Type> > ());
		for (sparse_matrix_element i : A) B[i.x].emplace_back(i.y, i.key);
		init(_n, B, A.size());
	}

	void calc_inv_K() {
		if (inv_K.empty()) {
			inv_K.resize(n);
			for (int i = 0; i < n; i++) inv_K[i] = 1. / K[i];
		}
	}

	void fma(std::vector <double> &x, std::vector <double> &y) {
		#pragma omp parallel
		{
			int k = omp_get_thread_num();
			int all = omp_get_num_threads();
			int l = n * k / all, r = n * (k + 1) / all;
			std::vector <unsigned> ::iterator it1 = V.begin() + T[l];
			typename std::vector <Type> ::iterator it2 = C.begin() + T[l];
			for (int i = l, j = T[l]; i != r; i++) {
				x[i] += K[i] * y[i];
				while (j != T[i + 1]) x[i] += (*it2) * y[*it1], it1++, it2++, j++;
			}
		}
	}

	void fms(std::vector <double> &x, std::vector <double> &y) {
		#pragma omp parallel
		{
			int k = omp_get_thread_num();
			int all = omp_get_num_threads();
			int l = n * k / all, r = n * (k + 1) / all;
			std::vector <unsigned> ::iterator it1 = V.begin() + T[l];
			typename std::vector <Type> ::iterator it2 = C.begin() + T[l];
			for (int i = l, j = T[l]; i != r; i++) {
				x[i] -= K[i] * y[i];
				while (j != T[i + 1]) x[i] -= (*it2) * y[*it1], it1++, it2++, j++;
			}
		}
	}

	std::vector <double> multiply(std::vector <double> &x, bool neg = false) {
		std::vector <double> y = std::vector <double> (n, 0);
		return (neg ? fms(y, x) : fma(y, x)), y;
	}
};

typedef compressed_matrix_basic <double> compressed_matrix;
typedef compressed_matrix_basic <float> compressed_matrix_f;

typedef std::vector <std::vector <int> > sparse_structure;

#ifdef _OPENMP

int parallel_num_threads() {
	int result;
	#pragma omp parallel
		result = omp_get_num_threads();
	return result;
}

struct Timer {
	double tl;
	void init() {tl = omp_get_wtime();}
	double time() {return omp_get_wtime() - tl;}
};

#else

int parallel_num_threads() {return 1;}

struct Timer {
	double tl;
	void init() {tl = clock();}
	double time() {return (clock() - tl) / CLOCKS_PER_SEC;}
};

#endif

typedef std::vector <std::vector <int> > vector_int_2d;
typedef std::vector <std::vector <double> > vector_double_2d;
typedef std::vector <vector_double_2d> vector_double_3d;
typedef std::vector <vector_int_2d> vector_int_3d;
typedef std::vector <vector_double_3d> vector_double_4d;
typedef std::vector <vector_int_3d> vector_int_4d;

#endif

