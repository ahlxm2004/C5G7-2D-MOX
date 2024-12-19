struct Order {
	int n;
	std::vector <int> order, iorder;

	Order(std::vector <int> _order = {}) {
		order = _order, n = _order.size(), iorder.resize(n);
		#pragma omp parallel for
		for (int i = 0; i < n; i++) iorder[order[i]] = i;
	}

	std::vector <double> X(std::vector <double> A) {
		std::vector <double> B(A.size());
		#pragma omp parallel for
		for (int i = 0; i < A.size(); i++) B[i] = A[order[i]];
		return B;
	}

	sparse_structure X(sparse_structure A) {
		sparse_structure B(A.size(), std::vector <int> ());
		for (auto i = 0; i < A.size(); i++)
			for (int j : A[i])
				B[iorder[j]].push_back(iorder[i]);
		return B;
	}

	sparse_matrix X(sparse_matrix A) {
		#pragma omp parallel for
		for (int i = 0; i < A.size(); i++)
			A[i].x = iorder[A[i].x], A[i].y = iorder[A[i].y];
		return A;
	}

	std::vector <double> Y(std::vector <double> A) {
		std::vector <double> B(A.size());
		#pragma omp parallel for
		for (int i = 0; i < A.size(); i++) B[order[i]] = A[i];
		return B;
	}

	sparse_matrix Y(sparse_matrix A) {
		#pragma omp parallel for
		for (int i = 0; i < A.size(); i++)
			A[i].x = order[A[i].x], A[i].y = order[A[i].y];
		return A;
	}

	Order duplicate(int k) {
		std::vector <int> result(k * n);
		for (int i = 0; i < n; i++)
			for (int j = 0; j < k; j++)
				result[i * k + j] = j * n + order[i];
		return Order(result); 
	}
};


namespace order_2D_square {
	int L1, L2;

	std::vector <int> ans;

	int id(int x, int y) {return L2 * x + y;}

	void divide(int l1, int r1, int l2, int r2) {
		if (r1 - l1 < 3 && r2 - l2 < 3) {
			for (int i = l1; i <= r1; i++)
				for (int j = l2; j <= r2; j++)
					ans.push_back(id(i, j));
		}
		else if (r1 - l1 < r2 - l2) {
			int mid2 = (l2 + r2) / 2;
			divide(l1, r1, l2, mid2 - 1);
			divide(l1, r1, mid2 + 1, r2);
			divide(l1, r1, mid2, mid2);
		}
		else {
			int mid1 = (l1 + r1) / 2;
			divide(l1, mid1 - 1, l2, r2);
			divide(mid1 + 1, r1, l2, r2);
			divide(mid1, mid1, l2, r2);
		}
	}

	std::vector <int> get_order(int _L1, int _L2) {
		L1 = _L1, L2 = _L2;
		ans.clear(); ans.reserve(L1 * L2);
		/*for (int i = 0; i < L1; i++)
			for (int j = 0; j < i; j++)
				ans.push_back(id(i, j));*/
		divide(0, L1 - 1, 0, L2 - 1);
		return ans;
	}
}

namespace order_3D_square {
	int L1, L2, L3;

	std::vector <int> ans;

	int id(int x, int y, int z) {return (x * L2 + y) * L3 + z;}

	void divide(int l1, int r1, int l2, int r2, int l3, int r3) {
		if (r1 - l1 < 2 && r2 - l2 < 2 && r3 - l3 < 2) {
			for (int i = l1; i <= r1; i++)
				for (int j = l2; j <= r2; j++)
					for (int k = l3; k <= r3; k++)
						if (i <= j) ans.push_back(id(i, j, k));
		}
		else if (r1 - l1 > r2 - l2 && r1 - l1 > r3 - l3) {
			int mid1 = (l1 + r1) / 2;
			divide(l1, mid1 - 1, l2, r2, l3, r3);
			divide(mid1 + 1, r1, l2, r2, l3, r3);
			divide(mid1, mid1, l2, r2, l3, r3);
		}
		else if (r2 - l2 >= r1 - l1 && r2 - l2 > r3 - l3) {
			int mid2 = (l2 + r2) / 2;
			divide(l1, r1, l2, mid2 - 1, l3, r3);
			divide(l1, r1, mid2 + 1, r2, l3, r3);
			divide(l1, r1, mid2, mid2, l3, r3);
		}
		else {
			int mid3 = (l3 + r3) / 2;
			divide(l1, r1, l2, r2, l3, mid3 - 1);
			divide(l1, r1, l2, r2, mid3 + 1, r3);
			divide(l1, r1, l2, r2, mid3, mid3);
		}
	}

	std::vector <int> get_order(int _L1, int _L2, int _L3) {
		L1 = _L1, L2 = _L2, L3 = _L3;
		ans.clear(); ans.reserve(L1 * L2 * L3);
		for (int i = 0; i < L1; i++)
			for (int j = 0; j < i; j++)
				for (int k = 0; k < L3; k++) ans.push_back(id(i, j, k));
		divide(0, L1 - 1, 0, L2 - 1, 0, L3 - 1);
		return ans;
	}
}

