#include <vector>
#include <algorithm>
#include <set>
#include <numeric>
#include <ctime>

struct jacobian_query {
	bool is_row;
	std::vector <int> list;
	std::vector <std::pair <int, int> > grid;
};

namespace jacobian_coloring {
	int N1, N2, E;
	int c1, c2;
	double t1, t2;

	std::vector <std::vector <int> > link_x, link_y;

	void init(int _N1, int _N2, sparse_structure G, bool binary = true) {
		N1 = _N1, N2 = _N2, E = 0; 
		link_x = std::vector <std::vector <int> > (N2, std::vector <int> ());
		for (int i = 0; i < N1; i++) E += G[i].size();
		for (int i = 0; i < N1; i++)
			for (int j : G[i]) link_x[j].push_back(i);
		if (binary) link_y = G;
	}

	std::vector <jacobian_query> solve1() {
		std::vector <int> V(N2), val2(N2);
		std::vector <bool> colored(N2), cover(N1);
		for (int i = 0; i < N2; i++)
			colored[i] = link_x[i].empty();
		std::vector <jacobian_query> result;
		for (int i = 0; i < N2; i++) val2[i] = link_x[i].size(); 
		std::iota(V.begin(), V.end(), 0);
		std::sort(V.begin(), V.end(), [&] (int u, int v) {
			return val2[u] > val2[v];
		});
		for (int num = std::count(colored.begin(), colored.end(), true); num < N2; ) {
			jacobian_query current; current.is_row = false;
			std::fill(cover.begin(), cover.end(), false);
			for (int i = 0; i < N2; i++) {
				if (colored[V[i]]) continue;
				for (int j : link_x[V[i]])
					if (cover[j]) goto stop;
				colored[V[i]] = true, num++;
				current.list.push_back(V[i]);
				for (int j : link_x[V[i]])
					current.grid.emplace_back(j, V[i]);
				for (int j : link_x[V[i]]) cover[j] = true;
				stop : ;
			}
			result.push_back(current);
		}
		return result;
	}

	std::vector <jacobian_query> solve2() {
		std::vector <bool> colored_1(N1), colored_2(N2), cover_1(N2), cover_2(N1);
		std::vector <int> V1(N1), V2(N2), val1(N1), val2(N2);
		for (int i = 0; i < N1; i++)
			for (int j : link_x[i]) val1[i] += link_y[j].size();
		for (int i = 0; i < N2; i++)
			for (int j : link_y[i]) val2[i] += link_x[j].size();
		std::iota(V1.begin(), V1.end(), 0);
		std::iota(V2.begin(), V2.end(), 0);
		std::sort(V1.begin(), V1.end(), [&] (int u, int v) {
			return val1[u] > val1[v];
		});
		std::sort(V2.begin(), V2.end(), [&] (int u, int v) {
			return val2[u] > val2[v];
		});
		for (int i = 0; i < N1; i++) colored_1[i] = link_x[i].empty();
		for (int i = 0; i < N2; i++) colored_2[i] = link_y[i].empty();
		std::vector <jacobian_query> result;
		for (int p1 = 0, p2 = 0, _ = 0; _ < E; ) {
			while (p1 < N1 && colored_1[V1[p1]]) p1++;
			while (p2 < N2 && colored_2[V2[p2]]) p2++; 
			jacobian_query current;
			if (p1 < N1 && (p2 == N2 || val1[V1[p1]] > val2[V2[p2]])) {
				current.is_row = true;
				std::fill(cover_1.begin(), cover_1.end(), false);
				for (int i = p1; i < N1; i++) {
					if (colored_1[V1[i]]) continue;
					for (int j : link_y[V1[i]])
						if (cover_1[j] && !colored_2[j]) goto stop1;
					colored_1[V1[i]] = true;
					current.list.push_back(V1[i]);
					for (int j : link_y[V1[i]])
						if (!colored_2[j]) current.grid.emplace_back(V1[i], j), _++;
					for (int j : link_y[V1[i]]) cover_1[j] = true;
					stop1 : ;
				}
				result.push_back(current);
			}
			else {
				current.is_row = false;
				std::fill(cover_2.begin(), cover_2.end(), false);
				for (int i = p2; i < N2; i++) {
					if (colored_2[V2[i]]) continue;
					for (int j : link_x[V2[i]])
						if (cover_2[j] && !colored_1[j]) goto stop2;
					colored_2[V2[i]] = true;
					current.list.push_back(V2[i]);
					for (int j : link_x[V2[i]])
						if (!colored_1[j]) current.grid.emplace_back(j, V2[i]), _++;
					for (int j : link_x[V2[i]]) cover_2[j] = true;
					stop2 : ;
				}
				result.push_back(current);
			}
		}
		return result;
	}

	std::vector <jacobian_query> solve(bool binary = true) {
		clock_t t = clock();
		std::vector <jacobian_query> X = solve1(); c1 = X.size();
		t1 = (double)(clock() - t) / CLOCKS_PER_SEC; t = clock();
		if (binary) {
			std::vector <jacobian_query> Y = solve2(); c2 = Y.size();
			t2 = (double)(clock() - t) / CLOCKS_PER_SEC;
			if (X.size() > Y.size()) X.swap(Y);
		}
		std::vector <std::vector <int> > ().swap(link_x);
		std::vector <std::vector <int> > ().swap(link_y);
		return X;
	}

	std::vector <jacobian_query> solve(int _N1, int _N2, std::vector <std::vector <int> > G, bool binary = true) {
		return init(_N1, _N2, G, binary), solve(binary);
	}
}

