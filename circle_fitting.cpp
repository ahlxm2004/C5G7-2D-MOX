/*
reference:
	A. Pautz, S. Langenbuch, A. Seubert, and W. Zwermann, ¡°Results on the OECD/NEA
	C5G7-MOX benchmark obtained with the discrete ordinates DORT code,¡± Progress in
	Nuclear Energy, vol. 45, no. 2, pp. 153¨C168, 2004.
*/

#include <bits/stdc++.h>

const double C[10][10] = {
	{0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000},
	{0.886227, 0.536748, 0.400733, 0.321452, 0.268550, 0.230491, 0.201761, 0.179313},
	{0.000000, 1.000000, 0.780889, 0.628202, 0.526487, 0.453302, 0.397889, 0.354413},
	{0.000000, 0.000000, 1.000000, 0.875520, 0.753324, 0.656709, 0.580791, 0.520017},
	{0.000000, 0.000000, 0.000000, 1.000000, 0.918891, 0.825076, 0.741269, 0.670152},
	{0.000000, 0.000000, 0.000000, 0.000000, 1.000000, 0.942572, 0.869406, 0.798369},
	{0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 1.000000, 0.957052, 0.898665},
	{0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 1.000000, 0.966606},
	{0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 1.000000}
};

double t1 = 0.54, t2 = 0.63;

int n;
double k;

std::vector <double> tmp;

void push_back(double x) {
	if (x < 2 * k) tmp.push_back(x);
	else {
		int s = (int)(x / k);
		for (int i = 0; i < s; i++) tmp.push_back(x / s);
	}
}

int main() {
	scanf("%d%lf", &n, &k);
	push_back(C[1][n - 1] * t1);
	for (int m = 1; m < n; m++)
		push_back((C[m + 1][n - 1] - C[m][n - 1]) * t1);
	push_back(t2 - t1 + (1 - C[n][n - 1]) * t1);
	printf("N = %d\n", (int)tmp.size() * 2);
	for (int i = (int)tmp.size() - 1; i >= 0; i--)
		printf("%.8lf\n", tmp[i]);
	for (int i = 0; i < (int)tmp.size(); i++)
		printf("%.8lf\n", tmp[i]);
	return 0;
}

