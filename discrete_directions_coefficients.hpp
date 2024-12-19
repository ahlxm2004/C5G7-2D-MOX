/*
reference :
	discrete ordinates angular quadrature of the neutron transport equation,
	LOS ALAMOS SCIENTIFIC LABORATORY OF THE UNIVERSITY OF CALIFORNIA
	website : https://www.osti.gov/servlets/purl/4666281
	DOI : https://doi.org/10.2172/4666281
*/ 

#ifndef _DISCRETE_DIRECTIONS_COEFFICIENTS_H_
#define _DISCRETE_DIRECTIONS_COEFFICIENTS_H_

#include <vector>
#include <algorithm>
#include <cassert>
#include <numeric>

std::vector <int> discrete_directions_N = {2, 4, 6, 8, 12, 16};

vector_double_2d discrete_directions_mu = {
	{0.5773503},
	{0.3500212, 0.8688903},
	{0.2666355, 0.6815076, 0.9261808},
	{0.2182179, 0.5773503, 0.7867958, 0.9511897},
	{0.1672126, 0.4595476, 0.6280191, 0.7600210, 0.8722706, 0.9716377},
	{0.1389568, 0.3922893, 0.5370966, 0.6504264, 0.7467506, 0.8319966, 0.9092855, 0.9805009}
};

vector_double_2d discrete_directions_weight = {
	{1.0000000},
	{0.3333333},
	{0.1761263, 0.1572071},
	{0.1209877, 0.0907407, 0.0925926},
	{0.0707626, 0.0558811, 0.0373377, 0.0502819, 0.0258513},
	{0.0489872, 0.0413296, 0.0212326, 0.0256207, 0.0360486, 0.0144589, 0.0344958, 0.0085179}
};

vector_double_2d get_discrete_directions_coefficients(int n) {
	int _ = std::find(discrete_directions_N.begin(),
		discrete_directions_N.end(), n) - discrete_directions_N.begin();
	assert(_ != discrete_directions_N.size());
	vector_int_2d index(n);
	for (int i = 0; i < n; i++) index[i].resize(i + 1);
	for (int i = 0, j = n / 2 - 1, u = 0; i <= j; i += 2, j--)
		for (int offset = 0; i + offset <= j - offset; offset++, u++)
			index[i + offset][i / 2] = u,
			index[i + offset][i / 2 + offset] = u,
			index[j - offset][i / 2] = u,
			index[j - offset][j - offset - i / 2] = u,
			index[j][i / 2 + offset] = u,
			index[j][j - offset - i / 2] = u;
	vector_double_2d result(4, std::vector <double> ());
	for (int i = 0; i < n / 2; i++)
		for (int j = 0; j <= i; j++)
			result[0].push_back(discrete_directions_mu[_][n / 2 - i - 1]),
			result[1].push_back(discrete_directions_mu[_][j]),
			result[2].push_back(discrete_directions_mu[_][i - j]),
			result[3].push_back(discrete_directions_weight[_][index[i][j]]);
	double sum = std::accumulate(result[3].begin(), result[3].end(), 0.);
	for (double &i : result[3]) i /= sum;
	return result;
}

#endif

