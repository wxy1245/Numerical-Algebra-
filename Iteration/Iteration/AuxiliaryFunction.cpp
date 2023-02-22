#include"auxiliary.h"
#include <cmath>
//#include <random>

vector<double> vector_minus(vector<double> x, vector<double> y)
{
	int n = x.size();
	vector<double> result(n);
	for (int i = 0; i < n; i++)
	{
		result[i] = x[i] - y[i];
	}
	return result;
}

double vector_norm_2(vector<double> z)
{

	int n = z.size();
	double sum = 0;
	for (int i = 0; i < n; i++)
	{
		sum += z[i] * z[i];
	}
	return sqrt(sum);
}