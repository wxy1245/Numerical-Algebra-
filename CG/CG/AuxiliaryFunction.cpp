#include"auxiliary.h"
#include <cmath>
//#include <random>

vector<double> matrix_vector(vector<vector<double>> A, vector<double> x)
{
	int m = A.size();	//m取A的行数
	int n = A[0].size();//n取A的列数
	vector<double> result(m, 0);
	for (int i = 0; i < m; i++)
	{
		for (int j = 0; j < n; j++)
		{
			result[i] += A[i][j] * x[j];
		}
	}
	return result;
}

double inner_product(vector<double> x, vector<double> y)
{
	int n = x.size();
	double result = 0;
	for (int i = 0; i < n; i++)
	{
		result += x[i] * y[i];
	}
	return result;
}

vector<double> vector_plus(vector<double> x, vector<double> y)
{
	int n = x.size();
	vector<double> result(n);
	for (int i = 0; i < n; i++)
	{
		result[i] = x[i] + y[i];
	}
	return result;
}

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

vector<double> scale_vector(double alpha, vector<double> x)
{
	int n = x.size();
	for (int i = 0; i < n; i++)
	{
		x[i] *= alpha;
	}
	return x;
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

double vector_norm_infinity(vector<double> z)
{
	int n = z.size();
	double maximum = 0;
	for (int i = 1; i < n; i++)
	{
		if (abs(z[i]) > maximum)
		{
			maximum = abs(z[i]);
		}
	}
	return maximum;	//返回的是取得最大模的下标
}