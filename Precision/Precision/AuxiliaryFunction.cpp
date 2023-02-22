#include"auxiliary.h"
#include <cmath>
#include <random>

void col_select(vector<vector<double>>& A, int k, int &p)
{//从矩阵第k列的k行及之后选取绝对值最大的那个，用引用参数带回元素的行标
	int n = A.size();
	p = k;
	for (int i = k; i < n; i++)
	{
		if (abs(A[i][k]) > abs(A[p][k]))
		{
			p = i;//记录这个较大值的行标
		}
	}
}

void transpose(vector<vector<double>>& A)
{//矩阵转置要注意，不能所有元素都交换一次，那样最后又回到初始的样子了；只操作一个下三角
	double temp;
	int n = A.size();
	for (int i = 1; i < n; i++)
	{
		for (int j = 0; j < i; j++)
		{
			temp = A[i][j]; A[i][j] = A[j][i]; A[j][i] = temp;
		}
	}
}

vector<double> sign_vector(vector<double> w)
{//虽然符号向量为整值的，但是反正最后都是用于计算的，浮点也无妨
	int n = w.size();
	vector<double> v(n);
	for (int i = 0; i < n; i++)
	{
		v[i] = w[i] > 0 ? 1 : -1;
	}
	return v;
}

double vector_norm_1(vector<double> w)
{
	double sum = 0;
	int n = w.size();
	for (int i = 0; i < n; i++)
	{
		sum += abs(w[i]);
	}
	return sum;
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

double vector_dotproduct(vector<double> z, vector<double> x)
{
	double sum = 0;
	int n = z.size();
	for (int i = 0; i < n; i++)
	{
		sum += z[i] * x[i];
	}
	return sum;
}

vector<double> orthonormal_basis(int j, int n)
{
	vector<double> e(n, 0);
	e[j] = 1;
	return e;
}

void random_value(vector<double>& b)
{
	std::random_device rd;  //Will be used to obtain a seed for the random number engine
	std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
	std::uniform_real_distribution<> dis(0.0, 3.0);
	int N = b.size();
	for (int i = 0; i < N; i++)
	{
		b[i] = dis(gen);
	}
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

vector<double> matrix_vector(vector<vector<double>> A, vector<double> x)
{
	int n = A.size();
	vector<double> result(n, 0);
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			result[i] += A[i][j] * x[j];
		}
	}
	return result;
}