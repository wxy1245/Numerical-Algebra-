#include"auxiliary.h"
#include <cmath>
#include <random>

void col_select(vector<vector<double>>& A, int k, int &p)
{//�Ӿ����k�е�k�м�֮��ѡȡ����ֵ�����Ǹ��������ò�������Ԫ�ص��б�
	int n = A.size();
	p = k;
	for (int i = k; i < n; i++)
	{
		if (abs(A[i][k]) > abs(A[p][k]))
		{
			p = i;//��¼����ϴ�ֵ���б�
		}
	}
}

void transpose(vector<vector<double>>& A)
{//����ת��Ҫע�⣬��������Ԫ�ض�����һ�Σ���������ֻص���ʼ�������ˣ�ֻ����һ��������
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
{//��Ȼ��������Ϊ��ֵ�ģ����Ƿ�����������ڼ���ģ�����Ҳ�޷�
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
	return maximum;	//���ص���ȡ�����ģ���±�
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