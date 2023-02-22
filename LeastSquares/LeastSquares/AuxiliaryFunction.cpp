#include"auxiliary.h"
#include <cmath>
#include <random>

void Ex1_Init1(vector<vector<double>>& A, vector<double>& b, int N)
{
	//��ʼ��A��b
	for (int i = 0; i < N - 1; i++)
	{
		A[i][i] = 6;
		A[i + 1][i] = 8;
		A[i][i + 1] = 1;
		b[i] = 15;
	}
	A[N - 1][N - 1] = 8;
	b[0] = 7;		//b������7
	b[N - 1] = 14;	//b��β��14
}

void Ex1_Init2(vector<vector<double>>& A, vector<double>& b, int N)
{
	for (int i = 0; i < N - 1; i++)
	{
		A[i][i] = 10;
		A[i + 1][i] = A[i][i + 1] = 1;
	}
	A[N - 1][N - 1] = 10;
	random_value(b);
}

void Ex1_Init3(vector<vector<double>>& A, vector<double>& b, int N)
{
	for (int i = 0; i < N; i++)
	{//����ϣ�����ؾ���
		for (int j = 0; j < N; j++)
		{
			A[i][j] = 1.0 / (i + j + 1);	//ע�⣺�������д��1.0��������������ֵ�����Ĭ�ϵõ���ֵ
			b[i] += A[i][j];
		}
	}
}

void Ex2_Init(vector<vector<double>>& A, vector<double> t)
{//A������Ϊ3
	int m = A.size();
	for (int i = 0; i < m; i++)
	{
		A[i][0] = t[i] * t[i]; A[i][1] = t[i]; A[i][2] = 1;
	}
}

void random_value(vector<double>& b)
{
	int N = b.size();
	//srand((unsigned)time(NULL));
	cout << "�������bΪ:" << "\n";
	for (int i = 0; i < N; i++)
	{
		b[i] = rand() % 10;
		cout << b[i] << " ";
		if (i % 10 == 9)cout << "\n";
	}
	cout << "\n";
}

double ResultError(vector<double> y, vector<double> y0)
{//��������1���������
	double error = 0;
	int N = y.size();
	for (int i = 0; i < N; i++)
	{
		error += y[i] > y0[i] ? (y[i] - y0[i]) : (y0[i] - y[i]);
	}
	return error;
}

double MSEError(vector<vector<double>> A, vector<double> b, vector<double> x)
{//����b-Ax��2����
	vector<double> y, z;
	y = matrix_vector(A, x);
	z = vector_minus(b, y);
	return vector_norm_2(z);
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
	int m = A.size();	//mȡA������
	int n = A[0].size();//nȡA������
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