#include"auxiliary.h"
#include <cmath>
#include <ctime>
#include <cstdlib>

void Ex1_Init(vector<vector<double>>& A,vector<double>& b,int N)
{
	//初始化A和b
	for (int i = 0; i < N - 1; i++)
	{
		A[i][i] = 6;
		A[i + 1][i] = 8;
		A[i][i + 1] = 1;
		b[i] = 15;
	}
	A[N - 1][N - 1] = 8;
	b[0] = 7;		//b的首行7
	b[N - 1] = 14;	//b的尾行14
}

void Ex2_Init1(vector<vector<double>>& A, vector<double>& b, int N)
{
	for (int i = 0; i < N - 1; i++)
	{
		A[i][i] = 10;
		A[i + 1][i] = A[i][i + 1] = 1;
	}
	A[N - 1][N - 1] = 10;
	random_value(b);
}

void Ex2_Init2(vector<vector<double>>& A, vector<double>& b, int N)
{
	for (int i = 0; i < N; i++)
	{//构造希尔伯特矩阵
		for (int j = 0; j < N; j++)
		{
			A[i][j] = 1.0 / (i + j + 1);	//注意：这里必须写成1.0，否则两个整型值运算后默认得到整值
			b[i] += A[i][j];
		}
	}
}

void full_select(vector<vector<double>>& A, int k,int &p,int &q)
{//从矩阵k行k列以后的元素中选取绝对值最大那个，用引用参数带回此元素的下标
	int n = A.size();
	p = q = k;
	for (int i = k; i < n; i++)
	{
		for (int j = k; j < n; j++)
		{
			if (abs(A[i][j]) > abs(A[p][q]))
			{
				p = i;//记录这个较大值的行标
				q = j;//记录这个较大值的列标
			}
		}
	}
}

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

void random_value(vector<double>& b)
{
	int N = b.size();
	//srand((unsigned)time(NULL));
	cout << "随机向量b为:" << "\n";
	for (int i = 0; i < N; i++)
	{
		b[i] = rand() % 10;
		cout << b[i] << " ";
		if (i % 10 == 9)cout << "\n";
	}
	cout << "\n";
}

double ResultError(vector<double> y, vector<double> y0)
{//用向量的1范数求误差
	double error = 0;
	int N = y.size();
	for (int i = 0; i < N; i++)
	{
		error += y[i] > y0[i] ? (y[i] - y0[i]) : (y0[i] - y[i]);
	}
	return error;
}