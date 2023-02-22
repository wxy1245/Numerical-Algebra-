#include"auxiliary.h"
#include <cmath>
#include <ctime>
#include <cstdlib>

void Ex1_Init(vector<vector<double>>& A,vector<double>& b,int N)
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
	{//����ϣ�����ؾ���
		for (int j = 0; j < N; j++)
		{
			A[i][j] = 1.0 / (i + j + 1);	//ע�⣺�������д��1.0��������������ֵ�����Ĭ�ϵõ���ֵ
			b[i] += A[i][j];
		}
	}
}

void full_select(vector<vector<double>>& A, int k,int &p,int &q)
{//�Ӿ���k��k���Ժ��Ԫ����ѡȡ����ֵ����Ǹ��������ò������ش�Ԫ�ص��±�
	int n = A.size();
	p = q = k;
	for (int i = k; i < n; i++)
	{
		for (int j = k; j < n; j++)
		{
			if (abs(A[i][j]) > abs(A[p][q]))
			{
				p = i;//��¼����ϴ�ֵ���б�
				q = j;//��¼����ϴ�ֵ���б�
			}
		}
	}
}

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