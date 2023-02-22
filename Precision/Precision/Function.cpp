#include "Function.h"
#include"auxiliary.h"
#include <cmath>

void forward_subs(vector<vector<double>>& L, vector<double>& b)
{//ǰ�����������Ǿ���
	int j;
	int n = L.size();
	for (j = 0; j < n - 1; j++)
	{
		b[j] = b[j] / L[j][j];
		for (int k = j + 1; k < n; k++)
		{
			b[k] -= b[j] * L[k][j];
		}
	}
	b[n - 1] = b[n - 1] / L[n - 1][n - 1];
}

void forward_subs1(vector<vector<double>>& L, vector<double>& b)
{//�Խ�ԪΪ1ʱ��ǰ����
	int j;
	int n = L.size();
	for (j = 0; j < n - 1; j++)
	{
		for (int k = j + 1; k < n; k++)
		{
			b[k] -= b[j] * L[k][j];
		}
	}
}

void back_subs(vector<vector<double>>& U, vector<double>& y)
{//�ش���
	int j;
	int n = U.size();
	for (j = n - 1; j > 0; j--)
	{
		y[j] = y[j] / U[j][j];
		for (int k = 0; k < j; k++)
		{
			y[k] -= y[j] * U[k][j];
		}
	}
	y[0] = y[0] / U[0][0];
}

void back_subs1(vector<vector<double>>& U, vector<double>& y)
{//�Խ�ԪΪ1�Ļش���
	int j;
	int n = U.size();
	for (j = n - 1; j > 0; j--)
	{
		for (int k = 0; k < j; k++)
		{
			y[k] -= y[j] * U[k][j];
		}
	}
}

void gauss_elim_col_pivoting(vector<vector<double>>& A, vector<int>& u)
{//����ԪGauss��ȥ��
	double temp;
	int n = A.size();
	int p;
	for (int k = 0; k < n - 1; k++)
	{
		col_select(A, k, p);
		for (int j = 0; j < n; j++)
		{
			temp = A[p][j]; A[p][j] = A[k][j]; A[k][j] = temp;	//����k,p��
		}
		u.push_back(p);	//��¼�û�����
		if (A[k][k] != 0)
		{
			for (int i = k + 1; i < n; i++)
			{
				A[i][k] = A[i][k] / A[k][k];	//A(k)���Խ����²���=0���ɴ��Ltiuta��Ԫ
				for (int j = k + 1; j < n; j++)
				{
					A[i][j] -= A[i][k] * A[k][j];
				}
			}
		}
		else
		{
			cout << "�������죬�ֽ��޷���������" << "\n";
			exit(0);
		}
	}
}

void vector_pb(vector<int>& u, vector<double>& b)
{//��������P*b����ѡ��;P=Pr...P1
	double temp;
	int r = u.size();
	for (int i = 0; i < r; i++)
	{
		temp = b[i]; b[i] = b[u[i]]; b[u[i]] = temp;	//����b��i��u(i)����
	}
}

double infinity_norm(vector<vector<double>> A)
{
	int n = A.size();
	double maxline = 0;
	for (int i = 0; i < n; i++)
	{
		double row = 0;
		for (int j = 0; j < n; j++)
		{
			row += abs(A[i][j]);
		}
		if (row > maxline)
		{
			maxline = row;
		}
	}
	return maxline;
}

double inverse_infinity_norm(vector<vector<double>> A, vector<double> x)
{//���Ѿ������A������Ԫ���Ƿֽ�(�����ҵ����۷��������о���P��Ӧ���ǲ���Ҫ��),�����Ax=b�ļ���������,����A��������
	vector<double> v;
	vector<double> w;
	while (1)
	{
		//�������A^t w = x;����������о�������¾���U^t L^t w = x
		transpose(A);	//���A��ת��(!!!����Ϊ��,���������)��������ΪU^t,������ΪL^t
		w = x;
		forward_subs(A, w);	//�����U^t y = x
		back_subs1(A, w);	//�����L^t w = y, ����ǵõ�������w
		v = sign_vector(w);
		//֮�����A z = v;����L U z = v, Ϊ��ʡ�洢������v�洢������z;
		transpose(A);	//!!!��ʱҪ��Aת�û�ȥ
		forward_subs1(A, v);	//�����L y = v
		back_subs(A, v);	//�����U z = y, ����z�Ѿ���ŵ�v������

		if (vector_norm_infinity(v) <= vector_dotproduct(v, x))
		{
			return vector_norm_1(w);
		}
		else
		{
			//���ȡ�����ģ���±�
			int k = 0, n = A.size();
			for (int i = 1; i < n; i++)
			{
				if (abs(v[i]) > abs(v[k]))
				{
					k = i;
				}
			}
			x = orthonormal_basis(k, x.size());
		}
	}
}