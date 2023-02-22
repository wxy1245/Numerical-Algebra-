#include "Function.h"
#include"auxiliary.h"
#include <cmath>

void forward_subs(vector<vector<double>>& L, vector<double>& b)
{//前代法（下三角矩阵）
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
{//对角元为1时的前代法
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
{//回代法
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
{//对角元为1的回代法
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
{//列主元Gauss消去法
	double temp;
	int n = A.size();
	int p;
	for (int k = 0; k < n - 1; k++)
	{
		col_select(A, k, p);
		for (int j = 0; j < n; j++)
		{
			temp = A[p][j]; A[p][j] = A[k][j]; A[k][j] = temp;	//交换k,p行
		}
		u.push_back(p);	//记录置换矩阵
		if (A[k][k] != 0)
		{
			for (int i = k + 1; i < n; i++)
			{
				A[i][k] = A[i][k] / A[k][k];	//A(k)主对角以下部分=0，可存放Ltiuta中元
				for (int j = k + 1; j < n; j++)
				{
					A[i][j] -= A[i][k] * A[k][j];
				}
			}
		}
		else
		{
			cout << "矩阵奇异，分解无法继续进行" << "\n";
			exit(0);
		}
	}
}

void vector_pb(vector<int>& u, vector<double>& b)
{//计算向量P*b【可选】;P=Pr...P1
	double temp;
	int r = u.size();
	for (int i = 0; i < r; i++)
	{
		temp = b[i]; b[i] = b[u[i]]; b[u[i]] = temp;	//交换b的i与u(i)分量
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
{//在已经计算好A的列主元三角分解(经过我的理论分析，排列矩阵P是应该是不需要的),且求出Ax=b的计算解情况下,估计A逆的无穷范数
	vector<double> v;
	vector<double> w;
	while (1)
	{
		//首先求解A^t w = x;无需计入排列矩阵情况下就是U^t L^t w = x
		transpose(A);	//求出A的转置(!!!别以为简单,这里大意了)，下三角为U^t,上三角为L^t
		w = x;
		forward_subs(A, w);	//先求解U^t y = x
		back_subs1(A, w);	//再求解L^t w = y, 这才是得到真正的w
		v = sign_vector(w);
		//之后求解A z = v;即：L U z = v, 为节省存储可以用v存储计算后的z;
		transpose(A);	//!!!此时要将A转置回去
		forward_subs1(A, v);	//先求解L y = v
		back_subs(A, v);	//再求解U z = y, 这样z已经存放到v里面了

		if (vector_norm_infinity(v) <= vector_dotproduct(v, x))
		{
			return vector_norm_1(w);
		}
		else
		{
			//求出取得最大模的下标
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