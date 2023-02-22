#include "Function.h"
#include "auxiliary.h"
#include<cmath>
#define ERROR 1e-5

void forward_subs(vector<vector<double>> L, vector<double>& b)
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

void forward_subs1(vector<vector<double>> L, vector<double>& b)
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

void back_subs(vector<vector<double>> U, vector<double>& y)
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

void back_subs1(vector<vector<double>> U, vector<double>& y)
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

void vector_pb(vector<int> u, vector<double>& b)
{//计算向量P*b【可选】;P=Pr...P1
	double temp;
	int r = u.size();
	for (int i = 0; i < r; i++)
	{
		temp = b[i]; b[i] = b[u[i]]; b[u[i]] = temp;	//交换b的i与u(i)分量
	}
}

void Jacobi_simplify_once(vector<vector<double>>& A, int p, int q, double c, double s)
{
	int n = A.size();
	double a_pp = A[p][p], a_pq = A[p][q], a_qp = A[q][p], a_qq = A[q][q];
	A[p][p] = c * c * a_pp - 2 * s * c * a_pq + s * s * a_qq;
	A[q][q] = s * s * a_pp + 2 * s * c * a_pq + c * c * a_qq;
	A[p][q] = A[q][p] = (c*c - s * s)*a_pq + s * c*(a_pp - a_qq);
	for (int i = 0; i < n; i++)
	{
		if (i != p && i != q)
		{
			A[i][p] = c * A[p][i] - s * A[q][i];
			A[i][q] = s * A[p][i] + c * A[q][i];
			A[p][i] = A[i][p];
			A[q][i] = A[i][q];
		}
	}
}

void update_Q(vector<vector<double>> &A, int p, int q, double c, double s)
{//无需进行矩阵乘法，只需计算p，q列（有变化的两列）即可
	int n = A.size();
	double a_ip, a_iq;
	for (int i = 0; i < n; i++)
	{
		a_ip = A[i][p]; a_iq = A[i][q];
		A[i][p] = c * a_ip - s * a_iq;
		A[i][q] = s * a_ip + c * a_iq;
	}
}

unsigned one_scan(vector<vector<double>> &A, vector<vector<double>> &Q, double delta)
{
	int n = A.size();
	unsigned one_scan_count = 0;
	double c, s, t, tau;
	for (int i = 0; i < n; i++)
	{
		for (int j = i + 1; j < n; j++)
		{
			if (abs(A[i][j]) > delta)
			{
				tau = (A[j][j] - A[i][i]) / (2.0*A[i][j]);
				t = tau >= 0 ? (1.0 / (tau + sqrt(1 + tau * tau))) : (-1.0 / (-tau + sqrt(1 + tau * tau)));
				c = 1.0 / sqrt(1 + t * t); s = t * c;	//计算出Jacobi矩阵的参数
				Jacobi_simplify_once(A, i, j, c, s);	//Jacobi一次约化
				update_Q(Q, i, j, c, s);	//更新变换矩阵
				one_scan_count++;	//迭代次数+1
			}
		}
	}
	return one_scan_count;
}

unsigned Gate_Jacobi(vector<vector<double>>& A, vector<vector<double>> &Q, double sigma)
{
	int n = A.size();
	double delta;
	unsigned count = 0;
	unsigned increment;
	delta = NotDiagonalNorm(A);
	Q = id_matrix(n);

	while (delta > ERROR)	//当关值充分小的时候可以跳出
	{
		increment = one_scan(A, Q, delta);
		count += increment;
		while (increment > 0)	//当所有非对角元没有全部过关(一次扫描后仍有增量)时，继续扫描
		{
			increment = one_scan(A, Q, delta);
			count += increment;
		}
		delta /= sigma;	//都过关之后，减小关值后重新扫描
	}

	return count;
}

unsigned sign_changes(vector<vector<double>> T, double x)
{
	int n = T.size();
	unsigned negatives = 0;
	double p_pre_pre = 1, p_pre = T[0][0] - x, p_cur;
	if (p_pre < 0)negatives++;

	for (int i = 2; i <= n; i++)
	{
		p_cur = (T[i - 1][i - 1] - x)*p_pre - T[i - 2][i - 1] * T[i - 2][i - 1] * p_pre_pre;
		if (p_cur / p_pre < 0)negatives++;
		p_pre_pre = p_pre;
		p_pre = p_cur;
	}

	return negatives;
}

double Bisection(vector<vector<double>> A, unsigned m)
{
	double left = -matrix_norm_infinity(A);
	double right = -left;
	double middle;
	unsigned s;
	while (right - left > ERROR)
	{
		middle = (left + right) / 2.0;
		s = sign_changes(A, middle);
		if (s >= m)right = middle;
		else left = middle;
	}
	middle = (left + right) / 2.0;
	return middle;
}

void Inverse_power_method(vector<vector<double>> A, double mu, vector<double> &z, int &times)
{
	int n = A.size();
	vector<int> u;
	vector<double> z_pre(n, 0);
	int count = 0;
	for (int i = 0; i < n; i++)	//求A-μI
	{
		A[i][i] -= mu;
	}
	gauss_elim_col_pivoting(A, u);	//对A-μI进行PA=LU分解（只需进行一次）

	z = scale_vector(1.0 / vector_norm_2(z), z);	//先将z做一下标准化

	while (count < times && abs(vector_norm_1(z)- vector_norm_1(z_pre)) > ERROR)	//当迭代次数足够多或者特征向量几乎不变时就可以停止迭代
	{
		z_pre = z;
		vector_pb(u, z);
		forward_subs1(A, z);
		back_subs(A, z);
		z = scale_vector(1.0 / vector_norm_2(z), z);
		count++;
	}

	times = count;	//times最后换成真实迭代次数
}



