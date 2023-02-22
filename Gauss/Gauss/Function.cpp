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

void gauss_elim(vector<vector<double>>& A)
{//Gauss消去法
	int n = A.size();
	for (int k = 0; k < n - 1; k++)
	{//依次求解A(k)
		for (int i = k + 1; i < n; i++)
		{
			A[i][k] = A[i][k] / A[k][k];//k列主对角以下0存放lk非零元
			for (int j = k + 1; j < n; j++)
			{
				A[i][j] -= A[i][k] * A[k][j];
			}
		}
	}
}

void gauss_elim_full_pivoting(vector<vector<double>>& A, vector<int>& u, vector<int>& v)
{//全主元Gauss消去法
	double temp;	//注意：这里要用double类型，不要随手就是int！
	int n = A.size();
	int p, q;
	for (int k = 0; k < n - 1; k++)
	{	
		full_select(A, k, p, q);
		for (int j = 0; j < n; j++)
		{
			temp = A[p][j]; A[p][j] = A[k][j]; A[k][j] = temp;	//交换k,p行
		}
		for (int i = 0; i < n; i++)
		{
			temp = A[i][q]; A[i][q] = A[i][k]; A[i][k] = temp;	//交换k,q列
		}
		//u[k] = p; v[k] = q;	//记录交换矩阵
		u.push_back(p); v.push_back(q);
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
			cout << "消去结束" << "\n";
			break;
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

void vector_qb(vector<int>& v, vector<double>& b)
{//计算向量Q*b【可选】;Q=Q1...Qr
	double temp;
	int r = v.size();
	for (int i = r - 1; i >= 0; i--)
	{
		temp = b[i]; b[i] = b[v[i]]; b[v[i]] = temp;	//交换b的i与v(i)分量
	}
}

void cholesky_decomp(vector<vector<double>>& A)
{//对称正定阵标准Cholesky分解
	int n = A.size();
	for (int k = 0; k < n; k++)
	{
		A[k][k] = sqrt(A[k][k]);
		for (int i = k + 1; i < n; i++)
		{
			A[i][k] = A[i][k] / A[k][k];
		}
		for (int j = k + 1; j < n; j++)
		{
			for (int i = j; i < n; i++)
			{
				A[i][j] = A[i][j] - A[i][k] * A[j][k];
			}
		}
	}
}

void modified_cholesky_decomp(vector<vector<double>>& A)
{//改进的平方根法
	int n = A.size();
	vector<double> v(n);	//设置double类型的临时向量v
	for (int j = 0; j < n; j++)
	{
		for (int k = 0; k <= j - 1; k++)
		{
			v[k] = A[j][k] * A[k][k];
		}
		for (int k = 0; k <= j - 1; k++)
		{
			A[j][j] -= A[j][k]*v[k];	//主对角线上存储对角元
		}
		for (int i = j + 1; i < n; i++)
		{
			for (int k = 0; k <= j - 1; k++)
			{
				A[i][j] -= A[i][k] * v[k];	//这里逻辑要注意:先进行求和，最后除以对角元
			}
			A[i][j] /= A[j][j];
		}
	}
}

void matrix_DLT(vector<vector<double>>& A)
{//计算矩阵D*L^T【可选】,置于A的上三角部分；其中D源于A的对角，L则源于A的主对角线以下
	int N = A.size();
	for (int i = 0; i < N - 1; i++)
	{
		for (int j = i + 1; j < N; j++)
		{
			A[i][j] = A[i][i] * A[j][i];
		}
	}
}

void GaussSolve(vector<vector<double>>& A, vector<double>& b)
{//Gauss消去法求解线性方程组
	int N = A.size();
	gauss_elim(A);	//A完成LU分解之后，就可以先后求解两个方程组Ly=b,Ux=y来得到最终解x
	forward_subs1(A, b);	//注意到L为单位下三角阵，可用对角元为1的前代法
	back_subs(A, b);		//注意前代法求解之后，b被y所取代;而回代法过后，b存储的就是x(最终解)
}

void FullGaussSolve(vector<vector<double>>& A, vector<double>& b)
{//全主元Gauss消去法求解线性方程组
	int N = A.size();
	vector<int> u, v;
	gauss_elim_full_pivoting(A, u, v);	//PAQ=LU完成分解之后，先后求解Ly=Pb,Uz=y,最后由x=Qz得到结果
	vector_pb(u, b);
	forward_subs1(A, b);
	back_subs(A, b);
	vector_qb(v, b);
}

void ColGaussSolve(vector<vector<double>>& A, vector<double>& b)
{//列主元Gauss消去法求解线性方程组
	int N = A.size();
	vector<int> u;
	gauss_elim_col_pivoting(A, u);	//PA=LU完成分解之后，先后求解Ly=Pb,Ux=y
	vector_pb(u, b);
	forward_subs1(A, b);
	back_subs(A, b);
}

void CholeskySolve(vector<vector<double>>& A, vector<double>& b)
{
	int N = A.size();
	cholesky_decomp(A);	//计算A的cholesky分解
	for (int i = 0; i < N - 1; i++)	//把A的上三角部分改造为L的转置
	{
		for (int j = i + 1; j < N; j++)
		{
			A[i][j] = A[j][i];
		}
	}
	forward_subs(A, b);	//求解Ly=b
	back_subs(A, b);	//求解L^t x=y
}

void ModifiedCholeskySolve(vector<vector<double>>& A, vector<double>& b)
{
	int N = A.size();
	modified_cholesky_decomp(A);	//先进行改进的cholesky分解
	matrix_DLT(A);	//A的上三角部分存放D*L^T
	forward_subs1(A, b);	//求解Ly=b
	back_subs(A, b);	//求解Ux=y
}
