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

void householder(vector<double> x, vector<double> &v, double &beta)
{
	int n = x.size();
	double eta = vector_norm_infinity(x);

	for (int i = 0; i < n; i++)
	{
		x[i] /= eta;
	}
	double sigma = 0;
	for (int i = 1; i < n; i++)
	{
		sigma += x[i] * x[i];
		v[i] = x[i];
	}

	if (sigma < 1.0e-10)
	{
		beta = 0;
	}
	else
	{
		double alpha = sqrt(x[0] * x[0] + sigma);
		if (x[0] <= 0)
		{
			v[0] = x[0] - alpha;
		}
		else
		{
			v[0] = -sigma / (x[0] + alpha);
		}
		beta = 2 * v[0] * v[0] / (sigma + v[0] * v[0]);
		for (int i = 1; i < n; i++)
		{
			v[i] /= v[0];
		}
		v[0] = 1;	//!!!注意：v进行规格化时，先用v[0]规格化其他元素，最后才把v[0]置1
	}
}

void matrix_householder_simplify(vector<vector<double>> &A, vector<double> v, double beta)
{
	int m = A.size();	//行数
	int n = A[0].size();	//列数
	vector<double> w(n, 0);

	//先计算w=beta*A^t v
	for (int i = 0; i < n; i++)
	{
		for (int k = 0; k < m; k++)
		{
			w[i] += A[k][i] * v[k];
		}
		w[i] *= beta;
	}

	//H A = A - v w^T,直接就保存在A里面
	for (int i = 0; i < m; i++)
	{
		for (int j = 0; j < n; j++)
		{
			A[i][j] -= v[i] * w[j];
		}
	}
}

void matrix_householder_simplify_j(vector<vector<double>> &A, vector<double> v, double beta, int j)
{//对A的j行j列以后的部分进行一步householder化简
	int m = A.size();	//总行数
	int n = A[0].size();	//总列数
	vector<double> w(n - j, 0);	//w的列数只需要n-j

	//先计算w=beta*A^t v
	for (int i = j; i < n; i++)
	{
		for (int k = j; k < m; k++) 
		{
			w[i - j] += A[k][i] * v[k - j];
		}
		w[i - j] *= beta;
	}

	//H A = A - v w^T,直接就保存在A里面
	for (int k = j; k < m; k++)
	{
		for (int i = j; i < n; i++)
		{
			A[k][i] -= v[k - j] * w[i - j];
		}
	}
}

vector<double> QR_decomposition(vector<vector<double>> &A)
{
	int m = A.size();	//行数
	int n = A[0].size();	//列数
	vector<double> v, temp, d;	//temp:用于取A的第j列(j行以后)计算householder变换的v,β;d:存放各次变换的β
	double beta;
	for (int j = 0; j < n; j++)
	{
		if (j < m - 1)	//j=m-1时就已经到达行底，则无需操作
		{
			v.resize(m - j); temp.resize(m - j);
			for (int k = j; k < m; k++)
			{
				temp[k - j] = A[k][j];
			}
			householder(temp, v, beta);
			matrix_householder_simplify_j(A, v, beta, j);
			d.push_back(beta);
			for (int k = j + 1; k < m; k++)
			{
				A[k][j] = v[k - j];
			}
		}
	}
	return d;
}

void Qt_b(vector<vector<double>> A, vector<double>& b, vector<double> d)
{//计算Q^T b = Hn...H1 b,n:变换次数
	int m = A.size();
	int n = d.size();	//!!!注意n取总变换的次数(d的维数,即β个数),而非A的列数(如果A行数≤A列数,则变换次数是<列数,故取变换次数较有一般性)
	vector<double> v;
	double tempsum;
	vector<double> temp;
	for (int j = 0; j < n; j++)	//依次计算b1=H1b,...,bn=Hnbn-1
	{
		//先把这一次的householder变换用的vj取出来(大小为m-j,首元为1)
		v.resize(m - j);
		v[0] = 1;
		for (int i = 1; i < m - j; i++)
		{
			v[i] = A[j + i][j];
		}

		//对b的后m-j项进行householder变换
		temp.resize(m - j);	//用于临时存放b变化后的值
		for (int i = 0; i < m - j; i++)	//把b临时存放进去
		{
			temp[i] = b[i + j];
		}
		for (int i = j; i < m; i++)	//先求出temp,即变换后的b
		{
			tempsum = 0;
			for (int k = j; k < m; k++)
			{
				tempsum += v[i - j] * v[k - j] * b[k];
			}
			temp[i - j] -= d[j] * tempsum;	//注意：这里逐个求b的时候，由于b在每次求和的时候都要用，过程中不可改变。所以必须用另外一个临时变量temp对其保存
		}
		for (int i = j; i < m; i++)	//最后b赋值temp完成一次操作
		{
			b[i] = temp[i - j];
		}
	}
}

void QRSolve(vector<vector<double>> &A, vector<double> &b)
{//给出A为可逆方阵时的QR分解求解Ax=b的方法
	vector<double> d;
	d = QR_decomposition(A);	//！！！注意：这里d只有n-1维; QR分解过后,Q=H1...Hn-1,且记录H变换的v保存于A的下三角,β保存于d；R保存于A的上三角；
	//求解Q R x = b,即R x = Q^T b = Hn-1...H1 b；
	Qt_b(A, b, d);	//先计算Q^T b,保存在b里面
	back_subs(A, b);	//再求解Rx = Q^t b,最终解保存在b里面
}

vector<double> QRLeastSquares(vector<vector<double>> A, vector<double> b)
{//超定情形(m>n),求最小二乘问题的解;另外,上述QRSolve可以看做此的特例
	int n = A[0].size();	//取A的列数
	vector<double> d;

	//求解步骤: 1.QR分解 2.c1 = Q1^T b(注意c1通过取Q^T b前n行就可以得到了) 3.R x = c1
	d = QR_decomposition(A);

	Qt_b(A, b, d);
	vector<double>c1(b.begin(), b.begin() + n);	//切片操作,取前n行

	vector<vector<double>> A1(n, vector<double>(n));	//为了能够应用回代法，需要取A的方阵那块(前n行n列,上三角部分即为R)
	for (int i = 0; i < n; i++)	
	{
		for (int j = 0; j < n; j++)
		{
			A1[i][j] = A[i][j];
		}
	}
	back_subs(A1, c1);
	return c1;	//返回解(n维向量)
}