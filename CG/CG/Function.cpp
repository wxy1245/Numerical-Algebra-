#include "Function.h"
#include "auxiliary.h"
#include<cmath>

unsigned Jacobi(vector<vector<double>> A, vector<double> b, vector<double> &x, double epsilon)
{//Jacobi迭代：传入的x是初值，最终带回的x是产生的近似解
	int n = x.size();
	vector<double> diagonal(n);
	vector<double> x_pre(n);
	unsigned count = 0;

	//首先求迭代矩阵 B = D^(-1) (L+U) = I - D^(-1) A ,就保存到A里面;以及迭代常数项 g = D^(-1) b ,就保存到b里面
	for (int i = 0; i < n; i++) diagonal[i] = A[i][i];	//用一个数组向量保存矩阵D
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			A[i][j] = -(1.0 / diagonal[i])*A[i][j];
		}
		A[i][i] += 1;
	}
	for (int i = 0; i < n; i++) b[i] = (1.0 / diagonal[i])*b[i];

	//逐次迭代：x(k) = B x(k-1) + g; 注意要保存前一步x的值，为了最后能够终止迭代
	while (1)
	{
		x_pre = x;
		for (int i = 0; i < n; i++)
		{
			x[i] = b[i];
			for (int j = 0; j < n; j++)
			{
				x[i] += A[i][j] * x_pre[j];
			}
		}
		count += 1;
		if (vector_norm_2(vector_minus(x, x_pre)) < 1e-6) break;	//当新得到的x和前一步的x_pre差距很小的时候结束迭代
	}

	return count;
}

unsigned GS(vector<vector<double>> A, vector<double> b, vector<double> &x, double epsilon)
{//Gauss-Seidel迭代：传入的x是初值，最终带回的x是产生的近似解
	int n = x.size();
	vector<double> diagonal(n);
	vector<double> x_pre(n);
	unsigned count = 0;

	for (int i = 0; i < n; i++) diagonal[i] = A[i][i];	//用一个数组向量保存矩阵D
	for (int i = 0; i < n; i++) b[i] = (1.0 / diagonal[i])*b[i];	//b可以用来存放g

	//逐次迭代：x(k) = D^(-1) L x(k) + D^(-1) U x(k-1) + g;
	while (1)
	{
		x_pre = x;
		for (int i = 0; i < n; i++)	//求新x的第i个分量时，用到了已经求出的新x的值前i-1个分量，而i+1到n分量仍用前一步x的值
		{
			x[i] = b[i];
			for (int j = 0; j <= i - 1; j++)
			{
				x[i] += (1.0 / diagonal[i])*(-A[i][j])*x[j];
			}
			for (int j = i + 1; j < n; j++)
			{
				x[i] += (1.0 / diagonal[i])*(-A[i][j])*x_pre[j];
			}
		}
		count += 1;
		if (vector_norm_2(vector_minus(x, x_pre)) < 1e-6) break;	//当新得到的x和前一步的x_pre差距很小的时候结束迭代
	}

	return count;
}

unsigned SOR(vector<vector<double>> A, vector<double> b, double omega, vector<double> &x, double epsilon)
{
	cout << "松弛因子：" << omega << "\n";
	int n = x.size();
	vector<double> diagonal(n);
	vector<double> x_pre(n);
	unsigned count = 0;	//记录迭代次数

	for (int i = 0; i < n; i++) diagonal[i] = A[i][i];	//用一个数组向量保存矩阵D
	for (int i = 0; i < n; i++) b[i] = (1.0 / diagonal[i])*b[i];	//b可以用来存放g

	//逐次迭代：x(k) = ω(D^(-1) L x(k) + D^(-1) U x(k-1) + g) + (1-ω)x(k-1);
	while (1)
	{
		x_pre = x;
		for (int i = 0; i < n; i++)
		{
			//先求出按GS迭代方法得到的分量值
			x[i] = b[i];
			for (int j = 0; j <= i - 1; j++)
			{
				x[i] += (1.0 / diagonal[i])*(-A[i][j])*x[j];
			}
			for (int j = i + 1; j < n; j++)
			{
				x[i] += (1.0 / diagonal[i])*(-A[i][j])*x_pre[j];
			}
			//然后再利用它计算出松弛迭代下的分量值
			x[i] = omega * x[i] + (1 - omega)*x_pre[i];
		}
		count += 1;
		if (vector_norm_2(vector_minus(x, x_pre)) < epsilon) break;	//当新得到的x和前一步的x_pre差距很小的时候结束迭代
	}

	return count;
}

unsigned CG(vector<vector<double>> A, vector<double> b, vector<double> &x, double epsilon)
{//共轭梯度法：求解对称正定方程组；假定x已经赋有初值
	int k = 0, kmax = 1000;	//给一个最大迭代次数限制
	vector<double> r, p, w, x_pre;
	double rho, beta, alpha, rho_pre;
	r = vector_minus(b, matrix_vector(A, x));
	rho = inner_product(r, r);

	while (k < kmax)
	{
		x_pre = x;	//先记录前一个x，方便后面终止迭代
		k = k + 1;
		if (k == 1)
		{
			p = r;
		}
		else
		{
			beta = rho / rho_pre;
			p = vector_plus(r, scale_vector(beta, p));
		}
		w = matrix_vector(A, p);
		alpha = rho / inner_product(p, w);
		x = vector_plus(x, scale_vector(alpha, p));
		if (vector_norm_infinity(vector_minus(x, x_pre)) < epsilon)
		{
			if (vector_norm_infinity(p) < epsilon) k -= 1;	//如果下山方向太小的话，其实最后一次就没有更新，因而迭代次数回退1
			break;
		}
		r = vector_minus(r, scale_vector(alpha, w));
		rho_pre = rho;
		rho = inner_product(r, r);
	}

	return k;	//返回迭代次数(k)
}

unsigned CGTest(vector<vector<double>> A, vector<double> b, vector<double> &x, double epsilon)
{//共轭梯度法：求解对称正定方程组；假定x已经赋有初值
	int k = 0, kmax = 1000;	//给一个最大迭代次数限制
	int n = x.size();
	vector<double> r, p, w, x_pre;
	double rho, beta, alpha, rho_pre;
	r = vector_minus(b, matrix_vector(A, x));
	rho = inner_product(r, r);
	cout << "初始的x值:";
	for (int i = 0; i < n; i++)
	{
		cout.precision(3);
		cout << x[i] << "\t";
	}
	cout << "\n";
	cout << "初始的r值:";
	for (int i = 0; i < n; i++)
	{
		cout.precision(3);
		cout << r[i] << "\t";
	}
	cout << "\n\n";

	while (k < kmax)
	{
		x_pre = x;	//先记录前一个x，方便后面终止迭代
		k = k + 1;
		if (k == 1)
		{
			p = r;
		}
		else
		{
			beta = rho / rho_pre;
			p = vector_plus(r, scale_vector(beta, p));
		}
		w = matrix_vector(A, p);
		alpha = rho / inner_product(p, w);
		x = vector_plus(x, scale_vector(alpha, p));
		if (vector_norm_infinity(vector_minus(x, x_pre)) < epsilon)
		{
			if (vector_norm_infinity(p) < epsilon) k -= 1;	//如果下山方向太小的话，其实最后一次就没有更新，因而迭代次数回退1
			else
			{
				cout << "最后一次产生的结果：" << "\n";
				cout << "下山方向p:";
				for (int i = 0; i < n; i++)
				{
					cout.precision(3);
					cout << p[i] << "\t";
				}
				cout << "\n";
				cout << "更新后的x:";
				for (int i = 0; i < n; i++)
				{
					cout.precision(3);
					cout << x[i] << "\t";
				}
				cout << "\n";
			}
			cout << "\n";
			break;
		}
		r = vector_minus(r, scale_vector(alpha, w));
		rho_pre = rho;
		rho = inner_product(r, r);

		cout << "第" << k << "次迭代后产生的结果:" << "\n";
		cout << "下山方向p:";
		for (int i = 0; i < n; i++)
		{
			cout.precision(3);
			cout << p[i] << "\t";
		}
		cout << "\n";
		cout << "更新后的x:";
		for (int i = 0; i < n; i++)
		{
			cout.precision(3);
			cout << x[i] << "\t";
		}
		cout << "\n";
		cout << "新负梯向r:";
		for (int i = 0; i < n; i++)
		{
			cout.precision(3);
			cout << r[i] << "\t";
		}
		cout << "\n";	
	}

	return k;	//返回迭代次数(k)
}