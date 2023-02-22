#include "Function.h"
#include "auxiliary.h"
#include<cmath>

unsigned Jacobi(vector<vector<double>> A, vector<double> b, vector<double> &x)
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

unsigned GS(vector<vector<double>> A, vector<double> b, vector<double> &x)
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

unsigned SOR(vector<vector<double>> A, vector<double> b, double omega, vector<double> &x)
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
		if (vector_norm_2(vector_minus(x, x_pre)) < 1e-6) break;	//当新得到的x和前一步的x_pre差距很小的时候结束迭代
	}

	return count;
}