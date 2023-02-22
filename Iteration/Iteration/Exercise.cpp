#include"Exercise.h"
#include<cmath>

void exercise_1()
{
	double epsilon = 0.0001, a = 0.5;	//epsilon可取的值为1，0.1，0.01, 0.0001；
	int n = 100, N = n - 1;
	double h = 1.0 / n, x;
	unsigned count;
	vector<vector<double>> A(N, vector<double>(N));
	vector<double> b(N);
	vector<double> y(N, 0);	//解的初值置为零
	vector<double> sol(N);	//精确解

	//对矩阵A，右端项b初始化
	for (int i = 0; i < N - 1; i++)
	{
		A[i][i] = -(2 * epsilon + h);
		A[i][i + 1] = epsilon + h;
		A[i + 1][i] = epsilon;
		b[i] = a * h*h;
	}
	A[N - 1][N - 1] = -(2 * epsilon + h);
	b[N - 1] = a * h*h - h - epsilon;

	//求出精确解：sol=(1-a)(1-exp(-x/epsilon))/(1-exp(-1/epsilon))+ax
	for (int i = 0; i < N; i++)
	{
		x = (i + 1)*h;
		sol[i] = (1 - a)*(1 - exp(-x / epsilon)) / (1 - exp(-1 / epsilon)) + a * x;
	}

	//分别采用三种不同迭代方法对其计算求解并且记录运行时间
	double t1 = GetTickCount();
	//count = Jacobi(A, b, y);
	//count = GS(A, b, y);
	count = SOR(A, b, 1.1, y);	//这里面的常数就是选取的松弛迭代因子
	double t2 = GetTickCount();
	cout << "迭代次数：" << count << "\n";
	cout << "求解用时：" << ((t2 - t1)*1.0 / 1000) << "\n";

	cout << "计算解为:" << "\n";
	for (int i = 0; i < N; i++)
	{
		cout.precision(4);
		cout << y[i] << "\t";
		if (i % 10 == 9)cout << "\n";
	}
	cout << "\n";
	cout <<"与精确解误差为:"<< vector_norm_2(vector_minus(y, sol)) << "\n";

}

void exercise_2_1()
{//Jacobi迭代法求解问题
	int n = 60, N = n + 1;	//扩充后的到外围一圈全为1的矩阵
	vector<vector<double>> U(N, vector<double>(N,1));	//解矩阵的初值就选作1(满足边界条件)
	vector<vector<double>> U_pre;
	double h = 1.0 / n;
	double err;	//每次迭代更新过后分析误差是否已经落入接受范围
	unsigned count = 0;	//迭代次数

	double t1 = GetTickCount();
	while (1)
	{
		U_pre = U;	//先保存前一次的结果

		for (int i = 1; i < n; i++)	//计算出新的解矩阵
		{
			for (int j = 1; j < n; j++)
			{
				U[i][j] = (1.0 / (4 + h * h * exp(i*h*j*h))) * (h*h*(i*h + j * h) + U_pre[i - 1][j] + U_pre[i][j - 1] + U_pre[i][j + 1] + U_pre[i + 1][j]);
			}
		}

		count += 1;	//迭代次数加一

		err = 0;	//求出相邻迭代解的误差(解矩阵拉直后看成是向量)并判断是否应该退出
		for (int i = 1; i < n; i++)
		{
			for (int j = 1; j < n; j++)
			{
				err += (U[i][j] - U_pre[i][j])*(U[i][j] - U_pre[i][j]);
			}
		}
		err = sqrt(err);
		if (err < 1e-7) break;
	}
	double t2 = GetTickCount();
	cout << "迭代次数：" << count << "\n";
	cout << "求解用时：" << ((t2 - t1)*1.0 / 1000) << "\n";

	cout << "计算解为" << "\n";
	for (int i = 1; i < n; i++)
	{
		for (int j = 1; j < n; j++)
		{
			cout.precision(4);
			cout << U[i][j] << "\t";
		}
		cout << "\n";
	}
}

void exercise_2_2()
{//Gauss-Seidel迭代法求解问题
	int n = 60, N = n + 1;	//扩充后的到外围一圈全为1的矩阵
	vector<vector<double>> U(N, vector<double>(N, 1));	//解矩阵的初值就选作1(满足边界条件)
	vector<vector<double>> U_pre;
	double h = 1.0 / n;
	double err;	//每次迭代更新过后分析误差是否已经落入接受范围
	unsigned count = 0;	//迭代次数

	double t1 = GetTickCount();
	while (1)
	{
		U_pre = U;	//先保存前一次的结果

		for (int i = 1; i < n; i++)	//计算出新的解矩阵
		{
			for (int j = 1; j < n; j++)
			{
				U[i][j] = (1.0 / (4 + h * h * exp(i*h*j*h))) * (h*h*(i*h + j * h) + U[i - 1][j] + U[i][j - 1] + U_pre[i][j + 1] + U_pre[i + 1][j]);		//只是把Jacobi稍作修改即可:部分项就用已经计算出的新的解就可以了
			}
		}

		count += 1;	//迭代次数加一

		err = 0;	//求出相邻迭代解的误差(解矩阵拉直后看成是向量)并判断是否应该退出
		for (int i = 1; i < n; i++)
		{
			for (int j = 1; j < n; j++)
			{
				err += (U[i][j] - U_pre[i][j])*(U[i][j] - U_pre[i][j]);
			}
		}
		err = sqrt(err);
		if (err < 1e-7) break;
	}
	double t2 = GetTickCount();
	cout << "迭代次数：" << count << "\n";
	cout << "求解用时：" << ((t2 - t1)*1.0 / 1000) << "\n";

	cout << "计算解为" << "\n";
	for (int i = 1; i < n; i++)
	{
		for (int j = 1; j < n; j++)
		{
			cout.precision(4);
			cout << U[i][j] << "\t";
		}
		cout << "\n";
	}
}

void exercise_2_3()
{
	int n = 60, N = n + 1;	//扩充后的到外围一圈全为1的矩阵
	vector<vector<double>> U(N, vector<double>(N, 1));	//解矩阵的初值就选作1(满足边界条件)
	vector<vector<double>> U_pre;
	double h = 1.0 / n;
	double err;	//每次迭代更新过后分析误差是否已经落入接受范围
	unsigned count = 0;	//迭代次数
	double omega = 1.90;	//松弛迭代因子(20阶时：1.75, 40阶时:1.855, 60阶时:1.90)

	double t1 = GetTickCount();
	while (1)
	{
		U_pre = U;	//先保存前一次的结果

		for (int i = 1; i < n; i++)	//计算出新的解矩阵
		{
			for (int j = 1; j < n; j++)
			{
				U[i][j] = (omega / (4 + h * h * exp(i*h*j*h)))*(h*h*(i*h + j * h) + U[i - 1][j] + U[i][j - 1] + U_pre[i][j + 1] + U_pre[i + 1][j]) - (omega - 1)*U_pre[i][j];	//松弛迭代格式
			}
		}

		count += 1;	//迭代次数加一

		err = 0;	//求出相邻迭代解的误差(解矩阵拉直后看成是向量)并判断是否应该退出
		for (int i = 1; i < n; i++)
		{
			for (int j = 1; j < n; j++)
			{
				err += (U[i][j] - U_pre[i][j])*(U[i][j] - U_pre[i][j]);
			}
		}
		err = sqrt(err);
		if (err < 1e-7) break;
	}
	double t2 = GetTickCount();
	cout << "松弛因子：" << omega << "\n";
	cout << "迭代次数：" << count << "\n";
	cout << "求解用时：" << ((t2 - t1)*1.0 / 1000) << "\n";

	cout << "计算解为" << "\n";
	for (int i = 1; i < n; i++)
	{
		for (int j = 1; j < n; j++)
		{
			cout.precision(4);
			cout << U[i][j] << "\t";
		}
		cout << "\n";
	}
}