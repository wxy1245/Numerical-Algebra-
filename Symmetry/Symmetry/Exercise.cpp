#include"Function.h"
#include"auxiliary.h"
#include"Exercise.h"
#include"Function.h"
#include<cmath>

void exercise_1()
{
	int min_order = 50, max_order = 100;
	vector<double> vals;
	for (int N = min_order; N <= max_order; N += 10)
	{
		cout << "矩阵阶数为" << N << "时：" << endl;

		vector<vector<double>> A(N, vector<double>(N, 0)), Q;
		unsigned k;
		double sigma = N;
		for (int i = 0; i < N - 1; i++)
		{
			A[i][i] = 4;
			A[i + 1][i] = A[i][i + 1] = 1;
		}
		A[N - 1][N - 1] = 4;

		double t1 = GetTickCount();
		k = Gate_Jacobi(A, Q, sigma);
		double t2 = GetTickCount();
		cout << "迭代次数：" << k << endl;
		cout << "求解用时(s)：" << ((t2 - t1)*1.0 / 1000) << "\n";
		cout << "最终矩阵Ak如下" << endl;
		matrix_print(A);
		cout << "变换矩阵Qk如下" << endl;
		matrix_print(Q);
		//cout << endl;
		//vals = get_sort_eigenvalues(A);
		//cout << "全体特征值为：";
		//vector_print(vals);
		cout << endl << endl << endl;
	}

}

void exercise_2()
{
	int N = 100;
	int times = 1000;
	vector<vector<double>> A(N, vector<double>(N, 0));
	double lambda;
	vector<double> z(N, 1);
	double t1, t2;

	for (int i = 0; i < N - 1; i++)
	{
		A[i][i] = 2;
		A[i + 1][i] = A[i][i + 1] = -1;
	}
	A[N - 1][N - 1] = 2;

	lambda = Bisection(A, 1);
	cout <<  "最小的特征值：" << lambda << endl;
	t1 = GetTickCount();
	Inverse_power_method(A, lambda, z, times);
	t2 = GetTickCount();
	cout << "反幂法求相应特征向量迭代次数：" << times << endl;
	cout << "求解用时(ms)：" << ((t2 - t1)*1.0) << "\n";
	cout << "此特征向量为:";
	for (int i = 0; i < N; i++)
	{
		cout.precision(3);
		cout << z[i] << " ";
	}
	cout << endl << endl;

	lambda = Bisection(A, N);
	cout << "最大的特征值：" << lambda << endl;

	t1 = GetTickCount();
	Inverse_power_method(A, lambda, z, times);
	t2 = GetTickCount();
	cout << "反幂法求相应特征向量迭代次数：" << times << endl;
	cout << "求解用时(ms)：" << ((t2 - t1)*1.0) << "\n";

	cout << "此特征向量为:";
	vector_print(z);
	cout << endl << endl;
}
