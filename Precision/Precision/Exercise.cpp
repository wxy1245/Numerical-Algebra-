#include"Exercise.h"
#include"auxiliary.h"

void ColGaussSolve(vector<vector<double>>& A, vector<double>& b)
{//列主元Gauss消去法求解线性方程组
	int N = A.size();
	vector<int> u;
	gauss_elim_col_pivoting(A, u);	//PA=LU完成分解之后，先后求解Ly=Pb,Ux=y
	vector_pb(u, b);
	forward_subs1(A, b);
	back_subs(A, b);
}

double Col_Precision(vector<vector<double>>& A, vector<double>& b)
{
	double mu = infinity_norm(A);
	double beta = vector_norm_infinity(b);
	double gamma, nu;
	vector<vector<double>>A0 = A;	//r的求解需要保留原来的矩阵与向量
	vector<double>b0 = b, r;

	ColGaussSolve(A, b);	//A放的是三角分解，b放的是计算解
	r = vector_minus(b0, matrix_vector(A0, b));
	gamma = vector_norm_infinity(r);
	nu = inverse_infinity_norm(A, b);
	return nu * mu * gamma / beta;
}

void exercise_1()
{
	int N;	//希尔伯特矩阵阶数
	for (N = 2; N <= 20; N++)
	{
		vector<vector<double>> A(N, vector<double>(N, 0));
		vector<double> b(N, 0);

		for (int i = 0; i < N; i++)
		{
			for (int j = 0; j < N; j++)
			{
				A[i][j] = 1.0 / (i + j + 1);
				b[i] += A[i][j];
			}
		}

		double mu = infinity_norm(A);	//!!!注意A的无穷范数用的原来的希尔伯特矩阵来求，而逆的无穷范数则是用到列主元消去得到的矩阵
		ColGaussSolve(A, b);
		cout << N << "阶Hilbert矩阵的无穷范数条件数为:" << mu * inverse_infinity_norm(A, b) << "\n";
	}
}

void exercise_2()
{
	int N;	//矩阵A的阶数
	for (N = 5; N <= 30; N++)
	{
		vector<vector<double>> A(N, vector<double>(N, 1));
		vector<double> x(N), b(N, 0);
		double accuracy, estimation;

		for (int i = 1; i < N - 1; i++)	//生成矩阵A
		{
			for (int j = 0; j < i; j++)
			{
				A[i][j] = -1; A[j][i] = 0;
			}
		}
		for (int j = 0; j < N - 1; j++)
		{
			A[N - 1][j] = -1;
		}
		random_value(x);	//随机选取一个x,它也是后面的真实解

		for (int i = 0; i < N; i++)	//因为知道矩阵A的特点，可以简单计算b
		{
			for (int j = 0; j < i; j++)
			{
				b[i] -= x[j];
			}
			b[i] += x[i] + x[N - 1];
		}
		b[N - 1] -= x[N - 1];

		estimation = Col_Precision(A, b);	//列主元Gauss消去法及解的精度估计,注意此步骤后A，b都已经改变
		accuracy = vector_norm_infinity(vector_minus(x, b)) / vector_norm_infinity(x);
		cout.setf(ios::scientific);
		cout << N << "阶矩阵,解的相对误差为：" << accuracy << ", 估计误差界为：" << estimation << "\n";
	}
}