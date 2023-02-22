#include "Exercise.h"
#include "auxiliary.h"

void exercise_1_1()
{
	int N = 84; //矩阵大小
	vector<vector<double>> A(N, vector<double>(N, 0));		//创建一个N行N列的数组，元素类型是double
	vector<double> b(N);	//创建一个N维列向量，元素类型也是double	
	vector<double> sol(N, 1);	//真实解

	Ex1_Init1(A, b, N);//初始化

	double t1 = GetTickCount();
	QRSolve(A, b);	//!!!注意A是方阵的时候，householder变换最后d是只有n-1列的
	double t2 = GetTickCount();
	cout << "求解用时：" << ((t2 - t1)*1.0 / 1000) << "\n";

	cout << "方程组的解如下:" << "\n";
	for (int i = 0; i < N; i++)
	{
		cout.precision(3);
		cout << b[i] << "\t";
		if (i % 10 == 9)cout << "\n";
	}
	cout << "\n";

	cout << "与真实解之差:" << ResultError(b, sol);
}

void exercise_1_2()
{
	int N = 100;
	vector<vector<double>> A(N, vector<double>(N, 0)), A0(N, vector<double>(N, 0));
	vector<double> b(N), b0(N);
	vector<double> y0(N);

	Ex1_Init2(A, b, N);
	A0 = A; b0 = b;	

	double t1 = GetTickCount();
	QRSolve(A, b);	
	double t2 = GetTickCount();
	cout << "求解用时：" << ((t2 - t1)*1.0 / 1000) << "\n";

	cout << "方程组的解如下:" << "\n";
	for (int i = 0; i < N; i++)
	{
		cout.precision(3);
		cout << b[i] << "\t";
		if (i % 10 == 9)cout << "\n";
	}
	cout << "\n";

	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++)
		{
			y0[i] += A0[i][j] * b[j];			//计算y=Ax(x已存放在b中),结果放在y0中
		}
	}
	cout << "Ax-b误差:" << ResultError(y0, b0) << "\n";

}

void exercise_1_3()
{
	int N = 40;
	vector<vector<double>> A(N, vector<double>(N, 0));
	vector<double> b(N, 0);
	vector<double> sol(N, 1);	//真实解

	Ex1_Init3(A, b, N);

	double t1 = GetTickCount();
	QRSolve(A, b);
	double t2 = GetTickCount();
	cout << "求解用时：" << ((t2 - t1)*1.0 / 1000) << "\n";

	cout << "方程组的解如下:" << "\n";
	for (int i = 0; i < N; i++)
	{
		cout.precision(3);
		cout << b[i] << "\t";
		if (i % 10 == 9)cout << "\n";
	}
	cout << "\n";

	cout << "与真实解之差:" << ResultError(b, sol);
}

void exercise_2()
{
	vector<double> t = { -1,-0.75,-0.5,0,0.25,0.5,0.75 };
	vector<double> y = { 1,0.8125,0.75,1,1.3125,1.75,2.3125 };
	int M = t.size();
	vector<vector<double>> A(M, vector<double>(3));		//创建一个M行3列的数组
	vector<double> c;
	
	Ex2_Init(A, t);	//系数矩阵A进行初始化(y就是右端项)

	double t1 = GetTickCount();
	c = QRLeastSquares(A, y);	//求出最小二乘解
	double t2 = GetTickCount();
	cout << "求解用时：" << ((t2 - t1)*1.0 / 1000) << "\n";

	cout << "方程组的解如下:" << "\n";
	for (int i = 0; i < 3; i++)
	{
		cout.precision(3);
		cout << c[i] << "\t";
	}
	cout << "\n";

	cout << "均方误差：" << MSEError(A, y, c);
}

void exercise_3()
{
	vector<vector<double>> A =
	{ {1,4.9176, 1, 3.472, 0.998, 1, 7, 4, 42, 3, 1, 0},
	{1,5.0208, 1, 3.531, 1.5, 2, 7, 4, 62, 1, 1, 0},
	{1,4.5429, 1, 2.275, 1.175, 1, 6, 3, 40,  2, 1, 0},
	{1,4.5573, 1, 4.05, 1.232, 1, 6, 3, 54, 4, 1, 0},
	{1,5.0597, 1, 4.455, 1.121, 1, 6, 3, 42, 3, 1, 0},
	{1,3.891, 1, 4.455, 0.988, 1, 6, 3, 56, 2, 1, 0},
	{1,5.898, 1, 5.85, 1.24, 1, 7, 3, 51, 2, 1,  1},
	{1,5.6039, 1, 9.52, 1.501, 0, 6, 3, 32, 1, 1, 0},
	{1,15.4202, 2.5,  9.8, 3.42, 2, 10, 5, 42, 2, 1, 1},
	{1,14.4598, 2.5, 12.8, 3, 2, 9, 5, 14, 4, 1, 1},
	{1,5.8282, 1, 6.435, 1.225, 2, 6, 3, 32, 1, 1, 0},
	{1,5.3003, 1, 4.9883, 1.552, 1, 6, 3, 30, 1, 2, 0},
	{1,6.2712, 1, 5.52, 0.975, 1, 5, 2, 30, 1, 2, 0},
	{1,5.9592, 1, 6.666, 1.121, 2, 6, 3, 32, 2, 1, 0},
	{1,5.05, 1, 5, 1.02, 0, 5, 2, 46, 4, 1, 1},
	{1,5.6039, 1, 9.52, 1.501, 0, 6, 3, 32, 1, 1, 0},
	{1,8.2464, 1.5, 5.15, 1.664, 2, 8, 4, 50, 4, 1, 0},
	{1,6.6969, 1.5, 6.092, 1.488, 1.5, 7, 3, 22, 1, 1, 1},
	{1,7.7841, 1.5, 7.102, 1.376, 1, 6, 3, 17, 2, 1, 0},
	{1,9.0384, 1, 7.8, 1.5, 1.5, 7, 3, 23, 3, 3, 0},
	{1,5.9894, 1, 5.52, 1.256, 2, 6, 3, 40, 4, 1, 1},
	{1,7.5422, 1.5, 4, 1.69, 1, 6, 3, 22, 1, 1, 0},
	{1,8.7951, 1.5, 9.89, 1.82, 2, 8, 4, 50, 1, 1, 1},
	{1,6.0931, 1.5, 6.7265, 1.652, 1, 6, 3, 44, 4, 1, 0},
	{1,8.3607, 1.5, 9.15, 1.777, 2., 8, 4, 48, 1, 1, 1},
	{1,8.14, 1, 8, 1.504, 2, 7, 3, 3, 1, 3, 0},
	{1,9.1416, 1.5, 7.3262, 1.831, 1.5, 8, 4, 31, 4, 1, 0},
	{1,12, 1.5, 5, 1.2, 2, 6, 3, 30, 3, 1, 1} };
	vector<double> b =
	{ 25.9, 29.5, 27.9, 25.9, 29.9, 29.9, 30.9,
	28.9, 84.9, 82.9, 35.9, 31.5, 31.0, 30.9,
	30.0, 28.9, 36.9, 41.9, 40.5, 43.9, 37.5,
	37.9, 44.5, 37.9, 38.9, 36.9, 45.8, 41.0 };
	int m = A.size();
	int n = A[0].size();
	vector<double> c;

	double t1 = GetTickCount();
	c = QRLeastSquares(A, b);
	double t2 = GetTickCount();
	cout << "求解用时：" << ((t2 - t1)*1.0 / 1000) << "\n";

	cout << "方程组的解如下:" << "\n";
	for (int i = 0; i < n; i++)
	{
		cout.precision(3);
		cout << c[i] << "\t";
	}
	cout << "\n";

	cout << "均方误差：" << MSEError(A, b, c);
}
