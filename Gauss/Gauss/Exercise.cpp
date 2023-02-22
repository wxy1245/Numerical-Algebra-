#include"Exercise.h"
#include"auxiliary.h"

//!!!注意：对于2.1,3.1才需要用Ax-b，而且A也需要保留才行
//!!!已经简洁又完善的程序(如GaussSolve)可以减少在内部添加不必要功能。

void exercise_1()
{
	int N = 84; //矩阵大小
	vector<vector<double>> A(N, vector<double>(N, 0));		//创建一个N行N列的数组，元素类型是double
	vector<double> b(N);	//创建一个N维列向量，元素类型也是double	
	vector<double> sol(N, 1);	//真实解

	Ex1_Init(A, b, N);//初始化

	double t1 = GetTickCount();
	//GaussSolve(A, b);
	//FullGaussSolve(A, b);
	ColGaussSolve(A, b);
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

void exercise_2_1()
{
	int N = 100; 
	vector<vector<double>> A(N, vector<double>(N, 0)), A0(N, vector<double>(N, 0));
	vector<double> b(N), b0(N);
	vector<double> y0(N);	

	Ex2_Init1(A, b, N);
	A0 = A; b0 = b;	//记录原始的A，b后续计算Ax-b的误差(注意:vector向量是可以整体赋值的!)
	
	double t1 = GetTickCount();
	//CholeskySolve(A, b);	//平方根法求解
	ModifiedCholeskySolve(A, b);	//改进的平方根法求解
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
	cout<< "Ax-b误差:" << ResultError(y0, b0) << "\n";
}

void exercise_2_2()
{
	int N = 40;
	vector<vector<double>> A(N, vector<double>(N, 0));
	vector<double> b(N,0);
	vector<double> sol(N, 1);	//真实解

	Ex2_Init2(A, b, N);

	double t1 = GetTickCount();
	//CholeskySolve(A, b);
	ModifiedCholeskySolve(A, b);
	double t2 = GetTickCount();
	cout << "求解用时：" << ((t2 - t1)*1.0 / 1000) << "\n";

	cout << "方程组的解如下:" << "\n";
	for (int i = 0; i < N; i++)
	{
		cout.precision(3);
		cout << b[i] << "\t";
	}
	cout << "\n";

	cout << "与真实解之差:" << ResultError(b, sol);
}

void exercise_3_1()
{//用第一题的方式求解第二题的第一个方程组
	int N = 100;
	vector<vector<double>> A(N, vector<double>(N, 0)), A0(N, vector<double>(N, 0));
	vector<double> b(N), b0(N);
	vector<double> y0(N);
	
	Ex2_Init1(A, b, N);
	A0 = A; b0 = b;	//记录原始的A，b后续计算Ax-b的误差

	double t1 = GetTickCount();
	//GaussSolve(A, b);
	//FullGaussSolve(A, b);
	ColGaussSolve(A, b);
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

void exercise_3_2()
{//用第一题的方式求解第二题的第二个方程组
	int N = 15;
	vector<vector<double>> A(N, vector<double>(N, 0));
	vector<double> b(N);
	vector<double> sol(N, 1);	//真实解

	Ex2_Init2(A, b, N);
	double t1 = GetTickCount();
	//GaussSolve(A, b);
	//FullGaussSolve(A, b);
	ColGaussSolve(A, b);
	double t2 = GetTickCount();
	cout << "求解用时：" << ((t2 - t1)*1.0 / 1000) << "\n";

	cout << "方程组的解如下:" << "\n";
	for (int i = 0; i < N; i++)
	{
		cout.precision(5);
		cout << b[i] << "\t";
		if (i % 10 == 9)cout << "\n";
	}
	cout << "\n";

	cout << "与真实解之差:" << ResultError(b, sol);
}