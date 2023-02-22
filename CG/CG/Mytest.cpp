#include<iostream>
#include<vector>
#include"Function.h"
#include"auxiliary.h"

/*

int main()
{
	int N = 4; //矩阵大小
	vector<vector<double>> A(N, vector<double>(N));		//创建一个N行N列的数组，元素类型是double
	vector<double> b(N, 0);
	vector<double> x(N, 0);
	//vector<double> b0(N);
	//b0 = b;	//可以整体赋值

	for (int i = 0; i < N; i++)
	{
		A[i][i] = i + 1;
		b[i] = i + 2;
	}
	A[0][2] = A[2][0] = 1;
	A[1][3] = A[3][1] = 1;

	cout << "迭代次数：" << CG(A, b, x) << "\n";

	for (int i = 0; i < N; i++)
	{
		cout << x[i] << "\t";
	}
}

*/