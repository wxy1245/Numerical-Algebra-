#include <iostream>
#include <vector>
#include"Exercise.h"
#include"Function.h"
#include"auxiliary.h"

/*

int main()
{
	int N = 4; //矩阵大小
	vector<vector<double>> A(N, vector<double>(N));		//创建一个N行N列的数组，元素类型是double
	vector<double> b(N,0);
	vector<double> b0(N);
	b0 = b;	//可以整体赋值
	vector<int> u, v;

	for (int i = 0; i < N; i++)
	{
		A[i][i] = i+1;
		b[i] = i;
	}
	A[0][2] = A[2][0] = 1;
	A[1][3] = A[3][1] = 2;
	ModifiedCholeskySolve(A, b);
	for (int i = 0; i < N; i++)
	{
		cout << b0[i] << " ";
	}
	cout << "\n";
	for (int i = 0; i < N; i++)
	{
		cout << b[i] << " ";
	}
}

*/

/*
	for (int i = 0; i < N - 1; i++)
	{
		for (int j = 0; j < N - 1; j++)
		{
			cout << A[i][j] << "\t";
		}
		cout << "\n";
	}
*/

	/*
		gauss_elim_col_pivoting(A, u);
		for (int i = 0; i < N; i++)
		{
			for (int j = 0; j < N; j++)
			{
				cout << A[i][j]<<" ";
			}
			cout << "\n";
		}
	*/

	/*
		//希尔伯特矩阵的打印
		for (int i = 0; i < N; i++)
		{
			for (int j = 0; j < N; j++)
			{
				cout.precision(2);
				cout << A[i][j] << "\t";
			}
			cout << "\n";
		}
		for (int i = 0; i < N; i++)
		{
			cout << b[i] << "\n";
		}
		*/