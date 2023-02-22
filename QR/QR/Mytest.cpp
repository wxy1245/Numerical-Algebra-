#include<iostream>
#include<vector>
#include<complex>
#include"Function.h"
#include"auxiliary.h"

/*
int main()
{
	int N = 6; //矩阵大小
	vector<vector<double>> A(N, vector<double>(N)), P;	//创建一个N行N列的数组，元素类型是double
	vector<complex<double>> z;
	A[0][0] = 1.1908; A[0][1] = -1.0565; A[0][2] = -2.1707; A[0][3] = 0.5913; A[0][4] = 0; A[0][5] = 0.7310;
	A[1][0] = -1.2025; A[1][1] = 1.4151; A[1][2] = -0.0592; A[1][3] = -0.6436; A[1][4] = -0.3179; A[1][5] = 0.5779;
	A[2][0] = -0.0198; A[2][1] = -0.8051; A[2][2] = -1.0106; A[2][3] = 0.3803; A[2][4] = 1.0950; A[2][5] = 0.0403;
	A[3][0] = -0.1567; A[3][1] = 0.5287; A[3][2] = 0.6145; A[3][3] = -1.0091; A[3][4] = -1.8740; A[3][5] = 0.6771;
	A[4][0] = -1.6041; A[4][1] = 0.2193; A[4][2] = 0.5077; A[4][3] = -0.0195; A[4][4] = 0.4282; A[4][5] = 0.5689;
	A[5][0] = 0.2573; A[5][1] = -0.9219; A[5][2] = 1.6924; A[5][3] = -0.0482; A[5][4] = 0.8956; A[5][5] = -0.2556;
	eigenvalues(A, z);
	for (int i = 0; i < N; i++)
	{
		cout << z[i] << "\n";
	}
}
*/

/*
hessenberg_decomp(A, V, b);
for (int i = 0; i < N; i++)
{
	for (int j = 0; j < N; j++)
	{
		cout << A[i][j] << "\t";
	}
	cout << "\n";
}

cout << "\n" << "QR迭代" << "\n";

for (int i = 0; i < 3; i++)
{
	P = two_step_displacement_qr(A);
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < i; j++)
		{
			if (A[i][j] < 1e-8)A[i][j] = 0;
		}
	}

	cout << "变换矩阵:" << "\n";
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++)
		{
			cout.precision(3);
			cout << P[i][j] << "\t";
		}
		cout << "\n";
	}

	cout << "一次迭代过后的矩阵" << "\n";
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++)
		{
			cout.precision(3);
			cout << A[i][j] << "\t";
		}
		cout << "\n";
	}
}
*/


	/*
		int N = 2; //矩阵大小
		vector<vector<double>> A(N, vector<double>(N)), P;	//创建一个N行N列的数组，元素类型是double
		vector<vector<double>> V;
		vector<double> b;

		A[0][0] = 4.48; A[0][1] = -6.25;
		A[1][0] = -4.61; A[1][1] = 11.2;

		for (int i = 0; i < 1; i++)
		{
			P = two_step_displacement_qr(A);
			cout << "变换矩阵:" << "\n";
			for (int i = 0; i < N; i++)
			{
				for (int j = 0; j < N; j++)
				{
					cout.precision(3);
					cout << P[i][j] << "\t";
				}
				cout << "\n";
			}
			cout << "QR得到的矩阵" << "\n";
			for (int i = 0; i < N; i++)
			{
				for (int j = 0; j < N; j++)
				{
					cout.precision(3);
					cout << A[i][j] << "\t";
				}
				cout << "\n";
			}
			matrix_print(matrix_product(matrix_product(P, A), transpose(P)));
		}
		*/


		/*
			P = implicit_qr(A);
			cout << "变换矩阵:" << "\n";
			for (int i = 0; i < N; i++)
			{
				for (int j = 0; j < N; j++)
				{
					cout.precision(3);
					cout << P[i][j] << "\t";
				}
				cout << "\n";
			}
			cout << "QR得到的矩阵" << "\n";
			for (int i = 0; i < N; i++)
			{
				for (int j = 0; j < N; j++)
				{
					cout.precision(3);
					cout << A[i][j] << "\t";
				}
				cout << "\n";
			}

			cout << "重构原矩阵" << "\n";
			matrix_print(matrix_product(matrix_product(P, A), transpose(P)));
			*/





	