#include<iostream>
#include <vector>
#include"Function.h"
#include"auxiliary.h"
using namespace std;

/*

int main()
{
	int M = 3, N = 2;
	//vector<double> x(M), v(M);
	//double beta;
	vector<double>d;


	
	x[0] = 0; x[1] = 3; x[2] = 4;
	householder(x, v, beta);
	for (int i = 0; i < M; i++)
	{
		cout << v[i] << "\t";
	}
	cout << beta << "\n";
	


	vector<vector<double>> A(M);
	for (int i = 0; i < M; i++)
	{
		A[i].resize(N);
	}

	A[0][0] = 0; A[0][1] = 1;
	A[1][0] = 3; A[1][1] = 1;
	A[2][0] = 4; A[2][1] = 1;

	//matrix_householder_simplify_j(A, v, beta, 0);

	d = QR_decomposition(A);

	for (int i = 0; i < M; i++)
	{
		for (int j = 0; j < N; j++)
		{
			cout << A[i][j] << "\t";
		}
		cout << "\n";
	}
	for (int j = 0; j < N; j++)
	{
		cout << d[j] << "\t";
	}
	cout << "\n";

	vector<double>b(M, 1);
	Qt_b(A, b, d);

	for (int i = 0; i < M; i++)
	{
		cout << b[i] << "\t";
	}
}

*/