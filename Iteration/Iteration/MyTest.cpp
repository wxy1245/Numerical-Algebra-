#include<iostream>
#include <vector>
#include"Function.h"
#include"auxiliary.h"
using namespace std;

/*

int main()
{
	int N = 3;
	unsigned count;
	vector<vector<double>> A(N, vector<double>(N, 0));
	vector<double> b(N, 0);
	vector<double> x(N, 0);	//解的迭代初值置为零
	A[0][0] = 1; A[0][1] = 2; A[0][2] = -2; b[0] = 1;
	A[1][0] = 1; A[1][1] = 1; A[1][2] = 1; b[1] = 3;
	A[2][0] = 2; A[2][1] = 2; A[2][2] = 1; b[2] = 5;
	
	count = Jacobi(A, b, x);

	for (int i = 0; i < N; i++)
	{	
		cout << x[i] << "\t";
	}
	cout << "\n";
	
	cout << count;
}

*/