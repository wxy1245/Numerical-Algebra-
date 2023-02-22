#include"Function.h"
#include"auxiliary.h"

/*
int main()
{
	int M = 5, N = 4; //矩阵大小
	vector<vector<double>> A(M, vector<double>(N)), U, V;	//创建一个N行N列的数组，元素类型是double
	U = id_matrix(M); V = id_matrix(N);

	A[0][0] = 1; A[0][1] = 1; A[0][2] = 0.5; A[0][3] = -1;
	A[1][0] = 1; A[1][1] = 2; A[1][2] = 3; A[1][3] = 6;
	A[2][0] = 1; A[2][1] = 3; A[2][2] = 6; A[2][3] = 10;
	A[3][0] = -1; A[3][1] = 6; A[3][2] = 10; A[3][3] = 20;
	A[4][0] = -2; A[4][1] = 3.1; A[4][2] = 9.2; A[4][3] = -2;
	
	SVD(A, U, V);

	cout << "变换后的矩阵：" << endl;
	matrix_print(A);

	cout << "正交变换矩阵U：" << endl;
	matrix_print(U);
	cout << "正交变换矩阵V：" << endl;
	matrix_print(V);

	cout << "验证" << endl;
	matrix_print(matrix_product(matrix_product(U, A), transpose(V)));
}
*/


/*
	int M = 4, N = 4; //矩阵大小
	vector<vector<double>> A(M, vector<double>(N)), U, V;	//创建一个N行N列的数组，元素类型是double
	U = id_matrix(M); V = id_matrix(N);

	A[0][0] = 1; A[0][1] = 1; A[0][2] = 1; A[0][3] = 1;
	A[1][0] = 1; A[1][1] = 2; A[1][2] = 3; A[1][3] = 6;
	A[2][0] = 1; A[2][1] = 3; A[2][2] = 6; A[2][3] = 10;
	A[3][0] = 1; A[3][1] = 6; A[3][2] = 10; A[3][3] = 20;
*/

/*

cout << "原矩阵" << endl;
matrix_print(A);

two_diag(A, U, V);
matrix_print(A);

for (int i = 0; i < 3; i++)
{
	Wilkinson_step(A, U, V);
}

*/