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
		cout << "�������Ϊ" << N << "ʱ��" << endl;

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
		cout << "����������" << k << endl;
		cout << "�����ʱ(s)��" << ((t2 - t1)*1.0 / 1000) << "\n";
		cout << "���վ���Ak����" << endl;
		matrix_print(A);
		cout << "�任����Qk����" << endl;
		matrix_print(Q);
		//cout << endl;
		//vals = get_sort_eigenvalues(A);
		//cout << "ȫ������ֵΪ��";
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
	cout <<  "��С������ֵ��" << lambda << endl;
	t1 = GetTickCount();
	Inverse_power_method(A, lambda, z, times);
	t2 = GetTickCount();
	cout << "���ݷ�����Ӧ������������������" << times << endl;
	cout << "�����ʱ(ms)��" << ((t2 - t1)*1.0) << "\n";
	cout << "����������Ϊ:";
	for (int i = 0; i < N; i++)
	{
		cout.precision(3);
		cout << z[i] << " ";
	}
	cout << endl << endl;

	lambda = Bisection(A, N);
	cout << "��������ֵ��" << lambda << endl;

	t1 = GetTickCount();
	Inverse_power_method(A, lambda, z, times);
	t2 = GetTickCount();
	cout << "���ݷ�����Ӧ������������������" << times << endl;
	cout << "�����ʱ(ms)��" << ((t2 - t1)*1.0) << "\n";

	cout << "����������Ϊ:";
	vector_print(z);
	cout << endl << endl;
}
