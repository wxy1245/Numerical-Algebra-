#include"Exercise.h"
#include"auxiliary.h"
#include"Function.h"
#include<cmath>

void exercise_1_1()
{
	cout << "Exercise_1_1" << endl << endl;

	vector<double> a = { -1,5,-3 };
	int N = 3;
	vector<double> u0 = { 0,0,1 };	//���������һ�����ȫΪ1��������ʼ���õ��Ľ����1

	vector<vector<double>> A(N, vector<double>(N,0));
	A[0] = a;
	for (int i = 1; i < N; i++)A[i][i - 1] = 1;

	int times = 1000;

	double t1 = GetTickCount();
	double x = power_method(A, u0, times);
	double t2 = GetTickCount();
	cout << "�����ʱ��" << ((t2 - t1)*1.0 / 1000) << "\n";
	cout << "����������" << times << "\n";
	cout << "���ģ����ֵ��" << x << "\n";
}

void exercise_1_2()
{
	cout << "Exercise_1_2" << endl << endl;

	vector<double> a = { 0,3,1 };
	int N = 3;
	vector<double>u0(N, 1);

	vector<vector<double>> A(N, vector<double>(N, 0));
	A[0] = a;
	for (int i = 1; i < N; i++)A[i][i - 1] = 1;

	int times = 1000;

	double t1 = GetTickCount();
	double x = power_method(A, u0, times);
	double t2 = GetTickCount();
	cout << "�����ʱ��" << ((t2 - t1)*1.0 / 1000) << "\n";
	cout << "����������" << times << "\n";
	cout << "���ģ����ֵ��" << x << "\n";
}

void exercise_1_3()
{
	cout << "Exercise_1_3" << endl << endl;

	vector<double> a = { -101,-208.01,-10891.01, -9802.08,-79108.9,99902,-790,1000 };
	int N = 8;
	vector<double>u0 = { 0,0,0,0,0,0,0,1 };	//���������һ�����ȫΪ1��������ʼ���õ��Ľ����1

	vector<vector<double>> A(N, vector<double>(N, 0));
	A[0] = a;
	for (int i = 1; i < N; i++)A[i][i - 1] = 1;

	int times = 1000;

	double t1 = GetTickCount();
	double x = power_method(A, u0, times);
	double t2 = GetTickCount();
	cout << "�����ʱ��" << ((t2 - t1)*1.0 / 1000) << "\n";
	cout << "����������" << times << "\n";
	cout << "���ģ����ֵ��" << x << "\n";
}

void exercise_2_1()
{
	cout << "Exercise_2_1" << endl << endl;
	int N = 41;
	vector<double> a(N, 0);
	a[N - 4] = -1; a[N - 1] = -1;
	vector<double>u0(N, 1);
	unsigned count;

	vector<vector<double>> A(N, vector<double>(N, 0));
	A[0] = a;
	for (int i = 1; i < N; i++)A[i][i - 1] = 1;
	vector<complex<double>> z;

	count = eigenvalues(A, z);

	cout << "����������" << count << "\n";
	cout << "����ֵ���£�" << "\n";
	for (int i = 0; i < N; i++)
	{
		cout << z[i] << "\n";
	}
}

void exercise_2_2()
{
	cout << "Exercise_2_2" << endl << endl;
	int N = 4;
	double x = 0.9;
	unsigned count;

	while (x <= 1.1)
	{
		cout << "x=" << x << "ʱ" << endl;
		vector<vector<double>> A = { {9.1,3.0,2.6,4.0},{4.2,5.3,4.7,1.6},{3.2,1.7,9.4,x},{6.1,4.9,3.5,6.2} };
		vector<complex<double>> z;

		count = eigenvalues(A, z);

		cout << "����������" << count << "\n";
		cout << "����ֵ���£�" << "\n";
		for (int i = 0; i < N; i++)
		{
			cout << z[i] << "\t";
		}
		cout << endl << endl;
		x += 0.1;
	}
}
