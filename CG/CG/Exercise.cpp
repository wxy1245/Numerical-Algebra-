#include"Exercise.h"
#include"auxiliary.h"
#include<cmath>

void exercise_1()
{
	int n = 20, N = (n - 1)*(n - 1);	//N:ϵ������A�Ľ�����Ҳ���Ҷ���b�ͽ�x��ά��
	vector<vector<double>> A(N, vector<double>(N, 0));	
	vector<double> b(N);
	vector<double> x(N, 0);	//��ֵ�򵥾�ȡ0
	double h = 1.0 / n;
	unsigned count;
	double omega = 1.754;	//�ɳ�����


	//��ϵ������A���Ҷ���b���г�ʼ��(A�ֿ��ʼ��,�ܹ���n-1��;b����ʱҪע��߽�����)
	for (int k = 0; k < n - 2; k++)	//A�ĳ�ʼ��
	{
		for (int i = k * (n - 1); i < k * (n - 1) + (n - 2); i++)
		{
			A[i][i] = 1 + h * h / 4;
			A[i][i + 1] = A[i + 1][i] = -1.0 / 4;
			A[i][i + n - 1] = A[i + n - 1][i] = -1.0 / 4;
		}
		A[k * (n - 1) + (n - 2)][k * (n - 1) + (n - 2)] = 1 + h * h / 4;
		A[k * (n - 1) + (n - 2)][k * (n - 1) + (n - 2) + (n - 1)] = A[k * (n - 1) + (n - 2) + (n - 1)][k * (n - 1) + (n - 2)] = -1.0 / 4;
	}
	for (int i = (n - 2)*(n - 1); i < (n - 2)*(n - 1) + (n - 2); i++)
	{
		A[i][i] = 1 + h * h / 4;
		A[i][i + 1] = A[i + 1][i] = -1.0 / 4;
	}
	A[N - 1][N - 1] = 1 + h * h / 4;
	for (int m = 0; m < N; m++)	//b�ĳ�ʼ��
	{
		int i = 1 + m / (n - 1);
		int j = (m + 1) - (i - 1)*(n - 1);
		b[m] = (h*h / 4)*sin(i*h*j*h);
		if (i - 1 == 0)b[m] += (1.0 / 4)*(j * j*h*h);
		if (j - 1 == 0)b[m] += (1.0 / 4)*(i * i*h*h);
		if (i + 1 == n)b[m] += (1.0 / 4)*(1 + j * j*h*h);
		if (j + 1 == n)b[m] += (1.0 / 4)*(i * i*h*h + 1);
	}

	double t1 = GetTickCount();
	//count = CG(A, b, x, 1.0e-7);	//�����ݶȷ����
	count = SOR(A, b, omega, x, 1.0e-7);	//SOR�������
	double t2 = GetTickCount();
	
	cout << "��������Ϊ��" << count << "\n";
	cout << "�����ʱ��" << ((t2 - t1)*1.0 / 1000) << "\n";
	cout << "������ƽ����£�" << "\n";
	for (int i = 0; i < N; i++)
	{
		cout.precision(4);
		cout << x[i] << '\t';
		if ((i + 1) % (n - 1) == 0)cout << "\n";
	}
}

void exercise_2()
{
	int N = 80;
	vector<vector<double>> A(N, vector<double>(N, 0));
	vector<double> b(N);
	vector<double> x(N, 0);	//��ֵ�򵥾�ȡ0
	unsigned count;

	//��ʼ��������ϣ�����ؾ���A���Ҷ���b
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++)
		{
			A[i][j] = 1.0 / (i + j + 1);	//ע�⣺�������д��1.0��������������ֵ�����Ĭ�ϵõ���ֵ
			b[i] += A[i][j];
		}
		b[i] /= 3.0;
	}

	double t1 = GetTickCount();
	count = CG(A, b, x, 1.0e-6);	//�����ݶȷ����
	double t2 = GetTickCount();

	cout << "��������Ϊ��" << count << "\n";
	cout << "�����ʱ��" << ((t2 - t1)*1.0 / 1000) << "\n";
	cout << "������ƽ����£�" << "\n";
	for (int i = 0; i < N; i++)
	{
		cout.precision(4);
		cout << x[i] << '\t';
		if ((i + 1) % 10 == 0)cout << "\n";
	}
}

void exercise_3()
{
	int N = 5; //�������
	vector<vector<double>> A(N, vector<double>(N));
	vector<double> b(N);
	vector<double> x(N, 0);
	unsigned count;

	A[0][0] = 10; A[0][1] = 1; A[0][2] = 2; A[0][3] = 3; A[0][4] = 4; b[0] = 12;
	A[1][0] = 1; A[1][1] = 9; A[1][2] = -1; A[1][3] = 2; A[1][4] = -3; b[1] = -27;
	A[2][0] = 2; A[2][1] = -1; A[2][2] = 7; A[2][3] = 3; A[2][4] = -5; b[2] = 14;
	A[3][0] = 3; A[3][1] = 2; A[3][2] = 3; A[3][3] = 12; A[3][4] = -1; b[3] = -17;
	A[4][0] = 4; A[4][1] = -3; A[4][2] = -5; A[4][3] = -1; A[4][4] = 15; b[4] = 12;

	double t1 = GetTickCount();
	//count = Jacobi(A, b, x, 1.0e-9);
	//count = GS(A, b, x, 1.0e-9);
	count = CGTest(A, b, x, 1.0e-9);
	double t2 = GetTickCount();

	cout << "��������Ϊ��" << count << "\n";
	cout << "�����ʱ��" << ((t2 - t1)*1.0 / 1000) << "\n";
	cout << "������ƽ����£�" << "\n";
	for (int i = 0; i < N; i++)
	{
		cout.precision(4);
		cout << x[i] << '\t';
	}
}

/*
for (int i = 0; i < N; i++)
{
	for (int j = 0; j < N; j++)
	{
		cout.precision(3);
		cout << A[i][j] << "\t";
	}
	cout << "\n";
}
*/

/*
for (int i = 0; i < N; i++)
{
	cout.precision(3);
	cout << b[i] << "\t";
}
*/