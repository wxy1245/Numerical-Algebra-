#include"Exercise.h"
#include"auxiliary.h"

//!!!ע�⣺����2.1,3.1����Ҫ��Ax-b������AҲ��Ҫ��������
//!!!�Ѿ���������Ƶĳ���(��GaussSolve)���Լ������ڲ���Ӳ���Ҫ���ܡ�

void exercise_1()
{
	int N = 84; //�����С
	vector<vector<double>> A(N, vector<double>(N, 0));		//����һ��N��N�е����飬Ԫ��������double
	vector<double> b(N);	//����һ��Nά��������Ԫ������Ҳ��double	
	vector<double> sol(N, 1);	//��ʵ��

	Ex1_Init(A, b, N);//��ʼ��

	double t1 = GetTickCount();
	//GaussSolve(A, b);
	//FullGaussSolve(A, b);
	ColGaussSolve(A, b);
	double t2 = GetTickCount();
	cout << "�����ʱ��" << ((t2 - t1)*1.0 / 1000) << "\n";

	cout << "������Ľ�����:" << "\n";
	for (int i = 0; i < N; i++)
	{
		cout.precision(3);
		cout << b[i] << "\t";
		if (i % 10 == 9)cout << "\n";
	}
	cout << "\n";

	cout << "����ʵ��֮��:" << ResultError(b, sol);
}

void exercise_2_1()
{
	int N = 100; 
	vector<vector<double>> A(N, vector<double>(N, 0)), A0(N, vector<double>(N, 0));
	vector<double> b(N), b0(N);
	vector<double> y0(N);	

	Ex2_Init1(A, b, N);
	A0 = A; b0 = b;	//��¼ԭʼ��A��b��������Ax-b�����(ע��:vector�����ǿ������帳ֵ��!)
	
	double t1 = GetTickCount();
	//CholeskySolve(A, b);	//ƽ���������
	ModifiedCholeskySolve(A, b);	//�Ľ���ƽ���������
	double t2 = GetTickCount();
	cout << "�����ʱ��" << ((t2 - t1)*1.0 / 1000) << "\n";

	cout << "������Ľ�����:" << "\n";
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
			y0[i] += A0[i][j] * b[j];			//����y=Ax(x�Ѵ����b��),�������y0��
		}
	}
	cout<< "Ax-b���:" << ResultError(y0, b0) << "\n";
}

void exercise_2_2()
{
	int N = 40;
	vector<vector<double>> A(N, vector<double>(N, 0));
	vector<double> b(N,0);
	vector<double> sol(N, 1);	//��ʵ��

	Ex2_Init2(A, b, N);

	double t1 = GetTickCount();
	//CholeskySolve(A, b);
	ModifiedCholeskySolve(A, b);
	double t2 = GetTickCount();
	cout << "�����ʱ��" << ((t2 - t1)*1.0 / 1000) << "\n";

	cout << "������Ľ�����:" << "\n";
	for (int i = 0; i < N; i++)
	{
		cout.precision(3);
		cout << b[i] << "\t";
	}
	cout << "\n";

	cout << "����ʵ��֮��:" << ResultError(b, sol);
}

void exercise_3_1()
{//�õ�һ��ķ�ʽ���ڶ���ĵ�һ��������
	int N = 100;
	vector<vector<double>> A(N, vector<double>(N, 0)), A0(N, vector<double>(N, 0));
	vector<double> b(N), b0(N);
	vector<double> y0(N);
	
	Ex2_Init1(A, b, N);
	A0 = A; b0 = b;	//��¼ԭʼ��A��b��������Ax-b�����

	double t1 = GetTickCount();
	//GaussSolve(A, b);
	//FullGaussSolve(A, b);
	ColGaussSolve(A, b);
	double t2 = GetTickCount();
	cout << "�����ʱ��" << ((t2 - t1)*1.0 / 1000) << "\n";

	cout << "������Ľ�����:" << "\n";
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
			y0[i] += A0[i][j] * b[j];			//����y=Ax(x�Ѵ����b��),�������y0��
		}
	}
	cout << "Ax-b���:" << ResultError(y0, b0) << "\n";
}

void exercise_3_2()
{//�õ�һ��ķ�ʽ���ڶ���ĵڶ���������
	int N = 15;
	vector<vector<double>> A(N, vector<double>(N, 0));
	vector<double> b(N);
	vector<double> sol(N, 1);	//��ʵ��

	Ex2_Init2(A, b, N);
	double t1 = GetTickCount();
	//GaussSolve(A, b);
	//FullGaussSolve(A, b);
	ColGaussSolve(A, b);
	double t2 = GetTickCount();
	cout << "�����ʱ��" << ((t2 - t1)*1.0 / 1000) << "\n";

	cout << "������Ľ�����:" << "\n";
	for (int i = 0; i < N; i++)
	{
		cout.precision(5);
		cout << b[i] << "\t";
		if (i % 10 == 9)cout << "\n";
	}
	cout << "\n";

	cout << "����ʵ��֮��:" << ResultError(b, sol);
}