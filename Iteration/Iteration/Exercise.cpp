#include"Exercise.h"
#include<cmath>

void exercise_1()
{
	double epsilon = 0.0001, a = 0.5;	//epsilon��ȡ��ֵΪ1��0.1��0.01, 0.0001��
	int n = 100, N = n - 1;
	double h = 1.0 / n, x;
	unsigned count;
	vector<vector<double>> A(N, vector<double>(N));
	vector<double> b(N);
	vector<double> y(N, 0);	//��ĳ�ֵ��Ϊ��
	vector<double> sol(N);	//��ȷ��

	//�Ծ���A���Ҷ���b��ʼ��
	for (int i = 0; i < N - 1; i++)
	{
		A[i][i] = -(2 * epsilon + h);
		A[i][i + 1] = epsilon + h;
		A[i + 1][i] = epsilon;
		b[i] = a * h*h;
	}
	A[N - 1][N - 1] = -(2 * epsilon + h);
	b[N - 1] = a * h*h - h - epsilon;

	//�����ȷ�⣺sol=(1-a)(1-exp(-x/epsilon))/(1-exp(-1/epsilon))+ax
	for (int i = 0; i < N; i++)
	{
		x = (i + 1)*h;
		sol[i] = (1 - a)*(1 - exp(-x / epsilon)) / (1 - exp(-1 / epsilon)) + a * x;
	}

	//�ֱ�������ֲ�ͬ�����������������Ⲣ�Ҽ�¼����ʱ��
	double t1 = GetTickCount();
	//count = Jacobi(A, b, y);
	//count = GS(A, b, y);
	count = SOR(A, b, 1.1, y);	//������ĳ�������ѡȡ���ɳڵ�������
	double t2 = GetTickCount();
	cout << "����������" << count << "\n";
	cout << "�����ʱ��" << ((t2 - t1)*1.0 / 1000) << "\n";

	cout << "�����Ϊ:" << "\n";
	for (int i = 0; i < N; i++)
	{
		cout.precision(4);
		cout << y[i] << "\t";
		if (i % 10 == 9)cout << "\n";
	}
	cout << "\n";
	cout <<"�뾫ȷ�����Ϊ:"<< vector_norm_2(vector_minus(y, sol)) << "\n";

}

void exercise_2_1()
{//Jacobi�������������
	int n = 60, N = n + 1;	//�����ĵ���ΧһȦȫΪ1�ľ���
	vector<vector<double>> U(N, vector<double>(N,1));	//�����ĳ�ֵ��ѡ��1(����߽�����)
	vector<vector<double>> U_pre;
	double h = 1.0 / n;
	double err;	//ÿ�ε������¹����������Ƿ��Ѿ�������ܷ�Χ
	unsigned count = 0;	//��������

	double t1 = GetTickCount();
	while (1)
	{
		U_pre = U;	//�ȱ���ǰһ�εĽ��

		for (int i = 1; i < n; i++)	//������µĽ����
		{
			for (int j = 1; j < n; j++)
			{
				U[i][j] = (1.0 / (4 + h * h * exp(i*h*j*h))) * (h*h*(i*h + j * h) + U_pre[i - 1][j] + U_pre[i][j - 1] + U_pre[i][j + 1] + U_pre[i + 1][j]);
			}
		}

		count += 1;	//����������һ

		err = 0;	//������ڵ���������(�������ֱ�󿴳�������)���ж��Ƿ�Ӧ���˳�
		for (int i = 1; i < n; i++)
		{
			for (int j = 1; j < n; j++)
			{
				err += (U[i][j] - U_pre[i][j])*(U[i][j] - U_pre[i][j]);
			}
		}
		err = sqrt(err);
		if (err < 1e-7) break;
	}
	double t2 = GetTickCount();
	cout << "����������" << count << "\n";
	cout << "�����ʱ��" << ((t2 - t1)*1.0 / 1000) << "\n";

	cout << "�����Ϊ" << "\n";
	for (int i = 1; i < n; i++)
	{
		for (int j = 1; j < n; j++)
		{
			cout.precision(4);
			cout << U[i][j] << "\t";
		}
		cout << "\n";
	}
}

void exercise_2_2()
{//Gauss-Seidel�������������
	int n = 60, N = n + 1;	//�����ĵ���ΧһȦȫΪ1�ľ���
	vector<vector<double>> U(N, vector<double>(N, 1));	//�����ĳ�ֵ��ѡ��1(����߽�����)
	vector<vector<double>> U_pre;
	double h = 1.0 / n;
	double err;	//ÿ�ε������¹����������Ƿ��Ѿ�������ܷ�Χ
	unsigned count = 0;	//��������

	double t1 = GetTickCount();
	while (1)
	{
		U_pre = U;	//�ȱ���ǰһ�εĽ��

		for (int i = 1; i < n; i++)	//������µĽ����
		{
			for (int j = 1; j < n; j++)
			{
				U[i][j] = (1.0 / (4 + h * h * exp(i*h*j*h))) * (h*h*(i*h + j * h) + U[i - 1][j] + U[i][j - 1] + U_pre[i][j + 1] + U_pre[i + 1][j]);		//ֻ�ǰ�Jacobi�����޸ļ���:����������Ѿ���������µĽ�Ϳ�����
			}
		}

		count += 1;	//����������һ

		err = 0;	//������ڵ���������(�������ֱ�󿴳�������)���ж��Ƿ�Ӧ���˳�
		for (int i = 1; i < n; i++)
		{
			for (int j = 1; j < n; j++)
			{
				err += (U[i][j] - U_pre[i][j])*(U[i][j] - U_pre[i][j]);
			}
		}
		err = sqrt(err);
		if (err < 1e-7) break;
	}
	double t2 = GetTickCount();
	cout << "����������" << count << "\n";
	cout << "�����ʱ��" << ((t2 - t1)*1.0 / 1000) << "\n";

	cout << "�����Ϊ" << "\n";
	for (int i = 1; i < n; i++)
	{
		for (int j = 1; j < n; j++)
		{
			cout.precision(4);
			cout << U[i][j] << "\t";
		}
		cout << "\n";
	}
}

void exercise_2_3()
{
	int n = 60, N = n + 1;	//�����ĵ���ΧһȦȫΪ1�ľ���
	vector<vector<double>> U(N, vector<double>(N, 1));	//�����ĳ�ֵ��ѡ��1(����߽�����)
	vector<vector<double>> U_pre;
	double h = 1.0 / n;
	double err;	//ÿ�ε������¹����������Ƿ��Ѿ�������ܷ�Χ
	unsigned count = 0;	//��������
	double omega = 1.90;	//�ɳڵ�������(20��ʱ��1.75, 40��ʱ:1.855, 60��ʱ:1.90)

	double t1 = GetTickCount();
	while (1)
	{
		U_pre = U;	//�ȱ���ǰһ�εĽ��

		for (int i = 1; i < n; i++)	//������µĽ����
		{
			for (int j = 1; j < n; j++)
			{
				U[i][j] = (omega / (4 + h * h * exp(i*h*j*h)))*(h*h*(i*h + j * h) + U[i - 1][j] + U[i][j - 1] + U_pre[i][j + 1] + U_pre[i + 1][j]) - (omega - 1)*U_pre[i][j];	//�ɳڵ�����ʽ
			}
		}

		count += 1;	//����������һ

		err = 0;	//������ڵ���������(�������ֱ�󿴳�������)���ж��Ƿ�Ӧ���˳�
		for (int i = 1; i < n; i++)
		{
			for (int j = 1; j < n; j++)
			{
				err += (U[i][j] - U_pre[i][j])*(U[i][j] - U_pre[i][j]);
			}
		}
		err = sqrt(err);
		if (err < 1e-7) break;
	}
	double t2 = GetTickCount();
	cout << "�ɳ����ӣ�" << omega << "\n";
	cout << "����������" << count << "\n";
	cout << "�����ʱ��" << ((t2 - t1)*1.0 / 1000) << "\n";

	cout << "�����Ϊ" << "\n";
	for (int i = 1; i < n; i++)
	{
		for (int j = 1; j < n; j++)
		{
			cout.precision(4);
			cout << U[i][j] << "\t";
		}
		cout << "\n";
	}
}