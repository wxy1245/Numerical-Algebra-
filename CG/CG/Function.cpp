#include "Function.h"
#include "auxiliary.h"
#include<cmath>

unsigned Jacobi(vector<vector<double>> A, vector<double> b, vector<double> &x, double epsilon)
{//Jacobi�����������x�ǳ�ֵ�����մ��ص�x�ǲ����Ľ��ƽ�
	int n = x.size();
	vector<double> diagonal(n);
	vector<double> x_pre(n);
	unsigned count = 0;

	//������������� B = D^(-1) (L+U) = I - D^(-1) A ,�ͱ��浽A����;�Լ����������� g = D^(-1) b ,�ͱ��浽b����
	for (int i = 0; i < n; i++) diagonal[i] = A[i][i];	//��һ�����������������D
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			A[i][j] = -(1.0 / diagonal[i])*A[i][j];
		}
		A[i][i] += 1;
	}
	for (int i = 0; i < n; i++) b[i] = (1.0 / diagonal[i])*b[i];

	//��ε�����x(k) = B x(k-1) + g; ע��Ҫ����ǰһ��x��ֵ��Ϊ������ܹ���ֹ����
	while (1)
	{
		x_pre = x;
		for (int i = 0; i < n; i++)
		{
			x[i] = b[i];
			for (int j = 0; j < n; j++)
			{
				x[i] += A[i][j] * x_pre[j];
			}
		}
		count += 1;
		if (vector_norm_2(vector_minus(x, x_pre)) < 1e-6) break;	//���µõ���x��ǰһ����x_pre����С��ʱ���������
	}

	return count;
}

unsigned GS(vector<vector<double>> A, vector<double> b, vector<double> &x, double epsilon)
{//Gauss-Seidel�����������x�ǳ�ֵ�����մ��ص�x�ǲ����Ľ��ƽ�
	int n = x.size();
	vector<double> diagonal(n);
	vector<double> x_pre(n);
	unsigned count = 0;

	for (int i = 0; i < n; i++) diagonal[i] = A[i][i];	//��һ�����������������D
	for (int i = 0; i < n; i++) b[i] = (1.0 / diagonal[i])*b[i];	//b�����������g

	//��ε�����x(k) = D^(-1) L x(k) + D^(-1) U x(k-1) + g;
	while (1)
	{
		x_pre = x;
		for (int i = 0; i < n; i++)	//����x�ĵ�i������ʱ���õ����Ѿ��������x��ֵǰi-1����������i+1��n��������ǰһ��x��ֵ
		{
			x[i] = b[i];
			for (int j = 0; j <= i - 1; j++)
			{
				x[i] += (1.0 / diagonal[i])*(-A[i][j])*x[j];
			}
			for (int j = i + 1; j < n; j++)
			{
				x[i] += (1.0 / diagonal[i])*(-A[i][j])*x_pre[j];
			}
		}
		count += 1;
		if (vector_norm_2(vector_minus(x, x_pre)) < 1e-6) break;	//���µõ���x��ǰһ����x_pre����С��ʱ���������
	}

	return count;
}

unsigned SOR(vector<vector<double>> A, vector<double> b, double omega, vector<double> &x, double epsilon)
{
	cout << "�ɳ����ӣ�" << omega << "\n";
	int n = x.size();
	vector<double> diagonal(n);
	vector<double> x_pre(n);
	unsigned count = 0;	//��¼��������

	for (int i = 0; i < n; i++) diagonal[i] = A[i][i];	//��һ�����������������D
	for (int i = 0; i < n; i++) b[i] = (1.0 / diagonal[i])*b[i];	//b�����������g

	//��ε�����x(k) = ��(D^(-1) L x(k) + D^(-1) U x(k-1) + g) + (1-��)x(k-1);
	while (1)
	{
		x_pre = x;
		for (int i = 0; i < n; i++)
		{
			//�������GS���������õ��ķ���ֵ
			x[i] = b[i];
			for (int j = 0; j <= i - 1; j++)
			{
				x[i] += (1.0 / diagonal[i])*(-A[i][j])*x[j];
			}
			for (int j = i + 1; j < n; j++)
			{
				x[i] += (1.0 / diagonal[i])*(-A[i][j])*x_pre[j];
			}
			//Ȼ����������������ɳڵ����µķ���ֵ
			x[i] = omega * x[i] + (1 - omega)*x_pre[i];
		}
		count += 1;
		if (vector_norm_2(vector_minus(x, x_pre)) < epsilon) break;	//���µõ���x��ǰһ����x_pre����С��ʱ���������
	}

	return count;
}

unsigned CG(vector<vector<double>> A, vector<double> b, vector<double> &x, double epsilon)
{//�����ݶȷ������Գ����������飻�ٶ�x�Ѿ����г�ֵ
	int k = 0, kmax = 1000;	//��һ����������������
	vector<double> r, p, w, x_pre;
	double rho, beta, alpha, rho_pre;
	r = vector_minus(b, matrix_vector(A, x));
	rho = inner_product(r, r);

	while (k < kmax)
	{
		x_pre = x;	//�ȼ�¼ǰһ��x�����������ֹ����
		k = k + 1;
		if (k == 1)
		{
			p = r;
		}
		else
		{
			beta = rho / rho_pre;
			p = vector_plus(r, scale_vector(beta, p));
		}
		w = matrix_vector(A, p);
		alpha = rho / inner_product(p, w);
		x = vector_plus(x, scale_vector(alpha, p));
		if (vector_norm_infinity(vector_minus(x, x_pre)) < epsilon)
		{
			if (vector_norm_infinity(p) < epsilon) k -= 1;	//�����ɽ����̫С�Ļ�����ʵ���һ�ξ�û�и��£����������������1
			break;
		}
		r = vector_minus(r, scale_vector(alpha, w));
		rho_pre = rho;
		rho = inner_product(r, r);
	}

	return k;	//���ص�������(k)
}

unsigned CGTest(vector<vector<double>> A, vector<double> b, vector<double> &x, double epsilon)
{//�����ݶȷ������Գ����������飻�ٶ�x�Ѿ����г�ֵ
	int k = 0, kmax = 1000;	//��һ����������������
	int n = x.size();
	vector<double> r, p, w, x_pre;
	double rho, beta, alpha, rho_pre;
	r = vector_minus(b, matrix_vector(A, x));
	rho = inner_product(r, r);
	cout << "��ʼ��xֵ:";
	for (int i = 0; i < n; i++)
	{
		cout.precision(3);
		cout << x[i] << "\t";
	}
	cout << "\n";
	cout << "��ʼ��rֵ:";
	for (int i = 0; i < n; i++)
	{
		cout.precision(3);
		cout << r[i] << "\t";
	}
	cout << "\n\n";

	while (k < kmax)
	{
		x_pre = x;	//�ȼ�¼ǰһ��x�����������ֹ����
		k = k + 1;
		if (k == 1)
		{
			p = r;
		}
		else
		{
			beta = rho / rho_pre;
			p = vector_plus(r, scale_vector(beta, p));
		}
		w = matrix_vector(A, p);
		alpha = rho / inner_product(p, w);
		x = vector_plus(x, scale_vector(alpha, p));
		if (vector_norm_infinity(vector_minus(x, x_pre)) < epsilon)
		{
			if (vector_norm_infinity(p) < epsilon) k -= 1;	//�����ɽ����̫С�Ļ�����ʵ���һ�ξ�û�и��£����������������1
			else
			{
				cout << "���һ�β����Ľ����" << "\n";
				cout << "��ɽ����p:";
				for (int i = 0; i < n; i++)
				{
					cout.precision(3);
					cout << p[i] << "\t";
				}
				cout << "\n";
				cout << "���º��x:";
				for (int i = 0; i < n; i++)
				{
					cout.precision(3);
					cout << x[i] << "\t";
				}
				cout << "\n";
			}
			cout << "\n";
			break;
		}
		r = vector_minus(r, scale_vector(alpha, w));
		rho_pre = rho;
		rho = inner_product(r, r);

		cout << "��" << k << "�ε���������Ľ��:" << "\n";
		cout << "��ɽ����p:";
		for (int i = 0; i < n; i++)
		{
			cout.precision(3);
			cout << p[i] << "\t";
		}
		cout << "\n";
		cout << "���º��x:";
		for (int i = 0; i < n; i++)
		{
			cout.precision(3);
			cout << x[i] << "\t";
		}
		cout << "\n";
		cout << "�¸�����r:";
		for (int i = 0; i < n; i++)
		{
			cout.precision(3);
			cout << r[i] << "\t";
		}
		cout << "\n";	
	}

	return k;	//���ص�������(k)
}