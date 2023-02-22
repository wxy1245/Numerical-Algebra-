#include"Exercise.h"
#include"auxiliary.h"

void ColGaussSolve(vector<vector<double>>& A, vector<double>& b)
{//����ԪGauss��ȥ��������Է�����
	int N = A.size();
	vector<int> u;
	gauss_elim_col_pivoting(A, u);	//PA=LU��ɷֽ�֮���Ⱥ����Ly=Pb,Ux=y
	vector_pb(u, b);
	forward_subs1(A, b);
	back_subs(A, b);
}

double Col_Precision(vector<vector<double>>& A, vector<double>& b)
{
	double mu = infinity_norm(A);
	double beta = vector_norm_infinity(b);
	double gamma, nu;
	vector<vector<double>>A0 = A;	//r�������Ҫ����ԭ���ľ���������
	vector<double>b0 = b, r;

	ColGaussSolve(A, b);	//A�ŵ������Ƿֽ⣬b�ŵ��Ǽ����
	r = vector_minus(b0, matrix_vector(A0, b));
	gamma = vector_norm_infinity(r);
	nu = inverse_infinity_norm(A, b);
	return nu * mu * gamma / beta;
}

void exercise_1()
{
	int N;	//ϣ�����ؾ������
	for (N = 2; N <= 20; N++)
	{
		vector<vector<double>> A(N, vector<double>(N, 0));
		vector<double> b(N, 0);

		for (int i = 0; i < N; i++)
		{
			for (int j = 0; j < N; j++)
			{
				A[i][j] = 1.0 / (i + j + 1);
				b[i] += A[i][j];
			}
		}

		double mu = infinity_norm(A);	//!!!ע��A��������õ�ԭ����ϣ�����ؾ������󣬶��������������õ�����Ԫ��ȥ�õ��ľ���
		ColGaussSolve(A, b);
		cout << N << "��Hilbert����������������Ϊ:" << mu * inverse_infinity_norm(A, b) << "\n";
	}
}

void exercise_2()
{
	int N;	//����A�Ľ���
	for (N = 5; N <= 30; N++)
	{
		vector<vector<double>> A(N, vector<double>(N, 1));
		vector<double> x(N), b(N, 0);
		double accuracy, estimation;

		for (int i = 1; i < N - 1; i++)	//���ɾ���A
		{
			for (int j = 0; j < i; j++)
			{
				A[i][j] = -1; A[j][i] = 0;
			}
		}
		for (int j = 0; j < N - 1; j++)
		{
			A[N - 1][j] = -1;
		}
		random_value(x);	//���ѡȡһ��x,��Ҳ�Ǻ������ʵ��

		for (int i = 0; i < N; i++)	//��Ϊ֪������A���ص㣬���Լ򵥼���b
		{
			for (int j = 0; j < i; j++)
			{
				b[i] -= x[j];
			}
			b[i] += x[i] + x[N - 1];
		}
		b[N - 1] -= x[N - 1];

		estimation = Col_Precision(A, b);	//����ԪGauss��ȥ������ľ��ȹ���,ע��˲����A��b���Ѿ��ı�
		accuracy = vector_norm_infinity(vector_minus(x, b)) / vector_norm_infinity(x);
		cout.setf(ios::scientific);
		cout << N << "�׾���,���������Ϊ��" << accuracy << ", ��������Ϊ��" << estimation << "\n";
	}
}