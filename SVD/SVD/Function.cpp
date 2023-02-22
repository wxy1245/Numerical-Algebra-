#include "Function.h"
#include "auxiliary.h"
#include <cmath>
#define ERROR 1e-6

void householder(vector<double> x, vector<double> &v, double &beta)
{
	int n = x.size();
	double eta = vector_norm_infinity(x);
	if (abs(eta) < 1e-5)
	{
		for (int i = 0; i < n; i++)
		{
			v[i] = 0;
		}
		beta = 0;
		return;
	}

	for (int i = 0; i < n; i++)
	{
		x[i] /= eta;
	}
	double sigma = 0;
	for (int i = 1; i < n; i++)
	{
		sigma += x[i] * x[i];
		v[i] = x[i];
	}

	if (sigma < 1.0e-10)
	{
		beta = 0;
	}
	else
	{
		double alpha = sqrt(x[0] * x[0] + sigma);
		if (x[0] <= 0)
		{
			v[0] = x[0] - alpha;
		}
		else
		{
			v[0] = -sigma / (x[0] + alpha);
		}
		beta = 2 * v[0] * v[0] / (sigma + v[0] * v[0]);
		for (int i = 1; i < n; i++)
		{
			v[i] /= v[0];
		}
		v[0] = 1;	//!!!ע�⣺v���й��ʱ������v[0]�������Ԫ�أ����Ű�v[0]��1
	}
}

void matrix_householder_rowsimplify_p_q(vector<vector<double>> &A, vector<double> v, double beta, int p, int q)
{//��A��p��q���Ժ�Ĳ��ֽ���һ��householder���� A(m-p)*(n-q) v:m-p
	int m = A.size();	//������
	int n = A[0].size();	//������
	vector<double> w(n - q, 0);	//w������ֻ��Ҫn-q

	//�ȼ���w=beta*A^t v
	for (int i = q; i < n; i++)
	{
		for (int k = p; k < m; k++)
		{
			w[i - q] += A[k][i] * v[k - p];
		}
		w[i - q] *= beta;
	}

	//H A = A - v w^T,ֱ�Ӿͱ�����A����
	for (int k = p; k < m; k++)
	{
		for (int i = q; i < n; i++)
		{
			A[k][i] -= v[k - p] * w[i - q];
		}
	}
}

void matrix_householder_colsimplify_p_q(vector<vector<double>> &A, vector<double> v, double beta, int p, int q)
{//��A��p��q���Ժ�Ĳ��ֽ���һ��householder���� A(m-p)*(n-q) v:n-q
	int m = A.size();	//������
	int n = A[0].size();	//������
	vector<double> w(m - p, 0);	//w������ֻ��Ҫm-p

	for (int i = p; i < m; i++)
	{
		for (int k = q; k < n; k++)
		{
			w[i - p] += A[i][k] * v[k - q];
		}
		w[i - p] *= beta;
	}

	for (int i = p; i < m; i++)
	{
		for (int k = q; k < n; k++)
		{
			A[i][k] -= w[i - p] * v[k - q];
		}
	}
}

void matrix_householder_rowsimplify_p1q1p2q2(vector<vector<double>> &A, vector<double> v, double beta, int p1, int q1, int p2, int q2)
{//��A��p1����p2�У�q1����q2�н�����householder�任����
	int m = p2;	//������
	int n = q2;	//������
	vector<double> w(n - q1, 0);	//w������ֻ��Ҫn-q1

	//�ȼ���w=beta*A^t v
	for (int i = q1; i < n; i++)
	{
		for (int k = p1; k < m; k++)
		{
			w[i - q1] += A[k][i] * v[k - p1];
		}
		w[i - q1] *= beta;
	}

	//H A = A - v w^T,ֱ�Ӿͱ�����A����
	for (int k = p1; k < m; k++)
	{
		for (int i = q1; i < n; i++)
		{
			A[k][i] -= v[k - p1] * w[i - q1];
		}
	}
}

void matrix_householder_colsimplify_p1q1p2q2(vector<vector<double>> &A, vector<double> v, double beta, int p1, int q1, int p2, int q2)
{//��A��p1����p2�У�q1����q2�н�����householder�任����
	int m = p2;	//������
	int n = q2;	//������
	vector<double> w(m - p1, 0);	//w������ֻ��Ҫm-p1

	//�ȼ���w=beta*A v
	for (int i = p1; i < m; i++)
	{
		for (int k = q1; k < n; k++)
		{
			w[i - p1] += A[i][k] * v[k - q1];
		}
		w[i - p1] *= beta;
	}

	//H A = A - v w^T,ֱ�Ӿͱ�����A����
	for (int i = p1; i < m; i++)
	{
		for (int k = q1; k < n; k++)
		{
			A[i][k] -= w[i - p1] * v[k - q1];
		}
	}
}

void two_diag(vector<vector<double>> &A, vector<vector<double>> &U, vector<vector<double>> &V)
{//ע��UΪm����������,VΪn���������󣬳�ʼ���߾�Ϊ��λ����
	int m = A.size();
	int n = A[0].size();
	double beta;
	for (int k = 0; k < n; k++)	//������ע���������n��
	{
		vector<double> x1(m - k);
		vector<double> v1(m - k), u1(n - k, 0);
		for (int i = k; i < m; i++)
		{
			x1[i - k] = A[i][k];
		}
		householder(x1, v1, beta);
		for (int i = 0; i < n - k; i++)
		{
			for (int j = 0; j < m - k; j++)
			{
				u1[i] += A[j + k][i + k] * v1[j];
			}
			u1[i] *= beta;
		}
		for (int i = k; i < m; i++)
		{
			for (int j = k; j < n; j++)
			{
				A[i][j] -= v1[i - k] * u1[j - k];
			}
		}
		matrix_householder_colsimplify_p_q(U, v1, beta, 0, k);	//���������任����U:������ע��U�ӵ�k��֮����ж���任
		if (k < n - 1)
		{
			vector<double> x2(n - 1 - k);
			vector<double> v2(n - 1 - k), u2(m - k, 0);
			for (int i = k + 1; i < n; i++)
			{
				x2[i - (k + 1)] = A[k][i];
			}
			householder(x2, v2, beta);
			for (int i = 0; i < m - k; i++)
			{
				for (int j = 0; j < n - 1 - k; j++)
				{
					u2[i] += A[i + k][j + (k + 1)] * v2[j];	//������ע��A��Ԫ�ص��б�Ϊj+(k+1)
				}
				u2[i] *= beta;
			}
			for (int i = k; i < m; i++)
			{
				for (int j = k + 1; j < n; j++)
				{
					A[i][j] -= u2[i - k] * v2[j - (k + 1)];
				}
			}
			matrix_householder_colsimplify_p_q(V, v2, beta, 0, k + 1);	//���������任����V:�������ӵ�k+1��֮�󣬸��ж��任
		}
	}
}

void Givens_row(vector<vector<double>>& A, int p, int q, double c, double s)
{//ֻ��p,q�з����仯
	int n = A[0].size();
	double a_pi, a_qi;
	for (int i = 0; i < n; i++)
	{
		a_pi = A[p][i]; a_qi = A[q][i];
		A[p][i] = c * a_pi + s * a_qi;
		A[q][i] = -s * a_pi + c * a_qi;
	}
}

void Givens_column(vector<vector<double>>& A, int p, int q, double c, double s)
{//ֻ��p,q�з����仯
	int n = A.size();
	double a_ip, a_iq;
	for (int i = 0; i < n; i++)
	{
		a_ip = A[i][p]; a_iq = A[i][q];
		A[i][p] = c * a_ip - s * a_iq;
		A[i][q] = s * a_ip + c * a_iq;
	}
}

void Wilkinson_step(vector<vector<double>>& B, vector<vector<double>>& U, vector<vector<double>>& V)
{
	int n = B.size();
	double alpha = B[n - 1][n - 1] * B[n - 1][n - 1] + B[n - 2][n - 1] * B[n - 2][n - 1];
	double delta = (B[n - 2][n - 2] * B[n - 2][n - 2] + (n > 2 ? B[n - 3][n - 2] * B[n - 3][n - 2] : 0) - alpha) / 2.0;	//��ע��n==2�����
	double beta = B[n - 2][n - 2] * B[n - 2][n - 1];
	double mu = alpha - beta * beta / (delta + (delta > 0 ? 1 : -1)*sqrt(delta*delta + beta * beta));
	double y = B[0][0] * B[0][0] - mu;
	double z = B[0][0] * B[0][1];
	double c, s;
	c = -y / sqrt(y * y + z * z); s = z / sqrt(y * y + z * z);
	Givens_column(B, 0, 1, c, s); //cout << endl; matrix_print(B);
	Givens_column(V, 0, 1, c, s); //cout << endl; matrix_print(V);
	//cout << endl;
	//matrix_print(matrix_product(matrix_product(U, B), transpose(V)));

	for (int k = 0; k < n - 1; k++)
	{
		y = B[k][k]; z = B[k + 1][k];
		c = y / sqrt(y * y + z * z); s = z / sqrt(y * y + z * z);
		Givens_row(B, k, k + 1, c, s); 
		Givens_column(U, k, k + 1, c, -s); //����U��������ע�⣬�ҳ�ʱ���۵���ת�þ����൱����ԭGivens����s��Ϊ-s

		if (k < n - 2)
		{
			y = B[k][k + 1]; z = B[k][k + 2];
			c = -y / sqrt(y * y + z * z); s = z / sqrt(y * y + z * z);
			Givens_column(B, k + 1, k + 2, c, s);
			Givens_column(V, k + 1, k + 2, c, s);	//����V
		}
	}
}

unsigned SVD(vector<vector<double>>& A, vector<vector<double>>& U, vector<vector<double>>& V)
{
	int m = A.size();
	int n = A[0].size();
	int q, p, pos;
	vector<vector<double>> B, B22, P, Q, Pnew, Qnew;
	double y, z, c, s;
	unsigned count = 0;

	//1.���Խǻ�
	two_diag(A, U, V);
	
	for (int i = n; i < m; i++)	//n���Ժ�Aȫ������
	{
		for (int j = 0; j < n; j++)
		{
			A[i][j] = 0;	
		}
	}

	B = sub_matrix(A, 0, n, 0, n);

	//2.�������ж�
	
	for (int i = 0; i < n - 2; i++)
	{
		for (int j = i + 2; j < n; j++)
		{
			B[i][j] = B[j][i] = 0;
		}
	}
	
	for (int i = 0; i < n - 1; i++)
	{
		B[i + 1][i] = 0;
		if (abs(B[i][i + 1]) <= (abs(B[i][i]) + abs(B[i + 1][i + 1]))*ERROR)B[i][i + 1] = 0;
	}
	for (int i = 0; i < n; i++)
	{
		if (abs(B[i][i]) <= matrix_norm_infinity(B)*ERROR)B[i][i] = 0;
	}
	//��Ҫȷ������q����С��p(���¶���)
	q = 0;
	while (q < n)
	{
		if (q == n - 1 || B[n - 2 - q][n - 1 - q] == 0)q++;
		else break;
	}
	p = n - q;
	if (p > 0)	//p����Ϊ2
	{
		p -= 2;
		while (p > 0)
		{
			if (B[p - 1][p] != 0)p--;
			else break;
		}
	}

	while (q < n)	//SVD����,q=nʱ�˳�;��������һ�ε���
	{
		B22 = sub_matrix(B, p, n - q, p, n - q);

		for (pos = n - q - p - 2; pos >= 0; pos--)
		{
			if (B22[pos][pos] == 0)break;
		}

		if (pos >= 0)
		{	
			P = id_matrix(n - q - p);
			for (int k = pos + 1; k < n - q - p; k++)
			{
				y = B22[pos][k]; z = B22[k][k];
				c = -z / sqrt(y * y + z * z); s = y / sqrt(y * y + z * z);
				Givens_row(B22, pos, k, c, s);
				Givens_column(P, pos, k, c, -s);
			}
			matrix_part_replace(B, p, n - q, p, n - q, B22);
			Pnew = id_matrix(m);
			matrix_part_replace(Pnew, p, n - q, p, n - q, P);
			U = matrix_product(U, Pnew);
		}
		else
		{
			P = id_matrix(n - q - p); Q = id_matrix(n - q - p);
			Wilkinson_step(B22, P, Q);
			matrix_part_replace(B, p, n - q, p, n - q, B22);
			Pnew = id_matrix(m); Qnew = id_matrix(n);
			matrix_part_replace(Pnew, p, n - q, p, n - q, P); matrix_part_replace(Qnew, p, n - q, p, n - q, Q);
			U = matrix_product(U, Pnew); V = matrix_product(V, Qnew);	//����U��V
		}

		//�ص��������ж�
		
		for (int i = 0; i < n - 2; i++)
		{
			for (int j = i + 2; j < n; j++)
			{
				B[i][j] = B[j][i] = 0;
			}
		}
		
		for (int i = 0; i < n - 1; i++)
		{
			B[i + 1][i] = 0;
			if (abs(B[i][i + 1]) <= (abs(B[i][i]) + abs(B[i + 1][i + 1]))*ERROR)B[i][i + 1] = 0;
		}
		for (int i = 0; i < n; i++)
		{
			if (abs(B[i][i]) <= matrix_norm_infinity(B)*ERROR)B[i][i] = 0;
		}
		//��Ҫȷ������q����С��p(���¶���)
		q = 0;
		while (q < n)
		{
			if (q == n - 1 || B[n - 2 - q][n - 1 - q] == 0)q++;
			else break;
		}
		p = n - q;
		if (p > 0)	//p����Ϊ2
		{
			p -= 2;
			while (p > 0)
			{
				if (B[p - 1][p] != 0)p--;
				else break;
			}
		}
		count += 1;
	}
	
	matrix_part_replace(A, 0, n, 0, n, B);	//�������Bȥ�滻A����Ӧ�Ĳ��ּ���
	return count;
}

vector<double> singular_values_sort(vector<vector<double>> D)
{
	int n = D[0].size();	//����
	vector<double> v;
	for (int j = 0; j < n; j++)
	{
		if (D[j][j] != 0)v.push_back(D[j][j]);
	}
	sort(v.begin(), v.end());
	return v;
}
