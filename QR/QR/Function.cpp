#include "Function.h"
#include"auxiliary.h"
#include <cmath>
#define ERROR 1e-6

vector<double> QR_decomposition(vector<vector<double>> &A)
{
	int m = A.size();	//����
	int n = A[0].size();	//����
	vector<double> v, temp, d;	//temp:����ȡA�ĵ�j��(j���Ժ�)����householder�任��v,��;d:��Ÿ��α任�Ħ�
	double beta;
	for (int j = 0; j < n; j++)
	{
		if (j < m - 1)	//j=m-1ʱ���Ѿ������еף����������
		{
			v.resize(m - j); temp.resize(m - j);
			for (int k = j; k < m; k++)
			{
				temp[k - j] = A[k][j];
			}
			householder(temp, v, beta);
			matrix_householder_simplify_j(A, v, beta, j);
			d.push_back(beta);
			for (int k = j + 1; k < m; k++)
			{
				A[k][j] = v[k - j];
			}
		}
	}
	return d;
}

double power_method(vector<vector<double>> A, vector<double> &u0, int &times)
{//�ݷ������ģ����ֵ�����أ�����������Ӧ������������������u0�У�,times����ʱ��Ϊ�������������������ʵ��������
	vector<double> y = u0;
	double mu_pre, mu;
	int count = 2;
	
	y = matrix_vector(A, y);
	mu_pre = y[vector_norm_maxindex(y)];
	y = scale_vector(1.0 / mu_pre, y);

	y = matrix_vector(A, y);
	mu = y[vector_norm_maxindex(y)];
	y = scale_vector(1.0 / mu, y);

	while(abs(mu-mu_pre)>1e-5 && count<times)	//count < times
	{
		y = matrix_vector(A, y);
		mu_pre = mu;
		mu = y[vector_norm_maxindex(y)];
		y = scale_vector(1.0 / mu, y);
		count++;
	}

	times = count;
	u0 = y;
	return mu;
}

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

void matrix_householder_simplify(vector<vector<double>> &A, vector<double> v, double beta)
{
	int m = A.size();	//����
	int n = A[0].size();	//����
	vector<double> w(n, 0);

	//�ȼ���w=beta*A^t v
	for (int i = 0; i < n; i++)
	{
		for (int k = 0; k < m; k++)
		{
			w[i] += A[k][i] * v[k];
		}
		w[i] *= beta;
	}

	//H A = A - v w^T,ֱ�Ӿͱ�����A����
	for (int i = 0; i < m; i++)
	{
		for (int j = 0; j < n; j++)
		{
			A[i][j] -= v[i] * w[j];
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

void matrix_householder_simplify_j(vector<vector<double>> &A, vector<double> v, double beta, int j)
{//��A��j��j���Ժ�Ĳ��ֽ���һ��householder����
	int m = A.size();	//������
	int n = A[0].size();	//������
	vector<double> w(n - j, 0);	//w������ֻ��Ҫn-j

	//�ȼ���w=beta*A^t v
	for (int i = j; i < n; i++)
	{
		for (int k = j; k < m; k++)
		{
			w[i - j] += A[k][i] * v[k - j];
		}
		w[i - j] *= beta;
	}

	//H A = A - v w^T,ֱ�Ӿͱ�����A����
	for (int k = j; k < m; k++)
	{
		for (int i = j; i < n; i++)
		{
			A[k][i] -= v[k - j] * w[i - j];
		}
	}
}

void matrix_householder_rowsimplify_p_q(vector<vector<double>> &A, vector<double> v, double beta,int p, int q)
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

	//�ȼ���w=beta*A v
	for (int i = p; i < m; i++)
	{
		for (int k = q; k < n; k++)
		{
			w[i - p] += A[i][k] * v[k - q];
		}
		w[i - p] *= beta;
	}

	//H A = A - v w^T,ֱ�Ӿͱ�����A����
	for (int i = p; i < m; i++)
	{
		for (int k = q; k < n; k++)
		{
			A[i][k] -= w[i - p] * v[k - q];
		}
	}
}

void hessenberg_decomp(vector<vector<double>>&A, vector<vector<double>>&V, vector<double> &b)
{
	double beta;
	int n = A.size();
	for (int k = 0; k < n - 2; k++)
	{
		//����householder�任����(��k+1�п�ʼ)
		vector<double> x(n - k - 1), v(n - k - 1, 0);
		for (int i = k + 1; i < n; i++)
		{	
			x[i - k - 1] = A[i][k];
		}
		householder(x, v, beta);

		//��A���ֽ���household�任
		matrix_householder_rowsimplify_p_q(A, v, beta, k + 1, k);
		matrix_householder_colsimplify_p_q(A, v, beta, 0, k + 1);
		
		//��¼�任��v�Լ�beta
		V.push_back(v); b.push_back(beta);
	}

}

vector<vector<double>> two_step_displacement_qr(vector<vector<double>> &H)
{//ע���������������ΪQR������Ҫ��
	int n = H.size();
	vector<double> v(3);
	double beta, x, y, z;
	vector<vector<double>> P, Pnew;
	vector<vector<double>> Psmall;
	P = id_matrix(n);
	
	if (n == 1)return P;
	else if (n == 2)
	{
		//����n=2�������������õ�����������ʱ����ʵ����ֵ�����Կ��Ǵ�ԭ��λ�Ƶ�QR����
		double mu = H[1][1];
		//H- mu I =QR
		H[0][0] -= mu; H[1][1] -= mu;
		vector<double> beta = QR_decomposition(H);
		vector<double> v1 = { 1,H[1][0] };
		double beta1 = beta[0];
		vector<vector<double>>Q(n, vector<double>(n, 0));
		for (int i = 0; i < n; i++)
		{
			for (int j = 0; j < n; j++)
			{
				Q[i][j] -= beta1 * v1[i] * v1[j];
				if (i == j)Q[i][j] += 1;
			}
		}
		//H1=RQ + mu I
		vector<vector<double>>H1(n, vector<double>(n, 0));
		for (int i = 0; i < n; i++)
		{
			for (int j = 0; j < n; j++)
			{
				for (int k = i; k < n; k++)
				{
					H1[i][j] += H[i][k] * Q[k][j];
				}
				if (i == j)H1[i][j] += mu;
			}
		}
		H = H1;
		P = Q;			//ע��Q���������任����
		return P;
	}
	else
	{
		int m = n - 1;
		double s = H[m - 1][m - 1] + H[n - 1][n - 1];
		double t = H[m - 1][m - 1] * H[n - 1][n - 1] - H[m - 1][n - 1] * H[n - 1][m - 1];
		x = H[0][0] * H[0][0] + H[0][1] * H[1][0] - s * H[0][0] + t;
		y = H[1][0] * (H[0][0] + H[1][1] - s);
		z = H[1][0] * H[2][1];
		for (int k = 0; k <= n - 3; k++)
		{
			vector<double> temp = { x,y,z };
			householder(temp, v, beta);

			//���ɲ����������任����
			//1.����
			Pnew = id_matrix(n);
			Psmall = household_matrix(v, beta);
			for (int i = k; i < k + 3; i++)
			{
				for (int j = k; j < k + 3; j++)
				{
					Pnew[i][j] = Psmall[i - k][j - k];
				}
			}
			//2.����(ע�⣺�µõ��ľ���Ӧ�����ҳ���ȥ)
			P = matrix_product(P, Pnew);

			int q = k >= 1 ? k : 1;
			matrix_householder_rowsimplify_p1q1p2q2(H, v, beta, k, q - 1, k + 3, n);
			int r = k + 4 <= n ? k + 4 : n;
			matrix_householder_colsimplify_p1q1p2q2(H, v, beta, 0, k, r, k + 3);
			x = H[k + 1][k];
			y = H[k + 2][k];
			if (k < n - 3)z = H[k + 3][k];
		}
		vector<double> temp = { x,y };
		vector<double> v(2);
		householder(temp, v, beta);

		//���ɲ����������任����
		//1.����
		Pnew = id_matrix(n);
		Psmall = household_matrix(v, beta);
		for (int i = n - 2; i < n; i++)
		{
			for (int j = n - 2; j < n; j++)
			{
				Pnew[i][j] = Psmall[i - (n - 2)][j - (n - 2)];
			}
		}
		//2.����(ע�⣺�µõ��ľ���Ӧ�����ҳ���ȥ)
		P = matrix_product(P, Pnew);

		matrix_householder_rowsimplify_p1q1p2q2(H, v, beta, n - 2, n - 3, n, n);
		matrix_householder_colsimplify_p1q1p2q2(H, v, beta, 0, n - 2, n, n);
		return P;	//��Ҫ���Ƿ���ֵ������
	}
}

unsigned implicit_qr(vector<vector<double>> &A)
{
	unsigned count = 0, maxcount = 100;
	int n = A.size();
	vector<vector<double>> Q, V, P;
	vector<vector<double>>A_pre(n, vector<double>(n, 0));
	Q = id_matrix(n);
	vector<double>b;

	//1.��Hessenberg��
	hessenberg_decomp(A, V, b);

	//����������hessenberg�ֽ�ı任����Q(Q^tAQ=H)��������V��b
	for (int i = 0; i < n - 2; i++)
	{
		vector<double> v = V[i];
		double beta = b[i];
		vector<vector<double>> Hnew = id_matrix(n);
		for (int p = i + 1; p < n; p++)
		{
			for (int q = i + 1; q < n; q++)
			{
				Hnew[p][q] -= beta * v[p - (i + 1)] * v[q - (i + 1)];
			}
		}
		Q = matrix_product(Q, Hnew);
	}

	//2.�������ж�
	//2.0�ζԽ�������Ԫ����Ϊ0
	for (int i = 2; i < n; i++)
	{
		for (int j = 0; j < i - 1; j++)
		{
			A[i][j] = 0;
		}
	}
	//2.1����С�ĴζԽ�����
	for (int i = 1; i < n; i++)
	{
		if (abs(A[i][i - 1]) <= (abs(A[i][i]) + abs(A[i - 1][i - 1]))*ERROR)
		{
			A[i][i - 1] = 0;
		}
	}
	//2.2���¶���ȷ�����ķǸ�����m����С�Ǹ�����k
	int m = 0;
	while (m < n)
	{
		if (m == n - 1 || A[n - m - 1][n - m - 2] == 0)
		{
			m += 1;
		}
		else
		{
			if (!judge_conjugate_root(A[n - m - 2][n - m - 2], A[n - m - 2][n - m - 1], A[n - m - 1][n - m - 2], A[n - m - 1][n - m - 1]))
			{
				break;
			}
			if (m <= n - 2 && A[n - m - 2][n - m - 3] != 0)
			{
				break;
			}
			m += 2;
		}
	}
	int k = n - m;
	if (k > 0)
	{
		k -= 1;
		while (k > 0)
		{
			if (A[k][k - 1] != 0)k--;
			else break;
		}
	}

	while (m != n && Frobenius_norm(matrix_minus(A, A_pre)) > ERROR && count < maxcount)	//���m=n�����,����ʱ����QR���������ٴν����ж�
	{
		A_pre = A;

		//3.QR����
		vector<vector<double>> H22 = sub_matrix(A, k, n - m, k, n - m);
		P = two_step_displacement_qr(H22);
		matrix_part_replace(A, k, n - m, k, n - m, H22);	//������ע�⣺Ҫ��ԭ�������H22�滻��

		//4.����任�����Լ����µ���hessenberg����
		vector<vector<double>> Pnew = id_matrix(n);
		matrix_part_replace(Pnew, k, n - m, k, n - m, P);
		Q = matrix_product(Q, Pnew);
		if (k > 0)
		{
			vector<vector<double>> H12 = sub_matrix(A, 0, k, k, n - m);
			H12 = matrix_product(H12, P);
			matrix_part_replace(A, 0, k, k, n - m, H12);
		}
		if (m > 0)
		{
			vector<vector<double>> H23 = sub_matrix(A, k, n - m, n - m, n);
			H23 = matrix_product(transpose(P), H23);
			matrix_part_replace(A, k, n - m, n - m, n, H23);
		}

		count++;

		//�ص��������ж�(2)
		//2.0�ζԽ�������Ԫ����Ϊ0
		for (int i = 2; i < n; i++)
		{
			for (int j = 0; j < i - 1; j++)
			{
				A[i][j] = 0;
			}
		}
		//2.1����С�ĴζԽ�����	
		for (int i = 1; i < n; i++)
		{
			if (abs(A[i][i - 1]) <= (abs(A[i][i]) + abs(A[i - 1][i - 1]))*ERROR)
			{
				A[i][i - 1] = 0;
			}
		}
		//2.2���¶���ȷ�����ķǸ�����m����С�Ǹ�����k
		while (m < n)
		{
			if (m == n - 1 || A[n - m - 1][n - m - 2] == 0)
			{
				m += 1;
			}
			else
			{
				if (m < n - 2 && A[n - m - 2][n - m - 3] != 0)
				{
					break;
				}
				if(judge_conjugate_root(A[n - m - 2][n - m - 2], A[n - m - 2][n - m - 1], A[n - m - 1][n - m - 2], A[n - m - 1][n - m - 1]))
				{
					m += 2;
					break;
				}
				//�������ر�ģ�����2*2������ʵ���ľ�������Ҫ�ر���2*2��QR������
				vector<vector<double>> H33 = sub_matrix(A, n - m - 2, n - m, n - m - 2, n - m);
				while (abs(H33[1][0]) > ERROR)
				{
					P = matrix_product(two_step_displacement_qr(H33), P);
				}
				vector<vector<double>> Pnew = id_matrix(n);
				matrix_part_replace(Pnew, n - m - 2, n - m, n - m - 2, n - m, P);
				Q = matrix_product(Q, Pnew);
				A = matrix_product(matrix_product(transpose(Pnew), A), Pnew);
				m += 2;
			}
		}
		int k = n - m;
		if (k > 0)
		{
			k -= 1;
			while (k > 0)
			{
				if (A[k][k - 1] != 0)k--;
				else break;
			}
		}
	}
	

	for (int i = 2; i < n; i++)
	{
		for (int j = 0; j < i - 1; j++)
		{
			A[i][j] = 0;
		}
	}
	for (int i = 1; i < n; i++)
	{
		if (abs(A[i][i - 1]) <= (abs(A[i][i]) + abs(A[i - 1][i - 1]))*ERROR)
		{
			A[i][i - 1] = 0;
		}
	}

	return count;
}

unsigned eigenvalues(vector<vector<double>> A, vector<complex<double>> &z)
{
	unsigned count = implicit_qr(A);
	int n = A.size();
	int k = 0;
	double delta, re, img;
	while (k < n)
	{
		if (k == n - 1 || A[k + 1][k] == 0)
		{
			z.push_back({ A[k][k],0 });
			k++;
		}
		else
		{
			delta = (A[k][k] - A[k + 1][k + 1])*(A[k][k] - A[k + 1][k + 1]) + 4 * A[k][k + 1] * A[k + 1][k];
			re = (A[k][k] + A[k + 1][k + 1]) / 2.0;
			img = sqrt(-delta) / 2.0;
			z.push_back({ re,img });
			z.push_back({ re,-img });
			k += 2;
		}
	}
	return count;
}