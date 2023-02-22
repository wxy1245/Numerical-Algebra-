#include "Function.h"
#include"auxiliary.h"
#include <cmath>

void forward_subs(vector<vector<double>>& L, vector<double>& b)
{//ǰ�����������Ǿ���
	int j;
	int n = L.size();
	for (j = 0; j < n - 1; j++)
	{
		b[j] = b[j] / L[j][j];
		for (int k = j + 1; k < n; k++)
		{
			b[k] -= b[j] * L[k][j];
		}
	}
	b[n - 1] = b[n - 1] / L[n - 1][n - 1];
}

void forward_subs1(vector<vector<double>>& L, vector<double>& b)
{//�Խ�ԪΪ1ʱ��ǰ����
	int j;
	int n = L.size();
	for (j = 0; j < n - 1; j++)
	{
		for (int k = j + 1; k < n; k++)
		{
			b[k] -= b[j] * L[k][j];
		}
	}
}

void back_subs(vector<vector<double>>& U, vector<double>& y)
{//�ش���
	int j;
	int n = U.size();
	for (j = n - 1; j > 0; j--)
	{
		y[j] = y[j] / U[j][j];
		for (int k = 0; k < j; k++)
		{
			y[k] -= y[j] * U[k][j];
		}
	}
	y[0] = y[0] / U[0][0];
}

void back_subs1(vector<vector<double>>& U, vector<double>& y)
{//�Խ�ԪΪ1�Ļش���
	int j;
	int n = U.size();
	for (j = n - 1; j > 0; j--)
	{
		for (int k = 0; k < j; k++)
		{
			y[k] -= y[j] * U[k][j];
		}
	}
}

void householder(vector<double> x, vector<double> &v, double &beta)
{
	int n = x.size();
	double eta = vector_norm_infinity(x);

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

void Qt_b(vector<vector<double>> A, vector<double>& b, vector<double> d)
{//����Q^T b = Hn...H1 b,n:�任����
	int m = A.size();
	int n = d.size();	//!!!ע��nȡ�ܱ任�Ĵ���(d��ά��,���¸���),����A������(���A������A����,��任������<����,��ȡ�任��������һ����)
	vector<double> v;
	double tempsum;
	vector<double> temp;
	for (int j = 0; j < n; j++)	//���μ���b1=H1b,...,bn=Hnbn-1
	{
		//�Ȱ���һ�ε�householder�任�õ�vjȡ����(��СΪm-j,��ԪΪ1)
		v.resize(m - j);
		v[0] = 1;
		for (int i = 1; i < m - j; i++)
		{
			v[i] = A[j + i][j];
		}

		//��b�ĺ�m-j�����householder�任
		temp.resize(m - j);	//������ʱ���b�仯���ֵ
		for (int i = 0; i < m - j; i++)	//��b��ʱ��Ž�ȥ
		{
			temp[i] = b[i + j];
		}
		for (int i = j; i < m; i++)	//�����temp,���任���b
		{
			tempsum = 0;
			for (int k = j; k < m; k++)
			{
				tempsum += v[i - j] * v[k - j] * b[k];
			}
			temp[i - j] -= d[j] * tempsum;	//ע�⣺���������b��ʱ������b��ÿ����͵�ʱ��Ҫ�ã������в��ɸı䡣���Ա���������һ����ʱ����temp���䱣��
		}
		for (int i = j; i < m; i++)	//���b��ֵtemp���һ�β���
		{
			b[i] = temp[i - j];
		}
	}
}

void QRSolve(vector<vector<double>> &A, vector<double> &b)
{//����AΪ���淽��ʱ��QR�ֽ����Ax=b�ķ���
	vector<double> d;
	d = QR_decomposition(A);	//������ע�⣺����dֻ��n-1ά; QR�ֽ����,Q=H1...Hn-1,�Ҽ�¼H�任��v������A��������,�±�����d��R������A�������ǣ�
	//���Q R x = b,��R x = Q^T b = Hn-1...H1 b��
	Qt_b(A, b, d);	//�ȼ���Q^T b,������b����
	back_subs(A, b);	//�����Rx = Q^t b,���սⱣ����b����
}

vector<double> QRLeastSquares(vector<vector<double>> A, vector<double> b)
{//��������(m>n),����С��������Ľ�;����,����QRSolve���Կ����˵�����
	int n = A[0].size();	//ȡA������
	vector<double> d;

	//��ⲽ��: 1.QR�ֽ� 2.c1 = Q1^T b(ע��c1ͨ��ȡQ^T bǰn�оͿ��Եõ���) 3.R x = c1
	d = QR_decomposition(A);

	Qt_b(A, b, d);
	vector<double>c1(b.begin(), b.begin() + n);	//��Ƭ����,ȡǰn��

	vector<vector<double>> A1(n, vector<double>(n));	//Ϊ���ܹ�Ӧ�ûش�������ҪȡA�ķ����ǿ�(ǰn��n��,�����ǲ��ּ�ΪR)
	for (int i = 0; i < n; i++)	
	{
		for (int j = 0; j < n; j++)
		{
			A1[i][j] = A[i][j];
		}
	}
	back_subs(A1, c1);
	return c1;	//���ؽ�(nά����)
}