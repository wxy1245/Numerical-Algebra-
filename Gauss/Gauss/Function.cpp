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

void gauss_elim(vector<vector<double>>& A)
{//Gauss��ȥ��
	int n = A.size();
	for (int k = 0; k < n - 1; k++)
	{//�������A(k)
		for (int i = k + 1; i < n; i++)
		{
			A[i][k] = A[i][k] / A[k][k];//k�����Խ�����0���lk����Ԫ
			for (int j = k + 1; j < n; j++)
			{
				A[i][j] -= A[i][k] * A[k][j];
			}
		}
	}
}

void gauss_elim_full_pivoting(vector<vector<double>>& A, vector<int>& u, vector<int>& v)
{//ȫ��ԪGauss��ȥ��
	double temp;	//ע�⣺����Ҫ��double���ͣ���Ҫ���־���int��
	int n = A.size();
	int p, q;
	for (int k = 0; k < n - 1; k++)
	{	
		full_select(A, k, p, q);
		for (int j = 0; j < n; j++)
		{
			temp = A[p][j]; A[p][j] = A[k][j]; A[k][j] = temp;	//����k,p��
		}
		for (int i = 0; i < n; i++)
		{
			temp = A[i][q]; A[i][q] = A[i][k]; A[i][k] = temp;	//����k,q��
		}
		//u[k] = p; v[k] = q;	//��¼��������
		u.push_back(p); v.push_back(q);
		if (A[k][k] != 0)
		{
			for (int i = k + 1; i < n; i++)
			{
				A[i][k] = A[i][k] / A[k][k];	//A(k)���Խ����²���=0���ɴ��Ltiuta��Ԫ
				for (int j = k + 1; j < n; j++)
				{
					A[i][j] -= A[i][k] * A[k][j];
				}
			}
		}
		else
		{
			cout << "��ȥ����" << "\n";
			break;
		}
	}	
}

void gauss_elim_col_pivoting(vector<vector<double>>& A, vector<int>& u)
{//����ԪGauss��ȥ��
	double temp;
	int n = A.size();
	int p;
	for (int k = 0; k < n - 1; k++)
	{
		col_select(A, k, p);
		for (int j = 0; j < n; j++)
		{
			temp = A[p][j]; A[p][j] = A[k][j]; A[k][j] = temp;	//����k,p��
		}
		u.push_back(p);	//��¼�û�����
		if (A[k][k] != 0)
		{
			for (int i = k + 1; i < n; i++)
			{
				A[i][k] = A[i][k] / A[k][k];	//A(k)���Խ����²���=0���ɴ��Ltiuta��Ԫ
				for (int j = k + 1; j < n; j++)
				{
					A[i][j] -= A[i][k] * A[k][j];
				}
			}
		}
		else
		{
			cout << "�������죬�ֽ��޷���������" << "\n";
			exit(0);
		}
	}
}

void vector_pb(vector<int>& u, vector<double>& b)
{//��������P*b����ѡ��;P=Pr...P1
	double temp;
	int r = u.size();
	for (int i = 0; i < r; i++)
	{
		temp = b[i]; b[i] = b[u[i]]; b[u[i]] = temp;	//����b��i��u(i)����
	}
}

void vector_qb(vector<int>& v, vector<double>& b)
{//��������Q*b����ѡ��;Q=Q1...Qr
	double temp;
	int r = v.size();
	for (int i = r - 1; i >= 0; i--)
	{
		temp = b[i]; b[i] = b[v[i]]; b[v[i]] = temp;	//����b��i��v(i)����
	}
}

void cholesky_decomp(vector<vector<double>>& A)
{//�Գ��������׼Cholesky�ֽ�
	int n = A.size();
	for (int k = 0; k < n; k++)
	{
		A[k][k] = sqrt(A[k][k]);
		for (int i = k + 1; i < n; i++)
		{
			A[i][k] = A[i][k] / A[k][k];
		}
		for (int j = k + 1; j < n; j++)
		{
			for (int i = j; i < n; i++)
			{
				A[i][j] = A[i][j] - A[i][k] * A[j][k];
			}
		}
	}
}

void modified_cholesky_decomp(vector<vector<double>>& A)
{//�Ľ���ƽ������
	int n = A.size();
	vector<double> v(n);	//����double���͵���ʱ����v
	for (int j = 0; j < n; j++)
	{
		for (int k = 0; k <= j - 1; k++)
		{
			v[k] = A[j][k] * A[k][k];
		}
		for (int k = 0; k <= j - 1; k++)
		{
			A[j][j] -= A[j][k]*v[k];	//���Խ����ϴ洢�Խ�Ԫ
		}
		for (int i = j + 1; i < n; i++)
		{
			for (int k = 0; k <= j - 1; k++)
			{
				A[i][j] -= A[i][k] * v[k];	//�����߼�Ҫע��:�Ƚ�����ͣ������ԶԽ�Ԫ
			}
			A[i][j] /= A[j][j];
		}
	}
}

void matrix_DLT(vector<vector<double>>& A)
{//�������D*L^T����ѡ��,����A�������ǲ��֣�����DԴ��A�ĶԽǣ�L��Դ��A�����Խ�������
	int N = A.size();
	for (int i = 0; i < N - 1; i++)
	{
		for (int j = i + 1; j < N; j++)
		{
			A[i][j] = A[i][i] * A[j][i];
		}
	}
}

void GaussSolve(vector<vector<double>>& A, vector<double>& b)
{//Gauss��ȥ��������Է�����
	int N = A.size();
	gauss_elim(A);	//A���LU�ֽ�֮�󣬾Ϳ����Ⱥ��������������Ly=b,Ux=y���õ����ս�x
	forward_subs1(A, b);	//ע�⵽LΪ��λ�������󣬿��öԽ�ԪΪ1��ǰ����
	back_subs(A, b);		//ע��ǰ�������֮��b��y��ȡ��;���ش�������b�洢�ľ���x(���ս�)
}

void FullGaussSolve(vector<vector<double>>& A, vector<double>& b)
{//ȫ��ԪGauss��ȥ��������Է�����
	int N = A.size();
	vector<int> u, v;
	gauss_elim_full_pivoting(A, u, v);	//PAQ=LU��ɷֽ�֮���Ⱥ����Ly=Pb,Uz=y,�����x=Qz�õ����
	vector_pb(u, b);
	forward_subs1(A, b);
	back_subs(A, b);
	vector_qb(v, b);
}

void ColGaussSolve(vector<vector<double>>& A, vector<double>& b)
{//����ԪGauss��ȥ��������Է�����
	int N = A.size();
	vector<int> u;
	gauss_elim_col_pivoting(A, u);	//PA=LU��ɷֽ�֮���Ⱥ����Ly=Pb,Ux=y
	vector_pb(u, b);
	forward_subs1(A, b);
	back_subs(A, b);
}

void CholeskySolve(vector<vector<double>>& A, vector<double>& b)
{
	int N = A.size();
	cholesky_decomp(A);	//����A��cholesky�ֽ�
	for (int i = 0; i < N - 1; i++)	//��A�������ǲ��ָ���ΪL��ת��
	{
		for (int j = i + 1; j < N; j++)
		{
			A[i][j] = A[j][i];
		}
	}
	forward_subs(A, b);	//���Ly=b
	back_subs(A, b);	//���L^t x=y
}

void ModifiedCholeskySolve(vector<vector<double>>& A, vector<double>& b)
{
	int N = A.size();
	modified_cholesky_decomp(A);	//�Ƚ��иĽ���cholesky�ֽ�
	matrix_DLT(A);	//A�������ǲ��ִ��D*L^T
	forward_subs1(A, b);	//���Ly=b
	back_subs(A, b);	//���Ux=y
}
