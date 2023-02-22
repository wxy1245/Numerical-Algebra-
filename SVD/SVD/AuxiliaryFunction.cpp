#include"auxiliary.h"
#include <cmath>

double vector_norm_1(vector<double> w)
{
	double sum = 0;
	int n = w.size();
	for (int i = 0; i < n; i++)
	{
		sum += abs(w[i]);
	}
	return sum;
}

double vector_norm_infinity(vector<double> z)
{
	int n = z.size();
	double maximum = 0;
	for (int i = 1; i < n; i++)
	{
		if (abs(z[i]) > maximum)
		{
			maximum = abs(z[i]);
		}
	}
	return maximum;	//返回的是取得最大模的下标
}

int vector_norm_maxindex(vector<double> z)
{
	int n = z.size();
	int k = 0;
	for (int i = 1; i < n; i++)
	{
		if (abs(z[i]) > abs(z[k]))
		{
			k = i;
		}
	}
	return k;	//返回的是取得最大模的下标
}

double vector_norm_2(vector<double> z)
{

	int n = z.size();
	double sum = 0;
	for (int i = 0; i < n; i++)
	{
		sum += z[i] * z[i];
	}
	return sqrt(sum);
}

vector<double> vector_minus(vector<double> x, vector<double> y)
{
	int n = x.size();
	vector<double> result(n);
	for (int i = 0; i < n; i++)
	{
		result[i] = x[i] - y[i];
	}
	return result;
}

vector<double> matrix_vector(vector<vector<double>> A, vector<double> x)
{
	int m = A.size();	//m取A的行数
	int n = A[0].size();//n取A的列数
	vector<double> result(m, 0);
	for (int i = 0; i < m; i++)
	{
		for (int j = 0; j < n; j++)
		{
			result[i] += A[i][j] * x[j];
		}
	}
	return result;
}

vector<vector<double>> transpose(vector<vector<double>> A)
{
	int m = A.size();
	int n = A[0].size();

	vector<vector<double>> At(n, vector<double>(m));
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < m; j++)
		{
			At[i][j] = A[j][i];
		}
	}
	return At;
}

double inner_product(vector<double> x, vector<double> y)
{
	int n = x.size();
	double result = 0;
	for (int i = 0; i < n; i++)
	{
		result += x[i] * y[i];
	}
	return result;
}

vector<double> scale_vector(double alpha, vector<double> x)
{
	int n = x.size();
	for (int i = 0; i < n; i++)
	{
		x[i] *= alpha;
	}
	return x;
}

vector<double> vector_plus(vector<double> x, vector<double> y)
{
	int n = x.size();
	vector<double> result(n);
	for (int i = 0; i < n; i++)
	{
		result[i] = x[i] + y[i];
	}
	return result;
}

void vector_print(vector<double> v)
{
	int n = v.size();
	for (int i = 0; i < n; i++)
	{
		cout.precision(3);
		cout << v[i] << " ";
	}
}

vector<vector<double>> matrix_product(vector<vector<double>> A, vector<vector<double>> B)
{
	int m = A.size();
	int n = A[0].size();
	int p = B[0].size();

	vector<vector<double>> C(m, vector<double>(p, 0));
	for (int i = 0; i < m; i++)
	{
		for (int j = 0; j < p; j++)
		{
			for (int k = 0; k < n; k++)
			{
				C[i][j] += A[i][k] * B[k][j];
			}
		}
	}

	return C;
}

vector<vector<double>> household_matrix(vector<double> v, double beta)
{
	int n = v.size();
	vector<vector<double>> H(n, vector<double>(n, 0));
	for (int i = 0; i < n; i++)H[i][i] = 1;
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			H[i][j] -= beta * v[i] * v[j];
		}
	}
	return H;
}

vector<vector<double>> id_matrix(int n)
{
	vector<vector<double>> I(n, vector<double>(n, 0));
	for (int i = 0; i < n; i++) I[i][i] = 1;
	return I;
}

bool judge_conjugate_root(double a, double b, double c, double d)
{
	double delta = (a - d)*(a - d) + 4 * b*c;//判别式
	if (delta < 0)return true;
	else return false;
}

vector<vector<double>> sub_matrix(vector<vector<double>> A, int i1, int i2, int j1, int j2)
{
	int m = i2 - i1;
	int n = j2 - j1;
	vector<vector<double>> subA(m, vector<double>(n));
	for (int i = 0; i < m; i++)
	{
		for (int j = 0; j < n; j++)
		{
			subA[i][j] = A[i + i1][j + j1];
		}
	}
	return subA;
}

void matrix_part_replace(vector<vector<double>>& A, int i1, int i2, int j1, int j2, vector<vector<double>> P)
{
	for (int i = i1; i < i2; i++)
	{
		for (int j = j1; j < j2; j++)
		{
			A[i][j] = P[i - i1][j - j1];
		}
	}
}

void matrix_print(vector<vector<double>> A)
{
	int m = A.size();
	if (m == 0)
	{
		cout << "A不存在" << "\n";
		return;
	}
	int n = A[0].size();
	for (int i = 0; i < m; i++)
	{
		for (int j = 0; j < n; j++)
		{
			cout.precision(3);
			cout << A[i][j] << "\t";
		}
		cout << "\n";
	}
}

double Frobenius_norm(vector<vector<double>> A)
{
	int n = A.size();
	double sumsquare = 0;
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			sumsquare += A[i][j] * A[i][j];
		}
	}
	return sqrt(sumsquare);
}

vector<vector<double>> matrix_minus(vector<vector<double>> A, vector<vector<double>> B)
{
	int m = A.size();
	int n = A[0].size();
	vector<vector<double>> C(m, vector<double>(n));
	for (int i = 0; i < m; i++)
	{
		for (int j = 0; j < n; j++)
		{
			C[i][j] = A[i][j] - B[i][j];
		}
	}
	return C;
}

double matrix_norm_infinity(vector<vector<double>> A)
{
	int n = A.size();
	double max_row_norm = 0, t;
	for (int i = 0; i < n; i++)
	{
		t = vector_norm_1(A[i]);
		if (t > max_row_norm)
		{
			max_row_norm = t;
		}
	}
	return max_row_norm;
}

double matrix_norm_max_elem(vector<vector<double>> A)
{
	int m = A.size();
	int n = A[0].size();
	double max_norm_elem = 0;
	for (int i = 0; i < m; i++)
	{
		for (int j = 0; j < n; j++)
		{
			if (abs(A[i][j]) > max_norm_elem)max_norm_elem = A[i][j];
		}
	}
	return max_norm_elem;
}
