#pragma once
#include<iostream>
#include<vector>
#include<complex>
using namespace std;

void householder(vector<double> x, vector<double> &v, double &beta);//计算householder变换

void matrix_householder_rowsimplify_p_q(vector<vector<double>> &A, vector<double> v, double beta, int p, int q);//对A的p行q列以后的部分进行一步行变换householder化简 A:(m-p)*(n-q) v:m-p

void matrix_householder_colsimplify_p_q(vector<vector<double>> &A, vector<double> v, double beta, int p, int q);//对A的p行q列以后的部分进行一步列变换householder化简 A(m-p)*(n-q) v:n-q

void matrix_householder_rowsimplify_p1q1p2q2(vector<vector<double>> &A, vector<double> v, double beta, int p1, int q1, int p2, int q2);

void matrix_householder_colsimplify_p1q1p2q2(vector<vector<double>> &A, vector<double> v, double beta, int p1, int q1, int p2, int q2);

void two_diag(vector<vector<double>> &A, vector<vector<double>> &U, vector<vector<double>> &V);//二对角化,即U^t A V为二对角阵

void Givens_row(vector<vector<double>>& A, int p, int q, double c, double s);//Givnes变换(左乘)

void Givens_column(vector<vector<double>> &A, int p, int q, double c, double s);//Givens变换(右乘)

void Wilkinson_step(vector<vector<double>> &B, vector<vector<double>> &U, vector<vector<double>> &V);//1次Wilkinson位移SVD迭代

unsigned SVD(vector<vector<double>> &A, vector<vector<double>> &U, vector<vector<double>> &V);//SVD算法,U,V记录正交变换矩阵

vector<double> singular_values_sort(vector<vector<double>> D);//获取全部奇异值并从小到大排序
