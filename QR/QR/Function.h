#pragma once
#include<iostream>
#include<vector>
#include<complex>
using namespace std;

vector<double> QR_decomposition(vector<vector<double>> &A);//计算矩阵A的QR分解,返回记录β的列向量d，v的内容则已经存放在A的主对角线以下部分

double power_method(vector<vector<double>> A, vector<double> &u0, int &times);//幂法计算方阵的最大模特征值（需要输入初始向量u0）

void householder(vector<double> x, vector<double> &v, double &beta);//计算householder变换

void matrix_householder_simplify(vector<vector<double>> &A, vector<double> v, double beta);//用householder变换化简矩阵

void matrix_householder_rowsimplify_p1q1p2q2(vector<vector<double>> &A, vector<double> v, double beta, int p1, int q1, int p2, int q2);

void matrix_householder_colsimplify_p1q1p2q2(vector<vector<double>> &A, vector<double> v, double beta, int p1, int q1, int p2, int q2);

void matrix_householder_simplify_j(vector<vector<double>> &A, vector<double> v, double beta, int j);//对A的j行j列以后的部分进行一步householder化简，j=0就是上面的特例

void matrix_householder_rowsimplify_p_q(vector<vector<double>> &A, vector<double> v, double beta, int p, int q);//对A的p行q列以后的部分进行一步行变换householder化简 A:(m-p)*(n-q) v:m-p

void matrix_householder_colsimplify_p_q(vector<vector<double>> &A, vector<double> v, double beta, int p, int q);//对A的p行q列以后的部分进行一步列变换householder化简 A(m-p)*(n-q) v:n-q

void hessenberg_decomp(vector<vector<double>>&A, vector<vector<double>>&V, vector<double> &b);//上hessenberg分解，Hessenberg阵存放在A中，householder变换v和beta存放在V与b中

vector<vector<double>> two_step_displacement_qr(vector<vector<double>> &H);//双重步位移的QR迭代，R存储在H中，返回正交变换矩阵

unsigned implicit_qr(vector<vector<double>> &A);//隐式QR算法,返回迭代次数

unsigned eigenvalues(vector<vector<double>> A, vector<complex<double>> &z);//求方阵A的特征值，返回QR迭代次数