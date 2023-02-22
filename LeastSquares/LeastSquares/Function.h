#pragma once
#include<iostream>
#include<vector>

using namespace std;

void forward_subs(vector<vector<double>>& L, vector<double>& b);//前代法

void forward_subs1(vector<vector<double>>& L, vector<double>& b);//对角元为1的前代法

void back_subs(vector<vector<double>>& U, vector<double>& b);//回代法

void back_subs1(vector<vector<double>>& U, vector<double>& y);//对角元为1的回代法

void householder(vector<double> x, vector<double> &v, double &beta);//计算householder变换

void matrix_householder_simplify(vector<vector<double>> &A, vector<double> v, double beta);//用householder变换化简矩阵

void matrix_householder_simplify_j(vector<vector<double>> &A, vector<double> v, double beta, int j);//对A的j行j列以后的部分进行一步householder化简，j=0就是上面的特例

vector<double> QR_decomposition(vector<vector<double>> &A);//计算矩阵A的QR分解,返回记录β的列向量d，v的内容则已经存放在A的主对角线以下部分

void Qt_b(vector<vector<double>> A, vector<double> &b, vector<double> d);//计算Q^T b = Hn...H1 b

void QRSolve(vector<vector<double>> &A, vector<double> &b);//(A:可逆方阵)用QR分解求解线性方程组

vector<double> QRLeastSquares(vector<vector<double>> A, vector<double> b);//用QR分解求解最小二乘问题
