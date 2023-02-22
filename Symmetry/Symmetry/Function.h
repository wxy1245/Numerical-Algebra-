#pragma once
#include<iostream>
#include<vector>
#include<complex>
using namespace std;

void forward_subs(vector<vector<double>> L, vector<double>& b);//前代法

void forward_subs1(vector<vector<double>> L, vector<double>& b);//对角元为1的前代法

void back_subs(vector<vector<double>> U, vector<double>& b);//回代法

void back_subs1(vector<vector<double>> U, vector<double>& y);//对角元为1的回代法

void gauss_elim_col_pivoting(vector<vector<double>>& A, vector<int>& u);//列主元Gauss消去法

void vector_pb(vector<int> u, vector<double>&b);//计算向量P*b【可选】

void Jacobi_simplify_once(vector<vector<double>> &A, int p, int q, double c, double s);//一次Jacobi约化，其中变换矩阵只需要用三个参数来记录即可

void update_Q(vector<vector<double>> &A, int p, int q, double c, double s);//更新变换矩阵:Jacobi矩阵右乘Q

unsigned one_scan(vector<vector<double>> &A, vector<vector<double>> &Q, double delta);//一次扫描，返回变换的次数

unsigned Gate_Jacobi(vector<vector<double>> &A, vector<vector<double>> &Q, double sigma);//过关Jacobi方法,返回迭代次数. 其中Q是积累的变换矩阵（各列即为特征向量）,sigma是过关Jacobi选取的关值缩减因子.

unsigned sign_changes(vector<vector<double>> T, double x);//计算变号数

double Bisection(vector<vector<double>> A, unsigned m);//求实对称三对角阵的第m个特征值（从小到大）

void Inverse_power_method(vector<vector<double>> A, double mu, vector<double> &z, int &times);//再已经求得特征值的情况下，用反幂法求相应的特征向量,times:输入最大迭代次数，带回真实迭代次数