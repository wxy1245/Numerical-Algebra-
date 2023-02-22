#pragma once
#include<iostream>
#include<vector>

using namespace std;

void forward_subs(vector<vector<double>>& L, vector<double>& b);//前代法

void forward_subs1(vector<vector<double>>& L, vector<double>& b);//对角元为1的前代法

void back_subs(vector<vector<double>>& U, vector<double>& b);//回代法

void back_subs1(vector<vector<double>>& U, vector<double>& y);//对角元为1的回代法

void gauss_elim_col_pivoting(vector<vector<double>>& A, vector<int>& u);//列主元Gauss消去法

void vector_pb(vector<int>&u, vector<double>&b);//计算向量P*b【可选】

double infinity_norm(vector<vector<double>> A);//计算矩阵A的无穷范数(行和范数)

double inverse_infinity_norm(vector<vector<double>> A, vector<double> x);//计算矩阵A逆的无穷范数,需要用到计算解x