#pragma once
#include<iostream>
#include<vector>
#include <random>
using namespace std;

void col_select(vector<vector<double>>& A, int k, int &p);//列主元的选取
void transpose(vector<vector<double>>& A);//方阵的转置
vector<double> sign_vector(vector<double> w);//返回w向量的符号向量
double vector_norm_1(vector<double> w);//计算列向量的1范数
double vector_norm_infinity(vector<double> z);//计算列向量最大模
double vector_dotproduct(vector<double> z, vector<double> x);//计算向量内积
vector<double> orthonormal_basis(int j, int n);//第j个分量为1，其余为0的一个标准基向量
void random_value(vector<double>& b);//随机向量的生成
vector<double> vector_minus(vector<double> x, vector<double> y);//向量减法
vector<double> matrix_vector(vector<vector<double>> A, vector<double> x);//矩阵与向量相乘