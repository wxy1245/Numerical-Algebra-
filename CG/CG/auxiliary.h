#pragma once
#include<iostream>
#include<vector>
using namespace std;

vector<double> matrix_vector(vector<vector<double>> A, vector<double> x);//m*n阶矩阵与n维列向量相乘，返回一个m维列向量

double inner_product(vector<double> x, vector<double> y);//向量内积

vector<double> scale_vector(double alpha, vector<double> x);//数乘向量

vector<double> vector_plus(vector<double> x, vector<double> y);//向量加法

vector<double> vector_minus(vector<double> x, vector<double> y);//向量减法

double vector_norm_2(vector<double> z);//计算列向量的2范数

double vector_norm_infinity(vector<double> z);//计算列向量最大模
