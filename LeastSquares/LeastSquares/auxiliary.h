#pragma once
#include<iostream>
#include<vector>
#include <random>
using namespace std;

void Ex1_Init1(vector<vector<double>>& A, vector<double>& b, int N);
void Ex1_Init2(vector<vector<double>>& A, vector<double>& b, int N);
void Ex1_Init3(vector<vector<double>>& A, vector<double>& b, int N);
void Ex2_Init(vector<vector<double>>& A, vector<double> t);
void random_value(vector<double>& b);//随机向量的生成
double ResultError(vector<double> y, vector<double> y0);//给出Ax,b之间的误差,用向量的1范数
double MSEError(vector<vector<double>> A, vector<double> b, vector<double> x);//均方误差的计算

double vector_norm_1(vector<double> w);//计算列向量的1范数
double vector_norm_infinity(vector<double> z);//计算列向量最大模
double vector_norm_2(vector<double> z);//计算列向量的2范数
vector<double> vector_minus(vector<double> x, vector<double> y);//向量减法
vector<double> matrix_vector(vector<vector<double>> A, vector<double> x);//m*n阶矩阵与n维列向量相乘，返回一个m维列向量
