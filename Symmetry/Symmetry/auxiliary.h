#pragma once
#include<iostream>
#include<vector>
#include <random>
using namespace std;

void col_select(vector<vector<double>>& A, int k, int &p);//列主元的选取
double vector_norm_1(vector<double> w);//计算列向量的1范数
double vector_norm_infinity(vector<double> z);//计算列向量最大模
int vector_norm_maxindex(vector<double> z);//计算列向量最大模下标 !返回值类型为int
double vector_norm_2(vector<double> z);//计算列向量的2范数
vector<double> vector_minus(vector<double> x, vector<double> y);//向量减法
vector<double> matrix_vector(vector<vector<double>> A, vector<double> x);//m*n阶矩阵与n维列向量相乘，返回一个m维列向量
vector<vector<double>> transpose(vector<vector<double>> A);//矩阵转置
double inner_product(vector<double> x, vector<double> y);//向量内积
vector<double> scale_vector(double alpha, vector<double> x);//数乘向量
vector<double> vector_plus(vector<double> x, vector<double> y);//向量加法
void vector_print(vector<double> v);//打印向量
vector<vector<double>> matrix_product(vector<vector<double>> A, vector<vector<double>> B);//矩阵乘法
vector<vector<double>> household_matrix(vector<double> v, double beta);//计算household变换矩阵
vector<vector<double>> id_matrix(int n);//生成单位矩阵
bool judge_conjugate_root(double a, double b, double c, double d);//计算矩阵A=((a,b),(c,d))的特征值是否为共轭复数
vector<vector<double>> sub_matrix(vector<vector<double>> A, int i1, int i2, int j1, int j2);
void matrix_part_replace(vector<vector<double>> &A, int i1, int i2, int j1, int j2, vector<vector<double>> P);//矩阵的copy
void matrix_print(vector<vector<double>> A);//打印矩阵
double Frobenius_norm(vector<vector<double>> A);// 方阵的Frobenius范数
vector<vector<double>> matrix_minus(vector<vector<double>> A, vector<vector<double>> B);//矩阵减法

double NotDiagonalNorm(vector<vector<double>> A);//非对角范数
double matrix_norm_infinity(vector<vector<double>> A);//矩阵的无穷范数（行和范数）
vector<double> get_sort_eigenvalues(vector<vector<double>> D);