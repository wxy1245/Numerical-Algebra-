#pragma once
#include<iostream>
#include<vector>
using namespace std;

vector<double> matrix_vector(vector<vector<double>> A, vector<double> x);//m*n�׾�����nά��������ˣ�����һ��mά������

double inner_product(vector<double> x, vector<double> y);//�����ڻ�

vector<double> scale_vector(double alpha, vector<double> x);//��������

vector<double> vector_plus(vector<double> x, vector<double> y);//�����ӷ�

vector<double> vector_minus(vector<double> x, vector<double> y);//��������

double vector_norm_2(vector<double> z);//������������2����

double vector_norm_infinity(vector<double> z);//�������������ģ
