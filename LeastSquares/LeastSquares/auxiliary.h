#pragma once
#include<iostream>
#include<vector>
#include <random>
using namespace std;

void Ex1_Init1(vector<vector<double>>& A, vector<double>& b, int N);
void Ex1_Init2(vector<vector<double>>& A, vector<double>& b, int N);
void Ex1_Init3(vector<vector<double>>& A, vector<double>& b, int N);
void Ex2_Init(vector<vector<double>>& A, vector<double> t);
void random_value(vector<double>& b);//�������������
double ResultError(vector<double> y, vector<double> y0);//����Ax,b֮������,��������1����
double MSEError(vector<vector<double>> A, vector<double> b, vector<double> x);//�������ļ���

double vector_norm_1(vector<double> w);//������������1����
double vector_norm_infinity(vector<double> z);//�������������ģ
double vector_norm_2(vector<double> z);//������������2����
vector<double> vector_minus(vector<double> x, vector<double> y);//��������
vector<double> matrix_vector(vector<vector<double>> A, vector<double> x);//m*n�׾�����nά��������ˣ�����һ��mά������
