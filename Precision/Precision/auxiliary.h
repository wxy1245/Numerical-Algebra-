#pragma once
#include<iostream>
#include<vector>
#include <random>
using namespace std;

void col_select(vector<vector<double>>& A, int k, int &p);//����Ԫ��ѡȡ
void transpose(vector<vector<double>>& A);//�����ת��
vector<double> sign_vector(vector<double> w);//����w�����ķ�������
double vector_norm_1(vector<double> w);//������������1����
double vector_norm_infinity(vector<double> z);//�������������ģ
double vector_dotproduct(vector<double> z, vector<double> x);//���������ڻ�
vector<double> orthonormal_basis(int j, int n);//��j������Ϊ1������Ϊ0��һ����׼������
void random_value(vector<double>& b);//�������������
vector<double> vector_minus(vector<double> x, vector<double> y);//��������
vector<double> matrix_vector(vector<vector<double>> A, vector<double> x);//�������������