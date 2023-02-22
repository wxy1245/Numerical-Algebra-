#pragma once
#include<iostream>
#include<vector>
#include <random>
using namespace std;

void col_select(vector<vector<double>>& A, int k, int &p);//����Ԫ��ѡȡ
double vector_norm_1(vector<double> w);//������������1����
double vector_norm_infinity(vector<double> z);//�������������ģ
int vector_norm_maxindex(vector<double> z);//�������������ģ�±� !����ֵ����Ϊint
double vector_norm_2(vector<double> z);//������������2����
vector<double> vector_minus(vector<double> x, vector<double> y);//��������
vector<double> matrix_vector(vector<vector<double>> A, vector<double> x);//m*n�׾�����nά��������ˣ�����һ��mά������
vector<vector<double>> transpose(vector<vector<double>> A);//����ת��
double inner_product(vector<double> x, vector<double> y);//�����ڻ�
vector<double> scale_vector(double alpha, vector<double> x);//��������
vector<double> vector_plus(vector<double> x, vector<double> y);//�����ӷ�
void vector_print(vector<double> v);//��ӡ����
vector<vector<double>> matrix_product(vector<vector<double>> A, vector<vector<double>> B);//����˷�
vector<vector<double>> household_matrix(vector<double> v, double beta);//����household�任����
vector<vector<double>> id_matrix(int n);//���ɵ�λ����
bool judge_conjugate_root(double a, double b, double c, double d);//�������A=((a,b),(c,d))������ֵ�Ƿ�Ϊ�����
vector<vector<double>> sub_matrix(vector<vector<double>> A, int i1, int i2, int j1, int j2);
void matrix_part_replace(vector<vector<double>> &A, int i1, int i2, int j1, int j2, vector<vector<double>> P);//�����copy
void matrix_print(vector<vector<double>> A);//��ӡ����
double Frobenius_norm(vector<vector<double>> A);// �����Frobenius����
vector<vector<double>> matrix_minus(vector<vector<double>> A, vector<vector<double>> B);//�������

double NotDiagonalNorm(vector<vector<double>> A);//�ǶԽǷ���
double matrix_norm_infinity(vector<vector<double>> A);//�������������кͷ�����
vector<double> get_sort_eigenvalues(vector<vector<double>> D);