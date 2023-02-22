#pragma once
#include<iostream>
#include<vector>
#include<complex>
using namespace std;

void householder(vector<double> x, vector<double> &v, double &beta);//����householder�任

void matrix_householder_rowsimplify_p_q(vector<vector<double>> &A, vector<double> v, double beta, int p, int q);//��A��p��q���Ժ�Ĳ��ֽ���һ���б任householder���� A:(m-p)*(n-q) v:m-p

void matrix_householder_colsimplify_p_q(vector<vector<double>> &A, vector<double> v, double beta, int p, int q);//��A��p��q���Ժ�Ĳ��ֽ���һ���б任householder���� A(m-p)*(n-q) v:n-q

void matrix_householder_rowsimplify_p1q1p2q2(vector<vector<double>> &A, vector<double> v, double beta, int p1, int q1, int p2, int q2);

void matrix_householder_colsimplify_p1q1p2q2(vector<vector<double>> &A, vector<double> v, double beta, int p1, int q1, int p2, int q2);

void two_diag(vector<vector<double>> &A, vector<vector<double>> &U, vector<vector<double>> &V);//���Խǻ�,��U^t A VΪ���Խ���

void Givens_row(vector<vector<double>>& A, int p, int q, double c, double s);//Givnes�任(���)

void Givens_column(vector<vector<double>> &A, int p, int q, double c, double s);//Givens�任(�ҳ�)

void Wilkinson_step(vector<vector<double>> &B, vector<vector<double>> &U, vector<vector<double>> &V);//1��Wilkinsonλ��SVD����

unsigned SVD(vector<vector<double>> &A, vector<vector<double>> &U, vector<vector<double>> &V);//SVD�㷨,U,V��¼�����任����

vector<double> singular_values_sort(vector<vector<double>> D);//��ȡȫ������ֵ����С��������
