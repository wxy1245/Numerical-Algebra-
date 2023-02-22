#pragma once
#include<iostream>
#include<vector>
#include<complex>
using namespace std;

vector<double> QR_decomposition(vector<vector<double>> &A);//�������A��QR�ֽ�,���ؼ�¼�µ�������d��v���������Ѿ������A�����Խ������²���

double power_method(vector<vector<double>> A, vector<double> &u0, int &times);//�ݷ����㷽������ģ����ֵ����Ҫ�����ʼ����u0��

void householder(vector<double> x, vector<double> &v, double &beta);//����householder�任

void matrix_householder_simplify(vector<vector<double>> &A, vector<double> v, double beta);//��householder�任�������

void matrix_householder_rowsimplify_p1q1p2q2(vector<vector<double>> &A, vector<double> v, double beta, int p1, int q1, int p2, int q2);

void matrix_householder_colsimplify_p1q1p2q2(vector<vector<double>> &A, vector<double> v, double beta, int p1, int q1, int p2, int q2);

void matrix_householder_simplify_j(vector<vector<double>> &A, vector<double> v, double beta, int j);//��A��j��j���Ժ�Ĳ��ֽ���һ��householder����j=0�������������

void matrix_householder_rowsimplify_p_q(vector<vector<double>> &A, vector<double> v, double beta, int p, int q);//��A��p��q���Ժ�Ĳ��ֽ���һ���б任householder���� A:(m-p)*(n-q) v:m-p

void matrix_householder_colsimplify_p_q(vector<vector<double>> &A, vector<double> v, double beta, int p, int q);//��A��p��q���Ժ�Ĳ��ֽ���һ���б任householder���� A(m-p)*(n-q) v:n-q

void hessenberg_decomp(vector<vector<double>>&A, vector<vector<double>>&V, vector<double> &b);//��hessenberg�ֽ⣬Hessenberg������A�У�householder�任v��beta�����V��b��

vector<vector<double>> two_step_displacement_qr(vector<vector<double>> &H);//˫�ز�λ�Ƶ�QR������R�洢��H�У����������任����

unsigned implicit_qr(vector<vector<double>> &A);//��ʽQR�㷨,���ص�������

unsigned eigenvalues(vector<vector<double>> A, vector<complex<double>> &z);//����A������ֵ������QR��������