#pragma once
#include<iostream>
#include<vector>

using namespace std;

void forward_subs(vector<vector<double>>& L, vector<double>& b);//ǰ����

void forward_subs1(vector<vector<double>>& L, vector<double>& b);//�Խ�ԪΪ1��ǰ����

void back_subs(vector<vector<double>>& U, vector<double>& b);//�ش���

void back_subs1(vector<vector<double>>& U, vector<double>& y);//�Խ�ԪΪ1�Ļش���

void householder(vector<double> x, vector<double> &v, double &beta);//����householder�任

void matrix_householder_simplify(vector<vector<double>> &A, vector<double> v, double beta);//��householder�任�������

void matrix_householder_simplify_j(vector<vector<double>> &A, vector<double> v, double beta, int j);//��A��j��j���Ժ�Ĳ��ֽ���һ��householder����j=0�������������

vector<double> QR_decomposition(vector<vector<double>> &A);//�������A��QR�ֽ�,���ؼ�¼�µ�������d��v���������Ѿ������A�����Խ������²���

void Qt_b(vector<vector<double>> A, vector<double> &b, vector<double> d);//����Q^T b = Hn...H1 b

void QRSolve(vector<vector<double>> &A, vector<double> &b);//(A:���淽��)��QR�ֽ�������Է�����

vector<double> QRLeastSquares(vector<vector<double>> A, vector<double> b);//��QR�ֽ������С��������
