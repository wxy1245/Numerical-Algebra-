#pragma once
#include<iostream>
#include<vector>
#include<complex>
using namespace std;

void forward_subs(vector<vector<double>> L, vector<double>& b);//ǰ����

void forward_subs1(vector<vector<double>> L, vector<double>& b);//�Խ�ԪΪ1��ǰ����

void back_subs(vector<vector<double>> U, vector<double>& b);//�ش���

void back_subs1(vector<vector<double>> U, vector<double>& y);//�Խ�ԪΪ1�Ļش���

void gauss_elim_col_pivoting(vector<vector<double>>& A, vector<int>& u);//����ԪGauss��ȥ��

void vector_pb(vector<int> u, vector<double>&b);//��������P*b����ѡ��

void Jacobi_simplify_once(vector<vector<double>> &A, int p, int q, double c, double s);//һ��JacobiԼ�������б任����ֻ��Ҫ��������������¼����

void update_Q(vector<vector<double>> &A, int p, int q, double c, double s);//���±任����:Jacobi�����ҳ�Q

unsigned one_scan(vector<vector<double>> &A, vector<vector<double>> &Q, double delta);//һ��ɨ�裬���ر任�Ĵ���

unsigned Gate_Jacobi(vector<vector<double>> &A, vector<vector<double>> &Q, double sigma);//����Jacobi����,���ص�������. ����Q�ǻ��۵ı任���󣨸��м�Ϊ����������,sigma�ǹ���Jacobiѡȡ�Ĺ�ֵ��������.

unsigned sign_changes(vector<vector<double>> T, double x);//��������

double Bisection(vector<vector<double>> A, unsigned m);//��ʵ�Գ����Խ���ĵ�m������ֵ����С����

void Inverse_power_method(vector<vector<double>> A, double mu, vector<double> &z, int &times);//���Ѿ��������ֵ������£��÷��ݷ�����Ӧ����������,times:����������������������ʵ��������