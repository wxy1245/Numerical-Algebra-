#pragma once
#include<iostream>
#include<vector>

using namespace std;

void forward_subs(vector<vector<double>>& L, vector<double>& b);//ǰ����

void forward_subs1(vector<vector<double>>& L, vector<double>& b);//�Խ�ԪΪ1��ǰ����

void back_subs(vector<vector<double>>& U, vector<double>& b);//�ش���

void back_subs1(vector<vector<double>>& U, vector<double>& y);//�Խ�ԪΪ1�Ļش���

void gauss_elim_col_pivoting(vector<vector<double>>& A, vector<int>& u);//����ԪGauss��ȥ��

void vector_pb(vector<int>&u, vector<double>&b);//��������P*b����ѡ��

double infinity_norm(vector<vector<double>> A);//�������A�������(�кͷ���)

double inverse_infinity_norm(vector<vector<double>> A, vector<double> x);//�������A��������,��Ҫ�õ������x