#pragma once
#include<iostream>
#include<vector>
using namespace std;

//һЩ��������

void Ex1_Init(vector<vector<double>>& A, vector<double>& b, int N);
void Ex2_Init1(vector<vector<double>>& A, vector<double>& b, int N);
void Ex2_Init2(vector<vector<double>>& A, vector<double>& b, int N);
void random_value(vector<double>& b);//�������������
void full_select(vector<vector<double>>& A, int k, int &p, int &q);//ȫ��Ԫ��ѡȡ
void col_select(vector<vector<double>>& A, int k, int &p);//����Ԫ��ѡȡ
double ResultError(vector<double> y, vector<double> y0);//����Ax,b֮������,��������1����
