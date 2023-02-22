#pragma once
#include<iostream>
#include<vector>
using namespace std;

//一些辅助函数

void Ex1_Init(vector<vector<double>>& A, vector<double>& b, int N);
void Ex2_Init1(vector<vector<double>>& A, vector<double>& b, int N);
void Ex2_Init2(vector<vector<double>>& A, vector<double>& b, int N);
void random_value(vector<double>& b);//随机向量的生成
void full_select(vector<vector<double>>& A, int k, int &p, int &q);//全主元的选取
void col_select(vector<vector<double>>& A, int k, int &p);//列主元的选取
double ResultError(vector<double> y, vector<double> y0);//给出Ax,b之间的误差,用向量的1范数
