#pragma once
#include<iostream>
#include<vector>
using namespace std;

unsigned Jacobi(vector<vector<double>> A, vector<double> b, vector<double> &x, double epsilon);//Jacobi����

unsigned GS(vector<vector<double>> A, vector<double> b, vector<double> &x, double epsilon);//Gauss-Seidel����

unsigned SOR(vector<vector<double>> A, vector<double> b, double omega, vector<double> &x, double epsilon);//SOR����

unsigned CG(vector<vector<double>> A, vector<double> b, vector<double> &x, double epsilon);//CG��

unsigned CGTest(vector<vector<double>> A, vector<double> b, vector<double> &x, double epsilon);//CG���ԣ������м�����������
