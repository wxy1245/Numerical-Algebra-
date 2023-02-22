#pragma once
#include<iostream>
#include<vector>
using namespace std;

unsigned Jacobi(vector<vector<double>> A, vector<double> b, vector<double> &x, double epsilon);//Jacobi迭代

unsigned GS(vector<vector<double>> A, vector<double> b, vector<double> &x, double epsilon);//Gauss-Seidel迭代

unsigned SOR(vector<vector<double>> A, vector<double> b, double omega, vector<double> &x, double epsilon);//SOR迭代

unsigned CG(vector<vector<double>> A, vector<double> b, vector<double> &x, double epsilon);//CG法

unsigned CGTest(vector<vector<double>> A, vector<double> b, vector<double> &x, double epsilon);//CG测试（带有中间变量的输出）
