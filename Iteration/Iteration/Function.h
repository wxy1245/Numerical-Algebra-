#pragma once
#include<iostream>
#include<vector>
using namespace std;

unsigned Jacobi(vector<vector<double>> A, vector<double> b, vector<double> &x);//Jacobi迭代

unsigned GS(vector<vector<double>> A, vector<double> b, vector<double> &x);//Gauss-Seidel迭代

unsigned SOR(vector<vector<double>> A, vector<double> b, double omega, vector<double> &x);//SOR迭代

/*

unsigned ModelJacobi(vector<vector<double>> &U)//模型问题的Jacobi迭代

unsigned ModelGS(vector<vector<double>> &U)//模型问题的GS迭代

unsigned ModelSOR(vector<vector<double>> &U)//模型问题的SOR迭代

*/

