#pragma once
#include<iostream>
#include<vector>
using namespace std;

unsigned Jacobi(vector<vector<double>> A, vector<double> b, vector<double> &x);//Jacobi����

unsigned GS(vector<vector<double>> A, vector<double> b, vector<double> &x);//Gauss-Seidel����

unsigned SOR(vector<vector<double>> A, vector<double> b, double omega, vector<double> &x);//SOR����

/*

unsigned ModelJacobi(vector<vector<double>> &U)//ģ�������Jacobi����

unsigned ModelGS(vector<vector<double>> &U)//ģ�������GS����

unsigned ModelSOR(vector<vector<double>> &U)//ģ�������SOR����

*/

