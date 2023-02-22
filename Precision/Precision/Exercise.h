#pragma once
#include <iostream>
#include <vector>
//#include<Windows.h>
#include"Function.h"
//#include<stdlib.h>

void ColGaussSolve(vector<vector<double>>& A, vector<double>& b);//列主元Gauss消去法求解线性方程组

double Col_Precision(vector<vector<double>>& A, vector<double>& b);//列主元消去法给出精度估计(注意：传入原矩阵和向量)

void exercise_1();
void exercise_2();