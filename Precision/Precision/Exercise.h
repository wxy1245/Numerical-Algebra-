#pragma once
#include <iostream>
#include <vector>
//#include<Windows.h>
#include"Function.h"
//#include<stdlib.h>

void ColGaussSolve(vector<vector<double>>& A, vector<double>& b);//����ԪGauss��ȥ��������Է�����

double Col_Precision(vector<vector<double>>& A, vector<double>& b);//����Ԫ��ȥ���������ȹ���(ע�⣺����ԭ���������)

void exercise_1();
void exercise_2();