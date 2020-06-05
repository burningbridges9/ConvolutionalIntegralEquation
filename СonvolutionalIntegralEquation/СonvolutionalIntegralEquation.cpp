// СonvolutionalIntegralEquation.cpp : Этот файл содержит функцию "main". Здесь начинается и заканчивается выполнение программы.
//

#include <iostream>
#include "IntegralCalculator.h"
#include "Well.h"
#include <cmath>
#include "СonvolutionalIntegralEquation.h"
#include "SequentialSolver.h"
#define _USE_MATH_DEFINES
#define M_PI 3.14159265358979323846
using namespace std;


int * GetIndexes(Well * wells)
{
	auto *indexes = new int[3];
	for (int i = 0; i < 3; i++)
	{
		int tempVal = 0;
		for (int j = 0; j <= i; j++)
		{
			tempVal += wells[j].N;
		}
		indexes[i] = tempVal - i;
	}
	return indexes;
}

int main()
{
	auto * wells = new Well[3];
	for (int i = 1; i <= 3; i++)
	{
		Well * well = new Well(5 * i, 5 * i, 3, 5 * (i - 1), 5 * i, 1, 1, 0.1, 10, 4, 0.3, 0, 100);
		wells[i - 1] = *well;
	}
	auto * indexes = GetIndexes(wells);
	auto sqSolver = SequentialSolver(wells, indexes);
	sqSolver.Solve();
	return 0;
}
