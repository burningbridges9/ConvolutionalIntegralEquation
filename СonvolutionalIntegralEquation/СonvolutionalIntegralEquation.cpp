// СonvolutionalIntegralEquation.cpp : Этот файл содержит функцию "main". Здесь начинается и заканчивается выполнение программы.
//

#include <iostream>
#include "IntegralCalculator.h"
#include "Well.h"
#include <cmath>
#include "СonvolutionalIntegralEquation.h"
#include "SequentialSolver.h"
#include "ParallelSolver.h"
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

void TestSequential()
{
	Logger * seqLogger = new Logger("SequentialSolver");

	for (int j = 100; j <= 1000; j += 200)
	{
		auto * wells = new Well[3];
		for (int i = 1; i <= 3; i++)
		{
			Well * well = new Well(5 * i, 5 * i, 3, 5 * (i - 1), 5 * i, 1, 1, 0.1, 10, 4, 0.3, 0, j);
			wells[i - 1] = *well;
		}
		auto * indexes = GetIndexes(wells);

		double * solve;
		double * error = new double;
		double * elapsedOnSlaeSolve = new double;
		double * elapsedOnWholeSolve = new double;
		auto sqSolver = new SequentialSolver(seqLogger, wells, indexes);
		solve = sqSolver->Solve(error, elapsedOnSlaeSolve, elapsedOnWholeSolve);
		for (int i = 0; i < wells[0].N * 3 - 3; i++)
		{
			//solve[i] = wells[0].Time1 + i * step;
			cout << "solve[" << i << "] = " << solve[i] << endl;
		}
		cout << "N = " << j * 3 << endl;
		cout << "error = " << *error << endl;
		cout << "elapsedOnSlaeSolve  = " << *elapsedOnSlaeSolve << endl;
		cout << "elapsedOnWholeSolve = " << *elapsedOnWholeSolve << endl;
	}
}

void TestParallel(Well * wells, int * indexes)
{
	Logger * plLogger = new Logger("ParallelSolver");

	double * solve;
	double * error = new double;
	double * elapsedOnSlaeSolve = new double;
	double * elapsedOnWholeSolve = new double;
	auto plSolver = new ParallelSolver(plLogger, wells, indexes);
	solve = plSolver->Solve(error, elapsedOnSlaeSolve, elapsedOnWholeSolve);
	for (int i = 0; i < wells[0].N * 3 - 3; i++)
	{
		//solve[i] = wells[0].Time1 + i * step;
		cout << "solve[" << i << "] = " << solve[i] << endl;
	}
	cout << "error = " << *error << endl;

}


int main()
{
	
	TestSequential();

	
	return 0;
}


