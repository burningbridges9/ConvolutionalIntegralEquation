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


void SaveConsumptionsAndTimes(string solverType, string solveType, string nNum, double * consumption, double * times, int N)
{
	string baseFilePath = "C:\\Users\\Rustam\\Documents\\Visual Studio 2017\\Projects\\СonvolutionalIntegralEquation\\СonvolutionalIntegralEquation\\Debug\\";
	string ext = ".txt";

	string fileNameCons = baseFilePath + solverType + solveType + nNum + "CONSUMPTION" + ext;
	string fileNameTimes = baseFilePath + solverType + solveType + nNum + "TIMES" + ext;
	ofstream writerCons;
	ofstream writerTimes;
	writerCons.open(fileNameCons);
	writerTimes.open(fileNameTimes);
	for (size_t i = 0; i < N; i++)
	{
		writerCons << consumption[i] << " ";
		writerTimes << times[i] << " ";
	}
	writerCons.close();
	writerTimes.close();
}

void SaveErrorTimeAndN(string solverType, string solveType, double * error, double * times, int * ns, int N)
{
	string baseFilePath = "C:\\Users\\Rustam\\Documents\\Visual Studio 2017\\Projects\\СonvolutionalIntegralEquation\\СonvolutionalIntegralEquation\\Debug\\";
	string ext = ".txt";

	string fileNameEr = baseFilePath + solverType + solveType + "ERROR" + ext;
	string fileNameTimes = baseFilePath + solverType + solveType + "TIMES" + ext;
	string fileNameNs = baseFilePath + solverType + solveType + "NS" + ext;

	ofstream writerEr;
	ofstream writerTimes;
	ofstream writerNs;
	writerEr.open(fileNameEr);
	writerTimes.open(fileNameTimes);
	writerNs.open(fileNameNs);
	for (size_t i = 0; i < N; i++)
	{
		writerEr << error[i] << " ";
		writerTimes << times[i] << " ";
		writerNs << ns[i] << " ";
	}
	writerEr.close();
	writerTimes.close();
	writerNs.close();
}


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
	cout << "Sequential solver is starting with Gauss - Seidel" << endl;
	int st = 100, fn = 800;
	int arNum = fn / st;
	double * errors = new double[arNum];
	double * times = new double[arNum];
	int * ns = new int[arNum];
	int k = 0;
	for (int j = st; j <= fn; j += 100)
	{
		auto * wells = new Well[3];
		for (int i = 1; i <= 3; i++)
		{
			Well * well = new Well(5 * i, 5 * i, 3, 5 * (i - 1), 5 * i, 1, 1, 0.1, 10, 4, 0.3, 0, j);
			wells[i - 1] = *well;
		}
		auto * indexes = GetIndexes(wells);

		double * solve;
		double * timesForSolve = new double[wells[0].N * 3 - 3];
		double * error = new double;
		double * elapsedOnSlaeSolve = new double;
		double * elapsedOnWholeSolve = new double;
		auto sqSolver = new SequentialSolver(seqLogger, wells, indexes);
		solve = sqSolver->Solve(1, error, elapsedOnSlaeSolve, elapsedOnWholeSolve, timesForSolve);
		cout << "\nN = " << j * 3 << endl;
		cout << "error = " << *error << endl;
		cout << "elapsedOnSlaeSolve  = " << *elapsedOnSlaeSolve << endl;
		cout << "elapsedOnWholeSolve = " << *elapsedOnWholeSolve << endl;
		SaveConsumptionsAndTimes("Sequential", "Seidel", to_string(j * 3), solve, timesForSolve, wells[0].N * 3 - 3);

		errors[k] = *error;
		times[k] = *elapsedOnWholeSolve;
		ns[k] = j * 3;
		k++;
	}
	SaveErrorTimeAndN("Sequential", "Seidel", errors, times, ns, arNum);

	k = 0;
	for (int j = st; j <= fn; j += 100)
	{
		auto * wells = new Well[3];
		for (int i = 1; i <= 3; i++)
		{
			Well * well = new Well(5 * i, 5 * i, 3, 5 * (i - 1), 5 * i, 1, 1, 0.1, 10, 4, 0.3, 0, j);
			wells[i - 1] = *well;
		}
		auto * indexes = GetIndexes(wells);

		double * solve;
		double * timesForSolve = new double[wells[0].N * 3 - 3];
		double * error = new double;
		double * elapsedOnSlaeSolve = new double;
		double * elapsedOnWholeSolve = new double;
		auto sqSolver = new SequentialSolver(seqLogger, wells, indexes);
		solve = sqSolver->Solve(2, error, elapsedOnSlaeSolve, elapsedOnWholeSolve, timesForSolve);
		cout << "\nN = " << j * 3 << endl;
		cout << "error = " << *error << endl;
		cout << "elapsedOnSlaeSolve  = " << *elapsedOnSlaeSolve << endl;
		cout << "elapsedOnWholeSolve = " << *elapsedOnWholeSolve << endl;
		SaveConsumptionsAndTimes("Sequential", "Reverse", to_string(j * 3), solve, timesForSolve, wells[0].N * 3 - 3);

		errors[k] = *error;
		times[k] = *elapsedOnWholeSolve;
		ns[k] = j * 3;
		k++;
	}
	SaveErrorTimeAndN("Sequential", "Reverse", errors, times, ns, arNum);
}

void TestParallel()
{
	Logger * plLogger = new Logger("ParallelSolver");
	cout << endl << "Parallel solver is starting with Gauss - Seidel" << endl;
	int st = 100, fn = 800;
	int arNum = fn / st;
	double * errors = new double[arNum];
	double * times = new double[arNum];
	int * ns = new int[arNum];
	int k = 0;
	for (int j = fn; j <= 800; j += 100)
	{
		auto * wells = new Well[3];
		for (int i = 1; i <= 3; i++)
		{
			Well * well = new Well(5 * i, 5 * i, 3, 5 * (i - 1), 5 * i, 1, 1, 0.1, 10, 4, 0.3, 0, j);
			wells[i - 1] = *well;
		}
		auto * indexes = GetIndexes(wells);

		double * solve;
		double * timesForSolve = new double[wells[0].N*3-3];
		double * error = new double;
		double * elapsedOnSlaeSolve = new double;
		double * elapsedOnWholeSolve = new double;
		auto plSolver = new ParallelSolver(plLogger, wells, indexes);
		solve = plSolver->Solve(1, error, elapsedOnSlaeSolve, elapsedOnWholeSolve, timesForSolve);
#pragma region dbg
		//for (int i = 0; i < wells[0].N * 3 - 3; i++)
//{
//	//solve[i] = wells[0].Time1 + i * step;
//	cout << "solve[" << i << "] = " << solve[i] << endl;
//}  
#pragma endregion
		cout << "\nN = " << j * 3 << endl;
		cout << "error = " << *error << endl;
		cout << "elapsedOnSlaeSolve  = " << *elapsedOnSlaeSolve << endl;
		cout << "elapsedOnWholeSolve  = " << *elapsedOnWholeSolve << endl;
		//SaveConsumptionsAndTimes("Parallel", "Seidel", to_string(j * 3), solve, timesForSolve, wells[0].N * 3 - 3);
		
		errors[k] = *error;
		times[k] = *elapsedOnWholeSolve;
		ns[k] = j * 3;
		k++;
		//delete solve, error, elapsedOnSlaeSolve, elapsedOnWholeSolve, indexes, wells;
	}
	SaveErrorTimeAndN("Parallel", "Seidel", errors, times, ns, arNum);

	k = 0;
	for (int j = st; j <= fn; j += 100)
	{
		auto * wells = new Well[3];
		for (int i = 1; i <= 3; i++)
		{
			Well * well = new Well(5 * i, 5 * i, 3, 5 * (i - 1), 5 * i, 1, 1, 0.1, 10, 4, 0.3, 0, j);
			wells[i - 1] = *well;
		}
		auto * indexes = GetIndexes(wells);

		double * solve;
		double * timesForSolve = new double[wells[0].N * 3 - 3];
		double * error = new double;
		double * elapsedOnSlaeSolve = new double;
		double * elapsedOnWholeSolve = new double;
		auto plSolver = new ParallelSolver(plLogger, wells, indexes);
		solve = plSolver->Solve(1, error, elapsedOnSlaeSolve, elapsedOnWholeSolve, timesForSolve);
#pragma region dbg
		//for (int i = 0; i < wells[0].N * 3 - 3; i++)
//{
//	//solve[i] = wells[0].Time1 + i * step;
//	cout << "solve[" << i << "] = " << solve[i] << endl;
//}  
#pragma endregion
		cout << "\nN = " << j * 3 << endl;
		cout << "error = " << *error << endl;
		cout << "elapsedOnSlaeSolve  = " << *elapsedOnSlaeSolve << endl;
		cout << "elapsedOnWholeSolve  = " << *elapsedOnWholeSolve << endl;
		//SaveConsumptionsAndTimes("Parallel", "Reverse", to_string(j * 3), solve, timesForSolve, wells[0].N * 3 - 3);

		errors[k] = *error;
		times[k] = *elapsedOnWholeSolve;
		ns[k] = j * 3;
		k++;
		//delete solve, error, elapsedOnSlaeSolve, elapsedOnWholeSolve, indexes, wells;
	}
	SaveErrorTimeAndN("Parallel", "Reverse", errors, times, ns, arNum);

	delete plLogger;
}

int main()
{
	//TestSequential();
	TestParallel();
	cout << "\nDone" << endl;
	cin.get();
	return 0;
}


