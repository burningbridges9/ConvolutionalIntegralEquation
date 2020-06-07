#pragma once

#include <iostream>
#include "IntegralCalculator.h"
#include "Well.h"
#include <cmath>
#include "СonvolutionalIntegralEquation.h"
#include <iostream>
#include <fstream>
#include <ctime>
#include "Logger.h"
#include <omp.h>
#define _USE_MATH_DEFINES
#define M_PI 3.14159265358979323846
using namespace std;

class Solver
{
public:
	Well * Wells;
	int * Indexes;
	Solver(Well * wells, int * indexes);
	double * Solve(double * error, double * elapsedOnSlaeSolve, double * elapsedOnWholeSolve);

protected:
	virtual double * GetTimesForPressures(Well * wells) = 0;
	void GetPressures(Well * wells);
	double Pressure(Well W, double q, double t);
	virtual double * GetTimesForCoefs(int allN, Well * wells) = 0;
	virtual double * GetPressuresForRightSide(int allN, int * indexes, Well * wells) = 0;
	virtual double ** PrepareCoefs(double * times, Well* wells, int N) = 0;
	virtual double * GaussSeidel(double** A, double* B, int N) = 0;
	virtual double * GaussReverse(double** A, double* B, int N) = 0;
	virtual bool Converge(double * xk, double * xkp, int N) = 0;
	void PrintMatrix(double** m, int N);
	double Error(double ** coefs, double * q, double *p, int N);
};

double * Solver::Solve(double * error, double * elapsedOnSlaeSolve, double * elapsedOnWholeSolve)
{
	double startTimeWhole = omp_get_wtime();
	GetPressures(Wells);

	double startTimeSlae = omp_get_wtime();
	int allN = Wells[0].N + Wells[1].N + Wells[2].N - 2;
	auto times = GetTimesForCoefs(allN, Wells);
	auto pressures = GetPressuresForRightSide(allN, Indexes, Wells);
	auto coefs = PrepareCoefs(times, this->Wells, allN);
	PrintMatrix(coefs, allN - 1);
	auto consumptions = GaussSeidel(coefs, pressures, allN - 1);
	 double endtimeSlae = omp_get_wtime();
	 *elapsedOnSlaeSolve = (endtimeSlae - startTimeSlae) / (CLOCKS_PER_SEC / 1000);
	Wells[0].CalculatedQ = consumptions[Indexes[0] - 2];
	Wells[1].CalculatedQ = consumptions[Indexes[1] - 1];
	Wells[2].CalculatedQ = consumptions[allN - 2];;
#pragma region dbg 

	//for (size_t i = 0; i < allN - 1; i++)
	//{
	//	cout << "consumptions[" << i << "] = " << consumptions[i] << endl;
	//	//cout << "pressures[" << i << "] = " << pressures[i] << endl;
	//}  
#pragma endregion

	*error = Error(coefs, consumptions, pressures, allN - 1);
	double endtimeWhole = omp_get_wtime(); 
	*elapsedOnWholeSolve = (endtimeWhole - startTimeWhole) / (CLOCKS_PER_SEC / 1000);
	return consumptions;
}


Solver::Solver(Well * wells, int * indexes)
{
	this->Wells = wells;
	this->Indexes = indexes;
}

inline void Solver::PrintMatrix(double ** m, int N)
{
	std::ofstream out;          // поток для записи
	out.open("D:\\coefs.txt"); // окрываем файл для записи
	if (out.is_open())
	{
		for (size_t i = 0; i < N; i++)
		{
			for (size_t j = 0; j < N; j++)
			{
				out << m[i][j] << " ";
			}
			out << endl;
		}
	}
	out.close();
}

inline double Solver::Error(double ** coefs, double * q, double * p, int N)
{
	auto Fmin = pow((Wells[0].Q - Wells[0].CalculatedQ), 2)
		+ pow((Wells[1].Q - Wells[1].CalculatedQ), 2)
		+ pow((Wells[2].Q - Wells[2].CalculatedQ), 2);
	Fmin = Fmin / (pow(Wells[0].Q, 2) + pow(Wells[1].Q, 2) + pow(Wells[2].Q, 2));
	return Fmin;
}

double Solver::Pressure(Well W, double q, double t)
{
	auto arg = pow(W.Rs, 2) / (4.0 * W.Kappa * t);
	return t == 0.0 ? 0 : (q * W.Mu) / (4.0 * M_PI * W.H0 * W.K) * (W.Ksi + IntegralCalculator::EIntegral(arg));
}

void Solver::GetPressures(Well* wells)
{
	double * times = GetTimesForPressures(wells);
	double Q1 = wells[0].Q;
	double Q2 = wells[1].Q;
	double Q3 = wells[2].Q;
	wells[0].CalculatedP = (Pressure(wells[0], Q1, times[wells[0].N - 1]) + wells[0].P0);
	wells[1].CalculatedP = (Pressure(wells[0], Q1, times[wells[0].N * 2 - 2])
		+ Pressure(wells[1], Q2 - Q1, times[wells[0].N * 2 - 2] - wells[1].Time1)
		+ wells[0].P0);
	wells[2].CalculatedP = (Pressure(wells[0], Q1, times[wells[0].N * 3 - 3])
		+ Pressure(wells[1], Q2 - Q1, times[wells[0].N * 3 - 3] - wells[1].Time1)
		+ Pressure(wells[1], Q3 - Q2, times[wells[0].N * 3 - 3] - wells[2].Time1)
		+ wells[0].P0);

}