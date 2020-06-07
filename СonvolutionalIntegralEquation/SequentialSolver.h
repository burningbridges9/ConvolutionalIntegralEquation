#pragma once

#include "IntegralCalculator.h"
#include "Well.h"
#include <cmath>
#include "ÑonvolutionalIntegralEquation.h"
#include <iostream>
#include <fstream>
#include <ctime>
#include <omp.h>
#include "Solver.h"
#define _USE_MATH_DEFINES
#define M_PI 3.14159265358979323846
using namespace std;

class SequentialSolver : public Solver
{
public:
	SequentialSolver(Logger *log, Well * wells, int * indexes);

protected:
	double * GetTimesForPressures(Well * wells);
	double * GetTimesForCoefs(int allN, Well * wells);
	double * GetPressuresForRightSide(int allN, int * indexes, Well * wells);
	double ** PrepareCoefs(double * times, Well* wells, int N);
	double * GaussSeidel(double** A, double* B, int N);
	bool Converge(double * xk, double * xkp, int N);
	double * GaussReverse(double** A, double* B, int N);
	Logger logger;
};


double * SequentialSolver::GetTimesForPressures(Well * wells)
{
	double startTime = omp_get_wtime();
	int allN = Wells[0].N + Wells[1].N + Wells[2].N;
	auto * times = new double[allN]; // to avoid the same rows in slae
	double step = (wells[0].Time2 - wells[0].Time1) / (wells[0].N - 1);
	for (int i = 0; i < allN; i++)
	{
		times[i] = wells[0].Time1 + i * step;
	}
	double endtime = omp_get_wtime();
	//printf("SEQUENTIAL: GetTimesForPressures elapsed time = %f \n", (endtime - startTime) / (CLOCKS_PER_SEC / 1000));
	return times;
}

SequentialSolver::SequentialSolver(Logger *log, Well * wells, int * indexes) : Solver(wells, indexes)
{
	logger = *log;
}

double * SequentialSolver::GetTimesForCoefs(int allN, Well * wells)
{
	double startTime = omp_get_wtime();
	auto * times = new double[allN]; // to avoid same rows in slae
	double step;
	step = (wells[0].Time2 - wells[0].Time1) / (wells[0].N - 1);

	for (int i = 0; i < allN; i++)
	{
		times[i] = wells[0].Time1 + i * step;
		//cout << "times[" << i << "] = " << times[i] << endl;
	}
	double endtime = omp_get_wtime();
	printf("SEQUENTIAL: GetTimesForPressures elapsed time = %f \n", (endtime - startTime) / (CLOCKS_PER_SEC / 1000));
	return times;
}

double * SequentialSolver::GetPressuresForRightSide(int allN, int * indexes, Well * wells)
{
	double startTime = omp_get_wtime();
	auto * eqPressures = new double[allN - 1];
	for (int i = 0; i < indexes[0] - 1; i++)
	{
		eqPressures[i] = wells[0].CalculatedP - wells[0].P0;
		//cout << "eqPressures[" << i << "] = " << eqPressures[i] << endl;
	}
	for (int i = indexes[0] - 1; i <= indexes[1]; i++)
	{
		eqPressures[i] = wells[1].CalculatedP - wells[0].P0;
		//cout << "eqPressures[" << i << "] = " << eqPressures[i] << endl;
	}
	for (int i = indexes[1]; i < indexes[2] - 1; i++)
	{
		eqPressures[i] = wells[2].CalculatedP - wells[0].P0;
		//cout << "eqPressures[" << i << "] = " << eqPressures[i] << endl;
	}
	double endtime = omp_get_wtime();
	printf("SEQUENTIAL: GetPressuresForRightSide elapsed time = %f \n", (endtime - startTime) / (CLOCKS_PER_SEC / 1000));
	return eqPressures;
}

double ** SequentialSolver::PrepareCoefs(double * times, Well* wells, int N)
{
	double ** coefs = new double *[N - 1];
	for (int i = 0; i < N - 1; i++)
		coefs[i] = new double[N - 1];


	for (int i = 0; i != N - 1; i++)
	{
		int n = i + 1;
		for (int j = 0; j != N - 1; j++)
		{
			int k = j + 1;
			if (k <= n)
			{
				double E1, E2, arg1, arg2;

				arg1 = pow(wells[0].Rs, 2) * 1.0 / (4 * wells[0].Kappa * (times[n] - times[k - 1]));
				E1 = arg1 < 1 ? IntegralCalculator::PolyApproxExpIntegral2(arg1) : IntegralCalculator::PolyApproxExpIntegral1(arg1);

				if (k == n)
				{
					E1 = E1 + wells[0].Ksi;
					arg2 = 0.0;
				}
				else
				{
					arg2 = pow(wells[0].Rs, 2) * 1.0 / (4 * wells[0].Kappa * (times[n] - times[k]));
				}

				E2 = (arg2 < 1) && (arg2 > 0.0) ? IntegralCalculator::PolyApproxExpIntegral2(arg2) : 0.0;
				auto coefBefore = wells[0].Mu / (4.0 * M_PI * wells[0].K * wells[0].H0); // 39788735772.973831
				coefs[i][j] = coefBefore * (E1 - E2);
			}
			else
			{
				coefs[i][j] = 0.0;
			}
			//cout << "coefs[" << i << "][" << j << "] = " << coefs[i][j] << endl;
		}
	}

	return coefs;
}

double * SequentialSolver::GaussSeidel(double** A, double* B, int N)
{
	double * prev = new double[N];
	double * x = new double[N];
	for (size_t i = 0; i < N; i++)
	{
		x[i] = 0;
	}
	int m = 1;
	do
	{
		for (int i = 0; i < N; i++)
			prev[i] = x[i];
		for (int i = 0; i < N; i++)
		{
			double var = 0;
			for (int j = 0; j < i; j++)
				var += (A[i][j] * x[j]);
			for (int j = i + 1; j < N; j++)
				var += (A[i][j] * prev[j]);
			x[i] = (B[i] - var) / A[i][i];
		}
		m++;
	} while (!Converge(x, prev, N));
	return x;
}

double * SequentialSolver::GaussReverse(double** A, double* B, int N)
{
	double * x = new double[N];
	for (size_t i = 0; i < N; i++)
	{
		auto sum = 0.0;
		for (size_t j = 0; j < i; j++)
		{
			sum += A[i][j] * x[j];
		}
		x[i] = (B[i] - sum) / A[i][i];
	}
	for (size_t i = 0; i < N; i++)
	{
		cout << "x[" << i << "] = " << x[i] << endl;
	}
	return x;
}

bool SequentialSolver::Converge(double * xk, double * xkp, int N)
{
	double eps = 0.000001;
	double norm = 0;
	for (int i = 0; i < N; i++)
		norm += (xk[i] - xkp[i]) * (xk[i] - xkp[i]);
	return (sqrt(norm) < eps);
}

