#pragma once

#include <iostream>
#include "IntegralCalculator.h"
#include "Well.h"
#include <cmath>
#include "ÑonvolutionalIntegralEquation.h"
#include <iostream>
#include <fstream>
#include <ctime>
#include <omp.h>
#include "Solver.h"
#include "Logger.h"
#define _USE_MATH_DEFINES
#define M_PI 3.14159265358979323846
using namespace std;

class ParallelSolver : public Solver
{
public:
	ParallelSolver(Logger *log, Well * wells, int * indexes);
	~ParallelSolver();
protected:
	double * GetTimesForPressures(Well * wells);
	double * GetTimesForCoefs(int allN, Well * wells);
	double * GetPressuresForRightSide(int allN, int * indexes, Well * wells);
	double ** PrepareCoefs(double * times, Well* wells, int N);
	double * GaussSeidel(double** A, double* B, int N);
	double * GaussReverse(double** A, double* B, int N);
	bool Converge(double * xk, double * xkp, int N);
};




double * ParallelSolver::GetTimesForPressures(Well * wells)
{
	double startTime = omp_get_wtime();
	int allN = Wells[0].N + Wells[1].N + Wells[2].N;
	auto * times = new double[allN]; // to avoid the same rows in slae

	double step = (wells[0].Time2 - wells[0].Time1) / (wells[0].N - 1);
#pragma omp parallel for schedule(guided) 
	for (int i = 0; i < allN; i++)
	{
		times[i] = wells[0].Time1 + i * step;
		//cout << "times[" << i << "] = " << times[i] << endl;
		//cout << "thread = " << omp_get_thread_num() << endl;
	}
	double endtime = omp_get_wtime();
	//printf("PARALLEL: GetTimesForPressures elapsed time = %f \n", (endtime - startTime) / (CLOCKS_PER_SEC / 1000));
	/*for (int i = 0; i < allN; i++)
	{
		cout << "times[" << i << "] = " << times[i] << endl;
	}*/

	logger->Log("GetTimesForPressures", (endtime - startTime) / (CLOCKS_PER_SEC / 1000));
	return times;
}

ParallelSolver::ParallelSolver(Logger *log, Well * wells, int * indexes) : Solver(log, wells, indexes)
{
	string msg = "Starting ParallelSolver N = " + to_string(wells[0].N * 3) + '\n';
	this->logger->Log(msg);
}

ParallelSolver::~ParallelSolver()
{
	this->logger->Log("Finishing ParallelSolver");
}

double * ParallelSolver::GetTimesForCoefs(int allN, Well * wells)
{
	double startTime = omp_get_wtime();
	auto * times = new double[allN]; // to avoid same rows in slae
	double step = (wells[0].Time2 - wells[0].Time1) / (wells[0].N - 1);

#pragma omp parallel for schedule(guided) 
	for (int i = 0; i < allN; i++)
	{
		times[i] = wells[0].Time1 + i * step;
		//cout << "times[" << i << "] = " << times[i] << endl;
	}
	double endtime = omp_get_wtime();
	//printf("PARALLEL: GetTimesForCoefs elapsed time = %f \n", (endtime - startTime) / (CLOCKS_PER_SEC / 1000)); 
#pragma region dbg
	/*for (int i = 0; i < allN; i++)
{
	cout << "times[" << i << "] = " << times[i] << endl;
}*/
#pragma endregion

	logger->Log("GetTimesForCoefs", (endtime - startTime) / (CLOCKS_PER_SEC / 1000));
	return times;
}

double * ParallelSolver::GetPressuresForRightSide(int allN, int * indexes, Well * wells)
{
	double startTime = omp_get_wtime();
	auto * eqPressures = new double[allN - 1];
#pragma omp parallel for schedule(guided) 
		for (int i = 0; i < indexes[0] - 1; i++)
		{
			eqPressures[i] = wells[0].CalculatedP - wells[0].P0;
		}

#pragma omp parallel for schedule(guided) 
		for (int i = indexes[0] - 1; i <= indexes[1]; i++)
		{
			eqPressures[i] = wells[1].CalculatedP - wells[0].P0;
		}

#pragma omp parallel for schedule(guided) 
		for (int i = indexes[1]; i < indexes[2] - 1; i++)
		{
			eqPressures[i] = wells[2].CalculatedP - wells[0].P0;
		}

	double endtime = omp_get_wtime();
	//printf("PARALLEL: GetPressuresForRightSide elapsed time = %f \n", (endtime - startTime) / (CLOCKS_PER_SEC / 1000));
#pragma region dbg
	/*for (int i = 0; i < allN - 1; i++)
{
	cout << "eqPressures[" << i << "] = " << eqPressures[i] << endl;
}*/
#pragma endregion

	logger->Log("GetPressuresForRightSide", (endtime - startTime) / (CLOCKS_PER_SEC / 1000));
	return eqPressures;
}

double ** ParallelSolver::PrepareCoefs(double * times, Well* wells, int N)
{
	double ** coefs = new double *[N - 1];

#pragma omp parallel for schedule(guided) 
	for (int i = 0; i < N - 1; i++)
		coefs[i] = new double[N - 1];

	int n, k;
#pragma region Collapse is not working
	//#pragma omp parallel for collapse(2)
		//for (int i = 0; i < N - 1; i++)
		//{
		//	n = i + 1;
		//	for (int j = 0; j < N - 1; j++)
		//	{
		//		k = j + 1;
		//		if (k <= n)
		//		{
		//			double E1, E2, arg1, arg2;

		//			arg1 = pow(wells[0].Rs, 2) * 1.0 / (4 * wells[0].Kappa * (times[n] - times[k - 1]));
		//			E1 = arg1 < 1 ? IntegralCalculator::PolyApproxExpIntegral2(arg1) : IntegralCalculator::PolyApproxExpIntegral1(arg1);

		//			if (k == n)
		//			{
		//				E1 = E1 + wells[0].Ksi;
		//				arg2 = 0.0;
		//			}
		//			else
		//			{
		//				arg2 = pow(wells[0].Rs, 2) * 1.0 / (4 * wells[0].Kappa * (times[n] - times[k]));
		//			}

		//			E2 = (arg2 < 1) && (arg2 > 0.0) ? IntegralCalculator::PolyApproxExpIntegral2(arg2) : 0.0;
		//			auto coefBefore = wells[0].Mu / (4.0 * M_PI * wells[0].K * wells[0].H0); // 39788735772.973831
		//			coefs[i][j] = coefBefore * (E1 - E2);
		//		}
		//		else
		//		{
		//			coefs[i][j] = 0.0;
		//		}
		//	}
		//}
#pragma endregion


	double startTime = omp_get_wtime();
	//(dynamic, N/omp_get_num_threads()) 
#pragma region test only
#pragma omp parallel for schedule(guided) private (n,k) // dynamic and etc // debug & release
	for (int ij = 0; ij < (N - 1)*(N - 1); ij++)
	{
		int i = ij / (N - 1);
		int j = ij % (N - 1);

		n = i + 1;
		k = j + 1;
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
	}

#pragma endregion


	double endtime = omp_get_wtime();
	logger->Log("PrepareCoefs", (endtime - startTime) / (CLOCKS_PER_SEC / 1000));
	return coefs;
}

double * ParallelSolver::GaussSeidel(double** A, double* B, int N)
{
	double startTime = omp_get_wtime();
	double * prev = new double[N];
	double * x = new double[N];

#pragma omp parallel for schedule(guided)
	for (int i = 0; i < N; i++)
		x[i] = 0;

	do
	{
#pragma omp parallel for schedule(guided)
		for (int i = 0; i < N; i++)
			prev[i] = x[i];

#pragma omp parallel for schedule(guided)
		for (int i = 0; i < N; i++)
		{
			double sum = 0;
			for (int j = 0; j < i; j++)
				sum += (A[i][j] * x[j]);

			for (int j = i + 1; j < N; j++)
				sum += (A[i][j] * prev[j]);
			x[i] = (B[i] - sum) / A[i][i];
		}

#pragma region second choice
//		for (int i = 0; i < N; i++)
//		{
//			double sum = 0;
//#pragma omp parallel for schedule(guided) reduction(+:sum)
//			for (int j = 0; j < i; j++)
//				sum += (A[i][j] * x[j]);
//
//#pragma omp parallel for schedule(guided) reduction(+:sum)
//			for (int j = i + 1; j < N; j++)
//				sum += (A[i][j] * prev[j]);
//
//			x[i] = (B[i] - sum) / A[i][i];
//		}
//
//		// slower 2x
#pragma endregion

	} while (!Converge(x, prev, N));

	double endtime = omp_get_wtime();
	logger->Log("GaussSeidel", (endtime - startTime) / (CLOCKS_PER_SEC / 1000));
	return x;
}

bool ParallelSolver::Converge(double * xk, double * xkp, int N)
{
	double eps = 0.000001;
	double norm = 0;
#pragma omp parallel for schedule(guided) reduction(+:norm)
	for (int i = 0; i < N; i++)
		norm += (xk[i] - xkp[i]) * (xk[i] - xkp[i]);
	return (sqrt(norm) < eps);
}

double * ParallelSolver::GaussReverse(double** A, double* B, int N)
{
	double startTime = omp_get_wtime();
	double * x = new double[N];
	for (int i = 0; i < N; i++)
	{
		double sum = 0.0;
#pragma omp parallel for schedule(guided) reduction(+:sum)
		for (int j = 0; j < i; j++)
		{
			sum += A[i][j] * x[j];
		}
		x[i] = (B[i] - sum) / A[i][i];
	}
	double endtime = omp_get_wtime();
	logger->Log("GaussReverse", (endtime - startTime) / (CLOCKS_PER_SEC / 1000));
	return x;
}
