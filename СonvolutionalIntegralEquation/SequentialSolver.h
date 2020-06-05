#pragma once

#include <iostream>
#include "IntegralCalculator.h"
#include "Well.h"
#include <cmath>
#include "ÑonvolutionalIntegralEquation.h"
#include <iostream>
#include <fstream>
#define _USE_MATH_DEFINES
#define M_PI 3.14159265358979323846
using namespace std;

class SequentialSolver
{
public:
	Well * Wells;
	int * Indexes;
	SequentialSolver(Well * wells, int * indexes);
	void Solve();

private:
	double * GetTimesForPressures(Well * wells);
	void GetPressures(Well * wells);
	double Pressure(Well W, double q, double t);
	double * GetTimesForCoefs(int allN, Well * wells);
	double * GetPressuresForRightSide(int allN, int * indexes, Well * wells);
	double ** PrepareCoefs(double * times, Well* wells, int N);
	double * GaussSeidel(double** A, double* B, int N);
	bool Converge(double * xk, double * xkp, int N);
	void PrintMatrix(double** m, int N);
	void Error(double ** coefs, double * q, double *p, int N);
};

void SequentialSolver::Solve()
{
	GetPressures(Wells);
	int allN = Wells[0].N + Wells[1].N + Wells[2].N - 2;
	auto times = GetTimesForCoefs(allN, Wells);

	auto pressures = GetPressuresForRightSide(allN, Indexes, Wells);
	auto coefs = PrepareCoefs(times, this->Wells, allN);
	PrintMatrix(coefs, allN - 1);
	auto consumptions = GaussSeidel(coefs, pressures, allN-1);
	Wells[0].CalculatedQ = consumptions[Indexes[0] - 2];
	Wells[1].CalculatedQ = consumptions[Indexes[1] - 1];
	Wells[2].CalculatedQ = consumptions[allN - 2];;

	//for (size_t i = 0; i < allN - 1; i++)
	//{
	//	cout << "consumptions[" << i << "] = " << consumptions[i] << endl;
	//	//cout << "pressures[" << i << "] = " << pressures[i] << endl;
	//}
	Error(coefs, consumptions, pressures, allN-1);
}



inline double * SequentialSolver::GetTimesForPressures(Well * wells)
{
	int allN = Wells[0].N + Wells[1].N + Wells[2].N;
	auto * times = new double[allN]; // to avoid the same rows in slae

	double step;
	step = (wells[0].Time2 - wells[0].Time1) / (wells[0].N - 1);

	for (int i = 0; i < allN; i++)
	{
		times[i] = wells[0].Time1 + i * step;
	}
	return times;
}

SequentialSolver::SequentialSolver(Well * wells, int * indexes)
{
	this->Wells = wells;
	this->Indexes = indexes;
}

double * SequentialSolver::GetTimesForCoefs(int allN, Well * wells)
{
	auto * times = new double[allN]; // to avoid same rows in slae

	double step;
	step = (wells[0].Time2 - wells[0].Time1) / (wells[0].N - 1);

	for (int i = 0; i < allN; i++)
	{
		times[i] = wells[0].Time1 + i * step;
		//cout << "times[" << i << "] = " << times[i] << endl;
	}
	return times;
}

double * SequentialSolver::GetPressuresForRightSide(int allN, int * indexes, Well * wells)
{
	auto * eqPressures = new double[allN-1];
	for (int i = 0; i < indexes[0]-1; i++)
	{
		eqPressures[i] = wells[0].CalculatedP - wells[0].P0;
		//cout << "eqPressures[" << i << "] = " << eqPressures[i] << endl;
	}
	for (int i = indexes[0]-1; i <= indexes[1]; i++)
	{
		eqPressures[i] = wells[1].CalculatedP - wells[0].P0;
		//cout << "eqPressures[" << i << "] = " << eqPressures[i] << endl;
	}
	for (int i = indexes[1]; i < indexes[2] - 1; i++)
	{
		eqPressures[i] = wells[2].CalculatedP - wells[0].P0;
		//cout << "eqPressures[" << i << "] = " << eqPressures[i] << endl;
	}
	return eqPressures;
}

double ** SequentialSolver::PrepareCoefs(double * times, Well* wells, int N)
{
	double ** coefs = new double *[N - 1];
	for (int i = 0; i < N - 1; i++)
		coefs[i] = new double[N - 1];

	int n = 1;
	for (int i = 0; i != N - 1; i++)
	{
		int k = 1;
		for (int j = 0; j != N - 1; j++)
		{
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
			k++;
		}
		n++;
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

bool SequentialSolver::Converge(double * xk, double * xkp, int N)
{
	double eps = 0.000001;
	double norm = 0;
	for (int i = 0; i < N; i++)
		norm += (xk[i] - xkp[i]) * (xk[i] - xkp[i]);
	return (sqrt(norm) < eps);
}

inline void SequentialSolver::PrintMatrix(double ** m, int N)
{
	std::ofstream out;          // ïîòîê äëÿ çàïèñè
	out.open("D:\\coefs.txt"); // îêðûâàåì ôàéë äëÿ çàïèñè
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

inline void SequentialSolver::Error(double ** coefs, double * q, double * p, int N)
{
	auto Fmin = pow((Wells[0].Q - Wells[0].CalculatedQ), 2)
		+  pow((Wells[1].Q - Wells[1].CalculatedQ), 2)
		+  pow((Wells[2].Q - Wells[2].CalculatedQ), 2);
	Fmin = Fmin / (pow(Wells[0].Q, 2) + pow(Wells[1].Q, 2) + pow(Wells[2].Q, 2));
}

double SequentialSolver::Pressure(Well W, double q, double t)
{
	double P, arg;
	arg = pow(W.Rs, 2) / (4.0 * W.Kappa * t);
	if (t == 0.0)
	{
		return 0;
	}
	return P = (q * W.Mu) / (4.0 * M_PI * W.H0 * W.K) * (W.Ksi + IntegralCalculator::EIntegral(arg));
}

void SequentialSolver::GetPressures(Well* wells)
{
	double * times = GetTimesForPressures(wells);
	int allN = Wells[0].N + Wells[1].N + Wells[2].N;
	double * pressures = new double[allN];
	/*double * pOne = new double[Wells[0].N];
	double * pTwo = new double[Wells[0].N];
	double * pThree = new double[Wells[0].N];*/
	
	double Q1 = wells[0].Q;
	double Q2 = wells[1].Q;
	double Q3 = wells[2].Q;

	bool entered1 = false;
	bool entered2 = false;

	for (int i = 0; i != allN; i++)
	{
		if (times[i] == 0.0)
		{
			pressures[i] = (0.0 + wells[0].P0);
		}
		else
		{
			if (!entered1 && !entered2)
			{
				pressures[i] = (Pressure(wells[0], Q1, times[i]) + wells[0].P0);
				cout << "pressures[" << i << "] = " << pressures[i] << endl;
			}
			if (times[i] > wells[1].Time1 || i > wells[0].N - 1)
			{
				entered1 = true;
				if (entered1 && !entered2)
				{
					pressures[i] = (Pressure(wells[0], Q1, times[i])
						+ Pressure(wells[1], Q2 - Q1, times[i] - wells[1].Time1)
						+ wells[0].P0);
					cout << "pressures[" << i << "] = " << pressures[i] << endl;
				}
			}
			if ((times[i] > wells[2].Time1) || (i > wells[0].N + wells[1].N - 1))
			{
				entered2 = true;
				if (entered1 && entered2)
				{
					pressures[i] = (Pressure(wells[0], Q1, times[i])
						+ Pressure(wells[1], Q2 - Q1, times[i] - wells[1].Time1)
						+ Pressure(wells[1], Q3 - Q2, times[i] - wells[2].Time1)
						+ wells[0].P0);
					cout << "pressures[" << i << "] = " << pressures[i] << endl;
				}
			}

		}
	}

	wells[0].CalculatedP = pressures[wells[0].N - 1];
	wells[1].CalculatedP = pressures[wells[0].N*2 - 2];
	wells[2].CalculatedP = pressures[wells[0].N*3 - 3];


}