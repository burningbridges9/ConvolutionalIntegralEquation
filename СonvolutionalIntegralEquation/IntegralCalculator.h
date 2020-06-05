#pragma once
#include <cmath>
#define _USE_MATH_DEFINES
#define M_E 2.7182818284590452354
using namespace std;

class IntegralCalculator
{
public:
	static double PolyApproxExpIntegral1(double arg);
	static double PolyApproxExpIntegral2(double arg);
	static double EIntegral(double arg);
};

double IntegralCalculator::PolyApproxExpIntegral1(double arg) //1 < ARG < INF
{
	double a1, a2, a3, a4, b1, b2, b3, b4, sum1, sum2;
	a1 = 8.5733287401;
	a2 = 18.0590169730;
	a3 = 8.6347608925;
	a4 = 0.2677737343;
	b1 = 9.5733223454;
	b2 = 25.6329561486;
	b3 = 21.0996530827;
	b4 = 3.9584969228;
	sum1 = pow(arg, 4) + a1 * pow(arg, 3) + a2 * pow(arg, 2) + a3 * arg + a4;
	sum2 = pow(arg, 4) + b1 * pow(arg, 3) + b2 * pow(arg, 2) + b3 * arg + b4;
	return (sum1 / (sum2 * arg * pow(M_E, arg)));
}

double IntegralCalculator::PolyApproxExpIntegral2(double arg) // 0 < ARG < 1
{
	double a0, a1, a2, a3, a4, a5, sum, E;
	a0 = -0.57721566;
	a1 = 0.99999193;
	a2 = -0.24991055;
	a3 = 0.05519968;
	a4 = -0.00976004;
	a5 = 0.00107857;
	sum = a0 + a1 * arg + a2 * pow(arg, 2) + a3 * pow(arg, 3) + a4 * pow(arg, 4) + a5 * pow(arg, 5);
	E = sum - log(arg);
	return E;
}

double IntegralCalculator::EIntegral(double arg)
{
	return  arg > 1 ? PolyApproxExpIntegral1(arg) : PolyApproxExpIntegral2(arg);
}