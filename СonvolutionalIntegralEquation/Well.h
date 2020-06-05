#pragma once
class Well
{
public:
	double Q;

	double P;

	double P0;

	double Time1;

	double Time2;

	double H0;

	double Mu;

	double Rw;

	double K;

	double Kappa;
	
	double Rs;

	double Ksi;

	int N;
	
	double CalculatedP;
	
	double CalculatedQ;

	Well(double Q, double P, double P0, double time1, double time2, double H0, double Mu, double Rw, double K, double kappa, double Rs, double Ksi, int N);
	Well();
};

Well::Well(double Q, double P, double P0, double time1, double time2, double H0, double Mu, double Rw, double K, double kappa, double Rs, double Ksi, int N)
{
	this->Q = 1.0 / (24.0 * 3600.0) * Q;
	this->P = pow(10.0, 6) * P;
	this->P0 = pow(10.0, 6) * P0;
	this->Time1 = 3600.0 * time1;
	this->Time2 = 3600.0 * time2;
	this->H0 = H0;
	this->Mu = pow(10.0, -3) * Mu;
	this->Rw = Rw;
	this->K = pow(10.0, -15) * K;
	this->Kappa = (1.0 / 3600.0) * kappa;
	this->Rs = Rs;
	this->Ksi = Ksi;
	this->N = N;
}

inline Well::Well()
{
}
