#pragma once
#include <iostream>
#include <fstream>
#include <string>
using namespace std;

//Time logger
class Logger
{
public:
	string SolverName;
	string Path;
	const string baseFilePath = "C:\\Users\\Rustam\\Documents\\Visual Studio 2017\\Projects\\ÑonvolutionalIntegralEquation\\ÑonvolutionalIntegralEquation\\Debug\\";
	const string ext = ".txt";
	ofstream writer;

	Logger(string solverName)
	{
		SolverName = solverName;
		Path = baseFilePath + SolverName + ext;

		writer.open(Path, std::ios::app);
		writer << ":::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::";
		writer << "New session" << endl;
		writer.close();
	}

	Logger()
	{

	}

	void Log(string methodName, double elapsedTime)
	{
		string message = methodName + " elapsed time = " + to_string(elapsedTime)+'\n';
		writer.open(Path, std::ios::app);
		writer << message;
		writer.close();
	}

	Logger& operator= (const Logger &logger)
	{
		SolverName = logger.SolverName;

		return *this;
	}
};

