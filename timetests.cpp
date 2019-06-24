#include <string>
#include <vector>
#include <ctime>
#include <iostream>
#include <stdlib.h>
#include "matrix.hpp"
#include "GaloisField.hpp"
#include "randmatrix.hpp"
#include "regmatrix.hpp"
#include "permmatrix.hpp"
#include "attacker.hpp"
#include "cryptosystem.hpp"
#include "polynomial.hpp"
#include "grsmatrix.hpp"
#include "timetests.hpp"


void TimeTests::makeParameters()
{
	//I calculated these using SageMath as Wieschebrink suggested in [12]
	parameters Pars;
	
	Pars.m = 7; // (q = 128)
	Pars.n = 95;
	Pars.k = 61;
	Pars.r = 12;
	params.push_back(Pars);

	Pars.m = 7; // (q = 128)
	Pars.n = 127;
	Pars.k = 79;
	Pars.r = 17;
	params.push_back(Pars);
	
	Pars.m = 8; // (q = 256)	
	Pars.n = 191;
	Pars.k = 119;
	Pars.r = 26;
	params.push_back(Pars);

	Pars.m = 8; // (q = 256)
	Pars.n = 255;
	Pars.k = 161;
	Pars.r = 34;
	params.push_back(Pars);

	Pars.m = 9; // (q = 512)
	Pars.n = 383;
	Pars.k = 243;
	Pars.r = 51;
	params.push_back(Pars);

	Pars.m = 9; // (q = 512)
	Pars.n = 511;
	Pars.k = 325;
	Pars.r = 68;
	params.push_back(Pars);

	Pars.m = 10; // (q = 1024)
	Pars.n = 767;
	Pars.k = 475;
	Pars.r = 106;
	params.push_back(Pars);
	
	Pars.m = 10; // (q = 1024)
	Pars.n = 1023;
	Pars.k = 639;
	Pars.r = 140;
	params.push_back(Pars);
}

TimeTests::TimeTests()
{
	id = "TimeTests";
	fileName = "ttest";
	makeFields();
	makeParameters();
}

std::string TimeTests::runTest(parameters pam)
{
	Matrix::field = fields[pam.m-2];
	GRSMatrix::field = fields[pam.m-2];
	Polynomial::field = fields[pam.m-2];
	Cryptosystem::field = fields[pam.m-2];
	Attacker::field = fields[pam.m-2];
	double time;
	bool res=false;
	
	Attacker eve(makeSystem(pam));

	std::clock_t start = std::clock();
	for (int i = 0; i < numOfRuns; i++)
	{
		std::clock_t wastedStart = std::clock();
		eve.newCryptosystem();
		std::clock_t wastedFinish = std::clock();
		time = time - ((double)(wastedFinish - wastedStart) / CLOCKS_PER_SEC);
		res = eve.ripntear();
		std::cout << i << ':';
		if (! res) std::cout << "error" << '\n';
		else std::cout << "ok" << '\n';
	}
	std::clock_t finish = std::clock();
	time = time + (double)(finish - start) / (CLOCKS_PER_SEC * numOfRuns);
	return "avgTime= " + to_string(time) + ";";	
}

Cryptosystem TimeTests::makeSystem(parameters pam)
{
	std::vector<galois::GFSymbol> locators(pam.n);
	std::vector<galois::GFSymbol> multiplicators(pam.n);
	
	for (int i = 0; i < pam.n; i++)
	{
		multiplicators[i] = (galois::GFSymbol) (rand() % Matrix::field->size())+1; //nonzero, not necessarily distinct
	}
	
	int ran;
	std::vector<int> pos(Matrix::field->size(),0);
	for (int i = 0; i < Matrix::field->size(); i++)
	{
		pos[i]=i+1;
	}
	for (int i = 0; i < pam.n; i++)
	{
		ran = rand() % pos.size();
		locators[i] = pos[ran]; 
		pos.erase(pos.begin() + ran); 
	}
	GRSMatrix code(pam.n-pam.k,pam.n,locators, multiplicators);
	Cryptosystem sys(code,pam.r);
	return sys;
}
