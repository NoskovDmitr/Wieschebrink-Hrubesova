#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <sys/stat.h>
#include <ctime>
#include "GaloisField.hpp"
#include "testin.hpp"

void Testin::makeFields()
{
	//I found these on the Internet somewhere
	unsigned int poly2[3] = {1,1,1};
	unsigned int poly3[4] = {1,1,0,1};
	unsigned int poly4[5] = {1,1,0,0,1};
	unsigned int poly5[6] = {1,1,1,1,0,1};
	unsigned int poly6[7] = {1,1,0,0,1,1,1};
	unsigned int poly7[8] = {1,0,1,1,1,0,0,1};	
	/*
	p(x) = 1x^8+1x^7+0x^6+0x^5+0x^4+0x^3+1x^2+1x^1+1x^0
			1    1    0    0    0    0    1    1    1
	*/
	unsigned int poly8[9] =   {1,1,1,0,0,0,0,1,1};
	unsigned int poly9[10] =  {1,0,0,1,0,1,1,0,0,1};
	unsigned int poly10[11] = {1,0,1,1,0,0,0,0,1,0,1};
	
	fields.resize(9);
	fields[0] = new galois::GaloisField (2,poly2);
	fields[1] = new galois::GaloisField (3,poly3);
	fields[2] = new galois::GaloisField (4,poly4);
	fields[3] = new galois::GaloisField (5,poly5);
	fields[4] = new galois::GaloisField (6,poly6);
	fields[5] = new galois::GaloisField (7,poly7);
	fields[6] = new galois::GaloisField (8,poly8);
	fields[7] = new galois::GaloisField (9,poly9);
	fields[8] = new galois::GaloisField (10,poly10);
	
}

std::string Testin::to_string(int a)
{
	std::stringstream stream;
	stream << a;
	return stream.str();
}

std::string Testin::to_string(double a)
{
	std::stringstream stream;
	stream << std::fixed << a;
	return stream.str();
}

void Testin::writeResults(parameters pam, std::string testResult)
{
	std::string fileOutName= id + "/" + fileName 
		+ to_string(pam.m) +  "_" + to_string(pam.n) 
		+ "_" +to_string(pam.k)+ "_" +to_string(pam.r) + ".txt";
		
	std::cout << fileOutName << '\n';
	std::ofstream fileOut;
	fileOut.open (fileOutName.c_str()); 
	fileOut << "m= " << pam.m << ";\n";
	fileOut << "n= " << pam.n << ";\n";
	fileOut << "k= " << pam.k << ";\n";
	fileOut << "r= " << pam.r << ";\n";
	fileOut << testResult << "\n";
	fileOut.close();
}

bool Testin::fileExists(parameters pam)
{
	std::string name = id + "/" + fileName 
		+ to_string(pam.m) +  "_" + to_string(pam.n) 
		+ "_" +to_string(pam.k)+ "_" +to_string(pam.r) + ".txt";
		
	std::ifstream f(name.c_str());
    if (f.good()) {
        f.close();
        return true;
    } else {
        f.close();
        return false;
    }   
}

void Testin::runTests()
{
	std::string result;
	for (int i = 0; i < params.size(); i++)
	{
		if (fileExists(params[i])) continue;
		std::clock_t start = std::clock();
		result = runTest(params[i]);
		std::clock_t finish = std::clock();
		writeResults(params[i], result);
		std::cout << "Test complete!" << '\n';
		std::cout << "Time elapsed: " << std::fixed << (double)(finish - start) / CLOCKS_PER_SEC << '\n';
	}
	std::cout << "Done " << id << " tests!" << '\n';
}
