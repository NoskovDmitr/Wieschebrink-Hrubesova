#include <string>
#include <vector>
#include <stdlib.h>
#include "GaloisField.hpp"
#include "matrix.hpp"
#include "grsmatrix.hpp"
#include "randmatrix.hpp"
#include "dimtests.hpp"

void DimTests::makeParameters()
{
	parameters Params;
	int m = 2;
	int k;
	for (int n = 4; n <= 64; n=n*2)
	{
		for (int r = 1; r <= (n-1)/4 ; r++)
		{
			k = (n-1-r)/2 - 1;
			Params.m = m;
			Params.k = k;
			Params.r = r;
			Params.n = n-1;
			params.push_back(Params);
		}
		m++;
	}
}

std::string DimTests::runTest(parameters pam)
{
	Matrix::field = fields[pam.m-2];
	int ok = 0;
	int dim;
	
	for (int j = 0; j < numOfMatrixes; j++)
	{
		Matrix temp(makeMatrix(pam.n, pam.r, pam.k));
		dim = temp.squareStarCodeDim();
		if (dim == 2*pam.k + pam.r -1) ok++;
		//else cout << dim << '\n';
	}
	return "ok= " + to_string(ok) + ";\ntotal= " + to_string(numOfMatrixes) + ";\nrate= " + to_string((double)ok / numOfMatrixes)+ ";";
}

Matrix DimTests::makeMatrix(int n, int r, int k)
{
	Matrix result(n+r,k);
	std::vector<galois::GFSymbol> locators(n);
	std::vector<galois::GFSymbol> multiplicators(n);
	
	for (int i = 0; i < n; i++)
	{
		multiplicators[i] = (galois::GFSymbol) (rand() % Matrix::field->size())+1; 
	}
	
	int ran;
	std::vector<int> pos(Matrix::field->size(),0);
	for (int i = 0; i < Matrix::field->size(); i++)
	{
		pos[i]=i+1;
	}
	for (int i = 0; i < n; i++)
	{
		ran = rand() % pos.size();
		locators[i] = pos[ran]; 
		pos.erase(pos.begin() + ran); 
	}
	GRSMatrix code(k,n,locators, multiplicators);
	RandMatrix rvectors(r,k);
	rvectors.generate();
	
	for (int i = 0; i < n+r; i++)
	{
		for (int j = 0; j < k; j++)
		{
			if ( i < n) result.setXY(i,j, code.getXY(i,j)); //column permutation does not affect dimension
			else result.setXY(i,j, rvectors.getXY(i-n,j));
		}
	}
	
	return result;
}

DimTests::DimTests()
{
	id = "DimTests";
	fileName = "dtest";
	makeFields();
	makeParameters();
}
