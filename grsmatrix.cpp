#include <vector>
#include "GaloisField.hpp"
#include "grsmatrix.hpp"
#include "regmatrix.hpp"
#include <stdlib.h>

class Matrix;
		
GRSMatrix::GRSMatrix(int rows, int cols, std::vector<galois::GFSymbol> locs,std::vector<galois::GFSymbol> multips)
{
	locators = locs;
	multiplicators = multips;
	
	n = cols;
	m = rows;
	
	array.clear();
	array.resize(this->m*this->n);
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < m; j++)
		{
			this->setXY(i,j,field->mul(multiplicators[i],field->exp(locators[i],j)));
		}
	}
}

GRSMatrix::GRSMatrix(int rows, int cols)
{
	n = cols;
	m = rows;
	
	array.clear();
	array.resize(this->m*this->n);
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < m; j++)
		{
			this->setXY(i,j,(galois::GFSymbol)0);
		}
	}	
}

GRSMatrix::GRSMatrix()
{
	
}

Matrix GRSMatrix::encode(Matrix information)
{
	Matrix infWord(n-k,1);
	for (int i = 0; i < n-k; i++)
	{
		infWord.setXY(i,0,information.getXY(i,0));
	}
	Matrix gen(findKer().transpose());
	return Matrix::multiplication(infWord, gen).transpose(); //should have made encode output compatible with decode input
}

Matrix GRSMatrix::decode(Matrix encoded)
{
	//syndrome
	Matrix syndrome(Matrix::multiplication(*this,encoded));
	//Euclid
	Polynomial synd(syndrome.getNumRows());
	for (int i = 0; i < syndrome.getNumRows(); i++)
	{
		synd.setA(i,syndrome.getXY(0,i));
	}
	Polynomial lambda(2);
	Polynomial gamma(2);
	synd.EuclidForDecoding(syndrome.getNumRows(),lambda, gamma);
	
	//Forney
	Matrix result(1,n);
	for (int i = 0; i < n; i++)
	{
		if (lambda.eval(field->inverse(locators[i])) == 0) 
		{
			galois::GFSymbol first = field->sub((galois::GFSymbol)0,field->div(locators[i],multiplicators[i]));
			galois::GFSymbol gammaEval = gamma.eval(field->inverse(locators[i]));
			Polynomial der = lambda.derivative();
			galois::GFSymbol lambdaEval = der.eval(field->inverse(locators[i]));
			result.setXY(0,i,field->mul(first, field->div(gammaEval,lambdaEval)));
		}
		else result.setXY(0,i,(galois::GFSymbol)0);
	}
	return result; 
}

void GRSMatrix::generate()
{
	//generate multiplicators (nonzero)
	for (int i = 0; i < n; i++)
	{
		multiplicators[i] = (galois::GFSymbol) (rand() % field->size())+1; //nonzero, not necessarily distinct
	}
	
	//generate locators (nonzero, distinct)
	int ran;
	std::vector<int> pos(field->size(),0);
	for (int i = 0; i < field->size(); i++)
	{
		pos[i]=i+1;
	}
	for (int i = 0; i < n; i++)
	{
		ran = rand() % pos.size();
		locators[i] = pos[ran]; //nonzero, distinct
		pos.erase(pos.begin() + ran); //iterator erase tricks
	}
	//make the matrix
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < m; j++)
		{
			this->setXY(i,j,field->mul(multiplicators[i],field->exp(locators[i],j)));
		}
	}
}

