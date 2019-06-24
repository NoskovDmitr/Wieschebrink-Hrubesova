#include "matrix.hpp"
#include "randmatrix.hpp"
#include "stdlib.h"

void RandMatrix::generate()
{
	int ran;
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < m; j++)
		{
			ran = rand() % field->size();
			setXY(i,j,(galois::GFSymbol)ran);
		}
		
	}
	
}

RandMatrix::RandMatrix(int col, int row)
{
	n = col;
	m = row;
	
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

RandMatrix::RandMatrix()
{
}

RandMatrix::RandMatrix(const Matrix& mat)
{
	this->n = mat.getNumCols();
	this->m = mat.getNumRows();
	array.clear();
	array.resize(this->m*this->n);
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < m; j++)
		{
			this->setXY(i,j,mat.getXY(i,j));
		}
	}
}
