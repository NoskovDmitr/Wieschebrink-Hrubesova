#include "randmatrix.hpp"
#include "regmatrix.hpp"
#include "stdlib.h"

class Matrix;

void RegMatrix::generate()
{
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < m; j++)
		{
			this->setXY(i,j,(galois::GFSymbol)0);
		}
	}
	int ran=1;
	int rank=0;
	for (int i = 0; i < n; i++)
	{
		while (rank != (i+1)) {
			for (int j = 0; j < m; j++)
			{
				ran = rand() % field->size();
				setXY(i,j,(galois::GFSymbol)ran);
			}
			rank = compRank();
		}
	}
	
}

RegMatrix::RegMatrix(int dim)
{
	n = dim;
	m = dim;
	
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

RegMatrix::RegMatrix()
{
}

RegMatrix::RegMatrix(const RegMatrix& mat)
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

RegMatrix RegMatrix::inverse()
{
	Matrix comp(2*n,n);
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			comp.setXY(i,j,getXY(i,j));
		}
		comp.setXY(n+i,i,(galois::GFSymbol)1);
	}
	Matrix comp2 = comp.GaussElimination();
	for (int i = n-1; i > 0; i--)
	{
		for (int j = i-1; j >= 0; j--)
		{
			comp2.plusKTimesRow(field->sub((galois::GFSymbol)0,comp2.getXY(i,j)),i,j);
		}
	}
	RegMatrix result(n);
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
				result.setXY(i,j,comp2.getXY(i+n,j));
		}
	}
	return result;
}
