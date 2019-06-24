#include "regmatrix.hpp"
#include "permmatrix.hpp"
#include "stdlib.h"
#include <vector>

class Matrix;

void PermMatrix::generate()
{
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < m; j++)
		{
			this->setXY(i,j,(galois::GFSymbol)0);
		}
	}
	int ran;
	std::vector<int> pos(n,0);
	for (int i = 0; i < n; i++)
	{
		pos[i]=i;
	}
	//std::cout << n << '\n';
	for (int i = 0; i < n; i++)
	{
		ran = rand() % pos.size();
		setXY(i,pos[ran],(galois::GFSymbol)1);
		pos.erase(pos.begin() + ran); 
	}
}

PermMatrix::PermMatrix(int dim)
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

PermMatrix::PermMatrix()
{

}

PermMatrix::PermMatrix(const PermMatrix& mat)
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

int PermMatrix::ChangePos(int pos)
{
	for (int i = 0; i < n; i++)
	{
		if (getXY(i,pos) == (galois::GFSymbol)1) return i;
	}
	return 0;
}
