#include "matrix.hpp"
#include "GaloisField.hpp"
#include <iostream>
#include <exception>
 
galois::GaloisField* Matrix::field = 0;

class noSolutionException: public std::exception
{
	virtual const char* what() const throw()
	{
		return "No solution for equation system.";
	}
} noSolutionEx;

bool Matrix::equals(Matrix b)
{
	if ((getNumCols() != b.getNumCols()) || (getNumRows() != b.getNumRows())) return false;
	else 
	{
		for (int i = 0; i < getNumCols(); i++)
		{
			for (int j = 0; j < getNumRows(); j++)
			{
				if (getXY(i,j) != b.getXY(i,j)) return false;
			}		
		}
	}
	return true;
}

int Matrix::compRank()
{
	Matrix temp;
	if (n < m) temp = transpose().GaussElimination();
	else temp = GaussElimination();
	int result = 0;
	bool emptyLine;
	for (int i = 0; i < m; i++)
	{
		emptyLine = true;
		for (int j = 0; j < n; j++)
		{
			if (temp.getXY(j,i) != (galois::GFSymbol)0)
			{
				emptyLine = false;
				result++;
				break;
			}
		}
		if (emptyLine) break;
	}
	return result;
}


int Matrix::squareStarCodeDim()
{
	Matrix comp((m*(m-1))/2+m,n);
	int counter = 0;
	for (int i = 0; i < m; i++)
	{
		for (int j = 0; j <= i; j++)
		{
			for (int k = 0; k < n; k++)
			{
				comp.setXY(counter,k,field->mul(getXY(k,i),getXY(k,j)));
			}	
			counter++;
		}
	}
	return comp.compRank();
}

Matrix Matrix::transpose()
{
	Matrix result(m,n);
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < m; j++)
		{
				result.setXY(j,i,getXY(i,j));
		}
	}
	return result;
}

Matrix Matrix::findKer()
{
	Matrix comp = GaussElimination();
	int pos = 0;
	std::vector<int> basis(n,0);
	int bcount = 0;
	std::vector<int> fre(n,0);
	int fcount = 0;
	for (int i = 0; i < n; i++)
	{
		if (comp.getXY(i,pos) != (galois::GFSymbol)0)
		{
			if (m<=i) 
			{
				fre[fcount] = i; 
				fcount++;
			} else {
				basis[bcount] = i;
				bcount++;
				pos++;
			}
		} 
		else 
		{
			fre[fcount] = i;
			fcount++;
		}
	}
	Matrix result(fcount,n);
	for (int i = 0; i < fcount; i++)
	{
		result.setXY(i,fre[i], (galois::GFSymbol)1);
	}
	galois::GFSymbol sum;
	for (int i = 0; i < fcount; i++)
	{
		for (int j = bcount-1; j >= 0 ; j--)
		{
			if (fre[i] > basis[j])
			{
				sum = 0;
				for (int k = basis[j]; k < n; k++)
				{
					sum = field->add(sum,field->mul(result.getXY(i,k),comp.getXY(k,j)));
				}
				result.setXY(i,basis[j],field->sub( (galois::GFSymbol)0, field->div(sum, comp.getXY(basis[j],j))));
			}
		}
	}	
	return result;
}

Matrix Matrix::puncture(std::vector<int>& positions)
{
	int vec = 0;
	Matrix result(n-positions.size(),m);
	for (int i = 0; i < n; i++)
	{
		if ((vec < (int)positions.size()) && (positions[vec] == i))
		{
			vec++;
			continue;
		}
		for (int j = 0; j < m; j++)
		{
			result.setXY(i-vec,j,getXY(i,j));
		}
	}
	return result;
}

Matrix Matrix::punctureInv(std::vector<int>& positions)
{
	Matrix result(positions.size(),m);
	for (int i = 0; i < positions.size(); i++)
	{
		for (int j = 0; j < m; j++)
		{
			result.setXY(i,j,getXY(positions[i],j));
		}
	}
	return result;
}

Matrix::Matrix(int col, int row)
{
	n = col;
	m = row;
	
	array.clear();
	array.resize(m*n,0);
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < m; j++)
		{
			this->setXY(i,j,(galois::GFSymbol)0);
		}
	}
}

Matrix::Matrix()
{
	n = 0;
	m = 0;
}

Matrix::Matrix(const Matrix& mat)
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

void Matrix::plusKTimesRow(galois::GFSymbol k, int rowFrom, int rowTo)
{
	for (int i = 0; i < n; i++)
	{		
		setXY(i,rowTo,field->add(getXY(i, rowTo),field->mul(k, getXY(i, rowFrom))));
	}
}
	
Matrix Matrix::GaussElimination()
{
	Matrix result(getNumCols(),getNumRows());
	for (int i = 0; i < getNumCols(); i++)
	{
		for (int j = 0; j < getNumRows(); j++)
		{
			result.setXY(i,j,getXY(i,j));
		}
		
	}
	
	bool ok = false;
	int toFirst;
	galois::GFSymbol swap,pivot;
	
	int pos = 0; //leftmost nonzero position
	int posMem = 0;
	for (int i = 0; i < m; i++)
	{
		ok = false;
		posMem = pos;
		bool cont = false;
		while (! ok) {
			if (pos == n) {pos = posMem; cont = true; break;}
			for (int j = i; j < m; j++)
			{
				if (result.getXY(pos,j) != (galois::GFSymbol)0) {
					toFirst = j;
					pivot = result.getXY(pos,j);
					ok = true;
					break;
				}
			}
			if (! ok ) pos++;
		}
		if (cont) continue;
		for (int j = pos; j < n; j++)
		{
			swap = result.getXY(j,toFirst);
			result.setXY(j,toFirst,result.getXY(j, i));
			result.setXY(j, i, field->mul(swap,field->inverse(pivot)));
		}
		
		galois::GFSymbol minus;
		if (i!=m-1) for (int j = i+1; j < m; j++)
		{
			if (result.getXY(pos,j) != (galois::GFSymbol)0) 
			{
				minus = field->sub((galois::GFSymbol)0,result.getXY(pos,j));		
				result.plusKTimesRow(minus,i,j);
			}
		}
		pos++;	
		if (pos == n) break;	
	}
	return result;
}

Matrix Matrix::solve(Matrix& rightSide)
{
	Matrix comp(n+1,m);
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < m; j++)
		{
			comp.setXY(i,j,getXY(i,j));
		}
	}
	if (rightSide.getNumCols() == 0) throw noSolutionEx; //this is plain wrong.
	for (int i = 0; i < m; i++)
	{
		comp.setXY(n,i,rightSide.getXY(0,i));
	}
	
	Matrix temp(comp.findKer());
	if (temp.getNumCols() == 0) /*st. terrible happened*/ throw noSolutionEx;
	for (int i = 0; i < temp.getNumCols()-1; i++)
	{
		for (int j = 0; j < temp.getNumRows(); j++)
		{
			temp.setXY(temp.getNumCols()-1,j,field->add(temp.getXY(temp.getNumCols()-1,j),temp.getXY(i,j)));
		}
	}
	
	Matrix result(1,n);
	for (int i = 0; i < n; i++)
	{
		result.setXY(0,i,temp.getXY(temp.getNumCols()-1,i));
	}
	
	return result;
}

Matrix Matrix::multiplication(Matrix& A, Matrix& B)
{
	galois::GFSymbol sum = (galois::GFSymbol)0;
	if (A.getNumCols() != B.getNumRows()) throw;
	Matrix result(B.getNumCols(),A.getNumRows());
	for (int i = 0; i < A.getNumRows(); i++)
	{
		for (int j = 0; j < B.getNumCols(); j++)
		{
			for (int k = 0; k < A.getNumCols(); k++)
			{
				sum = field->add(sum, (field->mul(A.getXY(k,i),B.getXY(j,k))));
			}
			result.setXY(j,i,sum);
			sum = (galois::GFSymbol)0;
		}
	}
	result = result;
	return result;
}

void Matrix::WriteMatrix()
{
		for (int i = 0; i < getNumRows(); i++)
		{
			for (int j = 0; j < getNumCols(); j++)
			{
				if ((int) getXY(j,i) >= 100) std::cout << getXY(j,i) << ' ';
				else if ((int) getXY(j,i) >= 10) std::cout << ' ' << getXY(j,i) << ' ';
				else std::cout << "  " << getXY(j,i) << ' ';
			}
			std::cout << '\n';
		}
		std::cout << "----------------------\n";
		
}
