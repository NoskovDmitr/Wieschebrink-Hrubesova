#include <vector>
#include <string>
#include "stdlib.h"
#include <algorithm>
#include "matrix.hpp"
#include "regmatrix.hpp"
#include "grsmatrix.hpp"
#include "randmatrix.hpp"
#include "permmatrix.hpp"
#include "cryptosystem.hpp"

galois::GaloisField* Cryptosystem::field = 0;


Cryptosystem::Cryptosystem(GRSMatrix H, int randCount)
{
	Hgrs = H;
	G = Hgrs.findKer().transpose();
	r = randCount;
	t = (H.getNumRows()-1) / 2; //check matrix
	RegMatrix em(r + H.getNumRows());
	M = em;
	PermMatrix pe(H.getNumCols()+r);
	P = pe;
	RandMatrix ce(r,G.getNumRows());
	C = ce;
	std::vector<int> positions(G.getNumCols()+r);
	pos = positions;
	generate();
}

void Cryptosystem::regenerate()
{
	Hgrs.generate();
	G = Hgrs.findKer().transpose();
	generate();
}

void Cryptosystem::generate()
{
	M.generate();
	P.generate();	
	C.generate();
	pos.resize(G.getNumCols()+r);
	int ran = 0;
	for (int i = 0; i < G.getNumCols()+r; i++)
	{
		pos[i]=i;
	}
	
	for (int i = 0; i < G.getNumCols(); i++)
	{
		ran = rand() % pos.size();
		pos.erase(pos.begin() + ran);
	} 
	
	Matrix tempG(G.getNumCols()+r,G.getNumRows());
	int point = 0;
	for (int i = 0; i < G.getNumCols()+r; i++)
	{
		if ((point < r) && (i == pos[point]))
		{
			for (int j = 0; j < G.getNumRows(); j++)
			{
				tempG.setXY(i,j,C.getXY(point,j));
			}
			point++;
		}
		else 
		{
			for (int j = 0; j < G.getNumRows(); j++)
			{
				tempG.setXY(i,j,G.getXY(i-point,j));
			}
		}
	}
	H = tempG.findKer().transpose();
	tempG = Matrix::multiplication(M,H);
	Matrix Pub(Matrix::multiplication(tempG,P));
	this->Pub = Pub;	
}

Matrix Cryptosystem::encrypt(std::string plain)
{
	Matrix transformedPlain(stringToVectors(plain));
	return Matrix::multiplication(Pub, transformedPlain);
}

std::string Cryptosystem::decrypt(Matrix cipher)
{
	Matrix Minv = M.inverse();
	Matrix Pinv = P.inverse();
	std::vector<int> punc(G.getNumCols()- G.getNumRows(),0);
	for (int i = 0; i < punc.size(); i++)
	{
		punc[i] = i;
	}
	Matrix Gx = H.findKer().transpose().puncture(pos).puncture(punc).transpose();
	Matrix Gprime = H.findKer().transpose();
	
	Matrix stepOne = Matrix::multiplication(Minv,cipher);
	Matrix result(cipher.getNumCols(),Pub.getNumCols());
	Matrix rightSide(1,stepOne.getNumRows());
	
	for (int j = 0; j < stepOne.getNumCols(); j++)
	{
		for (int i = 0; i < stepOne.getNumRows(); i++)
		{
			rightSide.setXY(0,i,stepOne.getXY(j,i));
		}
		Matrix eprime(H.solve(rightSide));
		Matrix east(eprime.transpose().puncture(pos).transpose());
		Matrix err(Hgrs.decode(east));
		Matrix s(1,err.getNumRows());
		for (int i = 0; i < s.getNumRows(); i++)
		{
			s.setXY(0,i,field->sub(east.getXY(0,i),err.getXY(0,i)));
		}
		Matrix sPunct(s.transpose().puncture(punc).transpose());
		Matrix a((Gx.solve(sPunct)).transpose());
		Matrix fin(Matrix::multiplication(a,Gprime));
		Matrix realFin(1,eprime.getNumRows());
		for (int i = 0; i < eprime.getNumRows(); i++)
		{
			realFin.setXY(0,i,field->sub(eprime.getXY(0,i),fin.getXY(i,0)));
		}
		realFin = Matrix::multiplication(Pinv, realFin);
		for (int i = 0; i < eprime.getNumRows(); i++)
		{
			result.setXY(j,i,realFin.getXY(0,i));
		}
		
	}
	return vectorsToString(result);
}

std::string Cryptosystem::vectorsToString(Matrix vectors)
{
	int m = vectors.getNumCols();
	int n = vectors.getNumRows();
	InsertBuffer buffer;
	for (int j = 0; j < m; j++)
	{
		for (int i = 0; i < n; i++)
		{
			if (vectors.getXY(j,i)!=(galois::GFSymbol)0)
			{
				for (int k = 1; k < field->pwr(); k++)
				{
					buffer.setNextBit(((int)vectors.getXY(j,i) >> (field->pwr()-1-k)) & 1);
				}
			}
		}
	}
	
	return buffer.getData();
}

Matrix Cryptosystem::stringToVectors(std::string text)
{
	int m = (text.length()*8)/((field->pwr()-1)*t)+1;
	int n = Hgrs.getNumCols()+r;
	std::vector<int> positions(n,0);
	int x = 0;
	int ran = 0;
	Buffer buffer(text);

	Matrix result(m,n);
	for (int j = 0; j < m; j++)
	{
		positions.resize(n);
		for (int i = 0; i < n; i++)
		{
			positions[i] = i;
		}
		for (int i = 0; i < n-t; i++)
		{
			ran = rand() % positions.size();
			positions.erase(positions.begin() + ran); 
		}
		for (int i = 0; i < t; i++)
		{
			if (buffer.isEmpty()) break; //done
			x = 1;
			for (int k = 0; k < field->pwr()-1; k++)
			{
				if (! buffer.isEmpty()) x = (x << 1) ^ buffer.getNextBit();
				else x = (x << 1)^ 0 ;
			}
			result.setXY(j,positions[i],(galois::GFSymbol)x);
		}
	}
	return result;
}

Buffer::Buffer(std::string data)
{
	this->data = data;
	counter = 0;
}

bool Buffer::isEmpty()
{
	return (counter >= data.length()*8);
}

int Buffer::getNextBit()
{
	int subCount = counter % 8;
	int result = (((int)data[counter / 8]) >> (7-subCount)) & 1;
	counter++;
	return result;
}

InsertBuffer::InsertBuffer()
{
	counter = 0;
	minBuf = (char) 0;
}

void InsertBuffer::setNextBit(int bit)
{
	minBuf = ((int)minBuf << 1) ^ (bit & 1);
	counter++;
	if ((counter % 8) == 0) {
		data.push_back(minBuf);
		minBuf = 0;
	}
}

std::string InsertBuffer::getData()
{
	int a = counter % 8;
	data.push_back((char)((int)minBuf << (8-a)));
	return data;
}
