#include <vector>
#include <algorithm>
#include <exception>
#include "cryptosystem.hpp"
#include <string>
#include "GaloisField.hpp"
#include "matrix.hpp"
#include "regmatrix.hpp"
#include "grsmatrix.hpp"
#include "attacker.hpp"

galois::GaloisField* Attacker::field = 0;

class starSquareDimException: public std::exception
{
	virtual const char* what() const throw()
	{
		return "Wrong square code dimension.";
	}
} starSquareDimEx;

class duplicateCodeLocatorsException: public starSquareDimException
{
	virtual const char* what() const throw()
	{
		return "Duplicate code locators."; 
	}
} duplicateCodeLocatorsEx;

class tooManyRandomVectorsException: public starSquareDimException
{
	virtual const char* what() const throw()
	{
		return "Too many random vectors."; 
	}
} tooManyRandomVectorsEx;

void Attacker::randomPositions1()
{
	std::vector<int> punc(1,0);
	Matrix Gen = CS.getPubMatrix().findKer().transpose();
	int d = Gen.squareStarCodeDim();
	
	for (int i = 0; i < Gen.getNumCols(); i++)
	{
		punc[0] = i;
		int dim = Gen.puncture(punc).squareStarCodeDim();
		if (dim == d-1) 
		{
			rPos.push_back(i);
		}
		else if (dim == d) ; //return false
		else throw starSquareDimEx; //return null
	}	
}

void Attacker::randomPositions2()
{
	std::vector<int> punc(1,0);
	int n = CS.getPubMatrix().getNumCols();
	int a = CS.getPubMatrix().getNumCols() - CS.getPubMatrix().getNumRows() - CS.getPubMatrix().getNumRows()/4;
	int trig;
	std::vector<int> A(a,0);
	
	for (int i = 0; i < n; i++)
	{
		trig = 0;
		for (int j = 0; j < a; j++)
		{
			if (j==i) trig++;
			A[j] = j + trig;
		}
		Matrix shortened = CS.getPubMatrix().puncture(A).findKer().transpose(); 
		int d = shortened.squareStarCodeDim();
		if (i <= a) punc[0] = 0;
		else punc[0] = i-a;
		Matrix temp = shortened.puncture(punc);
		int dim = temp.squareStarCodeDim();
		if (dim == d-1) 
		{
			rPos.push_back(i);
		}
		else if (dim == d) ; //return false
		else throw starSquareDimEx; //return null		
	}	
}

std::string Attacker::decrypt(Matrix cipher)
{
	Matrix Minv = M.inverse();
	std::vector<int> punc(Gen.getNumCols()- Gen.getNumRows(),0);
	for (int i = 0; i < punc.size(); i++)
	{
		punc[i] = i;
	}
	Matrix Gx = Hprime.findKer().transpose().puncture(rPos).puncture(punc).transpose();
	Matrix Gprime = Hprime.findKer().transpose();
	
	Matrix stepOne = Matrix::multiplication(Minv,cipher);
	Matrix result(cipher.getNumCols(),CS.getPubMatrix().getNumCols());
	Matrix rightSide(1,stepOne.getNumRows());
	
	for (int j = 0; j < stepOne.getNumCols(); j++)
	{
		for (int i = 0; i < stepOne.getNumRows(); i++)
		{
			rightSide.setXY(0,i,stepOne.getXY(j,i));
		}
		Matrix eprime(Hprime.solve(rightSide));
		Matrix east(eprime.transpose().puncture(rPos).transpose());
		Matrix err(Hgrs.decode(east));
		Matrix s(1,err.getNumRows());
		for (int i = 0; i < s.getNumRows(); i++)
		{
			s.setXY(0,i,field->sub(east.getXY(0,i),err.getXY(0,i)));
		}
		Matrix sPunct(s.transpose().puncture(punc).transpose());
		Matrix a((Gx.solve(sPunct)).transpose());
		Matrix fin(Matrix::multiplication(a,Gprime));
		for (int i = 0; i < eprime.getNumRows(); i++)
		{
			result.setXY(j,i,field->sub(eprime.getXY(0,i),fin.getXY(i,0)));
		}
	}
	return cutString(CS.vectorsToString(result));
}

void Attacker::matrixM()
{
	Matrix Hpub = CS.getPubMatrix();
	RegMatrix result(Hpub.getNumRows());
	int rankHprime = Hprime.compRank();
	std::vector<int> punc(Hprime.getNumCols()-rankHprime,0);
	for (int i = 0; i < punc.size(); i++)
	{
		punc[i] = i;
	}
	
	Matrix tempPub(Hpub.puncture(punc).transpose());
	Matrix tempPrime(Hprime.puncture(punc).transpose());
	Matrix sol;
	Matrix rightSide(1,Hpub.getNumRows());
	for (int j = 0; j < Hpub.getNumRows(); j++)
	{
		for (int i = 0; i < Hpub.getNumRows(); i++)
		{
			rightSide.setXY(0,i,tempPub.getXY(j,i));
		}
		sol = tempPrime.solve(rightSide);
		for (int i = 0; i < Hpub.getNumRows(); i++)
		{
			result.setXY(i,j,sol.getXY(0,i));
		}
	}
	M = result;
}

void Attacker::GRScode(Matrix Gen2)
{
	//last SegFault: too many random vectors
	if (CS.getPubMatrix().getNumCols() - rPos.size() < CS.getPubMatrix().getNumRows()) throw tooManyRandomVectorsEx;
	
	std::vector<galois::GFSymbol> muls(1,0);
	locators = SidelnikovSestakov1(Gen2,false);
	
	RegMatrix M2(1);
	muls = SidelnikovSestakov2(Gen2,M2,false);
	RegMatrix Minv = M2.inverse();
	Matrix nakedG (Matrix::multiplication(Minv,G));
	GRSMatrix toH(Gen2.getNumCols()-1,Gen2.getNumCols(),locators,muls);
	Matrix temp(toH.findKer());
	std::vector<galois::GFSymbol> mulsToH(Gen2.getNumCols(),0);
	for (int i = 0; i < Gen2.getNumCols(); i++)
	{
		mulsToH[i] = temp.getXY(0,i); 
	}
	GRSMatrix H(Gen2.getNumCols()-Gen2.getNumRows(),Gen2.getNumCols(),locators,mulsToH);
	Hgrs = H;
	Hprime = (nakedG.findKer().transpose());
}
	
bool Attacker::ripntear()
{
	rPos.clear();
	bool iDidIt = false;
	std::string test = "Testin, one, two...";
	try 
	{
		int t =CS.getT();
		int k =CS.getPubMatrix().getNumCols() - CS.getPubMatrix().getNumRows();
		int nr = CS.getPubMatrix().getNumCols();
		int n1 = 2*t+k;
		int r1 = nr - n1;
		int n2 = 2*t+k+1;
		int r2 = nr - n2;
		if (nr < 4*t+1)
		{
			randomPositions1();
		}
		else if (nr >= 4*t+3)
		{
			randomPositions2();		
		}
		else
		{
			randomPositions1();
			if ((rPos.size() != r1) && (rPos.size() != r2)) 
			{
				rPos.clear();
				randomPositions2();
				if ((rPos.size() != r1) && (rPos.size() != r2)) 
					throw starSquareDimEx;
			}
		}
		
		G = CS.getPubMatrix().findKer().transpose();
		Gen = G.puncture(rPos);
		GRScode(Gen);
		matrixM();
			
		Matrix enc = CS.encrypt(test);			
		if (test == (decrypt(enc))) 
		{	
			iDidIt = true;
		}
	}
	catch (std::exception& e)
	{
		std::cout << "Exception: ";
		std::cout << e.what() << '\n';
		iDidIt = false;
	}
	return iDidIt;
}


std::vector<std::string> displayResult()
{
	std::vector<std::string> result;
	return result;
}

bool Attacker::iHaveDataForAllAlpha(std::vector<bool> checklist)
{
	for (int i = 0; i < checklist.size(); i++)
	{
		if (checklist[i] == false) return false;
	}
	return true;
}

std::vector<galois::GFSymbol> Attacker::SidelnikovSestakov1(Matrix & T, bool parityCheck)
{
	int n = T.getNumCols();
	int k = T.getNumRows();
	std::vector<int> alpha(n,0);
	alpha[0] = 1;
	alpha[1] = 0;
	alpha[2] = -1;
	std::vector<bool> iHaveIt(n,false);
	iHaveIt[0] = true;
	iHaveIt[1] = true;
	iHaveIt[2] = true;
	std::vector<int> punc(1,0);
	Matrix u1(k,1);
	Matrix u2(k,1);
	
	bool gotOne;
	int goodCounter;
	int badCounter;
	while (! iHaveDataForAllAlpha(iHaveIt))
	{
		punc.clear();
		punc.push_back(0);
		goodCounter = 3;
		badCounter = 3;
		for (int i = 0; i < k-2; i++)
		{
			gotOne = false;
			while (! gotOne)
			{
				if ((iHaveIt[goodCounter] == true))
				{
					gotOne = true;
					punc.push_back(goodCounter);
					goodCounter++;
				}
				else if (goodCounter < n)
				{
					goodCounter++;
				}
				else break;
			}
			if (! gotOne)
			{
				while (! gotOne)
				{
					if ((iHaveIt[badCounter] == false))
					{
						gotOne = true;
						punc.push_back(badCounter);
						badCounter++;
					}
					else 
					{
						badCounter++;
					}
				}
			}	
		}
		std::sort(punc.begin(), punc.end());
		Matrix temp = T.punctureInv(punc).transpose();
		Matrix ker = temp.findKer();
		for (int i = 0; i < k; i++)
		{
			u1.setXY(i,0,ker.getXY(0,i));
		}
		punc[0] = 1;
		
		temp = T.punctureInv(punc).transpose();
		ker = temp.findKer();
		for (int i = 0; i < k; i++)
		{
			u2.setXY(i,0,ker.getXY(0,i));
		}
		Matrix c1 = Matrix::multiplication(u1, T);
		Matrix c2 = Matrix::multiplication(u2, T);
		
		galois::GFSymbol t3 = field->div(c1.getXY(2,0),c2.getXY(2,0));
		for (int i = 3; i < n; i++)
		{
			if ((iHaveIt[i] == false) && (c1.getXY(i,0) != 0) && (c2.getXY(i,0) != 0))
			{
				alpha[i] = field->div(t3, field->sub(t3, field->div(c1.getXY(i,0),c2.getXY(i,0))));
				iHaveIt[i] = true;
			}
		}
	}
	
	//lastSegFault: nondistinct code locators
	for (int i = 0; i < n-1; i++)
	{
		for (int j = i+1; j < n; j++)
		{
			if (alpha[i] == alpha[j]) throw duplicateCodeLocatorsEx;
		}
	}
	
	//transformation
	galois::GFSymbol a = 0;
	int f = (int)field->size();
	std::vector<bool> whatIsTheNumber(f,true);
	for (int i = 0; i < n; i++)
	{
		if (alpha[i] != -1) whatIsTheNumber[alpha[i]] = false;
	}
	for (int i = 0; i < n; i++)
	{
		if (whatIsTheNumber[i] == true) { 
			a = (galois::GFSymbol) i ;
			break;
		}
	}
	alpha[2] = 0;
	for (int i = 0; i < n; i++)
	{
		if (i != 2) alpha[i] = (int)field->div((galois::GFSymbol)1,field->sub(a,(galois::GFSymbol)alpha[i]));
	}
	
	//transformation2
	a = 0;
	for (int i = 0; i < f; i++)
	{
		whatIsTheNumber[i] = true;
	}
	for (int i = 0; i < n; i++)
	{
		if (alpha[i] != -1) whatIsTheNumber[alpha[i]] = false;
	}
	for (int i = 0; i < n; i++)
	{
		if (whatIsTheNumber[i] == true) { 
			a = (galois::GFSymbol) i ;
			break;
		}
	}
	for (int i = 0; i < n; i++)
	{
		alpha[i] = (int)field->sub((galois::GFSymbol)alpha[i],a);
	}	
	
	std::vector<galois::GFSymbol> result(n,0);
	for (int i = 0; i < n; i++)
	{
		result[i] = (galois::GFSymbol)alpha[i];
		
	}
	locators = result;
	
	return result;
}

std::vector<galois::GFSymbol> Attacker::SidelnikovSestakov2(Matrix & T, RegMatrix & M, bool parityCheck) 
{
	int n = T.getNumCols();
	int k;
	k = T.getNumCols() - T.getNumRows(); //algorithm works with $(n-k) \times n$ matrix
	
	
	std::vector<galois::GFSymbol> result(n,0); //v_i
	
	
	//step one: find a_i
	std::vector<int> toPuncture(k-1,0);
	for (int i = 0; i < k-1; i++)
	{
		toPuncture[i]= n-k+1+i;
	}
	Matrix comp(T.puncture(toPuncture));
	Matrix ker(comp.findKer());
	Matrix c(n-k+1,1);
	for (int i = 0; i < n-k+1; i++)
	{
		c.setXY(i,0,ker.getXY(0,i));
	}
	Matrix trans(c.transpose());
	
	//step two: find v_i (up to n-k)
	result[0] = (galois::GFSymbol)1;
	
	Matrix right(1,n-k);
	galois::GFSymbol a = field->sub((galois::GFSymbol)0,c.getXY(0,0));
	galois::GFSymbol alpha = locators[0];
	Matrix comp2(n-k,n-k);
	for (int i = 0; i < n-k; i++)
	{
		for (int j = 0; j < n-k; j++)
		{
			comp2.setXY(j,i,field->mul(c.getXY(j+1,0),field->exp(locators[j+1],i)));
		}
	}
	
	for (int i = 0; i < n-k; i++)
	{		
		right.setXY(0,i,a);
		a = field->mul(a,alpha);
	}
	Matrix sol(comp2.solve(right));
	for (int i = 0; i < n-k; i++)
	{
		result[i+1] = sol.getXY(0,i);
	}
	//OK
	
	//step three: find M
	Matrix moreRight(n-k,n-k);
	for (int i = 0; i < n-k; i++)
	{
		for (int j = 0; j < n-k; j++)
		{
			moreRight.setXY(i,j,field->mul(field->inverse(result[j]),T.getXY(j,i)));
		}
	}
	
	Matrix bigOne(2*(n-k),n-k);
	for (int i = 0; i < n-k; i++)
	{
		for (int j = 0; j < n-k; j++)
		{
			bigOne.setXY(j,i, field->exp(locators[i],j));
		}
		for (int j = 0; j < n-k; j++)
		{
			bigOne.setXY(j+n-k,i, moreRight.getXY(j,i));
		}
	}
	Matrix biggerOne(bigOne.GaussElimination());
	toPuncture.resize(n-k);
	for (int i = 0; i < n-k; i++)
	{
		toPuncture[i] = n-k+i;
	}
	
	Matrix leftSide(biggerOne.puncture(toPuncture));
	RegMatrix M2(n-k);
	M = M2;
	for (int i = 0; i < n-k; i++)
	{
		for (int j = 0; j < n-k; j++)
		{
			right.setXY(0,j,biggerOne.getXY(i+n-k,j));
		}
		sol = leftSide.solve(right);
		for (int j = 0; j < n-k; j++)
		{
			M.setXY(j,i,sol.getXY(0,j));
		}
	}
	
	//step four, finish
	Matrix bait(M.inverse());
	Matrix lastOne(Matrix::multiplication(bait,T));
	for (int i = n-k; i < n; i++)
	{
		result[i] = lastOne.getXY(i,0);
	}
	return result;
}

