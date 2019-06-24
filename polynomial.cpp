#include "polynomial.hpp"
#include <iostream>

galois::GaloisField* Polynomial::field = 0;

Polynomial Polynomial::divBy(const Polynomial& b) const
{
	Polynomial q(0);
	Polynomial r(*this);
	if (b.getDeg() > getDeg()) return q;
	else q.setDeg(getDeg() - b.getDeg());
	for (int i = getDeg() - b.getDeg(); i >=0 ; i--)
	{
		q.setA(i,field->div(r.getA(i+b.getDeg()),b.getA(b.getDeg())));
		r = r.minus((b.constMul(q.getA(i))).shift(i));
	}
	return q;
}

Polynomial Polynomial::modBy(const Polynomial& b) const 
{
	Polynomial q(0);
	Polynomial r (*this);
	if (b.getDeg() > getDeg()) return r;
	else q.setDeg(getDeg() - b.getDeg());
	for (int i = getDeg() - b.getDeg(); i >=0 ; i--)
	{
		q.setA(i,field->div(r.getA(i+b.getDeg()),b.getA(b.getDeg())));
		r = r.minus((b.constMul(q.getA(i))).shift(i));
	}
	return r;
}

Polynomial Polynomial::derivative() const
{
	galois::GFSymbol sum;
	Polynomial result(getDeg()-1);
	for (int i = 1; i <= getDeg(); i++)
	{
		sum = (galois::GFSymbol) 0;
		for (int j = 0; j < i; j++)
		{
			sum = field->add(sum,getA(i));
		}
		result.setA(i-1,sum);
	}
	return result;
}
	
void Polynomial::EuclidForDecoding(int d, Polynomial & lambda, Polynomial & gamma) const 
{
	Polynomial a(d);
	int count = 1;
	a.setA(d,1);
	Polynomial b(*this);
	Polynomial q(0);
	std::vector<Polynomial> t(2, Polynomial(0));
	t[0].setA(0,1);
	t[1].setA(0,0);
	std::vector<Polynomial> r(2, Polynomial(0));
	r[0] = b;
	r[1] = a;
	for (int i = 1; r[(i-1) % 2].getDeg() >= d/2; i++)
	{
		q = r[i % 2].divBy(r[(i+1) % 2]);
		r[i % 2] = r[i % 2].minus(q.mul(r[(i+1) % 2]));
		t[i % 2] = t[i % 2].minus(q.mul(t[(i+1) % 2]));
		count = (count + 1) % 2;
	}
	lambda = t[(count+1) % 2];
	gamma = r[(count+1) % 2];	
}

galois::GFSymbol Polynomial::eval(galois::GFSymbol x) const
{
	galois::GFSymbol result = 0;
	int d = getDeg();
	for (int i = d; i > 0; i--)
	{
		galois::GFSymbol z = field->add(result,getA(i));
		galois::GFSymbol y = field->mul(z,x);
		result = y;
	}
	result = field->add(result, getA(0));
	return result;
}

galois::GFSymbol Polynomial::lameEval(galois::GFSymbol x) const
{
	galois::GFSymbol result = 0;
	for (int i = 0; i <= getDeg(); i++)
	{
		result = field->add(result, field->mul(getA(i),field->exp(x,i)));
	}
	return result;
}

Polynomial Polynomial::mul(const Polynomial& b) const
{
	galois::GFSymbol sum;
	Polynomial result(getDeg() + b.getDeg());
	for (int i = 0; i <= getDeg() + b.getDeg(); i++)
	{
		sum = (galois::GFSymbol)0;
		for (int j = 0; j <= i; j++)
		{
			sum = field->add(sum, field->mul(getA(j),b.getA(i-j)));
		}
		result.setA(i,sum);
	}
	return result;
}

Polynomial::Polynomial(int n)
{
	array.resize(n+1);
	for (int i = 0; i < n+1; i++)
	{
		array[i] = 0;
	}
	
}

Polynomial::Polynomial(const Polynomial& pol)
{
	array.resize(pol.getDeg()+1);
	for (int i = 0; i <= pol.getDeg(); i++)
	{
		array[i] = pol.getA(i);
	}
	
}

void Polynomial::write() const
{
	for (int i = getDeg(); i > 0; i--)
	{
		std::cout << getA(i) << "x^" << i << " + ";
	}
	std::cout << getA(0) << "x^0";
	std::cout << '\n';
}
