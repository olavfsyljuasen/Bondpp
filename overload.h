#ifndef OVERLOAD_H
#define OVERLOAD_H

#include<iostream>
#include<complex>

using namespace std;

#ifdef FLOAT
complex<float> operator+(const int& l,const complex<float>& r){return complex<float>(l+real(r),imag(r));}
complex<float> operator*(const int& l,const complex<float>& r){return complex<float>(l*real(r),l*imag(r));}
#elif defined LONGDOUBLE
complex<long double> operator+(const int& l,   const complex<long double>& r){return complex<long double>(l+real(r),imag(r));}
complex<long double> operator*(const int& l,   const complex<long double>& r){return complex<long double>(l*real(r),l*imag(r));}
#else
complex<double> operator+(const int& l,  const complex<double>& r){return complex<double>(l+real(r),imag(r));}
complex<double> operator*(const int& l,  const complex<double>& r){return complex<double>(l*real(r),l*imag(r));}
#endif



#endif //OVERLOAD_H
