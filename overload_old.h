#ifndef OVERLOAD_H
#define OVERLOAD_H

#include<iostream>
#include<complex>

using namespace std;

#ifdef FLOAT
complex<float> operator+(const int& l,const complex<float>& r){return complex<float>(l+real(r),imag(r));}
complex<float> operator*(const int& l,const complex<float>& r){return complex<float>(l*real(r),l*imag(r));}
complex<float> expi(const float phi){return polar(1.,phi);}

#elif defined LONGDOUBLE
complex<long double> operator+(const int& l,   const complex<long double>& r){return complex<long double>(l+real(r),imag(r));}
complex<long double> operator*(const int& l,   const complex<long double>& r){return complex<long double>(l*real(r),l*imag(r));}
complex<long double> operator+(const float& l, const complex<long double>& r){return complex<long double>(l+real(r),imag(r));}
complex<long double> operator*(const float& l, const complex<long double>& r){return complex<long double>(l*real(r),l*imag(r));}
complex<long double> operator+(const double& l,const complex<long double>& r){return complex<long double>(l+real(r),imag(r));}
complex<long double> operator*(const double& l,const complex<long double>& r){return complex<long double>(l*real(r),l*imag(r));}

complex<long double> expi(const long double phi){return polar(1.,phi);}
#else
complex<double> operator+(const int& l,  const complex<double>& r){return complex<double>(l+real(r),imag(r));}
complex<double> operator*(const int& l,  const complex<double>& r){return complex<double>(l*real(r),l*imag(r));}
complex<double> operator+(const float& l,const complex<double>& r){return complex<double>(l+real(r),imag(r));}
complex<double> operator*(const float& l,const complex<double>& r){return complex<double>(l*real(r),l*imag(r));}

complex<double> expi(const double phi){return polar(1.,phi);}
//complex<double> expi(const double phi){return complex<double>(cos(phi),sin(phi));}
#endif



#endif //OVERLOAD_H
