#ifndef OVERLOAD_H
#define OVERLOAD_H

#include<iostream>
#include<complex>

using namespace std;

complextype operator+(const int& l, const complextype& r){ return complextype(l+real(r),imag(r));}
complextype operator*(const int& l, const complextype& r){ return complextype(l*real(r),l*imag(r));}

#endif //OVERLOAD_H
