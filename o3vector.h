#ifndef O3VECTOR_H
#define O3VECTOR_H

//#include<fstream>
#include<iostream>
//#include<stdio.h>
//#include<stdlib.h>
#include"math.h"

using namespace std;

inline int         conj(int a){return a;}
inline float       conj(float a){return a;}
inline double      conj(double a){return a;}
inline long double conj(long double a){return a;}


template<class T>
class O3vector
{  friend ostream& operator<<(ostream& os,O3vector& c)
    { os << c.x << " " << c.y << " " << c.z; return os; }
 public:
 O3vector(T xx=0,T yy=0,T zz=0):x(xx),y(yy),z(zz){}
 O3vector(Coord c):x(c.x),y(c.y),z(c.z){}
  T x;
  T y;
  T z;
  void clear(){x=0.;y=0.;z=0.;}
  realtype Norm(){return sqrt(conj(x)*x+conj(y)*y+conj(z)*z);}
};

template<class T>
inline bool operator==(const O3vector<T>& a,const O3vector<T>& b)
{ return bool( (a.x == b.x) && (a.y==b.y) && (a.z==b.z));}

template<class T>
O3vector<T> operator-(const O3vector<T>& l,const O3vector<T>& r){
  return O3vector<T>(l.x-r.x,l.y-r.y,l.z-r.z);
}

template<class T>
O3vector<T> operator+(const O3vector<T>& l,const O3vector<T>& r){
  return O3vector<T>(l.x+r.x,l.y+r.y,l.z+r.z);
}

template<class T>
O3vector<T> operator*(const O3vector<T>& l,const O3vector<T>& r){
  return O3vector<T>(l.x*r.x+l.y*r.y+l.z*r.z);
}

template<class T>
O3vector<T> operator*(const int& k,const O3vector<T>& a){
  return O3vector<T>(k*a.x,k*a.y,k*a.z);
}

template<class T>
O3vector<T> operator*(const O3vector<T>& a,const int& k){
  return O3vector<T>(k*a.x,k*a.y,k*a.z);
}

template<class T>
O3vector<T> operator*(const realtype& k,const O3vector<T>& a){
  return O3vector<T>(k*a.x,k*a.y,k*a.z);
}

template<class T>
O3vector<T> operator*(const O3vector<T>& a,const realtype& k){
  return O3vector<T>(k*a.x,k*a.y,k*a.z);
}

template<class T>
O3vector<T> operator*(const complex<realtype>& k,const O3vector<T>& a){
  return O3vector<T>(k*a.x,k*a.y,k*a.z);
}

template<class T>
O3vector<T> operator*(const O3vector<T>& a,const complex<realtype>& k){
  return O3vector<T>(k*a.x,k*a.y,k*a.z);
}

template<class T>
O3vector<T> operator*(const Coord& k,const O3vector<T>& a){
  return O3vector<T>(k.x*a.x,k.y*a.y,k.z*a.z);
}

template<class T>
O3vector<T> operator*(const O3vector<T>& a,const Coord& k){
  return O3vector<T>(k.x*a.x,k.y*a.y,k.z*a.z);
}

template<class T>
O3vector<T>& operator+=(O3vector<T>& a,const O3vector<T>& b){
  a.x += b.x;   a.y += b.y;   a.z += b.z; return a;
}

template<class T>
O3vector<T>& operator-=(O3vector<T>& a,const O3vector<T>& b){
  a.x -= b.x;   a.y -= b.y;   a.z -= b.z; return a;
}

template<class T>
T scalarproduct(const O3vector<T>& a,const O3vector<T>& b){
  return conj(a.x)*b.x+conj(a.y)*b.y+conj(a.z)*b.z;
}

template<class T>
O3vector<T> crossproduct(const O3vector<T>& a,const O3vector<T>& b){
  O3vector<T> result(conj(a.y)*b.z-conj(a.z)*b.y,conj(a.z)*b.x-conj(a.x)*b.z,conj(a.x)*b.y-conj(a.y)*b.x); 
  return result;
}




#endif
