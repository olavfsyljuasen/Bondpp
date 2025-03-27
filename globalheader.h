#ifndef GLOBALHEADER_H
#define GLOBALHEADER_H

#include<fstream>
#include<iostream>
#include<iomanip>
#include<stdio.h>
#include<stdlib.h>
#include<string>
#include<vector>
#include<time.h>
#include"math.h"

using namespace std;





string bitstring(int s,int cycle)
{
  //  string str=(s==0 ? "0":"");
  string str="";
  //  for(int i=0; (s!=0 || i<cycle); s>>=1, i++){if(s & 1){str='1'+str;}else{str='0'+str;}}
  for(int i=0; i<cycle; s>>=1, i++){if(s & 1){str='1'+str;}else{str='0'+str;}}
  return str;
}


class Timer{
  friend ostream&
    operator<<(ostream& os, Timer timer)
    { time_t t;
      if( time(&t) == time_t(-1)){exit(1);}
      os << ctime(&t); return os;
    }  
 public:
  Timer(){};
  ~Timer(){if(TRACE) cout << "deleting timer\n";}
  void Start(){t1=time(0);}
  void Stop(){t2=time(0);}
  double GetTimeElapsed(){return difftime(t2,t1);}
 private:
  time_t t1;
  time_t t2;
};


class Coord
{  friend ostream& operator<<(ostream& os,Coord& c)
    { os << c.x << " " << c.y << " " << c.z; return os; }
 public:
  Coord(realtype xx=0.,realtype yy=0.,realtype zz=0.):x(xx),y(yy),z(zz){}
  realtype x;
  realtype y;
  realtype z;
  realtype operator[](int i) const {return (i==0 ? x: (i==1 ? y: z));}
  void clear(){x=0.;y=0.;z=0.;}
  realtype Norm(){return sqrt(x*x+y*y+z*z);}
  realtype Norm2(){return x*x+y*y+z*z;}
  void Normalize(){realtype g=mysqrt(x*x+y*y+z*z); x /= g; y /= g; z /= g;}
  int size() const {return 3;}
};

inline bool operator==(const Coord& a,const Coord& b)
{ return bool( (a.x == b.x) && (a.y==b.y) && (a.z==b.z));}

Coord operator-(const Coord& l,const Coord& r){
  return Coord(l.x-r.x,l.y-r.y,l.z-r.z);
}

Coord operator+(const Coord& l,const Coord& r){
  return Coord(l.x+r.x,l.y+r.y,l.z+r.z);
}

realtype operator*(const Coord& l,const Coord& r){
  return realtype(l.x*r.x+l.y*r.y+l.z*r.z);
}

Coord operator*(const int& k,const Coord& a){
  return Coord(k*a.x,k*a.y,k*a.z);
}

Coord operator*(const Coord& a,const int& k){
  return Coord(k*a.x,k*a.y,k*a.z);
}

Coord operator*(const realtype& k,const Coord& a){
  return Coord(k*a.x,k*a.y,k*a.z);
}

Coord operator*(const Coord& a,const realtype& k){
  return Coord(k*a.x,k*a.y,k*a.z);
}

Coord& operator+=(Coord& a,const Coord& b){
  a.x += b.x;   a.y += b.y;   a.z += b.z; return a;
}

Coord& operator-=(Coord& a,const Coord& b){
  a.x -= b.x;   a.y -= b.y;   a.z -= b.z; return a;
}

realtype scalarproduct(const Coord& a,const Coord& b){
return a.x*b.x+a.y*b.y+a.z*b.z;
}

Coord crossproduct(const Coord& a,const Coord& b){
  Coord result(a.y*b.z-a.z*b.y,a.z*b.x-a.x*b.z,a.x*b.y-a.y*b.x); return result;
}




// a routine that converts integers to strings
const int MAXDIGITS = 20;
const int ASCII0 = 48;
string int2string(int a)
{
  int digit[MAXDIGITS];
  int d=0;
  string s("");
  if (a == 0){ s+=ASCII0; return s;}
  while( a > 0){ digit[d++] = a-10*(a/10); a /=10 ;}
  for(int i=d-1; i>=0; i--) s+= digit[i]+ASCII0;
  return s;
}


// this routine also handles leading zeros
string int2string(int a,int nsiffer)
{
  int digit[MAXDIGITS];
  string s("");
  for(int d=0; d<MAXDIGITS; d++){ digit[d] = a-10*(a/10); a /=10 ;}
  for(int i=nsiffer-1; i>=0; i--) s+= digit[i]+ASCII0;
  return s;
}

const int NELEMENTSTOPRINT=4;

// Routine for printing complex matrices
ostream& operator<<(ostream& os,vector<complex<realtype> >& M)
{
  bool printdots=false;
  int end=M.size();
  if(M.size() > NELEMENTSTOPRINT){end=NELEMENTSTOPRINT; printdots=true;}
  os << setprecision(20);
  for(int i=0; i<end; i++) 
    os << setprecision(20) << M[i] << " ";
  if(printdots) os << " ... " << M[M.size()-1];
  return os;
}

ostream& operator<<(ostream& os,vector<realtype >& M)
{
  bool printdots=false;
  int end=M.size();
  if(M.size() > NELEMENTSTOPRINT){end=NELEMENTSTOPRINT; printdots=true;}
  os << setprecision(20);
  for(int i=0; i<end; i++) 
    os << setprecision(20) << M[i] << " ";
  if(printdots) os << " ... " << M[M.size()-1];
  return os;
}



struct NumberList
{
  NumberList(int n=0,realtype a=0.):v(n,a){}
  //  NumberList(NumberList r):v(r.size(),0){for(unsigned int i=0; i<r.size(); i++){v[i]=r[i];}}
  friend ostream& operator<<(ostream& os,NumberList& d){for(unsigned int i=0; i<d.size(); i++){ os << d[i] << " ";} return os;}
  unsigned int size(){return v.size();}
  realtype& operator[](const int i){return v[i];}

  //  NumberList& operator=(NumberList& r) { for(unsigned int i=0; i<v.size(); i++){v[i]=r[i];} return *this;}
  NumberList& operator=(NumberList r) { for(unsigned int i=0; i<v.size(); i++){v[i]=r[i];} return *this;}
  NumberList& operator+=(NumberList& r){ for(unsigned int i=0; i<v.size(); i++){v[i]+=r[i];} return *this;}
  NumberList& operator*=(realtype& k)  { for(unsigned int i=0; i<v.size(); i++){v[i]*=k;} return *this;}
private:
  vector<realtype> v;
}; 

bool operator==(NumberList& l,NumberList& r)
{
  bool isequal=true;
  for(unsigned int i=0; i<l.size(); i++){isequal &= (l[i]==r[i]); if(!isequal){break;}}
  return isequal;
}

NumberList operator-(NumberList& l,NumberList& r)
{
  NumberList a=l;
  for(unsigned int i=0; i<a.size(); i++){ a[i] -= r[i];}
  return a;
}

NumberList operator+(NumberList& l,NumberList& r)
{
  NumberList a=l;
  for(unsigned int i=0; i<a.size(); i++){ a[i] += r[i];}
  return a;
}

NumberList operator*(realtype& k,NumberList& r)
{
  NumberList a=r;
  for(unsigned int i=0; i<a.size(); i++){ a[i] *= k;}
  return a;
}

NumberList operator*(NumberList& l, realtype& k)
{
  NumberList a=l;
  for(unsigned int i=0; i<a.size(); i++){ a[i] *= k;}
  return a;
}


realtype maxabs(NumberList& a)
{
  realtype maxval=abs(a[0]);
  if(a.size() >1)
    {
      for(unsigned int i=1; i<a.size(); i++)
	{
	  realtype thisabs = abs(a[i]);
	  if( thisabs > maxval){maxval = thisabs;}
	}
    }
  return maxval;
}

#endif
