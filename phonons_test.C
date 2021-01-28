#include<iomanip>
#include<vector>
//#include<complex>
#include<iostream>
#include<string>
#include <fstream>

using namespace std;

typedef double realtype;


#ifndef M_PI
#define M_PI 3.14159265358979323846264338328
#define M_PIl 3.141592653589793238462643383279502884L 
#endif

#ifndef SQRTTHREEOVERTWO
#ifdef LONGDOUBLE
#define SQRTTHREEOVERTWO 0.8660254037844386467637231707529362L
#else
#define SQRTTHREEOVERTWO 0.866025403784438646763723170753
#endif
#endif


// universal definitions
#ifdef LONGDOUBLE
const realtype PI=M_PIl;
#else
const realtype PI=M_PI;
#endif
const realtype TWOPI=2.*PI;


#define QSPACEDIRECTLATTICE
#define SIMPLECUBICBRAVAISLATTICE

#include "overload.h"

const bool TRACE=false;

const int MAXRUNS=10;
const int MAXPAR=200;
const string READIN="read.in";

const string SITERPTS="siteRpts.dat";
const string SITEQPTS="siteQpts.dat";
const bool PRINTRPTS=true;
const bool PRINTQPTS=true;

const int MAXNSELECTEDQPTS=20;
const string SELECTEDQPTSFILENAME="selectedqpts.in";

const int NPARAMS = 8;
enum params{LINEID,ALPHA1,ALPHA2,NX,NY,NZ,NBINS,EQFLAG};

// compile using g++ -std=c++11 -I/scratch/sylju/eigen-3.3.7 -I/scratch/sylju/boost_1_73_0 phonons_test.C; ./a.out


#include "globalheader.h"

ofstream logfile("log.txt",ios::app);


#include "RunParameter.h"

#include "bravaislattices.h"

#include "phonons.h"



int main()
{
  RunParameters rp;
  logfile << "parameters: " << rp;

  double* par = rp.GetPars(1);

  cout << "nx= " << par[NX] << endl;
  BravaisLattice la(par);

  cout << la.SiterVol() << " " << la.SiteqVol() << endl;

  Phonons ph(par,la);

  



}
