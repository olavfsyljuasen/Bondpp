#ifndef OBSERVABLES_H
#define OBSERVABLES_H

#include <functional>
using namespace std;

typedef std::function<complex<realtype>(int)> qfunc;

// a function for use in calculating observables 
class KernelFunction
{
 public:
  KernelFunction(){}
  virtual complex<realtype> operator()(int qi,int m1,int m2) =0;
};


class FixedIndxObs: public KernelFunction
{
 public:
 FixedIndxObs(const int in_m1,const int in_m2,qfunc in_f):mym1(in_m1),mym2(in_m2),myfunc(in_f){} 
 FixedIndxObs(const int in_s1,const int in_l1,const int in_s2, const int in_l2,qfunc in_f):mym1(mindx(in_s1,in_l1)),mym2(mindx(in_s2,in_l2)),myfunc(in_f){}
  complex<realtype> operator()(int qi,int m1,int m2)
    {
      if(m1==mym1 && m2==mym2){return myfunc(qi);}
      else{ return 0.;}
    }
 private:
  const int mym1;
  const int mym2;
  qfunc myfunc;
};


// q-functions

static const Triplet Tx={1,0,0};
static const Triplet Ty={0,1,0};

complex<realtype> cosaxes(int qi)
{
  return 0.5*(cos( la.qr(qi,Tx))-cos( la.qr(qi,Ty)));
}

static const Triplet Td= {1, 1,0};
static const Triplet Tmd={1,-1,0};

complex<realtype> cosdiag(int qi)
{
  return 0.5*(cos( la.qr(qi,Td))-cos( la.qr(qi,Tmd)));
}

static const Triplet TO={0,0,0};
static const Triplet Ta1={1, 0,0};
static const Triplet Ta2={0, 1,0};
static const Triplet Ta3={-1,1,0};
static const complex<realtype> w(-0.5,SQRTTHREEOVERTWO);

complex<realtype> honey01(int qi)
{
  return 0.5*(1.+w*expi(-la.qr(qi,Ta3))+w*w*expi(-la.qr(qi,Ta2)));
}

complex<realtype> honey00(int qi)
{
  return 0.5*(cos(la.qr(qi,Ta1))+w*cos(la.qr(qi,Ta2))+w*w*cos(la.qr(qi,Ta3)));
}



typedef complex<realtype> obstype; 

//OBSERVABLES:

#if defined SQUARELATTICE
const int NOBSERVABLES=2;
//enum observables{SA,SD};
string NAMESOFOBSERVABLES[NOBSERVABLES]={"sa","sd"};


const int NSPINOBSERVABLES=2;
enum spinobservables{SA,SD};
string NAMESOFSPINOBSERVABLES[NSPINOBSERVABLES]={"sa","sd"};


FixedIndxObs saobs(0,0,0,0,cosaxes);
FixedIndxObs sdobs(0,0,0,0,cosdiag);

KernelFunction* spinobservables[NSPINOBSERVABLES]={&saobs,&sdobs};


#elif defined HEXAGONALLATTICE
const int NOBSERVABLES=2;
//enum observables{SA,SD};
string NAMESOFOBSERVABLES[NOBSERVABLES]={"s00","s01"};


const int NSPINOBSERVABLES=2;
enum spinobservables{S00,S01};
string NAMESOFSPINOBSERVABLES[NSPINOBSERVABLES]={"s00","s01"};

FixedIndxObs s00obs(0,0,0,0,honey00);
FixedIndxObs s01obs(0,0,0,1,honey01);

KernelFunction* spinobservables[NSPINOBSERVABLES]={&s00obs,&s01obs};



#elif defined CUBICLATTICE
const int NOBSERVABLES=2;
//enum observables{SA,SD};
string NAMESOFOBSERVABLES[NOBSERVABLES]={"sa","sd"};


const int NSPINOBSERVABLES=2;
enum spinobservables{SA,SD};
string NAMESOFSPINOBSERVABLES[NSPINOBSERVABLES]={"sa","sd"};


FixedIndxObs saobs(0,0,0,0,cosaxes);
FixedIndxObs sdobs(0,0,0,0,cosdiag);

KernelFunction* spinobservables[NSPINOBSERVABLES]={&saobs,&sdobs};

//const int NOBSERVABLES=3;
//enum observables{A1M,A2M,A2P};
//string NAMESOFOBSERVABLES[NOBSERVABLES]={"a1m","a2m","a2p"};

#elif defined FCCLATTICE
const int NOBSERVABLES=2;
//enum observables{SA,SD};
string NAMESOFOBSERVABLES[NOBSERVABLES]={"sa","sd"};


const int NSPINOBSERVABLES=2;
enum spinobservables{SA,SD};
string NAMESOFSPINOBSERVABLES[NSPINOBSERVABLES]={"sa","sd"};


FixedIndxObs saobs(0,0,0,0,cosaxes);
FixedIndxObs sdobs(0,0,0,0,cosdiag);

KernelFunction* spinobservables[NSPINOBSERVABLES]={&saobs,&sdobs};

//const int NOBSERVABLES=3;
//enum observables{A1M,A2M,A2P};
//string NAMESOFOBSERVABLES[NOBSERVABLES]={"a1m","a2m","a2p"};

#elif defined BCCLATTICE
const int NOBSERVABLES=2;
//enum observables{SA,SD};
string NAMESOFOBSERVABLES[NOBSERVABLES]={"sa","sd"};


const int NSPINOBSERVABLES=2;
enum spinobservables{SA,SD};
string NAMESOFSPINOBSERVABLES[NSPINOBSERVABLES]={"sa","sd"};


FixedIndxObs saobs(0,0,0,0,cosaxes);
FixedIndxObs sdobs(0,0,0,0,cosdiag);

KernelFunction* spinobservables[NSPINOBSERVABLES]={&saobs,&sdobs};

//const int NOBSERVABLES=3;
//enum observables{A1M,A2M,A2P};
//string NAMESOFOBSERVABLES[NOBSERVABLES]={"a1m","a2m","a2p"};

#elif defined HEXAGONALBRAVAISLATTICE
const int NOBSERVABLES=5;
enum observables{ROT,PHASEII,PHASEIII,SIGMATA,SIGMATD};
string NAMESOFOBSERVABLES[NOBSERVABLES]={"rot","II","III","SA","SD"};
#else
const int NOBSERVABLES=4;
enum observables{SIGMAA,SIGMAD,SIGMAMX,SIGMAMY};
string NAMESOFOBSERVABLES[NOBSERVABLES]={"sa","sd","smx","smy"};
string NAMESOFALPHAS[NOBSERVABLES]={"sa_alpha.dat","sd_alpha.dat","smx_alpha.dat","smy_alpha.dat"};
#endif































#endif // OBSERVABLES_H
