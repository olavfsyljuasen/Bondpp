#include<iostream>
#include<complex>

typedef double realtype ;
const int NELEMENTSTOPRINT=10;


#include "vecmat.h"

using namespace std;



main()
{
  int NMAT=3;
  int NQ=2;

  VecMat<complex<realtype> > K(NQ,NMAT,NMAT);

  K(0,1,1)=complex<realtype>(2,0);
  K(0,0,2)=complex<realtype>(1,-1);
  cout << K[0] << endl;

  VecMat<complex<realtype> > G(K);

  G *= K;

  G += K;

  cout << tr(G[0]) << endl;

  cout << G[0] << endl;
  cout << G[1] << endl;


}
