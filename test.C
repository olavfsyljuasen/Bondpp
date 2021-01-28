#include<iostream>
#include<complex>
#include<vector>


using namespace std;

#include "matrixinverse.h"

int main()
{
  cout << "Start" << endl;
  const int n=5e7;
  
  vector<complex<double> > A(n*4);
  
#ifdef QMAJORORDER
  for(int i=0; i<n; i++)
    {
      const int istart=4*i;
      A[istart  ] = complex<double>(1.,0.);
      A[istart+1] = complex<double>(0.,0.);
      A[istart+2] = complex<double>(0.,0.);
      A[istart+3] = complex<double>(0.,1.);
    }
#else
  const int istart11=0;
  const int istart12=n;
  const int istart21=2*n;
  const int istart22=3*n;
  for(int i=0; i<n; i++)
    {
      A[istart11+i] = complex<double>(1.,0.);
      A[istart12+i] = complex<double>(0.,0.);
      A[istart21+i] = complex<double>(0.,0.);
      A[istart22+i] = complex<double>(0.,1.);
    }
#endif
  
  MatrixInverse(A);

  /*  
  for(int i=0; i<A.size(); i++)
    {
      cout << A[i] << endl;
    }
  */
  cout << "End" << endl;
  return 0;
}
