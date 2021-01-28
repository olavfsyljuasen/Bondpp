#include<iostream>
#include<complex>

#include <fftw3.h>

using namespace std;



extern "C"
{
  //#g++ -I/scratch/sylju/include -L/scratch/sylju/lib FFTWtest.C -lfftw3 -lm;

}

const bool TRACE = false;

const int NSUBL=1;
const int NSUBL2=NSUBL*NSUBL;

#include "globalheader.h"

#include "qsmatrix_new.h"

#include "matrixroutines.h"

#include "rnddef.h"

int main()
{
int dim=3;
vector<int> dims(dim);

for(int L=1; L<102; L++)
  {

cout << "L= " << L;
dims[0]=L;
dims[1]=1;
dims[2]=1;

 int Vq=dims[0]*dims[1]*dims[2];
 
 int Vtot=Vq*NSUBL2;

 QSmatrix A(Vtot);
 QSmatrix B(Vtot);
 
 fftw_complex* m=A.start();
 
 fftw_plan Aq_to_Ar = fftw_plan_many_dft(dim,&dims[0],NSUBL2,m,0,1,Vq,m,0,1,Vq,1,FFTW_MEASURE);
 
 fftw_plan Ar_to_Aq = fftw_plan_many_dft(dim,&dims[0],NSUBL2,m,0,1,Vq,m,0,1,Vq,-1,FFTW_MEASURE);
 
 for(int i=0; i<Vtot; i++)
   {
     double r=RAN();
     A[i]=complex<double>(r,0.);
   }
 
 //cout << "Aq= " << A << endl;
 
 fftw_execute(Aq_to_Ar);
 
 //cout << "Ar= " << A << endl;
 
 //for(int i=0; i<Vq; i++) A[i]*=(1./Vq);
 
 fftw_execute(Ar_to_Aq);
 
 //cout << "Aq= " << A << endl;
 cout << "Aq max imag part " << FindMaxImag(A) << endl; 

 cout << A(1,2)[0] << endl;
 
 fftw_destroy_plan(Aq_to_Ar);
 fftw_destroy_plan(Ar_to_Aq);
 
  }
}

