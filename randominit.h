#ifndef RANDOMINIT_H
#define RANDOMINIT_H

//routine to make a random self-energy that respects inversion symmetry.
//It works by making random real real-space entries. This is Fourier-transformed,
//and the real part is kept. This ensures invariance under q->-q.

#include "rnddef.h"
  
#include<vector>
#include<complex>
#include<iostream>

extern "C"
{
  //#g++ -I/data/sylju/include -L/data/sylju/lib FFTWtest.C -lfftw3 -lm;
#include <fftw3.h>
}


void MakeRandomMatrixWithInversionSymmetry(QSmatrix& Kinit,int dim, vector<int>& dims,double da)
{
  if(TRACE) cout << "Starting MakeRandomMatrixWithInversionSymmetry" << endl;

  int Vq=1;
  int Vf=1;
  for(int i=0; i<(dim-1); i++){Vq*=dims[i]; Vf*=dims[i];}
  Vq*= dims[dim-1];
  Vf*= dims[dim-1]/2+1;


  vector<complex<double> > RandomRealSpaceEntries(Vf);
  //  vector<double> Kinit(Vq);

  fftw_plan p1 = fftw_plan_dft_c2r(dim,&dims[0],reinterpret_cast<fftw_complex*>(&RandomRealSpaceEntries[0]),&Kinit[0],FFTW_ESTIMATE);

  //Assign random elements
  const double factor=2.*da/Vq;
  for(int i=0; i<Vf; i++){RandomRealSpaceEntries[i]=complex<double>(factor*(RAN()-0.5),0.);}

  if(TRACE)
    { 
      cout << "RandomRealSpaceEntries" << endl;
      //      FourierPrint(Vq,RandomRealSpaceEntries);
    }


  fftw_execute(p1); // RandomRealSpaceEntries -> Kinit

  if(TRACE)
    { 
      cout << "Kinit" << endl;
      cout << Kinit << endl;
    }

  if(TRACE) cout << "Done with MakeRandomMatrixWithInversionSymmetry " << endl;
}


void MakeRandomMatrix(vector<double>& Kinit,int dim, vector<int>& dims,double da)
{
  if(TRACE) cout << "Starting MakeRandomMatrix" << endl;

  int Vq=1;
  int Vf=1;
  for(int i=0; i<(dim-1); i++){Vq*=dims[i]; Vf*=dims[i];}
  Vq*= dims[dim-1];
  Vf*= dims[dim-1]/2+1;

  // forget about inversion symmetry and just do a random thing
  for(int i=0; i<Vq; i++){Kinit[i]=da*(RAN()-0.5);}

  if(TRACE)
    { 
      cout << "Kinit" << endl;
      cout << Kinit << endl;
    }

  if(TRACE) cout << "Done with MakeRandomMatrix " << endl;
}


#endif
