#ifndef MATRIXROUTINES_H
#define MATRIXROUTINES_H

#include<vector>
#include<complex>
#include<iostream>
#include<limits>
#include<cassert>
#include "lowlevelmatrixroutines.h"

const bool Checkevals=true;




template<class T,int Nrows,int Ncols>
void Setqzerotozero(VecMat<T,Nrows,Ncols>& A)
{
  for(int s1=0; s1<Nrows; s1++)
    for(int s2=0; s2<Ncols; s2++)
      {
	A(0,s1,s2)=0.; // assuming q=0 is first
      }
}

// A routine to sum over the ln of sublattice determinants for all q.
template<class T,int Nrows,int Ncols>
realtype SumLogDet(VecMat<T,Nrows,Ncols>& A)
{
  if(TRACELEVEL>0) cout << spaces(ir++) << "Starting SumLogDet" << endl;

  realtype sum=0.;
  const int n    = A.Nvecs(); // the number of q values

  vector<complextype> a(Nrows*Ncols);
  
  for(int k=0; k<n; k++)
    {
      for(int i=0; i<Nrows; i++)
	for(int j=0; j<Ncols; j++)
	  a[i+j*Nrows]=A(k,i,j);
      	    
      realtype value=SmallHermitianMatrixDeterminant(Nrows,a); // col-major order
      if(TRACELEVEL>4 ) cout << spaces(ir) << "k=" << k << " logdet:" << setprecision(DEBUGPRECISION) << value << " " << mylog(value) << endl;
      if(value > 0.){ sum += mylog(value);} 
    }
  if(TRACELEVEL>0) cout << spaces(--ir) << "Finished SumLogDet()" << endl;
  return sum;
}



// Routine to subtract the minimum of an array from all its element. 
// The routine modifies the array, and returns the value subtracted. 
template<class T,int Nrows,int Ncols>
  vector<realtype> SubtractMinimum(VecMat<T,Nrows,Ncols>& K)
{
  if(TRACELEVEL>0) cout << spaces(ir++) << "Start SubtractMinimum" << endl;

  const int nq=K.Nvecs(); // the number of q components.

  vector<realtype> min(Nrows);
  
  for(int s=0; s<Nrows; s++)
    {
      //      complex<realtype>* Kstart=K(s,s); // diagonal elements
      
      double tmin=real(K(0,s,s)); // the first element 
      
      for(int i=0; i<nq; i++){if( real(K(i,s,s))<tmin){tmin=real(K(i,s,s));}} // find minimum
      for(int i=0; i<nq; i++){K(i,s,s)-=tmin;}   // subtract it
      min[s]=tmin; // store it.
    }

  if(TRACELEVEL>0) cout << spaces(--ir) << "Done SubtractMinimum" << endl;

  return min;
}




// Routine to subtract the minimum of an array from all its diagonal elements. 
// The routine modifies the array, and returns the value subtracted. 
// This routine also returns the element for which the minimum occurred
template<class T,int Nrows,int Ncols>
  vector<realtype> SubtractMinimum(VecMat<T,Nrows,Ncols>& K,vector<int>& element)
{
  if(TRACELEVEL>0) cout << spaces(ir++) << "Start SubtractMinimum" << endl;
  const int nq=K.Nvecs(); // the number of q components.

  vector<realtype> min(Nrows);
  
  for(int s=0; s<Nrows; s++)
    {
      //      complex<realtype>* Kstart=K(s,s); // diagonal
      
      realtype tmin=real(K(0,s,s)); // the first element 
      
      int telement=0;
      
      for(int i=0; i<nq; i++){if(real(K(i,s,s))<tmin){tmin=real(K(i,s,s)); telement=i;}} // find minimum
      for(int i=0; i<nq; i++){K(i,s,s)-=tmin;}   // subtract it
      
      min[s]=tmin;
      element[s]=telement; // store it
    }
  
  if(TRACELEVEL>0) cout <<  spaces(--ir) << "Done SubtractMinimum" << endl;
  return min;
}


// clean matrix, set entries less than epsilon to 0



template<class T,int Nrows,int Ncols>
  void AddDelta(VecMat<T,Nrows,Ncols>& K,NumberList& delta)
{
  for(int q=0; q<K.Nvecs(); q++)
    {
      for(int m=0; m<Nrows; m++)
	K(q,m,m) += delta[subl(m)];
    }
}



template<class T,int Nrows,int Ncols>
  realtype FindMinimumEigenvalue(VecMat<T,Nrows,Ncols>& A)
{
  if(TRACELEVEL>0) cout << spaces(ir++) << "Starting FindMinimumEigenvalue()" << endl;

  const int n=A.Nvecs();   // the number of q values

  realtype global_min=numeric_limits<realtype>::max();

  vector<complextype> a(Nrows*Ncols);
  
  for(int k=0; k<n; k++)
    {
      for(int i=0; i<Nrows; i++)
	for(int j=0; j<Ncols; j++)
	  a[i+j*Nrows]=A(k,i,j);

      
      realtype eval=SmallHermitianMatrixMinEigenvalue(Nrows,a); // col-major order
      if(eval < global_min){global_min=eval;}
    }
  if(TRACELEVEL>2) cout << spaces(ir) << "Minimum eval=" << global_min << endl;
  if(TRACELEVEL>0) cout << spaces(--ir) << "Finished FindMinimumEigenvalue()" << endl;
  return global_min;
}

template<class T,int Nrows,int Ncols>
realtype SubtractMinimumEigenvalue(VecMat<T,Nrows,Ncols>& A)
{
  realtype emin=FindMinimumEigenvalue(A);
  if(TRACELEVEL>0) cout << spaces(ir) << "having found MinimumEigenvalue: " << emin << endl;
  
  SubtractFromDiagonal(A,complex<realtype>(emin,0.));
  
  return emin;
}



// A routine to perform the matrix inversion, replacing the input data with output values.
template<class T,int Nrows,int Ncols>
  void MatrixInverse(VecMat<T,Nrows,Ncols>& A)
{
  if(TRACELEVEL>0) cout << spaces(ir++) << "Starting MatrixInverse" << endl;
  const int n=A.Nvecs();   // the number of q values

  assert(Nrows == Ncols);
  
  vector<complextype> a(Nrows*Ncols);
  
  for(int k=0; k<n; k++)
    {
      for(int i=0; i<Nrows; i++)
	for(int j=0; j<Ncols; j++)
	  a[i+j*Nrows]=A(k,i,j);
      	    
      SmallMatrixInverse(Nrows,a); // col-major order

      for(int i=0; i<Nrows; i++)
	for(int j=0; j<Ncols; j++)
	  A(k,i,j) = a[i+j*Nrows];
    }

  if(TRACELEVEL>0) cout << spaces(--ir) << "Finished  MatrixInverse" << endl;
}

/*
template<class T,int Nrows,int Ncols>
  void MatrixPseudoInverse(VecMat<T,Nrows,Ncols>& A)
{
  const int n=A.Nvecs();   // the number of q values
  
  if(Nrows==1 && Ncols==1) // just take the inverse of each element
    {
      for(int k=0; k<n; k++){ A(k,0,0)=static_cast<realtype>(1.)/A(k,0,0);}
    }
  else 
    {
      //#ifdef EIGENBOOST
      const realtype tolerance=Nrows*std::numeric_limits<realtype>::epsilon();
      
      
      Matrix<eigen_complex_type,Dynamic,Dynamic> M(Nrows,Ncols);
      Matrix<eigen_complex_type,Dynamic,Dynamic> Minv(Nrows,Ncols);
      Matrix<eigen_complex_type,Dynamic,Dynamic> singularValuesInv(Nrows,Ncols);
      
      for(int k=0; k<n; k++)
	{
	  
	  for(int i=0; i<Nrows; i++)
	    for(int j=0; j<Ncols; j++)
	      M(i,j)=static_cast<eigen_complex_type>(A(k,i,j)); // row-major order
	  
	  auto svd = M.jacobiSvd(Eigen::ComputeFullU | Eigen::ComputeFullV);
	  const auto &singularValues = svd.singularValues();
	  
	  singularValuesInv.setZero();
	  
	  for (unsigned int i = 0; i < singularValues.size(); ++i) {
	    if (singularValues(i) > tolerance)
	      {
		singularValuesInv(i, i) = 1. / singularValues(i);
	      }
	    else
	      {
		singularValuesInv(i, i) = 0;
	      }
	  }
	  Minv= svd.matrixV() * singularValuesInv * svd.matrixU().adjoint();
	  
	  
	  for(int i=0; i<Nrows; i++)
	    for(int j=i+1; j<Ncols; j++)
	      {
		Minv(i,i).imag(0); // Hermitian: set diagonal real
		Minv(j,i) = conj(Minv(i,j)) ; // Hermitian: offd are c.c.
	      }

#ifdef EIGENBOOST	  
	  for(int i=0; i<Nrows; i++)
	    for(int j=0; j<Ncols; j++)
	      A[k](i,j) = Minv(i,j).convert_to<complex<realtype> >(); 
#else
	  for(int i=0; i<Nrows; i++)
	    for(int j=0; j<Ncols; j++)
	      A(k,i,j) = Minv(i,j);
#endif	  
	} 
      //#endif     
    }
}
*/

// A routine to sum over the trace of sublattice matrices for all q
// it computes sum_q Tr(Aq Bq)
template<class T,int Nrows,int Ncols>
  realtype SumTr(VecMat<T,Nrows,Ncols>& A,VecMat<T,Nrows,Ncols>& B)
{
  if(TRACELEVEL>0) cout << spaces(ir++) << "Starting SumTr" << endl;
  realtype sum=0.;

  for(int k=0; k<A.Nvecs(); k++)
    {
      SMatrix<complex<realtype>,Nrows,Ncols> temp(A[k]);
      temp *= B[k];
      sum += real( tr(temp) );
    }

  if(TRACELEVEL>0) cout <<  spaces(--ir) << "Ending SumTr" << endl;
  return sum;
}



#endif // MATRIXROUTINES_H
