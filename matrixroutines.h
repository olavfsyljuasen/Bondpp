#ifndef MATRIXROUTINES_H
#define MATRIXROUTINES_H

#include<vector>
#include<complex>
#include<iostream>
#include<limits>

const bool Checkevals=true;

//#define EIGENBOOST // use both EIGEN and BOOST

#ifdef EIGENBOOST

// Different types of multiprecision backends:

#include <boost/multiprecision/cpp_complex.hpp> // complex numbers
typedef boost::multiprecision::cpp_bin_float_quad   eigen_real_type;
typedef boost::multiprecision::cpp_complex_quad eigen_complex_type;


/*
#include <boost/multiprecision/gmp.hpp> // must download gmp.h
#include <boost/multiprecision/mpc.hpp>
typedef boost::multiprecision::mpfr_float_100   eigen_real_type;
typedef boost::multiprecision::mpc_complex_100 eigen_complex_type;
*/

/* // cannot get this to work yet
#include <boost/multiprecision/float128.hpp> // real numbers
#include <boost/multiprecision/complex128.hpp> // complex numbers
typedef boost::multiprecision::float128   eigen_real_type;
typedef boost::multiprecision::complex128 eigen_complex_type;
*/

#include <boost/multiprecision/eigen.hpp>

#else
typedef realtype   eigen_real_type;
typedef complex<realtype> eigen_complex_type;
#endif



#include<Eigen/Dense>
using namespace Eigen;







void Setqzerotozero(VecMat<complex<realtype>>& A)
{
  const int Nrows=A.Nrows;
  const int Ncols=A.Ncols;

  for(int s1=0; s1<Nrows; s1++)
    for(int s2=0; s2<Ncols; s2++)
      {
	A(0,s1,s2)=0.; // assuming q=0 is first
      }
}

// A routine to sum over the ln of sublattice determinants for all q.
realtype SumLogDet(VecMat<complex<realtype>>& A)
{
  if(TRACE) cout << "Starting SumLogDet" << endl;

  realtype sum=0.;
  const int Nrows= A.Nrows;
  const int Ncols= A.Ncols;
  const int n    = A.Nvec; // the number of q values
  const int N    = A.size(); // the total array size

  
  if(Nrows==1 && Ncols==1) // det is the entry itself
    {
      for(int i=0; i<N; i++)
	{ 
	  realtype value=real(A(i,0));
	  if(value != 0.)
	    {
	      //cout << "in SumLogDet value " << i << "= " << value << endl;
	      sum+=log(value);
	    }
	}
    }
  else // NSUBL>1
    {
//#ifdef EIGENBOOST
      Matrix<eigen_complex_type,Dynamic,Dynamic> M(Nrows,Ncols);
      
      for(int k=0; k<n; k++)
	{
	  
	  for(int i=0; i<Nrows; i++)
	    for(int j=0; j<Ncols; j++)
	      M(i,j)=static_cast<eigen_complex_type>(A(k,i,j)); // row-major order
	       
	  
	  eigen_complex_type boost_det=M.determinant();
	  //if(TRACE) cout << "boost_det=" << boost_det << endl;
#ifdef EIGENBOOST
	  complex<realtype> det=boost_det.convert_to<complex<realtype> >();
#else
	  complex<realtype> det=boost_det;
#endif

	  if(TRACE) cout << setprecision(35) << det << endl;

	  realtype realdet=real(det);
	  if(realdet > 0.){ sum+=log(realdet);}
	  else
	    {
	      cout << "Warning: determinant of q=" << k
		   << " component is not positive" << endl;
	    }	  
	} 
//#endif
    }    
    return sum;
}



// Routine to subtract the minimum of an array from all its element. 
// The routine modifies the array, and returns the value subtracted. 
vector<realtype> SubtractMinimum(VecMat<complex<realtype>>& K)
{
  if(TRACE) cout << "Start SubtractMinimum" << endl;

  const int nq=K.Nvec; // the number of q components.
  const int Nrows=K.Nrows;

  vector<realtype> min(Nrows);
  
  for(int s=0; s<NSUBL; s++)
    {
      //      complex<realtype>* Kstart=K(s,s); // diagonal elements
      
      double tmin=real(K(0,s,s)); // the first element 
      
      for(int i=0; i<nq; i++){if( real(K(i,s,s))<tmin){tmin=real(K(i,s,s));}} // find minimum
      for(int i=0; i<nq; i++){K(i,s,s)-=tmin;}   // subtract it
      min[s]=tmin; // store it.
      if(TRACE) cout << "Done SubtractMinimum" << endl;
    }
  return min;
}




// Routine to subtract the minimum of an array from all its diagonal elements. 
// The routine modifies the array, and returns the value subtracted. 
// This routine also returns the element for which the minimum occurred
vector<realtype> SubtractMinimum(VecMat<complex<realtype>>& K,vector<int>& element)
{
  if(TRACE) cout << "Start SubtractMinimum" << endl;
  const int nq=K.Nvec; // the number of q components.
  const int Nrows=K.Nrows;

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
  
  if(TRACE) cout << "Done SubtractMinimum" << endl;
  return min;
}


// clean matrix, set entries less than epsilon to 0
const realtype epsilon=1e-15;

void Chomp(SMatrix<complex<realtype>>& M)
{
  for(int i=0; i<M.size(); i++)
    {
      if(abs(real(M[i]))<epsilon){ M[i].real(0.);} 
      if(abs(imag(M[i]))<epsilon){ M[i].imag(0.);} 
    }
}


void Chomp(VecMat<complex<realtype>>& K)
{
  for(int i=0; i<K.size(); i++){Chomp(K[i]);}
}



void AddDelta(VecMat<complex<realtype>>& K,const vector<realtype>& delta)
{
  for(int i=0; i<K.size(); i++)
    {
      for(int m=0; m<K.Nrows; m++)
	K[i](m,m) += delta[subl(m)];
    }
}

/*
// Sum over all q values for the matrix block indx1,indx2
complex<realtype> SumoverQ(VecMat<complex<realtype>>& A,const int indx1,const int indx2)
{
  const int nq=A.Nvec; // the number of q components.
  
  complex<realtype> sum=0.;
  for(int i=0; i<nq; i++){sum+=A(i,indx1,indx2);}
  return sum;
}
*/


void AddToDiagonal(VecMat<complex<realtype>>& A,vector<complex<realtype> > value)
{
  const int N=A.size(); // the total array size
  const int n=A.Nvec; // the number of q values
  const int Nrows=A.Nrows;

  for(int s=0; s<Nrows; s++)
    {
      for(int i=0; i<n; i++){ A(i,s,s)+=value[s];}
    }
}




void AddToDiagonal(VecMat<complex<realtype>>& A,complex<realtype> value)
{
  const int N=A.size(); // the total array size
  const int n=A.Nvec; // the number of q values
  const int Nrows=A.Nrows;
  
  for(int s=0; s<Nrows; s++)
    {
      for(int i=0; i<n; i++)
	{
	  A(i,s,s)+=value;
	}
    }
}


void SubtractFromDiagonal(VecMat<complex<realtype>>& A,complex<realtype> value)
{
  AddToDiagonal(A,-value);
}

void SubtractFromDiagonal(VecMat<complex<realtype>>& A,vector<complex<realtype> > value)
{
  for(int s=0; s<value.size(); s++){value[s]=-value[s];}
  AddToDiagonal(A,value);
}


complex<realtype> FindMax(SMatrix<complex<realtype>>& M)
{
  complex<realtype> maxentry;
  realtype maxval=0;
  for(int i=0; i<M.size(); i++)
    {
      if( abs(M[i]) > maxval ){maxval = abs(M[i]); maxentry=M[i];}
    }
  return maxentry;
}


complex<realtype> FindMax(VecMat<complex<realtype>>& A)
{
  complex<realtype> maxentry;
  realtype maxval=0;
  for(int i=0; i<A.size(); i++)
    {
      complex<realtype> thisentry=FindMax(A[i]);
      if(abs(thisentry)>maxval){maxval=abs(thisentry); maxentry=thisentry;}
    }
  return maxentry;
}


complex<realtype> FindMaxImag(SMatrix<complex<realtype>>& M)
{
  complex<realtype> maxentry;
  realtype maxval=0.;
  for(int i=0; i<M.size(); i++)
    {
      if( abs(M[i].imag()) > maxval ){maxval = abs(M[i].imag()); maxentry=M[i];}
    }
  return maxentry;
}

complex<realtype> FindMaxImag(VecMat<complex<realtype>>& A)
{
  complex<realtype> maxentry;
  realtype maxval=0;
  for(int i=0; i<A.size(); i++)
    {
      complex<realtype> thisentry=FindMaxImag(A[i]);
      if(abs(thisentry.imag()) > maxval){maxval=abs(thisentry.imag()); maxentry=thisentry;}
    }
  return maxentry;
}


void MakeReal(SMatrix<complex<realtype>>& M)
{
  for(int i=0; i<M.size(); i++){M[i].imag(0.);}
}

void MakeReal(VecMat<complex<realtype>>& A)
{
  for(int i=0; i<A.size(); i++){MakeReal(A[i]);}
}


void ComplexConjugate(SMatrix<complex<realtype>>& M)
{
  for(int i=0; i<M.size(); i++)
    {
      M[i]=conj(M[i]);
    }
}

void ComplexConjugate(VecMat<complex<realtype>>& A)
{
  for(int i=0; i<A.size(); i++){ ComplexConjugate(A[i]);}
}


void MakeRealSymmetric(SMatrix<complex<realtype>>& M)
{
  const int Nrows=M.Nrows;
  const int Ncols=M.Ncols;
  
  for(int s=0; s<Nrows; s++)
    {
      M(s,s).imag(0.); // real diagonal
    }
  
  for(int s1=0; s1<Nrows; s1++)
    for(int s2=s1+1; s2<Ncols; s2++)
      {
	    M(s1,s2).imag(0.);
	    M(s2,s1)=M(s1,s2);
      }
}

void MakeRealSymmetric(VecMat<complex<realtype>>& A)
{
  for(int i=0; i<A.size(); i++){ MakeRealSymmetric(A[i]);}
}

void MakeHermitian(SMatrix<complex<realtype>>& M)
{
  const int Nrows=M.Nrows;
  const int Ncols=M.Ncols;
  
  for(int s=0; s<Nrows; s++)
    {
      M(s,s).imag(0.);
    }
  
  for(int s1=0; s1<Nrows; s1++)
    for(int s2=s1+1; s2<Ncols; s2++)
      {
	M(s2,s1) = conj(M(s1,s2)); // off-diagonals are cc of each other.
      }
}

void MakeHermitian(VecMat<complex<realtype>>& A)
{
  for(int i=0; i<A.size(); i++){ MakeHermitian(A[i]);}
}


bool IsHermitian(SMatrix<complex<realtype>>& M)
{
  if(TRACE) cout << "Checking if Hermitian: ";
  //  bool retval=true;
  const int Nrows=M.Nrows;
  const int Ncols=M.Ncols;
  
  for(int s=0; s<Nrows; s++)
      if(imag(M(s,s)) != 0.){return false;}
  
  for(int s1=0; s1<Nrows; s1++)
    for(int s2=s1+1; s2<Ncols; s2++)
      {
	if( M(s2,s1) != conj(M(s1,s2))){ return false;}
      }
  return true;
}


bool IsHermitian(VecMat<complex<realtype>>& A)
{
 for(int i=0; i<A.size(); i++)
   { 
     if(!IsHermitian(A[i])) return false;
   }
 return true;
}



realtype FindMinimumEigenvalue(VecMat<complex<realtype>>& A)
{
  if(TRACE) cout << "Starting FindMinimumEigenvalue " << endl;

  const int n=A.size();   // the number of q values
  const int Nrows=A[0].Nrows;
  const int Ncols=A[0].Ncols;  
  
  realtype lambda_min; // the min eigenvalue of a sublattice-matrix
  realtype global_min=numeric_limits<realtype>::max();
  
  
  if(Nrows==1 && Ncols==1) 
    {
      for(int i=0; i<n; i++){if(real(A(i,0,0)) < global_min){global_min=real(A(i,0,0));}}
    }
  else // NSUBL>1
    {
      //#ifdef EIGENBOOST
      Matrix<eigen_complex_type,Dynamic,Dynamic> M(Nrows,Ncols); // use the dynamic type h
      
      for(int k=0; k<n; k++)
	{
	  for(int i=0; i<Nrows; i++)
	    for(int j=0; j<Ncols; j++)
	      M(i,j)=static_cast<eigen_complex_type>(A[k](i,j)); // row-major order
	  
	  SelfAdjointEigenSolver<Matrix<eigen_complex_type,Dynamic,Dynamic> > es(M,EigenvaluesOnly);
	  
	  eigen_real_type val= es.eigenvalues()[0]; //eigen sorts eigenvalues,least first 

#ifdef EIGENBOOST	  
	  lambda_min = val.convert_to<realtype>();
#else
	  lambda_min = val;
#endif
	  if( lambda_min < global_min){global_min=lambda_min;} 
	}
      //#endif
    }

  if(TRACE) cout << "Ending FindMinimumEigenvalue: " << global_min << endl;
  return global_min;
}

realtype SubtractMinimumEigenvalue(VecMat<complex<realtype>>& A)
{
  if(TRACE) cout << "SubtractMinimumEigenvalue" << endl;
  if(TRACE) cout << "Starting SubtractMinimumEigenvalue: " << endl;
  if(TRACE) cout << A[0] << endl;
  
  realtype emin=FindMinimumEigenvalue(A);
  if(TRACE) cout << "having found MinimumEigenvalue: " << emin << endl;
  if(TRACE) cout << A[0] << endl;
  SubtractFromDiagonal(A,emin);

  if(TRACE) cout << "mininum eigenvalue: " << setprecision(17) << emin << endl;

  if(TRACE) cout << "After subtracting from diagonal: " << endl;
  if(TRACE) cout << A[0] << endl;
  

  return emin;
}



// A routine to perform the matrix inversion, replacing the input data with output values.
void MatrixInverse(VecMat<complex<realtype>>& A)
{
  const int n=A.Nvec;   // the number of q values
  const int Nrows=A[0].Nrows;
  const int Ncols=A[0].Ncols;  

  if(Nrows==1 && Ncols==1) // just take the inverse of each element
    {
      for(int k=0; k<n; k++){ A(k,0,0)=static_cast<realtype>(1.)/A(k,0,0);}
    }
  else 
    {
      //#ifdef EIGENBOOST
      Matrix<eigen_complex_type,Dynamic,Dynamic> M(Nrows,Ncols);
      Matrix<eigen_complex_type,Dynamic,Dynamic> Minv(Nrows,Ncols);

      for(int k=0; k<n; k++)
	{

	  for(int i=0; i<Nrows; i++)
	    for(int j=0; j<Ncols; j++)
	      M(i,j)=static_cast<eigen_complex_type>(A[k](i,j)); // row-major order
	  	  
	  Minv=M.inverse();
	  
	  for(int i=0; i<Nrows; i++)
	    for(int j=i+1; j<Ncols; j++)
	      {
		Minv(i,i).imag(0); // Hermitian: set diagonal real
		Minv(j,i) = conj(Minv(i,j)) ; // Hermitian: offd are c.c.
	      }
	  
	  for(int i=0; i<Nrows; i++)
	    for(int j=0; j<Ncols; j++)
#ifdef EIGENBOOST
	      A[k](i,j) = Minv(i,j).convert_to<complex<realtype> >(); 	  
#else
	      A[k](i,j) = Minv(i,j); 	  
#endif
	} 
      //#endif     
    }
}


void MatrixPseudoInverse(VecMat<complex<realtype>>& A)
{
  const int n=A.Nvec;   // the number of q values
  const int Nrows=A.Nrows;
  const int Ncols=A.Ncols;  
  
  
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
	      M(i,j)=static_cast<eigen_complex_type>(A[k](i,j)); // row-major order
	  
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
	      A[k](i,j) = Minv(i,j);
#endif	  
	} 
      //#endif     
    }
}


// A routine to sum over the trace of sublattice matrices for all q
// it computes sum_q Tr(Aq Bq)
realtype SumTr(VecMat<complex<realtype>>& A,VecMat<complex<realtype>>& B)
{
  if(TRACE) cout << "Starting SumTr" << endl;
  realtype sum=0.;

  for(int k=0; k<A.size(); k++)
    {
      SMatrix<complex<realtype>> temp(A[k]);
      temp *= B[k];
      sum += real( tr(temp) );
    }

  if(TRACE) cout << "Ending SumTr" << endl;
  return sum;
}


/*
// A routine to sum over all q values of diagonal matrix elements s
realtype Sumq(VecMat<complex<realtype>>& A,const int s)
{
  realtype sum=0.;

  for(int i=0; i<A.size(); i++)
    {
      sum += real(A[i](s,s));
    }
  return sum;
}
*/

#endif // MATRIXROUTINES_H
