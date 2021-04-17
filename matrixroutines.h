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
  const int n    = A.Nvecs; // the number of q values
  const int N    = A.size(); // the total array size

  
  if(Nrows==1 && Ncols==1) // det is the entry itself
    {
      for(int i=0; i<N; i++)
	{ 
	  realtype value=real(A(i,0));
	  if(value != 0.)
	    {
	      if(TRACE) cout << "in SumLogDet value " << i << "= " << value << endl;
	      sum+=log(value);
	    }
	}
    }
  else 
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

  const int nq=K.Nvecs; // the number of q components.
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
  const int nq=K.Nvecs; // the number of q components.
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




void AddDelta(VecMat<complex<realtype>>& K,NumberList& delta)
{
  for(int q=0; q<K.Nvecs; q++)
    {
      for(int m=0; m<K.Nrows; m++)
	K(q,m,m) += delta[subl(m)];
    }
}




realtype FindMinimumEigenvalue(VecMat<complex<realtype>>& A)
{
  if(TRACE) cout << "Starting FindMinimumEigenvalue " << endl;

  const int n=A.Nvecs;   // the number of q values
  const int Nrows=A.Nrows;
  const int Ncols=A.Ncols;  
  
  realtype lambda_min; // the min eigenvalue of a sublattice-matrix
  realtype global_min=numeric_limits<realtype>::max();
  
  
  if(Nrows==1 && Ncols==1) 
    {
      for(int i=0; i<n; i++){if(real(A(i,0,0)) < global_min){global_min=real(A(i,0,0));}}
    }
  else
    {
      //#ifdef EIGENBOOST
      Matrix<eigen_complex_type,Dynamic,Dynamic> M(Nrows,Ncols); // use the dynamic type h
      
      for(int k=0; k<n; k++)
	{
	  for(int i=0; i<Nrows; i++)
	    for(int j=0; j<Ncols; j++)
	      M(i,j)=static_cast<eigen_complex_type>(A(k,i,j)); // row-major order
	  
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
  realtype emin=FindMinimumEigenvalue(A);
  if(TRACE) cout << "having found MinimumEigenvalue: " << emin << endl;

  SubtractFromDiagonal(A,complex<realtype>(emin,0.));

  return emin;
}



// A routine to perform the matrix inversion, replacing the input data with output values.
void MatrixInverse(VecMat<complex<realtype>>& A)
{
  if(TRACE) cout << "Starting MatrixInverse" << endl;
  const int n=A.Nvecs;   // the number of q values
  const int Nrows=A.Nrows;
  const int Ncols=A.Ncols;  

  if(TRACE) cout << "Nrows= " << Nrows << " Ncols="<< Ncols << " n=" << n << endl;

  if(Nrows==1 && Ncols==1) // just take the inverse of each element
    {
      for(int k=0; k<n; k++){ A(k,0,0)=complex<realtype>(1./A(k,0,0).real(),0.);}
    }
  else 
    {
      // TEST to avoid matrix diag.
      /*
      for(int k=0; k<n; k++)
	for(int s=0; s<NSPIN; s++)
	  {
	    A(k,s,s)=complex<realtype>(1./A(k,s,s).real(),0.);
	  }
      */

      Matrix<eigen_complex_type,Dynamic,Dynamic> M(Nrows,Ncols);
      Matrix<eigen_complex_type,Dynamic,Dynamic> Minv(Nrows,Ncols);

      if(TRACE) cout << "Initialized M and Minv " << endl;


      for(int k=0; k<n; k++)
	{
	  for(int i=0; i<Nrows; i++)
	    for(int j=0; j<Ncols; j++)
	      {
		M(i,j)=static_cast<eigen_complex_type>(A(k,i,j)); // row-major order
	      }

	  Minv=M.inverse();

	  
	  for(int i=0; i<Nrows; i++)
	    for(int j=0; j<Ncols; j++)
	      {
#ifdef EIGENBOOST
		A(k,i,j) = Minv(i,j).convert_to<complex<realtype> >(); 	  
#else
		A(k,i,j) = Minv(i,j); 	  
#endif
	      }
	}

    }
  if(TRACE) cout << "Finished  MatrixInverse" << endl;
}


void MatrixPseudoInverse(VecMat<complex<realtype>>& A)
{
  const int n=A.Nvecs;   // the number of q values
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


// A routine to sum over the trace of sublattice matrices for all q
// it computes sum_q Tr(Aq Bq)
realtype SumTr(VecMat<complex<realtype>>& A,VecMat<complex<realtype>>& B)
{
  if(TRACE) cout << "Starting SumTr" << endl;
  realtype sum=0.;

  for(int k=0; k<A.Nvecs; k++)
    {
      SMatrix<complex<realtype> > temp(A.Nrows,A.Ncols,A[k]);
      temp *= B[k];
      sum += real( tr(temp) );
    }

  if(TRACE) cout << "Ending SumTr" << endl;
  return sum;
}



#endif // MATRIXROUTINES_H
