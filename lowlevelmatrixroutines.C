#include<iostream>
#include "lowlevelmatrixroutines.h"
//This routine must include explicit conversions between
// realtype    <-> matrix_real_type
// complextype <-> matrix_complex_type
// matrix_real_type and matrix_complex_type are possibly higher precision number types used in the matrix routines.
// they must be declared in the lowlevelmatrixroutines.h file

//const bool Checkevals=true;


//#define USEEIGEN
//#define EIGENBOOST // use both EIGEN and BOOST

//We want to use numbers with increased precision in the matrix routines
//we typedef their precision with the keywords matrix_real_type and matrix_complex_type
//for instance
//typedef std::double matrix_real_type;
//typedef std::complex<double> matrix_complex_type;



#ifdef USEEIGEN

#ifdef EIGENBOOST
// Different types of multiprecision backends:
#include <boost/multiprecision/cpp_complex.hpp> // complex numbers
typedef boost::multiprecision::cpp_bin_float_quad   matrix_real_type;
typedef boost::multiprecision::cpp_complex_quad matrix_complex_type;

#include <boost/multiprecision/eigen.hpp>

#else

typedef realtype matrix_real_type;
typedef complextype matrix_complex_type;
#endif

#include<Eigen/Dense>
using namespace Eigen;

#elif defined LAPACKE
//namespace mylapacke{
//#define LAPACK_DISABLE_NAN_CHECK
//#define LAPACK_COMPLEX_CUSTOM
typedef realtype matrix_real_type;
typedef complextype matrix_complex_type;
//#define lapack_complex_float std::complex<float>
#define lapack_complex_double complextype
#include<lapacke.h>
//typedef realtype double;
//typedef complextype lapack_complex_double;
//}

#elif defined MKL_LAPACKE
//#define LAPACK_DISABLE_NAN_CHECK
//#define LAPACK_COMPLEX_CUSTOM
//// #define lapack_complex_float std::complex<float>
// #define lapack_complex_double std::complex<double>
#define lapack_complex_double complextype
#include <mkl_lapacke.h>
typedef realtype matrix_real_type;
typedef complextype matrix_complex_type;
#elif defined MPLAPACK
// use the double-double type:
#include <mpblas_dd.h>
#include <mplapack_dd.h>
#include <mplapack_utils_dd.h>
typedef dd_real matrix_real_type;
typedef dd_complex matrix_complex_type;
#endif


void SmallMatrixMakeHermitian(const int n,vector<complextype>& A)
{
  if(n==1)
    {
      A[0].imag(0.);
    }
  else
    {
      for(int i=0; i<n; i++){ A[i+i*n].imag(0.); }
      
      for(int i=0; i<n; i++)
	for(int j=i+1; j<n; j++)
	  {
	    A[i+j*n] = conj( A[j+i*n] );
	  }
    }
}




#ifdef USEEIGEN
void SmallMatrixInverse(const int n,vector<complextype>& A)
{
  if(n==1) // just take the inverse of each element
    {
      A[0]=realtype(1)/A[0];
    }
  else
    {
      Matrix<matrix_complex_type,Dynamic,Dynamic> M(n,n);
      Matrix<matrix_complex_type,Dynamic,Dynamic> Minv(n,n);
      
      for(int i=0; i<n; i++)
	for(int j=0; j<n; j++)
	  M(i,j)=matrix_complex_type(A[i+j*n]); // col-major order

      Minv=M.inverse();
	  
      for(int i=0; i<n; i++)
	for(int j=0; j<n; j++)
	  A[i+j*n] = complextype(Minv(i,j));
    }  
}

#elif defined LAPACKE || defined MKL_LAPACKE

void SmallMatrixInverse(const int n,vector<complextype>& A)
{
  if(n==1) // just take the inverse of each element
    {
      A[0]=complextype(1.)/A[0];
    }
  else if(n==2)
    {
      matrix_complex_type a=matrix_complex_type(A[0]);
      matrix_complex_type b=matrix_complex_type(A[2]); // col-major order
      matrix_complex_type c=matrix_complex_type(A[1]);
      matrix_complex_type d=matrix_complex_type(A[3]); 

      matrix_complex_type invdet= matrix_complex_type(1.)/(a*d-b*c);
      A[0] = complextype( d*invdet);
      A[2] = complextype(-b*invdet);
      A[1] = complextype(-c*invdet);
      A[3] = complextype( a*invdet) ;
    }
  else
    {      
      int lda=n;
      vector<matrix_complex_type> a(n*n);
      vector<int> ipiv(n);
      int info;
      
      //lapack is fastest using col-major order.
      for(int i=0; i<n*n; i++) a[i]=matrix_complex_type(A[i]); // col-major order
      
      info=LAPACKE_zgetrf(LAPACK_COL_MAJOR, n, n, &a[0], lda, &ipiv[0]); // LU-decomposition
      if(info !=0){cout << "LAPACKE_zgetrf gives info=" << info << endl;}
      info=LAPACKE_zgetri(LAPACK_COL_MAJOR, n,    &a[0], lda, &ipiv[0]); // inversion
      if(info !=0){cout << "LAPACKE_zgetri gives info=" << info << endl;}
      
      for(int i=0; i<n*n; i++) A[i]=complextype(a[i]); // col-major order
    }
}

#elif defined MPLAPACK

void SmallMatrixInverse(const int n,vector<complextype>& A)
{
  if(n==1) // just take the inverse of each element
    {
      A[0]=complextype( matrix_complex_type(1.,0.)/matrix_complex_type(A[0]));
    }
  else if(n==2)
    {
      matrix_complex_type a=matrix_complex_type(A[0]);
      matrix_complex_type b=matrix_complex_type(A[2]); // col-major order
      matrix_complex_type c=matrix_complex_type(A[1]);
      matrix_complex_type d=matrix_complex_type(A[3]); 

      matrix_complex_type invdet= matrix_complex_type(1.)/(a*d-b*c);

      matrix_complex_type ainv=d*invdet;
      matrix_complex_type binv=-b*invdet;
      matrix_complex_type cinv=-c*invdet;
      matrix_complex_type dinv=a*invdet;
      
      A[0] = complextype(cast2double(ainv.real()),cast2double(ainv.imag()));
      A[2] = complextype(cast2double(binv.real()),cast2double(binv.imag()));
      A[1] = complextype(cast2double(cinv.real()),cast2double(cinv.imag()));
      A[3] = complextype(cast2double(dinv.real()),cast2double(dinv.imag()));
    }
  else
    {      
      mplapackint lda(n),lwork(n),info(0);
      vector<matrix_complex_type> a(n*n),work(n);
      vector<mplapackint> ipiv(n);

      for(int i=0; i<n*n; i++) a[i]=matrix_complex_type(A[i]); // col-major order

      Cgetrf(n, n, &a[0], lda, &ipiv[0],info); // LU-decomposition
      if(info !=0){cout << "Cgetrf gives info=" << info << endl;}

      //      int *n, singlecomplex *a, int *lda, int *ipiv, singlecomplex *work, int *lwork, int *info);
      
      Cgetri(n,    &a[0], lda, &ipiv[0],&work[0],lwork,info); // inversion
      if(info !=0){cout << "Cgetri gives info=" << info << endl;}
      
      for(int i=0; i<n*n; i++) A[i]=complextype(cast2double(a[i].real()),cast2double(a[i].imag())); 
    }
}
#endif	  


#ifdef USEEIGEN
realtype SmallHermitianMatrixMinEigenvalue(const int n,vector<complextype>& A)
{
  if(n==1) 
    {
      return real(A[0]);
    }
  else if(n==2)
    {
      matrix_complex_type a = matrix_complex_type(A[0]);
      matrix_complex_type b = matrix_complex_type(A[2]);
      matrix_complex_type c = matrix_complex_type(A[1]);
      matrix_complex_type d = matrix_complex_type(A[3]);
      
      matrix_real_type rad=real( (a-d)*(a-d)+matrix_complex_type(4)*b*c); 
      matrix_real_type lambda = matrix_real_type(0.5)*(real(a+d) - sqrt(rad));

      return realtype(lambda);      
    }
  else
    {
      Matrix<matrix_complex_type,Dynamic,Dynamic> M(n,n); // use the dynamic type h
      
      for(int i=0; i<n; i++)
	for(int j=0; j<n; j++)
	  M(i,j)=matrix_complex_type(A[i+j*n]); // col-major order
	  
      SelfAdjointEigenSolver<Matrix<matrix_complex_type,Dynamic,Dynamic> > es(M,EigenvaluesOnly);
      
      realtype lambda = realtype(es.eigenvalues()[0]);
      return lambda;
    }
}

#elif defined LAPACKE || defined MKL_LAPACKE      

realtype SmallHermitianMatrixMinEigenvalue(const int n,vector<complextype>& A)
{
  if(n==1) 
    {
      return real(A[0]);
    }
  else if(n==2)
    {
      matrix_complex_type a = matrix_complex_type(A[0]);
      matrix_complex_type b = matrix_complex_type(A[2]);
      matrix_complex_type c = matrix_complex_type(A[1]);
      matrix_complex_type d = matrix_complex_type(A[3]);

      matrix_real_type rad=real( (a-d)*(a-d)+matrix_complex_type(4)*b*c); 
      matrix_real_type lambda = matrix_real_type(0.5)*(real(a+d) - sqrt(rad));

      return realtype(lambda);      
    }
  else
    {
      char jobz='N'; // no left eigenvectors
      char uplo='L'; // upper triangular is uploaded
      vector<matrix_complex_type> a(n*n);
      lapack_int lda(n),info(0);
      vector<matrix_real_type> w(n); // the eigenvalues 

      for(int i=0; i<n*n; i++){ a[i]=matrix_complex_type(A[i]); }
      
      info=LAPACKE_zheev(LAPACK_COL_MAJOR,jobz, uplo,n,&a[0],lda, &w[0]);
      if(info !=0){ cout << "Warning LAPACKE_zheev exited with info=" << info << endl;}

      realtype lambda=realtype(w[0]);
      return lambda;
    }
}

#elif defined MPLAPACK
realtype SmallHermitianMatrixMinEigenvalue(const int n,vector<complextype>& A)
{
  if(n==1) 
    {
      return real(A[0]);
    }
  else if(n==2)
    {
      matrix_complex_type a = matrix_complex_type(A[0]);
      matrix_complex_type b = matrix_complex_type(A[2]);
      matrix_complex_type c = matrix_complex_type(A[1]);
      matrix_complex_type d = matrix_complex_type(A[3]);
      
      matrix_complex_type rad= (a-d)*(a-d)+matrix_complex_type(4)*b*c;
      matrix_complex_type minev= matrix_complex_type(0.5)*(a+d - sqrt(rad));
      realtype lambda = cast2double(minev.real());

      return lambda;      
    }
  else
    {
      char jobz='N'; // no left eigenvectors
      char uplo='L'; // upper triangular is uploaded
      vector<matrix_complex_type> a(n*n);
      mplapackint lda(n),info(0);
      vector<matrix_real_type> w(n); // the eigenvalues 
      mplapackint lwork(2*n-1);
      vector<matrix_complex_type> work(lwork);            
      vector<matrix_real_type> rwork(3*n-2);      
      
      // use column-major order and upload the lower triangular matrix only
      /*
      for(int j=0; j<n; j++)
	for(int i=j; i<n; i++)
	  { a[i+j*n]=matrix_complex_type(A[i+j*n]); }
      */
      for(int i=0; i<n*n; i++){ a[i]=matrix_complex_type(A[i]); }

      Cheev(&jobz, &uplo,n,&a[0],lda, &w[0], &work[0],lwork,&rwork[0],info);
      if(info !=0){ cout << "Warning Cheev exited with info=" << info << endl;}

      realtype lambda=cast2double(w[0]);
      return lambda;
    }
}




#endif	  




#ifdef USEEIGEN
realtype SmallHermitianMatrixDeterminant(const int n,vector<complextype>& A)
{
  if(n==1) 
    {
      return real(A[0]);
    }
  else if(n==2)
    {
      matrix_complex_type a = matrix_complex_type(A[0]);
      matrix_complex_type b = matrix_complex_type(A[2]);
      matrix_complex_type c = matrix_complex_type(A[1]);
      matrix_complex_type d = matrix_complex_type(A[3]);
      
      matrix_real_type det =real(a*d-b*c);
      
      return det;      
    }
  else
    {
      Matrix<matrix_complex_type,Dynamic,Dynamic> M(n,n); // use the dynamic type h
      
      for(int i=0; i<n; i++)
	for(int j=0; j<n; j++)
	  M(i,j)=matrix_complex_type(A[i+j*n]); // col-major order

      complextype det=complextype(M.determinant());
      //      cout << "determinant " << det << endl;
      return real(det);
    }
}

#elif defined LAPACKE || defined MKL_LAPACKE      

realtype SmallHermitianMatrixDeterminant(const int n,vector<complextype>& A)
{
  if(n==1) 
    {
      return real(A[0]);
    }
  else if(n==2)
    {
      matrix_complex_type a = matrix_complex_type(A[0]);
      matrix_complex_type b = matrix_complex_type(A[2]);
      matrix_complex_type c = matrix_complex_type(A[1]);
      matrix_complex_type d = matrix_complex_type(A[3]);
      
      matrix_real_type det =real(a*d-b*c);
      
      return realtype(det);      
    }
  else
    {
      char jobz='N'; // no left eigenvectors
      char uplo='L'; // upper triangular is uploaded
      vector<matrix_complex_type> a(n*n);
      lapack_int lda(n),info(0);
      vector<matrix_real_type> w(n); // the eigenvalues 

      for(int i=0; i<n*n; i++){ a[i]=matrix_complex_type(A[i]); }
      
      info=LAPACKE_zheev(LAPACK_COL_MAJOR,jobz, uplo,n,&a[0],lda, &w[0]);
      if(info !=0){ cout << "Warning LAPACKE_zheev exited with info=" << info << endl;}

      matrix_real_type det=matrix_real_type(1.);      
      for(int i=0; i<n; i++){det *= w[i];}

      return realtype(det);      

      /*
      char uplo='L';
      int lda(n);
      vector<matrix_complex_type> a(n*n);
      vector<int> ipiv(n); // permutations of rows

      for(int i=0; i<n*n; i++){ a[i]=matrix_complex_type(A[i]); }

      LAPACKE_zhetrf(LAPACK_COL_MAJOR, uplo, n, &a[0], lda, &ipiv[0]); // LU-decomposition of hermitian m.


      for(int i=0; i<n; i++){ det *= real(a[i*n+i]);}
      for(int i=0; i<n; i++){ if(ipiv[i] != i+1) det *= matrix_real_type(-1.);}

      return realtype(det);      
      */
    }
}

#elif defined MPLAPACK 
realtype SmallHermitianMatrixDeterminant(const int n,vector<complextype>& A)
{
  if(n==1) 
    {
      return real(A[0]);
    }
  else if(n==2)
    {
      matrix_complex_type a = matrix_complex_type(A[0]);
      matrix_complex_type b = matrix_complex_type(A[2]);
      matrix_complex_type c = matrix_complex_type(A[1]);
      matrix_complex_type d = matrix_complex_type(A[3]);
      
      matrix_complex_type det = a*d-b*c;
      
      return cast2double(det.real());      
    }
  else
    {
      char jobz='N'; // no left eigenvectors
      char uplo='L'; // upper triangular is uploaded
      vector<matrix_complex_type> a(n*n);
      mplapackint lda(n),info(0);
      vector<matrix_real_type> w(n); // the eigenvalues 
      mplapackint lwork(2*n-1);
      vector<matrix_complex_type> work(lwork);            
      vector<matrix_real_type> rwork(3*n-2);      
      
      // use column-major order and upload the lower triangular matrix only
      /*
      for(int j=0; j<n; j++)
	for(int i=j; i<n; i++)
	  { a[i+j*n]=matrix_complex_type(A[i+j*n]); }
      */
      for(int i=0; i<n*n; i++){ a[i]=matrix_complex_type(A[i]); }

      Cheev(&jobz, &uplo,n,&a[0],lda, &w[0], &work[0],lwork,&rwork[0],info);
      if(info !=0){ cout << "Warning Cheev exited with info=" << info << endl;}

      matrix_real_type det=1.;
      for(int i=0; i<n; i++){ det *= w[i];}
      
      return cast2double(det);      
    }
}
#endif



// solving the linear system of equations using Moore-Penrose inversion
#ifdef USEEIGEN

vector<realtype> SmallMatrixPseudoSolve(vector<realtype>& Ain, vector<realtype>& Bin,const int M, const int N)
{
  //  if(TRACE)  cout << "In SmallMatrixPseudoSolve" << endl;
  const realtype tolerance=1e-18;
    
  Matrix<matrix_real_type,Dynamic,Dynamic> A(M,N);
  Matrix<matrix_real_type,Dynamic,Dynamic> Ainv(N,M);
  Matrix<matrix_real_type,Dynamic,Dynamic> singularValuesInv(N,M);
  
  for(int i=0; i<M; i++)
    for(int j=0; j<N; j++)
      A(i,j)=matrix_real_type(Ain[i*N+j]); // row-major order
  
      
  auto svd = A.jacobiSvd(Eigen::ComputeFullU | Eigen::ComputeFullV);
  const auto &singularValues = svd.singularValues();
  
  singularValuesInv.setZero();
  
  for (unsigned int i = 0; i < singularValues.size(); ++i)
    {
      if (singularValues(i) > tolerance)
	{
	  singularValuesInv(i, i) = realtype(1.) / singularValues(i);
	}
      else
	{
	  cout << "Setting sing val " << i << " to 0" << endl;
	  singularValuesInv(i, i) = 0;
	}
    }

  
  Ainv = svd.matrixV() * singularValuesInv * svd.matrixU().transpose();
  
  Matrix<matrix_real_type,Dynamic,Dynamic> B(M,1);
  Matrix<matrix_real_type,Dynamic,Dynamic> ans(N,1);
  
  for(int i=0; i<M; i++) B(i)=matrix_real_type(Bin[i]);
  
  ans = Ainv * B;
  
  vector<realtype> answer(N);
  for(int i=0; i<N; i++){ answer[i]=realtype(ans(i)); }
  
  return answer;
}
#endif


/*
int main()
{
  // nothing to do.
}
*/



