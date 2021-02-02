#ifndef vecmat_h
#define vecmat_h

#include<iostream>
#include<iomanip>
#include<complex>
#include<vector>
#include<algorithm>

using namespace std;

// first define a small matrix class
/*
template <class T>
class SMatrix
{
  friend ostream& operator<<(ostream& os,SMatrix& m)
  { 
    for(int i=0; i<m.Nrows; i++)
      {
	for(int j=0; j<m.Ncols; j++){ os << setprecision(35) << m(i,j) << " ";}
	os << endl;
      }
    return os;
  }
 public:
 SMatrix(const int ncols,const int nrows):A(ncols*nrows),Ncols(ncols),Nrows(nrows),Nmat(ncols*nrows){}
 SMatrix(const SMatrix& c):Ncols(c.Ncols),Nrows(c.Nrows),Nmat(Ncols*Nrows),A(c.size())
  {
    copy(c.A.begin(), c.A.end(), A.begin());
  }
  int size() const {return A.size();}

  T* start(){return &A[0];} // probably not needed, but just in case

  void SetToZero(){ for(int i=0; i<A.size(); i++){ A[i]=0;}}

  T& operator()(const int i1,const int i2=0){ return A[smindx(i1,i2)];}  
  T& operator[](const int i){ return A[i];}  

  SMatrix& operator=(SMatrix<T>& rhs)
    {
      copy(rhs.A.begin(), rhs.A.end(), A.begin());
      return *this;
    }

  SMatrix& operator+=(SMatrix<T>& rhs) 
    {
      for(int i=0; i<A.size(); i++){A[i]+=rhs[i];}
      return *this;
    }

  SMatrix& operator-=(SMatrix<T>& rhs) 
    {
      for(int i=0; i<A.size(); i++){A[i]-=rhs[i];}
      return *this;
    }

  SMatrix<T>& operator*=(SMatrix<T> rhs) 
    {
      vector<T> temprow(Ncols); // vector for storing temporary result
      for(int i=0; i<Nrows; i++)
	{
	  for(int j=0; j<rhs.Ncols; j++)
	    {
	      T sum(0);
	      for(int k=0; k<Ncols; k++){ sum += (*this)(i,k)*rhs(k,j);}
	      temprow[j]=sum;
	    }
	  for(int l=0; l<Ncols; l++){ (*this)(i,l)=temprow[l];} // fill in matrix
	}
      return *this;
    }


  SMatrix& operator*=(const complex<realtype>& c)
    {
      for(int i=0; i<A.size(); i++){ A[i]*=c;} 
      return *this;
    }

  
  SMatrix& operator*=(const realtype& c)
    {
      for(int i=0; i<A.size(); i++){ A[i]*=c;} 
      return *this;
    }
  
  
  const int Ncols; // number of columns;
  const int Nrows; // number of rows;
  const int Nmat; // number of matrix elements = Ncols*Nrows;
  
  vector<T> A; // the actual storage area

 private:
  
  int smindx(const int m1,const int m2){return m1*Nrows+m2;} // row-first
  
};

// matrix product
template<class T>
SMatrix<T> operator*(SMatrix<T>& lhs,SMatrix<T>& rhs)
{
  SMatrix<T> res(lhs.Nrows,rhs.Ncols);
  
  for(int i=0; i<lhs.Nrows; i++)
    for(int j=0; j<rhs.Ncols; j++)
      {
	T sum(0);
	for(int k=0; k<lhs.Ncols; k++){ sum += lhs(i,k)*rhs(k,j);}
	res(i,j)=sum;
      }
  return res;
}

template<class T>
SMatrix<T> operator*(const T& c,const SMatrix<T>& rhs)
{
  SMatrix<T> res(rhs.Nrows,rhs.Ncols);
  for(int i=0; i<rhs.size(); i++){ res[i] = rhs[i]*c;}
  return res;
}

template<class T>
T tr(SMatrix<T>& M)
{
  T sum(0);
  for(int i=0; i<M.Nrows; i++){ sum+= M(i,i);}
  return sum;
}


template<class T>
SMatrix<T> operator+(const SMatrix<T>& lhs,const SMatrix<T>& rhs)
{
  SMatrix<T> res(lhs.Nrows,lhs.Ncols);
  
  for(int i=0; i<lhs.Nrows; i++)
    for(int j=0; j<lhs.Ncols; j++)
      {
	res(i,j)=lhs(i,j)+rhs(i,j);
      }
  return res;
}

template<class T>
SMatrix<T> operator-(const SMatrix<T>& lhs,const SMatrix<T>& rhs)
{
  SMatrix<T> res(lhs.Nrows,lhs.Ncols);
  
  for(int i=0; i<lhs.Nrows; i++)
    for(int j=0; j<lhs.Ncols; j++)
      {
	res(i,j)=lhs(i,j)-rhs(i,j);
      }
  return res;
}



*/



//MUST make this class with templates, where the quantity to store is templated.
// A class which is a vector of matrices. (or equivalently a matrix of vectors)
#include <boost/numeric/ublas/storage.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>

template <class T>
using SMatrix = boost::numeric::ublas::matrix<T, boost::numeric::ublas::row_major,std::vector<T> >;


// matrix product
template<class T>
SMatrix<T> operator*(SMatrix<T>& lhs,SMatrix<T>& rhs)
{
#ifndef NDEBUG
  if( lhs.size2() != rhs.size1() ){ cout << "Error: Matrix dimensions wrong for matrix product" << endl; exit(1);}
#endif

  SMatrix<T> res(lhs.size1(),rhs.size2());
  
  for(int i=0; i<lhs.size1(); i++)
    for(int j=0; j<rhs.size2(); j++)
      {
	T sum(0);
	for(int k=0; k<lhs.size2(); k++){ sum += lhs(i,k)*rhs(k,j);}
	res(i,j)=sum;
      }
  return res;
}

template<class T>
SMatrix<T>& operator*=(SMatrix<T>& lhs,SMatrix<T>& rhs) 
{
#ifndef NDEBUG
  if( lhs.size2() != rhs.size1() ){ cout << "Error: Matrix dimensions wrong for matrix product" << endl; exit(1);}
#endif

  vector<T> temprow(lhs.size2()); // vector for storing temporary result
  for(int i=0; i<lhs.size1(); i++)
    {
      for(int j=0; j<rhs.size2(); j++)
	{
	  T sum(0);
	  for(int k=0; k<lhs.size2(); k++){ sum += lhs(i,k)*rhs(k,j);}
	  temprow[j]=sum;
	}
      for(int l=0; l<lhs.size2(); l++){ lhs(i,l)=temprow[l];} // fill in matrix
    }
  return lhs;
}


template<class T>
SMatrix<T> conj(SMatrix<T>& M)
{
  SMatrix<T> res(M.Nrows,M.Ncols);
  
  for(int i=0; i<M.Nrows; i++)
    for(int j=0; j<M.Ncols; j++)
      {
	T temp=M(i,j);
	res(i,j)=conj(temp);
      }
  return res;
}





template <class T> 
class VecMat
{
  friend ostream& operator<<(ostream& os,VecMat& c)
  { 
    bool printdots=false;
    int end=c.size();
    if(c.size() > NELEMENTSTOPRINT){end=NELEMENTSTOPRINT; printdots=true;}
    T* start=c.start();
    for(int i=0; i<end; i++) 
      os << start[i] << " ";
    if(printdots) os << " ... " << start[end-1];
    return os;
  }
 public:
 VecMat(const int nvec,const int nrows,const int ncols=1):Nvec(nvec),Nrows(nrows),Ncols(ncols),A(nvec,SMatrix<T>(nrows,ncols)){}
 VecMat(VecMat& c):Nvec(c.Nvec),Nrows(c.Nrows),Ncols(c.Ncols),A(c.size(),SMatrix<T>(c.Nrows,c.Ncols))
    {
      for(int i=0; i<A.size(); i++){ A[i]=c[i];}
    }
  
  const int size() const {return Nvec;}
  //  int size() {return Nvec;}
  
  T* start(){return &A[0](0,0);} // needed for pointing to storage are
  
  T& operator()(const int v,const int i1,const int i2=0){ return A[v](i1,i2);}  

  void SetToZero(){for(int i=0; i<A.size(); i++){A[i].clear();}}
  
  SMatrix<T>& operator[](const int i){ return A[i];}  
  
  int Nvec; // size of vector
  int Nrows; // number of rows;
  int Ncols; // number of columns;


  vector<SMatrix<T>> A; // the actual storage area

  VecMat& operator=(VecMat<T>& rhs)
    {
      if(this == &rhs) return *this; // self-assignment
      for(int i=0; i<A.size(); i++){A[i]=rhs[i];}
      return *this;
    }

  VecMat& operator+=(VecMat<T>& rhs) 
    {
      for(int i=0; i<A.size(); i++){A[i]+=rhs[i];}
      return *this;
    }

  VecMat& operator-=(VecMat<T>& rhs) 
    {
      for(int i=0; i<A.size(); i++){A[i]-=rhs[i];}
      return *this;
    }

  // sublattice matrix product for each q
  VecMat& operator*=(VecMat<T>& rhs) 
    {
      for(int i=0; i<A.size(); i++){A[i]*=rhs[i];}
      return *this;
    }

  VecMat& operator*=(const complex<realtype>& c)
    {
      for(int i=0; i<A.size(); i++){ A[i]*=c;} 
      return *this;
    }

  
  VecMat& operator*=(const realtype& c)
    {
      for(int i=0; i<A.size(); i++){ A[i]*=c;} 
      return *this;
    }
 

  SMatrix<T> Sumq()
    {
      SMatrix<T> res(A[0]);
      for(int i=1; i<A.size(); i++){ res += A[i];} 
      return res;
    }

  T Sumq(const int m1,const int m2)
    {
      T res(A[0](m1,m2));
      for(int i=1; i<A.size(); i++){ res += A[i](m1,m2);} 
      return res;
    }
 
private:
};


template<class T>
VecMat<T> operator+(VecMat<T>& lhs,VecMat<T>& rhs)
{
  VecMat<T> res(lhs);
  res += rhs;
  return res;
}

template<class T>
VecMat<T> operator-(VecMat<T>& lhs,VecMat<T>& rhs)
{
  VecMat<T> res(lhs);
  res -= rhs;
  return res;
}

template<class T>
VecMat<T> operator*(VecMat<T>& lhs,VecMat<T>& rhs)
{
  VecMat<T> res(lhs);
  res *= rhs;
  return res;
}

template<class T>
VecMat<T> operator*(VecMat<T>& lhs,const complex<realtype>& c)
{
  VecMat<T> res(lhs);
  res *= c;
  return res;
}

template<class T>
VecMat<T> operator*(VecMat<T>& lhs,const realtype& c)
{
  VecMat<T> res(lhs);
  res *= c;
  return res;
}


template<class T>
T tr(SMatrix<T>& M)
{
  T sum(0);
  for(int i=0; i<M.size1(); i++){ sum+= M(i,i);}
  return sum;
}


template<class T>
T Tr(VecMat<T>& K)
{
  T sum(0);
  for(int i=0; i<K.size(); i++){ sum+= tr(K[i]);}
  return sum;
}


template<class T>
SMatrix<T>& Sumq(VecMat<T>& K)
{
  return K.Sumq();
}

template<class T>
T Sumq(VecMat<T>& K,const int m1,const int m2)
{
  return K.Sumq(m1,m2);
}


#endif // vecmat_h
