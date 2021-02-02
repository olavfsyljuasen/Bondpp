#ifndef vecmat_h
#define vecmat_h

#include<iostream>
#include<iomanip>
#include<complex>
#include<vector>
#include<algorithm>

using namespace std;

// first define a small matrix class. It should just be a pointer class unless defined specially to keep contents.

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
 SMatrix(const int nrows,const int ncols,T* in_ptr):Nrows(nrows),Ncols(ncols),ptr(in_ptr),A(0){}
 SMatrix(const int nrows,const int ncols):Nrows(nrows),Ncols(ncols),A(nrows*ncols),ptr(&A[0]){} 
  T& operator()(const int i1,const int i2=0){return *(ptr+smindx(i1,i2));}  
  T& operator[](const int i){return *(ptr+i);}  

  const int Nrows;
  const int Ncols;

  SMatrix& operator=(SMatrix<T>& rhs)
    {
      for(int i=0; i<Nrows*Ncols; i++) (*this)[i]=rhs[i];
      return *this;
    }

  SMatrix& operator=(T* rhs_ptr)
    {
      for(int i=0; i<Nrows*Ncols; i++) (*this)[i]=*(rhs_ptr+i);
      return *this;
    }

  SMatrix& operator+=(SMatrix<T>& rhs) 
    {
      for(int i=0; i<Nrows*Ncols; i++) (*this)[i]+=rhs[i];
      return *this;
    }

  SMatrix& operator-=(SMatrix<T>& rhs) 
    {
      for(int i=0; i<Nrows*Ncols; i++) (*this)[i]-=rhs[i];
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


  SMatrix<T>& operator*=(T* start_ptr) 
    {
      vector<T> temprow(Ncols); // vector for storing temporary result
      for(int i=0; i<Nrows; i++) // assuming Ncols==Nrows
	{
	  for(int j=0; j<Ncols; j++)
	    {
	      T sum(0);
	      for(int k=0; k<Ncols; k++){ sum += (*this)(i,k)*(*(start_ptr+smindx(k,j)));}
	      temprow[j]=sum;
	    }
	  for(int l=0; l<Ncols; l++){ (*this)(i,l)=temprow[l];} // fill in matrix
	}
      return *this;
    }


  SMatrix& operator*=(const complex<realtype>& c)
    {
      for(int i=0; i<Nrows*Ncols; i++) (*this)[i]*=c;
      return *this;
    }

  SMatrix& operator*=(const realtype& c)
    {
      for(int i=0; i<Nrows*Ncols; i++) (*this)[i]*=c;
      return *this;
    }

 private:
  vector<T> A;
  T* ptr;
  inline int smindx(const int i1,const int i2){return i1*Ncols+i2;}
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




//MUST make this class with templates, where the quantity to store is templated.

// A class which is a vector of matrices. (or equivalently a matrix of vectors)
template <class T> 
class VecMat
{
  /*
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
  */
  friend ostream& operator<<(ostream& os,VecMat& m)
  { 
    for(int q=0; q<m.Nvecs; q++)
      {
	os << q << ":" << endl;
	for(int i=0; i<m.Nrows; i++)
	  {
	    for(int j=0; j<m.Ncols; j++){ os << setprecision(16) << m(q,i,j) << " ";}
	    os << endl;
	  }
      }
    m.ContainsNaN("");

    return os;
  }

 public:
 VecMat(const int nvecs,const int nrows,const int ncols=1):Nvecs(nvecs),Nrows(nrows),Ncols(ncols),A(Nvecs*Nrows*Ncols){}
 VecMat(VecMat& c):Nvecs(c.Nvecs),Nrows(c.Nrows),Ncols(c.Ncols),A(Nvecs*Nrows*Ncols)
    {
      T* ptr=c.start();
      for(int i=0; i<A.size(); i++){ A[i]=*(ptr+i);}
    }
  
  const int size() const {return A.size();}
  
  T* start(){return &A[0];} // needed to access start of data-array
  
  T& operator()(const int v,const int i1,const int i2=0){ return A[indx(v,i1,i2)];}  
  
  void SetToZero(){for(int i=0; i<A.size(); i++){A[i]=0;}}
  
  T* operator[](const int i){ return &A[indx(i,0,0)];}  
  //const SMatrix<T>& constant(const int i) const { return A[i];}  
  
  const int Nvecs; // size of vector
  const int Nrows; // number of rows;
  const int Ncols; // number of columns;
  vector<T> A; // the actual storage area

  bool ContainsNaN(string message="")
  {
    for(int i=0; i<A.size(); i++)
      { 
	if( A[i] != A[i] )	
	  {
	    cout << "Nan encountered in " << message << endl;
	    exit(1);
	  }
      }
    return false;
  }


  VecMat& operator=(VecMat<T>& rhs)
    {
      T* ptr=rhs.start();
      if(this == &rhs) return *this; // self-assignment
      for(int i=0; i<A.size(); i++){A[i]= *(ptr+i);}
      return *this;
    }

  VecMat& operator+=(VecMat<T>& rhs) 
    {
      T* ptr=rhs.start();
      for(int i=0; i<A.size(); i++){A[i]+= *(ptr+i);}
      return *this;
    }

  VecMat& operator-=(VecMat<T>& rhs) 
    {
      T* ptr=rhs.start();
      for(int i=0; i<A.size(); i++){A[i]-= *(ptr+i);}
      return *this;
    }

  // sublattice matrix product for each q
  VecMat& operator*=(VecMat<T>& rhs) 
    {
      vector<T> temprow(Ncols); // vector for storing temporary result
      for(int q=0; q<Nvecs; q++)
	{
	  for(int i=0; i<Nrows; i++)
	    {
	      for(int j=0; j<rhs.Ncols; j++)
		{
		  T sum(0);
		  for(int k=0; k<Ncols; k++){ sum += (*this)(q,i,k)*rhs(q,k,j);}
		  temprow[j]=sum;
		}
	      for(int l=0; l<Ncols; l++){ (*this)(q,i,l)=temprow[l];} // fill in matrix
	    }
	}
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
      SMatrix<T> res(Nrows,Ncols);
      for(int q=0; q<Nvecs; q++)
	{
	  T* ptr=A[q];
	  for(int j=0; j<Nrows*Ncols; j++){ res[j] += *(ptr+j);}
	} 
      return res;
    }

  T Sumq(const int m1,const int m2)
    {
      T res(0);
      for(int q=0; q<Nvecs; q++){ res += A[indx(q,m1,m2)];} 
      return res;
    }

private:
inline int indx(const int j,const int m1,const int m2){return (j*Nrows+m1)*Ncols+m2;}
inline int smindx(const int m1,const int m2){return m1*Ncols+m2;}
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
T Tr(VecMat<T>& K)
{
  T sum(0);
  for(int q=0; q<K.Nvecs; q++)
    for(int i=0; i<K.Nrows; i++)
      sum += K(q,i,i);

  return sum;
}



template<class T>
SMatrix<T> Sumq(VecMat<T>& K)
{
  return K.Sumq();
}

template<class T>
T Sumq(VecMat<T>& K,const int m1,const int m2)
{
  return K.Sumq(m1,m2);
}


const realtype epsilon=1.e-15;

template<class T>
void Chomp(VecMat<T>& K)
{
  T* ptr=K.start();
  for(int i=0; i<K.size(); i++)
    {
      T& t(*(ptr+i));
      if(abs(real(t))<epsilon){ t.real(0.);} 
      if(abs(imag(t))<epsilon){ t.imag(0.);} 
    }
}

template<class T>
T FindMax(VecMat<T>& K)
{
  T maxentry;
  realtype maxval=0;

  T* ptr=K.start();
  for(int i=0; i<K.size(); i++)
    {
      T& thisentry(*(ptr+i));
      if(abs(thisentry)>maxval){maxval=abs(thisentry); maxentry=thisentry;}
    }
  return maxentry;
}

template<class T>
T FindMaxImag(VecMat<T>& K)
{
  T maxentry;
  realtype maxval=0;

  T* ptr=K.start();
  for(int i=0; i<K.size(); i++)
    {
      T& thisentry(*(ptr+i));
      if(abs(thisentry.imag()) > maxval){maxval=abs(thisentry.imag()); maxentry=thisentry;}
    }
  return maxentry;
}

template<class T>
void MakeReal(VecMat<T>& K)
{
  T* ptr=K.start();
  for(int i=0; i<K.size(); i++){ (*(ptr+i)).imag(0.);}
}

template<class T>
void ComplexConjugate(VecMat<T>& K)
{
  T* ptr=K.start();
  for(int i=0; i<K.size(); i++){ conj(*(ptr+i));}
}

template<class T>
void MakeRealSymmetric(VecMat<T>& M)
{
  const int Nrows=M.Nrows;
  const int Ncols=M.Ncols;
  
  for(int q=0; q<M.Nvecs; q++)
    {
      for(int s=0; s<Nrows; s++)
	{
	  M(q,s,s).imag(0.); // real diagonal
	}
      
      for(int s1=0; s1<Nrows; s1++)
	for(int s2=s1+1; s2<Ncols; s2++)
	  {
	    M(q,s1,s2).imag(0.);
	    M(q,s2,s1)=M(q,s1,s2);
	  }
    }
}

template<class T>
void MakeRealOnDiagonal(VecMat<T>& M)
{
  const int Nrows=M.Nrows;
  const int Ncols=M.Ncols;
  
  for(int q=0; q<M.Nvecs; q++)
    for(int s=0; s<Nrows; s++)
      M(q,s,s).imag(0.); // real diagonal
}      


template<class T>
void MakeAntiSymmetric(VecMat<T>& M)
{
  const int Nrows=M.Nrows;
  const int Ncols=M.Ncols;
  
  for(int q=0; q<M.Nvecs; q++)
    {
      for(int s=0; s<Nrows; s++)
	{
	  M(q,s,s)=0.;
	}
      
      for(int s1=0; s1<Nrows; s1++)
	for(int s2=s1+1; s2<Ncols; s2++)
	  {
	    M(q,s2,s1)=-M(q,s1,s2);
	  }
    }
}


template<class T>
void MakeHermitian(VecMat<T>& M)
{
  const int Nrows=M.Nrows;
  const int Ncols=M.Ncols;
  
  for(int q=0; q<M.Nvecs; q++)
    {
      for(int s=0; s<Nrows; s++)
	{
	  M(q,s,s).imag(0.); // real diagonal
	}
      
      for(int s1=0; s1<Nrows; s1++)
	for(int s2=s1+1; s2<Ncols; s2++)
	  {
	    M(q,s2,s1) = conj(M(q,s1,s2)); // off-diagonals are cc of each other.
	  }
    }
}

template<class T>
void MakeInversionTransposedSymmetric(BravaisLattice& la,VecMat<T>& M)
{
  const int Nrows=M.Nrows;
  const int Ncols=M.Ncols;
  
  for(int q=0; q<M.Nvecs; q++)
    {
      const int invq=la.GetInversionIndx(q); // the negative q
      
      for(int s1=0; s1<Nrows; s1++)
	for(int s2=0; s2<Ncols; s2++)
	  {
	    if( q==invq ){ M(q,s1,s2).imag(0.);}
	    M(invq,s2,s1) = M(q,s1,s2);
	  }
    }
}




bool IsHermitian(VecMat<complex<realtype>>& M)
{
  if(TRACE) cout << "Checking if Hermitian: ";
  //  bool retval=true;
  const int Nrows=M.Nrows;
  const int Ncols=M.Ncols;

  for(int q=0; q<M.Nvecs; q++)
    {
      
      for(int s=0; s<Nrows; s++)
	if(imag(M(q,s,s)) != 0.){return false;}
  
      for(int s1=0; s1<Nrows; s1++)
	for(int s2=s1+1; s2<Ncols; s2++)
	  {
	    if( M(q,s2,s1) != conj(M(q,s1,s2))){ return false;}
	  }
    }
  return true;

}

bool IsInversionTransposedSymmetric(BravaisLattice& la,VecMat<complex<realtype>>& M)
{
  if(TRACE) cout << "Checking if InversionTransposedSymmetric: ";

  const int Nrows=M.Nrows;
  const int Ncols=M.Ncols;
  
  for(int q=0; q<M.Nvecs; q++)
    {
      const int invq=la.GetInversionIndx(q); // the negative q
      
      for(int m1=0; m1<Nrows; m1++)
	for(int m2=m1; m2<Ncols; m2++)
	  {
	    if( q==invq && M(q,m1,m2).imag() != 0.){ return false;}
	    if( M(invq,m2,m1) != M(q,m1,m2)){ return false;}
	  }
    }
  return true;
}

template<class T>
void AddToDiagonal(VecMat<T>& M,T& value)
{
  for(int q=0; q<M.Nvecs; q++)
    for(int s=0; s<M.Nrows; s++)
      M(q,s,s)+=value;
}


template<class T>
void SubtractFromDiagonal(VecMat<T>& A,T value)
{
  T mvalue=-value;
  AddToDiagonal(A,mvalue);
}


template<class T>
void AddToDiagonal(VecMat<T>& M,vector<T>& value)
{
  for(int q=0; q<M.Nvecs; q++)
    for(int s=0; s<M.Nrows; s++)
      M(q,s,s)+=value[s];
}



template<class T>
void SubtractFromDiagonal(VecMat<T>& A,vector<T> value)
{
  for(int s=0; s<value.size(); s++){value[s]=-value[s];}
  AddToDiagonal(A,value);
}







#endif // vecmat_h
