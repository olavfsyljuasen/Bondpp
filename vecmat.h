#ifndef vecmat_h
#define vecmat_h

#include<string>
#include<iostream>
#include<iomanip>
#include<complex>
#include<vector>
#include<algorithm>

using namespace std;

const realtype epsilon=1.e-15;   // used to repair matrices, chomp.
const realtype sensitivity=1.e-10; // used in SanityCheck

// first define a small matrix class. It should just be a pointer class unless defined specially to keep contents.


/*
template <class T>
class SMatrix
{
  friend ostream& operator<<(ostream& os,SMatrix& m)
  { 
    for(int i=0; i<m.Nrows; i++)
      {
	for(int j=0; j<m.Ncols; j++){ os << setprecision(16) << m(i,j) << " ";}
	os << endl;
      }
    return os;
  }

  friend ostream& operator<<(ostream& os,SMatrix m)
  { 
    for(int i=0; i<m.Nrows; i++)
      {
	for(int j=0; j<m.Ncols; j++){ os << setprecision(16) << m(i,j) << " ";}
	os << endl;
      }
    return os;
  }

 public:
 SMatrix(const int nrows,const int  ncols,T* in_ptr):Nrows(nrows),Ncols(ncols),A(0),ptr(in_ptr){}
 SMatrix(const int nrows,const int  ncols):Nrows(nrows),Ncols(ncols),A(nrows*ncols),ptr(&A[0]){} 
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
      //      if(Nrows==1 && Ncols==1){ (*ptr) *= rhs(0,0); return *this;}

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

  SMatrix& Transpose()
    {
      for(int i=0; i<Nrows; i++)
	for(int j=i+1; j<Ncols; j++)
	  {
	    const int ij=smindx(i,j);
	    const int ji=smindx(j,i);
	    
	    T temp=A[ij];
	    A[ij]=A[ji];
	    A[ji]=temp;
	  }
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
*/


template <class T,int Nrows,int Ncols>
class SMatrix
{
  friend ostream& operator<<(ostream& os,SMatrix& m)
  { 
    for(int i=0; i<Nrows; i++)
      {
	for(int j=0; j<Ncols; j++){ os << setprecision(16) << m(i,j) << " ";}
	os << endl;
      }
    return os;
  }

  friend ostream& operator<<(ostream& os,SMatrix m)
  { 
    for(int i=0; i<Nrows; i++)
      {
	for(int j=0; j<Ncols; j++){ os << setprecision(16) << m(i,j) << " ";}
	os << endl;
      }
    return os;
  }

 public:
 SMatrix(T* in_ptr):A(0),ptr(in_ptr){}
 SMatrix():A(Nrows*Ncols),ptr(&A[0]){} 
  T& operator()(const int i1,const int i2=0){return *(ptr+smindx(i1,i2));}  
  T& operator[](const int i){return *(ptr+i);}  

  SMatrix& operator=(SMatrix& rhs)
    {
      for(int i=0; i<Nrows*Ncols; i++) (*this)[i]=rhs[i];
      return *this;
    }

  SMatrix& operator=(T* rhs_ptr)
    {
      for(int i=0; i<Nrows*Ncols; i++) (*this)[i]=*(rhs_ptr+i);
      return *this;
    }

  SMatrix& operator+=(SMatrix& rhs) 
    {
      for(int i=0; i<Nrows*Ncols; i++) (*this)[i]+=rhs[i];
      return *this;
    }

  SMatrix& operator-=(SMatrix& rhs) 
    {
      for(int i=0; i<Nrows*Ncols; i++) (*this)[i]-=rhs[i];
      return *this;
    }

  SMatrix& operator*=(SMatrix& rhs) // multiplication of square matrices 
    {
      //      if(Nrows==1 && Ncols==1){ (*ptr) *= rhs(0,0); return *this;}

      vector<T> temprow(Ncols); // vector for storing temporary result
      for(int i=0; i<Nrows; i++)
	{
	  for(int j=0; j<Ncols; j++)
	    {
	      T sum(0);
	      for(int k=0; k<Ncols; k++){ sum += (*this)(i,k)*rhs(k,j);}
	      temprow[j]=sum;
	    }
	  for(int l=0; l<Ncols; l++){ (*this)(i,l)=temprow[l];} // fill in matrix
	}
      return *this;
    }


  SMatrix& operator*=(T* start_ptr) 
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

  SMatrix& Transpose()
    {
      for(int i=0; i<Nrows; i++)
	for(int j=i+1; j<Ncols; j++)
	  {
	    const int ij=smindx(i,j);
	    const int ji=smindx(j,i);
	    
	    T temp=A[ij];
	    A[ij]=A[ji];
	    A[ji]=temp;
	  }
      return *this;
    }
  
 private:
  vector<T> A;
  T* ptr;
  inline int smindx(const int i1,const int i2){return i1*Ncols+i2;}
};
  

// matrix product of square matrices
template<class T,int Nrows,int Ncols>
  SMatrix<T,Nrows,Ncols> operator*(SMatrix<T,Nrows,Ncols>& lhs,SMatrix<T,Nrows,Ncols>& rhs)
{
  SMatrix<T,Nrows,Ncols> res;
  
  for(int i=0; i<Nrows; i++)
    for(int j=0; j<Ncols; j++)
      {
	T sum(0);
	for(int k=0; k<Ncols; k++){ sum += lhs(i,k)*rhs(k,j);}
	res(i,j)=sum;
      }
  return res;
}

template<class T,int Nrows,int Ncols>
  SMatrix<T,Nrows,Ncols> operator*(const T& c,const SMatrix<T,Nrows,Ncols>& rhs)
{
  SMatrix<T,Nrows,Ncols> res;
  for(int i=0; i<Nrows*Ncols; i++){ res[i] = rhs[i]*c;}
  return res;
}

template<class T,int Nrows,int Ncols>
  T tr(SMatrix<T,Nrows,Ncols>& M)
{
  T sum(0);
  for(int i=0; i<Nrows; i++){ sum+= M(i,i);}
  return sum;
}


template<class T,int Nrows,int Ncols>
  SMatrix<T,Nrows,Ncols> operator+(const SMatrix<T,Nrows,Ncols>& lhs,const SMatrix<T,Nrows,Ncols>& rhs)
{
  SMatrix<T,Nrows,Ncols> res;
  
  for(int i=0; i<Nrows; i++)
    for(int j=0; j<Ncols; j++)
      {
	res(i,j)=lhs(i,j)+rhs(i,j);
      }
  return res;
}

template<class T,int Nrows,int Ncols>
SMatrix<T,Nrows,Ncols> operator-(const SMatrix<T,Nrows,Ncols>& lhs,const SMatrix<T,Nrows,Ncols>& rhs)
{
  SMatrix<T,Nrows,Ncols> res;
  
  for(int i=0; i<Nrows; i++)
    for(int j=0; j<Ncols; j++)
      {
	res(i,j)=lhs(i,j)-rhs(i,j);
      }
  return res;
}


template<class T,int Nrows,int Ncols>
  SMatrix<T,Nrows,Ncols> conj(SMatrix<T,Nrows,Ncols>& M)
{
  SMatrix<T,Nrows,Ncols> res;
  
  for(int i=0; i<Nrows; i++)
    for(int j=0; j<Ncols; j++)
      {
	T temp=M(i,j);
	res(i,j)=conj(temp);
      }
  return res;
}




/*
//MUST make this class with templates, where the quantity to store is templated.

// A class which is a vector of matrices. (or equivalently a matrix of vectors)
template <class T> 
class VecMat
{
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
      for(unsigned int i=0; i<A.size(); i++){ A[i]=*(ptr+i);}
    }
  
  unsigned int size() const {return A.size();}
  
  T* start(){return &A[0];} // needed to access start of data-array
  
  T& operator()(const int v,const int i1,const int i2=0){ return A[indx(v,i1,i2)];}  
  
  void SetToZero(){for(unsigned int i=0; i<A.size(); i++){A[i]=0;}}
  
  T* operator[](const int i){ return &A[indx(i,0,0)];}  

  SMatrix<T> qcomp(const int q){return SMatrix<T>(Nrows,Ncols,&A[indx(q,0,0)]);} // only used for printing

  
  const int Nvecs; // size of vector
  const int Nrows; // number of rows;
  const int Ncols; // number of columns;
  vector<T> A; // the actual storage area

  bool ContainsNaN(string message="")
  {
    for(unsigned int i=0; i<A.size(); i++)
      { 
	if( A[i] != A[i] )	
	  {
	    cout << "Nan encountered in " << message << endl;
	    return true;
	  }
      }
    return false;
  }


  VecMat& operator=(VecMat<T>& rhs)
    {
      T* ptr=rhs.start();
      if(this == &rhs) return *this; // self-assignment
      for(unsigned int i=0; i<A.size(); i++){A[i]= *(ptr+i);}
      return *this;
    }

  VecMat& operator+=(VecMat<T>& rhs) 
    {
      T* ptr=rhs.start();
      for(unsigned int i=0; i<A.size(); i++){A[i]+= *(ptr+i);}
      return *this;
    }

  VecMat& operator-=(VecMat<T>& rhs) 
    {
      T* ptr=rhs.start();
      for(unsigned int i=0; i<A.size(); i++){A[i]-= *(ptr+i);}
      return *this;
    }

  // sublattice matrix product for each q
  VecMat& operator*=(VecMat<T>& rhs) 
    {
      vector<T> temprow(Ncols); // vector for storing temporary result
      for(int q=0; q<Nvecs; q++)
	{
	  for(unsigned int i=0; i<Nrows; i++)
	    {
	      for(unsigned int j=0; j<rhs.Ncols; j++)
		{
		  T sum(0);
		  for(int k=0; k<Ncols; k++){ sum += (*this)(q,i,k)*rhs(q,k,j);}
		  temprow[j]=sum;
		}
	      for(unsigned int l=0; l<Ncols; l++){ (*this)(q,i,l)=temprow[l];} // fill in matrix
	    }
	}
      return *this;
    }

  VecMat& operator*=(const complex<realtype>& c)
    {
      for(unsigned int i=0; i<A.size(); i++){ A[i]*=c;} 
      return *this;
    }

  
  VecMat& operator*=(const realtype& c)
    {
      for(unsigned int i=0; i<A.size(); i++){ A[i]*=c;} 
      return *this;
    }
 

  SMatrix<T> Sumq()
    {
      SMatrix<T> res(Nrows,Ncols);
      for(int q=0; q<Nvecs; q++)
	{
	  T* ptr=A[q];
	  for(unsigned int j=0; j<Nrows*Ncols; j++){ res[j] += *(ptr+j);}
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
  for(unsigned int q=0; q<K.Nvecs; q++)
    for(unsigned int i=0; i<K.Nrows; i++)
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




template<class T>
void Chomp(VecMat<T>& K)
{
  T* ptr=K.start();
  for(unsigned int i=0; i<K.size(); i++)
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
  for(unsigned int i=0; i<K.size(); i++)
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
  for(unsigned int i=0; i<K.size(); i++)
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
  for(unsigned int i=0; i<K.size(); i++)
    {
      if(TRACE)
	{
	  realtype ipart=imag(*(ptr+i));
	  if( fabs(ipart) > sensitivity)
	    {
	      cout << "Warning, MakeReal sets entry " << ipart << " to 0" << endl;
	    }
	}
      (*(ptr+i)).imag(0.);
    } 
}

template<class T>
void ComplexConjugate(VecMat<T>& K)
{
  T* ptr=K.start();
  for(unsigned int i=0; i<K.size(); i++){ conj(*(ptr+i));}
}

template<class T>
void MakeRealSymmetric(VecMat<T>& M)
{
  const int Nrows=M.Nrows;
  const int Ncols=M.Ncols;
  
  for(unsigned int q=0; q<M.Nvecs; q++)
    {
      for(unsigned int s=0; s<Nrows; s++)
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
  for(unsigned int s=0; s<value.size(); s++){value[s]=-value[s];}
  AddToDiagonal(A,value);
}
*/


//MUST make this class with templates, where the quantity to store is templated.

// A class which is a vector of matrices. (or equivalently a matrix of vectors)
template <class T,int Nrows,int Ncols> 
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
    for(int q=0; q< m.size()/(Nrows*Ncols); q++)
      {
	os << q << ":" << endl;
	for(int i=0; i<Nrows; i++)
	  {
	    for(int j=0; j<Ncols; j++){ os << setprecision(16) << m(q,i,j) << " ";}
	    os << endl;
	  }
      }
    m.ContainsNaN("");

    return os;
  }

 public:
  VecMat(const int nvecs):A(nvecs*Nrows*Ncols){}
  VecMat(VecMat& c):A(c.size())
  {
    T* ptr=c.start();
    for(unsigned int i=0; i<A.size(); i++){ A[i]=*(ptr+i);}
  }
  
  unsigned int size() const {return A.size();}

  int Nvecs(){ return A.size()/(Nrows*Ncols);}
  
  T* start(){return &A[0];} // needed to access start of data-array
  
  T& operator()(const int v,const int i1,const int i2=0){ return A[indx(v,i1,i2)];}  
  
  void SetToZero(){for(unsigned int i=0; i<A.size(); i++){A[i]=0;}}
  
  T* operator[](const int i){ return &A[indx(i,0,0)];}  
  
  SMatrix<T,Nrows,Ncols> qcomp(const int q){return SMatrix<T,Nrows,Ncols>(&A[indx(q,0,0)]);} // only used for printing

  
  vector<T> A; // the actual storage area
  
  bool ContainsNaN(string message="")
  {
    for(unsigned int i=0; i<A.size(); i++)
      { 
	if( A[i] != A[i] )	
	  {
	    cout << "Nan encountered in " << message << endl;
	    return true;
	  }
      }
    return false;
  }
  
  
  VecMat& operator=(VecMat& rhs)
  {
    T* ptr=rhs.start();
    if(this == &rhs) return *this; // self-assignment
    for(unsigned int i=0; i<A.size(); i++){A[i]= *(ptr+i);}
    return *this;
  }
  
  VecMat& operator+=(VecMat& rhs) 
  {
    T* ptr=rhs.start();
    for(unsigned int i=0; i<A.size(); i++){A[i]+= *(ptr+i);}
    return *this;
  }
  
  VecMat& operator-=(VecMat& rhs) 
  {
    T* ptr=rhs.start();
    for(unsigned int i=0; i<A.size(); i++){A[i]-= *(ptr+i);}
    return *this;
  }
  
  // sublattice matrix product for each q
  VecMat& operator*=(VecMat& rhs) 
  {
    int Nvecs=size()/(Nrows*Ncols);
    vector<T> temprow(Ncols); // vector for storing temporary result
    for(int q=0; q<Nvecs; q++)
      {
	  for(unsigned int i=0; i<Nrows; i++)
	    {
	      for(unsigned int j=0; j<Ncols; j++)
		{
		  T sum(0);
		  for(int k=0; k<Ncols; k++){ sum += (*this)(q,i,k)*rhs(q,k,j);}
		  temprow[j]=sum;
		}
	      for(unsigned int l=0; l<Ncols; l++){ (*this)(q,i,l)=temprow[l];} // fill in matrix
	    }
	}
      return *this;
    }

  VecMat& operator*=(const complex<realtype>& c)
    {
      for(unsigned int i=0; i<A.size(); i++){ A[i]*=c;} 
      return *this;
    }

  
  VecMat& operator*=(const realtype& c)
    {
      for(unsigned int i=0; i<A.size(); i++){ A[i]*=c;} 
      return *this;
    }
 

  SMatrix<T,Nrows,Ncols> Sumq()
    {
      SMatrix<T,Nrows,Ncols> res;
      for(int q=0; q< size()/(Nrows*Ncols); q++)
	{
	  T* ptr=A[q];
	  for(unsigned int j=0; j<Nrows*Ncols; j++){ res[j] += *(ptr+j);}
	} 
      return res;
    }

  T Sumq(const int m1,const int m2)
    {
      T res(0);
      for(int q=0; q< size()/(Ncols*Nrows); q++){ res += A[indx(q,m1,m2)];} 
      return res;
    }

private:
inline int indx(const int j,const int m1,const int m2){return (j*Nrows+m1)*Ncols+m2;}
inline int smindx(const int m1,const int m2){return m1*Ncols+m2;}
};


template<class T,int Nrows,int Ncols>
  VecMat<T,Nrows,Ncols> operator+(VecMat<T,Nrows,Ncols>& lhs,VecMat<T,Nrows,Ncols>& rhs)
{
  VecMat<T,Nrows,Ncols> res(lhs);
  res += rhs;
  return res;
}


template<class T,int Nrows,int Ncols>
VecMat<T,Nrows,Ncols> operator-(VecMat<T,Nrows,Ncols>& lhs,VecMat<T,Nrows,Ncols>& rhs)
{
  VecMat<T,Nrows,Ncols> res(lhs);
  res -= rhs;
  return res;
}

template<class T,int Nrows,int Ncols>
VecMat<T,Nrows,Ncols> operator*(VecMat<T,Nrows,Ncols>& lhs,VecMat<T,Nrows,Ncols>& rhs)
{
  VecMat<T,Nrows,Ncols> res(lhs);
  res *= rhs;
  return res;
}

template<class T,int Nrows,int Ncols>
VecMat<T,Nrows,Ncols> operator*(VecMat<T,Nrows,Ncols>& lhs,const complex<realtype>& c)
{
  VecMat<T,Nrows,Ncols> res(lhs);
  res *= c;
  return res;
}

template<class T,int Nrows,int Ncols>
VecMat<T,Nrows,Ncols> operator*(VecMat<T,Nrows,Ncols>& lhs,const realtype& c)
{
  VecMat<T,Nrows,Ncols> res(lhs);
  res *= c;
  return res;
}


template<class T,int Nrows,int Ncols>
T Tr(VecMat<T,Nrows,Ncols>& K)
{
  T sum(0);
  for(unsigned int q=0; q<K.Nvecs; q++)
    for(unsigned int i=0; i<Nrows; i++)
      sum += K(q,i,i);

  return sum;
}



template<class T,int Nrows,int Ncols>
SMatrix<T,Nrows,Ncols> Sumq(VecMat<T,Nrows,Ncols>& K)
{
  return K.Sumq();
}

template<class T,int Nrows,int Ncols>
T Sumq(VecMat<T,Nrows,Ncols>& K,const int m1,const int m2)
{
  return K.Sumq(m1,m2);
}




template<class T,int Nrows,int Ncols>
void Chomp(VecMat<T,Nrows,Ncols>& K)
{
  T* ptr=K.start();
  for(unsigned int i=0; i<K.size(); i++)
    {
      T& t(*(ptr+i));
      if(abs(real(t))<epsilon){ t.real(0.);} 
      if(abs(imag(t))<epsilon){ t.imag(0.);} 
    }
}

template<class T,int Nrows,int Ncols>
T FindMax(VecMat<T,Nrows,Ncols>& K)
{
  T maxentry;
  realtype maxval=0;

  T* ptr=K.start();
  for(unsigned int i=0; i<K.size(); i++)
    {
      T& thisentry(*(ptr+i));
      if(abs(thisentry)>maxval){maxval=abs(thisentry); maxentry=thisentry;}
    }
  return maxentry;
}

template<class T,int Nrows,int Ncols>
T FindMaxImag(VecMat<T,Nrows,Ncols>& K)
{
  T maxentry;
  realtype maxval=0;

  T* ptr=K.start();
  for(unsigned int i=0; i<K.size(); i++)
    {
      T& thisentry(*(ptr+i));
      if(abs(thisentry.imag()) > maxval){maxval=abs(thisentry.imag()); maxentry=thisentry;}
    }
  return maxentry;
}

template<class T,int Nrows,int Ncols>
void MakeReal(VecMat<T,Nrows,Ncols>& K)
{
  T* ptr=K.start();
  for(unsigned int i=0; i<K.size(); i++)
    {
      if(TRACE)
	{
	  realtype ipart=imag(*(ptr+i));
	  if( fabs(ipart) > sensitivity)
	    {
	      cout << "Warning, MakeReal sets entry " << ipart << " to 0" << endl;
	    }
	}
      (*(ptr+i)).imag(0.);
    } 
}

template<class T,int Nrows,int Ncols>
void ComplexConjugate(VecMat<T,Nrows,Ncols>& K)
{
  T* ptr=K.start();
  for(unsigned int i=0; i<K.size(); i++){ conj(*(ptr+i));}
}

template<class T,int Nrows,int Ncols>
void MakeRealSymmetric(VecMat<T,Nrows,Ncols>& M)
{
  for(unsigned int q=0; q<M.size()/(Nrows*Ncols); q++)
    {
      for(unsigned int s=0; s<Nrows; s++)
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



template<class T,int Nrows,int Ncols>
void MakeRealOnDiagonal(VecMat<T,Nrows,Ncols>& M)
{ 
  for(int q=0; q<M.size()/(Nrows*Ncols); q++)
    for(int s=0; s<Nrows; s++)
      M(q,s,s).imag(0.); // real diagonal
}      


template<class T,int Nrows,int Ncols>
void MakeAntiSymmetric(VecMat<T,Nrows,Ncols>& M)
{
  for(int q=0; q<M.size()/(Nrows*Ncols); q++)
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




template<class T,int Nrows,int Ncols>
void AddToDiagonal(VecMat<T,Nrows,Ncols>& M,T& value)
{
  for(int q=0; q<M.size()/(Nrows*Ncols); q++)
    for(int s=0; s<Nrows; s++)
      M(q,s,s)+=value;
}


template<class T,int Nrows,int Ncols>
void SubtractFromDiagonal(VecMat<T,Nrows,Ncols>& A,T value)
{
  T mvalue=-value;
  AddToDiagonal(A,mvalue);
}


template<class T,int Nrows,int Ncols>
void AddToDiagonal(VecMat<T,Nrows,Ncols>& M,vector<T>& value)
{
  for(int q=0; q<M.size()/(Nrows*Ncols); q++)
    for(int s=0; s<Nrows; s++)
      M(q,s,s)+=value[s];
}



template<class T,int Nrows,int Ncols>
void SubtractFromDiagonal(VecMat<T,Nrows,Ncols>& A,vector<T> value)
{
  for(unsigned int s=0; s<value.size(); s++){value[s]=-value[s];}
  AddToDiagonal(A,value);
}



#endif // vecmat_h


