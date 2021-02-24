#ifndef SYMMETRYROUTINES_H
#define SYMMETRYROUTINES_H


template<class T>
void MakeHermitian(VecMat<T>& M,bool warning=true)
{
  const int Nrows=M.Nrows;
  const int Ncols=M.Ncols;
  
  for(int q=0; q<M.Nvecs; q++)
    {
      for(int s=0; s<Nrows; s++)
	{
	  if(TRACE && warning)
	    {
	      if( fabs(imag(M(q,s,s))) > sensitivity){
		cout << "Warning: MakeHermitian sets an entry " << imag(M(q,s,s)) << " to 0" << endl;}
	    }
	  M(q,s,s).imag(0.); // real diagonal
	}
      
      for(int s1=0; s1<Nrows; s1++)
	for(int s2=s1+1; s2<Ncols; s2++)
	  {
	    if(TRACE && warning)
	      {
		complex<realtype> difference= M(q,s2,s1) - conj(M(q,s1,s2));
		if( abs(difference)  > sensitivity)
		  {
		    cout << "Warning: MakeHermitian makes two entries that are different by "
			 << difference << " equal" << endl;
		  }
	      }
	    M(q,s2,s1) = conj(M(q,s1,s2)); // off-diagonals are cc of each other.
	  }
    }
}

// make all diagonal spin entries equal 
template<class T>
void MakeEqualSpinDiagonalEntries(VecMat<T>& M,bool warning=true)
{
  for(int q=0; q<M.Nvecs; q++)
    {
      for(int l=0; l<NSUBL; l++)
	{
	  const int m0=mindx(0,l);
	  for(int s=0; s<NSPIN; s++)
	    {
	      const int m =mindx(s,l);

	      if(TRACE && warning)
	      {
		realtype diff= abs(M(q,m,m)-M(q,m0,m0));
		if( diff  > sensitivity)
		  {
		    cout << "Warning: MakeEqualSpinDiagonalEntries sets "
			 << M(q,m,m) << " equal to " << M(q,m0,m0) << " dev: " << diff << endl;
		  }
	      }
	      M(q,m,m)=M(q,m0,m0);
	    }
	}
    }
}

// erase all off-diagonal spin entries. 
template<class T>
void EraseSpinOffdiagonals(VecMat<T>& M,bool warning=true)
{
  for(int q=0; q<M.Nvecs; q++)
    {
      for(int l=0; l<NSUBL; l++)
	for(int s1=0; s1<NSPIN; s1++)
	  for(int s2=s1+1; s2<NSPIN; s2++)
	    {
	      int m1=mindx(s1,l);
	      int m2=mindx(s2,l);
	      
	      if(TRACE && warning)
		{
		  if( abs(M(q,m1,m2))  > sensitivity){
		    cout << "Warning: MakeSpinDiagonal sets an entry " << M(q,m1,m2) << " to 0" << endl;}
		  if( abs(M(q,m2,m1)) > sensitivity){
		    cout << "Warning: MakeSpinDiagonal sets an entry " << M(q,m2,m1) << " to 0" << endl;}
		}
	      M(q,m1,m2)=0.;
	      M(q,m2,m1)=0.;
	    }
    }
}




template<class T>
void MakeSpinSymmetric(VecMat<T>& M,bool warning=true)
{
  EraseSpinOffdiagonals(M,warning);
  MakeEqualSpinDiagonalEntries(M,warning);
}


template<class T>
bool IsSpinSymmetric(VecMat<T>& M)
{
  for(int q=0; q<M.Nvecs; q++)
    {
      for(int l=0; l<NSUBL; l++)
	for(int s1=0; s1<NSPIN; s1++)
	  {
	    int m0=mindx( 0,l);
	    int m1=mindx(s1,l);

	    realtype diff= abs( M(q,m1,m1)-M(q,m0,m0));
	    
	    if( diff > sensitivity)
	      {
		cout << "M(" << q << "," << m1 << "," << m1 << ")=" << M(q,m1,m1) << " deviates from "
		     << M(q,m0,m0) << " dev = " << diff << endl;   
		return false;
	      }
	 
	    for(int s2=s1+1; s2<NSPIN; s2++)
	      {
		int m2=mindx(s2,l);
		if( abs(M(q,m1,m2)) > sensitivity)
		  {
		    cout << "M(" << q << "," << m1 << "," << m2 << ")=" << M(q,m1,m2) << " is not zero " << endl;
		    return false;
		  }
		if( abs(M(q,m2,m1)) > sensitivity)
		  {
		    cout << "M(" << q << "," << m2 << "," << m1 << ")=" << M(q,m2,m1) << " is not zero " << endl;
		    return false;
		  }
	      }
	  }
    }
  return true;
}



// makes a matrix Hermitian in the block, but anti-Hermitian between the blocks
template<class T>
void MakeMixedHermitian(VecMat<T>& M,const int N1rows,const int N1cols)
{
  // N1rows and N1cols are the dimensions of the upper block matrix

  const int Nrows=M.Nrows;
  const int Ncols=M.Ncols;
  
  for(int q=0; q<M.Nvecs; q++)
    {
      for(int s=0; s<Nrows; s++)
	{
	  if(TRACE)
	    {
	      if( fabs(imag(M(q,s,s))) > sensitivity){
		cout << "Warning: MakeMixedHermitian sets the M(" << q << "," << s << "," << s <<").imag= " << imag(M(q,s,s)) << " to 0" << endl;}
	    }
	  M(q,s,s).imag(0.); // real diagonal
	}
      
      for(int s1=0; s1<Nrows; s1++)
	for(int s2=s1+1; s2<Ncols; s2++)
	  {
	    double sign=( ((s1>=N1rows && s2<N1cols) || (s1<N1rows && s2>=N1cols)) ? -1. : 1.);
	    if(TRACE)
	      {
		complex<realtype> difference= M(q,s2,s1) - sign*conj(M(q,s1,s2));
		if( abs(difference)  > sensitivity)
		  {
		    cout << "Warning: MakeMixedHermitian makes two entries that are different: "
			 << "q=" << q << " s1= " << s1 << " s2=" << s2 << " "
		         << M(q,s2,s1) << " and " << sign*conj(M(q,s1,s2)) << " equal. (difference="
			 << difference << ")" << endl;
		    cout << "M(" << q << "," << s1 << "," << s2 << ")=" << M(q,s1,s2)
			 << " M(" << q << "," << s2 << "," << s1 << ")=" << M(q,s2,s1) << endl;
		      
		  }
	      }
	    M(q,s2,s1) = sign*conj(M(q,s1,s2)); // off-diagonals are cc of each other.
	  }
    }
}


template<class T>
void MakeInversionTransposedSymmetric(VecMat<T>& M)
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



template<class T>
bool IsHermitian(VecMat<T>& M)
{
  //  if(TRACE) cout << "Checking if Hermitian: ";
  //  bool retval=true;
  const int Nrows=M.Nrows;
  const int Ncols=M.Ncols;

  for(int q=0; q<M.Nvecs; q++)
    {
      
      for(int s=0; s<Nrows; s++)
	if(abs(imag(M(q,s,s)))>sensitivity)
	  {
	    cout << "diagonal not real: M(" << q << "," << s << "," << s << ")=" << M(q,s,s) << endl;   
	    return false;
	  }
  
      for(int s1=0; s1<Nrows; s1++)
	for(int s2=s1+1; s2<Ncols; s2++)
	  {
	    complex<realtype> diff= M(q,s2,s1)-conj(M(q,s1,s2));
	    if( abs(real(diff))>sensitivity || abs(imag(diff))>sensitivity )
	      {
		cout << "M(" << q << "," << s1 << "," << s2 << ")^*=" << setprecision(20) << conj(M(q,s1,s2))
		     << " != " 
		     << "M(" << q << "," << s2 << "," << s1 << ")=" << setprecision(20) << M(q,s2,s1)
		     << " violation: " << conj(M(q,s1,s2))-M(q,s2,s1) << endl;
		return false;
	      }
	  }
    }
  //  if(TRACE) cout << "Hermitian it is!" << endl;
  return true;

}

template<class T>
bool IsMixedHermitian(VecMat<T>& M,const int N1rows,const int N1cols)
{
  //  if(TRACE) cout << "Checking if Hermitian: ";
  //  bool retval=true;
  const int Nrows=M.Nrows;
  const int Ncols=M.Ncols;

  for(int q=0; q<M.Nvecs; q++)
    {
      
      for(int s=0; s<Nrows; s++)
	if(abs(imag(M(q,s,s)))>sensitivity)
	  {
	    cout << "diagonal not real: M(" << q << "," << s << "," << s << ")=" << M(q,s,s) << endl;   
	    return false;
	  }
  
      for(int s1=0; s1<Nrows; s1++)
	for(int s2=s1+1; s2<Ncols; s2++)
	  {
	    double sign=( (s1>=N1rows && s2<N1cols) || (s1<N1rows && s2>=N1cols) ? -1. : 1.);

	    complex<realtype> diff= M(q,s2,s1)-sign*conj(M(q,s1,s2));
	    if( abs(real(diff))>sensitivity || abs(imag(diff))>sensitivity )
	      {
		cout << "M(" << q << "," << s1 << "," << s2 << ")^*=" << setprecision(20) << sign*conj(M(q,s1,s2))
		     << " != " 
		     << "M(" << q << "," << s2 << "," << s1 << ")=" << setprecision(20) << M(q,s2,s1)
		     << " violation: " << sign*conj(M(q,s1,s2))-M(q,s2,s1) << endl;
		return false;
	      }
	  }
    }
  //  if(TRACE) cout << "MixedHermitian it is!" << endl;
  return true;

}

template<class T>
bool IsInversionTransposedSymmetric(VecMat<T>& M)
{
  //  if(TRACE) cout << "Checking if InversionTransposedSymmetric: ";

  const int Nrows=M.Nrows;
  const int Ncols=M.Ncols;
  
  for(int q=0; q<M.Nvecs; q++)
    {
      const int invq=la.GetInversionIndx(q); // the negative q
      
      for(int m1=0; m1<Nrows; m1++)
	for(int m2=m1; m2<Ncols; m2++)
	  {
	    if( q==invq && abs(M(q,m1,m2).imag())>sensitivity)
	      {
		cout << "diagonal not real: M(" << q << "," << m1 << "," << m2 << ")=" << setprecision(20)
		     << M(q,m1,m2) << endl;  return false;
	      }
	    complex<realtype> diff= M(invq,m2,m1)-M(q,m1,m2);
	    if( abs(real(diff))>sensitivity || abs(imag(diff))>sensitivity )
	      {
		cout << "M(" << q << "," << m1 << "," << m2 << ")=" << setprecision(20) << M(q,m1,m2)
		     << " != " 
		     << "M(" << invq << "," << m2 << "," << m1 << ")=" << setprecision(20) << M(invq,m2,m1)
		     << " violation: " << M(q,m1,m2)-M(invq,m2,m1) << endl;
		     return false;
	      }
	  }
    }
  return true;
}

void SanityCheck(VecMat<complex<realtype>>& M,string message,bool checkHermitian=true)
{
  bool stop=false;
  if(M.ContainsNaN(message)){stop=true;}
  
  if(checkHermitian && !IsHermitian(M)){cout << "Matrix: " << message << " is not Hermitian."<< endl; stop=true;}
  if(!IsInversionTransposedSymmetric(M)){cout << "Matrix: " << message << " is not inv trans sym."; stop=true;}
  if(!IsSpinSymmetric(M)){cout << "WARNING Matrix: " << message << " is not spin symmetric" << endl; cout << M << endl; stop=false;}
  if(stop){ cout << "Sanity check failed! \n" << "Matrix " << message << "\n" << M << "\n SanityCheck FAILED" << endl; exit(1);}
  else
    {cout << "Sanity check passed for " << message << endl;}
}




















#endif // SYMMETRYROUTINES_H
