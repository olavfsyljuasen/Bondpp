#ifndef COUPLINGS_H
#define COUPLINGS_H

class Couplings
{
 public:
  Couplings(double* in_par,const int);
 private:
  double* par;
  const int nc;
  void InitializeJ();
#ifdef PHONONS
  void Initializeg();
#endif
 public:
  VecMat<complextype ,NMAT,NMAT> J;
#ifdef PHONONS
  VecMat<complextype ,NMAT,NMAT> g;
#endif
};

Couplings::Couplings(double* in_par,const int in_nc):par(in_par),nc(in_nc),J(in_nc)
#ifdef PHONONS
,g(in_nc)
#endif
{
  if(TRACE){cout << "Initializing J" << endl;}
  InitializeJ();
  if(TRACE)
    {
      cout << "Done initializing J" << endl;
      for(int i=0; i<nc; i++)
	cout << "i=" << i << "\n" << J[i] << endl;
    }

#ifdef PHONONS
  if(TRACE){cout << "Initializing g" << endl;}
  Initializeg();
  if(TRACE)
    {
      cout << "Done initializing g" << endl;
      for(int i=0; i<nc; i++)
	cout << "i=" << i << "\n" << g[i] << endl;
    }
#endif
}




#endif // COUPLINGS_H
