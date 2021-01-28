#ifndef COUPLINGS_H
#define COUPLINGS_H

class Couplings
{
 public:
  Couplings(double* in_par,const int,const int);
 private:
  double* par;
  const int nc;
  const int nmat;
  void InitializeJ();
#ifdef PHONONS
  void Initializeg();
#endif
 public:
  VecMat<complex<realtype>> J;
#ifdef PHONONS
  VecMat<complex<realtype>> g;
#endif
};

Couplings::Couplings(double* in_par,const int in_nc,const int in_nmat):par(in_par),nc(in_nc),nmat(in_nmat),J(in_nc,in_nmat,in_nmat)
#ifdef PHONONS
,g(in_nc,in_nmat,in_nmat)
#endif
{
  logfile << "Initializing J" << endl;
  InitializeJ();
  logfile << "Done initializing J" << endl;
  for(int i=0; i<nc; i++)
    logfile << "i=" << i << " " << J[i] << endl;

#ifdef PHONONS
  logfile << "Initializing g" << endl;
  Initializeg();
  logfile << "Done initializing g" << endl;
  for(int i=0; i<nc; i++)
    logfile << "i=" << i << " " << g[i] << endl;
#endif
}




#endif // COUPLINGS_H
