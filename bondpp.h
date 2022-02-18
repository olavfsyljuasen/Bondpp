#ifndef BOND_H
#define BOND_H

#include<vector>
#include<complex>
#include<iostream>
#include<fstream>
#include<algorithm>
#include<sstream>

#if defined PHONONS && !defined ELASTICONLY 
const auto NDMAT=NSUBL+NMODE;
#else
const auto NDMAT=NSUBL;
#endif
const auto NDMAT2=NDMAT*NDMAT;

  //#g++ -I/data/sylju/include -L/data/sylju/lib FFTWtest.C -lfftw3 -lm;
#include <fftw3.h>
#include "matrixroutines.h"

#ifdef LONGDOUBLE
typedef fftwl_plan FFTWPLAN;
typedef fftwl_complex FFTWCOMPLEX;
#define FFTWEXECUTE(x) fftwl_execute(x)
#define FFTWDESTROYPLAN(x) fftwl_destroy_plan(x)
#elif defined FLOAT
typedef fftwf_plan FFTWPLAN;
typedef fftwf_complex FFTWCOMPLEX;
#define FFTWEXECUTE(x) fftwf_execute(x)
#define FFTWDESTROYPLAN(x) fftwf_destroy_plan(x)
#elif defined QUAD
typedef fftwq_plan FFTWPLAN;
typedef fftwq_complex FFTWCOMPLEX;
#define FFTWEXECUTE(x) fftwq_execute(x)
#define FFTWDESTROYPLAN(x) fftwq_destroy_plan(x)
#else
typedef fftw_plan FFTWPLAN;
typedef fftw_complex FFTWCOMPLEX;
#define FFTWEXECUTE(x) fftw_execute(x)
#define FFTWDESTROYPLAN(x) fftw_destroy_plan(x)
#endif



class Driver
{
 public:
  Driver(Rule&);
  ~Driver()
    {
      FFTWDESTROYPLAN(A1q_to_A1r);
      FFTWDESTROYPLAN(A1r_to_A1q);
      FFTWDESTROYPLAN(A2q_to_A2r);
      FFTWDESTROYPLAN(A2r_to_A2q);
      FFTWDESTROYPLAN(Bq_to_Br);
      FFTWDESTROYPLAN(Br_to_Bq);
      FFTWDESTROYPLAN(F1pluss);
      FFTWDESTROYPLAN(F1minus);
      FFTWDESTROYPLAN(F2pluss);
      FFTWDESTROYPLAN(F2minus);
    }
  realtype CalculateT(int);
  void CalculateTs(vector<realtype>&);


#ifdef PHONONS
  NumberList CalculateEpsilonsOverT();
#endif

  realtype CalculateFreeEnergy(const realtype);
  vector<obstype> CalculateSpinOrderPars(realtype);
  vector<obstype> CalculateOrderPars(realtype,int,int);
  vector<obstype> CalculateAlphas(realtype);

  void FindMaxVals(SMatrix<int,NMAT,NMAT>&, SMatrix<realtype,NMAT,NMAT>&);

  //void Convolve(const bool);
  //  realtype CalculateSecondDerivative(const realtype T,const int k);

  void ConstructKinvq();
  void ComputeDq(const bool,const bool);  
  void ComputeSelfEnergy(const bool);
  
  void SetQsToZero(); // routine to set some q's to zero in self-energy
  void MakeRandomSigma();

  void MakeSymmetric(VecMat<complex<realtype>,NMAT,NMAT>&);
  
  void SolveSelfConsistentEquation(NumberList Delta); 

  void Solve(NumberList delta,NumberList thisepsilon,const bool pinfo)
  {
    Delta = delta;

#ifdef PHONONS
    if( USEPREVIOUSEPSILONS && converged && EpsilonInitialized)
      {
	// use the old epsilons if previous step converged, do not do anything
	logfile << "Use previous converged epsilons" << endl;
      }
    else
      {
	// set in new starting epsilon values
	epsilon=thisepsilon;
	EpsilonInitialized=true;
      }
#endif
    
    if(TRACE) cout << "Starting Solve with Delta= " << delta << " epsilon=" << epsilon << " Printinfo= " << pinfo << endl;
    logfile << "Starting Solver with Delta= " << delta << " epsilon=" << epsilon << " Printinfo= " << pinfo << endl;
    
    
    Printinfo=pinfo;
    //    Sigma = rule.GetInitialState(); // copy the Initial guess for Sigma
    SolveSelfConsistentEquation(Delta);
    if(TRACE) cout << "Done Solve" << endl;
  }
  
  
 private:
  Rule& rule;

  int dim;              // the number of dimensions
  vector<int> dims;     // the size of each dimension
  
  const int Nq; // number of q-space sites.
  const realtype invNq;
  const realtype invSqrtNq;

  bool converged; //  
  
  NumberList Delta; 
  
#ifdef PRESERVESYMMETRY
  int TransformationPeriod;
  vector<int> TransformationTable;
#endif
  
  bool Printinfo; // set this to print out correlation functions
  int lineid; 
  bool SigmaInitialized; // flag to indicate whether Sigma is initialized.  
  bool EpsilonInitialized; // flag to indicate whether Epsilonlist is initialized.  
  int MaxIterMultiplier; // multiplier for MaxIter, usually only used in the first step

  realtype mineigenvalue; // for storing the minimum SigmaE value
  realtype currT; // for storing the current value of the temperature 

  VecMat<complex<realtype>,NMAT,NMAT>& Jq;
  // The actual storage areas
  VecMat<complex<realtype>,NMAT,NMAT> A1;  // holds Kq,Kinvq
  VecMat<complex<realtype>,NMAT,NMAT> A2;  // holds Sigmaq,
  VecMat<complex<realtype>,NDMAT,NDMAT> B;  // holds Dq,Dinvq


  vector<complex<realtype>>  F1 ; // holds inplace intermediate Fourier-transform
  vector<complex<realtype>>  F2 ; // holds inplace intermediate Fourier-transform

#ifdef FFTS_INPLACE
  VecMat<complex<realtype>,NMAT,NMAT>& A1r; // holds Kinvr
  VecMat<complex<realtype>,NMAT,NMAT>& A2r; // holds Sigmar
  VecMat<complex<realtype>,NDMAT,NDMAT>& Br; // holds Dr
#else   
  VecMat<complex<realtype>,NMAT,NMAT> A1r;  // holds Kinvr
  VecMat<complex<realtype>,NMAT,NMAT> A2r;  // holds Sigmar
  VecMat<complex<realtype>,NDMAT,NDMAT> Br;  // holds Dr
#endif


  
  FFTWPLAN A1q_to_A1r;
  FFTWPLAN A1r_to_A1q;
  FFTWPLAN A2q_to_A2r;
  FFTWPLAN A2r_to_A2q;

  FFTWPLAN Bq_to_Br;
  FFTWPLAN Br_to_Bq;

  FFTWPLAN F1pluss;
  FFTWPLAN F1minus;

  FFTWPLAN F2pluss;
  FFTWPLAN F2minus;
 
  // the following references are used for readability of the code


  VecMat<complex<realtype>,NMAT,NMAT>& Kq;     // points to the A1 array
  VecMat<complex<realtype>,NMAT,NMAT>& Kinvq;  // points to the A1 array
  VecMat<complex<realtype>,NMAT,NMAT>& Kinvr;  // points to the A1 array

  VecMat<complex<realtype>,NMAT,NMAT>& Sigmar; // points to the A2 array
  VecMat<complex<realtype>,NMAT,NMAT>& Sigmaq; // points to the A2 array

  VecMat<complex<realtype>,NDMAT,NDMAT>& Dq;     // points to the B array
  VecMat<complex<realtype>,NDMAT,NDMAT>& Dinvq;  // points to the B array
  VecMat<complex<realtype>,NDMAT,NDMAT>& Dr;     // points to the B array

#ifdef PHONONS
  VecMat<complex<realtype>,NC,NMODE>& f; // holds the vertex information
  VecMat<complex<realtype>,NMAT,NMAT>& g; // points to rule
  
  vector<VecMat<complex<realtype>,NMAT,NMAT>* > gel;
#endif
  NumberList epsilon;  // amplitude of elastic modes


  
};


Driver::Driver(Rule& r):rule(r),dim(la.D()),dims(la.SiteqDims()),Nq(la.NqSites()),invNq(static_cast<realtype>(1.)/Nq),invSqrtNq(1./sqrt(Nq)),converged(false),Delta(NSUBL),Printinfo(false),lineid(0),SigmaInitialized(false),EpsilonInitialized(false),MaxIterMultiplier(1),Jq(r.Jq),  
  A1(Nq),A2(Nq),B(Nq),F1(Nq),F2(Nq),
#ifdef FFTS_INPLACE
  A1r(A1),A2r(A2),Br(B),
#else
  A1r(Nq),A2r(Nq),Br(Nq),
#endif
  Kq(A1),
  Kinvq(A1),Kinvr(A1r),Sigmar(A2r),Sigmaq(A2),
  Dq(B),Dinvq(B),Dr(Br)
#ifdef PHONONS
  ,f(r.Getf()),g(r.g)
  ,epsilon(NELASTIC)
#else
  ,epsilon(0)
#endif
{
  if(TRACE) cout << "Initializing solver " << endl;

  //Setting up fftw_plans, in-place ffts:
  /*
  fftw_plan fftw_plan_many_dft(int rank, const int *n, int howmany,
			       fftw_complex *in, const int *inembed,
			       int istride, int idist,
			       fftw_complex *out, const int *onembed,
			       int ostride, int odist,
			       int sign, unsigned flags);
  */
  //  unsigned flags=FFTW_MEASURE;
  unsigned flags=( Nq > 1000 ? FFTW_PATIENT: FFTW_ESTIMATE);
  //unsigned flags=FFTW_PATIENT;

  logfile << "Making FFT plans" << endl;    
  FFTWCOMPLEX* A1_ptr  =reinterpret_cast<FFTWCOMPLEX*>(A1.start());
  FFTWCOMPLEX* A1r_ptr =reinterpret_cast<FFTWCOMPLEX*>(A1r.start());
  FFTWCOMPLEX* A2_ptr  =reinterpret_cast<FFTWCOMPLEX*>(A2.start());
  FFTWCOMPLEX* A2r_ptr =reinterpret_cast<FFTWCOMPLEX*>(A2r.start());
  FFTWCOMPLEX* B_ptr  =reinterpret_cast<FFTWCOMPLEX*>(B.start());
  FFTWCOMPLEX* Br_ptr =reinterpret_cast<FFTWCOMPLEX*>(Br.start());

  FFTWCOMPLEX* F1_ptr=reinterpret_cast<FFTWCOMPLEX*>(&F1[0]);
  FFTWCOMPLEX* F2_ptr=reinterpret_cast<FFTWCOMPLEX*>(&F2[0]);


#ifdef LONGDOUBLE
  A1q_to_A1r = fftwl_plan_many_dft(dim,&dims[0],NMAT2 ,A1_ptr ,0,NMAT2 ,1,A1r_ptr,0,NMAT2 ,1,+1,flags);
  A1r_to_A1q = fftwl_plan_many_dft(dim,&dims[0],NMAT2 ,A1r_ptr,0,NMAT2 ,1,A1_ptr ,0,NMAT2 ,1,-1,flags);
  A2q_to_A2r = fftwl_plan_many_dft(dim,&dims[0],NMAT2 ,A2_ptr ,0,NMAT2 ,1,A2r_ptr,0,NMAT2 ,1,+1,flags);
  A2r_to_A2q = fftwl_plan_many_dft(dim,&dims[0],NMAT2 ,A2r_ptr,0,NMAT2 ,1,A2_ptr ,0,NMAT2 ,1,-1,flags);
  Bq_to_Br   = fftwl_plan_many_dft(dim,&dims[0],NDMAT2,B_ptr  ,0,NDMAT2,1,Br_ptr ,0,NDMAT2,1,+1,flags);  
  Br_to_Bq   = fftwl_plan_many_dft(dim,&dims[0],NDMAT2,Br_ptr ,0,NDMAT2,1,B_ptr  ,0,NDMAT2,1,-1,flags);
 
  F1pluss    = fftwl_plan_many_dft(dim,&dims[0],1,F1_ptr,0,1,1,F1_ptr,0,1,1,+1,flags);
  F1minus    = fftwl_plan_many_dft(dim,&dims[0],1,F1_ptr,0,1,1,F1_ptr,0,1,1,-1,flags);  
  F2pluss    = fftwl_plan_many_dft(dim,&dims[0],1,F2_ptr,0,1,1,F2_ptr,0,1,1,+1,flags);
  F2minus    = fftwl_plan_many_dft(dim,&dims[0],1,F2_ptr,0,1,1,F2_ptr,0,1,1,-1,flags);  

#elif defined FLOAT
  A1q_to_A1r = fftwf_plan_many_dft(dim,&dims[0],NMAT2 ,A1_ptr ,0,NMAT2 ,1,A1r_ptr,0,NMAT2 ,1,+1,flags);
  A1r_to_A1q = fftwf_plan_many_dft(dim,&dims[0],NMAT2 ,A1r_ptr,0,NMAT2 ,1,A1_ptr ,0,NMAT2 ,1,-1,flags);
  A2q_to_A2r = fftwf_plan_many_dft(dim,&dims[0],NMAT2 ,A2_ptr ,0,NMAT2 ,1,A2r_ptr,0,NMAT2 ,1,+1,flags);
  A2r_to_A2q = fftwf_plan_many_dft(dim,&dims[0],NMAT2 ,A2r_ptr,0,NMAT2 ,1,A2_ptr ,0,NMAT2 ,1,-1,flags);
  Bq_to_Br   = fftwf_plan_many_dft(dim,&dims[0],NDMAT2,B_ptr  ,0,NDMAT2,1,Br_ptr ,0,NDMAT2,1,+1,flags);  
  Br_to_Bq   = fftwf_plan_many_dft(dim,&dims[0],NDMAT2,Br_ptr ,0,NDMAT2,1,B_ptr  ,0,NDMAT2,1,-1,flags);

  F1pluss    = fftwf_plan_many_dft(dim,&dims[0],1,F1_ptr,0,1,1,F1_ptr,0,1,1,+1,flags);
  F1minus    = fftwf_plan_many_dft(dim,&dims[0],1,F1_ptr,0,1,1,F1_ptr,0,1,1,-1,flags);  
  F2pluss    = fftwf_plan_many_dft(dim,&dims[0],1,F2_ptr,0,1,1,F2_ptr,0,1,1,+1,flags);
  F2minus    = fftwf_plan_many_dft(dim,&dims[0],1,F2_ptr,0,1,1,F2_ptr,0,1,1,-1,flags);  
#else
  A1q_to_A1r = fftw_plan_many_dft(dim,&dims[0],NMAT2 ,A1_ptr ,0,NMAT2 ,1,A1r_ptr,0,NMAT2 ,1,+1,flags);
  A1r_to_A1q = fftw_plan_many_dft(dim,&dims[0],NMAT2 ,A1r_ptr,0,NMAT2 ,1,A1_ptr ,0,NMAT2 ,1,-1,flags);
  A2q_to_A2r = fftw_plan_many_dft(dim,&dims[0],NMAT2 ,A2_ptr ,0,NMAT2 ,1,A2r_ptr,0,NMAT2 ,1,+1,flags);
  A2r_to_A2q = fftw_plan_many_dft(dim,&dims[0],NMAT2 ,A2r_ptr,0,NMAT2 ,1,A2_ptr ,0,NMAT2 ,1,-1,flags);
  Bq_to_Br   = fftw_plan_many_dft(dim,&dims[0],NDMAT2,B_ptr  ,0,NDMAT2,1,Br_ptr ,0,NDMAT2,1,+1,flags);  
  Br_to_Bq   = fftw_plan_many_dft(dim,&dims[0],NDMAT2,Br_ptr ,0,NDMAT2,1,B_ptr  ,0,NDMAT2,1,-1,flags);

  F1pluss    = fftw_plan_many_dft(dim,&dims[0],1,F1_ptr,0,1,1,F1_ptr,0,1,1,+1,flags);
  F1minus    = fftw_plan_many_dft(dim,&dims[0],1,F1_ptr,0,1,1,F1_ptr,0,1,1,-1,flags);  
  F2pluss    = fftw_plan_many_dft(dim,&dims[0],1,F2_ptr,0,1,1,F2_ptr,0,1,1,+1,flags);
  F2minus    = fftw_plan_many_dft(dim,&dims[0],1,F2_ptr,0,1,1,F2_ptr,0,1,1,-1,flags);  
#endif

  logfile << "Done making FFT plans" << endl;


#ifdef PRESERVESYMMETRY
  TransformationPeriod=lattice.TransformationPeriod;
  TransformationTable=lattice.TransformationTable;
#endif

  
  if(TRACE) cout << "Done initializing solver " << endl;
};


realtype Driver::CalculateT(int sl)
{
  if(TRACE) cout << "Starting CalculateT for sublattice: " << sl << endl;
  realtype sumalpha=0.;
  vector<realtype> alpha(NSPIN);
  for(int s=0; s<NSPIN; s++) 
    {
      int m=mindx(s,sl); // make the composite index.
      alpha[s] = NFAKESPINTRACE*real(Sumq(Kinvq,m,m))/(2.*Nq);
      sumalpha += alpha[s];
    }

  if(TRACE)
    {
      cout << "T sums: "; 
      for(int s=0; s<NSPIN; s++){ cout << "alpha[" << s << "]=" << alpha[s] << " ";}
      cout << endl;
    }
  
  realtype T = 1./sumalpha;

  if(TRACE) cout << "Done CalculateT "  << T << endl;
  return T;
}


void Driver::CalculateTs(vector<realtype>& Ts)
{
  if(TRACE) cout << "Starting CalculateTs " << endl;
  for(int s=0; s<NSUBL; s++)
    {
      Ts[s]=CalculateT(s);
    }
}


#ifdef PHONONS
NumberList Driver::CalculateEpsilonsOverT()
{
  if(TRACE) cout << "Starting CalculateEpsilonsOverT" << endl;

  NumberList epsoverT(NELASTIC);
  for(int i=0; i<NELASTIC; i++)
    {
      realtype mui= rule.elasticeigenvalues[i];
      realtype sum=0.;
      if(mui != 0.)
	{
	  for(int qi=0; qi<Nq; qi++)
	    {
	      SMatrix<complex<realtype>,NMAT,NMAT> tmp;
	      tmp=Kinvq[qi];
	      tmp *= (*rule.gelptrs[i])[qi];
	      sum += NFAKESPINTRACE*real(tr(tmp));
	    }
	  sum *= -1./(2.*mui*Nq*la.UnitCellVolume());
	}
      epsoverT[i]= sum;
    }
 
  if(TRACE) cout << "Done CalculateEpsilonsOverT "  << endl;
  return epsoverT;
}
#endif


void Driver::FindMaxVals(SMatrix<int,NMAT,NMAT>& maxqs,SMatrix<realtype,NMAT,NMAT>& maxvals)
{
  if(TRACE) cout << "Starting FindMaxVals" << endl;

  for(int m1=0; m1<NMAT; m1++)
    for(int m2=0; m2<NMAT; m2++)
      {
	int maxq(-1);
	realtype maxv(0.);
	for(int qi=0; qi<Nq; qi++)
	  {
	    complex<realtype> val=(m1==m2 ? Kinvq(qi,m1,m2): 0.5*(Kinvq(qi,m1,m2)+Kinvq(qi,m2,m1)));
	    if( val.real() > maxv){ maxv=val.real(); maxq=qi;}
	  }
	maxvals(m1,m2) = maxv;
	maxqs(m1,m2)   = maxq;
      }

  if(TRACE) cout << "Done FindMaxVals" << endl;
}



vector<obstype> Driver::CalculateSpinOrderPars(realtype T)
{
  if(TRACE) cout << "Starting CalculateSpinOrderPars" << endl;

  vector<obstype> opars(NSPINOBSERVABLES);
  for(int j=0; j<NSPINOBSERVABLES; j++)
    {
      obstype sum=0.;
      KernelFunction* f=spinobservables[j];

      for(int qi=0; qi<Nq; qi++)
	{
	  for(int m1=0; m1<NMAT; m1++)
	    for(int m2=0; m2<NMAT; m2++)
	      sum+= (*f)(qi,m1,m2)*Kinvq(qi,m1,m2);
	}
      opars[j]=0.5*T*invNq*NFAKESPINTRACE*sum;
    }

  if(TRACE) cout << "Done CalculateSpinOrderPars" << endl;
  return opars;
}



vector<obstype> Driver::CalculateOrderPars(realtype T,int m1,int m2)
{
  if(TRACE) cout << "Starting CalculateOrderPars" << endl;

  vector<obstype> opars(NOBSERVABLES);
  for(int j=0; j<NOBSERVABLES; j++)
    {
      obstype sum=0.;
      vector<obstype>& f=rule.GetIrrep(j);

      for(int i=0; i<Nq; i++)
	{
	  sum+= f[i]*Kinvq(i,m1,m2);
	}
      opars[j]=0.5*T*invNq*sum;
    }

  if(TRACE) cout << "Done CalculateOrderPars" << endl;
  return opars;
}

#ifdef PRESERVESYMMETRY
void Driver::MakeSymmetric(VecMat<complex<realtype>,NMAT,NMAT>& m)
{
  if(TRACE) cout << "Starting MakeSymmetric" << endl;
  if( NSUBL !=1){ cout << "MakeSymmetric works only for NSUBL=1" << endl; exit(1);}
  
  complex<realtype>* mstart=m.start(); // This only works for NSUBL=1 
  
  //  cout << "m=" << m << endl;

  VecMat<complex<realtype>> temp(m);
  complex<realtype>* tempstart=temp.start(); // This only works for NSUBL=1 
  
  //  cout << "temp=" << temp << endl;

  int p=1; // p=0 is the identity transformation
  vector<int> ThisT(TransformationTable);
  
  //cout << "Period: " << TransformationPeriod << endl;
  //cout << "Nq=" << Nq << endl;


  while( p < TransformationPeriod)
    {
      for(int i=0; i<Nq; i++)
	{
	  int q=ThisT[i];
	  tempstart[i]+=mstart[q];
	  ThisT[i] = TransformationTable[q]; // update ThisT for the next period
	}
      p++;
    }
  
  // Average over all transformations and transfer the result to the input array
  realtype invperiod=1./TransformationPeriod;
  for(int i=0; i<Nq; i++){ mstart[i]=tempstart[i]*invperiod;}  

  if(TRACE) cout << "End of MakeSymmetric" << endl;

}
#endif



// Free energy per unit volume
// we use the convention \nu=\nup
realtype Driver::CalculateFreeEnergy(realtype T)
{
  if(TRACE) cout << "Starting CalculateFreeEnergy " << endl;

  const int Ns=NSPIN*NFAKESPINTRACE; 
  
  realtype f=0;
 
  if(TRACE) cout << "--- T  = " << T << endl;
  // constants:
  //realtype betaf_constants= -( 0.5*NSUBL*log(2.*Nq) + 0.5*NSUBL*(NS-1)*log(PI));
  realtype betaf_constants = -0.5*NSUBL*log(Nq) -0.5*NSUBL*(Ns-2)*log(PI) - 0.5*NSUBL*log(TWOPI);
    
    
  if(TRACE) cout << "betaf_constants  = " << betaf_constants << endl;
  
  f += T*betaf_constants;

  //temperature factors
  realtype betaf_Tdep = 0.5*NSUBL*(Ns-2)*log(1./T); 
  if(TRACE) cout << "betaf_Tdep  = " << betaf_Tdep << endl;
  
  f += T*betaf_Tdep;



#if defined PHONONS 
  //elastic modes
  realtype f_elastic = 0;
  for(int i=0; i<NELASTIC; i++){ f_elastic += 0.5*epsilon[i]*epsilon[i]*rule.elasticeigenvalues[i];}


#if !defined OMITELASTICCONT
  if(TRACE) cout << "f_elastic  = " << f_elastic << endl;
  f += f_elastic;
#endif

  //#if !defined(ELASTICONLY)  || !defined(OMITPHONONCONT)
#ifndef OMITPHONONCONT

  // the following term is subtracted from D, but added here, so as to count phonon
  // contrubutions to the free energy where it is appropriate. Remember,
  // the Dinv in this code contains a term beta on the phonon diagonal part, also
  // when the coupling to phonons g=0. Thus in order to count this as a contribution to
  // the phonons part of the free energy, we subtract it from D and add it here.    
  realtype betaf_phononTdep = -0.5*NMODE*log(T); 
  if(TRACE) cout << "betaf_phononTdep  = " << betaf_Tdep << endl;
  
  f += T*betaf_phononTdep;


  // constants:
  realtype betaf_phononconst = -0.5*2.*NMODE*log(TWOPI);
  if(TRACE) cout << "betaf_phononconst  = " << betaf_phononconst << endl;

  f += T*betaf_phononconst;

  //the phonon spectrum
  realtype betaf_phonons = 2.*rule.GetSumLogOmegaoverV(); // 2 both X^2 and P^2

  if(TRACE) cout << "betaf_phonons  = " << betaf_phonons << endl;
  
  f += T*betaf_phonons;

#endif
  
#endif  
  
  //Must correct the Delta values for the subtraction of the minimum from SigmaE
  for(int i=0; i<NSUBL; i++) f += -(Delta[i]-mineigenvalue);

  if(TRACE) cout << "betaf_delta      = " << -(Delta[0]-mineigenvalue)/T << " (Delta= " << Delta[0] << " mineig: " << mineigenvalue << ")" << endl;
  
  realtype betaf_logKinvq   = -0.5*invNq*NFAKESPINTRACE*SumLogDet(Kinvq);

  if(TRACE) cout << "betaf_logKinvq   =  " << betaf_logKinvq << endl;
  f += T*betaf_logKinvq;
  
  ComputeDq(false,true); // excludeqzero=false, preserveinput=true not to jeopardize Keff
  
  realtype betaf_logD       = -0.5*invNq*SumLogDet(Dq) + 0.5*NMODE*log(T);
  if(TRACE) cout << "betaf_logD       =  " << betaf_logD << endl;
  f += T*betaf_logD;

  ComputeSelfEnergy(true); // ensure that Kinvq is the same.

  realtype betaf_KinvqSigma = -0.5*invNq*NFAKESPINTRACE*SumTr(Kinvq,Sigmaq);
  if(TRACE) cout << "betaf_KinvqSigma = " << betaf_KinvqSigma << endl;
  f += T*betaf_KinvqSigma;

  // the correction to the saddle-point are already taken into account
  return f;
}


/* The old free energy routine
// Free energy per unit volume
// we use the convention \nu=\nup
realtype Driver::CalculateFreeEnergy(realtype T)
{
  if(TRACE) cout << "Starting CalculateFreeEnergy " << endl;

  realtype f=0;
 
  if(TRACE) cout << "--- T  = " << T << endl;
  // constants:
  realtype betaf_constants= -( 0.5*NSUBL*log(2.*Nq) + 0.5*NSUBL*(NS-1)*log(PI)); 
  if(TRACE) cout << "betaf_constants  = " << betaf_constants << endl;
  
  f += T*betaf_constants;
  
  //Must correct the Delta values for the subtraction of the minimum from SigmaE
#ifdef SOFTCONSTRAINT
  for(int i=0; i<NSUBL; i++) f += -(Delta[i]-mineigenvalue)*(Delta[i]-mineigenvalue)/(4.*g*T);
#endif
  for(int i=0; i<NSUBL; i++) f += -(Delta[i]-mineigenvalue);

  if(TRACE) cout << "betaf_delta      = " << -(Delta[0]-mineigenvalue)/T << " (Delta= " << Delta[0] << " mineig: " << mineigenvalue << ")" << endl;
  
  // NSUBL*log(T) is needed because T is not part of Kinvq in the program
  realtype betaf_logKinvq   = -0.5*NS*(invNq*SumLogDet(Kinvq) + NSUBL*log(T));
  //  cout << "invNqq*SumLogDet= " << invNq*SumLogDet(Kinvq) << " log(T)= " << log(T) << endl;
  if(TRACE) cout << "betaf_logKinvq   =  " << betaf_logKinvq << endl;					  f += T*betaf_logKinvq;
  
  ComputeDq(false,true); // excludeqzero=false, preserveinput=true not to jeopardize Keff
  
  // NSUBL*2.*log(T) needed because T is not part of Dq in the program
  realtype betaf_logD       = -0.5*invNq*( SumLogDet(Dq) - Nq*NSUBL*2.*log(T) );
  if(TRACE) cout << "betaf_logD       =  " << betaf_logD << endl;
  f += T*betaf_logD;

  ComputeSelfEnergy(true); // ensure that Kinvq is the same.

  realtype betaf_KinvqSigma = -0.5*NS*invNq*SumTr(Kinvq,Sigmaq);
  if(TRACE) cout << "betaf_KinvqSigma = " << betaf_KinvqSigma << endl;
  f += T*betaf_KinvqSigma;

  // the correction to the saddle-point are already taken into account
  return f;
}
*/

// ComputeDq() computes the renormalized constraint propagator
//
// Input: Kinvq (stored in A1).
// Output: Dq (stored in B)
// 
// The routine uses in-place FFTs so info contained in A and B are modified
// On output:
// A1 = Kinvr, unless preserveinput=true, then A1=Kinvq
// B = Dq
//
// FFTW omits volume prefactors in fourier transforms, we adopt the convention that the Fourier transform
// is without prefactors in going from q->r, and the prefactor 1/Nq is inserted on going from r->q.
// This means that using FFTW for transforming q->r->q, one should divide the result by Nq. 
void Driver::ComputeDq(const bool excludeqzero=true, const bool preserveinput=false)
{
#ifdef PHONONS
  if(NSUBL !=1){ cout << "Error:ComputeDq is not implemented for NSUBL !=1, exiting." << endl; exit(1);}
#endif

  if(TRACE) cout << "Starting ComputeDq, excludeqzero=" << excludeqzero << ",preserveinput=" << preserveinput << endl;

  if(TRACE){SanityCheck(Kinvq,"Kinvq, Initializing ComputeDq");}


  FFTWEXECUTE(A1q_to_A1r); // Kinvq->Kinvr, stored in Ar, after: Ar=Kinvr*Sqrt(Nq)
  Kinvr *= invSqrtNq; // Ar=Kinvr

  
#ifdef FORCEINVERSIONSYMMETRY
  MakeReal(Kinvr);  // is this needed here?
  MakeInversionTransposedSymmetric(Kinvr); 
#endif
  
  // Initialize Dinvq:
  Dinvq.SetToZero();

  // first compute the constraint block

  complex<realtype> tmp(0);
  
  for(int m1=0; m1<NSUBL; m1++)
    for(int m2=m1; m2<NSUBL; m2++)
      {
	for(int r=0; r<Nq; r++)
	  {
	    complex<realtype> tt(0.);
	    for(int s1=0; s1<NSPIN; s1++)
	      for(int s2=0; s2<NSPIN; s2++)
		{
		  int alpha=mindx(s1,m2); // alpha_s = m2
		  int delta=mindx(s2,m1); // delta_s = m1
		  
		  tmp=Kinvr(r,alpha,delta);
		  tt += tmp*tmp;   
		}
	    F1[r]=0.5*NFAKESPINTRACE*tt;
	  }

	FFTWEXECUTE(F1pluss); // 
	
	for(int qi=0; qi<Nq; qi++)
	  {
	    Dinvq(qi,m1,m2)=F1[qi];
	    if(m2 != m1) Dinvq(qi,m2,m1)=conj(F1[qi]); 
	  }
      }

  if(TRACE2)
    {
      cout << "after constraint block: " << endl;
      cout << Dinvq << endl;
    }

#ifdef PHONONS
  // the offdiagonal phonon-constraint part
#ifdef CPOSITIVE
  realtype multiplier=2.;
#else
  realtype multiplier=1.;
#endif  
    
  for(int m1=0; m1<NSUBL; m1++) // m1 must be 0 here because phonons NSUBL=1
    for(int c=0; c<NC; c++)
      {
	for(int r=0; r<Nq; r++)
	  {
	    SMatrix<complex<realtype>,NMAT,NMAT> my;
	    my  = g[c];
	    my *= Kinvr[r];
	    int mrpc = la.rAdd(la.GetInversionIndx(r),clist[c]); // index of -r+c	    
	    my *= Kinvr[mrpc];
	    
	    F1[r] = NFAKESPINTRACE*tr(my);  
	  }
	
	// do fourier-transform of theta
	FFTWEXECUTE(F1pluss); 
	
	for(int qi=0; qi<Nq; qi++)
	  {
	    for(int m2=NSUBL; m2<NDMAT; m2++)
	      {
		int n2=m2-NSUBL; 
		
		//	cout << "q=" << qi << " f=" << f(qi,c,n2) << endl;
		Dinvq(qi,m1,m2)+= multiplier*0.5*invSqrtNq*F1[qi]*f(qi,c,n2)*conj(expi(la.qr(qi,clist[c])));
	      }		   
	  }
      }
  
  // Set the constraint-phonon part.
  for(int qi=0; qi<Nq; qi++)
    for(int m1=0; m1<NSUBL; m1++)
      for(int m2=NSUBL; m2<NDMAT; m2++)
	{
	  Dinvq(la.GetInversionIndx(qi),m2,m1)=Dinvq(qi,m1,m2);
	}

  if(TRACE2)
    {
      cout << "after constraint-phonon block: " << endl;
      cout << Dinvq << endl;
    }


  if(TRACE2)
    {
      cout << "Dinvq[4]: " << Dinvq.qcomp(4) << endl;
    }

  // finally the phonon-phonon part
#ifdef CPOSITIVE
  /*
    for(int c2=0; c2<NC; c2++)
      for(int c4=0; c4<NC; c4++)
	{	    
	  if(TRACE2)
	    cout << "EXPERIMENTAL code: c2=" << c2 << " c4=" << c4 << endl;
	  
	  for(int r=0; r<Nq; r++)
	    {
	      SMatrix<complex<realtype>> my(NMAT,NMAT);
	      
	      int mrpc4=la.rAdd(r,clist[c4]);
	      
	      my  = g[c4];
	      my *= Kinvr[mrpc4];
	      my *= g[c2];
	      
	      int mrpc2=la.rAdd(la.GetInversionIndx(r),clist[c2]);
	      
	      my *= Kinvr[mrpc2];
	      
	      F1[r]=invNq*NFAKESPINTRACE*tr(my);
	    } 
	  FFTWEXECUTE(F1pluss); 
	  
	  
	  for(int m1=NSUBL; m1<NDMAT; m1++) 
	    for(int m2=m1; m2<NDMAT; m2++) 
	      {
		
		int n1=m1-NSUBL; 
		int n2=m2-NSUBL; 
		
		for(int qi=0; qi<Nq; qi++)
		  {
		    const complex<realtype> myF1 =( qi> la.qlimit ? F1[la.GetInversionIndx(qi)].conj():F1[qi]);
		    
		    //		    Dinvq(qi,m1,m2)+=-F1[qi]*conj(f(qi,c2,n1))*f(qi,c4,n2);
		    Dinvq(qi,m1,m2)+=-myF1*conj(f(qi,c2,n1))*f(qi,c4,n2);
		  }
		
		if(TRACE2)
		  {
		    if(m1==1 && m2==1)
		      {
			cout << "F1[4]=" << F1[4] 
			     << " f(4,"<< c2 << "," << n1 << ")^*=" << conj(f(4,c2,n1)) 
			     << " f(4," << c4 << "," << n2 << ")=" << f(4,c4,n2) 
			     << " expi(-q,c4)= " << conj(expi(la.qr(4,clist[c4]))) << endl;
			cout << "Dinvq(4,1,1) = " << Dinvq(4,m1,m2) << endl;
		      }
		  }
		
		
	      }
	  
	  if(TRACE2)
	    {
	      cout << "B: Dinvq[4]: " << Dinvq.qcomp(4) << endl;
	    }
	  
	  	  
	  myivec=clist[c2]-clist[c4];
	  
	  for(int r=0; r<Nq; r++)
	    {
	      SMatrix<complex<realtype>> my(NMAT,NMAT);
	      
	      my  = g[c4];
	      my.Transpose();
	      my *= Kinvr[r];
	      my *= g[c2];
	      
	      int mrc2mc4=la.rAdd(la.GetInversionIndx(r),myivec);
	      
	      my *= Kinvr[mrc2mc4];
	      
	      F1[r]=NFAKESPINTRACE*tr(my);
	    } 
	  FFTWEXECUTE(F1pluss); 
	  
	  
	  for(int m1=NSUBL; m1<NDMAT; m1++) 
	    for(int m2=m1; m2<NDMAT; m2++) 
	      {
		
		int n1=m1-NSUBL; 
		int n2=m2-NSUBL; 
		
		for(int qi=0; qi<Nq; qi++)
		  {
		    Dinvq(qi,m1,m2)+=-invNq*F1[qi]*conj(f(qi,c2,n1))*f(qi,c4,n2);
		  }

		if(TRACE2)
		  {
		    if(m1==1 && m2==1)
		      {
			cout << "F1[4]=" << F1[4] 
			     << " f(4,"<< c2 << "," << n1 << ")^*=" << conj(f(4,c2,n1)) 
			     << " f(4," << c4 << "," << n2 << ")=" << f(4,c4,n2) 
			     << " new cont=" << -invNq*F1[4]*conj(f(4,c2,n1))*f(4,c4,n2) << endl; 
			cout << "Dinvq(4,1,1) = " << Dinvq(4,m1,m2) << endl;
		      }
		  }

	      }


	  if(TRACE2)
	    {
	      //	      cout << "Dinvq=" << Dinvq << endl;
	      cout << "C: Dinvq[4]=" << Dinvq.qcomp(4) << endl;
	    }

  */
    for(int c2=0; c2<NC; c2++)
      for(int c4=0; c4<NC; c4++)
	{	    
	  if(TRACE2)
	    cout << "c2=" << c2 << " c4=" << c4 << endl;

	  	  	  
	  Triplet myivec=clist[c2]+clist[c4];
	  
	  for(int r=0; r<Nq; r++)
	    {
	      SMatrix<complex<realtype>,NMAT,NMAT> my;
	      
	      my  = g[c4];
	      my *= Kinvr[r];
	      my *= g[c2];
	      
	      int mrpc2c4=la.rAdd(la.GetInversionIndx(r),myivec);
	      
	      my *= Kinvr[mrpc2c4];
	      	      
	      //    F1[r]=NFAKESPINTRACE*tr(my);
	      F1[r]=invNq*NFAKESPINTRACE*tr(my);
	    } 
	  FFTWEXECUTE(F1pluss); 
	  
	  
	  for(int m1=NSUBL; m1<NDMAT; m1++) 
	    for(int m2=m1; m2<NDMAT; m2++) 
	      {
		
		int n1=m1-NSUBL; 
		int n2=m2-NSUBL; 
		
		for(int qi=0; qi<Nq; qi++)
		  {
		    //		    Dinvq(qi,m1,m2)+=-invNq*F1[qi]*conj(f(qi,c2,n1))*f(qi,c4,n2)*conj(expi(la.qr(qi,clist[c4])));
		    Dinvq(qi,m1,m2)+=-F1[qi]*conj(f(qi,c2,n1))*f(qi,c4,n2)*conj(expi(la.qr(qi,clist[c4])));
		  }

		if(TRACE2)
		  {
		    if(m1==1 && m2==1)
		      {
			cout << "F1[4]=" << F1[4] 
			     << " f(4,"<< c2 << "," << n1 << ")^*=" << conj(f(4,c2,n1)) 
			     << " f(4," << c4 << "," << n2 << ")=" << f(4,c4,n2) 
			     << " expi(-q,c4)= " << conj(expi(la.qr(4,clist[c4]))) << endl;
			cout << "Dinvq(4,1,1) = " << Dinvq(4,m1,m2) << endl;
		      }
		  }
                       

	      }

	  if(TRACE2)
	    {
	      cout << "B: Dinvq[4]: " << Dinvq.qcomp(4) << endl;
	    }


	  
	  myivec=clist[c2]-clist[c4];
	  
	  for(int r=0; r<Nq; r++)
	    {
	      SMatrix<complex<realtype>,NMAT,NMAT> my;
	      
	      my  = g[c4];
	      my.Transpose();
	      my *= Kinvr[r];
	      my *= g[c2];
	      
	      int mrc2mc4=la.rAdd(la.GetInversionIndx(r),myivec);
	      
	      my *= Kinvr[mrc2mc4];
	      
	      F1[r]=NFAKESPINTRACE*tr(my);
	    } 
	  FFTWEXECUTE(F1pluss); 
	  
	  
	  for(int m1=NSUBL; m1<NDMAT; m1++) 
	    for(int m2=m1; m2<NDMAT; m2++) 
	      {
		
		int n1=m1-NSUBL; 
		int n2=m2-NSUBL; 
		
		for(int qi=0; qi<Nq; qi++)
		  {
		    Dinvq(qi,m1,m2)+=-invNq*F1[qi]*conj(f(qi,c2,n1))*f(qi,c4,n2);
		  }

		if(TRACE2)
		  {
		    if(m1==1 && m2==1)
		      {
			cout << "F1[4]=" << F1[4] 
			     << " f(4,"<< c2 << "," << n1 << ")^*=" << conj(f(4,c2,n1)) 
			     << " f(4," << c4 << "," << n2 << ")=" << f(4,c4,n2) 
			     << " new cont=" << -invNq*F1[4]*conj(f(4,c2,n1))*f(4,c4,n2) << endl; 
			cout << "Dinvq(4,1,1) = " << Dinvq(4,m1,m2) << endl;
		      }
		  }

	      }

	  if(TRACE2)
	    {
	      //	      cout << "Dinvq=" << Dinvq << endl;
	      cout << "C: Dinvq[4]=" << Dinvq.qcomp(4) << endl;
	    }
	}
#else  
  for(int c2=0; c2<NC; c2++)
    for(int c4=0; c4<NC; c4++)
      {	    
	Triplet myivec=clist[c2]+clist[c4];
	
	for(int i=0; i<Nq; i++)
	  {
	    SMatrix<complex<realtype>> my(NMAT,NMAT);
	    
	    my  = g[c4];
	    my *= Kinvr[i];
	    my *= g[c2];
	    
	    int newindx=la.rAdd(la.GetInversionIndx(i),myivec);
	    
	      my *= Kinvr[newindx];
	      
	      F1[i]=0.5*NFAKESPINTRACE*tr(my);
	  } 
	FFTWEXECUTE(F1pluss); 
	
	
	for(int m1=NSUBL; m1<NDMAT; m1++) 
	  for(int m2=m1; m2<NDMAT; m2++) 
	    {
	      
	      int n1=m1-NSUBL; 
	      int n2=m2-NSUBL; 
	      
	      for(int qi=0; qi<Nq; qi++)
		{
		  Dinvq(qi,m1,m2)+=-invNq*F1[qi]*conj(f(qi,c2,n1))*f(qi,c4,n2)*conj(expi(la.qr(qi,clist[c4])));
		}
	    }
      }
#endif

 
  
  // Set the remainding phonon-phonon part.
  for(int qi=0; qi<Nq; qi++)
    for(int m1=NSUBL; m1<NDMAT; m1++) 
      for(int m2=m1; m2<NDMAT; m2++) 
	{
	  Dinvq(qi,m2,m1)=conj(Dinvq(qi,m1,m2));
	}


  if(TRACE2)
    {
      cout  << "Dinvq before adding bare phonons: " << Dinvq << endl;
    }

  // add the bare phonon part
  
  realtype barepart=1./currT; // the one-half is included in the def. of propagator 
  for(int qi=0; qi<Nq; qi++)
    {
      for(int n=0; n<NMODE; n++)
	{
	  int m=n+NSUBL;
	  Dinvq(qi,m,m) += barepart;
	}
    }
#endif

  if(TRACE2) 
    {
      cout << "Dinvq after adding bare phonon part: " << Dinvq << endl;
    }
  
  // FORCING
  MakeMixedHermitian(Dinvq,NSUBL,NSUBL);
  MakeInversionTransposedSymmetric(Dinvq);
  

  Chomp(Dinvq);
  
  if(TRACE){SanityCheck(Dinvq,"Dinvq, after constructing it",false);}
  
  if(TRACE2){ cout << "Dinvq before inversion=" << Dinvq << endl;}

  MatrixInverse(Dinvq); // B = Dq

  if(TRACE){SanityCheck(Dq,"Dq, after inverting Dinvq",false);}
  
  if(TRACE2){ cout << "Dq, after inverting Dinvq, Dq=" << Dq << endl;}

  if(excludeqzero) Setqzerotozero(Dq);


  if(preserveinput)
    {
      FFTWEXECUTE(A1r_to_A1q);  // Kinvr->Kinvq, so after transform: Aq=Kinvq*sqrt(Nq) 
      Kinvq *= invSqrtNq; // Aq=Kinvq;
      
#ifdef PRESERVESYMMETRY
      MakeSymmetric(Kinvq);
#endif
      MakeHermitian(Kinvq);
    }

  if(TRACE){SanityCheck(Dq,"Dq, at end of ComputeDq",false);}

  if(TRACE2){ cout << "At end of ComputDq, Dq=" << Dq << endl;}
  
  if(TRACE) cout << "Done with ComputeDq " << endl;
}


// ComputeSelfEnergy() computes the self-energy from the self-consistent equations
//
// Input: Kinv_q (stored in A1).
// Output: Sigma_q (stored in A2)
// 
// The routine uses in-place FFTs so info contained in A and B are modified
// On output:
// A1 = Kinvr, unless preserveinput=true, then A1=Kinvq
// A2 = Sigmaq 
//
// FFTW omits volume prefactors in fourier transforms, so it just carries out the sums
// Therefore volume factors must be included explicitly. 
void Driver::ComputeSelfEnergy(const bool preserveinput=false)
{
  if(TRACE) cout << "Starting ComputeSelfEnergy" << endl;

  if(TRACE){SanityCheck(Kinvq,"Kinvq, at start of ComputeSelfEnergy",true);}

  ComputeDq(true,false); // exclude zero mode, do not preserve input here, Ar must contain Kinvr for ComputeSelfEnergy to work.

  
  Sigmaq.SetToZero();

  for(int alpha=0; alpha<NMAT; alpha++)
    for(int delta=alpha; delta<NMAT; delta++)
      {
	//	const int alpha_s=spin(alpha);
	const int alpha_l=subl(alpha);
	
	//	const int delta_s=spin(delta);
	const int delta_l=subl(delta);

	// Term 1:
	for(int q=0; q<Nq; q++){F1[q]=Dq(q,delta_l,alpha_l);}
	FFTWEXECUTE(F1pluss); // we take the Fourier volume factor X1=0
	
	for(int r=0; r<Nq; r++){F2[r]=Kinvr(la.GetInversionIndx(r),alpha,delta)*F1[r];}
	FFTWEXECUTE(F2pluss);

	for(int k=0; k<Nq; k++){Sigmaq(k,alpha,delta) += invSqrtNq*F2[k];}

	//	cout << "Sigmaq after term 1: alpha=" << alpha << " delta=" << delta << endl;
	//cout << Sigmaq << endl;

	
#if defined PHONONS && !defined ELASTICONLY
	if(NSUBL != 1){cout << "Error: PHONONS are not implemented for NSUBL!=1, exiting" << endl; exit(1);}


#ifdef CPOSITIVE
	// Term 2
	for(int c=0; c<NC; c++)
	  {
 	    for(int q=0; q<Nq; q++)
	      {
		F1[q]=0;
		for(int n=0; n<NMODE; n++)
		  {
		    F1[q]+=Dq(q,NSUBL+n,0)*f(q,c,n);
		  }
	      }
	    FFTWEXECUTE(F1pluss);
	    
	    for(int r=0; r<Nq; r++)
	      {
		F2[r]=0;
		for(int s=0; s<NSPIN; s++)
		  F2[r]+=F1[r]*Kinvr(r,s,alpha)*g(c,s,delta); // conj(Kinv(r,s,alpha))= Kinv(r,s,alpha)) ; real		  
	      }
	    
	    FFTWEXECUTE(F2pluss);
	    
	    for(int k=0; k<Nq; k++){ Sigmaq(k,alpha,delta) += invNq*expi(la.qr(k,clist[c]))*F2[k];}

	    for(int r=0; r<Nq; r++)
	      {
		F2[r]=0;
		int rpc=la.rAdd(r,clist[c]);
		for(int s=0; s<NSPIN; s++)
		  F2[r]+=F1[r]*Kinvr(rpc,s,alpha)*g(c,delta,s); 
	      }
	    
	    FFTWEXECUTE(F2pluss);
	    
	    for(int k=0; k<Nq; k++){ Sigmaq(k,alpha,delta) += invNq*F2[k];}

	  }

	// Term 3 (as term 2, but switch alpha and delta and take complex conj)
	for(int c=0; c<NC; c++)
	  {
 	    for(int q=0; q<Nq; q++)
	      {
		F1[q]=0;
		for(int n=0; n<NMODE; n++)
		  {
		    F1[q]+=Dq(q,NSUBL+n,0)*f(q,c,n);
		  }
	      }
	    FFTWEXECUTE(F1pluss);


	    for(int r=0; r<Nq; r++)
	      {
		F2[r]=0;
		for(int s=0; s<NSPIN; s++)
		  F2[r]+=F1[r]*Kinvr(r,s,delta)*g(c,s,alpha); // conj(Kinv(r,s,alpha))= Kinv(r,s,alpha)) ; real		  
	      }
	    
	    FFTWEXECUTE(F2pluss);
	    
	    for(int k=0; k<Nq; k++){ Sigmaq(k,alpha,delta) += invNq*conj(expi(la.qr(k,clist[c]))*F2[k]);}

	    for(int r=0; r<Nq; r++)
	      {
		F2[r]=0;
		int rpc=la.rAdd(r,clist[c]);
		for(int s=0; s<NSPIN; s++)
		  F2[r]+=F1[r]*Kinvr(rpc,s,delta)*g(c,alpha,s); 
	      }
	    
	    FFTWEXECUTE(F2pluss);
	    
	    for(int k=0; k<Nq; k++){ Sigmaq(k,alpha,delta) += invNq*conj(F2[k]);}
	  }

	// Term4:
	for(int c1=0; c1<NC; c1++)
	  for(int c2=0; c2<NC; c2++)
	    {
	      for(int q=0; q<Nq; q++)
		{
		  F1[q]=0;
		  
		  for(int n1=0; n1<NMODE; n1++)
		    for(int n2=0; n2<NMODE; n2++)
		      { F1[q]+= conj(f(q,c1,n1))*f(q,c2,n2)*Dq(q,NSUBL+n2,NSUBL+n1)*expi(la.qr(q,clist[c1]));}
		}
	      FFTWEXECUTE(F1pluss);

	      
	      for(int r=0; r<Nq; r++)
		{
		  F2[r]=0;
		  for(int be=0; be<NSPIN; be++)
		    for(int ga=0; ga<NSPIN; ga++)
		      {
			F2[r]+=g(c1,alpha,be)*Kinvr(r,ga,be)*g(c2,ga,delta); // conj(Kinv(r,s2,s1))=Kinv(r,s2,s1)
		      }
		  F2[r] *= F1[r];
		}
	      
	      FFTWEXECUTE(F2pluss);
	      for(int k=0; k<Nq; k++){Sigmaq(k,alpha,delta)+= -invNq*invSqrtNq*F2[k]*expi(la.qr(k,clist[c1]))*expi(la.qr(k,clist[c2]));}

	      
	      for(int r=0; r<Nq; r++)
		{
		  int rpc2=la.rAdd(r,clist[c2]);
		  F2[r]=0;
		  for(int be=0; be<NSPIN; be++)
		    for(int ga=0; ga<NSPIN; ga++)
		      {
			F2[r]+=g(c1,alpha,be)*Kinvr(rpc2,ga,be)*g(c2,delta,ga); 
		      }
		  F2[r] *= F1[r];
		}
	      
	      FFTWEXECUTE(F2pluss);
	      for(int k=0; k<Nq; k++){Sigmaq(k,alpha,delta)+= -invNq*invSqrtNq*F2[k]*expi(la.qr(k,clist[c1]));}


	      for(int r=0; r<Nq; r++)
		{
		  int rpc1=la.rAdd(r,clist[c1]);
		  F2[r]=0;
		  for(int be=0; be<NSPIN; be++)
		    for(int ga=0; ga<NSPIN; ga++)
		      {
			F2[r]+=g(c1,be,alpha)*Kinvr(rpc1,ga,be)*g(c2,ga,delta); 
		      }
		  F2[r] *= F1[r];
		}
	      
	      FFTWEXECUTE(F2pluss);
	      for(int k=0; k<Nq; k++){Sigmaq(k,alpha,delta)+= -invNq*invSqrtNq*F2[k]*expi(la.qr(k,clist[c2]));}

	      for(int r=0; r<Nq; r++)
		{
		  int rpc1c2=la.rAdd(r,clist[c1]+clist[c2]);
		  F2[r]=0;
		  for(int be=0; be<NSPIN; be++)
		    for(int ga=0; ga<NSPIN; ga++)
		      {
			F2[r]+=g(c1,be,alpha)*Kinvr(rpc1c2,ga,be)*g(c2,delta,ga); 
		      }
		  F2[r] *= F1[r];
		}
	      
	      FFTWEXECUTE(F2pluss);
	      for(int k=0; k<Nq; k++){Sigmaq(k,alpha,delta)+= -invNq*invSqrtNq*F2[k];}
	    }
#else	
	// Term 2
	for(int c=0; c<NC; c++)
	  {
 	    for(int q=0; q<Nq; q++)
	      {
		F1[q]=0;
		for(int n=0; n<NMODE; n++)
		  {
		    F1[q]+=Dq(q,NSUBL+n,0)*f(q,c,n);
		  }
	      }
	    FFTWEXECUTE(F1pluss);

	    for(int r=0; r<Nq; r++)
	      {
		F2[r]=0;
		for(int s=0; s<NSPIN; s++)
		  F2[r]+=F1[r]*Kinvr(r,s,alpha)*g(c,s,delta); // conj(Kinv(r,s,alpha))= Kinv(r,s,alpha)) ; real		  
	      }
	    
	    FFTWEXECUTE(F2pluss);
	    
	    for(int k=0; k<Nq; k++){ Sigmaq(k,alpha,delta) += invNq*expi(la.qr(k,clist[c]))*F2[k];}
	  }

	// Term 3 (as term 2, but switch alpha and delta and take complex conj)
	for(int c=0; c<NC; c++)
	  {
 	    for(int q=0; q<Nq; q++)
	      {
		F1[q]=0;
		for(int n=0; n<NMODE; n++){ F1[q]+=Dq(q,NSUBL+n,0)*f(q,c,n);}
	      }
	    FFTWEXECUTE(F1pluss);
	    
	    for(int r=0; r<Nq; r++)
	      {
		F2[r]=0;
		for(int s=0; s<NSPIN; s++)
		  //		  F2r[r]+=F1r[r]*conj(Kinvr(r,s,delta))*g(c,s,alpha); // Kinv(-r,alpha,s) = Kinv(r,s,alpha)
		  F2[r]+=F1[r]*Kinvr(r,s,delta)*g(c,s,alpha); // conj(Kinv(r,s,delta))= Kinv(r,s,delta)
	      }
	    FFTWEXECUTE(F2pluss);
	    
	    for(int k=0; k<Nq; k++){ Sigmaq(k,alpha,delta) += invNq*conj(expi(la.qr(k,clist[c]))*F2[k]);}
	  }
	
	// Term 4:	
	for(int c1=0; c1<NC; c1++)
	  for(int c2=0; c2<NC; c2++)
	    {
	      for(int q=0; q<Nq; q++)
		{
		  F1[q]=0;
		  
		  for(int n1=0; n1<NMODE; n1++)
		    for(int n2=0; n2<NMODE; n2++)
		      { F1[q]+= conj(f(q,c1,n1))*f(q,c2,n2)*Dq(q,NSUBL+n2,NSUBL+n1)*expi(la.qr(q,clist[c1]));}
		}
	      FFTWEXECUTE(F1pluss);
	      
	      for(int r=0; r<Nq; r++)
		{
		  F2[r]=0;
		  for(int s1=0; s1<NSPIN; s1++)
		    for(int s2=0; s2<NSPIN; s2++)
		      {
			F2[r]+=g(c1,alpha,s1)*Kinvr(r,s2,s1)*g(c2,s2,delta); // conj(Kinv(r,s2,s1))=Kinv(r,s2,s1)
		      }
		  F2[r] *= F1[r];
		}
	      
	      FFTWEXECUTE(F2pluss);
	      for(int k=0; k<Nq; k++){Sigmaq(k,alpha,delta)+= -invNq*invSqrtNq*F2[k]*expi(la.qr(k,clist[c1]))*expi(la.qr(k,clist[c2]));}
	    }
#endif
#endif
      }
  
#ifdef PRESERVESYMMETRY
  MakeSymmetric(Sigmaq);
#endif  
  MakeHermitian(Sigmaq,false);
  //TRY
  MakeSpinDiagonal(Sigmaq,false);

  if(TRACE) SanityCheck(Sigmaq,"Sigmaq after explicitly constructing it");
  
  if(preserveinput)
    {
      FFTWEXECUTE(A1r_to_A1q); // Kinvr->Kinvq, so after transform: A1q=Kinvq*sqrt(Nq)
      Kinvq *= invSqrtNq;

#ifdef PRESERVESYMMETRY
      MakeSymmetric(Kinvq);
#endif
      MakeHermitian(Kinvq);
    }
  
  if(TRACE) cout << "Done with ComputeSelfEnergy " << endl;
}




// Routine to build Kinq from Jq and Sigmaq
//
//
void Driver::ConstructKinvq()
{
  if(TRACE) cout << "Starting ConstructKinvq" << endl;
  
  Kq =  Jq;

  if(TRACE) cout << "Initializing Kq" << endl;

  if(TRACE){SanityCheck(Kq,"Kq, after setting Jq");}

  Kq += Sigmaq;

  if(TRACE){SanityCheck(Kq,"Kq, after adding Sigmaq");}

  
#ifdef PHONONS
  if(TRACE) cout << "Adding elastic modes " << endl;
  for(int j=0; j<NELASTIC; j++)
    {
      VecMat<complex<realtype>,NMAT,NMAT> temp(*rule.gelptrs[j]);
      temp *= epsilon[j];
      Kq += temp;
    }
#endif

  

  if(TRACE){SanityCheck(Kq,"Kq, after adding Elastic modes");}
  
  mineigenvalue=SubtractMinimumEigenvalue(Kq);
  
  if(TRACE){SanityCheck(Kq,"Kq, after subtracting eigenvalues");}

  AddDelta(Kq,Delta);

  if(TRACE){SanityCheck(Kq,"Kq, after adding Delta");}
  
  // construct K1inv
  MatrixInverse(Kq); // K = Kinv_q 

  // Force correct properties on the matrix
  MakeHermitian(Kinvq); 
  MakeInversionTransposedSymmetric(Kinvq); 
  // TRY
  MakeSpinDiagonal(Kinvq);

  if(TRACE){SanityCheck(Kinvq,"Kinvq, after inverting");}
  
  if(TRACE) cout << "Done with ConstructKinvq" << endl;
}



void Driver::MakeRandomSigma()
{
  if(TRACE) cout << "Starting MakeRandomSigma()" << endl;
  logfile << "Making a random initialization of the self-energy" << endl;

  realtype da = par[DA];
  realtype de = Delta[0]; //scale randomness with Delta

  for(int m=0; m<NMAT; m++)
    {
      for(int i=0; i<Nq; i++)
	{
	  complex<realtype> c=da*de*complex<realtype>(RAN(),0.); // real positive value 
	  Sigmaq(i,m,m)=c;
	}
    }

  if(NMAT>1)
    {
      for(int m1=0; m1<NMAT; m1++)
	for(int m2=m1+1; m2<NMAT; m2++)
	  {
	    for(int i=0; i<Nq; i++)
	      {
		complex<realtype> c=da*de*complex<realtype>(RAN(),RAN());
		Sigmaq(i,m1,m2)=c;
	      }
	  }
    }


  MakeHermitian(Sigmaq,false); // no warnings
  //MakeHermitian(Sigmaq);

  //  cout << "Sigmaq after MakeHermitian..." << endl;
  // cout << Sigmaq << endl;

  //  MakeSpinSymmetric(Sigmaq,false);
  // TRY
  MakeSpinDiagonal(Sigmaq,false);
  //  MakeSpinDiagonalEqual(Sigmaq,false);

#ifdef FORCEINVERSIONSYMMETRY
  MakeInversionTransposedSymmetric(Sigmaq);
  /*
  FFTWEXECUTE(A2q_to_A2r);
  MakeReal(Sigmar);
  FFTWEXECUTE(A2r_to_A2q);
  Sigmaq*=invNq;  // do not change magnitude
  */
#endif

  //  cout << "Sigmaq after MakeInversionTrans..." << endl;
  //  cout << Sigmaq << endl;

#ifdef PRESERVESYMMETRY
  MakeSymmetric(Sigmaq);
#endif


  
  if(TRACE) SanityCheck(Sigmaq,"Sigmaq, at end of MakeRandomSigma");
  if(TRACE) cout << "Finished MakeRandomSigma()" << endl;
}


void Driver::SetQsToZero()
{
  logfile << "Setting the following qpts to zero in the self-energy" << endl;

  bool found=false;
  ifstream ifile("qstozero.in");
  
  while(ifile)
    {
      int qindx;
      ifile >> qindx;
      if(qindx >=0 && qindx <Nq && ifile)
	{
	  found=true;
	  Coord thisq=la.qPos(qindx);
	  logfile << "(" << thisq << ")" << endl;
	  Sigmaq(qindx,0,0)=complex<realtype>(0.,0.);
	}        
    }
  if(!found) logfile << "None" << endl;
}



// assumes initialized Sigma and mu
void Driver::SolveSelfConsistentEquation(NumberList Delta)
{
  if(TRACE) cout << "Starting SolveSelfConsistentEquation " << endl;

  //  lineid++;


#ifdef RANDOMINITIALIZATION
  MakeRandomSigma();
#elif defined USELASTSIGMA
  if(!SigmaInitialized){MakeRandomSigma(); SigmaInitialized=true; MaxIterMultiplier=20;}
  else{MaxIterMultiplier=1;}
#else
  rule.InitializeSigma(Sigmaq); // Get initial values of Sigmaq
#endif
  SetQsToZero();

#ifdef PRESERVESYMMETRY
  MakeSymmetric(Sigmaq);
#endif

#ifdef BIASQS  // set some q's to zero in 

#endif

  Chomp(Jq); // set very small entries to 0

  if(TRACE) SanityCheck(Jq,"Jq, input to SolveSelfConsistentEquation");

  
  //  if(TRACE) cout << "Min eigenvalue: " << FindMinimumEigenvalue(Jq) << endl;

  SubtractMinimumEigenvalue(Jq);

  Chomp(Jq); // set very small entries to 0

  if(TRACE) SanityCheck(Jq,"Jq, after subtracting minimum");

  /*
#ifdef PRINTTCONVERGENCE
  ostringstream sstrs; sstrs << "DELTA" << Delta << ".s.dat";
  ostringstream tstrs; tstrs << "DELTA" << Delta << ".t.dat";
  ostringstream ustrs; ustrs << "DELTA" << Delta << ".u.dat";
  string tfilename = tstrs.str();
  string sfilename = sstrs.str();
  string ufilename = ustrs.str();
  ofstream tfile(tfilename.c_str());
  ofstream sfile(sfilename.c_str());
  ofstream ufile(ufilename.c_str());
#endif
  */

  // Put together K 

  
  realtype newT=-1.;
  realtype m2=0.; // magnetic order parameter squared.
  vector<obstype> nobs(NOBSERVABLES); // nematic order parameters
  vector<obstype> nspinobs(NSPINOBSERVABLES); // different types of spin order parameters

  //  vector<obstype> nalphas(NOBSERVABLES); // alpha


  // construct Kinvq:
  
  ConstructKinvq();

  realtype oldT=CalculateT(0); 

  currT=oldT; // set the current operating temperature

  realtype absTdev(0.);
  realtype oldabsTdev(0.);
  int nincreases(0); // number of times that a new iteration gets a T that deviates more from the oldT than the previous iteration. Sign of no convergence.


  vector<realtype> Ts(NSUBL); // A list of Ts
  
  if(TRACE)
    {
      cout << "Initial temperatures: " << endl;
      CalculateTs(Ts); // calculate all the temperatures
      
      streamsize ss=cout.precision();
      cout.precision(17);
      for(int i=0; i<NSUBL; i++)
	cout << Ts[i] << " ";
      cout << endl;
      cout.precision(ss); // restore precision
    }



 
  int iter=0;
  converged=false;
  bool pconverged=true; // convergence in the previous iteration,
  bool done=false;

  bool reachedMAXITER=false;

 
  while(!done && !reachedMAXITER)  
    {
      if(TRACE) cout << "New iteration: " << iter << endl;
      iter++;
      if(converged) pconverged=true; // record if the previous iteration had converged

      if(TRACE) SanityCheck(Kinvq,"Kinvq, at start of new iteration");

      
      ComputeSelfEnergy(); // careful with this, it overwrites K

#ifdef NOSELFENERGY
      Sigmaq.SetToZero();
#endif

      if(TRACE) SanityCheck(Sigmaq,"Sigmaq, after ComputeSelfEnergy");
      
      ConstructKinvq();

      if(TRACE) SanityCheck(Kinvq,"Kinvq, in iteration loop after making sigmaq");

      // convergence checks
      CalculateTs(Ts); // calculate all the temperatures      
      newT=Ts[0];

      absTdev=fabs((newT-oldT)/oldT);
      if(absTdev > oldabsTdev)
	{
	  if(iter>40) nincreases++; // only increment this after 40 iterations.
	  if(nincreases > 0)
	    PRINTPROGRESSTICKLER=1; // monitor every step.
	  else
	    PRINTPROGRESSTICKLER=10;
	}
      oldabsTdev=absTdev;

            
#ifdef PHONONS
      NumberList epsoverT=CalculateEpsilonsOverT();
#endif

      // monitor maximum values of Kinv
      if(NSPIN>1)
	{
	  SMatrix<realtype,NMAT,NMAT> maxvals;
	  SMatrix<int,NMAT,NMAT> maxqs;
	  
	  FindMaxVals(maxqs,maxvals);

	  ofstream outfile(CONVERGENCEMONITORNAME,ios::app);
	  outfile << iter << " ";
	  for(int s=0; s<NSPIN; s++)
	    {
	      realtype outmax = maxvals(mindx(s,0),mindx(s,0));
	      outfile << outmax << " ";
	    }
	  outfile << endl;
	  outfile.close();
	}


      if(TRACE)
	{
	  cout << "iteration: " << iter << " T= " << newT << " oldT=" << oldT
	       << " dev: " << absTdev;
#ifdef PHONONS
	  cout << " epsilon/T= " << epsoverT;
#endif
	  cout << endl;
	}
      
      // write progress to logfile
      if(PRINTPROGRESS && iter%PRINTPROGRESSTICKLER==0)
	{
	  logfile.precision(9);
	  logfile << "iteration: " << iter << " T= "; 
	  for(int i=0; i<NSUBL; i++) logfile << Ts[i] << " ";
	  logfile << " oldT=" << oldT << " dev: " << absTdev << " ninc: " << nincreases;
#ifdef PHONONS
	  logfile << " epsilon/T= " << epsoverT;
#endif
	  logfile << endl;
	}
      
      if( absTdev < par[TOLERANCE]) converged=true;

      //      const realtype inertia=0.5; // how much to resist changes: [0,1] 
      const realtype inertia=0.2; // how much to resist changes: [0,1] 

      currT= (1.-inertia)*newT+inertia*oldT; 

      oldT=currT;

#ifdef PHONONS      
      for(int i=0; i<NELASTIC; i++){epsilon[i]= (1.-inertia)*currT*epsoverT[i]+inertia*epsilon[i];}
#endif       
      
      if(TRACE) cout << "converged= " << converged << " TOLERANCE:" << par[TOLERANCE] << endl;
      
      if(converged && pconverged){ done=true; continue;} // two iterations must fulfill conv. crit.

      if(nincreases > 4){ done =true; logfile << "Too many nincr=" << nincreases << ", Not converging, exiting" << endl; continue;} // to many increases in Tdiff, so exiting without converging.


      reachedMAXITER=(iter >= MaxIterMultiplier*par[MAXITER]);
    }
  
  //  if(TRACE) cout << "Final Kinv_q: " << Kinvq << endl;  

  //  nalphas=CalculateAlphas(newT); // calculate alphas

  m2=NFAKESPINTRACE*NSPIN*newT/(2.*Delta[0]*Nq); // calculate magnetic moment

  CalculateTs(Ts); // calculate the final temperatures
  
  // Print final Ts and epsilons
  logfile << "iteration: " << iter << " Ts: ";
  streamsize ss=logfile.precision();
  logfile.precision(17);
  for(int i=0; i<NSUBL; i++){logfile << Ts[i] << " ";}
#ifdef PHONONS
  logfile << " epsilon: " << epsilon;
#endif
  logfile << " converged: " << (converged ? "true": "false") << endl;
  logfile.precision(ss);

      
#ifdef PRINTTCONVERGENCE
      tfile << setprecision(17) << newT << endl;
      sfile << setprecision(17) << nobs[2].real() << endl;
#endif
      

  if(converged) logfile << "Convergence reached after " << iter << " steps." << endl;

  //Use MAXITER as convergence criterion itself
  if(par[TOLERANCE]==0.)
    {
      logfile << "Exiting after MAXITER steps, using end result" << endl;
      converged=true;
    } 


  if(!converged)
    {
      if(reachedMAXITER)
	logfile << "reached MAXITER=" << par[MAXITER] << " iterations without converging, increase MAXITER!" << endl;
      else
	logfile << "no uniform convergence of iterations, exiting!" << endl;
    }
  else
    {      
      lineid++; 

      SMatrix<realtype,NMAT,NMAT> maxvals;
      SMatrix<int,NMAT,NMAT> maxqs;

      FindMaxVals(maxqs,maxvals);

      //Writing diagonal magnetization for different spins and sublattices.
      for(int l1=0; l1<NSUBL; l1++)
	for(int s1=0; s1<NSPIN; s1++)
	  {
	    stringstream ss;
	    ss << "m" << l1 << l1 << "_" << s1 << s1 << ".dat";
	    
	    int m1=mindx(s1,l1);
	    int m2=mindx(s1,l1);
	    realtype mm=0.5*newT*invNq*maxvals(m1,m2);
	  
	    ofstream outfile(ss.str().c_str(),ios::app);
	    outfile << setprecision(16) << newT << " " << mm << endl;
	    outfile.close();
	  }




      ofstream peakfile(PEAKREPORTNAME.c_str(),ios::app);
      
      peakfile << "Delta=" << Delta[0] << " T=" << newT << " Max K^-1  @ q:" << endl; 
      for(int l1=0; l1<NSUBL; l1++)
	{
	  peakfile << "(l1=" << l1 << ",l2=" << l1 << ") ";

	  for(int s1=0; s1<NSPIN; s1++)
	    {
	      int m1=mindx(s1,l1);
	      int m2=mindx(s1,l1);
	      realtype mval=maxvals(m1,m2);
	      Coord q=la.qPos(maxqs(m1,m2));
	      if(mval > 0)
		peakfile << "(" << s1 << "," << s1 << "):" << mval << " @q=" << q << " ";
	    }
	  peakfile << endl;
	}

      for(int l1=0; l1<NSUBL; l1++)
	{
	  peakfile << "(l1=" << l1 << ",l2=" << l1 << ") ";

	  for(int s1=0; s1<NSPIN; s1++)
	    for(int s2=s1+1; s2<NSPIN; s2++)
	    {
	      int m1=mindx(s1,l1);
	      int m2=mindx(s2,l1);
	      realtype mval=maxvals(m1,m2);
	      Coord q=la.qPos(maxqs(m1,m2));
	      if(mval > 0)
		peakfile << "(" << s1 << "," << s2 << "):" << mval << " @q=" << q << " ";
	    }
	  peakfile << endl;
	}

      for(int l1=0; l1<NSUBL; l1++)
	for(int l2=l1+1; l2<NSUBL; l2++)
	  {
	    peakfile << "(l1=" << l1 << ",l2=" << l2 << ") ";
	    for(int s1=0; s1<NSPIN; s1++)
	      {
		int m1=mindx(s1,l1);
		int m2=mindx(s1,l2);
		realtype mval=maxvals(m1,m2);
		Coord q=la.qPos(maxqs(m1,m2));
		if(mval > 0)
		  peakfile << "(" << s1 << "," << s1 << "):" << mval << " @q=" << q << " ";
	      }
	    peakfile << endl;
	  }

      for(int l1=0; l1<NSUBL; l1++)
	for(int l2=l1+1; l2<NSUBL; l2++)
	  {
	    peakfile << "(l1=" << l1 << ",l2=" << l2 << ") ";
	    
	    for(int s1=0; s1<NSPIN; s1++)
	      for(int s2=s1+1; s2<NSPIN; s2++)
		{
		  int m1=mindx(s1,l1);
		  int m2=mindx(s2,l2);
		  realtype mval=maxvals(m1,m2);
		  Coord q=la.qPos(maxqs(m1,m2));
		  if(mval > 0)
		    peakfile << "(" << s1 << "," << s2 << "):" << mval << " @q=" << q << " ";
		}
	    peakfile << endl;
	  }

      peakfile << endl; // extra space in report
      peakfile.close();
      //      lineid++;



      //      if(Printinfo)
	{

	  /*
	  vector<record> lv=FindLowestValues(K1,nvals);
	  
	  logfile << "largest values" << endl;
	  for(int i=0; i<nvals; i++)
	    {
	      const int p=lv[i].pos;
	      Coord q=lattice.qPos(p);
	      Coord qoverpi=(1./PI)*q;
	      logfile << q << " (" << qoverpi << ")  : " << 1./lv[i].value << endl;  
	    }
	  
	  
	  logfile << "Printing info to " << MAXQNAME;
	  ofstream maxq(MAXQNAME.c_str(),ios::app);
	  Coord maximumq=lattice.qPos(lv[0].pos);
	  maxq << setprecision(16) << newT << " " << maximumq << " " << endl;
	  */


#ifdef REDUCEDOUTPUT

	  if(PRINTQCORRS)
	    {
	      stringstream ss;
	      ss << QCORRSFILENAME << "_" << lineid << ".dat";
	      
	      logfile << ss.str();

	      //	      ofstream qcorrfile(QCORRSFILENAME.c_str(),ios::app);
	      ofstream qcorrfile(ss.str().c_str());

	      realtype factor=newT*NSPIN*0.5;
	      complex<realtype>* start=Kinvq(0,0);
	      
	      if(BINARYOUTFILES)
		{
		  qcorrfile.write((char*) &lineid,sizeof(lineid)); // general format
		  qcorrfile.write((char*) &la.nindx_q,sizeof(la.nindx_q));
		  for(int i=0; i<la.nindx_q; i++)
		    {
		      complex<realtype> value=factor*start[la.indx_site_q[i]];
		      qcorrfile.write((char*) &value,sizeof(value));
		    }
		}
	      else
		{
		  //		  qcorrfile << setprecision(16) << lineid << " ";
		  for(int i=0; i<la.nindx_q; i++)
		    {
		      complex<realtype> value=factor*start[la.indx_site_q[i]];
		      qcorrfile << real(value) << " " << imag(value) << endl;;}
		}
	      qcorrfile.close();
	    }
#endif
	  /*
	  if(PRINTSIGMAE)
	    {
	      logfile << ", " << SIGMAEFILENAME;
	      ofstream sigmaefile(SIGMAEFILENAME.c_str(),ios::app);
	      
	      realtype min=SubtractMinimum(Sigma);

	      if(BINARYOUTFILES)
		{
		  sigmaefile.write((char*) &lineid,sizeof(lineid)); // general format
		  sigmaefile.write((char*) &lattice.nindx_q,sizeof(lattice.nindx_q));
		  for(int i=0; i<lattice.nindx_q; i++)
		    {
		      sigmaefile.write((char*) &Sigma[lattice.indx_site_q[i]],sizeof(Sigma[0]));
		    }
		}
	      else
		{
		  sigmaefile << setprecision(16) << lineid << " ";
		  for(int i=0; i<lattice.nindx_q; i++){sigmaefile << Sigma[lattice.indx_site_q[i]] << " ";}
		  sigmaefile << endl;
		}
	      sigmaefile.close();
	    }
	  	  
	  if(PRINTRCORRS)
	    {
	      logfile << ", " << RCORRSFILENAME;
	      vector<realtype> rcorr(Vf);
	      RealSpaceCorrelationFunction(K1,newT,NS,dim,dims,rcorr);
	      ofstream rcorrfile(RCORRSFILENAME.c_str(),ios::app);
	      if(BINARYOUTFILES)
		{
		  rcorrfile.write((char*) &lineid,sizeof(lineid)); // general format
		  rcorrfile.write((char*) &lattice.nindx_r,sizeof(lattice.nindx_r));
		  for(int i=0; i<lattice.nindx_r; i++)
		    {
		      rcorrfile.write((char*) &rcorr[lattice.indx_site_r[i]],sizeof(rcorr[0]));
		    }
		}
	      else
		{
		  rcorrfile << setprecision(16) << lineid << " ";
		  for(int i=0; i<lattice.nindx_r; i++){rcorrfile << rcorr[lattice.indx_site_r[i]] << " ";}
		}
	      rcorrfile.close();

	      ofstream mostremotefile(RCORRMOSTREMOTE.c_str(),ios::app);
	      mostremotefile << setprecision(16) << newT << " " << rcorr[lattice.indx_most_remote] << endl;
	      mostremotefile.close();

	    }
#else
      // Full output	  
	  if(PRINTQCORRS)
	    {
	      logfile << ", " << QCORRSFILENAME;
	      ofstream qcorrfile(QCORRSFILENAME.c_str(),ios::app);
	      realtype factor=newT*NS*0.5;
	      vector<realtype> qcorr(Nq);
	      //	      for(int i=0; i<Nq; i++) qcorr[i]=factor/K1[i];
	      
	      if(BINARYOUTFILES)
		{
		  qcorrfile.write((char*) &lineid,sizeof(lineid)); // general format
		  qcorrfile.write((char*) &Nq,sizeof(Nq));
		  qcorrfile.write((char*) &qcorr[0],Nq*sizeof(qcorr[0]));
		}
	      else
		{
		  qcorrfile << setprecision(16) << lineid << " ";
		  for(int i=0; i<Nq; i++){ qcorrfile << qcorr[i] << " ";}
		  qcorrfile << endl;
		}
	      qcorrfile.close();
	    }

	  if(PRINTSIGMAE)
	    {
	      logfile << ", " << SIGMAEFILENAME;
	      ofstream sigmaefile(SIGMAEFILENAME.c_str(),ios::app);
	      
	      //	      vector<realtype> min=SubtractMinimum(Sigma);

	      if(BINARYOUTFILES)
		{
		  sigmaefile.write((char*) &lineid,sizeof(lineid)); // general format
		  sigmaefile.write((char*) &Nq,sizeof(Nq));
		  sigmaefile.write((char*) &Sigma[0],Nq*sizeof(Sigma[0]));
		}
	      else
		{
		  sigmaefile << setprecision(16) << lineid << " ";
		  for(int i=0; i<Nq; i++){ sigmaefile << Sigma[i] << " ";}
		  sigmaefile << endl;
		}
	      sigmaefile.close();
	    }
	  
	  
	  if(PRINTRCORRS)
	    {
	      logfile << ", " << RCORRSFILENAME;
	      vector<realtype> rcorr(Vf);
	      //	      RealSpaceCorrelationFunction(K1,newT,NS,dim,dims,rcorr);
	      ofstream rcorrfile(RCORRSFILENAME.c_str(),ios::app);
	      if(BINARYOUTFILES)
		{
		  rcorrfile.write((char*) &lineid,sizeof(lineid)); // general format
		  rcorrfile.write((char*) &Vf,sizeof(Vf));
		  rcorrfile.write((char*) &rcorr[0],Vf*sizeof(rcorr[0]));
		}
	      else
		{
		  rcorrfile << setprecision(16) << lineid << " ";
		  for(int i=0; i<Vf; i++){ rcorrfile << rcorr[i] << " ";}
		  rcorrfile << endl;
		}
	      rcorrfile.close();

	      ofstream mostremotefile(RCORRMOSTREMOTE.c_str(),ios::app);
	      mostremotefile << setprecision(16) << newT << " " << rcorr[lattice.indx_most_remote] << endl;
	      mostremotefile.close();


	    }
#endif // REDUCEDOUTPUT 
	  logfile << endl;
	  */
	}

      

      for(int s=0; s<NSUBL; s++)
	{
	  stringstream ss;
	  ss << "td" << "_" << s << ".dat";
	  
	  ofstream outfile_a(ss.str().c_str(),ios::app);
	  outfile_a << setprecision(16) << Ts[s] << " " << Delta[s] << endl;
	  outfile_a.close();
	  
	  ss.str("");
	  ss << "dt" << "_" << s << ".dat";
	  
	  ofstream outfile_b(ss.str().c_str(),ios::app);
	  outfile_b << setprecision(16) << Delta[s] << " " << Ts[s] << endl;
	  outfile_b.close();
	}

#ifdef PHONONS
      for(int i=0; i<NELASTIC; i++)
	{
	  stringstream ss;
	  ss << "teps" << "_" << i << ".dat";
	  ofstream outfile_a(ss.str().c_str(),ios::app);
	  outfile_a << setprecision(16) << Ts[0] << " " << epsilon[i] << endl;
	  outfile_a.close();	
	}
#endif

      /*
      
      for(int s1=0; s1<NSUBL; s1++)
	for(int s2=0; s2<NSUBL; s2++)
	  {
	    nobs=CalculateOrderPars(newT,s1,s2); // calculate order pars.	    

	    stringstream ss;

	    for(int j=0; j<NOBSERVABLES; j++)
	      {
		ss.str("");
		ss << NAMESOFOBSERVABLES[j] << "_" << s1 << s2 << ".dat";

		ofstream outfile_a(ss.str().c_str(),ios::app);
		outfile_a << setprecision(16) << newT << " "
			  << real(nobs[j]) << " " << imag(nobs[j]) << endl;
		outfile_a.close();

		ss.str("");
		ss << NAMESOFOBSERVABLES[j] << "1_" << s1 << s2 << ".dat";

		ofstream outfile_b(ss.str().c_str(),ios::app);
		outfile_b << setprecision(16) << newT << " "
			  << abs(nobs[j]) << endl;
		outfile_b.close();

		ss.str("");
		ss << NAMESOFOBSERVABLES[j] << "2_" << s1 << s2 << ".dat";

		ofstream outfile_c(ss.str().c_str(),ios::app);
		outfile_c << setprecision(16) << newT << " "
			  << norm(nobs[j]) << endl;
		outfile_c.close();
	      }
	  }
  
      */

      {
	nspinobs=CalculateSpinOrderPars(newT); // calculate order pars.	    
	
	stringstream ss;
	
	for(int j=0; j<NSPINOBSERVABLES; j++)
	  {
	    ss.str("");
	    ss << NAMESOFSPINOBSERVABLES[j] << ".dat";
	    
	    ofstream outfile_a(ss.str().c_str(),ios::app);
	    outfile_a << setprecision(16) << newT << " "
		      << real(nspinobs[j]) << " " << imag(nspinobs[j]) << endl;
	    outfile_a.close();
	    
	    ss.str("");
	    ss << NAMESOFSPINOBSERVABLES[j] << ".abs.dat";
	    
	    ofstream outfile_b(ss.str().c_str(),ios::app);
	    outfile_b << setprecision(16) << newT << " "
		      << abs(nspinobs[j]) << endl;
	    outfile_b.close();
	    
	    ss.str("");
	    ss << NAMESOFSPINOBSERVABLES[j] << ".norm.dat";
	    
	    ofstream outfile_c(ss.str().c_str(),ios::app);
	    outfile_c << setprecision(16) << newT << " "
		      << norm(nspinobs[j]) << endl;
	    outfile_c.close();
	  }
      }
   
  


      
      logfile << "Magnetic moment: " << m2 << endl;
      
      stringstream ss;
      ss << "m2" << ".dat";
      
      ofstream outfilem(ss.str().c_str(),ios::app);
      outfilem << setprecision(16) << newT << " " << m2 << endl;
      outfilem.close();
      
      
      realtype f=CalculateFreeEnergy(newT);
      logfile << "Free energy: " << f << endl;
      
      ss.str("");
      ss << "tf" << ".dat";
      
      ofstream outfile(ss.str().c_str(),ios::app);
      outfile << setprecision(16) << newT << " " << f << endl;
      outfile.close();
      
      ss.str("");
      ss << "df" << ".dat";
      
      ofstream outfile2(ss.str().c_str(),ios::app);
      outfile2 << setprecision(16) << Delta[0] << " " << f << endl;
      outfile2.close();
    }
  
  if(TRACE) cout << "Done SolveSelfConsistentEquation " << endl;
}





class Simulation{
  friend ostream& operator<<(ostream& os,Simulation& s){
    os << endl; return os;}
 public:
  Simulation();
  void Run();
 private:
  Couplings couplings;
  Rule rule;

  Driver mysolver;
  vector<NumberList> Deltalist;
  vector<NumberList> Deltastoshowlist;
  vector<bool> Printinfolist;
  vector<NumberList> epsilonlist;
};


Simulation::Simulation(): couplings(par,NC),rule(couplings),mysolver(rule),Deltalist(0),Deltastoshowlist(0),Printinfolist(0),epsilonlist(0)
{
  if(TRACE) cout << "Initializing Simulation" << endl;

  ifstream parameterfile(PARAMETERFILENAME.c_str());
  if(!parameterfile)
    {
      if(TRACE)
	cout << "No file " << PARAMETERFILENAME << " found." << endl;
      logfile << "No file " << PARAMETERFILENAME << " found." << endl;

      bool logscale=int(par[DELTALOGSCALE])==1;

      int NDeltas=par[NDELTAS];
      
      realtype D2=(logscale ? log(par[DELTASLUTT]) : par[DELTASLUTT]);
      realtype D1=(logscale ? log(par[DELTASTART]) : par[DELTASTART]);

      realtype dDelta=(NDeltas>1 ? (D2-D1)/(NDeltas-1):0);
      
      if(TRACE) 
        cout << "Creating " << NDeltas << " on " << (logscale ? "log" : "linear") << " scale from " 
             << par[DELTASTART] << " to " << par[DELTASLUTT] << endl;
      logfile << "Creating " << NDeltas << " on " << (logscale ? "log" : "linear") << " scale from " 
	      << par[DELTASTART] << " to " << par[DELTASLUTT] << endl;

      for(int i=0; i<NDeltas; i++)
        {
          realtype this_Delta=(logscale ? exp(D1+dDelta*i): D1+dDelta*i);
          NumberList myval(NSUBL,this_Delta);
          Deltalist.push_back(myval); 
          Printinfolist.push_back(true);
        }
    }
  else
    {
      if(TRACE) cout << "Reading " << PARAMETERFILENAME << " from disk" << endl;
      logfile << "Reading " << PARAMETERFILENAME << " from disk" << endl;


      
      string line;
      while (getline(parameterfile, line))
	{
	  istringstream iss(line);
	  realtype newval;
	  if(!(iss >> newval)){ break;}
	  
	  NumberList newDelta(NSUBL,newval); // initialize with one value for all
	  
	  for(int i=1; i<NSUBL; i++)
	    {
	      if(!(iss >> newval)){ break;}
	      newDelta[i]=newval;
	    }
	  Deltalist.push_back(newDelta);
	  Printinfolist.push_back(false);
	}    
    }
  if(TRACE) cout << "Deltalist has " << Deltalist.size() << " entries" << endl;
  logfile << "Deltalist has " << Deltalist.size() << " entries" << endl;  
  
  
  ifstream parameterfile2(DELTASTOSHOWFILENAME.c_str());
  if(!parameterfile)
    {
      if(TRACE) 
	cout << "No file " << DELTASTOSHOWFILENAME << " found." 
	     << " Using Delta=" << par[DELTASTART] << endl;
      logfile << "No file " << DELTASTOSHOWFILENAME << " found." 
	      << " Using Delta=" << par[DELTASTART] << endl;
      
      NumberList myval(NSUBL,par[DELTASTART]);
      Deltastoshowlist.push_back(myval); 
    }
  else
    {
      if(TRACE) cout << "Reading " << DELTASTOSHOWFILENAME << " from disk" << endl;
      logfile << "Reading " << DELTASTOSHOWFILENAME << " from disk" << endl;
      
      string line;
      while (getline(parameterfile, line))
	{
	  istringstream iss(line);
	  realtype newval;
	  if(!(iss >> newval)){ break;}
	  
	  NumberList newDelta(NSUBL,newval); // initialize with one value for all
	  
	  for(int i=1; i<NSUBL; i++)
	    {
	      if(!(iss >> newval)){ break;}
	      newDelta[i]=newval;
	    }
	  Deltastoshowlist.push_back(newval);
	}
      
    }
  if(TRACE) cout << "Deltastoshowlist has " << Deltastoshowlist.size() << " entries" << endl;
  logfile << "Deltastoshowlist has " << Deltastoshowlist.size() << " entries" << endl;  
      
  
  //Search Deltastoshowlist and mark if it is present in Deltalist
  for(unsigned int j=0; j<Deltastoshowlist.size(); j++)
    {
      NumberList d=Deltastoshowlist[j];
      bool found=false;
      unsigned int indx=0;
      while( indx<Deltalist.size() && !found)
	{
	  if(TRACE) cout  << indx << " d=" << d << " Deltalist " << Deltalist[indx] << endl;
	  if(Deltalist[indx]==d){found=true;}
	  else{indx++;}
	}
      if(found){Printinfolist[indx]=true;}
    }
  
  if(TRACE)
    {
      cout << "Printinfolist" << endl;
      for(unsigned int i=0; i<Printinfolist.size(); i++)
	{
	  cout << i << " " << Deltalist[i] << " " << "info: " << Printinfolist[i] << endl;
	}
    }


  // Then seek initial value file for epsilons
  ifstream epsilonfile(EPSILONFILENAME.c_str());
  if(!epsilonfile)
    {
      // set default starting values
      NumberList newepsilon(NELASTIC);
      for(int i=0; i<NELASTIC; i++) newepsilon[i]=(RAN()-0.5)*par[EPSILONNULL]; // small values
      // put the same starting value for all entries:
      for(unsigned int i=0; i<Deltalist.size(); i++){epsilonlist.push_back(newepsilon);}
      
      if(TRACE)
	{
	  cout << "No file " << EPSILONFILENAME << " found." ;
	  cout << "Using epsilon=" << epsilon << endl;
	}
      
      if(TRACE)
	{
	  logfile << "No file " << EPSILONFILENAME << " found." ;
	  logfile << "Using epsilon=" << epsilon << endl;
	}
    }
  else
    {
      if(TRACE) cout << "Reading " << EPSILONFILENAME << " from disk" << endl;
      logfile << "Reading " << EPSILONFILENAME << " from disk" << endl;
      
      string line;
      while (getline(epsilonfile, line))
	{
	  istringstream iss(line);
	  realtype newval;
	  if(!(iss >> newval)){ break;}
	  
	  NumberList newepsilon(NELASTIC,newval); // initialize with one value for all
	  
	  for(int i=1; i<NELASTIC; i++)
	    {
	      if(!(iss >> newval)){ break;}
		newepsilon[i]=newval;
	    }
	  epsilonlist.push_back(newepsilon);
	}
      
      // duplicate the last entry until size is equal to Deltalist
      while(epsilonlist.size() < Deltalist.size())
	{
	  NumberList newepsilon=epsilonlist.back(); // repeat last element
	  epsilonlist.push_back(newepsilon);
	} 
    }
  if(TRACE) cout << "epsilonlist has " << epsilonlist.size() << " entries" << endl;
  logfile << "epsilonlist has " << epsilonlist.size() << " entries" << endl; 
  
  if(TRACE) cout << "Done Initializing Simulation" << endl;
}


void Simulation::Run()
{
  if(TRACE) cout << "Starting Run" << endl;
  for(unsigned int i=0; i<Deltalist.size(); i++)
    {
      mysolver.Solve(Deltalist[i],epsilonlist[i],Printinfolist[i]);
    }
  if(TRACE) cout << "Done Run" << endl;
}



#endif //BOND_H
