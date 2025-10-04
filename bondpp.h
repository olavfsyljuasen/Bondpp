#ifndef BOND_H
#define BOND_H

#include<vector>
#include<complex>
#include<iostream>
#include<fstream>
#include<algorithm>
#include<sstream>

#if defined LATTICEDISTORTIONS && defined PHONONS 
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


struct Qandvals
{
  Qandvals(Coord qin=Coord(0,0,0),realtype vin=0):q(qin),v(vin){}
  Coord q;
  realtype v;
};

bool operator<(const Qandvals& l,Qandvals r){return l.v > r.v;} // gives the sorted order, true if correct


struct IndxVal
{
  IndxVal(int in_indx=0,realtype in_val=0):indx(in_indx),val(in_val){}
  int indx;
  realtype val;
};

bool operator<(const IndxVal& l,IndxVal r){return l.val > r.val;} // gives the sorted order, true if correct



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

  void SaveState(string filename);
  bool LoadState(string filename);

  vector<IndxVal> FindMatrixMaxVals(VecMat<complextype,NMAT,NMAT>& A,const int nmax); // output the maxvalue and q indx of the nmax max values




  
  
  realtype CalculateT(int);
  void CalculateTs(vector<realtype>&);
  




#if defined LATTICEDISTORTIONS && defined ELASTIC
  NumberList CalculateEpsilonsOverT();
#endif

  realtype CalculateFreeEnergy(const realtype);
  realtype CalculatePhononFreeEnergy(const realtype);
  realtype CalculateElasticFreeEnergy(const realtype);
  vector<obstype> CalculateSpinOrderPars(realtype);
  vector<obstype> CalculateOrderPars(realtype,int,int);
  vector<obstype> CalculateAlphas(realtype);

  void FindMaxVals(SMatrix<int,NMAT,NMAT>&, SMatrix<realtype,NMAT,NMAT>&);

  //void Convolve(const bool);
  //  realtype CalculateSecondDerivative(const realtype T,const int k);

  void ConstructKinvq(const bool);
  void ComputeDq(const bool,const bool);  
  void ComputeSelfEnergy(const bool);
  
  void SetQsToZero(); // routine to set some q's to zero in self-energy
  void BiasSigma(); // routine to bias sigma by specifying Kinvq
  void MakeRandomSigma();

  void MakeSymmetric(VecMat<complextype,NMAT,NMAT>&);
  
  bool SolveSelfConsistentEquation(NumberList Delta,bool); 
  bool SolveSaddlePointEquations(realtype&,NumberList&);

  
  bool Solve(NumberList delta,NumberList thisepsilon,const bool pinfo,bool loadstate)
  {
    Delta = delta;
    RenormalizedDelta = delta; //  just as starting value
    
#if defined LATTICEDISTORTIONS && defined ELASTIC
    if( USEPREVIOUSEPSILONS && EpsilonInitialized)
      {
	// use the old epsilons if previous step converged, do not do anything
	logfile << "Use previous converged epsilons:" << epsilon << endl;
      }
    else
      {
	// set in new starting epsilon values
	epsilon=thisepsilon;
	EpsilonInitialized=true;
      }
#endif

    if(TRACE) cout << "Starting Solve with Delta= " << scientific << setprecision(LOGPRECISION) << delta << " epsilon=" << epsilon << " Printinfo= " << pinfo << endl;

    logfile << "Starting Solver with Delta= " << scientific << delta << " epsilon=" << epsilon << endl;

    
    Printinfo=pinfo;

    bool go_on=SolveSelfConsistentEquation(Delta,loadstate);
    if(TRACE) cout << "Done Solve" << endl;
    return go_on;
  }
  
  
 private:
  Rule& rule;

  int dim;              // the number of dimensions
  vector<int> dims;     // the size of each dimension
  
  const int Nq; // number of q-space sites.
  const realtype invNq;
  const realtype invSqrtNq;

  bool converged; //  
  bool pconverged;
  bool done;
  
  NumberList Delta; 
  NumberList RenormalizedDelta;

  vector<realtype> Ts; // holding the temperatures
  
#ifdef PRESERVESYMMETRY
  int TransformationPeriod;
  vector<int> TransformationTable;
#endif

  void WriteEpsilon(string filename);
  void WriteKinvq(string filename);
  void ReadKinvq(string filename);
  
  bool Printinfo; // set this to print out correlation functions
  int lineid; 
  bool SigmaInitialized; // flag to indicate whether Sigma is initialized.  
  bool EpsilonInitialized; // flag to indicate whether Epsilonlist is initialized.
  bool InitializeKinvqFromFile; // flag to indicate that one should try to initialize Kinvq from file if it exists.
  bool failed;           // flag to indicate failure to converge
  int MaxIterMultiplier; // multiplier for MaxIter, usually only used in the first step

  realtype mineigenvalue; // for storing the minimum SigmaE value
  //  realtype currT; // for storing the current value of the temperature
  realtype newT; // for storing the current value of the temperature 

  VecMat<complextype,NMAT,NMAT>& Jq;
  // The actual storage areas
  VecMat<complextype,NMAT,NMAT> A1;  // holds Kq,Kinvq
  VecMat<complextype,NMAT,NMAT> A2;  // holds Sigmaq,
  VecMat<complextype,NDMAT,NDMAT> B;  // holds Dq,Dinvq


  vector<complextype>  F1 ; // holds inplace intermediate Fourier-transform
  vector<complextype>  F2 ; // holds inplace intermediate Fourier-transform

#ifdef FFTS_INPLACE
  VecMat<complextype,NMAT,NMAT>& A1r; // holds Kinvr
  VecMat<complextype,NMAT,NMAT>& A2r; // holds Sigmar
  VecMat<complextype,NDMAT,NDMAT>& Br; // holds Dr
#else   
  VecMat<complextype,NMAT,NMAT> A1r;  // holds Kinvr
  VecMat<complextype,NMAT,NMAT> A2r;  // holds Sigmar
  VecMat<complextype,NDMAT,NDMAT> Br;  // holds Dr
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


  VecMat<complextype,NMAT,NMAT>& Kq;     // points to the A1 array
  VecMat<complextype,NMAT,NMAT>& Kinvq;  // points to the A1 array
  VecMat<complextype,NMAT,NMAT>& Kinvr;  // points to the A1 array

  VecMat<complextype,NMAT,NMAT>& Sigmar; // points to the A2 array
  VecMat<complextype,NMAT,NMAT>& Sigmaq; // points to the A2 array

  VecMat<complextype,NDMAT,NDMAT>& Dq;     // points to the B array
  VecMat<complextype,NDMAT,NDMAT>& Dinvq;  // points to the B array
  VecMat<complextype,NDMAT,NDMAT>& Dr;     // points to the B array

#ifdef LATTICEDISTORTIONS
  VecMat<complextype,NC,NMODE>& f; // holds the vertex information
  VecMat<complextype,NMAT,NMAT>& g; // points to rule
  
  vector<VecMat<complextype,NMAT,NMAT>* > gel;
#endif
  NumberList epsilon;  // amplitude of elastic modes


  
};


Driver::Driver(Rule& r):rule(r),dim(la.D()),dims(la.SiteqDims()),Nq(la.NqSites()),invNq(static_cast<realtype>(1.)/Nq),invSqrtNq(1./sqrt(Nq)),converged(false),pconverged(false),done(false),Delta(NSUBL),RenormalizedDelta(NSUBL),Ts(NSUBL),Printinfo(false),lineid(0),SigmaInitialized(false),EpsilonInitialized(false),InitializeKinvqFromFile(true),failed(false),MaxIterMultiplier(1),newT(-1.),Jq(r.Jq),  
  A1(Nq),A2(Nq),B(Nq),F1(Nq),F2(Nq),
#ifdef FFTS_INPLACE
  A1r(A1),A2r(A2),Br(B),
#else
  A1r(Nq),A2r(Nq),Br(Nq),
#endif
  Kq(A1),
  Kinvq(A1),Kinvr(A1r),Sigmar(A2r),Sigmaq(A2),
  Dq(B),Dinvq(B),Dr(Br)
#ifdef LATTICEDISTORTIONS
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


void Driver::SaveState(string filename)
{
  logfile << "SaveState to file " << filename << endl;
  
  if(TRACELEVEL>0) cout << spaces(ir++) << "Starting SaveState(" << filename << ")" << endl;

  ofstream outfile(filename.c_str());
  outfile.write((char*) &dim,sizeof(dim));
  outfile.write((char*) &dims[0],dims.size()*sizeof(dims[0]));
  outfile.write((char*) &Nq,sizeof(Nq));
  outfile.write((char*) &invNq,sizeof(invNq));
  outfile.write((char*) &invSqrtNq,sizeof(invSqrtNq));
  outfile.write((char*) &converged,sizeof(converged));
  outfile.write((char*) &pconverged,sizeof(pconverged));
  outfile.write((char*) &done,sizeof(done));
  outfile.write((char*) &Delta[0],Delta.size()*sizeof(Delta[0]));
  outfile.write((char*) &RenormalizedDelta[0],RenormalizedDelta.size()*sizeof(RenormalizedDelta[0]));
  outfile.write((char*) &Ts[0],Ts.size()*sizeof(Ts[0]));
#ifdef PRESERVESYMMETRY
  outfile.write((char*) &TransformationPeriod,sizeof(TransformationPeriod));
  outfile.write((char*) &TransformationTable[0],TransformationTable.size()*sizeof(TransformationTable[0]));
#endif
  outfile.write((char*) &Printinfo,sizeof(Printinfo));
  outfile.write((char*) &lineid,sizeof(lineid));
  outfile.write((char*) &SigmaInitialized,sizeof(SigmaInitialized));
  outfile.write((char*) &EpsilonInitialized,sizeof(EpsilonInitialized));
  outfile.write((char*) &InitializeKinvqFromFile,sizeof(InitializeKinvqFromFile));
  outfile.write((char*) &failed,sizeof(failed));
  outfile.write((char*) &MaxIterMultiplier,sizeof(MaxIterMultiplier));
  outfile.write((char*) &mineigenvalue,sizeof(mineigenvalue));
  outfile.write((char*) &newT,sizeof(newT));

  outfile.write((char*) A1.start(),A1.size()*sizeof(complextype));
  outfile.write((char*) A2.start(),A2.size()*sizeof(complextype));
  outfile.write((char*) B.start(),B.size()*sizeof(complextype));

  outfile.write((char*) &F1[0],F1.size()*sizeof(complextype));
  outfile.write((char*) &F2[0],F2.size()*sizeof(complextype));

#ifndef FFTS_INPLACE
  outfile.write((char*) A1r.start(),A1r.size()*sizeof(complextype));
  outfile.write((char*) A2r.start(),A2r.size()*sizeof(complextype));
  outfile.write((char*) Br.start(),Br.size()*sizeof(complextype));
#endif
  outfile.write((char*) &epsilon[0],epsilon.size()*sizeof(epsilon[0]));
  
  outfile.close();
  if(TRACELEVEL>0) cout << spaces(--ir) << "Done SaveState(" << filename << ")" << endl;
}



bool Driver::LoadState(string filename)
{  
  if(TRACELEVEL > 0) cout << spaces(ir++) << "Start LoadState(" << filename << ")" << endl;
  ifstream infile(filename.c_str());
  if(!infile)
    {
      if(TRACELEVEL > 0) cout << spaces(--ir) << "Done LoadState(" << filename << "), no file found" << endl;
      return false;
    }
  logfile << "LoadState from " << filename << endl;

  infile.read((char*) &dim,sizeof(dim));
  infile.read((char*) &dims[0],dims.size()*sizeof(dims[0]));
  infile.read((char*) &Nq,sizeof(Nq));
  infile.read((char*) &invNq,sizeof(invNq));
  infile.read((char*) &invSqrtNq,sizeof(invSqrtNq));
  infile.read((char*) &converged,sizeof(converged));
  infile.read((char*) &pconverged,sizeof(pconverged));
  infile.read((char*) &done,sizeof(done));
  infile.read((char*) &Delta[0],Delta.size()*sizeof(Delta[0]));
  infile.read((char*) &RenormalizedDelta[0],RenormalizedDelta.size()*sizeof(RenormalizedDelta[0]));
  infile.read((char*) &Ts[0],Ts.size()*sizeof(Ts[0]));
#ifdef PRESERVESYMMETRY
  infile.read((char*) &TransformationPeriod,sizeof(TransformationPeriod));
  infile.read((char*) &TransformationTable[0],TransformationTable.size()*sizeof(TransformationTable[0]));
#endif
  infile.read((char*) &Printinfo,sizeof(Printinfo));
  infile.read((char*) &lineid,sizeof(lineid));
  infile.read((char*) &SigmaInitialized,sizeof(SigmaInitialized));
  infile.read((char*) &EpsilonInitialized,sizeof(EpsilonInitialized));
  infile.read((char*) &InitializeKinvqFromFile,sizeof(InitializeKinvqFromFile));
  infile.read((char*) &failed,sizeof(failed));
  infile.read((char*) &MaxIterMultiplier,sizeof(MaxIterMultiplier));
  infile.read((char*) &mineigenvalue,sizeof(mineigenvalue));
  infile.read((char*) &newT,sizeof(newT));

  infile.read((char*) A1.start(),A1.size()*sizeof(complextype));
  infile.read((char*) A2.start(),A2.size()*sizeof(complextype));
  infile.read((char*) B.start(),B.size()*sizeof(complextype));

  infile.read((char*) &F1[0],F1.size()*sizeof(complextype));
  infile.read((char*) &F2[0],F2.size()*sizeof(complextype));

#ifndef FFTS_INPLACE
  infile.read((char*) A1r.start(),A1r.size()*sizeof(complextype));
  infile.read((char*) A2r.start(),A2r.size()*sizeof(complextype));
  infile.read((char*) Br.start(),Br.size()*sizeof(complextype));
#endif
  infile.read((char*) &epsilon[0],epsilon.size()*sizeof(epsilon[0]));

  
  infile.close();
  if(TRACELEVEL>0) cout << spaces(--ir) << "Done LoadState(" << filename << ")" << endl;
  return true;
}

  
// find the nmax largest values of the matrix mymatrix and output its q-indx
vector<IndxVal> Driver::FindMatrixMaxVals(VecMat<complextype,NMAT,NMAT>& mymatrix,const int nmax)
{
      // look for max among all matrix elements in first brillouinzone
      vector<IndxVal> maxlist(nmax); // sorted list of qs and their values, biggest first
      
      for(int qi=0; qi<la.NqSites(); qi++) 
	{
	  realtype matrixmaxval(0.);

	  for(int s1=0; s1<NSUBL; s1++)
	    for(int s2=0; s2<NSUBL; s2++)
	      {
		realtype val = abs(mymatrix(qi,s1,s2));
		if(val > matrixmaxval){ matrixmaxval=val;}
	      }
	  
	  if(matrixmaxval > maxlist[nmax-1].val) // compare with last element, insert if true
	    {
	      IndxVal newelem(qi,matrixmaxval);
	      vector<IndxVal>::iterator iter; // an iterator to the insertion point
	      iter=upper_bound(maxlist.begin(),maxlist.begin()+nmax,newelem);
	      maxlist.insert(iter,newelem); // insert and keep it sorted
	      maxlist.resize(nmax); 
	    }
	}
      return maxlist;
}


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


#if defined LATTICEDISTORTIONS && defined ELASTIC
NumberList Driver::CalculateEpsilonsOverT()
{
  if(TRACE) cout << "Starting CalculateEpsilonsOverT()" << endl;

  NumberList epsoverT(NELASTIC);
  for(int i=0; i<NELASTIC; i++)
    {
      realtype mui= rule.elasticeigenvalues[i];
      realtype sum=0.;
      if(mui != 0.)
	{
	  for(int qi=0; qi<Nq; qi++)
	    {
	      SMatrix<complextype,NMAT,NMAT> tmp;
	      tmp=Kinvq[qi];
	      tmp *= (*rule.gelptrs[i])[qi];
	      sum += NFAKESPINTRACE*real(tr(tmp));
	    }
	  //	  sum *= -1./(2.*mui*Nq*la.UnitCellVolume());
	  sum *= -1./(2.*mui*Nq); // define elastic moduli tensor per site instead of per unit volume, then no UnitCellVolume is necessary.
	}
      epsoverT[i]= sum;
    }
 
  if(TRACE) cout << "Done CalculateEpsilonsOverT "  << endl;
  return epsoverT;
}
#endif


void Driver::WriteEpsilon(string filename)
{
  ofstream outfile(filename);
  outfile << scientific << setprecision(OUTPUTPRECISION);
  outfile << setw(OUTPUTPRECISION+8) << epsilon << endl;
  outfile.close();
}

void Driver::WriteKinvq(string filename)
{
  ofstream outfile(filename);
  outfile << scientific << setprecision(OUTPUTPRECISION);
  for(int q=0; q<la.NqSites(); q++)
    {
      for(int i=0; i<NMAT; i++)
	for(int j=0; j<NMAT; j++)
	  outfile << setw(OUTPUTPRECISION+8) << Kinvq(q,i,j) << " ";
      outfile << endl;
    }
  outfile.close();
}

void Driver::ReadKinvq(string filename)
{
  ifstream infile(filename);
  if(infile)
    {
      logfile << "Reading in Kinvq from " << filename << endl;
      for(int q=0; q<la.NqSites(); q++)
	{
	  for(int i=0; i<NMAT; i++)
	    for(int j=0; j<NMAT; j++)
	      infile >> Kinvq(q,i,j);
	}
      infile.close();
    }
}


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
	    complextype val=(m1==m2 ? Kinvq(qi,m1,m2): 0.5*(Kinvq(qi,m1,m2)+Kinvq(qi,m2,m1)));
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
void Driver::MakeSymmetric(VecMat<complextype,NMAT,NMAT>& m)
{
  if(TRACE) cout << "Starting MakeSymmetric" << endl;
  if( NSUBL !=1){ cout << "MakeSymmetric works only for NSUBL=1" << endl; exit(1);}
  
  complextype* mstart=m.start(); // This only works for NSUBL=1 
  
  //  cout << "m=" << m << endl;

  VecMat<complextype> temp(m);
  complextype* tempstart=temp.start(); // This only works for NSUBL=1 
  
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
  realtype betaf_constants = realtype(-0.5*NSUBL)*(mylog(Nq)+(Ns-2)*mylog(PI)+mylog(TWOPI));

  if(TRACEFREEENERGY) cout << scientific << setprecision(6) << "T=" << T << " bf_const=" << betaf_constants;
    
  if(TRACE) cout << "betaf_constants  = " << betaf_constants << endl;
  
  f += T*betaf_constants;

  //temperature factors
  realtype betaf_Tdep = realtype(0.5*NSUBL*(Ns-2))*mylog(1./T); 
  if(TRACE) cout << "betaf_Tdep  = " << betaf_Tdep << endl;
  
  f += T*betaf_Tdep;

  if(TRACEFREEENERGY) cout << " bf_Tdep=" << betaf_Tdep;

#if defined LATTICEDISTORTIONS

#if defined ELASTIC
  //#if !(defined OMITELASTICCONT || defined NOELASTIC)
  //elastic modes  
  realtype f_elastic = 0;
  for(int i=0; i<NELASTIC; i++){ f_elastic += realtype(0.5)*epsilon[i]*epsilon[i]*rule.elasticeigenvalues[i];}

  if(TRACE) cout << "f_elastic  = " << f_elastic << endl;
  if(TRACEFREEENERGY) cout << " f_elastic=" << f_elastic;
  f += f_elastic;
#endif

  //  #if !defined(ELASTICONLY)  || !defined(OMITPHONONCONT)
#ifndef OMITPHONONCONT
  
  // the following term is subtracted from D, but added here, so as to count phonon
  // contributions to the free energy where it is appropriate. Remember,
  // the Dinv in this code contains a term beta on the phonon diagonal part, also
  // when the coupling to phonons g=0. Thus in order to count this as a contribution to
  // the phonons part of the free energy, we subtract it from D and add it here.    
  realtype betaf_phononTdep = realtype(-0.5)*NMODE*mylog(T); 
  if(TRACE) cout << "betaf_phononTdep  = " << betaf_phononTdep << endl;

  if(TRACEFREEENERGY) cout << " bf_phononTdep=" << betaf_phononTdep;
  
  f += T*betaf_phononTdep;

  // new in v1.95: in addition there is a similar term from P^2 as that
  // is alos multiplied by beta

  realtype betaf_phononTdep_fromP2 = realtype(-0.5)*NMODE*mylog(T); 
  if(TRACE) cout << "betaf_phononTdep_fromP2  = " << betaf_phononTdep_fromP2 << endl;

  if(TRACEFREEENERGY) cout << " bf_phononTdep_fromP2=" << betaf_phononTdep_fromP2;
  
  f += T*betaf_phononTdep_fromP2;
  
  
  // constants:
  realtype betaf_phononconst = realtype(-0.5)*2*NMODE*mylog(TWOPI); // 2 because both X^2 and P^2
  if(TRACE) cout << "betaf_phononconst  = " << betaf_phononconst << endl;
  if(TRACEFREEENERGY) cout << " bf_phononconst=" << betaf_phononconst;
  
  f += T*betaf_phononconst;

  //the phonon spectrum
  //The following line has been replaced from version 1.95 as only X^2 is multiplied by omega^2:
  //realtype betaf_phonons = 2*rule.GetSumLogOmegaoverV(); // 2 both X^2 and P^2
  realtype betaf_phonons = invNq*rule.GetSumLogOmega(); // from v1.95: only X^2

  if(TRACE) cout << "betaf_phonons  = " << betaf_phonons << endl;
  if(TRACEFREEENERGY) cout << " bf_phonons=" << betaf_phonons;
  
  f += T*betaf_phonons;

#endif  // not OMITPHONONCONT
  
#endif  // lattice distortions
  
  //Must correct the Delta values for the subtraction of the minimum from SigmaE
  //  for(int i=0; i<NSUBL; i++) f += -(Delta[i]-mineigenvalue);
  realtype f_delta(0.);
  for(int i=0; i<NSUBL; i++) f_delta += -RenormalizedDelta[i];

  if(TRACE) cout << "betaf_delta      = " << f_delta/T << endl;
  if(TRACEFREEENERGY) cout << " f_delta=" << f_delta;

  f += f_delta;
  
  realtype betaf_logKinvq   = realtype(-0.5*invNq)*NFAKESPINTRACE*SumLogDet(Kinvq);

  if(TRACE) cout << "betaf_logKinvq   =  " << betaf_logKinvq << endl;
  if(TRACEFREEENERGY) cout << " betaf_logKinvq=" << betaf_logKinvq;
  f += T*betaf_logKinvq;
  
  ComputeDq(false,true); // excludeqzero=false, preserveinput=true not to jeopardize Keff
  
  realtype betaf_logD       = realtype(-0.5*invNq)*SumLogDet(Dq);
  if(TRACEFREEENERGY) cout << " betaf_logD1=" << betaf_logD;
  //#if defined LATTICEDISTORTIONS && !defined ELASTICONLY
#if defined LATTICEDISTORTIONS && defined PHONONS 
  betaf_logD += realtype(0.5)*NMODE*mylog(T); //  -0.5*NMODE*log(1/T) 
#endif
  if(TRACE) cout << "betaf_logD       =  " << betaf_logD << endl;
  if(TRACEFREEENERGY) cout << " betaf_logD2=" << betaf_logD;
  
  f += T*betaf_logD;

  ComputeSelfEnergy(true); // ensure that Kinvq is the same.

  realtype betaf_KinvqSigma = -realtype(0.5*invNq)*NFAKESPINTRACE*SumTr(Kinvq,Sigmaq);
  if(TRACE) cout << "betaf_KinvqSigma = " << betaf_KinvqSigma << endl;
  if(TRACEFREEENERGY) cout << " betaf_KinvqSigma=" << betaf_KinvqSigma;
  f += T*betaf_KinvqSigma;

  if(TRACEFREEENERGY) cout << " f=" << f << endl;
  // the correction to the saddle-point are already taken into account
  return f;
}


// Phonon Free energy per unit volume
// we use the convention \nu=\nup
realtype Driver::CalculatePhononFreeEnergy(realtype T)
{
  if(TRACE) cout << "Starting CalculatePhononFreeEnergy " << endl;

  realtype f=0;

  realtype zeromodeomissionfactor= (Nq-1)*invNq;  // corrects for the missing phonon q=0 modes, which are elastic modes 


  //the bare phonon spectrum, from scaled out Xs
  realtype betaf_phonons = invNq*rule.GetSumLogOmega(); // only X^2 
  
  if(TRACE) cout << "betaf_phonons  = " << betaf_phonons << endl;
  if(TRACEFREEENERGY) cout << " bf_phonons=" << betaf_phonons;
  
  f += T*betaf_phonons;
  
  realtype betaf_phononTdep = -realtype(0.5)*zeromodeomissionfactor*NMODE*mylog(T); // from P^2, the X^2 contribution beta is in Dinvq 
  if(TRACE) cout << "betaf_phononTdep  = " << betaf_phononTdep << endl;

  if(TRACEFREEENERGY) cout << " bf_phononTdep=" << betaf_phononTdep;
  
  f += T*betaf_phononTdep;

  // constants:
  realtype betaf_phononconst = -realtype(0.5)*2*zeromodeomissionfactor*NMODE*mylog(TWOPI); // 2 as it comes from both X^2 and P^2 Gaussian integrals
  if(TRACE) cout << "betaf_phononconst  = " << betaf_phononconst << endl;
  if(TRACEFREEENERGY) cout << " bf_phononconst=" << betaf_phononconst;
  
  f += T*betaf_phononconst;
  
  ComputeDq(false,true); // excludeqzero=false, preserveinput=true not to jeopardize Keff

  MatrixInverse(Dq); // B=Dinvq
  realtype betaf_logD = realtype(0.5*invNq)*SumLogDet(Dinvq,1,true); // 1: only count phonons, excludeqzero=true
  MatrixInverse(Dinvq); // B=Dq, probably not necessary
    
  if(TRACE) cout << "betaf_logD       =  " << betaf_logD << endl;
  if(TRACEFREEENERGY) cout << " betaf_logD=" << betaf_logD;
  
  f += T*betaf_logD;

  if(TRACEFREEENERGY) cout << " f_phonons=" << f << endl;

  return f;
}

// Phonon Free energy per unit volume
// we use the convention \nu=\nup

realtype Driver::CalculateElasticFreeEnergy(realtype T)
{
  if(TRACE) cout << "Starting CalculateElasticFreeEnergy " << endl;

  realtype f=0;

  realtype f_elastic = 0;
  for(int i=0; i<NELASTIC; i++){ f_elastic += realtype(0.5)*epsilon[i]*epsilon[i]*rule.elasticeigenvalues[i];}

  if(TRACE) cout << "f_elastic  = " << f_elastic << endl;
  if(TRACEFREEENERGY) cout << " f_elastic=" << f_elastic;
  f += f_elastic;

  if(TRACEFREEENERGY) cout << " f=" << f << endl;

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
#ifdef LATTICEDISTORTIONS
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

  if(TRACELEVEL>=4) cout << "Kinvr=" << Kinvr << endl; 
  
  // Initialize Dinvq:
  Dinvq.SetToZero();

  // first compute the constraint block (this implementation is valid also for NSUBL >1)

  complextype tmp(0.,0.);
  
  for(int m1=0; m1<NSUBL; m1++)
    for(int m2=m1; m2<NSUBL; m2++)
      {
	for(int r=0; r<Nq; r++)
	  {
	    complextype tt(0.,0.);
	    for(int s1=0; s1<NSPIN; s1++)
	      for(int s2=0; s2<NSPIN; s2++)
		{
		  int alpha=mindx(s1,m2); // alpha_s = m2
		  int delta=mindx(s2,m1); // delta_s = m1
		  
		  tmp=Kinvr(r,alpha,delta);
		  tt += tmp*tmp;   
		}
	    F1[r]=complextype(0.5*NFAKESPINTRACE,0.)*tt;
	  }

	FFTWEXECUTE(F1pluss); // 
	
	for(int qi=0; qi<Nq; qi++)
	  {
	    Dinvq(qi,m1,m2)=F1[qi];
	    if(m2 != m1) Dinvq(qi,m2,m1)=conj(F1[qi]); 
	  }
      }

  if(TRACELEVEL>=4)
    {
      cout << "after constraint block, Dinvq: " << endl;
      cout << Dinvq << endl;
    }

#if defined LATTICEDISTORTIONS && defined PHONONS
  // the offdiagonal phonon-constraint part
#ifdef CPOSITIVE
  realtype multiplier=realtype(2.);
#else
  realtype multiplier=realtype(1.);
#endif  
    
  for(int m1=0; m1<NSUBL; m1++) // m1 must be 0 here because phonons NSUBL=1
    for(unsigned int cc=0; cc<nonzeroclist.size(); cc++)
      {
	int c = nonzeroclist[cc]; // loop thru only c values for which gc is non-zero
	for(int r=0; r<Nq; r++)
	  {
	    SMatrix<complextype,NMAT,NMAT> my;
	    my  = g[c];
	    my *= Kinvr[r];
	    int mrpc = la.rAdd(la.GetInversionIndx(r),clist[c]); // index of -r+c	    
	    my *= Kinvr[mrpc];
	    
	    F1[r] = complextype(NFAKESPINTRACE,0.)*tr(my);  
	  }
	
	// do fourier-transform of theta
	FFTWEXECUTE(F1pluss); 
	
	for(int qi=0; qi<Nq; qi++)
	  {
	    for(int m2=NSUBL; m2<NDMAT; m2++)
	      {
		int n2=m2-NSUBL; 
		
		//	cout << "q=" << qi << " f=" << f(qi,c,n2) << endl;
		Dinvq(qi,m1,m2)+= complextype(multiplier*realtype(0.5)*invSqrtNq,0.)*F1[qi]*f(qi,c,n2)*conj(expi(la.qr(qi,clist[c])));
	      }		   
	  }
      }
  
  // Set the constraint-phonon part.
  for(int qi=0; qi<Nq; qi++)
    for(int m1=0; m1<NSUBL; m1++)
      for(int m2=NSUBL; m2<NDMAT; m2++)
	{
	  //	  Dinvq(la.GetInversionIndx(qi),m2,m1)=Dinvq(qi,m1,m2);
	  Dinvq(qi,m2,m1)=-conj(Dinvq(qi,m1,m2));  // should be the same as the above, but probably faster, changed from v1.67
	}

  if(TRACELEVEL>=4)
    {
      cout << "after constraint-phonon block, Dinvq: " << endl;
      cout << Dinvq << endl;
    }


  // finally the phonon-phonon part
#ifdef CPOSITIVE
  /*
    for(unsigned int cc2=0; cc2< nonzeroclist.size(); cc2++)
      for(unsigned int cc4=0; cc4< nonzeroclist.size(); cc4++)
	{	    
	  int c2 = nonzeroclist[cc2];
	  int c4 = nonzeroclist[cc4];

	  if(TRACE2)
	    cout << "EXPERIMENTAL code: c2=" << c2 << " c4=" << c4 << endl;
	  
	  for(int r=0; r<Nq; r++)
	    {
	      SMatrix<complextype> my(NMAT,NMAT);
	      
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
		    const complextype myF1 =( qi> la.qlimit ? F1[la.GetInversionIndx(qi)].conj():F1[qi]);
		    
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
	      SMatrix<complextype> my(NMAT,NMAT);
	      
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
  for(unsigned int cc2=0; cc2<nonzeroclist.size(); cc2++)
    for(unsigned int cc4=0; cc4<nonzeroclist.size(); cc4++)
      {
	int c2 = nonzeroclist[cc2];
	int c4 = nonzeroclist[cc4];
	
	  if(TRACE2)
	    cout << "c2=" << c2 << " c4=" << c4 << endl;

	  	  	  
	  Triplet myivec=clist[c2]+clist[c4];
	  
	  for(int r=0; r<Nq; r++)
	    {
	      SMatrix<complextype,NMAT,NMAT> my;
	      
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
	      SMatrix<complextype,NMAT,NMAT> my;
	      
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
		    Dinvq(qi,m1,m2)+=-complextype(invNq,0.)*F1[qi]*conj(f(qi,c2,n1))*f(qi,c4,n2);
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
  for(unsigned int cc2=0; cc2<nonzeroclist.size(); cc2++)
    for(unsigned int cc4=0; cc4<nonzeroclist.size(); cc4++)
      {
	int c2 = nonzeroclist[cc2];
	int c4 = nonzeroclist[cc4];

	Triplet myivec=clist[c2]+clist[c4];
	
	for(int i=0; i<Nq; i++)
	  {
	    SMatrix<complextype> my(NMAT,NMAT);
	    
	    my  = g[c4];
	    my *= Kinvr[i];
	    my *= g[c2];
	    
	    int newindx=la.rAdd(la.GetInversionIndx(i),myivec);
	    
	      my *= Kinvr[newindx];
	      
	      F1[i]=complextype(0.5*NFAKESPINTRACE,0.)*tr(my);
	  } 
	FFTWEXECUTE(F1pluss); 
	
	
	for(int m1=NSUBL; m1<NDMAT; m1++) 
	  for(int m2=m1; m2<NDMAT; m2++) 
	    {
	      
	      int n1=m1-NSUBL; 
	      int n2=m2-NSUBL; 
	      
	      for(int qi=0; qi<Nq; qi++)
		{
		  Dinvq(qi,m1,m2)+=complextype(-invNq,0.)*F1[qi]*conj(f(qi,c2,n1))*f(qi,c4,n2)*conj(expi(la.qr(qi,clist[c4])));
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
  
  //  realtype barepart=1./currT; // the one-half is included in the def. of propagator
  complextype barepart(1./newT,0.); // the one-half is included in the def. of propagator 
  for(int qi=0; qi<Nq; qi++)
    {
      for(int n=0; n<NMODE; n++)
	{
	  int m=n+NSUBL;
	  Dinvq(qi,m,m) += barepart;
	}
    }
#endif

  if(TRACELEVEL>=4) 
    {
      cout << "Dinvq after adding bare phonon part: " << Dinvq << endl;
    }
  
  // FORCING
  MakeMixedHermitian(Dinvq,NSUBL,NSUBL);
  MakeInversionTransposedSymmetric(Dinvq);
  

  Chomp(Dinvq);
  
  if(TRACELEVEL>0){SanityCheck(Dinvq,"Dinvq, after constructing it",false);}
  
  if(TRACELEVEL>=4){ cout << "Dinvq before inversion=" << Dinvq << endl;}

  MatrixInverse(Dinvq); // B = Dq

  if(TRACELEVEL>0){SanityCheck(Dq,"Dq, after inverting Dinvq",false);}
  
  if(TRACELEVEL>=4){ cout << "Dq, after inverting Dinvq, Dq=" << Dq << endl;}

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

  if(TRACELEVEL>0){SanityCheck(Dq,"Dq, at end of ComputeDq",false);}

  if(TRACELEVEL>=4){ cout << "At end of ComputDq, Dq=" << Dq << endl;}
  
  if(TRACELEVEL>0) cout << "Done with ComputeDq " << endl;
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

	for(int k=0; k<Nq; k++){Sigmaq(k,alpha,delta) += complextype(invSqrtNq,0.)*F2[k];}

	//	cout << "Sigmaq after term 1: alpha=" << alpha << " delta=" << delta << endl;
	//cout << Sigmaq << endl;

	
#if defined LATTICEDISTORTIONS && defined PHONONS
	if(NSUBL != 1){cout << "Error: LATTICEDISTORTIONS are not implemented for NSUBL!=1, exiting" << endl; exit(1);}


#ifdef CPOSITIVE
	// Term 2
	for(unsigned int cc=0; cc<nonzeroclist.size(); cc++)
	  {
	    int c=nonzeroclist[cc];
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
	    
	    for(int k=0; k<Nq; k++){ Sigmaq(k,alpha,delta) += complextype(invNq,0.)*F2[k];}

	  }

	// Term 3 (as term 2, but switch alpha and delta and take complex conj)
	for(unsigned int cc=0; cc<nonzeroclist.size(); cc++)
	  {
	    int c=nonzeroclist[cc];
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
	    
	    for(int k=0; k<Nq; k++){ Sigmaq(k,alpha,delta) += complextype(invNq,0.)*conj(expi(la.qr(k,clist[c]))*F2[k]);}

	    for(int r=0; r<Nq; r++)
	      {
		F2[r]=0;
		int rpc=la.rAdd(r,clist[c]);
		for(int s=0; s<NSPIN; s++)
		  F2[r]+=F1[r]*Kinvr(rpc,s,delta)*g(c,alpha,s); 
	      }
	    
	    FFTWEXECUTE(F2pluss);
	    
	    for(int k=0; k<Nq; k++){ Sigmaq(k,alpha,delta) += complextype(invNq,0.)*conj(F2[k]);}
	  }

	// Term4:
	for(unsigned int cc1=0; cc1<nonzeroclist.size(); cc1++)
	  for(unsigned int cc2=0; cc2<nonzeroclist.size(); cc2++)
	    {
	      int c1 = nonzeroclist[cc1];
	      int c2 = nonzeroclist[cc2];
	      
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
	      for(int k=0; k<Nq; k++){Sigmaq(k,alpha,delta)+= complextype(-invNq,0.)*invSqrtNq*F2[k]*expi(la.qr(k,clist[c1]))*expi(la.qr(k,clist[c2]));}

	      
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
	      for(int k=0; k<Nq; k++){Sigmaq(k,alpha,delta)+= complextype(-invNq,0.)*invSqrtNq*F2[k]*expi(la.qr(k,clist[c1]));}


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
	      for(int k=0; k<Nq; k++){Sigmaq(k,alpha,delta)+= complextype(-invNq,0.)*invSqrtNq*F2[k]*expi(la.qr(k,clist[c2]));}

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
	      for(int k=0; k<Nq; k++){Sigmaq(k,alpha,delta)+= complextype(-invNq,0.)*invSqrtNq*F2[k];}
	    }
#else	
	// Term 2
	for(unsigned int cc=0; cc< nonzeroclist.size(); cc++)
	  {
	    int c = nonzeroclist[cc];
	    
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
	    
	    for(int k=0; k<Nq; k++){ Sigmaq(k,alpha,delta) += complextype(invNq,0.)*expi(la.qr(k,clist[c]))*F2[k];}
	  }

	// Term 3 (as term 2, but switch alpha and delta and take complex conj)
	for(unsigned int cc=0; cc<nonzeroclist.size(); cc++)
	  {
	    int c = nonzeroclist[cc];
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
	    
	    for(int k=0; k<Nq; k++){ Sigmaq(k,alpha,delta) += complextype(invNq,0.)*conj(expi(la.qr(k,clist[c]))*F2[k]);}
	  }
	
	// Term 4:	
	for(unsigned int cc1=0; cc1 < nonzeroclist.size(); cc1++)
	  for(unsigned int cc2=0; cc2 < nonzeroclist.size(); cc2++)
	    {
	      int c1 = nonzeroclist[cc1];
	      int c2 = nonzeroclist[cc2];
	      
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
	      for(int k=0; k<Nq; k++){Sigmaq(k,alpha,delta)+= complextype(-invNq*invSqrtNq,0.)*F2[k]*expi(la.qr(k,clist[c1]))*expi(la.qr(k,clist[c2]));}
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
      Kinvq *= complextype(invSqrtNq,0.);

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
void Driver::ConstructKinvq(bool firsttime=false)
{
  if(TRACE) cout << "Starting ConstructKinvq" << endl;
  
  Kq =  Jq;

  if(TRACE) cout << "Initializing Kq" << endl;

  if(TRACE){SanityCheck(Kq,"Kq, after setting Jq");}

  Kq += Sigmaq;

  if(TRACE){SanityCheck(Kq,"Kq, after adding Sigmaq");}

  
#if defined LATTICEDISTORTIONS && defined ELASTIC
  if(TRACE) cout << "Adding elastic modes " << endl;
  for(int j=0; j<NELASTIC; j++)
    {
      VecMat<complextype,NMAT,NMAT> temp(*rule.gelptrs[j]);
      temp *= epsilon[j];
      Kq += temp;
    }
#endif

 
  if(TRACE){SanityCheck(Kq,"Kq, after adding Elastic modes");}
  
  mineigenvalue=SubtractMinimumEigenvalue(Kq);

  for(int i=0; i<NSUBL; i++){ RenormalizedDelta[i] = Delta[i]-mineigenvalue;}
  
  if(TRACE){SanityCheck(Kq,"Kq, after subtracting eigenvalues");}

  AddDelta(Kq,Delta);

  if(TRACE){SanityCheck(Kq,"Kq, after adding Delta");}
  
  // construct Kinvq
  MatrixInverse(Kq); // K = Kinv_q 

  // Force correct properties on the matrix
  MakeHermitian(Kinvq); 
  MakeInversionTransposedSymmetric(Kinvq); 
  MakeSpinDiagonal(Kinvq);

  if(firsttime)
    {
      if(InitializeKinvqFromFile){ ReadKinvq(INPUTKINVQFILENAME); }
      
      InitializeKinvqFromFile=false;

      
      CalculateTs(Ts);

      //      currT=Ts[0];
      newT=Ts[0];
      
      logfile << "initial Deltas and Ts:";
      logfile << scientific << setprecision(LOGPRECISION);
      for(int i=0; i<NSUBL; i++)
	{
	  logfile << setw(LOGPRECISION+8) << Delta[i] << setw(LOGPRECISION+8) << Ts[i] << "  ";
	}
      logfile << endl;
      //      logfile << scientific << setprecision(0) << setw(7) << FindSpread(Delta) << setw(7) << FindSpread(Ts) << endl;
    }
  
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
	  complextype c=da*de*complextype(RAN(),0.); // real positive value 
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
		complextype c=da*de*complextype(RAN(),RAN());
		Sigmaq(i,m1,m2)=c;
	      }
	  }
    }

  MakeHermitian(Sigmaq,false); // no warnings

  MakeSpinDiagonal(Sigmaq,false);

#ifdef FORCEINVERSIONSYMMETRY
  MakeInversionTransposedSymmetric(Sigmaq);
#endif

#ifdef PRESERVESYMMETRY
  MakeSymmetric(Sigmaq);
#endif

  
  if(TRACE) SanityCheck(Sigmaq,"Sigmaq, at end of MakeRandomSigma");
  if(TRACE) cout << "Finished MakeRandomSigma()" << endl;
}



void Driver::SetQsToZero()
{
  if(TRACE) cout << "Starting SetQsToZero()" << endl;
  logfile << "Setting the following qpts to zero in the self-energy" << endl;

  bool found=false;
  ifstream ifile("qstozero.in");
  
  while(ifile)
    {
      realtype qx_val;
      realtype qy_val;
      realtype qz_val;
      ifile >> qx_val;
      ifile >> qy_val;
      ifile >> qz_val;
      Coord qval(qx_val,qy_val,qz_val);

      if(ifile)
	{
	  Coord fqval=la.TranslateToFundamentalDomain(qval);
	  int bestqindx=la.FindMatchingqindx(fqval);

	  found=true;
	  Coord matchingq=la.qPos(bestqindx);
	  logfile << bestqindx << " (" << matchingq << ")" << endl;
	  for(int i=0; i<NMAT; i++)
	    for(int j=0; j<NMAT; j++)
	      Sigmaq(bestqindx,i,j)=complextype(0.,0.);
	}        
    }
  if(!found) logfile << "None" << endl;
  if(TRACE) cout << "Finished SetQsToZero()" << endl;
}


// a routine to bias Sigma so as to select a particular Keff 
void Driver::BiasSigma()
{
  // read in q's to bias. and which matrix elements to bias
  if(TRACE) cout << "Starting BiasSigma()" << endl;
  logfile << "Biasing the self-energy.";
  
  ifstream ifile("qbias.in");

  if(!ifile)
    {
      logfile << " No bias file found. Returning." << endl;
      if(TRACE) cout << "No bias file found. Finishing BiasSigma()" << endl;
      return;
    }
  logfile << endl;

  
  ConstructKinvq(true);

  
  while(ifile)
    {
      realtype qx_val;
      realtype qy_val;
      realtype qz_val;
      ifile >> qx_val;
      ifile >> qy_val;
      ifile >> qz_val;
      Coord qval(qx_val,qy_val,qz_val);

      
      if(ifile)
	{
	  Coord fqval=la.TranslateToFundamentalDomain(qval);
	  int bestqindx=la.FindMatchingqindx(fqval);
	  //	  qslots.push_back(bestqindx); // so that one can monitor a q slot
	  Coord matchingq=la.qPos(bestqindx);
	  logfile << "Biasing q=(" << qval << ")" << " q#=" <<  bestqindx << " (" << matchingq << ")" << endl;
	  if(TRACE) cout << "Biasing q=(" << qval << ") funddomainq=(" << fqval << ") q#=" <<  bestqindx << " (" << matchingq << ")" << endl;
	  for(int i=0; i<NMAT; i++)
	    for(int j=0; j<NMAT; j++)
	      {
		realtype realpart; ifile >> realpart;
		realtype imagpart; ifile >> imagpart;
		
		Kinvq(bestqindx,i,j)=complextype(realpart,imagpart);
	      }
	}         
    }

  MakeHermitian(Kinvq); // force the inverse to be Hermitian

  if(TRACELEVEL >2 )
    {
      cout << "Kinvq=" << endl;
      cout << Kinvq << endl;
    }
  
  ComputeSelfEnergy(); // this is then the new self-energy to use

  if(TRACELEVEL >2 )
    {
      cout << "Sigmaq:" << endl;
      cout << Sigmaq << endl;
    }
  
  if(TRACE) cout << "Finished BiasSigma()" << endl;
}





/*
void Driver::SetQsToZero()
{
  logfile << "Setting the following qpts entries to zero in the self-energy:" << endl;

  bool found = false;
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
	  for(int i1=0; i1<NMAT; i1++)
	    for(int i2=0; i2<NMAT; i2++)
	      Sigmaq(qindx,i1,i2)=complextype(0.,0.);
	}        
    }
  if(!found){logfile << "None" << endl;}
}
*/

bool Driver::SolveSaddlePointEquations(realtype& thisT,NumberList& thisepsilon)
{
  if(TRACE) cout << "Starting SolveSaddlePointEquations(" << thisT << " " << thisepsilon << ") " << endl;

  realtype myoldT(thisT);
  NumberList myoldeps(thisepsilon);

  CalculateTs(Ts);
  realtype myT=Ts[0];  
  NumberList myeps(thisepsilon); // initial value
  
#if defined LATTICEDISTORTIONS && defined ELASTIC
  NumberList epsoverT=CalculateEpsilonsOverT();
  myeps=epsoverT*myT;
#if defined CLAMPEDMODES
  for(int k=0; k<Nclampedmodes; k++){ myeps[clampedmodes[k]]=realtype(0.);}
#endif      
#endif
  int counter = 0;
  bool TooManyIter=false;
  NumberList epsdiff = myeps - myoldeps;
  while( (abs(myT - myoldT) > par[TOLERANCE] || maxabs(epsdiff) > par[TOLERANCE] ) && !TooManyIter)
    {
      if(TRACE)
	{
	  cout << "begin: myT:" << myT << " eps:" << myeps << endl;
	  //	  cout << "begin: myT:" << myT << " myoldT:" << myoldT << " epsdiff:" << epsdiff << " TooManyIter:" << TooManyIter << endl;
	}
      
      myoldT   = myT;
      myoldeps = myeps;
      
      
      epsilon = myeps; // set Kinvq using the new value of epsilon
      ConstructKinvq(); // also subtract minimal value and renormalize Delta

      // Extra lines here for better accuracy ???
      ComputeSelfEnergy();
      ConstructKinvq(); // also subtract minimal value and renormalize Delta
      // end extra lines

      myT=CalculateT(0);  // only valid for one sublattice here
#if defined LATTICEDISTORTIONS && defined ELASTIC
      epsoverT=CalculateEpsilonsOverT();
      myeps=epsoverT*myT;
#if defined CLAMPEDMODES
      for(int k=0; k<Nclampedmodes; k++){ myeps[clampedmodes[k]]= realtype(0.);}
#endif      
#endif
      epsdiff = myeps - myoldeps;

      
      /*
      if(counter > 10) // write log output
	{
	  logfile << "T and eps: " << myT << " " << myeps << endl;
	}
      */
      TooManyIter = ( ++counter >= 10000);

      if(TRACE)
	{
	  //	  cout << "end: myT:" << myT << " myoldT:" << myoldT << " eps:" << myeps << " epsdiff:" << epsdiff << " TooManyIter:" << TooManyIter << endl;
	  cout << "end: myT:" << myT << " eps:" << myeps << endl;
	}
    }
  
  
  if(counter > 100){ logfile << "used " << counter << " steps in getting saddlepoint eqs. to converge" << endl;}
  
  if(TooManyIter)
    {
      logfile << "Too many iterations in solving saddlepoint equations" << endl;
      return false;
    }

  thisT = myT;
  thisepsilon = myeps;

  if(TRACE)
    {
      cout << "Doing a final check on the saddlepoints" << endl;
      // final check, just to check:
      epsilon = thisepsilon; // set Kinvq using the new value of epsilon
      ConstructKinvq(); // also subtract minimal value and renormalize Delta
      
      myT=CalculateT(0);  // only valid for one sublattice here
#if defined LATTICEDISTORTIONS && defined ELASTIC
      epsoverT=CalculateEpsilonsOverT();
      myeps=epsoverT*myT;
#if defined CLAMPEDMODES
      for(int k=0; k<Nclampedmodes; k++){ myeps[clampedmodes[k]]=realtype(0.);}
#endif      
#endif
      epsdiff=thisepsilon-myeps;
       
      cout << "final, T: LHS:" << thisT << " eps: LHS:" << thisepsilon << endl;
      cout << "final, T: RHS:" << myT   << " eps: RHS:" << myeps << " Tdiff:" << thisT-myT << " epsdiff:" << epsdiff << endl;
    }
  
  if(TRACE) cout << "Finished SolveSaddlePointEquations() " << endl;
  return true;
}


// assumes initialized Sigma and mu
bool Driver::SolveSelfConsistentEquation(NumberList Delta,bool load_state)
{
  if(TRACE) cout << "Starting SolveSelfConsistentEquation " << endl;

  //  if(TRACE) cout << "Starting Solve with Delta= " << Delta << " epsilon=" << epsilon << endl;
  //  logfile << "Starting Solver with Delta= " << Delta << " epsilon=" << epsilon << endl;

  
  bool loadstate(load_state);
  bool go_on=true; // set the return variable to true by default, it becomes false when converged if EXITWHENCONVERGED is set
  failed=false;

  Chomp(Jq); // set very small entries to 0

  if(TRACE) SanityCheck(Jq,"Jq, input to SolveSelfConsistentEquation");

  //  if(TRACELEVEL>3) cout << "Jq=" << Jq << endl;
  
  SubtractMinimumEigenvalue(Jq);

  Chomp(Jq); // set very small entries to 0

  if(TRACE) SanityCheck(Jq,"Jq, after subtracting minimum");
  
  //  realtype T=-1;

  realtype m2=0.; // magnetic order parameter squared.
  vector<obstype> nobs(NOBSERVABLES); // nematic order parameters
  vector<obstype> nspinobs(NSPINOBSERVABLES); // different types of spin order parameters


  
#ifdef RANDOMINITIALIZATION
  MakeRandomSigma();
  SetQsToZero();
  BiasSigma();
  InitializeKinvqFromFile=true;
#elif defined USELASTSIGMA
  if(!SigmaInitialized || failed){MakeRandomSigma(); SetQsToZero(); BiasSigma(); SigmaInitialized=true; InitializeKinvqFromFile=true; MaxIterMultiplier=10;}
  else{MaxIterMultiplier=1; InitializeKinvqFromFile=false;}
#else
  logfile << "Initializing Sigma from rules-file" << endl;
  rule.InitializeSigma(Sigmaq); // Get initial values of Sigmaq
  SetQsToZero();
  BiasSigma();
  InitializeKinvqFromFile=true;
#endif

  
#ifdef PRESERVESYMMETRY
  MakeSymmetric(Sigmaq);
#endif


  // construct Kinvq:
  
  ConstructKinvq(true);

  if(TRACE) SanityCheck(Kinvq,"Kinvq, first time construction");
  
  if(TRACELEVEL>=4) cout << "Kinvq=" << Kinvq << endl;

  
  realtype oldT=Ts[0]; 

  //  currT=oldT; // set the current operating temperature
  newT=oldT; // set the current operating temperature

  realtype absTdev(0.);
  realtype oldabsTdev(0.);
  int nincreases(0); // number of times that a new iteration gets a T that deviates more from the oldT than the previous iteration. Sign of no convergence.


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
  bool reachedMAXITER=false;
  bool saddlepts_ok=true; // an indicator to flag whether iterations of saddlept eqs are succesful or not.

  converged=false;
  pconverged=false; // convergence in the previous iteration,
  done=false;


  // Load State variables from file STATEFILENAME
  bool stateloaded=false;  
  if( loadstate && LoadState(STATEFILENAME) )
    {
      loadstate=false;
      stateloaded=true;
    }
      
  while(!done && !reachedMAXITER && saddlepts_ok && !failed)  
    {
      if(TRACE) cout << "New iteration: " << iter << endl;
      iter++;
      if(converged){pconverged=true;} // record if the previous iteration had converged
      else{ pconverged=false;}
      
      if(TRACE) SanityCheck(Kinvq,"Kinvq, at start of new iteration");

      if(TRACELEVEL>=4) cout << "Kinvq=" << Kinvq << endl;
      
      ComputeSelfEnergy(); // careful with this, it overwrites K

#ifdef NOSELFENERGY
      Sigmaq.SetToZero();
#endif

      if(TRACE) SanityCheck(Sigmaq,"Sigmaq, after ComputeSelfEnergy");

      if(TRACELEVEL>=4) cout << "Sigmaq=" << Sigmaq << endl;

      
      ConstructKinvq();

      if(TRACE) SanityCheck(Kinvq,"Kinvq, in iteration loop after making Sigmaq");

      if(TRACELEVEL>=4) cout << "Kinvq=" << Kinvq << endl;
      
      // This is where to solve the saddlepointequations
      //-------------------------------------------------------------

      NumberList neweps(epsilon);
      saddlepts_ok=SolveSaddlePointEquations(newT,neweps); // return true if converged, false if not, modifies newT

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


      if(TRACE)
	{
	  cout << "iteration: " << iter << " T= " << newT << " oldT=" << oldT
	       << " dev=" << absTdev;
#if defined LATTICEDISTORTIONS && defined ELASTIC
	  cout << " epsilon= " << neweps;
#endif
	  cout << endl;
	}
      
      // write progress to logfile
      if(PRINTPROGRESS && iter%PRINTPROGRESSTICKLER==0)
	{
	  
	  logfile << "iteration: " << iter << " T=";
	  logfile << scientific << setprecision(LOGPRECISION);
	  for(int i=0; i<NSUBL; i++)
	    logfile << setw(LOGPRECISION+7) << Ts[i] << " ";
	  logfile << "oldT=" << setw(LOGPRECISION+7) << oldT << " dev=" << setw(LOGPRECISION+7) << absTdev << " ninc=" << setw(2) << nincreases;
#if defined LATTICEDISTORTIONS && defined ELASTIC
	  logfile << " epsilon=" << setw(LOGPRECISION+8) << neweps;
#endif
	  logfile << endl;
	}
      
      if( absTdev < par[TOLERANCE]) converged=true;

      //      const realtype inertia=0.5; // how much to resist changes: [0,1] 
      const realtype inertia=realtype(0.2); // how much to resist changes: [0,1] 

#if defined LATTICEDISTORTIONS && defined ELASTIC     
      for(int i=0; i<NELASTIC; i++)
	{
	  //	  epsilon[i]= (1.-inertia)*currT*epsoverT[i]+inertia*epsilon[i];
	  epsilon[i]= (1.-inertia)*neweps[i]+inertia*epsilon[i]; // maybe change to this
	}
#if defined CLAMPEDMODES
      for(int k=0; k<Nclampedmodes; k++){ epsilon[clampedmodes[k]]=realtype(0.);}
#endif      
#endif       

      //currT= (1.-inertia)*newT+inertia*oldT; 
      //oldT=currT;
      newT= (realtype(1.)-inertia)*newT+inertia*oldT; 
      oldT=newT;

      
      if(TRACE) cout << "converged= " << converged << " TOLERANCE:" << par[TOLERANCE] << endl;
      
      if(converged && pconverged){ done=true; continue;} // two iterations must fulfill conv. crit.

      if(nincreases > MAXNINCREASES){ done =true; logfile << "Too many nincr=" << nincreases << ", Not converging, exiting" << endl; continue;} // to many increases in Tdiff, so exiting without converging.

      reachedMAXITER=(iter >= MaxIterMultiplier*par[MAXITER]);
    }

  //Use MAXITER as convergence criterion itself
  if(par[TOLERANCE]==0.)
    {
      logfile << "Exiting after MAXITER steps, using end result" << endl;
      converged=true;
    } 

  if(converged && !stateloaded)
    {
      SaveState(STATEFILENAME); // Save state if converged and the state was not loaded.

    }

  if(converged)
    {
      WriteEpsilon(EPSILONFILENAME);
      WriteKinvq(KINVQFILENAME);
    }
  
  realtype f(0.);
  if(converged){ f= CalculateFreeEnergy(newT);}
  
  realtype f_elastic(0.);
#if defined ELASTIC
  if(converged){ f_elastic= CalculateElasticFreeEnergy(newT);}
#endif
  
  realtype f_phonons(0.);
#if defined LATTICEDISTORTIONS 
  if(converged){ f_phonons= CalculatePhononFreeEnergy(newT);}
#endif

  
  //  if(TRACE) cout << "Final Kinv_q: " << Kinvq << endl;  

  //  nalphas=CalculateAlphas(newT); // calculate alphas

  m2=NFAKESPINTRACE*NSPIN*newT/(realtype(2.)*Delta[0]*Nq); // calculate magnetic moment

  //  CalculateTs(Ts); // calculate the final temperatures
  
  // Print final Ts and epsilons
  logfile << "iteration: " << iter << " Ts: ";
  streamsize ss=logfile.precision();
  logfile.precision(17);
  for(int i=0; i<NSUBL; i++){logfile << Ts[i] << " ";}
#if defined LATTICEDISTORTIONS && defined ELASTIC
  logfile << " epsilon: " << epsilon;
#endif
  logfile << " converged: " << (converged ? "true": "false") << endl;
  logfile.precision(ss);

      
#ifdef PRINTTCONVERGENCE
      tfile << setprecision(17) << newT << endl;
      sfile << setprecision(17) << nobs[2].real() << endl;
#endif
      
  if(converged)
    {
      logfile << "Convergence reached after " << iter << " steps." << endl;
      if(TRACE) cout << "Convergence reached after " << iter << " steps." << endl;
    }


  if(!converged)
    {
      if(reachedMAXITER)
	logfile << "reached MAXITER=" << par[MAXITER] << " iterations without converging, increase MAXITER!" << endl;
      else
	logfile << "no uniform convergence of iterations, exiting!" << endl;
      SigmaInitialized   = false;
      EpsilonInitialized = false;
      InitializeKinvqFromFile= true;
    }
  else
    {      
      lineid++; 

      //      T = newT;

      

      realtype factor=newT*NSPIN*NFAKESPINTRACE*realtype(0.5); // conversion factor from Kinvq to Spin-corr  
      
      vector<Qandvals> maxqs(NMAXQS); // sorted list of qs and their values, biggest first
      
      for(int qi=0; qi<la.NqSites(); qi++) 
	{
	  Coord q1bz=la.qPos(qi);
	  for(int Z=0; Z<la.Nextendedzones; Z++)
	    {
	      Coord bzorigin=la.extendedzonesorigin[Z];
	      
	      realtype val(0.);
	      Coord q=q1bz+bzorigin;
	      
	      for(int s1=0; s1<NSUBL; s1++)
		for(int s2=0; s2<NSUBL; s2++)
		  {
		    val+= real(Kinvq(qi,s1,s2)*rule.GetCorrFactor(s1,s2,q)); 
		  }
	      val *=factor;
	      
	      if(val > maxqs[NMAXQS-1].v) // compare with last element, insert if true
		{
		  Qandvals newelem(q,val);
		  vector<Qandvals>::iterator iter; // an iterator to the insertion point
		  iter=upper_bound(maxqs.begin(),maxqs.begin()+NMAXQS,newelem);
		  maxqs.insert(iter,newelem); // insert and keep it sorted
		  maxqs.resize(NMAXQS); 
		}
	    }
	}
      
      realtype maxval=maxqs[0].v; // the overall maximum value. To be used for picking the largest corrs.
      
      const int OUTTP = 6; // T digits in printout
      const int OUTVP = 3; // val digits in printout
      const int OUTQP = 2; // q digits in printout
      
      
      ofstream outfile_maxqs(MAXQNAME.c_str(),ios::app);
      outfile_maxqs << lineid << " T= " << scientific << setprecision(OUTTP) << setw(OUTTP+6) << newT
		    << " Delta: " << scientific << setprecision(OUTTP) << setw(OUTTP+6) << Delta[0] << endl;
      int colsperline=2;
      int colcounter= 0;
      realtype oldval=0;
      for(int i=0; i<NMAXQS; i++)
	{
	  realtype v=maxqs[i].v;
	  Coord    q=maxqs[i].q;
	  Coord    qb=la.GetUVW(q);
	  string   qid="";
	  Coord    qBz=la.TranslateToFirstBZ(q);
	  bool idfound=la.FindQid(q,qid);
	  if( (colcounter%colsperline == 0 || v != oldval) && i != 0 ){outfile_maxqs << "\n"; colcounter=0;} 
	  outfile_maxqs << scientific << setprecision(OUTVP) << setw(OUTVP+6) << v;
	  outfile_maxqs << fixed << setprecision(OUTQP)
			<< " (" << setw(OUTQP+4) <<  q.x  << "," << setw(OUTQP+4) <<  q.y   << "," << setw(OUTQP+4) <<  q.z << ")"
			<< " [" << setw(OUTQP+3) << qb.x  << "," << setw(OUTQP+3) <<  qb.y  << "," << setw(OUTQP+3) <<  qb.z << "]"
			<< " (" << setw(OUTQP+4) << qBz.x << "," << setw(OUTQP+4) <<  qBz.y << "," << setw(OUTQP+4) <<  qBz.z << ")"
			<< setw(5) << qid << "  ";
	  colcounter++;
	  oldval=v;
	  
	  if(i==0) // also write to logfile
	    {
	      logfile << "maxQ:" << fixed << setprecision(OUTQP)
		      << " (" << setw(OUTQP+4) <<  q.x  << "," << setw(OUTQP+4) << q.y   << "," << setw(OUTQP+4) << q.z << ")"
		      << " [" << setw(OUTQP+2) << qb.x  << "," << setw(OUTQP+2) << qb.y  << "," << setw(OUTQP+2) << qb.z << "]"
		      << " (" << setw(OUTQP+4) << qBz.x << "," << setw(OUTQP+4) << qBz.y << "," << setw(OUTQP+4) << qBz.z << ")";
	      if(idfound){ logfile << setw(4) << qid << endl;}
	      else{ logfile << " not identified" << endl;}
	    }
	}
      outfile_maxqs << endl;
      outfile_maxqs << endl;
      
      outfile_maxqs.close();
      

      // look for max among all matrix elements in first brillouinzone
      // ONLY OUTPUT Kinvq values, do not multiply by factor=T*NS*realtype(0.5)
      vector<IndxVal> maxlist = FindMatrixMaxVals(Kinvq,NMAXQSINFIRSTBZ);


      ofstream mfile(MAXQSINFIRSTBZNAME.c_str(),ios::app);
      mfile << lineid << " T= " << scientific << setprecision(OUTTP) << setw(OUTTP+6) << newT
	    << " Delta: " << scientific << setprecision(OUTTP) << setw(OUTTP+6) << Delta[0] << endl;

      mfile << "  epsilons:";
      for(int i=0; i<NELASTIC; i++) mfile << scientific << setprecision(OUTTP) << setw(OUTTP+6) << epsilon[i] << " ";
      mfile << endl;
      
      realtype absmax=maxlist[0].val; // the absolute biggest element
      
      const int OUTMP=2; // precision of matrix elements
      for(int i=0; i<NMAXQSINFIRSTBZ; i++)
	{
	  realtype v=maxlist[i].val;

	  if(v < absmax/10.) break; // do not write out matrix if max is less than 10% of absolute max
	  
	  int      qi=maxlist[i].indx;
	  Coord    q=la.qPos(qi);
	  Coord    qb=la.GetUVW(q);
	  string   qid="";
	  Coord    qfd=la.TranslateToFundamentalDomain(q);
	  Coord    qBz=la.TranslateToFirstBZ(q);

	  la.FindQid(q,qid);

	  mfile << scientific << setprecision(OUTVP) << setw(OUTVP+6) << v;
	  mfile << fixed << setprecision(OUTQP)
		<< " (" << setw(OUTQP+4) <<  q.x  << "," << setw(OUTQP+4) <<  q.y   << "," << setw(OUTQP+4) <<  q.z << ")"
		<< " [" << setw(OUTQP+3) << qb.x  << "," << setw(OUTQP+3) <<  qb.y  << "," << setw(OUTQP+3) <<  qb.z << "]"
		<< " (" << setw(OUTQP+4) << qBz.x << "," << setw(OUTQP+4) <<  qBz.y << "," << setw(OUTQP+4) <<  qBz.z << ")"
		<< setw(5) << qid << endl;

	  mfile << setprecision(10) << fixed
		<< setw(14) << qfd.x << " " << setw(14) <<  qfd.y << " " << setw(14) <<  qfd.z << endl;
	  
	  mfile << scientific << setprecision(OUTMP);
	  for(int s1=0; s1<NMAT; s1++)
	    {
	      for(int s2=0; s2<NMAT; s2++)
		{
		  complex<realtype> mval(Kinvq(qi,s1,s2));
		  realtype mval_real= (abs(real(mval)) < 1.e-6 ? 0.: real(mval));
		  realtype mval_imag= (abs(imag(mval)) < 1.e-6 ? 0.: imag(mval)); 
	      
		  mfile << setw(OUTMP+8) << mval_real << setw(OUTMP+8) << mval_imag << "   ";


		}
	      mfile << endl;
	    }
	  mfile << endl;	  
	}
      mfile << endl;	  
      mfile << endl;
      
      mfile.close();


      if( (PRINTQCORRS && PRINTLARGESTQCORRS && lineid%PRINTTICKLER==0) || (PRINTPHONONSPECTRUM && lineid%PRINTPHONONSPECTRUMTICKLER==0))
	{
#ifndef SEPARATEZONEFILES
	  stringstream ss;
	  ss << LARGESTQCORRSFILENAME << "_" << lineid << ".dat";
	  //	  cout << "Opening " << ss.str() << endl;
	  ofstream qcorrfile(ss.str().c_str());
#endif	  
	  for(int Z=0; Z<la.Nextendedzones; Z++)
	    {	      
#ifdef SEPARATEZONEFILES
	      stringstream ss;
	      ss << LARGESTQCORRSFILENAME << "_Z" << Z << "_" << lineid << ".dat";
	      ofstream qcorrfile(ss.str().c_str());
#endif
	      qcorrfile << "# lineid=" << lineid << " T=" << newT << endl;
	      Coord bzorigin=la.extendedzonesorigin[Z];
	      
	      for(int qi=0; qi<la.NqSites(); qi++) 
		{
		  Coord q1bz=la.TranslateToFirstBZ(la.qPos(qi));
		  
		  realtype val(0.);
		  Coord q=q1bz+bzorigin;
#ifdef REDUCEDOUTPUT
		  if(!InRegionOfInterest(q)){ continue;}
#endif
		  
		  for(int s1=0; s1<NMAT; s1++)
		    for(int s2=0; s2<NMAT; s2++)
		      {
			val+= real(Kinvq(qi,s1,s2)*rule.GetCorrFactor(s1,s2,q)); 
		      }
		  val *= factor;

		  //		  cout << val << " " << QCORRSTHRESHOLD*maxval;
		  if(val >= QCORRSTHRESHOLD*maxval)
		    {
		      qcorrfile << setprecision(10) << q << " " << val << endl;
		    }
		}
	    }
	}     





      /*      
      SMatrix<realtype,NMAT,NMAT> maxvals;
      SMatrix<int,NMAT,NMAT> maxqs;

      FindMaxVals(maxqs,maxvals);
      
      if(TRACE) cout << "Write diagonal magnetization" << endl;
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
      */

      /*
      if(TRACE) cout << "Write peakfile: " << PEAKREPORTNAME << endl;
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
      if(TRACE) cout << "Done writing peakreport" << endl;
      */
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
	      complextype* start=Kinvq(0,0);
	      
	      if(BINARYOUTFILES)
		{
		  qcorrfile.write((char*) &lineid,sizeof(lineid)); // general format
		  qcorrfile.write((char*) &la.nindx_q,sizeof(la.nindx_q));
		  for(int i=0; i<la.nindx_q; i++)
		    {
		      complextype value=factor*start[la.indx_site_q[i]];
		      qcorrfile.write((char*) &value,sizeof(value));
		    }
		}
	      else
		{
		  //		  qcorrfile << setprecision(16) << lineid << " ";
		  for(int i=0; i<la.nindx_q; i++)
		    {
		      complextype value=factor*start[la.indx_site_q[i]];
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
	    outfile_a << scientific << setprecision(OUTPUTPRECISION)
		      << setw(OUTPUTPRECISION+8) << newT << " "
		      << setw(OUTPUTPRECISION+8) << real(nspinobs[j]) << " "
 		      << setw(OUTPUTPRECISION+8) << imag(nspinobs[j]) << endl;
	    outfile_a.close();
	    
	    ss.str("");
	    ss << NAMESOFSPINOBSERVABLES[j] << "1.dat";
	    
	    ofstream outfile_b(ss.str().c_str(),ios::app);
	    outfile_b << scientific << setprecision(OUTPUTPRECISION)
		      << setw(OUTPUTPRECISION+8) << newT << " "
		      << setw(OUTPUTPRECISION+8) << abs(nspinobs[j]) << endl;
	    outfile_b.close();
	    
	    ss.str("");
	    ss << NAMESOFSPINOBSERVABLES[j] << "2.dat";
	    
	    ofstream outfile_c(ss.str().c_str(),ios::app);
	    outfile_c << scientific << setprecision(OUTPUTPRECISION)
		      << setw(OUTPUTPRECISION+8) << newT << " "
		      << setw(OUTPUTPRECISION+8) << norm(nspinobs[j]) << endl;
	    outfile_c.close();
	  }
      }
   
  


      
      logfile << "Magnetic moment: " << m2 << endl;
      
      stringstream ss;
      ss << "m2.dat";
      
      ofstream outfilem(ss.str().c_str(),ios::app);
      outfilem << scientific << setprecision(OUTPUTPRECISION)
	       << setw(OUTPUTPRECISION+8) << newT << " " << setw(OUTPUTPRECISION+8) << m2 << endl;
      outfilem.close();
      
     
      logfile << "Free energy: " << scientific << setprecision(LOGPRECISION) << setw(LOGPRECISION+8) << f << endl;
      
      ss.str("");
      ss << "tf.dat";
      
      ofstream outfile(ss.str().c_str(),ios::app);
      outfile << scientific << setprecision(OUTPUTPRECISION)
	      << setw(OUTPUTPRECISION+8) << newT << " " << setw(OUTPUTPRECISION+8) << f << endl;
      outfile.close();
      
      ss.str("");
      ss << "df.dat";
      
      ofstream outfile2(ss.str().c_str(),ios::app);
      outfile2 << scientific << setprecision(OUTPUTPRECISION)
	       << setw(OUTPUTPRECISION+8) << Delta[0] << " " << setw(OUTPUTPRECISION+8) << f << endl;
      outfile2.close();

      ss.str("");
      ss << "td.dat";
	  
      ofstream outfile_a(ss.str().c_str(),ios::app);
      outfile_a << scientific << setprecision(OUTPUTPRECISION)
		<< setw(OUTPUTPRECISION+8) << newT << " " << setw(OUTPUTPRECISION+8) << Delta[0] << endl;
      outfile_a.close();
      
      ss.str("");
      ss << "dt.dat";
      
      ofstream outfile_b(ss.str().c_str(),ios::app);
      outfile_b << scientific << setprecision(OUTPUTPRECISION)
		<< setw(OUTPUTPRECISION+8) << Delta[0] << " " << setw(OUTPUTPRECISION+8) << newT << endl;
      outfile_b.close();

      logfile << "Phonon Free energy: " << scientific << setprecision(LOGPRECISION) << setw(LOGPRECISION+8) << f_phonons << endl;
      
      ss.str("");
      ss << "tfphonons.dat";
      
      ofstream outfile_phonons(ss.str().c_str(),ios::app);
      outfile_phonons << scientific << setprecision(OUTPUTPRECISION)
	      << setw(OUTPUTPRECISION+8) << newT << " " << setw(OUTPUTPRECISION+8) << f_phonons << endl;
      outfile_phonons.close();

      logfile << "Elastic Free energy: " << scientific << setprecision(LOGPRECISION) << setw(LOGPRECISION+8) << f_elastic << endl;
      
      ss.str("");
      ss << "tfelastic.dat";
      
      ofstream outfile_elastic(ss.str().c_str(),ios::app);
      outfile_elastic << scientific << setprecision(OUTPUTPRECISION)
	      << setw(OUTPUTPRECISION+8) << newT << " " << setw(OUTPUTPRECISION+8) << f_elastic << endl;
      outfile_elastic.close();

      
#if defined LATTICEDISTORTIONS && defined PHONONS
      if( lineid % PRINTPHONONSPECTRUMTICKLER == 0 && PRINTPHONONSPECTRUM)
	{
	  MatrixInverse(Dq); // B=Dinvq
	  rule.phonons.PrintPhononModes("renormalizedphonons_"+int2string(lineid)+".dat",Dinvq,newT);
	  MatrixInverse(Dinvq); // B=Dq
	}
#endif
      
#if defined LATTICEDISTORTIONS && defined ELASTIC
      for(int i=0; i<NELASTIC; i++)
	{
	  stringstream ss;
	  ss << "teps" << "_" << i << ".dat";
	  ofstream outfile_a(ss.str().c_str(),ios::app);
	  outfile_a << setprecision(16) << newT << " " << epsilon[i] << endl;
	  outfile_a.close();	
	}
#endif



     
#if defined EXITWHENCONVERGED
      go_on = false;
#endif     
    }
  
  if(TRACE) cout << "Done SolveSelfConsistentEquation " << endl;
  return go_on;
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
      
      realtype D2=(logscale ? mylog(par[DELTASLUTT]) : par[DELTASLUTT]);
      realtype D1=(logscale ? mylog(par[DELTASTART]) : par[DELTASTART]);

      realtype dDelta=(NDeltas>1 ? (D2-D1)/(NDeltas-1):0);
      
      if(TRACE) 
        cout << "Creating " << NDeltas << " on " << (logscale ? "log" : "linear") << " scale from " 
             << par[DELTASTART] << " to " << par[DELTASLUTT] << endl;
      logfile << "Creating " << NDeltas << " on " << (logscale ? "log" : "linear") << " scale from " 
	      << par[DELTASTART] << " to " << par[DELTASLUTT] << endl;

      for(int i=0; i<NDeltas; i++)
        {
          realtype this_Delta=(logscale ? myexp(D1+dDelta*i): D1+dDelta*i);
          NumberList myval(NSUBL,this_Delta); // set them all equal
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
	{
	  cout << "No file " << DELTASTOSHOWFILENAME << " found." 
	       << " Using Delta=" << par[DELTASTART] << endl;
	}
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
  ifstream epsilonfile(INPUTEPSILONFILENAME.c_str());
  if(!epsilonfile)
    {
      // set default starting values
      NumberList newepsilon(NELASTIC);
      for(int i=0; i<NELASTIC; i++) newepsilon[i]=(RAN()-0.5)*par[EPSILONNULL]; // small values

#if defined CLAMPEDMODES
      for(int k=0; k<Nclampedmodes; k++){ newepsilon[clampedmodes[k]]=realtype(0.);}
#endif            

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
      bool loadstate=(i==0 ? true: false);
      bool go_on = mysolver.Solve(Deltalist[i],epsilonlist[i],Printinfolist[i],loadstate);
      if(!go_on){logfile << "Stopping run because go_on=false" << endl; break;}
    }
  if(TRACE) cout << "Done Run" << endl;
}



#endif //BOND_H

