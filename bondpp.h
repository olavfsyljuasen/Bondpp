#ifndef BOND_H
#define BOND_H

#include<vector>
#include<complex>
#include<iostream>
#include<fstream>
#include<algorithm>
#include<sstream>

#ifdef PHONONS
const int NDMAT=NSUBL+NMODE;
#else
const int NDMAT=NSUBL;
#endif
const int NDMAT2=NDMAT*NDMAT;

const int NS=1; // remove this 


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



struct NumberList
{
NumberList(int n=NSUBL,realtype a=0.):v(n,a){}
  friend ostream& operator<<(ostream& os,const NumberList& d){for(int i=0; i<d.v.size(); i++){ os << d.v[i] << " ";} return os;}
  vector<realtype> v;
}; 

bool operator==(NumberList& l,NumberList& r)
{
  bool isequal=true;
  for(int i=0; i<l.v.size(); i++){isequal &= (l.v[i]==r.v[i]); if(!isequal){break;}}
  return isequal;
}



class Driver
{
 public:
  Driver(realtype*,BravaisLattice&,Rule&);
  ~Driver()
    {
      FFTWDESTROYPLAN(A1q_to_A1r);
      FFTWDESTROYPLAN(A1r_to_A1q);
      FFTWDESTROYPLAN(A2q_to_A2r);
      FFTWDESTROYPLAN(A2r_to_A2q);
      FFTWDESTROYPLAN(Bq_to_Br);
      FFTWDESTROYPLAN(Br_to_Bq);
      FFTWDESTROYPLAN(F1r_to_F1q);
      FFTWDESTROYPLAN(F1q_to_F1r);
      FFTWDESTROYPLAN(F2r_to_F2q);
      FFTWDESTROYPLAN(F2q_to_F2r);
    }
  realtype CalculateT(int);
  void CalculateTs(vector<realtype>&);

#ifdef PHONONS
  realtype CalculateMus();
#endif

  realtype CalculateFreeEnergy(const realtype);
  vector<obstype> CalculateOrderPars(realtype,int,int);
  vector<obstype> CalculateAlphas(realtype);

  //void Convolve(const bool);
  //  realtype CalculateSecondDerivative(const realtype T,const int k);

  void ConstructKinvq();
  void ComputeDq(const bool,const bool);  
  void ComputeSelfEnergy(const bool);
  
  void SetQsToZero(); // routine to set some q's to zero in self-energy
  void MakeRandomSigma();

  void MakeSymmetric(VecMat<complex<realtype>>&);
  
  void SolveSelfConsistentEquation(vector<realtype> Delta); 

  void Solve(const NumberList delta,const bool pinfo)
  {
    if(TRACE) cout << "Starting Solve with Delta= " << delta << " Printinfo= " << pinfo << endl;
    logfile << "Starting Solver with Delta= " << delta << " Printinfo= " << pinfo << endl;
    
    for(int s=0; s<NSUBL; s++){Delta[s]=delta.v[s];} 

    for(int s=0; s<NELASTIC; s++){mu[s]=1.;} // default starting value 1


    Printinfo=pinfo;
    //    Sigma = rule.GetInitialState(); // copy the Initial guess for Sigma
    SolveSelfConsistentEquation(Delta);
    if(TRACE) cout << "Done Solve" << endl;
  }
  
  
 private:
  realtype* par;
  BravaisLattice& la; 
  Rule& rule;

  int dim;              // the number of dimensions
  vector<int> dims;     // the size of each dimension
  
  const int Vq; // number of q-space sites.
  const realtype invVq;
  
  vector<realtype> Delta; 
  vector<realtype> mu;  // elastic constants 
  
#ifdef PRESERVESYMMETRY
  int TransformationPeriod;
  vector<int> TransformationTable;
#endif
  
  bool Printinfo; // set this to print out correlation functions
  int lineid; 

  realtype mineigenvalue; // for storing the minimum SigmaE value
  //  realtype g; // the soft constraint variable
  realtype currT; // for storing the current value of the temperature 


  // The actual storage areas
  VecMat<complex<realtype>> A1;  // holds Kq,Kinvq
  VecMat<complex<realtype>> A2;  // holds Sigmaq,
  VecMat<complex<realtype>> B;  // holds Dq,Dinvq

  vector<complex<realtype>> F1q ; // holds intermediate Fourier-transform
  vector<complex<realtype>> F2q ; // holds intermediate Fourier-transform

#ifdef PHONONS
  VecMat<complex<realtype>>& f; // holds the vertex information
  VecMat<complex<realtype>>& g; // points to rule

  vector<VecMat<complex<realtype>>*> gel;
#endif

  VecMat<complex<realtype>>& Jq;


#ifdef FFTS_INPLACE
  VecMat<complex<realtype>>& A1r; // holds Kinvr
  VecMat<complex<realtype>>& A2r; // holds Sigmar
  VecMat<complex<realtype>>& Br; // holds Dr
  vector<complex<realtype>>& F1r; // 
  vector<complex<realtype>>& F2r; // 
#else   
  VecMat<complex<realtype>> A1r;  // holds Kinvr
  VecMat<complex<realtype>> A2r;  // holds Sigmar
  VecMat<complex<realtype>> Br;  // holds Dr
  vector<complex<realtype>> F1r; // 
  vector<complex<realtype>> F2r; // 
#endif


  // the following references are used for readability of the code


  VecMat<complex<realtype>>& Kq;     // points to the A1 array
  VecMat<complex<realtype>>& Kinvq;  // points to the A1 array
  VecMat<complex<realtype>>& Kinvr;  // points to the A1 array

  VecMat<complex<realtype>>& Sigmar; // points to the A2 array
  VecMat<complex<realtype>>& Sigmaq; // points to the A2 array

  VecMat<complex<realtype>>& Dq;     // points to the B array
  VecMat<complex<realtype>>& Dinvq;  // points to the B array
  VecMat<complex<realtype>>& Dr;     // points to the B array


  FFTWPLAN A1q_to_A1r;
  FFTWPLAN A1r_to_A1q;
  FFTWPLAN A2q_to_A2r;
  FFTWPLAN A2r_to_A2q;

  FFTWPLAN Bq_to_Br;
  FFTWPLAN Br_to_Bq;

  FFTWPLAN F1q_to_F1r;
  FFTWPLAN F1r_to_F1q;

  FFTWPLAN F2q_to_F2r;
  FFTWPLAN F2r_to_F2q;


};


Driver::Driver(realtype* in_par,BravaisLattice& in_la,Rule& r):par(in_par),la(in_la),rule(r),dim(la.D()),dims(la.SiteqDims()),Vq(la.SiteqVol()),invVq(static_cast<realtype>(1.)/Vq),Delta(NSUBL),Printinfo(false),lineid(0),
  A1(Vq,NMAT,NMAT),A2(Vq,NMAT,NMAT),B(Vq,NDMAT,NDMAT),F1q(Vq),F2q(Vq),
#ifdef FFTS_INPLACE
  A1r(A1),A2r(A2),Br(B),F1r(F1q),F2r(F2q),
#else
  A1r(Vq,NMAT,NMAT),A2r(Vq,NMAT,NMAT),Br(Vq,NDMAT,NDMAT),F1r(Vq),F2r(Vq),
#endif
  Kq(A1),
  Kinvq(A1),Kinvr(A1r),Sigmar(A2r),Sigmaq(A2),
  Dq(B),Dinvq(B),Dr(Br),Jq(r.Jq)  
#ifdef PHONONS
  ,f(r.Getf()),g(r.g)
  ,mu(NELASTIC)
#else
  ,mu(0)
#endif
{
  if(TRACE) cout << "Initializing solver " << endl;




#include "fourierplans.h"


#ifdef PRESERVESYMMETRY
  TransformationPeriod=lattice.TransformationPeriod;
  TransformationTable=lattice.TransformationTable;
#endif

  
  if(TRACE) cout << "Done initializing solver " << endl;
};


realtype Driver::CalculateT(int sl)
{
  if(TRACE) cout << "Starting CalculateT for sublattice: " << sl << endl;
  realtype alpha=0.;
  for(int spin=0; spin<NSPIN; spin++) 
    {
      int m=mindx(spin,sl); // make the composite index.
      alpha=real(Sumq(Kinvq,m,m));
    }
  alpha *= 1./(2.*Vq);

  realtype T = 1./alpha;

  /*
#ifdef SOFTCONSTRAINT
  // remember to use the correct value of Delta which is not the tilde{Delta}
  const realtype b = (Delta[s]-mineigenvalue)/(2.*g);
  const realtype fouralphab = 4.*b*alpha;
  if( fouralphab < -1.){ cout << " WARNING 4alphadelta= " << fouralphab << " < -1" << endl;} 
  realtype T= (1. + sqrt( 1. + fouralphab))/(2.*alpha);
#else
  realtype T = 1./alpha;
#endif
  */
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
realtype Driver::CalculateMus()
{
  if(TRACE) cout << "Starting CalculateMus" << endl;

  for(int mui=0; mui<rule.Nelastic; mui++)
    {
      realtype lambda= rule.elasticeigenvalues[mui];
      realtype sum=0.;
      if(lambda != 0.)
	{
	  for(int qi=0; qi<Vq; qi++)
	    {
	      SMatrix<complex<realtype>> tmp(NMAT,NMAT);
	      tmp=Kinvq[qi];
	      tmp *= (*rule.gelptrs[mui])[qi];
	      sum += real(tr(tmp));
	    }
	  sum *= currT/(2.*Vq*lambda);
	}
      mu[mui]= sum;
    }
 
  if(TRACE) cout << "Done CalculateMus "  << endl;
}
#endif


vector<obstype> Driver::CalculateOrderPars(realtype T,int m1,int m2)
{
  if(TRACE) cout << "Starting CalculateOrderPars" << endl;

  vector<obstype> opars(NOBSERVABLES);
  for(int j=0; j<NOBSERVABLES; j++)
    {
      obstype sum=0.;
      vector<obstype>& f=rule.GetIrrep(j);

      for(int i=0; i<Vq; i++)
	{
	  sum+= f[i]*Kinvq(i,m1,m2);
	}
      opars[j]=0.5*T*invVq*sum;
    }

  if(TRACE) cout << "Done CalculateOrderPars" << endl;
  return opars;
}

#ifdef PRESERVESYMMETRY
void Driver::MakeSymmetric(VecMat<complex<realtype>>& m)
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
  //cout << "Vq=" << Vq << endl;


  while( p < TransformationPeriod)
    {
      for(int i=0; i<Vq; i++)
	{
	  int q=ThisT[i];
	  tempstart[i]+=mstart[q];
	  ThisT[i] = TransformationTable[q]; // update ThisT for the next period
	}
      p++;
    }
  
  // Average over all transformations and transfer the result to the input array
  realtype invperiod=1./TransformationPeriod;
  for(int i=0; i<Vq; i++){ mstart[i]=tempstart[i]*invperiod;}  

  if(TRACE) cout << "End of MakeSymmetric" << endl;

}
#endif




// Free energy per unit volume
// we use the convention \nu=\nup
realtype Driver::CalculateFreeEnergy(realtype T)
{
  if(TRACE) cout << "Starting CalculateFreeEnergy " << endl;

  realtype f=0;
 
  if(TRACE) cout << "--- T  = " << T << endl;
  // constants:
  realtype betaf_constants= -( 0.5*NSUBL*log(2.*Vq) + 0.5*NSUBL*(NS-1)*log(PI)); 
  if(TRACE) cout << "betaf_constants  = " << betaf_constants << endl;
  
  f += T*betaf_constants;
  
  //Must correct the Delta values for the subtraction of the minimum from SigmaE
#ifdef SOFTCONSTRAINT
  for(int i=0; i<NSUBL; i++) f += -(Delta[i]-mineigenvalue)*(Delta[i]-mineigenvalue)/(4.*g*T);
#endif
  for(int i=0; i<NSUBL; i++) f += -(Delta[i]-mineigenvalue);

  if(TRACE) cout << "betaf_delta      = " << -(Delta[0]-mineigenvalue)/T << " (Delta= " << Delta[0] << " mineig: " << mineigenvalue << ")" << endl;
  
  // NSUBL*log(T) is needed because T is not part of Kinvq in the program
  realtype betaf_logKinvq   = -0.5*NS*(invVq*SumLogDet(Kinvq) + NSUBL*log(T));
  //  cout << "invVqq*SumLogDet= " << invVq*SumLogDet(Kinvq) << " log(T)= " << log(T) << endl;
  if(TRACE) cout << "betaf_logKinvq   =  " << betaf_logKinvq << endl;					  f += T*betaf_logKinvq;
  
  ComputeDq(false,true); // excludeqzero=false, preserveinput=true not to jeopardize Keff
  
  // NSUBL*2.*log(T) needed because T is not part of Dq in the program
  realtype betaf_logD       = -0.5*invVq*( SumLogDet(Dq) - Vq*NSUBL*2.*log(T) );
  if(TRACE) cout << "betaf_logD       =  " << betaf_logD << endl;
  f += T*betaf_logD;

  ComputeSelfEnergy(true); // ensure that Kinvq is the same.

  realtype betaf_KinvqSigma = -0.5*NS*invVq*SumTr(Kinvq,Sigmaq);
  if(TRACE) cout << "betaf_KinvqSigma = " << betaf_KinvqSigma << endl;
  f += T*betaf_KinvqSigma;

  // the correction to the saddle-point are already taken into account
  return f;
}


// ComputeDq() computes the renormalized constraint propagator
//
// Input: Kinvq (stored in A).
// Output: Dq (stored in B)
// 
// The routine uses in-place FFTs so info contained in A and B are modified
// On output:
// A = Kinvr/Vq, unless preserveinput=true, then A=Kinvq
// B = Dq
//
// FFTW omits volume prefactors in fourier transforms, we adopt the convention that the Fourier transform
// is without prefactors in going from q->r, and the prefactor 1/Vq is inserted on going from r->q.
// This means that using FFTW for transforming q->r->q, one should divide the result by Vq. 
void Driver::ComputeDq(const bool excludeqzero=true, const bool preserveinput=false)
{
#ifdef PHONONS
  if(NSUBL !=1){ cout << "Error:ComputeDq is not implemented for NSUBL !=1, exiting." << endl; exit(1);}
#endif

  if(TRACE) cout << "Starting ComputeDq, excludeqzero=" << excludeqzero << ",preserveinput=" << preserveinput << endl;

  if(TRACE)
    { 
      cout << "Kinv_q " << Kinvq << endl;
     }
  //  if(TRACE) cout << "SumLogDet(Kinvq)=" << SumLogDet(Kinvq) << endl; 

  if(TRACE) cout << "Kinv_q max imag value= " << FindMaxImag(Kinvq) << endl;

  // multiply by invVq so that the values in the intermediate steps is not too large
 
  if(TRACE) 
    {
      cout << "Hallo2, Kinvq=0" << endl;
      cout << Kinvq << endl;
   }
  

  //for sammenligning ta ut følgende linje:
  //  for(int i=0; i<Kinvq.size(); i++){Kinvq[i]*=invVq;}  // Aq =  Kinv_q / Vq
  Kinvq *= invVq;

  FFTWEXECUTE(A1q_to_A1r); // Kinvq->Kinvr, stored in Ar, so after transform: Ar=Kinvr/Vq

  if(TRACE) 
    {
      cout << "Hallo3, Kinv(r=0)=" << endl;
      cout << Kinvr << endl;
    }
  
#ifdef FORCEINVERSIONSYMMETRY
  MakeReal(Kinvr);  // is this needed here?
  MakeInversionTransposedSymmetric(la,Kinvr); 
#endif
  
  
  
  if(TRACE)
    { 
      cout << "after enforcing symmetry properties" << endl;
      cout << "Kinv_r /Vq " << Kinvr << endl;
    }
  
  
  
  //prepare Kinv Kinv kernel, store it in B to save space.
  //  const realtype Nsover2=NS/2.;
  //  for(int i=0; i<Vtot; i++){Br[i]=Nsover2*Kinvr[i]*Kinvr[i];} // Br = (NS/2) Kinv_r^2/(Vq*Vq)
  //  for(int i=0; i<Br.size(); i++){Br[i]=Nsover2*Kinvr[i]*conj(Kinvr[i]);} // Br = (NS/2) Kinv_r^2/(Vq*Vq)
  
  
  // first compute the constraint block

  complex<realtype> tmp(0);
  
  for(int l1=0; l1<NSUBL; l1++)
    for(int l2=l1; l2<NSUBL; l2++)
      {
	for(int i=0; i<Vq; i++)
	  {
	    complex<realtype> tt(0.);
	    for(int s1=0; s1<NSPIN; s1++)
	      for(int s2=0; s2<NSPIN; s2++)
		{
		  int m1=mindx(s1,l2); // lprime = l2
		  int m2=mindx(s2,l1); // l=l1
		  
		  tmp=Kinvr(i,m1,m2);
		  tt += tmp*tmp;   
		}
	    F1r[i]=0.5*tt;
	    F1r[i].imag(0.); // make sure the imaginary part is zero
	  }

	FFTWEXECUTE(F1r_to_F1q); // F1r->F1q 

	for(int i=0; i<Vq; i++)
	  {
	    Dinvq(i,l1,l2)=F1q[i];
	    Dinvq(i,l2,l1)=F1q[i]; 
	  }
      }
  
  
  if(TRACE) cout << "after constraint part Dinvq=" << endl;
  if(TRACE) cout << Dinvq << endl;

#ifdef PHONONS
  // the offdiagonal phonon-constraint part
  for(int l1=0; l1<NSUBL; l1++) // l1 must be 0 here because phonons NSUBL=1
    for(int l2=NSUBL; l2<NDMAT; l2++)
      {
	int n2=l2-NSUBL; 

	for(int c4=0; c4<NC; c4++)
	  {
	    for(int i=0; i<Vq; i++)
	      {
		SMatrix<complex<realtype>> my(NMAT,NMAT);
		my  = g[c4];
		my *= Kinvr[i];
		
		//		cout << "ph-co: i=" << i << " my=" << my << endl;
		
		F1r[i]=0;
		for(int s1=0; s1<NSPIN; s1++)
		  for(int s2=0; s2<NSPIN; s2++)
		    F1r[i] += 0.5*my(s1,s2)*my(s1,s2);  
	      }
	    
	    // do fourier-transform of theta
	    FFTWEXECUTE(F1r_to_F1q); 
	    
	    for(int qi=0; qi<Vq; qi++)
	      {
		if(qi==0)
		  {
		    // avoid using the zero mode of the phonons *\\/ *\/ */
		    Dinvq(0,l1,l2)=0.; 
		    Dinvq(0,l2,l1)=0.; 
		    continue; 
		  } 
		if(c4==0)
		  {	  
		    Dinvq(qi,l1,l2) = F1q[qi]*f(qi,c4,n2)*conj(expi(la.qr(qi,c4)));
		  } 
		else 
		  { 
		    Dinvq(qi,l1,l2)+= F1q[qi]*f(qi,c4,n2)*conj(expi(la.qr(qi,c4)));
		  }
		
		Dinvq(qi,l2,l1) =-Dinvq(qi,l1,l2);
	      }		   
 	  }
      }
  
  // finally the phonon-phonon part
  for(int l1=NSUBL; l1<NDMAT; l1++) 
    for(int l2=l1; l2<NDMAT; l2++) 
      {
	
	int n1=l1-NSUBL; 
	int n2=l2-NSUBL; 
	
	// should only choose positive c2 and c4 to sum over
 	for(int c2=0; c2<NC; c2++)
	  for(int c4=0; c4<NC; c4++)
	    {	    
	      Triplet myivec=clist[c2]+clist[c4];
	      
 	      for(int i=0; i<Vq; i++)
 		{
 		  SMatrix<complex<realtype>> my(NMAT,NMAT);
		  
 		  my  = g[c4];
 		  my *= Kinvr[i];
		  my *= g[c2];
		  
		  int newindx=la.rAdd(i,myivec);
		  
		  F1r[i]=0;
 		  for(int sd=0; sd < NSPIN; sd++) 
 		    for(int sg=0; sg < NSPIN; sg++)
 		      { 
 			F1r[i] += 0.5*my(sd,sg)*Kinvr(newindx,sd,sg);
 		      } 
 		} 
 	      FFTWEXECUTE(F1r_to_F1q); // F1r->F1q;
	      
 	      for(int qi=0; qi<Vq; qi++)
 		{
 		  if(qi==0) 
 		    {
 		      // avoid using the zero mode of the phonons
 		      Dinvq(0,l1,l2)=0.; 
 		      Dinvq(0,l2,l1)=0.; 
		      continue; 
		    }
		  
		  if(c2==0 && c4==0)
		    {
		      Dinvq(qi,l1,l2)=-F1q[qi]*conj(f(qi,c2,n1))*f(qi,c4,n2)*conj(expi(la.qr(qi,c4)));
		    }
		  else
		    {
		      Dinvq(qi,l1,l2)+=-F1q[qi]*conj(f(qi,c2,n1))*f(qi,c4,n2)*conj(expi(la.qr(qi,c4)));
		    }
		  Dinvq(qi,l2,l1)=Dinvq(qi,l1,l2);
		}
	    }
      }
  // add the bare phonon part

  realtype barepart=0.5/currT; 
  for(int qi=0; qi<Vq; qi++)
    {
      for(int n=0; n<NMODE; n++)
	{
	  int l=n+NSUBL;
	  Dinvq(qi,l,l) += barepart;
	}
    }
#endif
  
  if(TRACE) cout << "Dinvq=" << endl;
  if(TRACE) cout << Dinvq << endl;

  
  MakeHermitian(Dinvq);




  if(TRACE) cout << "Dinvq max imag value= " << FindMaxImag(Dinvq) << endl;

  Dinvq *= Vq;
  //  for(int i=0; i<Dinvq.size(); i++) Dinvq[i] *= Vq; // B = Dinv_q

  /*
#ifdef SOFTCONSTRAINT
  complex<realtype> softconstraintmodifier=complex<realtype>(0.5*Vq/g,0.);
  if(TRACE) cout << "Softconstraint: adding " <<  softconstraintmodifier << endl;
  AddToDiagonal(Dinvq,softconstraintmodifier);
#endif
  */
  //
  //  if(TRACE) cout << "SumLogDet(Dinvq)=" << SumLogDet(Dinvq) << endl; 
  //

  if(TRACE) cout << "Dinvq=" << endl;
  if(TRACE) cout << Dinvq << endl;

  MatrixInverse(Dinvq); // B = Dq

  //  if(TRACE) cout << "SumLogDet(Dq)=" << SumLogDet(Dq) << endl; 

  if(excludeqzero) Setqzerotozero(Dq);

  // experiemental code, added 24.6.2020:
  //MakeHermitian(Dq);

  if(TRACE)
    { 
      cout << "Dq=" << Dq << endl;
    }

  if(TRACE) cout << "Dq max imag value= " << FindMaxImag(Dq) << endl;

  //
  if(TRACE) cout << "Dq" << endl;
  if(TRACE) cout << Dq << endl;
  //

  if(preserveinput)
    {
      FFTWEXECUTE(A1r_to_A1q);  // Kinvr->Kinvq, so after transform: Aq=Kinvq 

#ifdef PRESERVESYMMETRY
      MakeSymmetric(Kinvq);
#endif
      MakeHermitian(Kinvq);
    }

  if(TRACE) cout << "Done with ComputeDq " << endl;
}


// ComputeSelfEnergy() computes the self-energy from the self-consistent equations
//
// Input: Kinv_q (stored in A).
// Output: Sigma_q (stored in B)
// 
// The routine uses in-place FFTs so info contained in A and B are modified
// On output:
// A = Kinvr/Vq, unless preserveinput=true, then A=Kinvq
// B = Sigmaq 
//
// FFTW omits volume prefactors in fourier transforms, we adopt the convention that the Fourier transform
// is without prefactors in going from q->r, and the prefactor 1/Vq is inserted on going from r->q.
// This means that using FFTW for transforming q->r->q, one should divide the result by Vq. 
void Driver::ComputeSelfEnergy(const bool preserveinput=false)
{
  if(TRACE) cout << "Starting ComputeSelfEnergy" << endl;

  if(TRACE)
    { 
      cout << "Kinv_q " << Kinvq << endl;
    }

  /*
  if(TRACE)
    {
      cout << "Hallo SumLogDet(Kinvr)=" << endl;
      cout << SumLogDet(Kinvr) << endl;
    }
  */




  ComputeDq(true,false); // exclude zero mode, do not preserve input here, Ar must contain Kinvr/Vq for ComputeSelfEnergy to work.

  //  FFTWEXECUTE(Bq_to_Br); // Dq -> Dr, after transform: B= Dr

  // MakeReal(Dr); // Dr is real as follows from hermiticity and D_{ij q} = D*_{ji,q} = D_{ji,-q} 

  //if(TRACE) cout << "Dr r=0" << endl;
  //if(TRACE) cout << Dr[0] << endl;

  // Ar contains Kinvr/Vq

  // need to reset sigmaq

  Sigmaq.SetToZero();


  for(int m1=0; m1<NMAT; m1++)
    for(int m2=m1; m2<NMAT; m2++)
      {
	const int s1=spin(m1);
	const int s2=spin(m2);

	const int l1=subl(m1);
	const int l2=subl(m2);

	for(int q=0; q<Vq; q++){F1q[q]=Dq(q,l1,l2);}
	FFTWEXECUTE(F1q_to_F1r);
	
	for(int r=0; r<Vq; r++){F2r[r]=F1r[r]*conj(Kinvr(r,m2,m1));}
	FFTWEXECUTE(F2r_to_F2q);

	for(int k=0; k<Vq; k++){Sigmaq(k,m1,m2)=F2q[k];}

#ifdef PHONONS
	if(NSUBL != 1){cout << "Error: PHONONS are not implemented for NSUBL!=1, exiting" << endl; exit(1);}

	for(int c=0; c<NC; c++)
	  {
	    	for(int q=0; q<Vq; q++)
		  {
		    F1q[q]=0;
		    for(int n=0; n<NMODE; n++){ F1q[q]+=Dq(q,l1+n,l2)*f(q,c,n);}
		  }
		FFTWEXECUTE(F1q_to_F1r);

		for(int r=0; r<Vq; r++)
		  {
		    F2r[r]=0;
		    for(int s=0; s<NSPIN; s++)
		      F2r[r]+=F1r[r]*conj(Kinvr(r,s,m1))*g(c,s,m1);
		  }
		FFTWEXECUTE(F2r_to_F2q);

		for(int k=0; k<Vq; k++){Sigmaq(k,m1,m2)+=2*real(expi(la.qr(k,c))*F2q[k]);}
	  }

	for(int c1=0; c1<NC; c1++)
	  for(int c2=0; c2<NC; c2++)
	  {
	    for(int q=0; q<Vq; q++)
		  {
		    F1q[q]=0;
		    
		    for(int n1=0; n1<NMODE; n1++)
		      for(int n2=0; n2<NMODE; n2++)
			{ F1q[q]+=-Dq(q,l1+n1,l2+n2)*f(q,c1,n1)*conj(f(q,c2,n2))*expi(la.qr(q,c1));}
		  }
		FFTWEXECUTE(F1q_to_F1r);

		for(int r=0; r<Vq; r++)
		  {
		    F2r[r]=0;
		    for(int s3=0; s3<NSPIN; s3++)
		      for(int s4=0; s4<NSPIN; s4++)
			F2r[r]+=F1r[r]*g(c1,s1,s3)*conj(Kinvr(r,s3,s4))*g(c2,s4,s2);
		  }
		FFTWEXECUTE(F2r_to_F2q);
		
		for(int k=0; k<Vq; k++){Sigmaq(k,m1,m2)+=F2q[k]*expi(la.qr(k,c1))*expi(la.qr(k,c2));}
	  }
#endif
      }


#ifdef PRESERVESYMMETRY
  MakeSymmetric(Sigmaq);
#endif  
  MakeHermitian(Sigmaq);

  if(TRACE)
    { 
      cout << "Sigmaq etter hermitian" << Sigmaq << endl;
      cout << "Sigmaq(q=0) etter hermitian" << endl;
    }



  if(preserveinput)
    {
      FFTWEXECUTE(A1r_to_A1q); // Kinvr->Kinvq, so after transform: Aq=Kinvq


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
  
  //  if(TRACE) cout << "SumLogDet(Sigmaq)=" << SumLogDet(Sigmaq) << endl;

  Kq =  Jq;

  if(TRACE) cout << "Initializing Kq=Jq = " << Kq << endl;

  Kq += Sigmaq;

  if(TRACE) cout << "Added Sigmaq, Kq = " << Kq << endl;

#ifdef PHONONS
  for(int j=0; j<NELASTIC; j++)
    {
      VecMat<complex<realtype>> temp(*rule.gelptrs[j]);
      temp *= mu[j];
      Kq += temp;
    }
#endif


  if(TRACE) cout << "Kq=" << Kq << endl; 
  //if(TRACE) cout << "SumLogDet(Kq)=" << SumLogDet(Kq) << endl;

  mineigenvalue=SubtractMinimumEigenvalue(Kq);
  
  if(TRACE) cout << "Kq=" << Kq << endl; 
  
  AddDelta(Kq,Delta);
  
  if(TRACE) cout << "Kq=" << Kq << endl; 

  if(TRACE && !IsHermitian(Kq))
    {
      cout << "Warning: Kq is not Hermitian" << endl;
    }

  if(TRACE) cout << "SumLogDet(Kq)=" << SumLogDet(Kq) << endl;

  // construct K1inv
  MatrixInverse(Kq); // K = Kinv_q 

  //Experimental code:
  if(TRACE && !IsHermitian(Kinvq))
    {
      cout << "Warning: Kinvq is not Hermitian" << endl;
    }
  // experimental code, added 24.6.2020:
  // MakeHermitian(Kinvq); // force the inverse to be Hemrmititan


  if(TRACE) cout << "SumLogDet(Kinvq)=" << SumLogDet(Kinvq) << endl;

  if(TRACE) cout << "Done with ConstructKinvq" << endl;
}

/*
// This should work for several sublattices
void Driver::MakeRandomSigma()
{
  logfile << "Making a random (hermitian) initialization of the self-energy" << endl;

  realtype da = param[DA];
  
  for(int s1=0; s1<NSUBL; s1++)
    {
      for(int i=0; i<Vq; i++)
	{
	  complex<realtype> c=da*complex<realtype>(RAN()-0.5,0.); // real value 
	  Sigmar(s1,s1)[i]=c;
	}
    }

  for(int s1=0; s1<NSUBL; s1++)
    for(int s2=s1+1; s2<NSUBL; s2++)
	  {
	    for(int i=0; i<Vq; i++)
	      {
		complex<realtype> c=da*complex<realtype>(RAN()-0.5,0.);
		Sigmar(s1,s2)[i]=c;
	      }
	  }
  fftw_execute(Br_to_Bq);

#ifdef PRESERVESYMMETRY
  MakeSymmetric(Sigmaq);
#endif

  MakeHermitian(Sigmaq);
}
*/


void Driver::MakeRandomSigma()
{
  if(TRACE) cout << "Starting MakeRandomSigma()" << endl;
  logfile << "Making a random initialization of the self-energy" << endl;

  realtype da = par[DA];

  for(int m=0; m<NMAT; m++)
    {
      for(int i=0; i<Vq; i++)
	{
	  complex<realtype> c=da*complex<realtype>(RAN(),0.); // real positive value 
	  Sigmaq(i,m,m)=c;
	}
    }

  if(NMAT>1)
    {
      for(int m1=0; m1<NMAT; m1++)
	for(int m2=m1+1; m2<NMAT; m2++)
	  {
	    for(int i=0; i<Vq; i++)
	      {
		complex<realtype> c=da*complex<realtype>(RAN(),RAN());
		Sigmaq(i,m1,m2)=c;
	      }
	  }
    }


  MakeHermitian(Sigmaq);

#ifdef FORCEINVERSIONSYMMETRY
  MakeInversionTransposedSymmetric(la,Sigmaq);
  /*
  FFTWEXECUTE(A2q_to_A2r);
  MakeReal(Sigmar);
  FFTWEXECUTE(A2r_to_A2q);
  Sigmaq*=invVq;  // do not change magnitude
  */
#endif

#ifdef PRESERVESYMMETRY
  MakeSymmetric(Sigmaq);
#endif

  MakeHermitian(Sigmaq);

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
      if(qindx >=0 && qindx <Vq && ifile)
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
void Driver::SolveSelfConsistentEquation(vector<realtype> Delta)
{
  if(TRACE) cout << "Starting SolveSelfConsistentEquation " << endl;

  lineid++;


#ifdef RANDOMINITIALIZATION
  MakeRandomSigma();
  if(TRACE) cout << "Sigmaq: " << Sigmaq << endl;
#else
  rule.InitializeSigma(Sigmaq); // Get initial values of Sigmaq
#endif
  SetQsToZero();

#ifdef PRESERVESYMMETRY
  MakeSymmetric(Sigmaq);
#endif

#ifdef BIASQS  // set some q's to zero in 

#endif


  if(TRACE) cout << "Sigma: " << Sigmaq << endl;

       
  if(TRACE)
    {
      cout << "Jq:" << Jq << endl;
      for(int i=0; i<Jq.Nvecs; i++)
	{
	  SMatrix<complex<realtype> > tmp(NMAT,NMAT,Jq[i]);
	  
	  cout << "Jq(q=" << i << ")=" << endl;
	  cout << tmp << endl;
	}
    }


  Chomp(Jq); // set very small entries to 0
  
  if(TRACE) cout << "Min eigenvalue: " << FindMinimumEigenvalue(Jq) << endl;

  SubtractMinimumEigenvalue(Jq);

  Chomp(Jq); // set very small entries to 0
 
  if(TRACE) cout << "Jq after chomp: " << Jq << endl;
	      


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

  //  vector<obstype> nalphas(NOBSERVABLES); // alpha


  // construct Kinvq:
  
  ConstructKinvq();

  realtype oldT=CalculateT(0); 

  currT=oldT; // set the current operating temperature

  vector<realtype> Ts(NSUBL); // A list of Ts
  
  if(TRACE)
    {
      cout << "Initial temperatures: " << endl;
      CalculateTs(Ts); // calculate all the temperatures
      
      
      for(int i=0; i<NSUBL; i++)
	cout << setprecision(17) << Ts[i] << " ";
      cout << endl;
    }



 
  int iter=0;
  bool converged=false;
  bool pconverged=true; // convergence in the previous iteration,
  bool done=false;
 
  while(iter< par[MAXITER] && !done)  
    {
      if(TRACE) cout << "New iteration: " << iter << endl;
      iter++;
      if(converged) pconverged=true; // record if the previous iteration had converged



      if(TRACE) cout << "Kinvq=" << Kinvq << endl;
      if(TRACE) cout << "Kinvq max imag value= " << FindMaxImag(Kinvq) << endl;


      ComputeSelfEnergy(); // careful with this, it overwrites K

#ifdef NOSELFENERGY
      for(int i=0; i<Vtot; i++){Sigmaq[i]=0.;}      
#endif

      if(TRACE) cout << "Sigmaq="<< Sigmaq << endl;
      if(TRACE) cout << "Sigmaq max imag value= " << FindMaxImag(Sigmaq) << endl;

      ConstructKinvq(); 

      if(TRACE) cout << "Kinvq=" << Kinvq << endl;      
      if(TRACE) cout << "Max element of Kinvq: " << FindMax(Kinvq) << endl;

      // These lines can be removed, here it is just for debugging
      if(TRACE)
	{
	  CalculateTs(Ts); // calculate all the temperatures
	  
	  
	  for(int i=0; i<NSUBL; i++)
	    cout << setprecision(17) << Ts[i] << " ";
	  cout << endl;
	}
      
      // convergence checks
      newT=CalculateT(0);
      currT=newT;

#ifdef PHONONS
      CalculateMus();
      if(TRACE)
	{
	  cout << "Elastic constants mu: ";
	  for(int i=0; i<NELASTIC; i++) cout << mu[i] << " ";
	  cout << endl;
	}
#endif

      if(TRACE) cout << "T= " << newT << " oldT= " << oldT << endl;

      if( fabs((newT-oldT)/oldT) < par[TOLERANCE]) converged=true;
      oldT=newT;

      if(TRACE) cout << "converged= " << converged << endl;
      
      if(converged && pconverged){ done=true; continue;} // two iterations must fulfill conv. crit.
    }



  
  if(TRACE) cout << "Final Kinv_q: " << Kinvq << endl;  

  //  nalphas=CalculateAlphas(newT); // calculate alphas

  m2=NS*newT/(2.*Delta[0]*Vq); // calculate magnetic moment

  CalculateTs(Ts); // calculate the final temperatures


  
  //  logfile << iter << " T=" << newT << " " << m2 << " ";
  logfile << iter << " Ts=" << Ts << " ";
  logfile << endl;

  for(int i=0; i<NSUBL; i++)
    {
      logfile << setprecision(35) << Ts[i] << endl;
    }
  /*
  for(int i=0; i<NOBSERVABLES; i++)
    {
      logfile << nobs[i] << " ";
    } 
  */

  logfile << " converged: " << (converged ? "true": "false") << endl;
  

  //  nobs=CalculateOrderPars(newT); // calculate nematic order pars.

  
  
      
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
      logfile << "reached MAXITER=" << par[MAXITER] << " iterations without converging, increase MAXITER!" << endl;
    }
  else
    {      

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
	      logfile << q << " ( " << qoverpi << ")  : " << 1./lv[i].value << endl;  
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

	      realtype factor=newT*NS*0.5;
	      complex<realtype>* start=Kinvq(0,0);
	      
	      if(BINARYOUTFILES)
		{
		  qcorrfile.write((char*) &lineid,sizeof(lineid)); // general format
		  qcorrfile.write((char*) &lattice.nindx_q,sizeof(lattice.nindx_q));
		  for(int i=0; i<lattice.nindx_q; i++)
		    {
		      complex<realtype> value=factor*start[lattice.indx_site_q[i]];
		      qcorrfile.write((char*) &value,sizeof(value));
		    }
		}
	      else
		{
		  //		  qcorrfile << setprecision(16) << lineid << " ";
		  for(int i=0; i<lattice.nindx_q; i++)
		    {
		      complex<realtype> value=factor*start[lattice.indx_site_q[i]];
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
	      vector<realtype> qcorr(Vq);
	      //	      for(int i=0; i<Vq; i++) qcorr[i]=factor/K1[i];
	      
	      if(BINARYOUTFILES)
		{
		  qcorrfile.write((char*) &lineid,sizeof(lineid)); // general format
		  qcorrfile.write((char*) &Vq,sizeof(Vq));
		  qcorrfile.write((char*) &qcorr[0],Vq*sizeof(qcorr[0]));
		}
	      else
		{
		  qcorrfile << setprecision(16) << lineid << " ";
		  for(int i=0; i<Vq; i++){ qcorrfile << qcorr[i] << " ";}
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
		  sigmaefile.write((char*) &Vq,sizeof(Vq));
		  sigmaefile.write((char*) &Sigma[0],Vq*sizeof(Sigma[0]));
		}
	      else
		{
		  sigmaefile << setprecision(16) << lineid << " ";
		  for(int i=0; i<Vq; i++){ sigmaefile << Sigma[i] << " ";}
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
  Simulation(realtype* pars,int i);
  void Run();
 private:
  realtype* param;
  int ic;
  BravaisLattice lattice;
  Couplings couplings;
  Rule rule;

  Driver mysolver;
  vector<NumberList> Deltalist;
  vector<NumberList> Deltastoshowlist;
  vector<bool> Printinfolist;
};


Simulation::Simulation(realtype* pars,int i): param(pars),ic(i),lattice(pars),couplings(pars,NC,NMAT),rule(pars,lattice,couplings),mysolver(pars,lattice,rule),Deltalist(0),Deltastoshowlist(0),Printinfolist(0)
{
  if(TRACE) cout << "Initializing Simulation" << endl;
  
  ifstream parameterfile(PARAMETERFILENAME.c_str());
  if(!parameterfile)
    {
      if(TRACE) 
	cout << "No file " << PARAMETERFILENAME << " found." 
	     << " Using Delta=" << param[DELTA] << endl;
      logfile << "No file " << PARAMETERFILENAME << " found." 
	      << " Using Delta=" << param[DELTA] << endl;
      
      NumberList myval(NSUBL,param[DELTA]);
      Deltalist.push_back(myval); 
      Printinfolist.push_back(true);
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
	      newDelta.v[i]=newval;
	    }
	  Deltalist.push_back(newDelta);
	  Printinfolist.push_back(false);
	}    
    }
  if(TRACE) cout << "Deltalist has " << Deltalist.size() << " entries" << endl;  logfile << "Deltalist has " << Deltalist.size() << " entries" << endl;  
  
  
  ifstream parameterfile2(DELTASTOSHOWFILENAME.c_str());
  if(!parameterfile)
    {
      if(TRACE) 
	cout << "No file " << DELTASTOSHOWFILENAME << " found." 
	     << " Using Delta=" << param[DELTA] << endl;
      logfile << "No file " << DELTASTOSHOWFILENAME << " found." 
	      << " Using Delta=" << param[DELTA] << endl;
      
      NumberList myval(NSUBL,param[DELTA]);
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
	      newDelta.v[i]=newval;
	    }
	  Deltastoshowlist.push_back(newval);
	}
      
    }
  if(TRACE) cout << "Deltastoshowlist has " << Deltastoshowlist.size() << " entries" << endl;  logfile << "Deltastoshowlist has " << Deltastoshowlist.size() << " entries" << endl;  
      
  
  //Search Deltastoshowlist and mark if it is present in Deltalist
  for(int j=0; j<Deltastoshowlist.size(); j++)
    {
      NumberList d=Deltastoshowlist[j];
      bool found=false;
      int indx=0;
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
      for(int i=0; i<Printinfolist.size(); i++)
	{
	  cout << i << " " << Deltalist[i] << " " << "info: " << Printinfolist[i] << endl;
	}
    }
  
  if(TRACE) cout << "Done Initializing Simulation" << endl;
}


void Simulation::Run()
{
  if(TRACE) cout << "Starting Run" << endl;
  for(int i=0; i< Deltalist.size(); i++)
    {
      mysolver.Solve(Deltalist[i],Printinfolist[i]);
    }
  if(TRACE) cout << "Done Run" << endl;
}



#endif //BOND_H
