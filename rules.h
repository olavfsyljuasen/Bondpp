#ifndef RULES_H
#define RULES_H




//#ifdef RANDOMINITIALIZATION
#include "rnddef.h"
//#endif

//#include <cstdlib>
//#include <time.h>
#include<complex>



using namespace std;


class Rule
{
  friend ostream& operator<<(ostream& os,Rule& c){ return c.Write(os);}
  friend istream& operator>>(istream& is,Rule& c){ c.Read(is); return is;}
 public:
  Rule(Couplings&);
  ~Rule()
    {
#ifdef PHONONS
      for(unsigned int i=0; i<gelptrs.size(); i++){ delete gelptrs[i];}
#endif
    }
  
  void MakeSigma(VecMat<complex<realtype>,NMAT,NMAT>& A);
  
  void InitializeJq();
  void InitializeSigma(VecMat<complex<realtype>,NMAT,NMAT>& A)
  {
    MakeSigma(A);
  }
  void InitializeObservables(){};
  
  void InitializeJ();
#ifdef PHONONS
  //  void Initializeg();
  
  void InitializeElasticMode(int i,VecMat<complex<realtype>,NMAT,NMAT>& gel);

  VecMat<complex<realtype>,NC,NMODE>& Getf(){return phonons.Getf();}
#endif

  
  
  vector<obstype>& GetIrrep(const int i){return irrep[i];}
  
  ostream& Write(ostream& os);
  istream& Read(istream& is);
  realtype* GetPars(){return par;}

  VecMat<complex<realtype>,NMAT,NMAT>& Jr; // defined in couplings.h 
  
  int Nq;
  
  vector<vector<obstype> > irrep;

  VecMat<complex<realtype>,NMAT,NMAT> Jq;

#ifdef PHONONS
  Phonons phonons;
  VecMat<complex<realtype>,NMAT,NMAT>& g; // defined in couplings.h

  vector<realtype> elasticeigenvalues;
  vector<VecMat<complex<realtype>,NMAT,NMAT>* > gelptrs;

  realtype GetSumLogOmegaoverV(){return phonons.GetSumLogOmegaoverV();}
#endif



 private:  
  void SetInteraction();  
  void SetIrreps();
  void SetInitialState(realtype);
  
};
 
 


Rule::Rule(Couplings& realspacecouplings):Jr(realspacecouplings.J),Nq(la.NqSites()),irrep(NOBSERVABLES,vector<obstype>(Nq)), Jq(Nq)
#ifdef PHONONS
  ,phonons()
  ,g(realspacecouplings.g)
  ,gelptrs(NELASTIC)
#endif
  {
  InitializeObservables();
  if(TRACE) cout << "Initializing Rule " << endl; 

#ifdef PHONONS

  for(int i=0; i<NELASTIC; i++)
    {
      VecMat<complex<realtype>,NMAT,NMAT>* gel_ptr=new VecMat<complex<realtype>,NMAT,NMAT>(Nq);
      gel_ptr->SetToZero();
      InitializeElasticMode(i,*gel_ptr);
      gelptrs[i]=gel_ptr;
      
      elasticeigenvalues.push_back(phonons.GetElasticeigenvalue(i));
    }
#endif

  if(TRACE)
    {
      cout << "Jr:" << Jr << endl;
    }


  InitializeJq();        
}
 

#ifdef CPOSITIVE

void Rule::InitializeJq()
{
  for(int ci=0; ci<NC; ci++)   
      for(int qj=0; qj<Nq; qj++)
	for(int i1=0; i1<NMAT; i1++)
	  for(int i2=0; i2<NMAT; i2++)
	    {
	      Jq(qj,i1,i2) += 0.5*(Jr(ci,i1,i2)*expi(la.qr(qj,clist[ci]))+Jr(ci,i2,i1)*conj(expi(la.qr(qj,clist[ci]))));
	    }
}

#else

void Rule::InitializeJq()
{
  for(int ci=0; ci<NC; ci++)   
      for(int qj=0; qj<Nq; qj++)
	for(int i1=0; i1<NMAT; i1++)
	  for(int i2=0; i2<NMAT; i2++)
	    {
	      Jq(qj,i1,i2) += 0.5*Jr(ci,i1,i2)*expi(la.qr(qj,clist[ci]));
	    }
}

#endif



#ifdef PHONONS

// called with eps={1,-1,0,0,0,0} gives xx-yy, etc.

#ifdef CPOSITIVE
void Rule::InitializeElasticMode(int i,VecMat<complex<realtype>,NMAT,NMAT>& gel)
{
  if(TRACE) cout << "Starting InitializeElasticMode " << i << endl;
  const int Nq=la.NqSites();
  
  voigtstring u=phonons.GetElasticMode(i);
  
  for(int ci=0; ci<NC; ci++)
    {      
      Coord c=la.rPos(clist[ci]);
      
      for(int qj=0; qj<Nq; qj++)
	{
	  for(int vi=0; vi<NELASTIC; vi++) // run thru all Voigt indices, same as NELASTIC
	    {
	      for(int i1=0; i1<NSUBL; i1++)
		for(int i2=0; i2<NSUBL; i2++)
		  {
		    Coord crr=c+roffset[i2]-roffset[i1];
		    //		    realtype crr2=crr[voigt1indx[vi]]*crr[voigt2indx[vi]]*(voigt1indx[vi] != voigt2indx[vi] ? 2:1);
		    realtype crr2=crr[voigt1indx[vi]]*crr[voigt2indx[vi]];
	      
		    for(int s1=0; s1<NSPIN; s1++)
		      for(int s2=0; s2<NSPIN; s2++)
			{
			  int m1=mindx(s1,i1);
			  int m2=mindx(s2,i2);		
			  
			  gel(qj,m1,m2) += 0.5*u[vi]*crr2*g(ci,m1,m2)*expi(la.qr(qj,clist[ci]));
			}
		  }

	      for(int i1=0; i1<NSUBL; i1++)
		for(int i2=0; i2<NSUBL; i2++)
		  {
		    Coord crr=c+roffset[i1]-roffset[i2];
		    //		    realtype crr2=crr[voigt1indx[vi]]*crr[voigt2indx[vi]]*(voigt1indx[vi] != voigt2indx[vi] ? 2:1);
		    realtype crr2=crr[voigt1indx[vi]]*crr[voigt2indx[vi]];
		    for(int s1=0; s1<NSPIN; s1++)
		      for(int s2=0; s2<NSPIN; s2++)
			{
			  int m1=mindx(s1,i1);
			  int m2=mindx(s2,i2);		
			  
			  gel(qj,m1,m2) += 0.5*u[vi]*crr2*g(ci,m2,m1)*conj(expi(la.qr(qj,clist[ci])));
			}
		  }
	    }
	}
    }
  if(TRACE) SanityCheck(gel,"gel, at end of InitializeElasticMode");
  if(TRACE) cout << "Finished InitializeElasticMode" << endl;
}
#else
void Rule::InitializeElasticMode(int i,VecMat<complex<realtype>>& gel)
{
  if(TRACE) cout << "Starting InitializeElasticMode " << i << endl;
  const int Nq=la.NqSites();
  
  voigtstring u=phonons.GetElasticMode(i);
  
  for(int ci=0; ci<NC; ci++)
    {      
      Coord c=la.rPos(clist[ci]);
      
      for(int qj=0; qj<Nq; qj++)
	{
	  for(int vi=0; vi<NELASTIC; vi++) // run thru all Voigt indices, same as NELASTIC
	    {
	      
	      for(int i1=0; i1<NSUBL; i1++)
		for(int i2=0; i2<NSUBL; i2++)
		  {
		    Coord crr=c+roffset[i2]-roffset[i1];
		    //		    realtype crr2=crr[voigt1indx[vi]]*crr[voigt2indx[vi]]*(voigt1indx[vi] != voigt2indx[vi] ? 2:1);
		    realtype crr2=crr[voigt1indx[vi]]*crr[voigt2indx[vi]]; // do not overcount
		    
		    for(int s1=0; s1<NSPIN; s1++)
		      for(int s2=0; s2<NSPIN; s2++)
			{
			  int m1=mindx(s1,i1);
			  int m2=mindx(s2,i2);		
			  
			  gel(qj,m1,m2) += 0.5*u[vi]*crr2*g(ci,m1,m2)*expi(la.qr(qj,clist[ci]));
			}
		  }
	    }
	}
    }
  if(TRACE) SanityCheck(gel,"gel, at end of InitializeElasticMode");
  if(TRACE) cout << "Finished InitializeElasticMode" << endl;
}  
#endif
#endif



ostream& Rule::Write(ostream& os)
{
  //  os.write( (char*) &Nbonds,sizeof(Nbonds));
  //  os.write( (char*) &lambda_total,sizeof(lambda_total));
  return os;
}


istream& Rule::Read(istream& is)
{
  //  is.read( (char*) &Nbonds,sizeof(Nbonds));
  //  is.read( (char*) &lambda_total,sizeof(lambda_total));
  return is;
}




#endif 

