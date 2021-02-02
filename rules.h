#ifndef RULES_H
#define RULES_H




//#ifdef RANDOMINITIALIZATION
#include "rnddef.h"
//#endif

//#include <cstdlib>
//#include <time.h>
#include<complex>



using namespace std;


typedef complex<realtype> obstype; 

//OBSERVABLES:

#if defined SQUARELATTICE
const int NOBSERVABLES=2;
enum observables{SA,SD};
string NAMESOFOBSERVABLES[NOBSERVABLES]={"sa","sd"};
//#elif defined SIMPLECUBICBRAVAISLATTICE
//const int NOBSERVABLES=3;
//enum observables{A1M,A2M,A2P};
//string NAMESOFOBSERVABLES[NOBSERVABLES]={"a1m","a2m","a2p"};
#elif defined HEXAGONALBRAVAISLATTICE
const int NOBSERVABLES=5;
enum observables{ROT,PHASEII,PHASEIII,SIGMATA,SIGMATD};
string NAMESOFOBSERVABLES[NOBSERVABLES]={"rot","II","III","SA","SD"};
#else
const int NOBSERVABLES=4;
enum observables{SIGMAA,SIGMAD,SIGMAMX,SIGMAMY};
string NAMESOFOBSERVABLES[NOBSERVABLES]={"sa","sd","smx","smy"};
string NAMESOFALPHAS[NOBSERVABLES]={"sa_alpha.dat","sd_alpha.dat","smx_alpha.dat","smy_alpha.dat"};
#endif







class Rule
{
  friend ostream& operator<<(ostream& os,Rule& c){ return c.Write(os);}
  friend istream& operator>>(istream& is,Rule& c){ c.Read(is); return is;}
 public:
  Rule(realtype*,BravaisLattice&,Couplings&);
  ~Rule(){for(int i=0; i<gelptrs.size(); i++){ delete gelptrs[i];}}
  
  void MakeSigma(VecMat<complex<realtype>>& A);
  
  void InitializeJq();
  void InitializeSigma(VecMat<complex<realtype>>& A)
  {
    MakeSigma(A);
  }
  void InitializeObservables(){};
  
  void InitializeJ();
#ifdef PHONONS
  //  void Initializeg();
  
  int GetNumberofElasticModes(){return phonons.GetNumberofElasticModes();}
  void InitializeElasticMode(int i,VecMat<complex<realtype>>& gel);

  VecMat<complex<realtype>>& Getf(){return phonons.Getf();}
#endif

  
  
  vector<obstype>& GetIrrep(const int i){return irrep[i];}
  
  ostream& Write(ostream& os);
  istream& Read(istream& is);
  realtype* GetPars(){return par;}

  realtype* par;
  BravaisLattice& la;

  VecMat<complex<realtype>>& Jr; // defined in couplings.h 
  
  int Vq;
  
  vector<vector<obstype> > irrep;

  VecMat<complex<realtype>> Jq;

#ifdef PHONONS
  Phonons phonons;
  VecMat<complex<realtype>>& g; // defined in couplings.h

  int Nelastic;
  vector<realtype> elasticeigenvalues;
  vector<VecMat<complex<realtype>>*> gelptrs;

#endif



 private:  
  void SetInteraction();  
  void SetIrreps();
  void SetInitialState(realtype);
  
};
 
 


Rule::Rule(realtype* in_par,BravaisLattice& in_la,Couplings& realspacecouplings):par(in_par),la(in_la),Jr(realspacecouplings.J),Vq(la.SiteqVol()),irrep(NOBSERVABLES,vector<obstype>(Vq)), Jq(Vq,NMAT,NMAT)
#ifdef PHONONS
    ,g(realspacecouplings.g),phonons(par,la)
    ,gelptrs(NELASTIC)
#endif
  {
  InitializeObservables();
  if(TRACE) cout << "Initializing Rule " << endl; 

#ifdef PHONONS

  //  Initializeg();
  Nelastic=phonons.GetNumberofElasticModes();

  for(int i=0; i<Nelastic; i++)
    {
      VecMat<complex<realtype>>* gel_ptr=new VecMat<complex<realtype>>(Vq,NMAT,NMAT);
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



  if(TRACE)
    {
      cout << "After initializing Jq:" << Jq << endl;
    }
        
}
 


void Rule::InitializeJq()
{
  for(int ci=0; ci<NC; ci++)   
      for(int qj=0; qj<Vq; qj++)
	for(int i1=0; i1<NMAT; i1++)
	  for(int i2=0; i2<NMAT; i2++)
	    {
	      Jq(qj,i1,i2) += 0.5*Jr(ci,i1,i2)*expi(la.qr(qj,clist[ci]));
	    }
}




#ifdef PHONONS

// voigt:                 xx,yy,zz,yz,xz,xy
vector<double> voigt1indx{0 ,1 ,2 ,1 ,0 ,0 };
vector<double> voigt2indx{ 0, 1, 2, 2, 2, 1};

// called with eps={1,-1,0,0,0,0} gives xx-yy, etc.

void Rule::InitializeElasticMode(int i,VecMat<complex<realtype>>& gel)
{
  const int Vq=la.SiteqVol();
  
  voigtstring eps=phonons.GetElasticMode(i);
  
  for(int vi=0; vi<6; vi++) // run thru Voigt indices
    {
      if(eps[vi]==0.) continue;
      double epsilon=eps[vi];
      
      for(int ci=0; ci<NC; ci++)
	{      
	  Coord c=la.rPos(clist[ci]);
	  
	  for(int qj=0; qj<Vq; qj++)
	    {
	      for(int i1=0; i1<NSUBL; i1++)
		for(int i2=0; i2<NSUBL; i2++)
		  {
		    Coord crr=c+roffset[i2]-roffset[i1];
		    realtype crr2=crr[voigt1indx[vi]]*crr[voigt2indx[vi]];
		    
		    
		    for(int s1=0; s1<NSPIN; s1++)
		      for(int s2=0; s2<NSPIN; s2++)
			{
			  int m1=mindx(s1,i1);
			  int m2=mindx(s2,i2);		
			  
			  gel(qj,m1,m2) += 0.5*epsilon*crr2*g(ci,m1,m2)*expi(la.qr(qj,ci));
			}
		  }
	    }
	}
    }
}  
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

