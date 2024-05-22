#ifndef BRAVAISLATTICES_H
#define BRAVAISLATTICES_H


#include<complex>
#include<iomanip>
#include<iostream>
#include<array>

using namespace std;

typedef std::array<int,3> Triplet; // vector of integers that specify positions on lattice in units of lattice vectors

Triplet operator+(const Triplet& lhs,const Triplet& rhs)
{
  Triplet res(lhs);
  for(unsigned int i=0; i<lhs.size(); i++){res[i]=lhs[i]+rhs[i];}
  return res;
}

Triplet operator-(const Triplet& lhs,const Triplet& rhs)
{
  Triplet res(lhs);
  for(unsigned int i=0; i<lhs.size(); i++){res[i]=lhs[i]-rhs[i];}
  return res;
}

Triplet& operator*=(Triplet& lhs,const int it)
{
  for(unsigned int i=0; i<lhs.size(); i++){lhs[i]*=it;}
  return lhs;
}

Triplet& operator+=(Triplet& lhs,const Triplet& rhs)
{
  for(unsigned int i=0; i<lhs.size(); i++){lhs[i]*=rhs[i];}
  return lhs;
}

realtype scalarproduct(const Triplet& a,const Triplet& b)
{
  realtype sum=0.;
  for(unsigned int i=0; i<a.size(); i++){sum += a[i]*b[i];}
  return sum;
}

realtype norm(const Triplet& a)
{
  return sqrt(scalarproduct(a,a));
}

ostream& operator<<(ostream& os,const Triplet& t){os << "(" << t[0] << "," << t[1] << "," << t[2] << ")"; return os;}



bool InRegionOfInterest(Coord c)
{
#if defined XYPLANE
  return c.z==0;
#elif defined XZPLANE
  return c.y==0;
#elif defined YZPLANE
  return c.x==0;
#elif defined XAXIS
  return c.y==0 && c.z==0;
#elif defined YAXIS
  return c.x==0 && c.z==0;
#elif defined ZAXIS
  return c.x==0 && c.y==0;
#elif defined YAXISONPIZEROZERO
  return c.x==PI && c.z==0;
#elif defined ZAXISONPIZEROZERO
  return c.x==PI && c.y==0;
#elif defined XAXISONZEROPIZERO
  return c.y==PI && c.z==0;
#elif defined ZAXISONZEROPIZERO
  return c.x==0 && c.y==PI;
#else
  return true;  // just pick all points 
#endif
}

// return true if a is an integer (closer to the integer than epsilon) 
bool IsInt(realtype a)
{
  realtype fa=fabs(a);
  return (fabs(fa-int(fa+0.5)) <= 1.e-10);
}




class BravaisLattice
{
 public:
  BravaisLattice(realtype* par);
  ~BravaisLattice(){if(TRACE) cout << "deleting BravaisLattice\n";}

  const int N1;
  const int N2;
  const int N3;

  const double invN1;
  const double invN2;
  const double invN3;

  
 private:
  const int d;

  const int Nr; // real space number of sites
  const int Nq; // q-space number of sites
  
 public:

  //  int Nsites() const {return v;} 
  int D() const {return d;}
  Coord rPos(const int i){return site_r[i];}
  Coord qPos(const int i){return site_q[i];}

  Coord rPos(const Triplet i){return i[0]*a1+i[1]*a2+i[2]*a3;}

  int NrSites() const {return Nr;} 
  vector<int> SiterDims() const {return site_rdims;}
  int NqSites() const {return Nq;}
  vector<int> SiteqDims() const {return site_qdims;}
  
  int qAdd(const int k1,const int k2);  // return matrix entry for k1+k2
  int qSub(const int k1,const int k2);  // return matrix entry for k1-k2


  int rAdd(const int s,const Triplet v);
  int rSub(const int s,const Triplet v);

  realtype UnitCellVolume(){return unitcellvolume;} // volume of unit cell
#ifdef PRESERVESYMMETRY
  int TransformationPeriod;
  vector<int> TransformationTable;
#endif

  Coord GetLatticeVector(int i){return ( i==1 ? a1 : i==2 ? a2: a3);}
  
 private:
  const Coord a1;
  const Coord a2;
  const Coord a3;

  const realtype unitcellvolume;
  const Coord origo;

  Coord b1;
  Coord b2;
  Coord b3;

  vector<int> site_rdims;
  vector<Coord> site_r;
  
  vector<int> site_qdims;
  vector<Coord> site_q;


  Triplet UsePBC(const Triplet); // Use PBC to put coordinates inside lattice
  Triplet GetTriplet(const int k);
  int GetIndx(const Triplet v);


 public:

  int nindx_r;
  int nindx_q;

  realtype qr(const Triplet&, const Triplet& );
  realtype qr(const int&    , const int&     );
  realtype qr(const int&    , const Triplet& );
  realtype qr(const Triplet&, const int&     );
  int GetInversionIndx(const int s);

#ifdef REDUCEDOUTPUT
  vector<int> indx_site_r;
  vector<int> indx_site_q;
#endif

  int nselectedqpts;
  vector<int> selectedqpts;

  int indx_most_remote;

  vector<int> ConstructTransformationTable(Coord,Coord,Coord);
};



BravaisLattice::BravaisLattice(realtype* par):
N1(static_cast<int>(par[NX])),
  N2(static_cast<int>(par[NY])),
  N3(static_cast<int>(par[NZ])),
  invN1(1./N1),
  invN2(1./N2),
  invN3(1./N3),
  d( (N1>1)+(N2>1)+(N3>1) ),
#if defined QSPACEDIRECTLATTICE
//  vr(N3*N2*(N1/2+1)),
  Nr(N3*N2*N1),
  Nq(N1*N2*N3),
#else
  Nr(N3*N2*N1),
  Nq(N3*N2*(N1/2+1)),
#endif    
  
// real space basis vectors:
  
#ifdef SIMPLECUBICBRAVAISLATTICE
  a1(1.,0.,0.),
  a2(0.,1.,0.),
  a3(0.,0.,1.),
#elif defined BCCBRAVAISLATTICE
  a1(-0.5, 0.5, 0.5),
  a2( 0.5,-0.5, 0.5),
  a3( 0.5, 0.5,-0.5),
#elif defined FCCBRAVAISLATTICE
  a1(0. ,0.5,0.5),
  a2(0.5,0.0,0.5),
  a3(0.5,0.5,0. ),
#elif defined HEXAGONALBRAVAISLATTICE
  a1(1.,0.,0.),
  a2(0.5,SQRTTHREEOVERTWO,0.),
  a3(0.,0.,1.),
#endif
  
  unitcellvolume(scalarproduct(a1,crossproduct(a2,a3))),
  origo(0.,0.,0.),   
  
  site_rdims(d),
  site_r(Nr),
  
  site_qdims(d),
  site_q(Nq),
  
  nindx_r(0),  
  nindx_q(0)

#ifdef REDUCEDOUTPUT
  ,indx_site_r(Nr)
  ,indx_site_q(Nq)
#endif

  ,nselectedqpts(0)
  ,selectedqpts(MAXNSELECTEDQPTS)
  ,indx_most_remote(0)
{
  logfile << "initializing BravaisLattice" << endl;
  
  if(TRACE) cout << "d: " << d << endl;
  if(TRACE) cout << "Nr: " << Nr << " Nq: " << Nq  << endl;
 
  logfile << "d: " << d << endl;
  logfile << "Nr: " << Nr << " Nq: " << Nq  << endl;

  // initialize all sizes with 1
  int N1r=1;
  int N2r=1;
  int N3r=1;
  int N1q=1;
  int N2q=1;
  int N3q=1;
  
#ifdef QSPACEDIRECTLATTICE
  /*
  if(d==1){N1q=N1;                       N1r=N1/2+1;                 }
  if(d==2){N1q=N1;     N2q=N2;           N1r=N1/2+1; N2r=N2;         }
  if(d==3){N1q=N1;     N2q=N2; N3q=N3;   N1r=N1/2+1; N2r=N2; N3r=N3; }
  */
  if(d==1){N1q=N1;                       N1r=N1;                     }
  if(d==2){N1q=N1;     N2q=N2;           N1r=N1; N2r=N2;             }
  if(d==3){N1q=N1;     N2q=N2; N3q=N3;   N1r=N1; N2r=N2; N3r=N3;     }
#else
  if(d==1){N1q=N1/2+1;                   N1r=N1;                     }
  if(d==2){N1q=N1/2+1; N2q=N2;           N1r=N1;     N2r=N2;         }
  if(d==3){N1q=N1/2+1; N2q=N2; N3q=N3;   N1r=N1;     N2r=N2; N3r=N3; }
#endif


  logfile << "dimensions: " << endl;
  logfile << "N1=" << N1 << " N2=" << N2 << " N3=" << N3 << endl;
  logfile << "rspace: " << N1r << "x" << N2r << "x" << N3r << endl;
  logfile << "qspace: " << N1q << "x" << N2q << "x" << N3q << endl;
  

  
  if(d==1){site_rdims[0]=N1r;}
  if(d==2){site_rdims[0]=N2r; site_rdims[1]=N1r;}
  if(d==3){site_rdims[0]=N3r; site_rdims[1]=N2r; site_rdims[2]=N1r;}

  // reciprocal basisvectors
  realtype invunitvolume=1./unitcellvolume;
  b1= TWOPI*invunitvolume*crossproduct(a2,a3);
  b2= TWOPI*invunitvolume*crossproduct(a3,a1);
  b3= TWOPI*invunitvolume*crossproduct(a1,a2);

  logfile << "q-space basis vectors: " << endl;
  logfile << "b1= " << b1 << endl;
  logfile << "b2= " << b2 << endl;
  logfile << "b3= " << b3 << endl;
  
  

  if(TRACE) cout << "Starting to record positions" << endl;
  if(TRACE) cout << "N1=" << N1 << " N2=" << N2 << " N3=" << N3 << endl;
  // set the coordinate positions and phase

  int s=-1;
  for(int i3=0; i3<N3r; i3++)
    for(int i2=0; i2<N2r; i2++)
      for(int i1=0; i1<N1r; i1++)
	{
	  Coord position=origo+i1*a1+i2*a2+i3*a3;

	  site_r[++s]=position;

	  if(TRACE) cout << "site: " << s 
			 << " int coord: (" 
			 << i1 << " " << i2 << " " << i3 << ") " << position << " \n"; 
	}

  indx_most_remote=((N3/2)*N2r+N2/2)*N1r+N1/2; 
  logfile << "setting indx_most_remote=" << indx_most_remote << endl;


#ifdef REDUCEDOUTPUT
  for(int i=0; i<Nr; i++)
    {
      if( InRegionOfInterest(site_r[i]) ){indx_site_r[nindx_r++]=i;}
    }
#endif


   if(PRINTRPTS)
    {
      ofstream file_siter( SITERPTS.c_str());
#ifdef REDUCEDOUTPUT
      for(int i=0; i<nindx_r; i++) file_siter << site_r[indx_site_r[i]] << endl; 
#else
      for(int i=0; i<Nr; i++) file_siter << site_r[i] << endl; 
#endif
    }

  // initializing the site q-coordinates
  if(d==1){site_qdims[0]=N1q;}
  if(d==2){site_qdims[0]=N2q; site_qdims[1]=N1q;}
  if(d==3){site_qdims[0]=N3q; site_qdims[1]=N2q; site_qdims[2]=N1q;}

  // specify the q-coordinates of the entries in the
  // fourier-transform matrix;

  int i=0;
  if(d==1)
    { 
      for(int q1=0; q1< N1q; q1++) 
	site_q[i++]=b1*static_cast<realtype>(q1*invN1);
    }
  if(d==2)
    {
      for(int q2=0; q2< N2q; q2++)
	for(int q1=0; q1< N1q; q1++) 
	  site_q[i++]=
	    b1*static_cast<realtype>(q1*invN1)+
	    b2*static_cast<realtype>(q2*invN2);
    }
  if(d==3)
    {
      for(int q3=0; q3< N3q; q3++)
	for(int q2=0; q2< N2q; q2++)
	  for(int q1=0; q1< N1q; q1++) 
	    site_q[i++]=
	    b1*static_cast<realtype>(q1*invN1)+
	    b2*static_cast<realtype>(q2*invN2)+
	    b3*static_cast<realtype>(q3*invN3);
    }
  

#ifdef REDUCEDOUTPUT
  for(int i=0; i<Nq; i++)
    {
      if( InRegionOfInterest(site_q[i]) ){indx_site_q[nindx_q++]=i;}
    }
#endif


  ifstream selectedqptsfile( SELECTEDQPTSFILENAME.c_str());
  
  while(selectedqptsfile)
    {
      realtype qx,qy,qz;
      selectedqptsfile >> qx >> qy >> qz;
      if(selectedqptsfile)
	{
	  cout << qx << " " << qy << " " << qz << endl;
	  
	  for(i=0; i<Nq; i++)
	    if( abs(site_q[i].x - qx)< 1.e-6 &&
		abs(site_q[i].y - qy)< 1.e-6 &&
		abs(site_q[i].z - qz)< 1.e-6){ selectedqpts[nselectedqpts++]=i;}
	}
    }
  logfile << "there are " << nselectedqpts << " selected q pts" << endl;


  if(PRINTQPTS)
    {
      ofstream file_siteq( SITEQPTS.c_str());
#ifdef REDUCEDOUTPUT
      for(int i=0; i<nindx_q; i++) file_siteq << setprecision(16) << site_q[indx_site_q[i]] << endl; 
#else
      for(int i=0; i<Nq; i++) file_siteq << setprecision(16) << site_q[i] << endl; 
#endif
    }


  // Transformations
  // specified by new transformed basis vectors and the maximal period 
#ifdef PRESERVESYMMETRY

#ifdef SIMPLECUBICBRAVAISLATTICE

#ifdef INVERSIONSYMMETRY
  Coord ap1(-1.,0.,0.);
  Coord ap2(0.,-1.,0.);
  Coord ap3(0.,0.,-1.);

  TransformationPeriod=2;
  TransformationTable=ConstructTransformationTable(ap1,ap2,ap3);

#elif defined C4ROTATIONSYMMETRYABOUTZAXIS
  Coord ap1(0.,1.,0.);
  Coord ap2(-1.,0,0.);
  Coord ap3(0.,0.,1.);

  TransformationPeriod=4;
  TransformationTable=ConstructTransformationTable(ap1,ap2,ap3);

#elif defined REFLECTIONSYMMETRYABOUTYAXIS
  Coord ap1(-1.,0.,0.);
  Coord ap2(0.,1.,0.);
  Coord ap3(0.,0.,1.);

  TransformationPeriod=2;
  TransformationTable=ConstructTransformationTable(ap1,ap2,ap3);

#elif defined REFLECTIONSYMMETRYABOUTXAXIS
  Coord ap1(1.,0.,0.);
  Coord ap2(0.,-1.,0.);
  Coord ap3(0.,0.,1.);

  TransformationPeriod=2;
  TransformationTable=ConstructTransformationTable(ap1,ap2,ap3);

#elif defined REFLECTIONSYMMETRYABOUTXYDIAGONAL
  Coord ap1(0.,1.,0.);
  Coord ap2(1.,0.,0.);
  Coord ap3(0.,0.,1.);

  TransformationPeriod=2;
  TransformationTable=ConstructTransformationTable(ap1,ap2,ap3);

#elif defined REFLECTIONSYMMETRYABOUTNEGATIVEXYDIAGONAL
  Coord ap1(0.,-1.,0.);
  Coord ap2(-1.,0.,0.);
  Coord ap3(0.,0.,1.);
  
  TransformationPeriod=2;
  TransformationTable=ConstructTransformationTable(ap1,ap2,ap3);  
  
#elif defined INVERSIONSYMMETRY
  Coord ap1(-1.,0.,0.);
  Coord ap2(0.,-1.,0.);
  Coord ap3(0.,0.,-1.);
  
  TransformationPeriod=2;
  TransformationTable=ConstructTransformationTable(ap1,ap2,ap3);  
  
#endif // SIMPLECUBICBRAVAISLATTICE

#elif defined HEXAGONALBRAVAISLATTICE

#ifdef INVERSIONSYMMETRY
  Coord ap1(-1.,0.,0.);
  Coord ap2(-0.5,-SQRTTHREEOVERTWO,0.);
  Coord ap3(0.,0.,-1.);

  TransformationPeriod=2;
  TransformationTable=ConstructTransformationTable(ap1,ap2,ap3);


#elif defined C6ROTATIONSYMMETRYABOUTZAXIS
  Coord ap1(0.5,SQRTTHREEOVERTWO,0.);
  Coord ap2(-0.5,SQRTTHREEOVERTWO,0.);
  Coord ap3(0.,0.,1.);

  TransformationPeriod=6;
  TransformationTable=ConstructTransformationTable(ap1,ap2,ap3);

#elif defined C3ROTATIONSYMMETRYABOUTZAXIS
  Coord ap1(-0.5,SQRTTHREEOVERTWO,0.);
  Coord ap2(-1.0,0,0.);
  Coord ap3(0.,0.,1.);

  TransformationPeriod=3;
  TransformationTable=ConstructTransformationTable(ap1,ap2,ap3);    

#elif defined REFLECTIONSYMMETRYABOUTXAXIS
  Coord ap1(1.0,0,0.);
  Coord ap2(0.5,-SQRTTHREEOVERTWO,0.);
  Coord ap3(0.,0.,1.);

  TransformationPeriod=2;
  TransformationTable=ConstructTransformationTable(ap1,ap2,ap3);    

#elif defined REFLECTIONSYMMETRYABOUTYAXIS
  Coord ap1(-1.0,0,0.);
  Coord ap2(-0.5,SQRTTHREEOVERTWO,0.);
  Coord ap3(0.,0.,1.);

  TransformationPeriod=2;
  TransformationTable=ConstructTransformationTable(ap1,ap2,ap3);    

#endif

#endif // HEXAGONALBRAVAISLATTICE  
  
#endif //PRESERVESYMMETRY
  
  
  if(TRACE)
    {
      for(int i=0; i<Nq; i++)
	{
	  Coord ic=qPos(i);
	  int ni=GetInversionIndx(i);
	  Coord nic=qPos(ni);
	  cout << "site: " << i 
	       << " pos=" << ic
	       << " inversion indx: " << ni 
	       << " pos=" << nic << endl;
	}
    }
}

inline realtype BravaisLattice::qr(const Triplet& qivec,const Triplet& rivec)
{
  return TWOPI*(   qivec[0]*rivec[0]*invN1 
		 + qivec[1]*rivec[1]*invN2
		 + qivec[2]*rivec[2]*invN3);
}

inline realtype BravaisLattice::qr(const int& qi,const int& ri)
{
  return qr(GetTriplet(qi),GetTriplet(ri));
}

inline realtype BravaisLattice::qr(const int& qi,const Triplet& rtri)
{
  return qr(GetTriplet(qi),rtri);
}

inline realtype BravaisLattice::qr(const Triplet& qtri,const int& ri)
{
  return qr(qtri,GetTriplet(ri));
}


// this needs to be tested. Perhaps a bit dangerous...
inline int BravaisLattice::GetInversionIndx(const int s)
{
  const int xc = s%N1;
  const int r  = s/N1;
  const int yc = r%N2;
  const int zc = r/N2;

  return (-xc+N1)%N1 + N1*( (-yc+N2)%N2 + N2*( (-zc+N3)%N3));
}


inline Triplet BravaisLattice::UsePBC(const Triplet v)
{
  Triplet w={(v[0]+N1)%N1,(v[1]+N2)%N2,(v[2]+N3)%N3};

  return w;
}


inline int BravaisLattice::GetIndx(const Triplet v)
{
  return v[0]+N1*(v[1]+N2*v[2]);
}


inline Triplet BravaisLattice::GetTriplet(const int s)
{
  const int r  = s/N1;
  return Triplet{s%N1,r%N2,r/N2};
}


int BravaisLattice::rAdd(const int s,const Triplet v)
{
  const int r  = s/N1;
  return (s%N1+v[0]+N1)%N1 + N1*( (r%N2+v[1]+N2)%N2 + N2*( (r/N2+v[2]+N3)%N3));
}


int BravaisLattice::rSub(const int s,const Triplet v)
{
  const int r  = s/N1;
  return (s%N1-v[0]+N1)%N1 + N1*( (r%N2-v[1]+N2)%N2 + N2*( (r/N2-v[2]+N3)%N3));
}




inline int BravaisLattice::qAdd(const int ka,const int kb)
{
  Triplet va=GetTriplet(ka);
  Triplet vb=GetTriplet(kb);
  
  return GetIndx(UsePBC(va+vb));
}

inline int BravaisLattice::qSub(const int ka,const int kb)
{
  Triplet va=GetTriplet(ka);
  Triplet vb=GetTriplet(kb);

  return GetIndx(UsePBC(va-vb));
}


vector<int> BravaisLattice::ConstructTransformationTable(Coord ap1,Coord ap2,Coord ap3)
{
  vector<int> T(Nq);
  // reciprocal basisvectors in transformed basis
  realtype invunitvolume=1./scalarproduct(ap1,crossproduct(ap2,ap3));
  Coord bp1= TWOPI*invunitvolume*crossproduct(ap2,ap3);
  Coord bp2= TWOPI*invunitvolume*crossproduct(ap3,ap1);
  Coord bp3= TWOPI*invunitvolume*crossproduct(ap1,ap2);  

  Coord Tq;

  for(int s=0; s<Nq; s++)
    {
      Triplet q=GetTriplet(s);

      Tq=(d>=1 ? bp1*static_cast<realtype>(q[0]*1./N1):0.) + (d>=2 ? bp2*static_cast<realtype>(q[1]*1./N2):0.)  + (d>=3 ? bp3*static_cast<realtype>(q[2]*1./N3):0.);
 
      //      cout << "s= " << s << " q=" << q[0] << " " << q[1] << " " << q[2] << " site: (" << site_q[s] << ")" << endl;
      // check Tq is equal to one of the original points mod a reciprocal lattice vector

      realtype a1dotTq= scalarproduct(a1,Tq)/TWOPI;
      realtype a2dotTq= scalarproduct(a2,Tq)/TWOPI;
      realtype a3dotTq= scalarproduct(a3,Tq)/TWOPI;

      //      cout << "Tq = (" << Tq << ") " << a1dotTq << " " << a2dotTq << " " << a3dotTq << endl;
      
      int j=0;
      bool done=false;
      
      while(!done && j<Nq)
	{
	  Triplet oq=GetTriplet(j);

	  //	  cout << oq[0] << " " << oq[1] << " " << oq[2] << endl;
	  
	  bool isinteger= (d>=1 ? IsInt(a1dotTq-oq[0]*1./N1):1) && (d>=2 ? IsInt(a2dotTq-oq[1]*1./N2):1) && (d>=3 ? IsInt(a3dotTq-oq[2]*1./N3):1);
	  
	  //	  cout << a1dotTq-oq[0]*1./N1 << " " << a2dotTq-oq[1]*1./N2 << " " << a3dotTq-oq[2]*1./N3 << " integer:" << isinteger << endl;
	  
	  if( isinteger)
	    {
	      T[s]=j;
	      done=true;
	    }
	  j++;
	}
    }
  ofstream tfile("symmetrytable.dat"); 
  for(int i=0; i<Nq; i++) tfile << i << " " << T[i] << endl;
  tfile.close();
  
  return T;
}  


// This routine should probabily be in the matrix section


#endif
