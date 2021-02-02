#ifndef PHONONS_H
#define PHONONS_H

using namespace std;

#include<complex>

// include eigenvalue routines



//#define EIGENBOOST // use both EIGEN and BOOST                                                                      
#ifdef EIGENBOOST
// Different types of multiprecision backends:                                                                     
#include <boost/multiprecision/cpp_complex.hpp> // complex numbers             
typedef boost::multiprecision::cpp_bin_float_quad   eigen_real_type;
typedef boost::multiprecision::cpp_complex_quad eigen_complex_type;
#include <boost/multiprecision/eigen.hpp>

#else
typedef realtype   eigen_real_type;
typedef complex<realtype> eigen_complex_type;
#endif

#include<Eigen/Dense>
using namespace Eigen;


#include "o3vector.h"

enum voigt{epsxx,epsyy,epszz,epsyz,epsxz,epsxy};

 
const int NDISP  =2; // the dimension of the displacements 
const int NMODE=NSUBL*NDISP; // the number of normal modes
const int NELASTIC=3; // the number of elastic constants.


typedef std::array<complex<realtype>,3> cmplxcoord;

typedef std::array<realtype,6> voigtstring;

class Phonons
{
 public:
  Phonons(double*,BravaisLattice&);
  ~Phonons(){for(int i=0; i<NSUBL; i++){delete normalmode[i]; }}
    
  double GetOmega(const int q,const int n){return omega(q,n);}
  cmplxcoord GetNormalMode(const int q,const int n,const int i=0){return (*normalmode[i])(q,n);}
  
  int GetNumberofElasticModes(){return elasticmode.size();}
  voigtstring GetElasticMode(int i){return elasticmode[i];}
  realtype GetElasticeigenvalue(int i){return elasticeigenvalue[i];}

  VecMat<complex<realtype>>& Getf(){return flist;}

 private:  
  double* par;
  BravaisLattice& la;


  const int Vq;
  
  VecMat<realtype> omega;
  
  vector<VecMat<cmplxcoord>* > normalmode;

  vector<voigtstring> elasticmode;
  vector<double> elasticeigenvalue;


  void Initializef();
  VecMat<complex<realtype>> flist;
};



Phonons::Phonons(double* inpar,BravaisLattice& in_la): par(inpar),la(in_la),Vq(la.SiteqVol()),omega(Vq,NMODE),normalmode(NSUBL),flist(Vq,NC,NMODE)
{
  if(TRACE) cout << "Initializing Phonons" << endl;
  // We start with the elastic constants
  // specifying the spings
  int Nsprings=4;
  vector<Coord> springs(Nsprings);
  springs[0]=Coord(1,0,0);
  springs[1]=Coord(0,1,0);
  springs[2]=Coord(1,1,0);
  springs[3]=Coord(1,-1,0);
  
  vector<double> couplings(Nsprings);
  couplings[0]=par[ALPHA1];
  couplings[1]=par[ALPHA1];
  couplings[2]=par[ALPHA2];
  couplings[3]=par[ALPHA2];
  
  
  Matrix<eigen_real_type,Dynamic,Dynamic> El(6,6); // the elastic matrix      
  El.setZero(6,6);
  
  for(int s=0; s<Nsprings; s++)
    {
      double c=couplings[s];
      Coord  sp=springs[s];
      double n2=sp.Norm()*sp.Norm();
      
      double spx=sp.x;
      double spy=sp.y;
      double spz=sp.z;
      
      El(0,0) += c * spx*spx * spx*spx / n2;
      El(0,1) += c * spx*spx * spy*spy / n2;
      El(0,2) += c * spx*spx * spz*spz / n2;
      El(0,3) += c * spx*spx * spy*spz / n2;
      El(0,4) += c * spx*spx * spx*spz / n2;
      El(0,5) += c * spx*spx * spx*spy / n2;
      
      El(1,0) += c * spy*spy * spx*spx / n2;
      El(1,1) += c * spy*spy * spy*spy / n2;
      El(1,2) += c * spy*spy * spz*spz / n2;
      El(1,3) += c * spy*spy * spy*spz / n2;
      El(1,4) += c * spy*spy * spx*spz / n2;
      El(1,5) += c * spy*spy * spx*spy / n2;
      
      El(2,0) += c * spz*spz * spx*spx / n2;
      El(2,1) += c * spz*spz * spy*spy / n2;
      El(2,2) += c * spz*spz * spz*spz / n2;
      El(2,3) += c * spz*spz * spy*spz / n2;
      El(2,4) += c * spz*spz * spx*spz / n2;
      El(2,5) += c * spz*spz * spx*spy / n2;
      
      El(3,0) += c * spy*spz * spx*spx / n2;
      El(3,1) += c * spy*spz * spy*spy / n2;
      El(3,2) += c * spy*spz * spz*spz / n2;
      El(3,3) += c * spy*spz * spy*spz / n2;
      El(3,4) += c * spy*spz * spx*spz / n2;
      El(3,5) += c * spy*spz * spx*spy / n2;
      
      El(4,0) += c * spx*spz * spx*spx / n2;
      El(4,1) += c * spx*spz * spy*spy / n2;
      El(4,2) += c * spx*spz * spz*spz / n2;
      El(4,3) += c * spx*spz * spy*spz / n2;
      El(4,4) += c * spx*spz * spx*spz / n2;
      El(4,5) += c * spx*spz * spx*spy / n2;
      
      El(5,0) += c * spx*spy * spx*spx / n2;
      El(5,1) += c * spx*spy * spy*spy / n2;
      El(5,2) += c * spx*spy * spz*spz / n2;
      El(5,3) += c * spx*spy * spy*spz / n2;
      El(5,4) += c * spx*spy * spx*spz / n2;
      El(5,5) += c * spx*spy * spx*spy / n2;
    }

  if(TRACE)
    {
      cout << "The elastic matrix" << endl;
      cout << El << endl;
    }

  // diagonalize it, and store the results
  SelfAdjointEigenSolver<Matrix<eigen_real_type,Dynamic,Dynamic> > elsolve(El);
  
  if(TRACE) cout << "Elastic modes:" << endl;

  ofstream ofile("elasticmodes.dat");  
  for(int n=0; n<6; n++)
    {
      eigen_real_type eval= elsolve.eigenvalues()[n]; //eigen sorts eigenvalues,least first                     

#ifdef EIGENBOOST
      realtype eigenvalue=eval.convert_to<realtype>();
#else
      realtype eigenvalue=eval;
#endif

      if(TRACE) cout << "eigenvalue[" << n << "]=" << eigenvalue << endl;

      VectorXd evec=elsolve.eigenvectors().col(n);

      if(TRACE)
	{
	  cout << "eigenvector:" << endl;
	  cout << evec << endl;
	}

      if(eigenvalue >0.)
	{
	  ofile << eigenvalue << " : " << endl;
	  ofile << evec << endl;
	  
	  voigtstring thisv;
	  for(int i=0; i<6; i++) thisv[i]=evec(i);
	  elasticmode.push_back(thisv);
	  elasticeigenvalue.push_back(eigenvalue);
	}
    } 
  ofile.close();  
  
    
  for(int i=0; i<NSUBL; i++)
    { 
      normalmode[i]=new VecMat<cmplxcoord>(Vq,NMODE,1);
    }
  
  
  for(int j=0; j<Vq; j++)
    {
      Coord q=la.qPos(j);
      
      
      // set up the dynamical matrix
      Matrix<eigen_complex_type,Dynamic,Dynamic> D(NMODE,NMODE);      
      D.setZero(NMODE,NMODE);

      // square lattice with horizontal and vertical bonds (alpha1), and diagonal bonds (alpha2). x and y displacements
      // Nsublattice=1, Ndimutslag=2.     
      double alpha1=par[ALPHA1];
      double alpha2=par[ALPHA2];

      D(0,0)=2*alpha1*(1-cos(q.x))+2*alpha2*(1-cos(q.x)*cos(q.y));
      D(0,1)=2*alpha2*sin(q.x)*sin(q.y);
      D(1,1)=2*alpha1*(1-cos(q.y))+2*alpha2*(1-cos(q.x)*cos(q.y));
      D(1,0)=2*alpha2*sin(q.x)*sin(q.y);


      // diagonalize it, and store the results
      SelfAdjointEigenSolver<Matrix<eigen_complex_type,Dynamic,Dynamic> > es(D);

      
      for(int n=0; n<NMODE; n++)
	{
	  eigen_real_type val= es.eigenvalues()[n]; //eigen sorts eigenvalues,least first 

	 
	  /*
#ifdef EIGENBOOST                                 
	  omega(j,n)=val.convert_to<realtype>();
#else
	  omega(j,n)=val;
#endif
	  */

	  omega(j,n)=val;
	  if(TRACE) cout << "eigenvalue[" << n << "]=" << omega(j,n) << endl;
	}
      
      for(int n=0; n<NMODE; n++)
	{
	  if(TRACE) cout << es.eigenvectors().col(n) << endl;
	 
	  for(int i=0; i<NSUBL; i++)
	    {
	      complex<realtype> u_x(0);
	      complex<realtype> u_y(0);
	      complex<realtype> u_z(0);

#ifdef EIGENBOOST
	      if(NDISP>0)
		{
		  eigen_complex_type valx = es.eigenvectors().col(n)[i*NSUBL+0];
		  u_x=valx.convert_to<complex<realtype>>();
		}

	      if(NDISP>1)
		{
		  eigen_complex_type valy = es.eigenvectors().col(n)[i*NSUBL+1];
		  u_y=valy.convert_to<complex<realtype>>();
		}

	      if(NDISP>2)
		{
		  eigen_complex_type valz = es.eigenvectors().col(n)[i*NSUBL+2];
		  u_z=valz.convert_to<complex<realtype>>();
		}
#else
	      if(NDISP>0)
		{
		  eigen_complex_type valx = es.eigenvectors().col(n)[i*NSUBL+0];
		  u_x=valx;
		}

	      if(NDISP>1)
		{
		  eigen_complex_type valy = es.eigenvectors().col(n)[i*NSUBL+1];
		  u_y=valy;
		}

	      if(NDISP>2)
		{
		  eigen_complex_type valz = es.eigenvectors().col(n)[i*NSUBL+2];
		  u_z=valz;
		}
#endif

	      cmplxcoord thisutslag={u_x,u_y,u_z};
	      (*normalmode[i])(j,n,0)=thisutslag;
	    } 
	}
    }
  if(TRACE) cout << "before f" << endl;

  Initializef();


  ofstream outfile("phononenergies.dat");
  for(int j=0; j<Vq; j++)
    {
      Coord q=la.qPos(j);
      outfile << q;
      for(int n=0; n<NMODE; n++)
	{
	  outfile << " " << GetOmega(j,n);
	}
      outfile << endl;  
    } 
  outfile.close();

  ofstream outfile2("normalmodes.dat");
  for(int j=0; j<Vq; j++)
    {
      Coord q=la.qPos(j);
      outfile2 << q;
      for(int n=0; n<NMODE; n++)
	{
	  outfile2 << " | ";
	  for(int i=0; i<NSUBL; i++)
	    {
	      cmplxcoord thismode=GetNormalMode(j,n,i);
	      for(int k=0; k<3; k++) outfile2 << thismode[k] << "  ";
	    }
	}
      outfile2 << endl;  
    } 
  outfile2.close();



  if(TRACE) cout << "Done Initializing Phonons" << endl;


}



void Phonons::Initializef()
{
  if(TRACE) cout << "Initializing f" << endl;
  if(NSUBL != 1){ cout << "ftensor is not implemented for NSUBL != 1, exiting." << endl; exit(1);}

  for(int qi=0; qi<Vq; qi++)
    {
      Coord q=la.qPos(qi);

      for(int n=0; n<NMODE; n++)
	{
	  double omega=GetOmega(qi,n);
	  cmplxcoord w=GetNormalMode(qi,n,0);


	  for(int ci=0; ci<NC; ci++)
	    {      
	      Coord c=la.rPos(clist[ci]);
	      cmplxcoord crr={c.x,c.y,c.z};	      

	      realtype qc=scalarproduct(q,c);


	      complex<realtype> tempf = 0;

	      for(int k=0; k<3; k++){tempf += crr[k]*w[k];}
	      tempf *= complex<realtype>(0,-0.5)*(expi(qc)-1.);
	      tempf *= invsqrtmasses[0]*(1./omega);

	      flist(qi,ci,n)=(qi==0 ? 0: tempf);
	    }
	}
    }
  if(TRACE) cout << "Done Initializing f" << endl;
}






#endif //PHONONS_H
