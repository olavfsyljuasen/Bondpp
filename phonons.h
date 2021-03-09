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

 





typedef std::array<complex<realtype>,3> cmplxcoord;

#if defined XDISPLACEMENTS

const int NDISP  =1; // the dimension of the displacements 
const int NELASTIC=1; // the number of elastic constants.
typedef std::array<realtype,1> voigtstring;
// voigt:                   xx
vector<double> voigt1indx{  0 };
vector<double> voigt2indx{   0};
vector<string> voigtnames{"Exx"};


#elif defined XYDISPLACEMENTS

const int NDISP  =2; // the dimension of the displacements 
const int NELASTIC=3; // the number of elastic constants.
typedef std::array<realtype,3> voigtstring;
// voigt:                   xx,   yy,   xy
vector<double> voigt1indx{  0 ,   1 ,   0 };
vector<double> voigt2indx{   0,    1,    1};
vector<string> voigtnames{"Exx","Eyy","Exy"};

#elif defined XYZDISPLACEMENTS

const int NDISP  =3; // the dimension of the displacements 
const int NELASTIC=6; // the number of elastic constants.
typedef std::array<realtype,6> voigtstring;
// voigt:                   xx,   yy,   zz,   yz,   xz,   xy
vector<double> voigt1indx{  0 ,   1 ,   2 ,   1 ,   0 ,   0 };
vector<double> voigt2indx{   0,    1,    2,    2,    2,    1};
vector<string> voigtnames{"Exx","Eyy","Ezz","Eyz","Exz","Exy"};
#endif

const int NMODE=NSUBL*NDISP; // the number of normal modes



class Phonons
{
 public:
  Phonons();
  ~Phonons(){for(int i=0; i<NSUBL; i++){delete normalmode[i]; }}
    
  double GetOmega(const int q,const int n){return omega(q,n);}
  cmplxcoord GetNormalMode(const int q,const int n,const int i=0){return (*normalmode[i])(q,n);}
  
  voigtstring GetElasticMode(int i){return elasticmode[i];}
  realtype GetElasticeigenvalue(int i){return elasticeigenvalue[i];}

  VecMat<complex<realtype>>& Getf(){return flist;}

  realtype GetSumLogOmegaoverV(){return sumlogomegaoverv;}
 private:  
  const int Vq;
  
  VecMat<realtype> omega;
  
  vector<VecMat<cmplxcoord>* > normalmode;

  vector<voigtstring> elasticmode;
  vector<double> elasticeigenvalue;

  realtype sumlogomegaoverv;
  
  void Initializef();
  VecMat<complex<realtype>> flist;
};



Phonons::Phonons(): Vq(la.SiteqVol()),omega(Vq,NMODE),normalmode(NSUBL),sumlogomegaoverv(0.),flist(Vq,NC,NMODE)
{
  if(TRACE) cout << "Initializing Phonons" << endl;
  // We start with the elastic constants
  // specifying the spings

#if defined SQUAREPHONONS
  int Nsprings=4;
  vector<Triplet> springs(Nsprings);
  springs[0]=Triplet{1,0,0};
  springs[1]=Triplet{0,1,0};
  springs[2]=Triplet{1,1,0};
  springs[3]=Triplet{1,-1,0};
  
  vector<double> couplings(Nsprings);
  couplings[0]=par[ALPHA1];
  couplings[1]=par[ALPHA1];
  couplings[2]=par[ALPHA2];
  couplings[3]=par[ALPHA2];

#elif defined TRIANGULARPHONONS
  int Nsprings=6;
  vector<Triplet> springs(Nsprings);
  springs[0]=Triplet{ 1, 0, 0};
  springs[1]=Triplet{ 0, 1, 0};
  springs[2]=Triplet{-1, 1, 0};
  springs[3]=Triplet{ 1, 1, 0};
  springs[4]=Triplet{-1, 2, 0};
  springs[5]=Triplet{-2, 1, 0};
 
  vector<double> couplings(Nsprings);
  couplings[0]=par[ALPHA1];
  couplings[1]=par[ALPHA1];
  couplings[2]=par[ALPHA1];
  couplings[3]=par[ALPHA2];
  couplings[4]=par[ALPHA2];
  couplings[5]=par[ALPHA2];
  
#elif defined CUBICPHONONS
  int Nsprings=9;
  vector<Triplet> springs(Nsprings);
  springs[0]=Triplet{ 1, 0, 0};
  springs[1]=Triplet{ 0, 1, 0};
  springs[2]=Triplet{ 0, 0, 1};
  springs[3]=Triplet{ 0, 1, 1};
  springs[4]=Triplet{ 0, 1,-1};
  springs[5]=Triplet{ 1, 0, 1};
  springs[6]=Triplet{ 1, 0,-1};
  springs[7]=Triplet{ 1, 1, 0};
  springs[8]=Triplet{ 1,-1, 0};
 
  vector<double> couplings(Nsprings);
  couplings[0]=par[ALPHA1];
  couplings[1]=par[ALPHA1];
  couplings[2]=par[ALPHA1];
  couplings[3]=par[ALPHA2];
  couplings[4]=par[ALPHA2];
  couplings[5]=par[ALPHA2];
  couplings[6]=par[ALPHA2];
  couplings[7]=par[ALPHA2];
  couplings[8]=par[ALPHA2];
#endif
  
  
  
  
  Matrix<eigen_real_type,Dynamic,Dynamic> El(NELASTIC,NELASTIC); // the elastic matrix      
  El.setZero(NELASTIC,NELASTIC);
  
  for(int s=0; s<Nsprings; s++)
    {
      double c=couplings[s];
      Coord  sp=la.rPos(springs[s]);
      double n2=sp.Norm()*sp.Norm();
      
      double spx=sp.x;
      double spy=sp.y;


#if defined XDISPLACEMENTS
      El(0,0) += c * spx*spx * spx*spx / n2;

#elif defined XYDISPLACEMENTS
      El(0,0) += c * spx*spx * spx*spx / n2;
      El(0,1) += c * spx*spx * spy*spy / n2;
      El(0,2) += c * spx*spx * spx*spy / n2;

      El(1,0) += c * spy*spy * spx*spx / n2;
      El(1,1) += c * spy*spy * spy*spy / n2;
      El(1,2) += c * spy*spy * spx*spy / n2;

      El(2,0) += c * spx*spy * spx*spx / n2;
      El(2,1) += c * spx*spy * spy*spy / n2;
      El(2,2) += c * spx*spy * spx*spy / n2;
      
#elif defined XYZDISPLACEMENTS            
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
#endif      
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
  for(int n=0; n<NELASTIC; n++)
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

      ofile << eigenvalue << " : " << endl;
      ofile << evec << endl;
      
      voigtstring thisv;
      for(int i=0; i<NELASTIC; i++) thisv[i]=evec(i);
      elasticmode.push_back(thisv);
      elasticeigenvalue.push_back(eigenvalue);       
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

#if defined SQUAREPHONONS
      /*
      D(0,0)=2*alpha1*(1-cos(q.x))+2*alpha2*(1-cos(q.x)*cos(q.y));
      D(1,1)=2*alpha1*(1-cos(q.y))+2*alpha2*(1-cos(q.x)*cos(q.y));
      D(0,1)=2*alpha2*sin(q.x)*sin(q.y);
      D(1,0)=D(0,1);
      */
      for(int p=0; p<Nsprings; p++)
	{
	  Coord R=la.rPos(springs[p]);
	  realtype R2=scalarproduct(R,R);

	  for(int i=0; i<NDISP; i++)
	    for(int j=0; j<NDISP; j++)
	      D(i,j)+= (R[i]*R[j]/R2)*2*couplings[p]*(1.-cos(q*R));
	}


#elif defined TRIANGULARPHONONS
      for(int p=0; p<Nsprings; p++)
	{
	  Coord R=la.rPos(springs[p]);
	  realtype R2=scalarproduct(R,R);

	  for(int i=0; i<NDISP; i++)
	    for(int j=0; j<NDISP; j++)
	      D(i,j)+= (R[i]*R[j]/R2)*2*couplings[p]*(1.-cos(q*R));
	}
#elif defined CUBICPHONONS
      D(0,0)=2*alpha1*(1-cos(q.x))+2*alpha2*(2-cos(q.x)*cos(q.y)-cos(q.x)*cos(q.z));
      D(1,1)=2*alpha1*(1-cos(q.y))+2*alpha2*(2-cos(q.y)*cos(q.z)-cos(q.y)*cos(q.x));
      D(2,2)=2*alpha1*(1-cos(q.z))+2*alpha2*(2-cos(q.z)*cos(q.x)-cos(q.z)*cos(q.y));
      D(0,1)=2*alpha2*sin(q.x)*sin(q.y);
      D(0,2)=2*alpha2*sin(q.x)*sin(q.z);
      D(1,2)=2*alpha2*sin(q.y)*sin(q.z);
      D(1,0)=D(0,1);
      D(2,0)=D(0,2);
      D(2,1)=D(1,2);
#else
      
#endif
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

	      // we are free to change the eigenvectors by a phase factor, so we choose a convention so
	      // that W_q = W_{-q}^*  THIS IS PROBABLY WRONG WAY TO DO THAT
	      realtype mysign=( real(u_x) < 0 || (real(u_x)==0. && real(u_y) < 0) ? -1.:1.);
	      
	      cmplxcoord thisutslag={mysign*u_x,mysign*u_y,mysign*u_z};

	      (*normalmode[i])(j,n,0)=thisutslag;
	    } 
	}
    }

  
  // enforce W_-q,n = W_q,n^*

  for(int q=0; q<Vq; q++)
    {
      const int mq=la.GetInversionIndx(q);
      if(mq==q) continue;
      
      for(int i=0; i<NSUBL; i++)
	for(int n=0; n<NMODE; n++)
	  {
	    cmplxcoord Wq=(*normalmode[i])( q,n);
	    
	    (*normalmode[i])(mq,n)[0]=conj(Wq[0]);
	    (*normalmode[i])(mq,n)[1]=conj(Wq[1]);
	    (*normalmode[i])(mq,n)[2]=conj(Wq[2]);
	  }
    }


  
  if(TRACE) cout << "Checking that normal/mode eigenvectors are complex conjugate of each other" << endl;
  for(int q=0; q<Vq; q++)
    {
      const int mq=la.GetInversionIndx(q);
      
      for(int i=0; i<NSUBL; i++)
	for(int n=0; n<NMODE; n++)
	  {
	    cmplxcoord Wq=(*normalmode[i])( q,n);
	    cmplxcoord Wmq=(*normalmode[i])(mq,n);
	    
	    if( abs(Wq[0]-conj(Wmq[0])) > 1e-14 || abs(Wq[1]-conj(Wmq[1])) > 1e-14 || abs(Wq[2]-conj(Wmq[2])) > 1e-14)
	      {
		cout << "Warning: q=" << q << " mq=" << mq << " n=" << n << " Wq != Wmq" << endl;
		for(int n=0; n<NMODE; n++)
		  {
		    cmplxcoord Wq=(*normalmode[i])( q,n);
		    cmplxcoord Wmq=(*normalmode[i])(mq,n);
		    cout << "n= " << n;
		    cout << " Wq =(" << Wq[0]  << "," << Wq[1]  << "," << Wq[2]  << ") : "
			 << " Wmq=(" << Wmq[0] << "," << Wmq[1] << "," << Wmq[2] << ")" << endl;
		  }
	      }
	  }
    }

  // compute sum log omega, not counting the q=0 mode.
  for(int q=1; q<Vq; q++)
    for(int n=0; n<NMODE; n++)
      {
	sumlogomegaoverv += log(omega(q,n));
      }
  sumlogomegaoverv /= Vq;
  if(TRACE) cout << "sumlogomega=" << sumlogomegaoverv << endl;
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

  ofstream outfile3("elasticmodes.dat");

  outfile3 << "mode" << " Eigenvalue"  << " Eigenvector( ";
  for(int j=0; j<NELASTIC; j++) outfile3 << voigtnames[j] << " ";
  outfile3 << ")" << endl;
  
  for(int i=0; i<NELASTIC; i++)
    {
      outfile3 << i << "    " << elasticeigenvalue[i] << "                    ( ";
      for(int j=0; j<NELASTIC; j++) outfile3 << elasticmode[i][j] << " ";
      outfile3 << ")" << endl;
    }
  outfile3.close();

  if(TRACE) cout << "Done Initializing Phonons" << endl;


}



void Phonons::Initializef()
{
  if(TRACE) cout << "Initializing f" << endl;
  if(NSUBL != 1){ cout << "ftensor is not implemented for NSUBL != 1, exiting." << endl; exit(1);}

  for(int qi=0; qi<Vq; qi++)
    {
      for(int n=0; n<NMODE; n++)
	{
	  double omega=GetOmega(qi,n);
	  cmplxcoord w=GetNormalMode(qi,n,0);


	  for(int ci=0; ci<NC; ci++)
	    {      
	      Coord c=la.rPos(clist[ci]);
	      cmplxcoord crr={c.x,c.y,c.z};	      

	      complex<realtype> tempf = 0;

	      for(int k=0; k<3; k++){tempf += crr[k]*w[k];}
	      tempf *= complex<realtype>(0,-0.5)*(expi(la.qr(qi,clist[ci]))-1.);
	      tempf *= invsqrtmasses[0]*(1./omega);

	      flist(qi,ci,n)=(qi==0 ? 0: tempf);
	    }
	}
    }
  if(TRACE) cout << "Done Initializing f" << endl;
}






#endif //PHONONS_H
