#ifndef PHONONVERTEX_H
#define PHONONVERTEX_H


class ftensor
{
 public:
  ftensor();
  void Initialize(Phonons& ph);

  VecMat<complex<realtype>> list; 

  complex<realtype> operator()(int j, int c,int n){return (list[j](c,n);}
  
  private:
    int Vq;
};

  ftensor::ftensor():Vq(ph.la.GetSiteqVol()),list(Vq,NC,NMODE){}


void ftensor::Initialize(Phonons& ph)
{
  if(NSUBL != 1){ cout << "ftensor is not implemented for NSUBL != 1, exiting." << endl; exit(1);}

  for(int qi=0; qi<Vq; qi++)
    {
      Coord q=ph.la.qPos(qi);

      for(int n=0; n<NMODE; n++)
	{
	  double omega=ph.GetOmega(j,n);
	  cmplcoord w=ph.GetNormalMode(qi,n,0);


	  for(int ci=0; ci<NC; ci++)
	    {      
	      Coord c=ph.la.rPos(clist[ci]);
	      cmplxcoord crr(c);	      

	      realtype qc=scalarproduct(q,c);


	      complex<realtype> tempf = scalarproduct(crr,w);

	      tempf *= complex<realtype>(0,-0.5)*(Expi(qc)-1.);
	      tempf *= invsqrtmasses[0]*(1./omega);

	      list[qi](ci,n)= tempf;
	    }
	}
    }
}

#endif
