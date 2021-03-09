#ifndef MODELDEF_H
#define MODELDEF_H

enum spindirections{SX,SY,SZ};

using namespace std;




#ifdef SQUAREPHONONS

#ifdef FAKEHEISENBERG
const int NSPIN= 1;
const int NFAKESPINTRACE=3; // a compensating factor for treating Heisenberg models with a single spin comp.
#else
const int NSPIN= 3;
const int NFAKESPINTRACE=1; 
#endif
const int NSUBL= 1;



#ifdef CPOSITIVE
const int NC=6;
const int NN=2;
const int NNN=4;
vector<Triplet> clist{{ 1,  0,  0 }, { 0,  1,  0 }, { 1,  1,  0 }, { 1, -1,  0 }, { 2, 0, 0}, { 0, 2, 0}};
#else
const int NC=12;
const int NN=4;
const int NNN=8;
vector<Triplet> clist{{ 1,  0,  0 }, {-1,  0,  0 }, { 0,  1,  0 }, { 0, -1,  0 }, { 1,  1,  0 }, {-1, -1, 0}, { 1, -1,  0 }, {-1,  1,  0 },{ 2,0,0},{-2,0,0},{0,2,0},{0,-2,0}};
#endif
//vector<int> minusc{1,0,3,2,5,4,7,6,9,8,11,10}; // refers to indx in clist
vector<Coord> roffset={Coord(0,0,0)};
vector<double> invsqrtmasses={1.};

// global routines to extract sublattice and spin indices.

int mindx(const int spin,const int subl=0){ return spin*NSUBL+subl;}
int spin(const int m){return m/NSUBL;}
int subl(const int m){return m%NSUBL;}

// specific coupling rules for squarephonons
void Couplings::InitializeJ()
{
#ifdef FAKEHEISENBERG
  for(int c=0; c<NN; c++)
    {
      J(c,mindx(SX,0),mindx(SX,0)) = par[J1X];
    } 

  for(int c=NN; c<NNN; c++)
    {
      J(c,mindx(SX,0),mindx(SX,0)) = par[J2X];
    } 

  for(int c=NNN; c<NC; c++)
    {
      J(c,mindx(SX,0),mindx(SX,0)) = par[J3X];
    } 
#else
  
  for(int c=0; c<NN; c++)
    {
      J(c,mindx(SX,0),mindx(SX,0)) = par[J1X];
      J(c,mindx(SY,0),mindx(SY,0)) = par[J1Y];
      J(c,mindx(SZ,0),mindx(SZ,0)) = par[J1Z]; 
    } 

  for(int c=NN; c<NNN; c++)
    {
      J(c,mindx(SX,0),mindx(SX,0)) = par[J2X];
      J(c,mindx(SY,0),mindx(SY,0)) = par[J2Y];
      J(c,mindx(SZ,0),mindx(SZ,0)) = par[J2Z];
    } 

  for(int c=NNN; c<NC; c++)
    {
      J(c,mindx(SX,0),mindx(SX,0)) = par[J3X];
      J(c,mindx(SY,0),mindx(SY,0)) = par[J3Y];
      J(c,mindx(SZ,0),mindx(SZ,0)) = par[J3Z];
    }
#endif  
}

#ifdef PHONONS
void Couplings::Initializeg()
{
#ifdef FAKEHEISENBERG
  for(int c=0; c<NN; c++)
    {
      g(c,mindx(SX,0),mindx(SX,0)) = par[g1X];
    } 

  for(int c=NN; c<NNN; c++)
    {
      g(c,mindx(SX,0),mindx(SX,0)) = par[g2X];
    } 
#else  
  for(int c=0; c<NN; c++)
    {
      g(c,mindx(SX,0),mindx(SX,0)) = par[g1X];
      g(c,mindx(SY,0),mindx(SY,0)) = par[g1Y];
      g(c,mindx(SZ,0),mindx(SZ,0)) = par[g1Z];
    } 

  for(int c=NN; c<NNN; c++)
    {
      g(c,mindx(SX,0),mindx(SX,0)) = par[g2X];
      g(c,mindx(SY,0),mindx(SY,0)) = par[g2Y];
      g(c,mindx(SZ,0),mindx(SZ,0)) = par[g2Z];
    }
#endif  
}
#endif


#elif defined TRIANGULARPHONONS

#ifdef FAKEHEISENBERG
const int NSPIN= 1;
const int NFAKESPINTRACE=3; // a compensating factor for treating Heisenberg models with a single spin comp.
#else
const int NSPIN= 3;
const int NFAKESPINTRACE=1; 
#endif
const int NSUBL= 1;


#ifdef CPOSITIVE
const int NC=9;
const int NN=3;
const int NNN=6;
vector<Triplet> clist{{ 1,  0,  0 }, { 0,  1,  0 }, {-1, 1,  0 }, { 1,  1,  0 }, { -1, 2, 0}, { -2, 1, 0}, { 2, 0, 0}, { 0, 2, 0}, { -2, 2, 0}};
#else
const int NC=18;
const int NN=6;
const int NNN=12;
vector<Triplet> clist{
{ 1,  0,  0 }, {-1,  0, 0 }, { 0, 1, 0 }, { 0, -1, 0 }, { -1, 1, 0 }, { 1, -1, 0 }, 
{ 1,  1,  0 }, {-1, -1, 0 }, {-1, 2, 0 }, { 1, -2, 0 }, { -2, 1, 0 }, { 2, -1, 0 }, 
{ 2,  0,  0 }, {-2,  0, 0 }, { 0, 2, 0 }, { 0, -2, 0 }, { -2, 2, 0 }, { 2, -2, 0 }};
#endif
vector<Coord> roffset={Coord(0,0,0)};
vector<double> invsqrtmasses={1.};

// global routines to extract sublattice and spin indices.

int mindx(const int spin,const int subl=0){ return spin*NSUBL+subl;}
int spin(const int m){return m/NSUBL;}
int subl(const int m){return m%NSUBL;}

// specific coupling rules for squarephonons
void Couplings::InitializeJ()
{
#ifdef FAKEHEISENBERG
  for(int c=0; c<NN; c++)
    {
      J(c,mindx(SX,0),mindx(SX,0)) = par[J1X];
    } 

  for(int c=NN; c<NNN; c++)
    {
      J(c,mindx(SX,0),mindx(SX,0)) = par[J2X];
    } 

  for(int c=NNN; c<NC; c++)
    {
      J(c,mindx(SX,0),mindx(SX,0)) = par[J3X];
    } 
#else
  
  for(int c=0; c<NN; c++)
    {
      J(c,mindx(SX,0),mindx(SX,0)) = par[J1X];
      J(c,mindx(SY,0),mindx(SY,0)) = par[J1Y];
      J(c,mindx(SZ,0),mindx(SZ,0)) = par[J1Z]; 
    } 

  for(int c=NN; c<NNN; c++)
    {
      J(c,mindx(SX,0),mindx(SX,0)) = par[J2X];
      J(c,mindx(SY,0),mindx(SY,0)) = par[J2Y];
      J(c,mindx(SZ,0),mindx(SZ,0)) = par[J2Z];
    } 

  for(int c=NNN; c<NC; c++)
    {
      J(c,mindx(SX,0),mindx(SX,0)) = par[J3X];
      J(c,mindx(SY,0),mindx(SY,0)) = par[J3Y];
      J(c,mindx(SZ,0),mindx(SZ,0)) = par[J3Z];
    }
#endif  
}

#ifdef PHONONS
void Couplings::Initializeg()
{
#ifdef FAKEHEISENBERG
  for(int c=0; c<NN; c++)
    {
      g(c,mindx(SX,0),mindx(SX,0)) = par[g1X];
    } 

  for(int c=NN; c<NNN; c++)
    {
      g(c,mindx(SX,0),mindx(SX,0)) = par[g2X];
    } 
#else  
  for(int c=0; c<NN; c++)
    {
      g(c,mindx(SX,0),mindx(SX,0)) = par[g1X];
      g(c,mindx(SY,0),mindx(SY,0)) = par[g1Y];
      g(c,mindx(SZ,0),mindx(SZ,0)) = par[g1Z];
    } 

  for(int c=NN; c<NNN; c++)
    {
      g(c,mindx(SX,0),mindx(SX,0)) = par[g2X];
      g(c,mindx(SY,0),mindx(SY,0)) = par[g2Y];
      g(c,mindx(SZ,0),mindx(SZ,0)) = par[g2Z];
    }
#endif  
}
#endif


#elif defined CUBICPHONONS

#ifdef FAKEHEISENBERG
const int NSPIN= 1;
const int NFAKESPINTRACE=3; // a compensating factor for treating Heisenberg models with a single spin comp.
#else
const int NSPIN= 3;
const int NFAKESPINTRACE=1; 
#endif
const int NSUBL= 1;


#ifdef CPOSITIVE
/*
const int NC=9;
const int NN=3;
vector<Triplet> clist{{ 1, 0, 0},{ 0, 1, 0},{ 0, 0, 1},{ 1, 1, 0},{ 1,-1, 0},{ 1, 0, 1},{ 1, 0,-1},{ 0, 1, 1},{ 0, 1,-1}};
#else
const int NC=18;
const int NN=6;
vector<Triplet> clist{{ 1, 0, 0},{-1, 0, 0},{ 0, 1, 0},{ 0,-1, 0},{ 0, 0, 1},{ 0, 0,-1},{ 1, 1, 0},{-1,-1, 0},{ 1,-1, 0},{-1, 1, 0},{ 1, 0, 1}, {-1, 0,-1},{ 1, 0,-1}, {-1, 0, 1},{ 0, 1, 1},{ 0,-1,-1},{ 0, 1,-1},{ 0,-1, 1}};
*/
/* nearest neighbor only */
const int NC=3;
const int NN=3;
vector<Triplet> clist{{ 1, 0, 0},{ 0, 1, 0},{ 0, 0, 1}};
#else
const int NC=6;
const int NN=6;
vector<Triplet> clist{{ 1, 0, 0},{-1, 0, 0},{ 0, 1, 0},{ 0,-1, 0},{ 0, 0, 1},{ 0, 0,-1}};

#endif
//vector<int> minusc{1,0,3,2,5,4,7,6,9,8,11,10,13,12,15,14,17,15}; // refers to indx in clist
vector<Coord> roffset={Coord(0,0,0)};
vector<double> invsqrtmasses={1.};

// global routines to extract sublattice and spin indices.

int mindx(const int spin,const int subl=0){ return spin*NSUBL+subl;}
int spin(const int m){return m/NSUBL;}
int subl(const int m){return m%NSUBL;}

// specific coupling rules for squarephonons
void Couplings::InitializeJ()
{
#ifdef FAKEHEISENBERG
  for(int c=0; c<NN; c++)
    {
      J(c,mindx(SX,0),mindx(SX,0)) = par[J1X];
    } 

  for(int c=NN; c<NC; c++)
    {
      J(c,mindx(SX,0),mindx(SX,0)) = par[J2X];
    } 
#else
  
  for(int c=0; c<NN; c++)
    {
      J(c,mindx(SX,0),mindx(SX,0)) = par[J1X];
      J(c,mindx(SY,0),mindx(SY,0)) = par[J1Y];
      J(c,mindx(SZ,0),mindx(SZ,0)) = par[J1Z]; 
    } 

  for(int c=NN; c<NC; c++)
    {
      J(c,mindx(SX,0),mindx(SX,0)) = par[J2X];
      J(c,mindx(SY,0),mindx(SY,0)) = par[J2Y];
      J(c,mindx(SZ,0),mindx(SZ,0)) = par[J2Z];
    } 
#endif  
}

#ifdef PHONONS
void Couplings::Initializeg()
{
#ifdef FAKEHEISENBERG
  for(int c=0; c<NN; c++)
    {
      g(c,mindx(SX,0),mindx(SX,0)) = par[g1X];
    } 

  for(int c=NN; c<NC; c++)
    {
      g(c,mindx(SX,0),mindx(SX,0)) = par[g2X];
    } 
#else  
  for(int c=0; c<NN; c++)
    {
      g(c,mindx(SX,0),mindx(SX,0)) = par[g1X];
      g(c,mindx(SY,0),mindx(SY,0)) = par[g1Y];
      g(c,mindx(SZ,0),mindx(SZ,0)) = par[g1Z];
    } 

  for(int c=NN; c<NC; c++)
    {
      g(c,mindx(SX,0),mindx(SX,0)) = par[g2X];
      g(c,mindx(SY,0),mindx(SY,0)) = par[g2Y];
      g(c,mindx(SZ,0),mindx(SZ,0)) = par[g2Z];
    }
#endif  
}
#endif


#elif defined CrI3

const int NSPIN= 3;
const int NSUBL= 2;

const int NC=7;

enum latticev{O,a1,a2,a3,ma1,ma2,ma3};
vector<Triplet> clist{{0, 0, 0},{ 1,  0,  0 }, {0,  1,  0 }, {-1,  1,  0 }, {-1, 0,  0 }, { 0, -1,  0 }, { 1, -1, 0}};

// global routines to extract sublattice and spin indices.

int mindx(const int spin,const int subl=0){ return spin*NSUBL+subl;}
int spin(const int m){return m/NSUBL;}
int subl(const int m){return m%NSUBL;}

// specific coupling rules for CrI3
void Couplings::InitializeJ()
{
  // 00 sublattices
  J(  O,mindx(SX,0),mindx(SX,0)) = par[AX];
  J(  O,mindx(SY,0),mindx(SY,0)) = par[AY];
  J(  O,mindx(SZ,0),mindx(SZ,0)) = par[AZ];

  J( a1,mindx(SX,0),mindx(SX,0)) = par[J2X];
  J( a1,mindx(SY,0),mindx(SY,0)) = par[J2Y];
  J( a1,mindx(SZ,0),mindx(SZ,0)) = par[J2Z];

  J( a2,mindx(SX,0),mindx(SX,0)) = par[J2X];
  J( a2,mindx(SY,0),mindx(SY,0)) = par[J2Y];
  J( a2,mindx(SZ,0),mindx(SZ,0)) = par[J2Z];

  J( a3,mindx(SX,0),mindx(SX,0)) = par[J2X];
  J( a3,mindx(SY,0),mindx(SY,0)) = par[J2Y];
  J( a3,mindx(SZ,0),mindx(SZ,0)) = par[J2Z];

  J(ma1,mindx(SX,0),mindx(SX,0)) = par[J2X];
  J(ma1,mindx(SY,0),mindx(SY,0)) = par[J2Y];
  J(ma1,mindx(SZ,0),mindx(SZ,0)) = par[J2Z];

  J(ma2,mindx(SX,0),mindx(SX,0)) = par[J2X];
  J(ma2,mindx(SY,0),mindx(SY,0)) = par[J2Y];
  J(ma2,mindx(SZ,0),mindx(SZ,0)) = par[J2Z];

  J(ma3,mindx(SX,0),mindx(SX,0)) = par[J2X];
  J(ma3,mindx(SY,0),mindx(SY,0)) = par[J2Y];
  J(ma3,mindx(SZ,0),mindx(SZ,0)) = par[J2Z];

  //11
  J(  O,mindx(SX,1),mindx(SX,1)) = par[AX];
  J(  O,mindx(SY,1),mindx(SY,1)) = par[AY];
  J(  O,mindx(SZ,1),mindx(SZ,1)) = par[AZ];

  J( a1,mindx(SX,1),mindx(SX,1)) = par[J2X];
  J( a1,mindx(SY,1),mindx(SY,1)) = par[J2Y];
  J( a1,mindx(SZ,1),mindx(SZ,1)) = par[J2Z];

  J( a2,mindx(SX,1),mindx(SX,1)) = par[J2X];
  J( a2,mindx(SY,1),mindx(SY,1)) = par[J2Y];
  J( a2,mindx(SZ,1),mindx(SZ,1)) = par[J2Z];

  J( a3,mindx(SX,1),mindx(SX,1)) = par[J2X];
  J( a3,mindx(SY,1),mindx(SY,1)) = par[J2Y];
  J( a3,mindx(SZ,1),mindx(SZ,1)) = par[J2Z];

  J(ma1,mindx(SX,1),mindx(SX,1)) = par[J2X];
  J(ma1,mindx(SY,1),mindx(SY,1)) = par[J2Y];
  J(ma1,mindx(SZ,1),mindx(SZ,1)) = par[J2Z];

  J(ma2,mindx(SX,1),mindx(SX,1)) = par[J2X];
  J(ma2,mindx(SY,1),mindx(SY,1)) = par[J2Y];
  J(ma2,mindx(SZ,1),mindx(SZ,1)) = par[J2Z];

  J(ma3,mindx(SX,1),mindx(SX,1)) = par[J2X];
  J(ma3,mindx(SY,1),mindx(SY,1)) = par[J2Y];
  J(ma3,mindx(SZ,1),mindx(SZ,1)) = par[J2Z];

  //01

#ifdef KITAEVMODEL
  J(ma2,mindx(SX,0),mindx(SX,1)) = par[J1X];
  J(ma3,mindx(SY,0),mindx(SY,1)) = par[J1Y];
  J(  O,mindx(SZ,0),mindx(SZ,1)) = par[J1Z];
#else
  J(  O,mindx(SX,0),mindx(SX,1)) = par[J1X];
  J(  O,mindx(SY,0),mindx(SY,1)) = par[J1Y];
  J(  O,mindx(SZ,0),mindx(SZ,1)) = par[J1Z];

  J(ma2,mindx(SX,0),mindx(SX,1)) = par[J1X];
  J(ma2,mindx(SY,0),mindx(SY,1)) = par[J1Y];
  J(ma2,mindx(SZ,0),mindx(SZ,1)) = par[J1Z];

  J(ma3,mindx(SX,0),mindx(SX,1)) = par[J1X];
  J(ma3,mindx(SY,0),mindx(SY,1)) = par[J1Y];
  J(ma3,mindx(SZ,0),mindx(SZ,1)) = par[J1Z];
#endif

  // DM-interactions
  J(  O,mindx(SY,0),mindx(SZ,1)) = par[D1X];
  J(  O,mindx(SZ,0),mindx(SY,1)) =-par[D1X];
  J(  O,mindx(SZ,0),mindx(SX,1)) = par[D1Y];
  J(  O,mindx(SX,0),mindx(SZ,1)) =-par[D1Y];
  J(  O,mindx(SX,0),mindx(SY,1)) = par[D1Z];
  J(  O,mindx(SY,0),mindx(SX,1)) =-par[D1Z];

  J(ma2,mindx(SY,0),mindx(SZ,1)) = par[D2X];
  J(ma2,mindx(SZ,0),mindx(SY,1)) =-par[D2X];
  J(ma2,mindx(SZ,0),mindx(SX,1)) = par[D2Y];
  J(ma2,mindx(SX,0),mindx(SZ,1)) =-par[D2Y];
  J(ma2,mindx(SX,0),mindx(SY,1)) = par[D2Z];
  J(ma2,mindx(SY,0),mindx(SX,1)) =-par[D2Z];

  J(ma3,mindx(SY,0),mindx(SZ,1)) = par[D3X];
  J(ma3,mindx(SZ,0),mindx(SY,1)) =-par[D3X];
  J(ma3,mindx(SZ,0),mindx(SX,1)) = par[D3Y];
  J(ma3,mindx(SX,0),mindx(SZ,1)) =-par[D3Y];
  J(ma3,mindx(SX,0),mindx(SY,1)) = par[D3Z];
  J(ma3,mindx(SY,0),mindx(SX,1)) =-par[D3Z];

  //10
#ifdef KITAEVMODEL
  J( a2,mindx(SX,1),mindx(SX,0)) = par[J1X];
  J( a3,mindx(SY,1),mindx(SY,0)) = par[J1Y];
  J(  O,mindx(SZ,1),mindx(SZ,0)) = par[J1Z];
#else
  J(  O,mindx(SX,1),mindx(SX,0)) = par[J1X];
  J(  O,mindx(SY,1),mindx(SY,0)) = par[J1Y];
  J(  O,mindx(SZ,1),mindx(SZ,0)) = par[J1Z];

  J( a2,mindx(SX,1),mindx(SX,0)) = par[J1X];
  J( a2,mindx(SY,1),mindx(SY,0)) = par[J1Y];
  J( a2,mindx(SZ,1),mindx(SZ,0)) = par[J1Z];

  J( a3,mindx(SX,1),mindx(SX,0)) = par[J1X];
  J( a3,mindx(SY,1),mindx(SY,0)) = par[J1Y];
  J( a3,mindx(SZ,1),mindx(SZ,0)) = par[J1Z];
#endif

     // DM-interactions
  J(  O,mindx(SZ,1),mindx(SY,0)) =-par[D1X];
  J(  O,mindx(SY,1),mindx(SZ,0)) = par[D1X];
  J(  O,mindx(SX,1),mindx(SZ,0)) =-par[D1Y];
  J(  O,mindx(SZ,1),mindx(SX,0)) = par[D1Y];
  J(  O,mindx(SY,1),mindx(SX,0)) =-par[D1Z];
  J(  O,mindx(SX,1),mindx(SY,0)) = par[D1Z];

  J( a2,mindx(SZ,1),mindx(SY,0)) =-par[D2X];
  J( a2,mindx(SY,1),mindx(SZ,0)) = par[D2X];
  J( a2,mindx(SX,1),mindx(SZ,0)) =-par[D2Y];
  J( a2,mindx(SZ,1),mindx(SX,0)) = par[D2Y];
  J( a2,mindx(SY,1),mindx(SX,0)) =-par[D2Z];
  J( a2,mindx(SX,1),mindx(SY,0)) = par[D2Z];

  J( a3,mindx(SZ,1),mindx(SY,0)) =-par[D3X];
  J( a3,mindx(SY,1),mindx(SZ,0)) = par[D3X];
  J( a3,mindx(SX,1),mindx(SZ,0)) =-par[D3Y];
  J( a3,mindx(SZ,1),mindx(SX,0)) = par[D3Y];
  J( a3,mindx(SY,1),mindx(SX,0)) =-par[D3Z];
  J( a3,mindx(SX,1),mindx(SY,0)) = par[D3Z];
}

#endif


const int NMAT = NSPIN*NSUBL;
const int NMAT2= NMAT*NMAT;

#endif // MODELDEF_H
