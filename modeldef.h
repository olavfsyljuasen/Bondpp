#ifndef MODELDEF_H
#define MODELDEF_H

enum spindirections{SX,SY,SZ};

using namespace std;




#ifdef SQUAREPHONONS


const int NSPIN= 3;
const int NSUBL= 1;
const int NFAKETRACE=1; // a compensating factor for treating Heisenberg models with a single spin comp.


const int NC=12;

vector<Triplet> clist{{ 1,  0,  0 }, {-1,  0,  0 }, { 0,  1,  0 }, { 0, -1,  0 }, { 1,  1,  0 }, {-1, -1, 0}, { 1, -1,  0 }, {-1,  1,  0 },{ 2,0,0},{-2,0,0},{0,2,0},{0,-2,0}};
vector<int> minusc{1,0,3,2,5,4,7,6,9,8,11,10}; // refers to indx in clist
vector<Coord> roffset={Coord(0,0,0)};
vector<double> invsqrtmasses={1.};

// global routines to extract sublattice and spin indices.

int mindx(const int spin,const int subl=0){ return spin*NSUBL+subl;}
int spin(const int m){return m/NSUBL;}
int subl(const int m){return m%NSUBL;}

// specific coupling rules for squarephonons
void Couplings::InitializeJ()
{
  for(int c=0; c<4; c++)
    {
      J(c,mindx(SX,0),mindx(SX,0)) = par[J1X];
      J(c,mindx(SY,0),mindx(SY,0)) = par[J1Y];
      J(c,mindx(SZ,0),mindx(SZ,0)) = par[J1Z];
    } 

  for(int c=4; c<8; c++)
    {
      J(c,mindx(SX,0),mindx(SX,0)) = par[J2X];
      J(c,mindx(SY,0),mindx(SY,0)) = par[J2Y];
      J(c,mindx(SZ,0),mindx(SZ,0)) = par[J2Z];
    } 

  for(int c=8; c<12; c++)
    {
      J(c,mindx(SX,0),mindx(SX,0)) = par[J3X];
      J(c,mindx(SY,0),mindx(SY,0)) = par[J3Y];
      J(c,mindx(SZ,0),mindx(SZ,0)) = par[J3Z];
    } 
}

#ifdef PHONONS
void Couplings::Initializeg()
{
  for(int c=0; c<4; c++)
    {
      g(c,mindx(SX,0),mindx(SX,0)) = par[g1X];
      g(c,mindx(SY,0),mindx(SY,0)) = par[g1Y];
      g(c,mindx(SZ,0),mindx(SZ,0)) = par[g1Z];
    } 

  for(int c=4; c<8; c++)
    {
      g(c,mindx(SX,0),mindx(SX,0)) = par[g2X];
      g(c,mindx(SY,0),mindx(SY,0)) = par[g2Y];
      g(c,mindx(SZ,0),mindx(SZ,0)) = par[g2Z];
    } 
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
