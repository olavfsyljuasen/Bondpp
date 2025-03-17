#ifndef MODELCOUPLINGS_H
#define MODELCOUPLINGS_H

enum spindirections{SX,SY,SZ};

using namespace std;

vector<int> nonzeroclist(0); // a list of c indices which have gc !=0. So it is only necessary to loop through those values.


#if defined SQUAREPHONONS

#if NBRRANGE==1

#ifdef CPOSITIVE
const auto NC=2;
const auto NN=2;
const auto NNN=2;
vector<Triplet> clist{{ 1,  0,  0 }, { 0,  1,  0 }};
#else
const auto NC=4;
const auto NN=4;
const auto NNN=4;
vector<Triplet> clist{{ 1,  0,  0 }, {-1,  0,  0 }, { 0,  1,  0 }, { 0, -1,  0 }};
#endif //CPOSITIVE


#elif NBRRANGE==2

#ifdef CPOSITIVE
const auto NC=4;
const auto NN=2;
const auto NNN=4;
vector<Triplet> clist{{ 1,  0,  0 }, { 0,  1,  0 }, { 1,  1,  0 }, { 1, -1,  0 }};
#else
const auto NC=8;
const auto NN=4;
const auto NNN=8;
vector<Triplet> clist{{ 1,  0,  0 }, {-1,  0,  0 }, { 0,  1,  0 }, { 0, -1,  0 }, { 1,  1,  0 }, {-1, -1, 0}, { 1, -1,  0 }, {-1,  1,  0 }};
#endif //CPOSITIVE

#elif NBRRANGE==3

#ifdef CPOSITIVE
const auto NC=6;
const auto NN=2;
const auto NNN=4;
vector<Triplet> clist{{ 1,  0,  0 }, { 0,  1,  0 }, { 1,  1,  0 }, { 1, -1,  0 }, {2, 0, 0}, { 0, 2, 0}};
#else
const auto NC=12;
const auto NN=4;
const auto NNN=8;
vector<Triplet> clist{{ 1,  0,  0 }, {-1,  0,  0 }, { 0,  1,  0 }, { 0, -1,  0 }, { 1,  1,  0 }, {-1, -1, 0}, { 1, -1,  0 }, {-1,  1,  0 }, { 2,  0,  0 }, {-2,  0,  0 }, { 0,  2,  0 }, { 0, -2,  0 }};
#endif //CPOSITIVE

#else

#ifdef CPOSITIVE
const auto NC=6;
const auto NN=2;
const auto NNN=4;
vector<Triplet> clist{{ 1,  0,  0 }, { 0,  1,  0 }, { 1,  1,  0 }, { 1, -1,  0 }, { 2, 0, 0}, { 0, 2, 0}};
#else
const auto NC=12;
const auto NN=4;
const auto NNN=8;
vector<Triplet> clist{{ 1,  0,  0 }, {-1,  0,  0 }, { 0,  1,  0 }, { 0, -1,  0 }, { 1,  1,  0 }, {-1, -1, 0}, { 1, -1,  0 }, {-1,  1,  0 },{ 2,0,0},{-2,0,0},{0,2,0},{0,-2,0}};
#endif //CPOSITIVE

#endif  // NBRRANGE





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
#endif // FAKEHEISENBERG 
}

#ifdef PHONONS
void Couplings::Initializeg()
{
#ifdef FAKEHEISENBERG
  for(int c=0; c<NN; c++)
    {
     //      g(c,mindx(SX,0),mindx(SX,0)) = par[g1X];  // use g =1/r dJ/dr
      g(c,mindx(SX,0),mindx(SX,0)) = par[g1X]/norm(clist[c]); // use g = dJ/dr
      if( g(c,mindx(SX,0),mindx(SX,0)) != 0.){ nonzeroclist.push_back(c);} // make a list of non-zero c's to loop around
    } 

  for(int c=NN; c<NNN; c++)
    {
      g(c,mindx(SX,0),mindx(SX,0)) = par[g2X]/norm(clist[c]);
      if( g(c,mindx(SX,0),mindx(SX,0)) != 0.){ nonzeroclist.push_back(c);} 
    } 
#else  
  for(int c=0; c<NN; c++)
    {
      g(c,mindx(SX,0),mindx(SX,0)) = par[g1X]/norm(clist[c]);
      g(c,mindx(SY,0),mindx(SY,0)) = par[g1Y]/norm(clist[c]);
      g(c,mindx(SZ,0),mindx(SZ,0)) = par[g1Z]/norm(clist[c]);
      if( g(c,mindx(SX,0),mindx(SX,0)) != 0. || g(c,mindx(SY,0),mindx(SY,0) != 0. || g(c,mindx(SZ,0),mindx(SZ,0) != 0. ){ nonzeroclist.push_back(c);} 
    } 

  for(int c=NN; c<NNN; c++)
    {
      g(c,mindx(SX,0),mindx(SX,0)) = par[g2X]/norm(clist[c]);
      g(c,mindx(SY,0),mindx(SY,0)) = par[g2Y]/norm(clist[c]);
      g(c,mindx(SZ,0),mindx(SZ,0)) = par[g2Z]/norm(clist[c]);
      if( g(c,mindx(SX,0),mindx(SX,0)) != 0. || g(c,mindx(SY,0),mindx(SY,0) != 0. || g(c,mindx(SZ,0),mindx(SZ,0) != 0. ){ nonzeroclist.push_back(c);} 
    }
#endif  // FAKEHEISENBERG
}
#endif //PHONONS


#elif defined TRIANGULARPHONONS

#if NBRRANGE==1

#ifdef CPOSITIVE
const auto NC=3;
const auto NN=3;
const auto NNN=3;
vector<Triplet> clist{{ 1,  0,  0 }, { 0,  1,  0 }, {-1, 1,  0 }};
#else
const auto NC=6;
const auto NN=6;
const auto NNN=6;
 vector<Triplet> clist{{ 1,  0,  0 }, {-1, 0, 0},{ 0,  1,  0 }, {0, -1, 0}, {-1, 1,  0 },{ 1,-1, 0}};
#endif //CPOSITIVE

#elif NBRRANGE==2

#ifdef CPOSITIVE
const auto NC=6;
const auto NN=3;
const auto NNN=6;
vector<Triplet> clist{{ 1,  0,  0 }, { 0,  1,  0 }, {-1, 1,  0 }, { 1,  1,  0 }, { -1, 2, 0}, { -2, 1, 0}, { 2, 0, 0}, { 0, 2, 0}, { -2, 2, 0}};
#else
const auto NC=12;
const auto NN=6;
const auto NNN=12;
vector<Triplet> clist{
{ 1,  0,  0 }, {-1,  0, 0 }, { 0, 1, 0 }, { 0, -1, 0 }, { -1, 1, 0 }, { 1, -1, 0 }, 
{ 1,  1,  0 }, {-1, -1, 0 }, {-1, 2, 0 }, { 1, -2, 0 }, { -2, 1, 0 }, { 2, -1, 0 }};
#endif //CPOSITIVE

#elif NBRRANGE==3

#ifdef CPOSITIVE
const auto NC=9;
const auto NN=3;
const auto NNN=6;
 vector<Triplet> clist{{ 1,  0,  0 }, { 0,  1,  0 }, {-1, 1,  0 }, { 1,  1,  0 }, { -1, 2, 0}, { -2, 1, 0}, { 2, 0, 0}, { 0, 2, 0}, { -2, 2, 0}};
#else
const auto NC=18;
const auto NN=6;
const auto NNN=12;
vector<Triplet> clist{
{ 1,  0,  0 }, {-1,  0, 0 }, { 0, 1, 0 }, { 0, -1, 0 }, { -1, 1, 0 }, { 1, -1, 0 }, 
{ 1,  1,  0 }, {-1, -1, 0 }, {-1, 2, 0 }, { 1, -2, 0 }, { -2, 1, 0 }, { 2, -1, 0 }, 
{ 2,  0,  0 }, {-2,  0, 0 }, { 0, 2, 0 }, { 0, -2, 0 }, { -2, 2, 0 }, { 2, -2, 0 }};
#endif // CPOSITIVE

#endif // NBRRANGE
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
#endif //FAKEHEISENBERG
}

#ifdef PHONONS
void Couplings::Initializeg()
{
#ifdef FAKEHEISENBERG
  for(int c=0; c<NN; c++)
    {
      //      g(c,mindx(SX,0),mindx(SX,0)) = par[g1X];
      g(c,mindx(SX,0),mindx(SX,0)) = par[g1X]/norm(clist[c]); // use g = dJ/dr
      if( g(c,mindx(SX,0),mindx(SX,0)) != 0.){ nonzeroclist.push_back(c);} // make a list of non-zero c's to loop around
      // if(g(c,mindx(SX,0),mindx(SX,0) != 0.){ nonzeroclist.pushback(c);} // make a list of non-zero c's to loop around
    } 

  for(int c=NN; c<NNN; c++)
    {
      g(c,mindx(SX,0),mindx(SX,0)) = par[g2X]/norm(clist[c]); // use g = dJ/dr;
      if( g(c,mindx(SX,0),mindx(SX,0)) != 0.){ nonzeroclist.push_back(c);} // make a list of non-zero c's to loop around
    }

  for(int c=NNN; c<NC; c++)
    {
      g(c,mindx(SX,0),mindx(SX,0)) = par[g3X]/norm(clist[c]); // use g = dJ/dr;
      if( g(c,mindx(SX,0),mindx(SX,0)) != 0.){ nonzeroclist.push_back(c);} // make a list of non-zero c's to loop around
    } 
#else  
  for(int c=0; c<NN; c++)
    {
      g(c,mindx(SX,0),mindx(SX,0)) = par[g1X]/norm(clist[c]); // use g = dJ/dr;
      g(c,mindx(SY,0),mindx(SY,0)) = par[g1Y]/norm(clist[c]); // use g = dJ/dr;
      g(c,mindx(SZ,0),mindx(SZ,0)) = par[g1Z]/norm(clist[c]); // use g = dJ/dr;
      if( g(c,mindx(SX,0),mindx(SX,0)) != 0. || g(c,mindx(SY,0),mindx(SY,0) != 0. || g(c,mindx(SZ,0),mindx(SZ,0) != 0. ){ nonzeroclist.push_back(c);} 
    } 

  for(int c=NN; c<NNN; c++)
    {
      g(c,mindx(SX,0),mindx(SX,0)) = par[g2X]/norm(clist[c]); // use g = dJ/dr;
      g(c,mindx(SY,0),mindx(SY,0)) = par[g2Y]/norm(clist[c]); // use g = dJ/dr;
      g(c,mindx(SZ,0),mindx(SZ,0)) = par[g2Z]/norm(clist[c]); // use g = dJ/dr;
      if( g(c,mindx(SX,0),mindx(SX,0)) != 0. || g(c,mindx(SY,0),mindx(SY,0) != 0. || g(c,mindx(SZ,0),mindx(SZ,0) != 0. ){ nonzeroclist.push_back(c);} 
    }

    for(int c=NNN; c<NC; c++)
    {
      g(c,mindx(SX,0),mindx(SX,0)) = par[g3X]/norm(clist[c]); // use g = dJ/dr;
      g(c,mindx(SY,0),mindx(SY,0)) = par[g3Y]/norm(clist[c]); // use g = dJ/dr;
      g(c,mindx(SZ,0),mindx(SZ,0)) = par[g3Z]/norm(clist[c]); // use g = dJ/dr;
      if( g(c,mindx(SX,0),mindx(SX,0)) != 0. || g(c,mindx(SY,0),mindx(SY,0) != 0. || g(c,mindx(SZ,0),mindx(SZ,0) != 0. ){ nonzeroclist.push_back(c);} 
    }
#endif //FAKEHEISENBERG 
}
#endif //PHONONS


#elif defined CUBICPHONONS

#ifdef CPOSITIVE
const auto NC=9;
const auto NN=3;
vector<Triplet> clist{{ 1, 0, 0}, { 0, 1, 0}, { 0, 0, 1}, { 0, 1, 1}, { 0, 1,-1}, { 1, 0, 1}, { 1, 0, -1}, { 1, 1, 0}, { 1,-1, 0}};
#else
const auto NC=18;
const auto NN=6;
vector<Triplet> clist{{ 1, 0, 0}, { 0, 1, 0}, { 0, 0, 1}, { 0, 1, 1}, { 0, 1,-1}, { 1, 0, 1}, { 1, 0, -1}, { 1, 1, 0}, { 1,-1, 0},{-1, 0, 0}, { 0,-1, 0}, { 0, 0,-1}, { 0,-1,-1}, { 0,-1, 1}, {-1, 0,-1}, {-1, 0,  1}, {-1,-1, 0}, {-1, 1, 0}}
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
#endif //FAKEHEISENBERG
}

#ifdef PHONONS
void Couplings::Initializeg()
{
#ifdef FAKEHEISENBERG
  for(int c=0; c<NN; c++)
    {
      g(c,mindx(SX,0),mindx(SX,0)) = par[g1X]/norm(clist[c]); // use g = dJ/dr;
      if( g(c,mindx(SX,0),mindx(SX,0)) != 0.){ nonzeroclist.push_back(c);} // make a list of non-zero c's to loop around
    } 

  for(int c=NN; c<NC; c++)
    {
      g(c,mindx(SX,0),mindx(SX,0)) = par[g2X]/norm(clist[c]); // use g = dJ/dr;
      if( g(c,mindx(SX,0),mindx(SX,0)) != 0.){ nonzeroclist.push_back(c);} // make a list of non-zero c's to loop around
    } 
#else  
  for(int c=0; c<NN; c++)
    {
      g(c,mindx(SX,0),mindx(SX,0)) = par[g1X]/norm(clist[c]); // use g = dJ/dr;
      g(c,mindx(SY,0),mindx(SY,0)) = par[g1Y]/norm(clist[c]); // use g = dJ/dr;
      g(c,mindx(SZ,0),mindx(SZ,0)) = par[g1Z]/norm(clist[c]); // use g = dJ/dr;
      if( g(c,mindx(SX,0),mindx(SX,0)) != 0. || g(c,mindx(SY,0),mindx(SY,0) != 0. || g(c,mindx(SZ,0),mindx(SZ,0) != 0. ){ nonzeroclist.push_back(c);} 
    } 

  for(int c=NN; c<NC; c++)
    {
      g(c,mindx(SX,0),mindx(SX,0)) = par[g2X]/norm(clist[c]); // use g = dJ/dr;
      g(c,mindx(SY,0),mindx(SY,0)) = par[g2Y]/norm(clist[c]); // use g = dJ/dr;
      g(c,mindx(SZ,0),mindx(SZ,0)) = par[g2Z]/norm(clist[c]); // use g = dJ/dr;
      if( g(c,mindx(SX,0),mindx(SX,0)) != 0. || g(c,mindx(SY,0),mindx(SY,0) != 0. || g(c,mindx(SZ,0),mindx(SZ,0) != 0. ){ nonzeroclist.push_back(c);} 
    }
#endif //FAKEHEISENBERG
}
#endif //PHONONS

#elif defined FCCPHONONS

#ifdef CPOSITIVE
const auto NC=9;
const auto NN=3;
vector<Triplet> clist{{ 1, 0, 0}, { 0, 1, 0}, { 0, 0, 1}, { 0, 1, 1}, { 0, 1,-1}, { 1, 0, 1}, { 1, 0, -1}, { 1, 1, 0}, { 1,-1, 0}};
#else
const auto NC=18;
const auto NN=6;
vector<Triplet> clist{{ 1, 0, 0}, { 0, 1, 0}, { 0, 0, 1}, { 0, 1, 1}, { 0, 1,-1}, { 1, 0, 1}, { 1, 0, -1}, { 1, 1, 0}, { 1,-1, 0},{-1, 0, 0}, { 0,-1, 0}, { 0, 0,-1}, { 0,-1,-1}, { 0,-1, 1}, {-1, 0,-1}, {-1, 0,  1}, {-1,-1, 0}, {-1, 1, 0}}
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
#endif //FAKEHEISENBERG 
}

#ifdef PHONONS
void Couplings::Initializeg()
{
#ifdef FAKEHEISENBERG
  for(int c=0; c<NN; c++)
    {
      g(c,mindx(SX,0),mindx(SX,0)) = par[g1X]/norm(clist[c]); // use g = dJ/dr;
      if( g(c,mindx(SX,0),mindx(SX,0)) != 0.){ nonzeroclist.push_back(c);} // make a list of non-zero c's to loop around
    } 

  for(int c=NN; c<NC; c++)
    {
      g(c,mindx(SX,0),mindx(SX,0)) = par[g2X]/norm(clist[c]); // use g = dJ/dr;
      if( g(c,mindx(SX,0),mindx(SX,0)) != 0.){ nonzeroclist.push_back(c);} // make a list of non-zero c's to loop around
    } 
#else  
  for(int c=0; c<NN; c++)
    {
      g(c,mindx(SX,0),mindx(SX,0)) = par[g1X]/norm(clist[c]); // use g = dJ/dr;
      g(c,mindx(SY,0),mindx(SY,0)) = par[g1Y]/norm(clist[c]); // use g = dJ/dr;
      g(c,mindx(SZ,0),mindx(SZ,0)) = par[g1Z]/norm(clist[c]); // use g = dJ/dr;
      if( g(c,mindx(SX,0),mindx(SX,0)) != 0. || g(c,mindx(SY,0),mindx(SY,0) != 0. || g(c,mindx(SZ,0),mindx(SZ,0) != 0. ){ nonzeroclist.push_back(c);} 
    } 

  for(int c=NN; c<NC; c++)
    {
      g(c,mindx(SX,0),mindx(SX,0)) = par[g2X]/norm(clist[c]); // use g = dJ/dr;
      g(c,mindx(SY,0),mindx(SY,0)) = par[g2Y]/norm(clist[c]); // use g = dJ/dr;
      g(c,mindx(SZ,0),mindx(SZ,0)) = par[g2Z]/norm(clist[c]); // use g = dJ/dr;
      if( g(c,mindx(SX,0),mindx(SX,0)) != 0. || g(c,mindx(SY,0),mindx(SY,0) != 0. || g(c,mindx(SZ,0),mindx(SZ,0) != 0. ){ nonzeroclist.push_back(c);} 
    }
#endif  //FAKEHEISENBERG
}
#endif //PHONONS

#elif defined BCCPHONONS

#ifdef CPOSITIVE
/* g's for nearest neighbors only, note that it is still possible to couple springs to next nearest neighbors*/
const auto NC=4;
const auto NN=4;
vector<Triplet> clist{{ 1, 0, 0}, { 0, 1, 0}, { 0, 0, 1}, { 1, 1, 1}};
#else
const auto NC=8;
const auto NN=8;
vector<Triplet> clist{{ 1, 0, 0}, { 0, 1, 0}, { 0, 0, 1}, { 1, 1, 1},{-1, 0, 0}, { 0,-1, 0}, { 0, 0,-1}, {-1,-1,-1}};
#endif // CPOSITIVE

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
#endif  //FAKEHEISENBERG
}

#ifdef PHONONS
void Couplings::Initializeg()
{
#ifdef FAKEHEISENBERG
  for(int c=0; c<NN; c++)
    {
      g(c,mindx(SX,0),mindx(SX,0)) = par[g1X]/norm(clist[c]); // use g = dJ/dr;
      if( g(c,mindx(SX,0),mindx(SX,0)) != 0.){ nonzeroclist.push_back(c);} // make a list of non-zero c's to loop around
    } 

  for(int c=NN; c<NC; c++)
    {
      g(c,mindx(SX,0),mindx(SX,0)) = par[g2X]/norm(clist[c]); // use g = dJ/dr;
      if( g(c,mindx(SX,0),mindx(SX,0)) != 0.){ nonzeroclist.push_back(c);} // make a list of non-zero c's to loop around
    } 
#else  
  for(int c=0; c<NN; c++)
    {
      g(c,mindx(SX,0),mindx(SX,0)) = par[g1X]/norm(clist[c]); // use g = dJ/dr;
      g(c,mindx(SY,0),mindx(SY,0)) = par[g1Y]/norm(clist[c]); // use g = dJ/dr;
      g(c,mindx(SZ,0),mindx(SZ,0)) = par[g1Z]/norm(clist[c]); // use g = dJ/dr;
      if( g(c,mindx(SX,0),mindx(SX,0)) != 0. || g(c,mindx(SY,0),mindx(SY,0) != 0. || g(c,mindx(SZ,0),mindx(SZ,0) != 0. ){ nonzeroclist.push_back(c);} 
    } 

  for(int c=NN; c<NC; c++)
    {
      g(c,mindx(SX,0),mindx(SX,0)) = par[g2X]/norm(clist[c]); // use g = dJ/dr;
      g(c,mindx(SY,0),mindx(SY,0)) = par[g2Y]/norm(clist[c]); // use g = dJ/dr;
      g(c,mindx(SZ,0),mindx(SZ,0)) = par[g2Z]/norm(clist[c]); // use g = dJ/dr;
      if( g(c,mindx(SX,0),mindx(SX,0)) != 0. || g(c,mindx(SY,0),mindx(SY,0) != 0. || g(c,mindx(SZ,0),mindx(SZ,0) != 0. ){ nonzeroclist.push_back(c);} 
    }
#endif  //FAKEHEISENBERG
}
#endif //PHONONS


#elif defined CrI3

enum latticev{O,a1,a2,a3,ma1,ma2,ma3};
const auto NC=7;
vector<Triplet> clist{{0, 0, 0},{ 1,  0,  0 }, {0,  1,  0 }, {-1,  1,  0 }, {-1, 0,  0 }, { 0, -1,  0 }, { 1, -1, 0}};


vector<Coord> roffset={Coord(0,0,0),Coord(0,1./SQRTTHREEOVERTWO,0)};
vector<double> invsqrtmasses={1.,1.};

// global routines to extract sublattice and spin indices.

int mindx(const int spin,const int subl=0){ return spin*NSUBL+subl;}
int spin(const int m){return m/NSUBL;}
int subl(const int m){return m%NSUBL;}

// specific coupling rules for CrI3
#ifdef SPINISOTROPIC
void Couplings::InitializeJ()
{
  double j1x= (par[USEPHI] >0 ? par[A]*cos(par[PHIOVERPI]*PI): par[J1X]);


  // 00 sublattices
  J( a1,mindx(SX,0),mindx(SX,0)) = par[J2X];
  J( a2,mindx(SX,0),mindx(SX,0)) = par[J2X];
  J( a3,mindx(SX,0),mindx(SX,0)) = par[J2X];
  J(ma1,mindx(SX,0),mindx(SX,0)) = par[J2X];
  J(ma2,mindx(SX,0),mindx(SX,0)) = par[J2X];
  J(ma3,mindx(SX,0),mindx(SX,0)) = par[J2X];

  //11
  J( a1,mindx(SX,1),mindx(SX,1)) = par[J2X];
  J( a2,mindx(SX,1),mindx(SX,1)) = par[J2X];
  J( a3,mindx(SX,1),mindx(SX,1)) = par[J2X];
  J(ma1,mindx(SX,1),mindx(SX,1)) = par[J2X];
  J(ma2,mindx(SX,1),mindx(SX,1)) = par[J2X];
  J(ma3,mindx(SX,1),mindx(SX,1)) = par[J2X];


  //01
  J(  O,mindx(SX,0),mindx(SX,1)) = j1x;
  J(ma2,mindx(SX,0),mindx(SX,1)) = j1x;
  J(ma3,mindx(SX,0),mindx(SX,1)) = j1x;

  //10
  J(  O,mindx(SX,1),mindx(SX,0)) = j1x;
  J( a2,mindx(SX,1),mindx(SX,0)) = j1x;
  J( a3,mindx(SX,1),mindx(SX,0)) = j1x;
}
#else
void Couplings::InitializeJ()
{
  double j1x= (par[USEPHI] >0 ? par[A]*cos(par[PHIOVERPI]*PI): par[J1X]);
  double j1y= (par[USEPHI] >0 ? par[A]*cos(par[PHIOVERPI]*PI): par[J1Y]);
  double j1z= (par[USEPHI] >0 ? par[A]*cos(par[PHIOVERPI]*PI): par[J1Z]);

  double k1x= (par[USEPHI] >0 ? 2*par[A]*sin(par[PHIOVERPI]*PI): par[K1X]);
  double k1y= (par[USEPHI] >0 ? 2*par[A]*sin(par[PHIOVERPI]*PI): par[K1Y]);
  double k1z= (par[USEPHI] >0 ? 2*par[A]*sin(par[PHIOVERPI]*PI): par[K1Z]);

  //DM interaction

  // specify DM vector on vertical bonds from 0 to 1
  double d1x_O  = (USEDMRPZ > 0 ? par[DM1BP] : par[DM1X]);
  double d1y_O  = (USEDMRPZ > 0 ? par[DM1BR] : par[DM1Y]);
  double d1z_O  = (USEDMRPZ > 0 ? par[DM1BZ] : par[DM1Z]);

  // rotate DM-vector with bond orientation if RPZ system, else uniform DM-vector
  double angle=( USEDMRPZ>0 ? TWOPI/3.: 0.);

  double d1x_ma2=  cos(angle)*d1x_O-sin(angle)*d1y_O;
  double d1y_ma2=  sin(angle)*d1x_O+cos(angle)*d1y_O;
  double d1z_ma2=  d1z_O;

  double d1x_ma3=  cos(2*angle)*d1x_O-sin(2*angle)*d1y_O;
  double d1y_ma3=  sin(2*angle)*d1x_O+cos(2*angle)*d1y_O;
  double d1z_ma3=  d1z_O;


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
  J(  O,mindx(SX,0),mindx(SX,1)) = j1x;
  J(  O,mindx(SY,0),mindx(SY,1)) = j1y;
  J(  O,mindx(SZ,0),mindx(SZ,1)) = j1z;

  J(ma2,mindx(SX,0),mindx(SX,1)) = j1x;
  J(ma2,mindx(SY,0),mindx(SY,1)) = j1y;
  J(ma2,mindx(SZ,0),mindx(SZ,1)) = j1z;

  J(ma3,mindx(SX,0),mindx(SX,1)) = j1x;
  J(ma3,mindx(SY,0),mindx(SY,1)) = j1y;
  J(ma3,mindx(SZ,0),mindx(SZ,1)) = j1z;

  J(ma2,mindx(SX,0),mindx(SX,1)) += k1x;
  J(  O,mindx(SZ,0),mindx(SZ,1)) += k1y;
  J(ma3,mindx(SY,0),mindx(SY,1)) += k1z;

  //10
  J(  O,mindx(SX,1),mindx(SX,0)) = j1x;
  J(  O,mindx(SY,1),mindx(SY,0)) = j1y;
  J(  O,mindx(SZ,1),mindx(SZ,0)) = j1z;

  J( a2,mindx(SX,1),mindx(SX,0)) = j1x;
  J( a2,mindx(SY,1),mindx(SY,0)) = j1y;
  J( a2,mindx(SZ,1),mindx(SZ,0)) = j1z;

  J( a3,mindx(SX,1),mindx(SX,0)) = j1x;
  J( a3,mindx(SY,1),mindx(SY,0)) = j1y;
  J( a3,mindx(SZ,1),mindx(SZ,0)) = j1z;

  J( a2,mindx(SX,1),mindx(SX,0)) += k1x;
  J(  O,mindx(SZ,1),mindx(SZ,0)) += k1y;
  J( a3,mindx(SY,1),mindx(SY,0)) += k1z;


  // DM-interactions
  J(  O,mindx(SY,0),mindx(SZ,1)) = d1x_O;
  J(  O,mindx(SZ,0),mindx(SY,1)) =-d1x_O;
  J(  O,mindx(SZ,0),mindx(SX,1)) = d1y_O;
  J(  O,mindx(SX,0),mindx(SZ,1)) =-d1y_O;
  J(  O,mindx(SX,0),mindx(SY,1)) = d1z_O;
  J(  O,mindx(SY,0),mindx(SX,1)) =-d1z_O;

  J(  O,mindx(SY,1),mindx(SZ,0)) =-d1x_O;
  J(  O,mindx(SZ,1),mindx(SY,0)) = d1x_O;
  J(  O,mindx(SZ,1),mindx(SX,0)) =-d1y_O;
  J(  O,mindx(SX,1),mindx(SZ,0)) = d1y_O;
  J(  O,mindx(SX,1),mindx(SY,0)) =-d1z_O;
  J(  O,mindx(SY,1),mindx(SX,0)) = d1z_O;

  J(ma2,mindx(SY,0),mindx(SZ,1)) = d1x_ma2;
  J(ma2,mindx(SZ,0),mindx(SY,1)) =-d1x_ma2;
  J(ma2,mindx(SZ,0),mindx(SX,1)) = d1y_ma2;
  J(ma2,mindx(SX,0),mindx(SZ,1)) =-d1y_ma2;
  J(ma2,mindx(SX,0),mindx(SY,1)) = d1z_ma2;
  J(ma2,mindx(SY,0),mindx(SX,1)) =-d1z_ma2;

  J( a2,mindx(SY,1),mindx(SZ,0)) =-d1x_ma2;
  J( a2,mindx(SZ,1),mindx(SY,0)) = d1x_ma2;
  J( a2,mindx(SZ,1),mindx(SX,0)) =-d1y_ma2;
  J( a2,mindx(SX,1),mindx(SZ,0)) = d1y_ma2;
  J( a2,mindx(SX,1),mindx(SY,0)) =-d1z_ma2;
  J( a2,mindx(SY,1),mindx(SX,0)) = d1z_ma2;

  J(ma3,mindx(SY,0),mindx(SZ,1)) = d1x_ma3;
  J(ma3,mindx(SZ,0),mindx(SY,1)) =-d1x_ma3;
  J(ma3,mindx(SZ,0),mindx(SX,1)) = d1y_ma3;
  J(ma3,mindx(SX,0),mindx(SZ,1)) =-d1y_ma3;
  J(ma3,mindx(SX,0),mindx(SY,1)) = d1z_ma3;
  J(ma3,mindx(SY,0),mindx(SX,1)) =-d1z_ma3;

  J( a3,mindx(SY,0),mindx(SZ,1)) =-d1x_ma3;
  J( a3,mindx(SZ,0),mindx(SY,1)) = d1x_ma3;
  J( a3,mindx(SZ,0),mindx(SX,1)) =-d1y_ma3;
  J( a3,mindx(SX,0),mindx(SZ,1)) = d1y_ma3;
  J( a3,mindx(SX,0),mindx(SY,1)) =-d1z_ma3;
  J( a3,mindx(SY,0),mindx(SX,1)) = d1z_ma3;
}
#endif //SPINISOTROPIC

#elif defined SQUARENOPHONONS
enum latticev{O,a1,a2,ma1,ma2,a1a2,ma1a2,ma1ma2,a1ma2};
const auto NC=9;
vector<Triplet> clist{{ 0, 0, 0}, { 1,  0,  0 }, {0,  1,  0 }, {-1, 0,  0 }, { 0, -1,  0 }, { 1, 1, 0}, {-1, 1, 0}, {-1,-1, 0},{ 1,-1, 0} };


vector<Coord> roffset={Coord(0,0,0)};
vector<double> invsqrtmasses={1.};

// global routines to extract sublattice and spin indices.

int mindx(const int spin,const int subl=0){ return spin*NSUBL+subl;}
int spin(const int m){return m/NSUBL;}
int subl(const int m){return m%NSUBL;}


#ifdef SPINISOTROPIC
void Couplings::InitializeJ()
{
  J( a1,mindx(SX,0),mindx(SX,0)) = par[J1X];
  J( a2,mindx(SX,0),mindx(SX,0)) = par[J1X];
  J(ma1,mindx(SX,0),mindx(SX,0)) = par[J1X];
  J(ma2,mindx(SX,0),mindx(SX,0)) = par[J1X];

  J(  a1a2,mindx(SX,0),mindx(SX,0)) = par[J2X];
  J( ma1a2,mindx(SX,0),mindx(SX,0)) = par[J2X];
  J(ma1ma2,mindx(SX,0),mindx(SX,0)) = par[J2X];
  J( a1ma2,mindx(SX,0),mindx(SX,0)) = par[J2X];
}
#else
void Couplings::InitializeJ()
{
  // 00 sublattices
  J(  O,mindx(SX,0),mindx(SX,0)) = par[AX];
  J(  O,mindx(SY,0),mindx(SY,0)) = par[AY];
  J(  O,mindx(SZ,0),mindx(SZ,0)) = par[AZ];

  J( a1,mindx(SX,0),mindx(SX,0)) = par[J1X];
  J( a1,mindx(SY,0),mindx(SY,0)) = par[J1Y];
  J( a1,mindx(SZ,0),mindx(SZ,0)) = par[J1Z];

  J( a2,mindx(SX,0),mindx(SX,0)) = par[J1X];
  J( a2,mindx(SY,0),mindx(SY,0)) = par[J1Y];
  J( a2,mindx(SZ,0),mindx(SZ,0)) = par[J1Z];

  J(ma1,mindx(SX,0),mindx(SX,0)) = par[J1X];
  J(ma1,mindx(SY,0),mindx(SY,0)) = par[J1Y];
  J(ma1,mindx(SZ,0),mindx(SZ,0)) = par[J1Z];

  J(ma2,mindx(SX,0),mindx(SX,0)) = par[J1X];
  J(ma2,mindx(SY,0),mindx(SY,0)) = par[J1Y];
  J(ma2,mindx(SZ,0),mindx(SZ,0)) = par[J1Z];

  J(  a1a2,mindx(SX,0),mindx(SX,0)) = par[J2X];
  J(  a1a2,mindx(SY,0),mindx(SY,0)) = par[J2Y];
  J(  a1a2,mindx(SZ,0),mindx(SZ,0)) = par[J2Z];

  J( ma1a2,mindx(SX,0),mindx(SX,0)) = par[J2X];
  J( ma1a2,mindx(SY,0),mindx(SY,0)) = par[J2Y];
  J( ma1a2,mindx(SZ,0),mindx(SZ,0)) = par[J2Z];

  J(ma1ma2,mindx(SX,0),mindx(SX,0)) = par[J2X];
  J(ma1ma2,mindx(SY,0),mindx(SY,0)) = par[J2Y];
  J(ma1ma2,mindx(SZ,0),mindx(SZ,0)) = par[J2Z];

  J( a1ma2,mindx(SX,0),mindx(SX,0)) = par[J2X];
  J( a1ma2,mindx(SY,0),mindx(SY,0)) = par[J2Y];
  J( a1ma2,mindx(SZ,0),mindx(SZ,0)) = par[J2Z];
}
#endif //SPINISOTROPIC

#elif defined SQUAREXY

enum latticev{O,a1,a2,ma1,ma2,a1a2,ma1a2,ma1ma2,a1ma2};
const auto NC=9;
vector<Triplet> clist{{ 0, 0, 0}, { 1,  0,  0 }, {0,  1,  0 }, {-1, 0,  0 }, { 0, -1,  0 }, { 1, 1, 0}, {-1, 1, 0}, {-1,-1, 0},{ 1,-1, 0} };


vector<Coord> roffset={Coord(0,0,0)};
vector<double> invsqrtmasses={1.};

// global routines to extract sublattice and spin indices.

int mindx(const int spin,const int subl=0){ return spin*NSUBL+subl;}
int spin(const int m){return m/NSUBL;}
int subl(const int m){return m%NSUBL;}


#ifdef SPINISOTROPIC
void Couplings::InitializeJ()
{
  J( a1,mindx(SX,0),mindx(SX,0)) = par[J1X];
  J( a2,mindx(SX,0),mindx(SX,0)) = par[J1X];
  J(ma1,mindx(SX,0),mindx(SX,0)) = par[J1X];
  J(ma2,mindx(SX,0),mindx(SX,0)) = par[J1X];

  J(  a1a2,mindx(SX,0),mindx(SX,0)) = par[J2X];
  J( ma1a2,mindx(SX,0),mindx(SX,0)) = par[J2X];
  J(ma1ma2,mindx(SX,0),mindx(SX,0)) = par[J2X];
  J( a1ma2,mindx(SX,0),mindx(SX,0)) = par[J2X];
}
#else
void Couplings::InitializeJ()
{
  // 00 sublattices
  J(  O,mindx(SX,0),mindx(SX,0)) = par[AX];
  J(  O,mindx(SY,0),mindx(SY,0)) = par[AY];

  J( a1,mindx(SX,0),mindx(SX,0)) = par[J1X];
  J( a1,mindx(SY,0),mindx(SY,0)) = par[J1Y];

  J( a2,mindx(SX,0),mindx(SX,0)) = par[J1X];
  J( a2,mindx(SY,0),mindx(SY,0)) = par[J1Y];

  J(ma1,mindx(SX,0),mindx(SX,0)) = par[J1X];
  J(ma1,mindx(SY,0),mindx(SY,0)) = par[J1Y];

  J(ma2,mindx(SX,0),mindx(SX,0)) = par[J1X];
  J(ma2,mindx(SY,0),mindx(SY,0)) = par[J1Y];

  J(  a1a2,mindx(SX,0),mindx(SX,0)) = par[J2X];
  J(  a1a2,mindx(SY,0),mindx(SY,0)) = par[J2Y];

  J( ma1a2,mindx(SX,0),mindx(SX,0)) = par[J2X];
  J( ma1a2,mindx(SY,0),mindx(SY,0)) = par[J2Y];

  J(ma1ma2,mindx(SX,0),mindx(SX,0)) = par[J2X];
  J(ma1ma2,mindx(SY,0),mindx(SY,0)) = par[J2Y];

  J( a1ma2,mindx(SX,0),mindx(SX,0)) = par[J2X];
  J( a1ma2,mindx(SY,0),mindx(SY,0)) = par[J2Y];
}
#endif //SPINISOTROPIC

#endif // 

#endif // MODELCOUPLINGS_H
