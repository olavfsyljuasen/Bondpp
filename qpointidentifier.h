#ifndef POINTIDENTIFIER_H
#define POINTIDENTIFIER_H
#include<iostream>
#include<string>


#include "globalheader.h"



class Qlocation
{
public:
  Qlocation(){}
  virtual ~Qlocation(){}
  virtual bool isinside(Coord q)=0;  
};


class Point : public Qlocation
{
public:
  Point(Coord p_in):p(p_in){}
  Point(const Point& p_in):p(p_in.p){}
  ~Point(){}
  bool isinside(Coord q){return( (q-p).Norm()< QTOL );}
  Coord p;
};

class Line : public Qlocation
{
public:
  Line(Coord p1_in,Coord p2_in):p1(p1_in),p2(p2_in),v(p2-p1){ v.Normalize();}
  Line(const Point& p1_in,const Point& p2_in):p1(p1_in.p),p2(p2_in.p),v(p2-p1){v.Normalize();}
  ~Line(){}
  bool isinside(Coord q) // q is inside if q-p1 is paralell to v;
  { 
    return ( crossproduct(q-p1,v).Norm()< QTOL);
  }
  
  Coord p1;
  Coord p2;
  Coord v;  // the vector defining the line p2-p1
};


class Surface : public Qlocation
{
public:
  Surface(Coord p1_in,Coord p2_in,Coord p3_in):p1(p1_in),p2(p2_in),p3(p3_in),nv(crossproduct(p2-p1,p3-p1)){nv.Normalize();}
  Surface(const Point& p1_in,const Point& p2_in,const Point& p3_in):p1(p1_in.p),p2(p2_in.p),p3(p3_in.p),nv(crossproduct(p2-p1,p3-p1)){ nv.Normalize();}
  ~Surface(){}
  bool isinside(Coord q) // q is inside if q-p1 is normal to nv
  {
    return ( abs(scalarproduct(q-p1,nv)) < QTOL );
  }
  Coord p1;
  Coord p2;
  Coord p3;
  Coord nv;  // the normal vector
};



class Qids
{
public:
  Qids();
  ~Qids(){for(long unsigned int i=0; i<sqpt.size(); i++) delete(sqpt[i]);}
  bool Findid(Coord q,string& qid)
  {
    //    cout << "q in=" << q << " " << endl;
    qid ="? ";
    for(long unsigned int i=0; i<sqpt.size(); i++)
      {
	if( sqpt[i]->isinside(symmetrytransform(q)) ){ qid=sqid[i]; return true; }
      }

    return false;
  }
  
private:
  vector<Qlocation*> sqpt;
  vector<string> sqid;
  Coord symmetrytransform(Coord);
};

#if defined FACECENTEREDCUBICBRAVAISLATTICE
Qids::Qids():sqpt(34),sqid(34)
{
  // Transform q to the positive octant first, then these are the points and lines to look for:
  // first points, then lines, and then surfaces last. sqid is a string of max 3 char.
  sqpt[ 0] = new Point(Coord(     0,     0,     0)); sqid[ 0] = "\u0393";  // Gamma
 
  sqpt[ 1] = new Point(Coord(    PI,    PI,    PI)); sqid[ 1] = "L ";

  sqpt[ 2] = new Point(Coord( TWOPI,     0,     0)); sqid[ 2] = "X ";
  sqpt[ 3] = new Point(Coord(     0, TWOPI,     0)); sqid[ 3] = "X ";
  sqpt[ 4] = new Point(Coord(     0,     0, TWOPI)); sqid[ 4] = "X ";

  sqpt[ 5] = new Point(Coord(    PI, TWOPI,     0)); sqid[ 5] = "W ";
  sqpt[ 6] = new Point(Coord(     0,    PI, TWOPI)); sqid[ 6] = "W ";
  sqpt[ 7] = new Point(Coord( TWOPI,     0,    PI)); sqid[ 7] = "W ";
  sqpt[ 8] = new Point(Coord( TWOPI,    PI,     0)); sqid[ 8] = "W ";
  sqpt[ 9] = new Point(Coord(     0, TWOPI,    PI)); sqid[ 9] = "W ";
  sqpt[10] = new Point(Coord(    PI,     0, TWOPI)); sqid[10] = "W ";

  sqpt[11] = new Point(Coord( HLFPI, HLFPI, TWOPI)); sqid[11] = "U ";
  sqpt[12] = new Point(Coord( HLFPI, TWOPI, HLFPI)); sqid[12] = "U ";
  sqpt[13] = new Point(Coord( TWOPI, HLFPI, HLFPI)); sqid[13] = "U ";

  sqpt[14] = new Point(Coord(     0, TRHPI, TRHPI)); sqid[14] = "K ";
  sqpt[15] = new Point(Coord( TRHPI,     0, TRHPI)); sqid[15] = "K ";
  sqpt[16] = new Point(Coord( TRHPI, TRHPI,     0)); sqid[16] = "K ";

  sqpt[17] = new Line( Coord(    0,    0,    0), Coord(   PI,   PI,   PI)); sqid[17] ="\u0393L";

  sqpt[18] = new Line( Coord(    0,    0,    0), Coord(TWOPI,    0,    0)); sqid[18] ="\u0393X";
  sqpt[19] = new Line( Coord(    0,    0,    0), Coord(    0,TWOPI,    0)); sqid[19] ="\u0393X";
  sqpt[20] = new Line( Coord(    0,    0,    0), Coord(    0,    0,TWOPI)); sqid[20] ="\u0393X";

  sqpt[21] = new Line( Coord(TWOPI,    0,    0), Coord(TWOPI,HLFPI,HLFPI)); sqid[21] ="XU ";
  sqpt[22] = new Line( Coord(    0,TWOPI,    0), Coord(HLFPI,TWOPI,HLFPI)); sqid[22] ="XU ";
  sqpt[23] = new Line( Coord(    0,    0,TWOPI), Coord(HLFPI,HLFPI,TWOPI)); sqid[23] ="XU ";

  sqpt[24] = new Line( Coord(TWOPI,    0,    0), Coord(TWOPI,   PI,    0)); sqid[24] ="XW ";
  sqpt[25] = new Line( Coord(TWOPI,    0,    0), Coord(TWOPI,    0,   PI)); sqid[25] ="XW ";
  sqpt[26] = new Line( Coord(    0,TWOPI,    0), Coord(   PI,TWOPI,    0)); sqid[26] ="XW ";
  sqpt[27] = new Line( Coord(    0,TWOPI,    0), Coord(    0,TWOPI,   PI)); sqid[27] ="XW ";
  sqpt[28] = new Line( Coord(    0,    0,TWOPI), Coord(   PI,    0,TWOPI)); sqid[28] ="XW ";
  sqpt[29] = new Line( Coord(    0,    0,TWOPI), Coord(    0,   PI,TWOPI)); sqid[29] ="XW ";

  sqpt[30] = new Surface(Coord(TWOPI,    0,    0), Coord(TWOPI,   PI,    0), Coord(TWOPI,HLFPI,HLFPI)); sqid[30] ="Sqr ";
  sqpt[31] = new Surface(Coord(    0,TWOPI,    0), Coord(    0,TWOPI,   PI), Coord(HLFPI,TWOPI,HLFPI)); sqid[31] ="Sqr ";
  sqpt[32] = new Surface(Coord(    0,    0,TWOPI), Coord(   PI,    0,TWOPI), Coord(HLFPI,HLFPI,TWOPI)); sqid[32] ="Sqr ";

  sqpt[33] = new Surface(Coord(   PI,   PI,   PI), Coord(TRHPI,TRHPI,    0), Coord(HLFPI,HLFPI,TWOPI)); sqid[33] ="Hex ";
}

Coord Qids::symmetrytransform(Coord q)
{
  return Coord(abs(q.x),abs(q.y),abs(q.z)); // transform into the positive octant.
}

#elif defined SIMPLECUBICBRAVAISLATTICE
Qids::Qids():sqpt(24),sqid(24)
{
  Point  G(Coord(    0,    0,   0));
  Point  R(Coord(   PI,   PI,  PI));
  Point X1(Coord(   PI,    0,   0));
  Point X2(Coord(    0,   PI,   0));
  Point X3(Coord(    0,    0,  PI));
  Point M1(Coord(    0,   PI,  PI));
  Point M2(Coord(   PI,    0,  PI));
  Point M3(Coord(   PI,   PI,   0));

  sqpt[ 0] = new Point(G);  sqid[ 0] = "\u0393";  // Gamma
  sqpt[ 1] = new Point(R);  sqid[ 1] = "R";
  sqpt[ 2] = new Point(X1); sqid[ 2] = "X";
  sqpt[ 3] = new Point(X2); sqid[ 3] = "X";
  sqpt[ 4] = new Point(X3); sqid[ 4] = "X";
  sqpt[ 5] = new Point(M1); sqid[ 5] = "M";
  sqpt[ 6] = new Point(M2); sqid[ 6] = "M";
  sqpt[ 7] = new Point(M3); sqid[ 7] = "M";

  sqpt[ 8] = new Line( G, R ); sqid[ 8] ="\u0393R"; // Lambda
  sqpt[ 9] = new Line( G, X1); sqid[ 9] ="\u0393X"; // Delta
  sqpt[10] = new Line( G, X2); sqid[10] ="\u0393X"; // Delta
  sqpt[11] = new Line( G, X3); sqid[11] ="\u0393X"; // Delta
  sqpt[12] = new Line( G, M1); sqid[12] ="\u0393M"; // Sigma
  sqpt[13] = new Line( G, M2); sqid[13] ="\u0393M"; // Sigma
  sqpt[14] = new Line( G, M3); sqid[14] ="\u0393M"; // Sigma
  sqpt[15] = new Line( R, X1); sqid[15] ="RX";      // S
  sqpt[16] = new Line( R, X2); sqid[16] ="RX";      // S
  sqpt[17] = new Line( R, X3); sqid[17] ="RX";      // S
  sqpt[18] = new Line( R, M1); sqid[18] ="RM";      // T
  sqpt[19] = new Line( R, M2); sqid[19] ="RM";      // T
  sqpt[20] = new Line( R, M3); sqid[20] ="RM";      // T
  
  sqpt[21] = new Surface( R, X1 , M2); sqid[21] = "Edg";
  sqpt[22] = new Surface( R, X2 , M3); sqid[22] = "Edg";
  sqpt[23] = new Surface( R, X3 , M1); sqid[23] = "Edg";
  
}

Coord Qids::symmetrytransform(Coord q)
{
  return Coord(abs(q.x),abs(q.y),abs(q.z)); // transform into the positive octant.
}
#elif defined HEXAGONALBRAVAISLATTICE
Qids::Qids():sqpt(29),sqid(29)
{
  // Transform q to the positive octant first, then these are the points and lines to look for:
  // first points, then lines, and then surfaces last. sqid is a string of max 3 char.

  const realtype FTHPI=4.*PI/3.;
  const realtype TTHPI=TWOPI/3.;
  const realtype TR3PI=TWOPI/sqrt(3.);
  const realtype  R3PI=   PI/sqrt(3.);

  Point  G(Coord(    0,    0,   0));
  Point  K(Coord(FTHPI,    0,   0));
  Point Kp(Coord(TTHPI,TR3PI,   0));
  Point Mp(Coord(   PI, R3PI,   0));
  Point  M(Coord(    0,TR3PI,   0));

  Point  A(Coord(    0,    0,  PI));
  Point  H(Coord(FTHPI,    0,  PI));
  Point Hp(Coord(TTHPI,TR3PI,  PI));
  Point Lp(Coord(   PI, R3PI,  PI));
  Point  L(Coord(    0,TR3PI,  PI));
  
  sqpt[ 0] = new Point(G ); sqid[ 0] = "\u0393";

  sqpt[ 1] = new Point(A ); sqid[ 1] = "A";

  sqpt[ 2] = new Point(K ); sqid[ 2] = "K";
  sqpt[ 3] = new Point(Kp); sqid[ 3] = "K";
  sqpt[ 4] = new Point(M ); sqid[ 4] = "M";
  sqpt[ 5] = new Point(Mp); sqid[ 5] = "M";

  sqpt[ 6] = new Point(H ); sqid[ 6] = "H";
  sqpt[ 7] = new Point(Hp); sqid[ 7] = "H";
  sqpt[ 8] = new Point(L ); sqid[ 8] = "L";
  sqpt[ 9] = new Point(Lp); sqid[ 9] = "L";

  sqpt[10] = new Line( G, K ); sqid[10] ="\u0393K";
  sqpt[11] = new Line( G, Kp); sqid[11] ="\u0393K";
  sqpt[12] = new Line( G, M ); sqid[12] ="\u0393M";
  sqpt[13] = new Line( G, Mp); sqid[13] ="\u0393M";

  sqpt[14] = new Line( A, H ); sqid[14] ="AH";
  sqpt[15] = new Line( A, Hp); sqid[15] ="AH";
  sqpt[16] = new Line( A, L ); sqid[16] ="AL";
  sqpt[17] = new Line( A, Lp); sqid[17] ="AL";

  sqpt[18] = new Line( K, Mp); sqid[18] ="KM";
  sqpt[19] = new Line(Mp, Kp); sqid[19] ="KM";
  sqpt[20] = new Line(Kp, M ); sqid[20] ="KM";
  sqpt[21] = new Line( H, Lp); sqid[21] ="HL";
  sqpt[22] = new Line(Lp, Hp); sqid[22] ="HL";
  sqpt[23] = new Line(Hp, L ); sqid[23] ="HL";

  sqpt[24] = new Surface( K,Mp, H); sqid[24] ="Zed"; // Zone edge surface.
  sqpt[25] = new Surface(Mp,Kp,Lp); sqid[25] ="Zed";
  sqpt[26] = new Surface(Kp, M,Hp); sqid[26] ="Zed";
  sqpt[27] = new Surface( G, K, M); sqid[27] ="Zbo"; // Zone planar bottom
  sqpt[28] = new Surface( A, H, L); sqid[28] ="Zto"; // Zone planar top  
}

Coord Qids::symmetrytransform(Coord q)
{
  return Coord(abs(q.x),abs(q.y),abs(q.z)); // transform into the positive octant.
}
#elif defined BODYCENTEREDCUBICBRAVAISLATTICE
Qids::Qids():sqpt(1),sqid(1)
{
 // NOT IMPLEMENTED
  sqpt[ 0] = new Point(Coord(     0,     0,     0)); sqid[ 0] = "\u0393";  // Gamma
}
Coord Qids::symmetrytransform(Coord q)
{
  return Coord(abs(q.x),abs(q.y),abs(q.z)); // transform into the positive octant.
}
#endif

#endif
