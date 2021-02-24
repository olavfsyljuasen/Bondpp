#ifndef INPUTPARAMETERS_H
#define INPUTPARAMETERS_H

using namespace std;

// The input parameter sequence

#ifdef SQUAREPHONONS

const int NPARAMS = 28;
enum params{LINEID,J1X,J1Y,J1Z,J2X,J2Y,J2Z,J3X,J3Y,J3Z,g1X,g1Y,g1Z,g2X,g2Y,g2Z,ALPHA1,ALPHA2,NX,NY,NZ,DA,DELTA,MU,MAXITER,TOLERANCE,NBINS,EQFLAG};

#elif defined CUBICPHONONS

const int NPARAMS = 28;
enum params{LINEID,J1X,J1Y,J1Z,J2X,J2Y,J2Z,J3X,J3Y,J3Z,g1X,g1Y,g1Z,g2X,g2Y,g2Z,ALPHA1,ALPHA2,NX,NY,NZ,DA,DELTA,MU,MAXITER,TOLERANCE,NBINS,EQFLAG};

#elif defined CrI3

const int NPARAMS = 32;
enum params{LINEID,J1X,J1Y,J1Z,J2X,J2Y,J2Z,J3x,J3Y,J3Z,AX,AY,AZ,D1X,D1Y,D1Z,D2X,D2Y,D2Z,D3X,D3Y,D3Z,NX,NY,NZ,DA,DELTA,MU,MAXITER,TOLERANCE,NBINS,EQFLAG};

#endif

















#endif // INPUTPARAMETERS_H
