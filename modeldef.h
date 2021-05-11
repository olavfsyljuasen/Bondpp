#ifndef MODELDEF_H
#define MODELDEF_H

using namespace std;


#if defined SQUAREPHONONS || defined TRIANGULARPHONONS || defined CUBICPHONONS || defined FCCPHONONS || defined BCCPHONONS

#ifdef FAKEHEISENBERG
const auto NSPIN= 1;
const auto NFAKESPINTRACE=3; // a compensating factor for treating Heisenberg models with a single spin comp.
#else
const auto NSPIN= 3;
const auto NFAKESPINTRACE=1; 
#endif // FAKEHEISENBERG
const auto NSUBL= 1;


#elif defined CrI3

const auto NELASTIC=0;

#ifdef SPINISOTROPIC
const auto NSPIN= 1;
const auto NFAKESPINTRACE=3; // a compensating factor for treating Heisenberg models with a single spin comp.
#else
const auto NSPIN= 3;
const auto NFAKESPINTRACE=1; // a compensating factor for treating Heisenberg models with a single spin comp.
#endif //SPINISOTROPIC

const auto NSUBL= NSUBLATTICES; // should be 2.


#elif defined SQUARENOPHONONS

const auto NELASTIC=0;

#ifdef SPINISOTROPIC
const auto NSPIN= 1;
const auto NFAKESPINTRACE=3; // a compensating factor for treating Heisenberg models with a single spin comp.
#else
const auto NSPIN= 3;
const auto NFAKESPINTRACE=1; // a compensating factor for treating Heisenberg models with a single spin comp.
#endif //SPINISOTROPIC

const auto NSUBL= NSUBLATTICES; // should be 1.


#elif defined SQUAREXY

const auto NELASTIC=0;

#ifdef SPINISOTROPIC
const auto NSPIN= 1;
const auto NFAKESPINTRACE=2; // a compensating factor for treating Heisenberg models with a single spin comp.
#else
const auto NSPIN= 2;
const auto NFAKESPINTRACE=1; // a compensating factor for treating Heisenberg models with a single spin comp.
#endif //SPINISOTROPIC

const auto NSUBL= NSUBLATTICES; // should be 1.


#endif // 

const auto NMAT = NSPIN*NSUBL;
const auto NMAT2= NMAT*NMAT;

#endif // MODELDEF_H
