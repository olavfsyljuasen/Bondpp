#include<string>
#include<iostream>
#include<complex>

//#define _USE_MATH_DEFINES
//#include<cmath>

#include<stdio.h>
#include<stdlib.h>
using namespace std;

#include "mynumbertypes.h"

#define _USE_MATH_DEFINES
#include<cmath>

#define KMAJORORDER   // first traverse k's, then next matrix element
//#undef  KMAJORORDER // first traverse matrix, then next k.

#ifdef KMAJORORDER   // first traverse k's, then next matrix element
constexpr int KSTRIDE=1;
#else                // first traverse matrix, then next k.
constexpr int KSTRIDE= NSUBL*NSUBL;
#endif

int ir(0); // inntrykk variabel, for setting innrykk i TRACELEVEL printouts
string spaces(const int n){ return string(n,' ');}


#include "overload.h"

#define watch(x) cout << (#x) << " is " << (x) << endl


// Options for this run
#if defined NDEBUG 
const bool TRACE= false; 
const bool TRACE2=false; // more debug info
const int TRACELEVEL=-1;
#else
const bool TRACE= true; // debugging option
const bool TRACE2=true; // more debug info
const int  TRACELEVEL= 4;     // 0: No info, 1: Location 2: Main results 3: Details, 4: Details that also takes a lot of space. 5: Unimportant details 
//const bool TRACE= false; // debugging option
#endif

const bool TRACEFREEENERGY = false; // print out individual terms to the free energy

// options regarding formatting in logfile, debugging and outputfiles
#if defined SETLOGPRECISION
const int LOGPRECISION  = SETLOGPRECISION; // number of digits in logfile info
#else
const int LOGPRECISION  = 2; // number of digits in logfile info
#endif

const int DEBUGPRECISION =6; // number of digits to use in debug-info
const int OUTPUTPRECISION=12; // number of digits to use in outfiles


const bool dcheck=false; // check if D is 0 or negative and output warning

const int  MAXRUNS = 200; // the maximum number of runs in the read.in file  
const int  MAXPAR = 100;  // the maximum allowed input parameters


const int Nlargestqvalues=6; // keep track of the locations of the ... susc peaks. 

const bool SHOWSUBTRACTED = true; // true to write subtractions to logfile when printino=true

const bool PRINTPROGRESS        = true; // true to write T and epsilons to logfile
int  PRINTPROGRESSTICKLER = 10;   // print every TICKLER steps, can be modified

#ifdef PRINTPHONONS
const int PRINTPHONONSPECTRUM = true;
#else
const int PRINTPHONONSPECTRUM = false;
#endif
const int PRINTPHONONSPECTRUMTICKLER = 10; // print renormalized phonon spectrum for every ..TICKLER steps.

const bool BINARYOUTFILES=false;

#ifdef USELASTSIGMA
const bool USEPREVIOUSEPSILONS = true; // use previous epsilon values whenever last run converged.
#else
const bool USEPREVIOUSEPSILONS = false;
#endif

const bool PRINTGQS = true; // output the couplings to files

const string STATEFILENAME ="conf"; 

const string RCORRSFILENAME="rcorrs.dat";
const string QCORRSFILENAME="qcorrs";
const string SIGMAEFILENAME="sigmaEq.dat";

const string LARGESTQCORRSFILENAME="qcorrs_sparse";

const string KINVQFILENAME="Kinvq.dat"; // the last converged Kinvq
const string INPUTKINVQFILENAME="Kinvq.in"; // initializing Kinvq


const string RCORRMOSTREMOTE="rcorr_L2.dat";

const string SITERPTS="siteRpts.dat";
const string SITEQPTS="siteQpts.dat";


const int MAXNSELECTEDQPTS=20;
const string SELECTEDQPTSFILENAME="selectedqpts.in";

const string DELTASTOSHOWFILENAME="Deltastoshow.in";
const string INPUTEPSILONFILENAME="epsilons.in";
const string EPSILONFILENAME="epsilons.dat";


const bool PRINTMAXQS=true; // print out the q values of the NMAXQS largest susceptibility peaks 
const int  NMAXQS = 32;
const int  NMAXQSINFIRSTBZ=16;
const realtype QTOL         = 1.e-7; // if two q-values are closer than this, they are regarded equal.


#if defined CLAMPEDMODES
const int clampedmodes[]={CLAMPEDMODES};
const int Nclampedmodes=sizeof(clampedmodes)/sizeof(clampedmodes[0]);
#endif

#ifdef SHOWCORRELATIONFUNCTIONS
const bool PRINTRCORRS=true;
const bool PRINTQCORRS=true;
const bool PRINTSIGMAE=false;
const bool PRINTRPTS=false;
const bool PRINTQPTS=true;
const bool PRINTLARGESTQCORRS=true;
const realtype QCORRSTHRESHOLD=realtype(0.); // print QCORRS that are above maxval*QCORRSTHRESHOLD
#if defined SETPRINTTICKLER
const int PRINTTICKLER=SETPRINTTICKLER;
#else
const int PRINTTICKLER=10;
#endif
#else
const bool PRINTRCORRS=false;
const bool PRINTQCORRS=false;
const bool PRINTSIGMAE=false;
const bool PRINTRPTS=false;
const bool PRINTQPTS=false;
const bool PRINTLARGESTQCORRS=false;
const realtype QCORRSTHRESHOLD=realtype(0.);
const int PRINTTICKLER=1; 
#endif


const int MAXNINCREASES = 10; // used to be very small = 5.
const int ITERMULTIPLIER = 10; // used to increase the number of iterations if convergence fails.

const bool PRINTNORMALMODES=true;
const bool PRINTPHONONENERGIES=true;

const string CONVERGENCEMONITORNAME="monitor.dat";  
const string RESULTSNAME= "res.dat";
const string RES1NAME= "res1.dat";
const string RES2NAME= "res2.dat";

const string RESULTS2NAME= "dt.dat";
const string RESULTS3NAME= "td.dat";
const string RESULTS4NAME= "m2.dat";
const string RESULTS5NAME= "chiq.dat";

const string MAXQNAME="maxq.dat";
const string MAXQSINFIRSTBZNAME="maxq1bz.dat";
const string PEAKREPORTNAME= "peaks.dat";

const string PARAMETERFILENAME = "Deltas.in";
const string ZONEORIGINSFILENAME= "zoneorigins.in"; // center of zones in units of reciprocal lattice vectors b1,b2,b3
const string READIN   = "read.in";

#include <fstream>
ofstream logfile("log.txt",ios::app);

//#include "rnddef.h" // the random number generator RAN


#include "globalheader.h"

#include "inputparameters.h" // should be modified for new systems

#include "RunParameter.h"

RunParameters rp; // declare parameters globally
realtype* par = rp.GetPars(1); // only the first line is read

#include "bravaislattices.h"

BravaisLattice la(par); // declare lattice globally 



#include "vecmat.h"

#include "modeldef.h" // set number of sublattices, and other compile-time parameters, modify it for new models

#include "couplings.h" // generic coupling class

#include "modelcouplings.h" // set couplings, modify for new models



#ifdef LATTICEDISTORTIONS
#include "phonons.h" // must set the springs, modify for different phonons
#endif

#include "symmetryroutines.h" // for enforcing and checking symmetries

#include "observables.h" // observables, modify it for new models

#include "rules.h" // generic rules


#include "bondpp.h"



int main()
{
  if(TRACE) cout << "Starting Main()" << endl;
  Timer mytimer;  
  logfile << "\nStarting at: " << mytimer << endl;
  mytimer.Start();
  
  //  RunParameters rp;
  logfile << "parameters: " << rp;

#ifdef PRESERVESYMMETRY
  logfile << " PRESERVESYMMETRY option is on " << endl;
#endif
  
  logfile << "Starting calculation " << endl;



  Simulation sim;

  logfile << "# of nonzero gc's= " << nonzeroclist.size() << endl;
  sim.Run();

            
  mytimer.Stop();
  double time_spent = mytimer.GetTimeElapsed();
  logfile << "Time elapsed:  " << time_spent << " sec. :" 
	  << time_spent/3600. << " hours." << endl;
  logfile << "Ending at: " << mytimer;
  logfile << "---Done--- " << endl;
  logfile.close();
  
  ofstream terminationsign("end.exec");
  terminationsign << "execution ended at " << mytimer;
  terminationsign.close();
  
  cout << "program ended \n";
  exit(0);
}











