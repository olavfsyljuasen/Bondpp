#ifndef MYNUMBERTYPES_H
#define MYNUMBERTYPES_H

#include<complex>
#include<array>
using namespace std;




#ifdef LONGDOUBLE
typedef         long double  realtype;
typedef complex<long double> complextype;
#elif defined FLOAT
typedef         float  realtype;
typedef complex<float> complextype;
#elif defined QUAD
typedef         __float128  realtype; // does not work                                                                                                                                               
typedef complex<__float128> complextype; // does not work                                                                                                                                            
#else
typedef         double  realtype;
typedef complex<double> complextype;
#endif

typedef array<complextype,3> cmplxcoord;

#ifdef LONGDOUBLE
#define mycos(x) cosl((x))
#define mysin(x) sinl((x))
#define mysqrt(x) sqrtl((x))
#define mylog(x) logl((x))
#define myexp(x) expl((x))
#define expi(x) complextype(cosl((x)),sinl((x)))
#else
#define mycos(x) cos((x))
#define mysin(x) sin((x))
#define mysqrt(x) sqrt((x))
#define mylog(x) log((x))
#define myexp(x) exp((x))
#define expi(x) complextype(cos((x)),sin((x)))
#endif



#ifndef M_PI
#define M_PI 3.14159265358979323846264338328
#define M_PIl 3.141592653589793238462643383279502884L 

#endif

#ifndef SQRTTHREEOVERTWO
#ifdef LONGDOUBLE
#define SQRTTHREEOVERTWO  0.8660254037844386467637231707529362L
#else
#define SQRTTHREEOVERTWO  0.866025403784438646763723170753
#endif
#endif

#ifdef LONGDOUBLE
#define ROOTTHREE         1.73205080756887729352744634150587236L
#define ROOTTHREEOVERTWO  0.8660254037844386467637231707529362L
#define ROOTTHREEOVERFOUR 0.43301270189221932338186158537646809L
#define ROOTSEVEN         2.64575131106459059050161575363926043L
#else
#define ROOTTHREE         1.7320508075688772935274463415058
#define ROOTTHREEOVERTWO  0.866025403784438646763723170753
#define ROOTTHREEOVERFOUR 0.4330127018922193233818615853765
#define ROOTSEVEN         2.6457513110645905905016157536393
#endif


// universal definitions
#ifdef LONGDOUBLE
const realtype PI=M_PIl;
#else
const realtype PI=M_PI;
#endif

const realtype TWOPI=2.*PI;
const realtype HLFPI=0.5*PI;
const realtype TRHPI=1.5*PI;





























#endif

