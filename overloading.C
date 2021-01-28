#include<iostream>
#include<complex>
#include<vector>
#include<algorithm>

using namespace std;

const int NSUBL=NSUBLATTICES;
const int NSUBL2=NSUBL*NSUBL;

#include "qsmatrix.h"

main()
{
  QSmatrix C(16);

  C[0]=complex<double>(1,1);
  C[13]=complex<double>(1,2);
  C[14]={3,5};

  cout << C << endl;

  cout << C(1,1)[0] << endl;
}




