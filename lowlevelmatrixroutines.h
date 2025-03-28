#include<vector>
#include<complex>
#include<iostream>
#include<limits>

#include "mynumbertypes.h"

void SmallMatrixMakeHermitian(const int n, vector<complextype>& A);
void SmallMatrixInverse(const int n,vector<complextype>& A);
realtype SmallHermitianMatrixMinEigenvalue(const int n, vector<complextype>& A);
realtype SmallHermitianMatrixDeterminant(const int n, vector<complextype>& A);
vector<realtype> SmallMatrixPseudoSolve(vector<realtype>& Ain, vector<realtype>& Bin,const int M, const int N);
