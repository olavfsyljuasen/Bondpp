#include<vector>
#include<complex>
#include<iostream>

using namespace std;

#include<Eigen/Dense>
using namespace Eigen;

// compile using: g++ -I/scratch/sylju/eigen-3.3.7 eigentest.C; ./a.out

int main()
{
  vector<complex<double>> a(4);

  a[0]=complex<double>(2,0);
  a[1]=complex<double>(0,-1);
  a[2]=complex<double>(0,+1);
  a[3]=complex<double>(1,0);

  cout << "Input array: " << endl;
  for(int i=0; i<4; i++)
    {
      cout << "a[" << i << "]="  << a[i] <<  endl;
    }

  // Note the difference between using Map and assigning matrix elements directly
  Map<Matrix<complex<double>,2,2>,RowMajor> M(&a[0]); // RowMAjor or ColMajor makes NO difference
  /* using indeces to initialize to RowMajor order:
    Matrix<complex<double>,2,2> M;
    M(0,0)=a[0];
    M(0,1)=a[1];
    M(1,0)=a[2];
    M(1,1)=a[3];
  */

  cout << "Storing input array on matrix M" << endl;
  cout << "M= \n" << M << endl;

  cout << endl;
  cout << "Writing out M using M.data(): " << endl;
  
  for(int i=0; i<4; i++)
    {
      cout << "M.data()[" << i << "]="  << M.data()[i] <<  endl;
    }

  cout << endl;

  cout << "Accessing M using matrix elements: " << endl;
  for(int i=0; i<2; i++)
    for(int j=0; j<2; j++)
      {
	cout << "M(" << i << "," << j << ")=" << M(i,j) << endl;
      }

  cout << "-------------------" << endl;
  
  Matrix<complex<double>,2,2> Minv=M.inverse();

  cout << "Minv= \n" << Minv <<  endl;

  cout << endl;
  for(int i=0; i<4; i++)
    {
      cout << "Minv.data()[" << i << "]="  << Minv.data()[i] <<  endl;
    }

  cout << endl;

  
  
  for(int i=0; i<2; i++)
    for(int j=0; j<2; j++)
      {
	cout << "Minv(" << i << "," << j << ")=" << Minv(i,j) << endl;
      }


  int array[8];
  for(int i = 0; i < 8; ++i) array[i] = i;
  cout << "Column-major:\n" << Map<Matrix<int,2,4> >(array) << endl;
  cout << "Row-major:\n" << Map<Matrix<int,2,4,RowMajor> >(array) << endl;
  cout << "Row-major using stride:\n" <<
    Map<Matrix<int,2,4>, Unaligned, Stride<1,4> >(array) << endl;
}
