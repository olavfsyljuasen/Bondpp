#include<iomanip>
#include<vector>
//#include<complex>
#include<iostream>

using namespace std;

//typedef long double realtype;
//typedef complex<double> complex_type; 
//typedef double realtype;
//typedef float realtype;

#include <boost/multiprecision/cpp_complex.hpp>
#include <boost/multiprecision/eigen.hpp>
#include <Eigen/Dense>
using namespace Eigen;


// compile using g++ -std=c++11 -I/scratch/sylju/eigen-3.3.7 -I/scratch/sylju/boost_1_73_0 dettest.C; ./a.out
typedef boost::multiprecision::cpp_complex_quad complex_type;
//typedef complex<double> complex_type;

//typedef boost::multiprecision::cpp_complex_oct complex_type;





int main()
{
  vector<complex_type> a(16);

  a[ 0] = complex_type(23.206819701248566190088240546174347,0);
  a[ 1] = complex_type(22.900115826543057551134552340954542,0);
  a[ 2] = complex_type(22.921520925654476741328835487365723,0);
  a[ 3] = complex_type(22.921520925654476741328835487365723,0);
  a[ 4] = complex_type(24.060767841830951851989084389060736,0);
  a[ 5] = complex_type(23.206819701248566190088240546174347,0);
  a[ 6] = complex_type(23.214035974890844471474338206462562,0);
  a[ 7] = complex_type(23.214035974890848024188017006963491,0);
  a[ 8] = complex_type(23.214035974890848024188017006963491,0);
  a[ 9] = complex_type(22.921520925654483846756193088367581,0);
  a[10] = complex_type(23.033257269740005312996800057590008,0);
  a[11] = complex_type(22.952954271614647296928524156101048,0);
  a[12] = complex_type(23.214035974890848024188017006963491,0); 
  a[13] = complex_type(22.921520925654483846756193088367581,0);
  a[14] = complex_type(22.952954271614643744214845355600119,0); 
  a[15]=  complex_type(23.033257269739983996714727254584432,0);


  //  Eigen::Matrix< complex_type,4,4 > M; 
  //Eigen::Matrix< complex_type,4,4 > Minv; 
  Eigen::Matrix<complex_type,Dynamic,Dynamic > M(4,4);  // This is better than the fixed size
  Eigen::Matrix<complex_type,Dynamic,Dynamic > Minv(4,4);  // This is better than the fixed size

  int s=0;
  for(int i=0; i<4; i++)
    for(int j=0; j<4; j++)
      {
	M(i,j)=a[s++];
      }

  complex_type det=M.determinant();


  complex_type a11=M(0,0);
  complex_type a12=M(0,1);
  complex_type a13=M(0,2);
  complex_type a14=M(0,3);
  complex_type a21=M(1,0);
  complex_type a22=M(1,1);
  complex_type a23=M(1,2);
  complex_type a24=M(1,3);
  complex_type a31=M(2,0);
  complex_type a32=M(2,1);
  complex_type a33=M(2,2);
  complex_type a34=M(2,3);
  complex_type a41=M(3,0);
  complex_type a42=M(3,1);
  complex_type a43=M(3,2);
  complex_type a44=M(3,3);


  complex_type  det_brute= 
            a14*a23*a32*a41- a13*a24*a32*a41- a14*a22*a33*a41+ a12*a24*a33*a41
           +a13*a22*a34*a41- a12*a23*a34*a41- a14*a23*a31*a42+ a13*a24*a31*a42
           +a14*a21*a33*a42- a11*a24*a33*a42- a13*a21*a34*a42+ a11*a23*a34*a42
           +a14*a22*a31*a43- a12*a24*a31*a43- a14*a21*a32*a43+ a11*a24*a32*a43
           +a12*a21*a34*a43- a11*a22*a34*a43- a13*a22*a31*a44+ a12*a23*a31*a44
            +a13*a21*a32*a44- a11*a23*a32*a44- a12*a21*a33*a44+ a11*a22*a33*a44;



  cout << "Storing input array on matrix M" << endl;
  cout << "M= \n" << setprecision(35) << M << endl;


  cout << "Determinant: " << endl;
  cout << setprecision(50) << det << " (using Eigen)" << endl;
  cout << "-0.099428409859066000615714468691597266326899171503636 (Mathematica)" << endl;
  cout << setprecision(50) << det_brute << " (Brute force method)" << endl;


  Minv=M.inverse();

  cout << "The inverse matrix:" << endl;
  cout << Minv(0,0) << " (Minv(0,0)) " << endl;
  cout << "-2.4162335384522813609672912614026554914457664306740 (Minv(0,0) Mathematica" << endl;
  cout << Minv(0,1) << " (Minv(0,1)) " << endl;
  cout << "1.8554284099166845044815163290756499357040235989045 (Minv(0,1) Mathematica)" << endl;

  



}
