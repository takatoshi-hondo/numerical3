#ifndef __N3_COMPLEX_NUMBER_HPP__
#define __N3_COMPLEX_NUMBER_HPP__

#include <cmath>
#include <iostream>

using namespace std;

#define N3R ("N3_REAL")
#define N3I ("N3_IMAG")

class N3Complex {
private:
  double re , im;
  int    cnt;
public:
  const static double ZERO_EPS=1e-13;
  N3Complex();
  double    getReal( void );
  double    getImag( void );
  bool      isreal( void );
  bool      isimag( void );
  N3Complex conj( void );
  double    abs( void );
  double    angle( void );
  void      polar2comp( double , double );
  void      operator()( double , double );
  N3Complex operator<<( double );
  double&   operator[]( const char * );
  N3Complex operator=( double );
  N3Complex operator+( N3Complex );
  N3Complex operator-( N3Complex );
  N3Complex operator*( N3Complex );
  N3Complex operator/( double );
  N3Complex operator/( N3Complex );
};

ostream&  operator<<( ostream& , N3Complex );
N3Complex sqrt( N3Complex );

#endif
