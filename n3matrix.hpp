#ifndef __N3_MATRIX_HPP__
#define __N3_MATRIX_HPP__

#include <iostream>
#include "n3dynamicarray.hpp"
#include <cmath>

using namespace std;

#define size_r ('R')
#define size_c ('C')

#define TOSTR(Var) # Var
#define pmat(Mat) cout << endl; cout << TOSTR(Mat) << " =" << endl; cout << Mat << endl

#define kroneckerDelta(i,j) (i==j?1.0:0.0)

class N3Matrix{
private:
  N3DynamicArray<double> elem;
  int             R , C , ptr;
  int    c2s( int , int );
public:
  const static double RANK_NORM_EPS=1e-13;
  N3Matrix();
  N3Matrix( int , int );
  void   allocate( int , int );
  int    size( char );
  void   set( int , int , double );
  double get( int , int );
  void   set( int , double );
  double get( int );
  bool   iszero( double );
  void   zero( void );
  void   eye( void );
  double norm( void );
  N3Matrix normalize( void );
  
  N3Matrix simp( N3Matrix* );
  N3Matrix simp( void );
  N3Matrix ge( N3Matrix& );
  int      rank( void );
  N3Matrix inv( void );
  N3Matrix pinv( void );
  N3Matrix opmat( void );
  
  N3Matrix gso( void );
  void     QRdec( N3Matrix * , N3Matrix * );

  N3Matrix getRowVec( int );
  N3Matrix getColVec( int );
  void setRowVec( int , N3Matrix );
  void setColVec( int , N3Matrix );
  
  //Operator ovarload
  N3Matrix  operator+( N3Matrix );
  N3Matrix  operator-( N3Matrix );
  N3Matrix  operator*( N3Matrix );
  N3Matrix& operator=( N3Matrix );
  N3Matrix& operator<<( double );
  N3Matrix& operator<<( int );
};

ostream& operator<<( ostream& , N3Matrix );
bool     operator==( N3Matrix& , N3Matrix& );
bool     operator!=( N3Matrix& , N3Matrix& );
N3Matrix operator*( double , N3Matrix );
N3Matrix operator/( N3Matrix , double );
N3Matrix operator+( N3Matrix );

#endif
