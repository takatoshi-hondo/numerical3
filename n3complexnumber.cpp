#include "n3complexnumber.hpp"

N3Complex::N3Complex(){
  cnt = 0;
  re = 0;
  im = 0;
}

double N3Complex::getReal( void ){
  return re;
}

double N3Complex::getImag( void ){
  return im;
}

bool N3Complex::isreal( void ){
  bool ret;
  if( im < ZERO_EPS && im > -ZERO_EPS ){
    ret = true;
  }else{
    ret = false;
  }
  return ret;
}

bool N3Complex::isimag( void ){
  bool ret;
  if( re < ZERO_EPS && re > -ZERO_EPS ){
    ret = true;
  }else{
    ret = false;
  }
  return ret;
}

N3Complex N3Complex::conj( void ){
  N3Complex ret;
  ret << re << -im;
  return ret;
}

double N3Complex::abs( void ){
  return sqrt( re*re + im*im );
}

double N3Complex::angle( void ){
  return atan2( im , re );
}

void N3Complex::polar2comp( double r , double angle ){
  re = r*cos(angle);
  im = r*sin(angle);
}

void N3Complex::operator()( double r , double i ){
  this->re = r;
  this->im = i;
}

N3Complex N3Complex::operator<<( double val ){
  if( cnt == 0 ){
    re = val;
    cnt = 1;
  }else{
    im = val;
    cnt = 0;
  }
  return *this;
}

double& N3Complex::operator[]( const char *arg ){
  if( arg == N3R ){
    return re;
  }else if( arg == N3I ){
    return im;
  }
}

N3Complex N3Complex::operator=( double val ){
  re = val;
  im = 0.0;
  cnt++;
  return *this;
}

N3Complex N3Complex::operator+( N3Complex q ){
  N3Complex ret;
  ret(q.getReal()+re,q.getImag()+im);
  return ret;
}

N3Complex N3Complex::operator-( N3Complex q ){
  N3Complex ret;
  ret(re-q.getReal(),im-q.getImag());
  return ret;
}

N3Complex N3Complex::operator*( N3Complex q ){
  N3Complex ret;
  ret(re*q.getReal()-im*q.getImag(),re*q.getImag()+im*q.getReal());
  return ret;
}

N3Complex N3Complex::operator/( double val ){
  N3Complex ret;
  ret(re/val,im/val);
  return ret;
}

N3Complex N3Complex::operator/( N3Complex q ){
  N3Complex ret;
  ret = (*this)*q.conj()/(q*q.conj()).getReal();
  return ret;
}

ostream& operator<<( ostream& ost , N3Complex q ){
  ost << q.getReal() << "+" << q.getImag() << "i";
  return ost;
}

N3Complex sqrt( N3Complex q ){
  N3Complex ret;
  double r , theta;
  r = q.abs();
  theta = q.angle();
  r = sqrt(r);
  theta = theta/2.0;
  ret.polar2comp(r,theta);
  return ret;
}
