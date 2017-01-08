#include "n3matrix.hpp"

N3Matrix::N3Matrix(){
  ptr = 0;
}

N3Matrix::N3Matrix( int aR , int aC ){
  R = aR;
  C = aC;
  elem.allocate( R*C );
  ptr = 0;
}

void N3Matrix::allocate( int aR , int aC ){
  R = aR;
  C = aC;
  elem.allocate( R*C );
  ptr = 0;
}

int N3Matrix::c2s( int r , int c ){
  int snum;
  snum = r*C + c;
  if( snum >= R*C || r >= R || c >= C ){
    return -1;
  }
  return snum;
}

void N3Matrix::set( int snum , double val ){
  elem[snum] = val;
}

void N3Matrix::set( int r , int c , double val ){
  elem[c2s(r,c)] = val;
}

double N3Matrix::get( int snum ){
  return elem[snum];
}

double N3Matrix::get( int r , int c ){
  return elem[c2s(r,c)];
}

int N3Matrix::size( char t ){
  int ret = -1;
  if( t == size_r ){
    ret = R;
  }else if( t == size_c ){
    ret = C;
  }
  return ret;
}

void N3Matrix::zero( void ){
  int i;
  for( i = 0 ; i < R*C ; i++ ){
    elem[i] = 0;
  }
}

void N3Matrix::eye( void ){
  int i , j;
  for( i = 0 ; i < R ; i++ ){
    for( j = 0 ; j < C ; j++ ){
      elem[c2s(i,j)] = kroneckerDelta(i,j);
    }
  }
}

double N3Matrix::norm( void ){
  double ret = 0;
  int i , j;
  for( i = 0 ; i < R ; i++ ){
    for( j = 0 ; j < C ; j++ ){
      ret += this->get(i,j)*this->get(i,j);
    }
  }
  return sqrt(ret);
}

N3Matrix N3Matrix::normalize( void ){
  return *this/this->norm();
}

bool N3Matrix::iszero( double val ){
  bool ret;
  if( val < RANK_NORM_EPS && val > -RANK_NORM_EPS ){
    ret = true;
  }else{
    ret = false;
  }
  return ret;
}

N3Matrix N3Matrix::simp( N3Matrix *M ){
  int i , j , k , l;
  N3Matrix RET = *this , Mtmp(R,R) , Mswap(R,R);

  Mtmp.zero();
  Mswap.zero();
  M->eye();
  
  for( k = 0 ; k < R ; k++ ){
    for( i = 0 ; i < R ; i++ ){
      for( j = 0 ; j < C ; j++ ){
	Mtmp.set(i,j,kroneckerDelta(i,j));
      }
    }
    if( iszero( RET.get(k,k) ) && k < R-1 ){
      for( l = k+1 ; l < R ; l++ ){
	if( !iszero( RET.get(l,k) ) ){
	  Mswap.eye();
	  Mswap.set(k,k,0.0);
	  Mswap.set(k,l,1.0);
	  Mswap.set(l,l,0.0);
	  Mswap.set(l,k,1.0);
	  RET = Mswap*RET;
	  *M = Mswap*(*M);
	  break;
	}
      }
    }
    if( !iszero( RET.get(k,k) ) ){
      //for( i = k+1 ; i < R ; i++ ){
      for( i = k+1 ; i < R ; i++ ){
	Mtmp.set(i,k,-RET.get(i,k)/RET.get(k,k));
      }
      RET = Mtmp*RET;
      *M = Mtmp*(*M);
    }
  }
  return RET;
}

N3Matrix N3Matrix::simp( void ){
  N3Matrix M(R,C);
  return simp( &M );
}

N3Matrix N3Matrix::ge( N3Matrix& b ){
  int i , j;
  N3Matrix RET = b , M(R,R) , A , bb;
  A  = this->simp( &M );
  bb = M*b;
  for( i = R-1 ; i >= 0 ; i-- ){
    RET.set(i,0,bb.get(i,0)/A.get(i,i));
    for( j = R-1 ; j > i ; j-- ){
      RET.set(i,0,RET.get(i,0)-RET.get(j,0)*A.get(i,j)/A.get(i,i));
    }
  }
  return RET;
}

int N3Matrix::rank( void ){
  N3Matrix A;
  int i , ret = 0;
  A = this->simp();
  for( i = 0 ; i < R ; i++ ){
    if( A.getRowVec(i).norm() > RANK_NORM_EPS ){
      ret++;
    }
  }
  return ret;
}

N3Matrix N3Matrix::inv( void ){
  N3Matrix RET = *this , I(R,C) , b;
  N3Matrix  ctmp(R,1);
  int i;
  I.eye();
  for( i = 0 ; i < R ; i++ ){
    b = I.getColVec(i);
    ctmp = this->ge(b);
    RET.setColVec( i , ctmp );
  }
  return RET;
}

N3Matrix N3Matrix::pinv( void ){
  N3Matrix RET;
  int i , rnk;
  rnk = this->rank();
  if( rnk == R ){
    RET = (+(*this))*((*this)*(+(*this))).inv();
  }else{
    N3Matrix Simp , M(R,R) , SL(rnk,C) , N , C;
    SL.eye();
    Simp = this->simp(&M);
    N = SL*M;
    C = SL*Simp;
    RET = (+C)*((C*(+C)).inv())*N;
  }
  return RET;
}

N3Matrix N3Matrix::opmat( void ){
  N3Matrix pinvThis , I(C,C);
  pinvThis = this->pinv();
  I.eye();
  return I-pinvThis*(*this);
}

N3Matrix N3Matrix::gso( void ){
  N3Matrix RET(R,C) , x , y , v;
  int i , j;
  for( i = 0 ; i < C ; i++ ){
    x = this->getColVec(i);
    RET.setColVec(i,x);
    for( j = i-1 ; j >= 0 ; j-- ){
      y = RET.getColVec(j);
      v = ((+y)*x).get(0,0)*y/((+y)*y).get(0,0);
      RET.setColVec(i,RET.getColVec(i)-v);
    }
  }
  for( i = 0 ; i < C ; i++ ){
    RET.setColVec(i,RET.getColVec(i).normalize());
  }
  return RET;
}

void N3Matrix::QRdec( N3Matrix *Q , N3Matrix *R ){
  *Q = this->gso();
  *R = Q->inv()*(*this);
}

N3Matrix N3Matrix::getRowVec( int r ){
  N3Matrix RET(1,C);
  int i;
  for( i = 0 ; i < C ; i++ ){
    RET.set(i,this->get(r,i));
  }
  return RET;
}

N3Matrix N3Matrix::getColVec( int c ){
  N3Matrix RET(R,1);
  int i;
  for( i = 0 ; i < R ; i++ ){
    RET.set(i,this->get(i,c));
  }
  return RET;
}

void N3Matrix::setRowVec( int r , N3Matrix v ){
  int i;
  for( i = 0 ; i < C ; i++ ){
    this->set(r,i,v.get(0,i));
  }
}

void N3Matrix::setColVec( int c , N3Matrix v ){
  int i;
  for( i = 0 ; i < R ; i++ ){
    this->set(i,c,v.get(i,0));
  }
}

N3Matrix& N3Matrix::operator<<( double val ){
  this->elem[ptr] = val;
  ptr++;
  if( ptr == R*C ){
    ptr = 0;
  }
  return *this;
}

N3Matrix& N3Matrix::operator<<( int val ){
  *this << (double)val;
  return *this;
}

N3Matrix N3Matrix::operator+( N3Matrix M ){
  N3Matrix RET(R,C);
  int i;
  for( i = 0 ; i < R*C ; i++ ){
    RET.set(i,this->elem[i]+M.get(i));
  }
  return RET;
}

N3Matrix N3Matrix::operator-( N3Matrix M ){
  N3Matrix RET(R,C);
  int i;
  for( i = 0 ; i < R*C ; i++ ){
    RET.set(i,this->elem[i]-M.get(i));
  }
  return RET;
}

N3Matrix N3Matrix::operator*( N3Matrix M ){
  N3Matrix RET(R,M.size(size_c));
  int i , j , k;
  double tmp;
  for( i = 0 ; i < R ; i++ ){
    for( j = 0 ; j < M.size(size_c) ; j++ ){
      tmp = 0;
      for( k = 0 ; k < C ; k++ ){
	tmp += this->get(i,k)*M.get(k,j);
      }
      RET.set(i,j,tmp);
    }
  }
  return RET;
}

N3Matrix& N3Matrix::operator=( N3Matrix M ){
  int i;
  if( !elem.isallocated() ){
    this->R = M.size(size_r);
    this->C = M.size(size_c);
    this->elem.allocate(M.size(size_r)*M.size(size_c));
  }
  for( i = 0 ; i < M.size(size_r)*M.size(size_c) ; i++ ){
    this->set(i,M.get(i));
  }
  return *this;
}

ostream& operator<<( ostream& ost , N3Matrix M ){
  int i , j;
  for( i = 0 ; i < M.size(size_r) ; i++ ){
    for( j = 0 ; j < M.size(size_c) ; j++ ){
      ost << M.get(i,j) << "\t";
    }
    if( i < M.size(size_r)-1 ){
      ost << endl;
    }
  }
  return ost;
}

N3Matrix operator*( double a , N3Matrix M ){
  int i;
  N3Matrix RET(M.size(size_r),M.size(size_c));
  for( i = 0 ; i < M.size(size_r)*M.size(size_c) ; i++ ){
    RET.set(i,a*M.get(i));
  }
  return RET;
}

N3Matrix operator/( N3Matrix M , double a ){
  int i;
  N3Matrix RET(M.size(size_r),M.size(size_c));
  for( i = 0 ; i < M.size(size_r)*M.size(size_c) ; i++ ){
    RET.set(i,M.get(i)/a);
  }
  return RET;
}

N3Matrix operator+( N3Matrix M ){
  int i , j;
  N3Matrix RET(M.size(size_c),M.size(size_r));
  for( i = 0 ; i < M.size(size_r) ; i++ ){
    for( j = 0 ; j < M.size(size_c) ; j++ ){
      RET.set(j,i,M.get(i,j));
    }
  }
  return RET;
}

bool  operator==( N3Matrix& M1 , N3Matrix& M2 ){
  bool ret;
  if( M1.size(size_r) == M2.size(size_r) &&
      M1.size(size_c) == M2.size(size_c) ){
    ret = true;
  }else{
    ret = false;
  }
  return ret;
}

bool operator!=( N3Matrix& M1 , N3Matrix& M2 ){
  return !(M1==M2);
}
