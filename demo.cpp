#include <numerical3/n3matrix.hpp>
#include <numerical3/n3complexnumber.hpp>

int main( void ){
  int i , j;
  N3Matrix M(3,3) , A(3,3) , C , D(3,3) , E(2,3) , b(3,1) , MM(3,3) , Ms , Msd(2,3) , bb , bbd(2,1) , SL(2,3) , I3(3,3);
  M << 1 << 2 << 3
    << 1 << 2 << 3
    << 2 << 4 << 6;
  A << 2 << 3 << 4
    << 5 << 6 << 7
    << 8 << 9 << 10;
  D << 1 << 0 << 0
    << 0 << 2 << 1
    << 3 << 3 << 0;
  b << 1
    << 2
    << 3;
  I3.eye();
  pmat(M);
  pmat(A);
  pmat(M+M);
  pmat(M-M);
  pmat(A-M);
  pmat(A*M);
  pmat(2.0*M);
  
  C = A*M;

  pmat(C);
  pmat(+C);

  pmat(M.simp());
  pmat(M.rank());
  pmat(A.simp());
  pmat(A.rank());
  pmat(D);
  pmat(b);
  pmat(D.rank());
  pmat(D.ge(b));
  D << 1 << 0 << 0
    << 0 << 2 << 1
    << 5 << 3 << 10;
  pmat(D);
  pmat(D.ge(b));
  pmat(D.inv());

  E << 1 << 2 << 3
    << 3 << 2 << 1;
  pmat(E.rank());
  pmat(E.pinv());
  pmat(M.pinv());
  pmat(M.pinv()*b);
  pmat(M.pinv()*M*M.pinv());
  pmat(M*M.pinv()*M);
  pmat(M.opmat());
  pmat(+(M*M.pinv()));
  pmat(M*M.pinv());

  M.simp(&MM);
  pmat(MM);
  pmat(MM.inv()*M.simp());
  pmat(+MM);

  pmat(D.gso());
  pmat(+D.gso());
  pmat((+D.gso())*D.gso());

  N3Matrix Q(3,3) , R(3,3);
  D.QRdec(&Q,&R);
  pmat(Q);
  pmat(R);
  pmat(Q.inv()*D);
  pmat(D);
  pmat(Q*R);
  pmat(+Q*Q);

  N3Complex q1 , q2;
  q1(0.0,1.0);
  q2(1,1);
  pmat(q1);
  pmat(q2);
  pmat(q1+q2);
  pmat(q1*q2);
  pmat(q1.angle());
}
