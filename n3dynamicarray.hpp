#ifndef __N3_DYNAMIC_ARRAY_HPP__
#define __N3_DYNAMIC_ARRAY_HPP__

#include <cmath>
#include <iostream>

using namespace std;

template <class T>
class N3DynamicArray{
private:
  T *val;
  T  null;
  int size;
  bool is_allocated;
public:
  N3DynamicArray(){
    val          = NULL;
    is_allocated = false;
  }
  N3DynamicArray( int a_size ){
    val           = NULL;
    allocate( a_size );
  }
  N3DynamicArray( const N3DynamicArray<T>& obj ){
    int i;
    this->allocate( obj.getsize() );
    for( i = 0 ; i < obj.getsize() ; i++ ){
      this->val[i] = obj.copy(i);
    }
  }
  bool isallocated( void ){
    return is_allocated;
  }
  int  getsize( void ) const{
    return size;
  }
  void allocate( int a_size ){
    size = a_size;
    val = new T[size];
    is_allocated = true;
  }
  T copy( int i ) const{
    return val[i];
  }
  T& operator[] ( int i ){
    if( i >= 0 && i < size ){
      return val[i];
    }
    return null;
  }
  ~N3DynamicArray(){
    delete [] val;
  }
};

#endif
