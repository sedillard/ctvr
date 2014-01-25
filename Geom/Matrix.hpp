#ifndef GEOM_MATRIX_HPP_INCLUDED
#define GEOM_MATRIX_HPP_INCLUDED

#include <algorithm> //for std::swap
#include <utility> //for std::pair
#include <iostream>
#include <cmath>

#include <Geom/Vec.hpp>

namespace Geom {

template <int m, int n, typename T>
struct Matrix : Vec< m, Vec<n,T> >
{
  typedef Vec< m, Vec<n,T> > base_type;
  typedef Vec<n,T> row_type;
  typedef Vec<m,T> col_type;

  Matrix () {}

  template <typename S>
  Matrix (S const* s) { for (int i=0; i<m; ++i) base_type::elems[i] = row_type(s+(i*n)); }

  template <typename S>
  Matrix (S (*s)[n]) { for (int i=0; i<m; ++i) base_type::elems[i] = row_type(s[i]); }

  Matrix( const T& d )
  {
    for (int i=0; i<m; ++i)
    for (int j=0; j<n; ++j) base_type::elems[i][j] = (i==j)?d:0;
  }

  static Matrix identity;

  const row_type& row(int i) const { return base_type::elems[i]; }
  col_type col(int j) const { return strided_vec<m,T,n>(&base_type::elems[0][j]); }

  void set_row( int i, const row_type& r ) { base_type::elems[i] = r; }
  void set_col( int j, const col_type& c ) { for (int i=0; i<m; ++i) base_type::elems[i][j] = c[i]; }

  void transpose()
  {
    for (int i=0; i<m; ++i)
    for (int j=(i+1); j<n; ++j)
      std::swap( base_type::elems[i][j], base_type::elems[j][i] );
  }

  void scale( const T& s ) { for (int i=0; i<m; ++i) base_type::elems[i][i] *= s; }
  void scale( const Vec<m,T>& s ) { for (int i=0; i<m; ++i) base_type::elems[i][i] *= s[i]; }


  template <int x, int y> static
  Matrix from_upper_left( Matrix<x,y,T> const& src ) {
    Matrix<m,n,T> dst;
    for (int i=0; i<m; ++i)
      for (int j=0; j<m; ++j)
        dst[i][j] = src[i][j];
    return dst;
  }
};

template <int m, int n, typename T>
Matrix<m,n,T> Matrix<m,n,T>::identity(1);

template <int m, int n, typename T> inline
Matrix<n,m,T> transpose(const Matrix<m,n,T>& a)
{
  Matrix<n,m,T> b;
  for (int i=0; i<m; ++i)
    for (int j=0; j<n; ++j)
      b[j][i] = a[i][j];
  return b;
}

template <int m, int n, int o, typename T> inline
Matrix<m,o,T> operator*( const Matrix<m,n,T>& a, const Matrix<n,o,T>& b )
{
  Matrix<m,o,T> c;
  for (int i=0; i<m; ++i)
  for (int j=0; j<o; ++j) {
    c[i][j] = dot( a.row(i), b.col(j) );
  }
  return c;
}

//left multiply: v is a column vector
template <int m, int n, typename T> inline
Vec<m,T> operator*( const Matrix<m,n,T>& a, const Vec<n,T>& v )
{
  Vec<m,T> c;
  for (int i=0; i<m; ++i) c[i] = dot(a.row(i),v);
  return c;
}

//right multiply: v is a row vector
template <int m, int n, typename T> inline
Vec<n,T> operator*( const Vec<n,T>& v,  const Matrix<m,n,T>& a )
{
  Vec<n,T> c;
  for (int j=0; j<n; ++j) c[j] = dot(a.col(j),v);
  return c;
}

//accumulate multiply
template <int n, typename T> inline
Matrix<n,n,T>& operator*=( Matrix<n,n,T>& a, const Matrix<n,n,T>& b )
{ return (a = a*b); }

template <int m, int n, typename T>
void gauss_elim( Matrix<m,n,T> & a )
{
  //U of LU decomp is put in upper part of a. L is left in lower part of a, but
  //entries of L still need to be divided by the diagonal of U
  for(int s=0; s<m-1;) {
    double pivot = fabs(a[s][s]);
    int r=s;
    for(int i=s+1; i<m; ++i) {
      double x = fabs(a[i][s]);
      if( pivot < x) {
        pivot=x;
        r=i;
      }
    }
    if(r != s)
      for(int j=0; j<n; ++j)
        std::swap(a[s][j],a[r][j]);
    for(int i=s+1; i<m; ++i) {
      double f = -a[i][s]/a[s][s];
      for(int j=s+1; j<n; ++j)
        a[i][j] += a[s][j]*f;
    }
    ++s;
  }
}

template <int m, int n, typename T> inline
void lu_decomp( Matrix<m,n,T> & a )
{
  gauss_elim(a);
  //finishes LU decomposition by dividing elements below the diagonal
  for (int j=0;   j<m; ++j)
  for (int i=j+1; i<m; ++i)
    a[i][j] /= a[j][j];
}

template <int n, typename T>
Matrix<n,n,T> inverse( const Matrix<n,n,T>& a)
{
  Matrix<n,n,T> r; //result
  Matrix<n,n+n,T> e; //extended Matrix
  for (int i=0; i<n; ++i)
  for (int j=0; j<n; ++j) {
    e[i][j] = a[i][j];
    e[i][j+n] = (i==j) ? 1 : 0;
  }
  gauss_elim(e); //put e in upper triangular
  //back substitution
  for(int i=n-1; i>=0; --i)
    for(int j=0; j<n; ++j) {
      for(int k=i+1; k<n; ++k)
        e[i][n+j] -= e[i][k]*e[k][n+j];
      e[i][n+j] /= e[i][i];
      r[i][j] = e[i][j+n]; //store result
  }
  return r;
}



template <int n, typename T>
T det( const Matrix<n,n,T>& a_)
{
  Matrix<n,n,T> a = a_;
  gauss_elim(a);
  T d = a[0][0];
  for (int i=1; i<n; ++i) d *= a[i][i];
  return d;
}

template <int n, typename T>
std::pair<Matrix<n,n,T>,T> inverse_and_det( const Matrix<n,n,T>& a)
{
  std::pair<Matrix<n,n,T>,T> r; //result
  Matrix<n,n+n,T> e; //extended matrix
  for (int i=0; i<n; ++i)
  for (int j=0; j<n; ++j) {
    e[i][j] = a[i][j];
    e[i][j+n] = (i==j) ? 1 : 0;
  }
  gauss_elim(e); //put e in upper triangular

  //get determinant
  r.second = e[0][0];
  for (int i=1; i<n; ++i) r.second *= e[i][i];

  //back substitution
  for(int i=n-1; i>=0; --i)
    for(int j=0; j<n; ++j) {
      for(int k=i+1; k<n; ++k)
        e[i][n+j] -= e[i][k]*e[k][n+j];
      e[i][n+j] /= e[i][i];
      r.first[i][j] = e[i][j+n]; //copy back into a
  }

  return r;
}


} //namespace vecmath

#endif

