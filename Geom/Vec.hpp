#ifndef GEOM_VEC_HPP_INCLUDED
#define GEOM_VEC_HPP_INCLUDED

#include <iostream>
#include <cmath>

namespace Geom {

template < int n, typename T >
struct Vec 
{
  T elems[n];

  T& operator[](int i) { return elems[i]; }
  const T& operator[](int i) const { return elems[i]; }

  Vec () {}
  Vec (const Vec& o) { for (int i=0; i<n; ++i) elems[i]=o.elems[i]; }

  Vec (const T& x)   { for (int i=0; i<n; ++i) elems[i]=x; }

  template <typename S>
  Vec (const S* a) { for (int i=0; i<n; ++i) elems[i]=T(a[i]); }

  template <typename U> 
  explicit 
  Vec( Vec<n,U> const& u ) { for (int i=0;i<n;++i) elems[i]=T(u[i]); }
 
  Vec operator+(const Vec& b) const { Vec a; for (int i=0; i<n; ++i) {a.elems[i] = elems[i]+b.elems[i];} return a; }
  Vec operator-(const Vec& b) const { Vec a; for (int i=0; i<n; ++i) {a.elems[i] = elems[i]-b.elems[i];} return a; }

  template<typename S>
  Vec operator*(const S& s) const { Vec a; for (int i=0; i<n; ++i) {a.elems[i] = elems[i]*T(s);} return a; }

  template<typename S>
  Vec operator/(const S& s) const { Vec a; for (int i=0; i<n; ++i) {a.elems[i] = elems[i]/T(s);} return a; } 

  Vec& operator= (const Vec& o) { for (int i=0; i<n; ++i) { elems[i] =o.elems[i]; } return *this; }
  Vec& operator+=(const Vec& o) { for (int i=0; i<n; ++i) { elems[i]+=o.elems[i]; } return *this; }
  Vec& operator-=(const Vec& o) { for (int i=0; i<n; ++i) { elems[i]-=o.elems[i]; } return *this; }
  Vec& operator*=(const T& x)   { for (int i=0; i<n; ++i) { elems[i]*=x; } return *this; }
  Vec& operator/=(const T& x)   { for (int i=0; i<n; ++i) { elems[i]/=x; } return *this; }

  typedef Vec< n+1, T > dim_higher;
  typedef Vec< n-1, T > dim_lower;

 
    //lexicographical ordering
  bool operator< ( const Vec& o ) const {
    for (int i=0; i<n; ++i) {
        if (elems[i] < o.elems[i]) return true;
        if (elems[i] > o.elems[i]) return false;
    }
    return false;
  }

    //lexicographical ordering
  bool operator> ( const Vec& o ) const {
    for (int i=0; i<n; ++i) {
        if (elems[i] > o.elems[i]) return true;
        if (elems[i] < o.elems[i]) return false;
    }
    return false;
  }

  bool operator== ( const Vec& o ) const {
    for (int i=0; i<n; ++i) if (elems[i]!=o.elems[i]) return false;
    return true;
  }
  bool operator!= ( const Vec& o ) const {
    for (int i=0; i<n; ++i) if (elems[i]!=o.elems[i]) return true;
    return false;
  }

  
  Vec& normalize() {
    T l = norm(*this);
    (*this) /= l;
    return *this;
  }

};


template <int n, typename T, int s> inline
Vec<n,T> strided_vec( const T* p ) { Vec<n,T> v; for (int i=0; i<n; ++i) v.elems[i]=p[i*s]; return v; }

template <int n, typename T, typename S> inline 
Vec<n,T> operator*(const S& s, const Vec<n,T>& v) 
{
  Vec<n,T> a;
  for (int i=0; i<n; ++i) a[i] = v[i]*T(s);
  return a;
};

template <int n, typename T> inline 
Vec<n,T> operator-(const Vec<n,T>& v) 
{ 
  Vec<n,T> a;
  for (int i=0; i<n; ++i) a[i] = -v[i];
  return a;
};

template <int n, typename T> inline
T dot(const Vec<n,T>& a, const Vec<n,T>& b) 
{
  T s=0;
  for (int i=0; i<n; ++i) s+=a[i]*b[i];
  return s;
};

template<typename T> inline
Vec<2,T> vec2( const T& x, const T& y ) 
{ 
  Vec<2,T> v; v[0]=x,v[1]=y; return v;
}

template<typename T> inline
Vec<3,T> vec3( const T& x, const T& y, const T& z ) 
{ 
  Vec<3,T> v; v[0]=x,v[1]=y,v[2]=z; return v; 
}

template<typename T> inline
Vec<4,T> vec4( const T& x, const T& y, const T& z, const T& w ) 
{ 
  Vec<4,T> v; v[0]=x,v[1]=y,v[2]=z,v[3]=w; return v;
}

template <typename T> inline
Vec<3,T> cross(const Vec<3,T>& a, const Vec<3,T>& b) 
{
  return vec3( a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0] );
}

template <int n, typename T> inline
T sqr_norm( const Vec<n,T>& v ) { return dot(v,v); }

template <int n, typename T> inline
T norm( const Vec<n,T>& v ) { return sqrt(dot(v,v)); }

template <int n, typename T> inline
Vec<n,T> normalize( Vec<n,T> const& v ) { return v/norm(v); }



template <int n, typename T> inline
typename Vec<n,T>::dim_higher hvec ( Vec<n,T> const& v ) 
{
  typename Vec<n,T>::dim_higher h;
  for (int i=0; i<n; ++i) h[i] = v[i];
  h[n] = 0;
  return h;
}

template <int n, typename T> inline
typename Vec<n,T>::homog hpoint ( Vec<n,T> const& v ) 
{
  typename Vec<n,T>::homog h;
  for (int i=0; i<n; ++i) h[i] = v[i];
  h[n] = 1;
  return h;
}

template <int n, typename T> inline
typename Vec<n,T>::dim_lower project( Vec<n,T> const& h )
{
  typename Vec<n,T>::dim_lower v(h.elems);
  v /= h[n];
  return v;
}


} //namespace geom

namespace std {
template <int n, typename T> inline
std::ostream& operator<< (std::ostream& o, const Geom::Vec<n,T>& a) 
{
  o << "[" << a.elems[0];
  for (int i=1; i<n; ++i) o << ',' << a.elems[i];
  o << ']';
  return o;
}
}


#endif

