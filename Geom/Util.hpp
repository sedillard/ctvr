#ifndef GEOM_UTIL_HPP_INCLUDED
#define GEOM_UTIL_HPP_INCLUDED

#include <utility> //for std::pair

namespace Geom {

template <typename float_t>
std::pair< Vec<3,float_t>, float_t >
trackball
  ( Vec<2,float_t> p0
  , Vec<2,float_t> p1
  )
{
  float_t d,z0,z1;
  d = 1-(dot(p0,p0));
  z0 = d < 0 ? 0 : sqrt(d);
  d = 1-(dot(p1,p1));
  z1 = d < 0 ? 0 : sqrt(d);
  Vec<3,float_t>
    v0 = vec3(p0[0],p0[1],z0),
    v1 = vec3(p1[0],p1[1],z1);

  return std::make_pair(
    normalize(cross(v0,v1)),  //axis
    acos(std::max(float_t(-1.0),std::min(float_t(1.0),dot(v0,v1)))) //angle
  );
}


template <int n, typename T>
T dist_to_line( Vec<n,T> const& p0, Vec<n,T> const& p1, Vec<n,T> const& q)
{
  Vec<n,T> l = p1-p0;
  l /= norm(l);
  return norm ( q - (dot(l,(q-p0))*l+p0) ) ;
}



} //namespace geom

#endif

