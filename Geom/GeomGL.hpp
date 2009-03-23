#ifndef GEOM_GL_HPP_INCLUDED
#define GEOM_GL_HPP_INCLUDED

#include <GL/gl.h>
#include <Geom/Matrix.hpp>

namespace Geom {


Matrix<4,4,double> get_gl_matrix (GLenum which)
{
  Matrix<4,4,double> m;
  glGetDoublev(which,reinterpret_cast<double*>(&m));
  return transpose(m);
}



};




#endif

