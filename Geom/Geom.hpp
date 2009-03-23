#ifndef GEOM_HPP
#define GEOM_HPP

#include <Geom/Vec.hpp>
#include <Geom/Matrix.hpp>
#include <Geom/Util.hpp>

namespace Geom  {

typedef Vec<2,double> Vec2d;
typedef Vec<2,float > Vec2f;
typedef Vec<3,double> Vec3d;
typedef Vec<3,float > Vec3f;
typedef Vec<4,double> Vec4d;
typedef Vec<4,float > Vec4f;

typedef Vec<2,double> Vec2d;
typedef Vec<2,float > Vec2f;
typedef Vec<3,double> Vec3d;
typedef Vec<3,float > Vec3f;
typedef Vec<4,double> Vec4d;
typedef Vec<4,float > Vec4f;

typedef Matrix<2,2,double> Matrix22d;
typedef Matrix<2,3,double> Matrix23d;
typedef Matrix<2,4,double> Matrix24d;
typedef Matrix<3,2,double> Matrix32d;
typedef Matrix<3,3,double> Matrix33d;
typedef Matrix<3,4,double> Matrix34d;
typedef Matrix<4,2,double> Matrix42d;
typedef Matrix<4,3,double> Matrix43d;
typedef Matrix<4,4,double> Matrix44d;

typedef Matrix<2,2,float> Matrix22f;
typedef Matrix<2,3,float> Matrix23f;
typedef Matrix<2,4,float> Matrix24f;
typedef Matrix<3,2,float> Matrix32f;
typedef Matrix<3,3,float> Matrix33f;
typedef Matrix<3,4,float> Matrix34f;
typedef Matrix<4,2,float> Matrix42f;
typedef Matrix<4,3,float> Matrix43f;
typedef Matrix<4,4,float> Matrix44f;
}

#endif

