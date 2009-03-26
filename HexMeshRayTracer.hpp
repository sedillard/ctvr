#ifndef HEX_MESH_RAY_TRACER_HPP_INCLUDED
#define HEX_MESH_RAY_TRACER_HPP_INCLUDED

#include <iostream>

namespace HexMeshRayTracer__private {

struct Box 
{
  double min[3], max[3];
  Box( const double pMin[3], const double pMax[3] ) 
  {
    for (int i=0; i<3; ++i) min[i]=pMin[i],max[i]=pMax[i];
  }
  
  bool intersectsPoint(const double p[3]) 
  {
    for(uint d = 0 ; d < 3 ; d++ ) {
      if (p[d] < min[d] || p[d] > max[d]) return false;
    }
    return true;
  }
  
  bool intersectsRay ( const double p[3], const double v[3], 
                       double & tNear, double & tFar ) 
  {

    double t0 = tNear;
    double t1 = tFar;
    
    double tmin, tmax, tymin, tymax, tzmin, tzmax;
    
    if (v[0] >= 0) {
      tmin = (min[0] - p[0]) / v[0];
      tmax = (max[0] - p[0]) / v[0];
    } else {
      tmin = (max[0] - p[0]) / v[0];
      tmax = (min[0] - p[0]) / v[0];
    }

    if (v[1] >= 0) {
      tymin = (min[1] - p[1]) / v[1];
      tymax = (max[1] - p[1]) / v[1];
    }
    else {
      tymin = (max[1] - p[1]) / v[1];
      tymax = (min[1] - p[1]) / v[1];
    }

    if ( (tmin > tymax) || (tymin > tmax) )
      return false;
    if (tymin > tmin)
      tmin = tymin;
    if (tymax < tmax)
      tmax = tymax;
    if (v[2] >= 0) {
      tzmin = (min[2] - p[2]) / v[2];
      tzmax = (max[2] - p[2]) / v[2];
    }
    else {
      tzmin = (max[2] - p[2]) / v[2];
      tzmax = (min[2] - p[2]) / v[2];
    }
    if ( (tmin > tzmax) || (tzmin > tmax) )
      return false;
    if (tzmin > tmin)
      tmin = tzmin;
    if (tzmax < tmax)
      tmax = tzmax;
    
    if ( (tmin < t1) && (tmax > t0) ) {
      tNear = tmin;
      tFar = tmax;
      return true;
    } else {
      return false;
    }
  }
};

} //namespace





template <typename Value>
struct HexMeshRayTracer 
{
  Value *values;
  uint32_t size[3];

  double p[3]; //ray origin 
  double v[3]; //ray direction

  uint cell[3];
  double tNear, tFar; //current ray positions, at entry and exit points of the current cell
  
  int stepDim; //dimension we will step along when we move to next cell
  int stepDelta[3]; //direction (+/-) we will step in that dimension
  
  double F[8]; //function values at vertices of current cell
  uint32_t verts[8]; //indcies of vertices of current cell
  
  double segmentQueue[10]; //queue of ray parameter values, to be fetched with * operator
  int segmentIndex; //index into parameter queue, increment this with ++ operator
  int numSegmentsThisCell; //when pqi == this, go to next cell

  operator bool() { 
    return segmentIndex < numSegmentsThisCell; 
  }

  void operator ++() { 
    segmentIndex++; if (segmentIndex >= numSegmentsThisCell) step(); 
    } //advance to next ray segment
  
  double tFront() { return segmentQueue[segmentIndex] + 1e-10; }
  double tBack() { return segmentQueue[segmentIndex+1] - 1e-10; }



  HexMeshRayTracer
  ( Value *values_, 
    const uint32_t size_[3], 
    const double p_[3], 
    const double v_[3] 
  ) 
  : values(values_),tNear(-HUGE_VAL),tFar(HUGE_VAL)
  { 
    p[0] = p_[0]; p[1] = p_[1]; p[2] = p_[2];
    v[0] = v_[0]; v[1] = v_[1]; v[2] = v_[2];
    size[0]=size_[0], size[1]=size_[1], size[2]=size_[2];
  
    HexMeshRayTracer__private::Box 
      box( (double[]){0,0,0} , (double[]){size[0]-1, size[1]-1, size[2]-1} );

    if ( box.intersectsPoint( p ) ) {//check to see if we're inside the mesh
      cell[0] = (uint32_t)( p[0] );
      cell[1] = (uint32_t)( p[1] );
      cell[2] = (uint32_t)( p[2] );
      tNear = 0;
      for (uint32_t d = 0; d < 3; d++) {
        stepDelta[d] = ( v[d] > 0 ) ? 1 : -1 ;
      }
    } else {
      if (!box.intersectsRay(p,v,tNear,tFar)) {
        //no intersection. our work is done here
        //cell[0] = cell[1] = cell[2] = 0xffffffff;
        segmentIndex = numSegmentsThisCell = 0;
        return;
      }
      
      //find cell containing intersection point
      double cp[3] = { p[0]+tNear*v[0], p[1]+tNear*v[1], p[2]+tNear*v[2] };     
      cell[0] = (uint32_t)( cp[0] );
      cell[1] = (uint32_t)( cp[1] );
      cell[2] = (uint32_t)( cp[2] );
      
      for (uint d = 0; d < 3; d++) {
        cell[d] = (uint32_t)( cp[d] );
        cell[d] = (cell[d] > size[d] - 2) ? size[d] - 2 : cell[d]; 
          //make sure it really is in the bounds of the data
        stepDelta[d] = ( v[d] > 0 ) ? 1 : -1 ;
      }
    }
    
    calcVertices();
    calcCellExitPoint();
    segmentQueue[0] = tNear;
    segmentIndex = 0;
    numSegmentsThisCell = 1;
    calcNextParameterValues();
    segmentQueue[numSegmentsThisCell] = tFar;
  }

  void calcVertices() //get vertex indices as function values
  {
    verts[0] = cell[0]+size[0]*((cell[1])+size[1]*cell[2]);
    verts[1] = verts[0]+1;
    verts[2] = verts[0]+size[0];
    verts[3] = verts[0]+size[0]+1;
    verts[4] = verts[0]+size[0]*size[1];
    verts[5] = verts[0]+size[0]*size[1]+1;
    verts[6] = verts[0]+size[0]*size[1]+size[0];
    verts[7] = verts[0]+size[0]*size[1]+size[0]+1;
    for (uint i = 0; i < 8; i++) F[i] = values[verts[i]];
  }


  void step() //advance to next cell
  {
    //step to next cell
    cell[stepDim] += stepDelta[stepDim];
    
    if ( cell[0] >= size[0]-1 || cell[1] >= size[1]-1 || cell[2] >= size[2]-1 ) {
      segmentIndex = 0;
      numSegmentsThisCell = 0;
      return; 
    }
    
    //look up function values for next cell
    calcVertices();
    
    //entry point to this cell is known (it was the exit point of the last cell)
    tNear = tFar;
    
    //compute exit point of this cell
    calcCellExitPoint();
    
    //q up list of criticalities along the ray
    segmentQueue[0] = tNear;
    segmentIndex = 0;
    numSegmentsThisCell = 1;
    calcNextParameterValues();
    
    segmentQueue[numSegmentsThisCell] = tFar;
  }


  
  void calcCellExitPoint() //compute cell exit point for that cell
  {
    double tMin = HUGE_VAL; //ray parameter value which yields closest cell wall intersection
    int dMin = 0; //dimension xyz which yields the closest cell wall intersection
    for( int d  = 0; d < 3; d++) {
      if ( v[d] == 0 ) continue;
      double plane = ( v[d] > 0 ) ? cell[d]+1 : cell[d] ;
      double t = (plane - p[d]) / v[d];
      if (t < tMin ) {
          tMin = t;
          dMin = d;
      }
    }
    tFar = tMin;
    stepDim = dMin;
  }


  /*
  A = 	3 * a * v[0] * v[1] * v[2];

  B = 	2 * (
            a*((p[0]*v[1]*v[2])+(v[0]*p[1]*v[2])+(v[0]*v[1]*p[2])) +
            b*v[1]*v[2] +
            c*v[0]+v[2] +
            d*v[0]*v[1]
          );

  C = 	a*((p[0]+p[1]+v[2])*(p[0]+v[1]+p[2])*(v[0]+p[1]+p[2])) +
          b*(v[1]*p[2] + p[1]*v[2]) +
          c*(v[0]*p[2] + p[0]*v[2]) +
          d*(v[0]*p[1] + p[0]*v[1]) +
          e*v[0] + f*v[1] + g*v[2];

  For reference, here's where my abcdefgh come from:
  F(x,y,z) = axyz + byz + cxz + dxy + ex + fy + gz + h
  */


  void calcNextParameterValues() //F[] is function values at vertices of cell
  { 

    double a = F[1]+F[2]+F[4]+F[7]-F[0]-F[6]-F[5]-F[3],
           b = F[0]-F[1]-F[2]+F[3], 
           c = F[0]-F[1]-F[4]+F[5], 
           d = F[0]-F[2]-F[4]+F[6], 
           e = F[1]-F[0],
           f = F[2]-F[0],
           g = F[4]-F[0];
        
    double A = 3*a*v[0]*v[1]*v[2];
    double B = 2*( b*v[0]*v[1] + c*v[0]*v[2] + d*v[1]*v[2] + a*(p[0]*v[1]*v[2] + p[1]*v[0]*v[2] + p[2]*v[0]*v[1] ) );
    double C = a*(p[0]*p[1]*v[2] + p[0]*p[2]*v[1] + p[1]*p[2]*v[0] )+ e*v[0] + f*v[1] + g*v[2] + b*(p[0]*v[1] + p[1]*v[0] ) + c*(p[0]*v[2] + p[2]*v[0] ) + d*(p[1]*v[2] + p[2]*v[1] );
    
    if (A != 0) {
      double discr = B*B - 4*A*C;
      if (discr >= 0) {
        discr = sqrt(discr);
        
        double roots[2];
        roots[0] = ( -B + discr ) / (2*A);
        roots[1] = ( -B - discr ) / (2*A);
        
        if (roots[0] > roots[1]) {
          double tmp = roots[1];
          roots[1] = roots[0];
          roots[0] = tmp;
        }
        
        if (roots[0] == roots[1]) roots[1] = HUGE_VAL;
        
        for( uint r = 0; r < 2; r++) {
          if (roots[r] > tNear + 1e-6 && roots[r] < tFar - 1e-6 ) {
            segmentQueue[ numSegmentsThisCell++ ] = roots[r];
          } 
        }
      } 
    } 
  }
	
};



#endif

