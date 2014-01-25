#ifndef TRILINEAR_HPP_INCLUDED
#define TRILINEAR_HPP_INCLUDED

#include <iostream>
#include <algorithm>
#include <cmath>

#include "Tourtre.hpp"


template<typename T>
void append_vector( std::vector<T> & dst, std::vector<T> & src )
{
  size_t s = dst.size();
  dst.resize(s+src.size());
  std::copy(src.begin(),src.end(),dst.begin()+s);
}

struct Saddle
{
  enum Type { YZFace, XZFace, XYFace, Body1, Body2 };
  enum Sort { SortNaturally, SortBefore, SortAfter };

  uint32_t vert; //the vertex number of this saddle;
  Saddle* next; //saddles are stored in a linked list in the voxel
  float value;
  uint32_t whom; //against whom must this saddle's sorting order be enforced?
  Tourtre::SweepComponent<uint32_t> *join_comp, *split_comp;

  //TODO: put these members into packed flags field
  int8_t type;
  int8_t symbolic;
  int8_t sort;
  int8_t orientation;
      // For pairs of body saddles:
      //   Indicates which cell vertex is reachable by a monotone path.
      // For face saddles:
      //  0 : vertices 0 and 3 are below
      //  1 : vertices 1 and 2 are below


  Saddle() :
    vert(-1),next(0),whom(-1),
    join_comp(0), split_comp(0),
    symbolic(false),
    sort(SortNaturally),
    orientation(-1)
    {}
};

template <typename Value>
struct Trilinear
{
  typedef Tourtre::SweepComponent<uint32_t> Component;
  uint32_t size[3]; //voxel grid size
  uint32_t stride[3]; //voxel grid strides
  uint32_t nvoxels; //product of sizes
  Value *values; //voxel values (the image)
  std::vector<Saddle*> saddle_map; //for each voxel, yields the linked list of saddles stored there
  std::vector<Component*> join_comps, split_comps;

  //reachable_max and reachable_min yield, for each voxel, a max or min that is
  //reachable by a monotone path
  std::vector<uint32_t> reachable_max,reachable_min;

  //maxes and mins (used for flooding out reachable regions)
  std::vector<uint32_t> maxes, mins;

  std::vector<Saddle*> nonvoxel_saddles;
  std::vector<uint32_t> voxel_saddles, place_holders;

  //About saddle indexes: A saddle is indexed by an integer. If there are n
  //voxels then an index i specifies the following..
  // 0  <= i < n  : a regular vertex
  // n  <= i < 2n : a yz face saddle
  // 2n <= i < 3n : an xz face saddle
  // 3n <= i < 4n : an xy face saddle
  // 4n <= i < 5n : a body 1-saddle
  // 5n <= i < 6n : a body 2-saddle
  //
  // NB: (i/n)-1 == enum Saddle::Type


  Trilinear( Value *values_, uint32_t ncols, uint32_t nrows, uint32_t nstacks )
  {
    values = values_; //does not own!
    size[0] = ncols;
    size[1] = nrows;
    size[2] = nstacks;
    nvoxels = nrows*ncols*nstacks;
    stride[0] = 1;
    stride[1] = ncols;
    stride[2] = ncols*nrows;
    saddle_map.resize(nvoxels,0);
    join_comps.resize(nvoxels,0);
    split_comps.resize(nvoxels,0);
  }

  ~Trilinear()
  {
    for ( uint32_t i=0; i<nvoxels; ++i ) {
      Saddle *s=saddle_map[i];
      while (s) {
        Saddle *t =s;
        s=s->next;
        delete t;
      }
    }
  }

  Saddle* saddle( uint32_t i ) const
    //here the i/nvoxels determines the type
  {
    Saddle::Type t = Saddle::Type( (i/nvoxels)-1 );
    Saddle *s = saddle_map[i%nvoxels];
    while ( s && s->type!=t ) s=s->next;
    return s;
  }

  Saddle* saddle( Saddle::Type t, uint32_t i ) const
    //here i < nvoxels
  {
    Saddle *s = saddle_map[i];
    while ( s && s->type!=t ) s=s->next;
    return s;
  }

  Saddle* saddle( Saddle::Type t, uint32_t x, uint32_t y, uint32_t z ) const
    //here i < nvoxels
  {
    return saddle(t, x+y*stride[1]+z*stride[2] );
  }

  void voxel_pos( uint32_t i, uint32_t & x, uint32_t & y, uint32_t & z ) const
  {
    assert( i < nvoxels );
    z = i / stride[2],
    y = (i-z*stride[2]) / size[0],
    x = i % size[0];
  }


  //Comparison functions
  static
  bool compare_ptrs( const Value *a, const Value *b )
  {
    return *a!=*b ? *a<*b : a<b;
  }

  struct compare_voxels
  {
    const Trilinear & t;
    compare_voxels( const Trilinear & t_ ) : t(t_) {}
    bool operator()( const uint32_t & a, const uint32_t b ) const
    {
      assert(a<t.nvoxels&&b<t.nvoxels);
      return compare_ptrs(t.values+a,t.values+b);
    }
  };

  struct compare_saddles
  {
    const Trilinear & t;
    compare_saddles( const Trilinear & t_ ) : t(t_) {}
    bool operator()( Saddle* sa, Saddle* sb ) const
    {
      return (sa->value!=sb->value)
                ? (sa->value < sb->value)
                : (sa->whom != sb->whom)
                    ? (sa->whom < sb->whom)
                    : (sa->sort < sb->sort);
    }
    bool operator()( const uint32_t & a, const uint32_t & b ) const
    {
      Saddle *sa = t.saddle(a);
      Saddle *sb = t.saddle(b);
      return (*this)(sa,sb);
    }
  };

  int voxel_neighbors(uint32_t x,uint32_t y,uint32_t z,uint32_t nbrs[]) const
  {
    int nnbrs=0;
    if ( x > 0 )         nbrs[nnbrs++] = (x-1)+y*stride[1]+z*stride[2];
    if ( y > 0 )         nbrs[nnbrs++] = x+(y-1)*stride[1]+z*stride[2];
    if ( z > 0 )         nbrs[nnbrs++] = x+y*stride[1]+(z-1)*stride[2];
    if ( x < size[0]-1 ) nbrs[nnbrs++] = (x+1)+y*stride[1]+z*stride[2];
    if ( y < size[1]-1 ) nbrs[nnbrs++] = x+(y+1)*stride[1]+z*stride[2];
    if ( z < size[2]-1 ) nbrs[nnbrs++] = x+y*stride[1]+(z+1)*stride[2];
    return nnbrs;
  }

  int voxel_neighbors( uint32_t i, uint32_t nbrs[] ) const
  {
    uint32_t x,y,z;
    voxel_pos(i,x,y,z);
    return voxel_neighbors(x,y,z,nbrs);
  }

  void push_if_less( Value *i, Value *j, int & nlink, uint32_t link[] ) const
  {
    if ( compare_ptrs(j,i) ) link[nlink++] = j-values;
  }

  int voxel_lower_link( uint32_t i, uint32_t link[] ) const
  {
    uint32_t x,y,z;
    voxel_pos(i,x,y,z);
    int nlink=0;
    Value *p = values+i;
    if (x>0)         push_if_less(p,p-stride[0],nlink,link);
    if (y>0)         push_if_less(p,p-stride[1],nlink,link);
    if (z>0)         push_if_less(p,p-stride[2],nlink,link);
    if (x<size[0]-1) push_if_less(p,p+stride[0],nlink,link);
    if (y<size[1]-1) push_if_less(p,p+stride[1],nlink,link);
    if (z<size[2]-1) push_if_less(p,p+stride[2],nlink,link);
    return nlink;
  }

  void push_if_greater( Value *i, Value *j, int & nlink, uint32_t link[] ) const
  {
    if ( compare_ptrs(i,j) ) link[nlink++] = j-values;
  }

  int voxel_upper_link( uint32_t i, uint32_t link[] ) const
  {
    uint32_t x,y,z;
    voxel_pos(i,x,y,z);
    int nlink=0;
    Value *p = values+i;
    if (x>0)         push_if_greater(p,p-stride[0],nlink,link);
    if (y>0)         push_if_greater(p,p-stride[1],nlink,link);
    if (z>0)         push_if_greater(p,p-stride[2],nlink,link);
    if (x<size[0]-1) push_if_greater(p,p+stride[0],nlink,link);
    if (y<size[1]-1) push_if_greater(p,p+stride[1],nlink,link);
    if (z<size[2]-1) push_if_greater(p,p+stride[2],nlink,link);
    return nlink;
  }



  int saddle_lower_link( Saddle::Type type, uint32_t i, uint32_t nbrs[] ) const
  {
    Saddle *s = saddle(type,i);
    switch(type) {
      case Saddle::YZFace:
      case Saddle::XZFace:
      case Saddle::XYFace: {
        uint32_t s1=stride[(type+1)%3], s2=stride[(type+2)%3];
        if (s->orientation == 0) {
          nbrs[0]=i, nbrs[1]=i+s1+s2;
        } else {
          nbrs[0]=i+s1, nbrs[1]=i+s2;
        }
        return 2;
      }
      case Saddle::Body1:
      case Saddle::Body2: {
        int not_connected = (s->orientation>=0) ? (~s->orientation&7) : -1 ;
        if (s->symbolic) {
          //perturbation saddle
          if (s->orientation==0) { //1-saddle
            nbrs[0] = i;
            nbrs[1] = i+stride[0]+stride[1];
            return 2;
          } else { //2-saddle
            nbrs[0] = i+stride[0]+stride[1];
            return 1;
          }
        } else {
          int nnbrs=0;
          for ( int j=0; j<8; ++j ) {
            if (j != not_connected) {
              uint32_t n = i + (j&1)*stride[0] +
                               ((j&2)>>1)*stride[1] +
                               ((j&4)>>2)*stride[2] ;
              if ( values[n] < s->value ) nbrs[nnbrs++] = n;
            }
          }
          return nnbrs;
        }
      }
    }
    assert(0&&"wtf?");
    return 0;
  }

  int saddle_upper_link( Saddle::Type type, uint32_t i, uint32_t nbrs[] ) const
  {
    Saddle *s = saddle(type,i);
    switch(type) {
      case Saddle::YZFace:
      case Saddle::XZFace:
      case Saddle::XYFace: {
        uint32_t s1=stride[(type+1)%3], s2=stride[(type+2)%3];
        if (s->orientation == 0) {
          nbrs[0]=i+s1, nbrs[1]=i+s2;
        } else {
          nbrs[0]=i, nbrs[1]=i+s1+s2;
        }
        return 2;
      }
      case Saddle::Body1:
      case Saddle::Body2: {
        int not_connected = (s->orientation>=0) ? (~s->orientation&7) : -1 ;
        if (s->symbolic) {
          //perturbation saddle
          if (s->orientation==0) { //1-saddle
            nbrs[0] = i+stride[0];
            return 1;
          } else { //2-saddle
            nbrs[0] = i+stride[0];
            nbrs[1] = i+stride[0]+stride[1]+stride[2];
            return 2;
          }
        } else {
          int nnbrs=0;
          for ( int j=0; j<8; ++j ) {
            if (j != not_connected) {
              uint32_t n = i + (j&1)*stride[0] +
                               ((j&2)>>1)*stride[1] +
                               ((j&4)>>2)*stride[2] ;
              if ( values[n] > s->value ) nbrs[nnbrs++] = n;
            }
          }
          return nnbrs;
        }
      }
    }
    assert(0&&"wtf?");
    return 0;
  }

  int saddle_lower_link( uint32_t i, uint32_t link[] ) const
  {
    return saddle_lower_link(Saddle::Type((i/nvoxels)-1),i%nvoxels,link);
  }

  int saddle_upper_link( uint32_t i, uint32_t link[] ) const
  {
    return saddle_upper_link(Saddle::Type((i/nvoxels)-1),i%nvoxels,link);
  }

  int lower_link( uint32_t i, uint32_t link[] ) const
  {
    if (i < nvoxels) {
      return voxel_lower_link(i,link);
    } else {
      return saddle_lower_link(i,link);
    }
  }

  int upper_link( uint32_t i, uint32_t link[] ) const
  {
    if (i < nvoxels) {
      return voxel_upper_link(i,link);
    } else {
      return saddle_upper_link(i,link);
    }
  }



  //returns a newly allocated saddle if there's one here, null otherwise
  Saddle* is_face_saddle ( Value* f[4] )  //vertices on the face here
  {
    //first check to see if theres a bilinear saddle in f
    bool e0 = compare_ptrs(f[0],f[1]) ,
         e1 = compare_ptrs(f[0],f[2]) ,
         e2 = compare_ptrs(f[3],f[1]) ,
         e3 = compare_ptrs(f[3],f[2]) ;
    if ( e0 == e1 && e1 == e2 && e2 == e3 ) {
      Saddle *s = new Saddle;
      s->orientation = e0 ? 0 : 1;

      //ok there's a candidate saddle in the face
      //check for real (numerical) saddle
      double a = (double)*f[0] - *f[1] - *f[2] + *f[3];
      if (a != 0) {

        if ( *f[0]==*f[1] ) { //saddle on bottom edge
          s->sort = Saddle::SortAfter;
          s->whom = f[0]-values;
          s->value = *f[0];
        } else if ( *f[0] == *f[2] ) { //saddle on left edge
          s->sort = Saddle::SortAfter;
          s->whom = f[0]-values;
          s->value = *f[0];
        } else if ( *f[1] == *f[3] ) { //saddle on right edge
          s->sort = Saddle::SortBefore;
          s->whom = f[3]-values;
          s->value = *f[3];
        } else if ( *f[2] == *f[3] ) { //saddle on top edge
          s->sort = Saddle::SortBefore;
          s->whom = f[3]-values;
          s->value = *f[3];
        } else { //interior saddle
          double
            b = double(*f[1]) - *f[0],
            c = double(*f[2]) - *f[0],
            x = -c/a,
            y = -b/a,
            v = a*x*y + b*x + c*y + *f[0];
          assert( x>0 && x<1 && y>0 && y<1 );
          assert(v != *f[0]);
          assert(v != *f[1]);
          assert(v != *f[2]);
          assert(v != *f[3]);
          s->value = v;
        }
      } else {
        if ( e0 ) {
          //   1 <---- 0
          //   ^       |  Saddle needs to sort
          //   |       |  after 0 and (1)
          //   | s     v
          //  (1)----> 1
          s->symbolic = true;
          s->sort = Saddle::SortAfter;
          s->value = *f[0];
          s->whom = f[0]-values;
        } else {
          //  0 ---->(0)
          //  |     s ^  Saddle needs to sort
          //  |       |  before 1 and (0)
          //  v       |
          //  1 <---- 0
          s->symbolic = true;
          s->sort = Saddle::SortBefore;
          s->value = *f[3];
          s->whom = f[3]-values;
        }
      }
      return s;
    } else {
      return 0;
    }
  }

  void find_saddles()
  {
    find_face_saddles();
    find_body_saddles();
  }

  void find_face_saddles()
  {
    std::cout << "face saddles" << std::endl;
    int32_t size[3] = { this->size[0], this->size[1], this->size[2] };
    uint32_t offs[3] = { stride[0], stride[1], stride[2] };
    for ( int dim=0; dim<3; ++dim ) {

      #pragma omp parallel for
      for ( int32_t k=0; k<size[2]-1; ++k ) {
        std::vector<uint32_t> place_holders;
        std::vector<Saddle*> saddles;
        for ( int32_t j=0; j<size[1]-1; ++j ) {
          Saddle::Type type = Saddle::Type(dim);

          Value* v[4] = { values + ( j )*offs[1] + ( k )*offs[2]
                        , values + (j+1)*offs[1] + ( k )*offs[2]
                        , values + ( j )*offs[1] + (k+1)*offs[2]
                        , values + (j+1)*offs[1] + (k+1)*offs[2] };

          uint32_t jk = j*offs[1] + k*offs[2];
          for ( int32_t i=0; i<size[0]; ++i ) {
            Saddle *s=0;
            if ( (s=is_face_saddle(v)) ) {
              uint32_t ind = i*offs[0]+jk;
              s->next = saddle_map[ind];
              saddle_map[ind] = s;
              s->vert = ind+(1+type)*nvoxels;
              s->type = type;
              saddles.push_back(s);
              if (s->whom != uint32_t(-1)) place_holders.push_back(s->whom);
              assert(saddle(s->vert) == s);
            }

            //move to next face
            for ( int i=0; i<4; ++i )  v[i] += offs[0];
          }
        }
        #pragma omp critical
        {
          append_vector(nonvoxel_saddles,saddles);
          append_vector(this->place_holders,place_holders);
        }
      }
      //rotate to next dimension
      { uint32_t t = size[0];
        size[0] = size[1];
        size[1] = size[2];
        size[2] = t;
        t = offs[0];
        offs[0] = offs[1];
        offs[1] = offs[2];
        offs[2] = t;
      }
    }
  }



  std::pair<Saddle*,Saddle*> are_body_saddles( Value v[8] )
  {
    double
      h = v[0],
      e = v[1] - v[0],
      f = v[2] - v[0],
      g = v[4] - v[0],
      b = v[0] - v[1] - v[2] + v[3],
      c = v[0] - v[2] - v[4] + v[6],
      d = v[0] - v[4] - v[1] + v[5],
      a = v[1] + v[2] + v[4] + v[7] - v[0] - v[6] - v[5] - v[3],
      ax = a*e - b*d,
      ay = a*f - b*c,
      az = a*g - c*d;

    Saddle *s1=0, *s2=0;
    if (a != 0) {
      if ( ax*ay*az < 0 ) {
        double discr = sqrt( -ax*ay*az );

        double x1 = -( c + discr / ax ) / a;
        double y1 = -( d + discr / ay ) / a;
        double z1 = -( b + discr / az ) / a;

        if (0<x1 && x1<1 && 0<y1 && y1<1 && 0<z1 && z1<1 ) {
          s1 = new Saddle;
          s1->value =
            a*x1*y1*z1 + b*x1*y1 + c*y1*z1 + d*x1*z1 + e*x1 + f*y1 + g*z1 + h;
          s1->type = Saddle::Body1;
        }

        double x2 =  -( c - discr / ax ) / a;
        double y2 =  -( d - discr / ay ) / a;
        double z2 =  -( b - discr / az ) / a;

        if (0<x2 && x2<1 && 0<y2 && y2<1 && 0<z2 && z2<1 ) {
          s2 = new Saddle;
          s2->value =
            a*x2*y2*z2 + b*x2*y2 + c*y2*z2 + d*x2*z2 + e*x2 + f*y2 + g*z2 + h;
          s2->type = Saddle::Body2;
        }

        if (s1 && s2) { //set the saddle orientations
          int o = (x1 < x2 ? 0 : 1) | (y1 < y2 ? 0 : 2) | (z1 < z2 ? 0 : 4);
          s1->orientation = o;
          s2->orientation = ~o&7;
        }
      }
    } else {
      if ( ax*ay*az != 0 ) {
        double x = (c*e - b*g - d*f ) / (2*b*d);
        double y = (d*f - b*g - c*e ) / (2*b*c);
        double z = (b*g - c*e - d*f ) / (2*c*d);
        if (0<x && x<1 && 0<y && y<1 && 0<z && z<1 ) {
          s1 = new Saddle;
          s1->type = Saddle::Body1;
          s1->value = b*x*y + c*y*z + d*x*z + e*x + f*y + g*z + h;
        }
      }
    }
    return std::make_pair(s1,s2);
  }

  void find_body_saddles ()
  {
    std::cout << "body saddles" << std::endl;
    #pragma omp parallel for
    for (int z = 0; z < (int)size[2]-1; ++z) {
      std::vector<Saddle*> saddles;
      for (uint32_t y = 0; y < size[1]-1; ++y) {
        uint32_t i=y*stride[1]+z*stride[2];
        for (uint32_t x = 0; x < size[0]-1; ++x,++i) {
          Value v[8] =
            { values[i]
            , values[i+1]
            , values[i+stride[1]]
            , values[i+1+stride[1]]
            , values[i+stride[2]]
            , values[i+1+stride[2]]
            , values[i+stride[1]+stride[2]]
            , values[i+1+stride[1]+stride[2]] };
          std::pair<Saddle*,Saddle*> s = are_body_saddles(v);
          if (s.first) {
            s.first->next = saddle_map[i];
            saddle_map[i] = s.first;
            s.first->vert = i+4*nvoxels;
            saddles.push_back( s.first );
            assert(saddle(s.first->vert) == s.first);
          }
          if (s.second) {
            s.second->next = saddle_map[i];
            saddle_map[i] = s.second;
            s.second->vert = i+5*nvoxels;
            saddles.push_back( s.second );
            assert(saddle(s.second->vert) == s.second);
          }
        }
      }
      #pragma omp critical
      {
        append_vector(nonvoxel_saddles,saddles);
      }
    }

    #pragma omp for
    for (int z = 0; z < (int)size[2]-1; ++z) {
      std::vector<uint32_t> place_holders;
      std::vector<Saddle*> saddles;
      for (uint32_t y = 0; y<size[1]-1; ++y) {
        uint32_t i=y*stride[1]+z*stride[2];
        for (uint32_t x = 0; x<size[0]-1; ++x,++i) {
          Saddle *s[3] = { saddle(Saddle::YZFace,i) ,
                           saddle(Saddle::XZFace,i) ,
                           saddle(Saddle::XYFace,i) };
          if ( s[0] && s[1] && s[2] &&
              s[0]->orientation == 0 &&
              s[1]->orientation == 0 &&
              s[2]->orientation == 0 &&
              s[0]->symbolic &&
              s[1]->symbolic &&
              s[2]->symbolic )
          {
            Saddle *bs = saddle(Saddle::Body1,i);
            if (!bs) bs = saddle(Saddle::Body2,i);
            if (!bs) {
              /*        0--------0
              *        /|       /|
              *       / |      / |
              *      1--------0  |
              *      ^  1-----|--0
              *      | /      | / Body saddle s needs
              *      |/'s     |/  to sort above (1) but
              *     (1)-->--- 1   below the other 1s
              */
              Saddle *s = new Saddle;
              s->type = Saddle::Body1;
              s->value = values[i];
              s->orientation = 0;
              s->symbolic = true;
              s->sort = Saddle::SortAfter;
              s->whom = i;
              s->next = saddle_map[i] ;
              saddle_map[i] = s;
              s->vert = i+4*nvoxels;
              saddles.push_back( s );
              if (s->whom != uint32_t(-1)) place_holders.push_back(s->whom);
              assert(saddle(s->vert)==s);
            }
          }
        }
      }
      #pragma omp critical
      {
        append_vector(nonvoxel_saddles,saddles);
        append_vector(this->place_holders,place_holders);
      }
    }
    #pragma omp parallel for
    for (int z = 1; z < (int)size[2]; ++z) {
      std::vector<uint32_t> place_holders;
      std::vector<Saddle*> saddles;
      for (size_t y = 1; y < size[1]; ++y) {
        uint32_t i=y*stride[1]+z*stride[2];
        for (size_t x = 1; x < size[0]; ++x,++i) {
          Saddle *s[3] = { saddle(Saddle::YZFace, x ,y-1,z-1) ,
                           saddle(Saddle::XZFace,x-1, y ,z-1) ,
                           saddle(Saddle::XYFace,x-1,y-1, z ) };
          if ( s[0] && s[1] && s[2] &&
               s[0]->orientation == 1 &&
               s[1]->orientation == 1 &&
               s[2]->orientation == 1 &&
               s[0]->symbolic &&
               s[1]->symbolic &&
               s[2]->symbolic )
          {
            Saddle *bs = saddle(Saddle::Body1,x-1,y-1,z-1);
            if (!bs) bs = saddle(Saddle::Body2,x-1,y-1,z-1);
            if (!bs)  {
              /*       0---->--(0)
              *      /|    s ,/|
              *     / |      / ^
              *    1--------0  |
              *    |  1-----|--0
              *    | /      | / Body saddle s needs
              *    |/       |/  to sort below (0) but
              *    1--------1   above the other 0's
              */
              Saddle *s = new Saddle;
              s->type = Saddle::Body2;
              s->value = values[i];
              s->orientation = 8+7;
              s->symbolic = true;
              s->sort = Saddle::SortBefore;
              s->whom = i+stride[0]+stride[1]+stride[2];
              s->next = saddle_map[i] ;
              saddle_map[i] = s;
              s->vert = i+5*nvoxels;
              saddles.push_back(s);
              if (s->whom != uint32_t(-1)) place_holders.push_back(s->whom);
              assert(saddle(s->vert)==s);
            }
          }
        }
      }
      #pragma omp critical
      {
        append_vector(nonvoxel_saddles,saddles);
        append_vector(this->place_holders,place_holders);
      }
    }
  }

  void find_critical_voxels()
  {
    std::cout << "find critical voxels" << std::endl;
    #pragma omp parallel for
    for ( int z=0; z<int(size[2]); ++z ) {
      std::vector<uint32_t> found_maxes,found_mins,found_saddles;
      for ( uint32_t y=0; y<size[1]; ++y ) {
        Value *i = values+y*stride[1]+z*stride[2];
        for ( uint32_t x=0; x<size[0]; ++x, i+=stride[0] ) {
          int32_t
            x0 = (x>0)         ? (compare_ptrs(i,i-stride[0])?1:-1) : 0,
            x1 = (x<size[0]-1) ? (compare_ptrs(i,i+stride[0])?1:-1) : 0,
            y0 = (y>0)         ? (compare_ptrs(i,i-stride[1])?1:-1) : 0,
            y1 = (y<size[1]-1) ? (compare_ptrs(i,i+stride[1])?1:-1) : 0,
            z0 = (z>0)         ? (compare_ptrs(i,i-stride[2])?1:-1) : 0,
            z1 = (z<int(size[2]-1)) ? (compare_ptrs(i,i+stride[2])?1:-1) : 0,
            ddx = x0+x1, //these are something like second derivatives
            ddy = y0+y1, //they are zero if i is not extremal in that direction
            ddz = z0+z1;

          bool critical = ddx!=0 && ddy!=0 && ddz!= 0;
          if (critical) {
            int index = 0-(ddx>>31)-(ddy>>31)-(ddz>>31); //arithmetic shift!
            switch(index) {
              case 0 :
                found_mins.push_back(i-values);
                break;
              case 1 :
              case 2 :
                found_saddles.push_back(i-values);
                break;
              case 3 :
                found_maxes.push_back(i-values);
                break;
              default: assert(0 && "wtf?");
            }
          }
        }
      }

      #pragma omp critical
      { //collect
        append_vector(maxes,found_maxes);
        append_vector(mins,found_mins);
        append_vector(voxel_saddles,found_saddles);
      }
    }
  }


  void mark_reachable_maxes()
  {
    std::cout << "mark reachable maxes (stack)" << std::endl;
    reachable_max.clear();
    reachable_max.resize(nvoxels,-1);
    int nmaxes = maxes.size();
    #pragma omp parallel for
    for ( int m=0; m<nmaxes; ++m ) {
      std::vector<uint32_t> stack;
      stack.push_back(maxes[m]);
      while (!stack.empty()) {
        uint32_t i = stack.back();
        Value *p = values+i;
        stack.pop_back();
        if ( reachable_max[i] == uint32_t(-1) ) {
          //at this point some other thread may have written to
          //reachable_max[i] but it doesn't really matter
          reachable_max[i] = maxes[m];
          uint32_t x,y,z;
          voxel_pos(i,x,y,z);
          if ( z>0         && compare_ptrs(p-stride[2],p) ) stack.push_back(i-stride[2]);
          if ( z<size[2]-1 && compare_ptrs(p+stride[2],p) ) stack.push_back(i+stride[2]);
          if ( y>0         && compare_ptrs(p-stride[1],p) ) stack.push_back(i-stride[1]);
          if ( y<size[1]-1 && compare_ptrs(p+stride[1],p) ) stack.push_back(i+stride[1]);
          if ( x>0         && compare_ptrs(p-stride[0],p) ) stack.push_back(i-stride[0]);
          if ( x<size[0]-1 && compare_ptrs(p+stride[0],p) ) stack.push_back(i+stride[0]);
        }
      }
    }
    for ( uint32_t i=0; i<nvoxels; ++i ) assert(reachable_max[i] != uint32_t(-1));
  }



  void mark_reachable_mins()
  {
    std::cout << "mark reachable mins (stack)" << std::endl;
    reachable_min.clear();
    reachable_min.resize(nvoxels,-1);
    int nmins = mins.size();
    #pragma omp parallel for
    for ( int m=0; m<nmins; ++m ) {
      std::vector<uint32_t> stack;
      stack.push_back(mins[m]);
      while (!stack.empty()) {
        uint32_t i = stack.back();
        Value *p = values+i;
        stack.pop_back();
        if ( reachable_min[i] == uint32_t(-1) ) {
          //at this point some other thread may have written to
          //reachable_min[i] but it doesn't really matter
          reachable_min[i] = mins[m];
          uint32_t x,y,z;
          voxel_pos(i,x,y,z);
          if ( z>0         && compare_ptrs(p,p-stride[2]) ) stack.push_back(i-stride[2]);
          if ( z<size[2]-1 && compare_ptrs(p,p+stride[2]) ) stack.push_back(i+stride[2]);
          if ( y>0         && compare_ptrs(p,p-stride[1]) ) stack.push_back(i-stride[1]);
          if ( y<size[1]-1 && compare_ptrs(p,p+stride[1]) ) stack.push_back(i+stride[1]);
          if ( x>0         && compare_ptrs(p,p-stride[0]) ) stack.push_back(i-stride[0]);
          if ( x<size[0]-1 && compare_ptrs(p,p+stride[0]) ) stack.push_back(i+stride[0]);
        }
      }
    }
    for ( uint32_t i=0; i<nvoxels; ++i ) assert(reachable_max[i] != uint32_t(-1));
  }

  void merge_sorted_voxels_and_saddles (
      std::vector<uint32_t> & order, //voxels already sorted.
      std::vector<Saddle*> & saddles_sorted ) //saddles sorted by compare_saddles
        //sorted_saddles will be merged into order
  {
    int vi = order.size()-1;
    int si = nonvoxel_saddles.size()-1;
    order.resize( order.size() + nonvoxel_saddles.size() );
    int end = order.size()-1;

    Saddle *s = saddles_sorted[si] ;
    while( vi >=0 && si >= 0 ) {

      switch(s->sort) {
        case Saddle::SortNaturally: {
          while( vi >= 0 && values[order[vi]] > s->value ) {
            order[end--] = order[vi--];
          }
          order[end--] = s->vert;
          --si;
          if (si>=0) s = saddles_sorted[si];
        }
        break;


        case Saddle::SortBefore: {
          while( vi >= 0 && order[vi] != s->whom ) {
            order[end--] = order[vi--];
          }
          uint32_t ovi = order[vi]; // save order[vi] as ovi, same as whom
          order[end--] = order[vi--];  // move whom into place
          while( s->whom == ovi && s->sort == Saddle::SortBefore) {
            // move all the saddles that sort before whom into place
            order[end--] = s->vert;
            if (si == 0) break;
            s = saddles_sorted[--si];
          }
        }
        break;


        case Saddle::SortAfter: {
          // this case is a bit complicated
          while ( vi >= 0 ) {
            // first scan down the list of grid verts until we find s->whom
            if ( order[vi] == s->whom ) {
              // due to careful choice of the initial sorting order
              // of saddles, we should find BBBBAAAAA at the end of
              // the saddle list, where the B's are the ones that
              // must come before whom, and the A's are those that
              // come after.
              while( s->whom == order[vi] && s->sort == Saddle::SortAfter ) {
                // move all the A's
                order[end--] = s->vert;
                if (si == 0) break;
                s = saddles_sorted[--si] ;
              }
              order[end--] = order[vi]; // move whom
              while( s->whom == order[vi] && s->sort == Saddle::SortBefore ) {
                // move all the B's
                order[end--] = s->vert;
                if (si == 0) break;
                s = saddles_sorted[--si] ;
              }
              --vi;

              break;
            } else {
              order[end--] = order[vi--];

            }
          }
        }
        break;

        default : abort();
      } //switch
    } //while
  }


  struct ComponentMap
  {
    typedef uint32_t key_type;
    typedef Component* value_type;
    Trilinear & t;
    Tourtre::SweepType type;
    ComponentMap( Trilinear & t_, Tourtre::SweepType type_ )
      : t(t_),type(type_) {}
    Component* & operator[]( uint32_t i )
    {
      if ( type == Tourtre::Join ) {
        if ( i < t.nvoxels ) return t.join_comps[i];
        else return t.saddle(i)->join_comp;
      } else { // type == Tourtre::Split
        if ( i < t.nvoxels ) return t.split_comps[i];
        else return t.saddle(i)->split_comp;
      }
    }
  };


  void clear_saddles()
  {
    for ( uint32_t i=0; i<nvoxels; ++i ) {
      Saddle *s = saddle_map[i];
      while(s) { Saddle *t=s; s=s->next; delete t; }
    }
    saddle_map = std::vector<Saddle*>();
  }

  float value( uint32_t i )
  {
    if ( i < nvoxels ) return values[i];
    else return saddle(i)->value;
  }
};

#endif


