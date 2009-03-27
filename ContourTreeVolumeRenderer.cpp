

#include <cstdlib>
#include <cstring>
#include <algorithm>
#include <set>
#include <sstream>

#define GL_GLEXT_PROTOTYPES 1
#include <GL/gl.h>
#include <GL/glu.h>

#include "Tourtre.hpp"
#include "Trilinear.hpp"
#include "ContourTree.hpp"
#include "ContourTreeVolumeRenderer.hpp"
#include "HexMeshRayTracer.hpp"

#include <boost/foreach.hpp>

#include <boost/pool/pool_alloc.hpp>

using namespace std;
using namespace boost;
using namespace Tourtre;
using namespace Geom;
#define foreach BOOST_FOREACH

static 
void
wrap_tex_size ( GLuint size, GLuint tf_res, GLuint & nrows, GLuint & ncols )
{
    assert(size < 4096*4096);
    nrows = ncols = 1;
    while (nrows < size) nrows *= 2;
    while (nrows > ncols || ncols < tf_res) {
        ncols *= 2;
        nrows /= 2;
    }
    if (nrows < 2) nrows = 2;
}

static 
void* 
malloc0 ( uint32_t s )
{
    void* p = malloc(s);
    memset(p,0,s);
    return p;
}


//create a rainbow transfer function in c
static 
void rainbow_colors( RGBA8 *c, uint32_t n )
{
    float r,g,b,t,h;
    for (uint32_t i=0; i<n; ++i,++c) {
        t = i / (float)(n-1);
        h = 90.0f + t*270.0f;
        //h = 360.f * t;
        hls_to_rgb( h,0.4 + 0.2*t,1, r,g,b );
        //hls_to_rgb( h,0.5,1, &r,&g,&b );
        c->r = (uint8_t)(r*255.0f);
        c->g = (uint8_t)(g*255.0f);
        c->b = (uint8_t)(b*255.0f);
        c->a = (uint8_t)(t*255.0f) ;
    }
}



ContourTreeVolumeRenderer::ContourTreeVolumeRenderer 
(
  ContourTree & ct_, 
  uint8_t *voxels_, 
  uint32_t ncols, 
  uint32_t nrows, 
  uint32_t nstacks 
) : 
  ct(ct_),
  tl(ct_.tl),
  nvoxels(ncols*nrows*nstacks),
  voxels(voxels_)
{
  vol_size[0] = ncols;
  vol_size[1] = nrows;
  vol_size[2] = nstacks;
  assert(ncols   == tl.size[0]);
  assert(nrows   == tl.size[1]);
  assert(nstacks == tl.size[2]);

  gradients = 0;

  //compute branch decomposition properties
  compute_branch_properties();
  compute_branch_map();
  compute_max_shader_itrs(); 

  tf_res = 256;
  tf_size = 256; 

  //create the global transfer function
  global_tf_tex = (RGBA8*)malloc0( sizeof(RGBA8)*tf_res );;
  rainbow_colors( global_tf_tex, tf_res );
  
  uint32_t nbranches = ct.branches.size();
  branch_tf_offset.resize( nbranches );
  fill(branch_tf_offset.begin(),branch_tf_offset.end(),-1);

  //assign branch texture function offsets
  for ( uint32_t b=0; b<nbranches; ++b ) {
    branch_tf_offset[b] = tf_size;
    tf_size += branch_tf_size(b);
  }


  //compute size of wrapped tf and branch textures
  wrap_tex_size( tf_size, tf_res, tf_tex_nrows, tf_tex_ncols);
  wrap_tex_size( nbranches, tf_res, br_tex_nrows, br_tex_ncols);

  //compute the sizes of the 3D textures, rounding up to powers of 2
  for (int i=0; i<3; ++i) {
      scalar_tex_size[i] = 1;
      while( (scalar_tex_size[i]*=2) < vol_size[i] ) {};
  }

  //create the per-branch tf texture
  tf_tex = (RGBA8*)malloc( sizeof(RGBA8)*tf_tex_nrows*tf_tex_ncols );
  //init to black to see errors
  for ( uint32_t i=0; i<tf_tex_nrows*tf_tex_ncols; ++i ) {
    tf_tex[i] = (RGBA8){0,0,0,255}; 
  }
  
  init_branch_textures();
  default_tf();

}



GLushort ContourTreeVolumeRenderer::br_tex_x( uint32_t b )
{
    return ( b % br_tex_ncols ) * 
            (0x10000/br_tex_ncols);// + (0x8000/br_tex_ncols) ;
}

GLushort ContourTreeVolumeRenderer::br_tex_y( uint32_t b )
{
    return ( b / br_tex_ncols ) * 
            (0x10000/br_tex_nrows);// + (0x8000/br_tex_nrows) ;
}

GLushort ContourTreeVolumeRenderer::tf_tex_x( uint32_t offset )
{
    return ( offset % tf_tex_ncols ) * 
            (0x10000/tf_tex_ncols) + (0x8000/tf_tex_ncols) ;
}

GLushort ContourTreeVolumeRenderer::tf_tex_y( uint32_t offset )
{
    return ( offset / tf_tex_ncols ) * 
            (0x10000L/tf_tex_nrows) + (0x8000/tf_tex_nrows) ;
}


void
ContourTreeVolumeRenderer::init_branch_textures()
{
  uint32_t br_tex_size = br_tex_nrows * br_tex_ncols;

  parent_tex   = (GLushort*)malloc0( 2*sizeof(GLushort)*br_tex_size );
  tf_index_tex = (GLushort*)malloc0( 2*sizeof(GLushort)*br_tex_size );
  depth_tex    = (GLubyte*) malloc0(   sizeof(GLubyte )*br_tex_size );
  saddle_val_tex = (float*)malloc0( sizeof(float)*br_tex_size );

  branch_map_tex = (GLushort*)malloc0( 2*sizeof(GLushort)*nvoxels );
  // set branch_map_tex values, splitting into 1D index into 2D tex coords.
  {   
    GLushort *p = branch_map_tex;

    for (uint32_t i=0; i<nvoxels; ++i, p+=2) { 
        assert(branch_map[i] != uint32_t(-1));
        p[0] = br_tex_x(branch_map[i]);
        p[1] = br_tex_y(branch_map[i]);
    }
  }


  // initialize the textures
  uint32_t nbranches = ct.branches.size();
  for ( uint32_t b=0; b<nbranches; ++b ) {
    uint32_t id0 = 2*b, id1 = id0+1;
    
    if (b == 0) { //root
        parent_tex[id0] = parent_tex[id1] = 0; 
    } else {
        parent_tex[id0] = br_tex_x(branch_parent[b]);
        parent_tex[id1] = br_tex_y(branch_parent[b]);
    }

    int32_t tf_index = branch_tf_offset[b] -
                                (int32_t)(floor(branch_lo_val[b]*(tf_res-1))) + 1 ;
    tf_index_tex[id0] = tf_tex_x(tf_index);
    tf_index_tex[id1] = tf_tex_y(tf_index);

    depth_tex[b] = branch_depth[b];
    saddle_val_tex[b] = branch_saddle_value[b]/255.0f;
  }

  vector<uint32_t> offs = branch_tf_offset;
  sort(offs.begin(),offs.end());
  assert( unique(offs.begin(),offs.end()) == offs.end() );
}



pair<uint32_t,uint32_t> 
ContourTreeVolumeRenderer::branch_tf_bounds( uint32_t b )
{
  pair<uint32_t,uint32_t> out;
  out.first  = (uint32_t)( floor(branch_lo_val[b] * tf_res) );
  out.second = (uint32_t)( ceil (branch_hi_val[b] * tf_res) );
  assert( out.first <= out.second );
  if (out.first > 0) --out.first;
  if (out.second < 255) ++out.second;
  return out;
}

uint32_t ContourTreeVolumeRenderer::branch_tf_size( uint32_t b )
{
  pair<uint32_t,uint32_t> bounds = branch_tf_bounds(b);
    return bounds.second-bounds.first+1;
}



void
ContourTreeVolumeRenderer::compile_shader()
{
    fshader_id = glCreateShader(GL_FRAGMENT_SHADER);
    program_id = glCreateProgram();

    //pass this info into the shader as a #define
    stringstream shader_defines; 
    shader_defines << "#define MAX_ITRS " << max_shader_itrs << endl  
                  << "#define TF_RES " << tf_res << endl;
    
    string src_string = shader_defines.str() + fshader_src;
    //load shader src
    const char *src[] = {src_string.c_str()};
    const int len[] = {src_string.size()};
    
    glShaderSource( fshader_id,1,src,len );

    //compile
    glCompileShader(fshader_id);

    //check compilation status
    GLint status;
    glGetShaderiv(fshader_id, GL_COMPILE_STATUS, &status);
    if (status == GL_FALSE) {
        //report error     
        int len; 
        glGetShaderiv(fshader_id,GL_INFO_LOG_LENGTH,&len);
        char *err = (char*)malloc(len);
        glGetShaderInfoLog(fshader_id,len,0,err);
        cerr << "shader_init_gl: compilation error\n " << err << endl;
        free(err);
        return;
    }

    //link
    glAttachShader(program_id,fshader_id);
    glLinkProgram(program_id);
    glGetProgramiv(program_id,GL_LINK_STATUS,&status);
    if (status == GL_FALSE) {
        //report error     
        int len; 
        glGetProgramiv(fshader_id,GL_INFO_LOG_LENGTH,&len);
        char *err = (char*)malloc(len);
        glGetProgramInfoLog(fshader_id,len,0,err);
        cerr << "shader_init_gl: compilation error\n " << err << endl;
        free(err);
        return;
    }
  
    //get uniform variable locations
#define TV_GET_UL(X) X##_ul = glGetUniformLocation(program_id,#X)

    TV_GET_UL(scalar_tex);
    TV_GET_UL(global_tf_tex);
    TV_GET_UL(tf_tex);
    TV_GET_UL(tf_index_tex);
    TV_GET_UL(parent_tex);
    TV_GET_UL(depth_tex);
    TV_GET_UL(saddle_val_tex);
    TV_GET_UL(branch_map_tex);
    TV_GET_UL(br_tex_size);
    TV_GET_UL(tf_tex_size);
    TV_GET_UL(x_size);
    TV_GET_UL(y_size);
    TV_GET_UL(z_size);
    TV_GET_UL(x_inch);
    TV_GET_UL(y_inch);
    TV_GET_UL(z_inch);
    TV_GET_UL(light_vec);
    TV_GET_UL(view_vec);

#undef TV_GET_UL

}


static 
void set_nn_tex_env()
{
    glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_BASE_LEVEL,0);
    glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MAX_LEVEL,0);
    glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_WRAP_S,GL_CLAMP);
    glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_WRAP_T,GL_CLAMP);
    glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MAG_FILTER,GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MIN_FILTER,GL_NEAREST);
    glTexEnvi(GL_TEXTURE_ENV,GL_TEXTURE_ENV_MODE,GL_REPLACE);
}

static 
void set_linear_tex_env()
{
    glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_BASE_LEVEL,0);
    glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MAX_LEVEL,0);
    glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_WRAP_S,GL_CLAMP);
    glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_WRAP_T,GL_CLAMP);
    glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MAG_FILTER,GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MIN_FILTER,GL_LINEAR);
    glTexEnvi(GL_TEXTURE_ENV,GL_TEXTURE_ENV_MODE,GL_REPLACE);
}

static
void set_pixel_store( GLuint row_len )
{
    glPixelStorei(GL_UNPACK_ALIGNMENT,1);
    glPixelStorei(GL_UNPACK_SKIP_PIXELS,0);
    glPixelStorei(GL_UNPACK_ROW_LENGTH, row_len );
    glPixelStorei(GL_UNPACK_SKIP_ROWS,0);
}





void
ContourTreeVolumeRenderer::load_textures()
{
    //load tf tex
    glGenTextures(1,&tf_tex_id);
    glBindTexture(GL_TEXTURE_2D,tf_tex_id);
    set_nn_tex_env();
    set_pixel_store( tf_tex_ncols );

    glTexImage2D(
        GL_TEXTURE_2D, 0, GL_RGBA8, 
        tf_tex_ncols, tf_tex_nrows, 0, 
        GL_RGBA, GL_UNSIGNED_BYTE, tf_tex );
    
    
    glGenTextures(1,&global_tf_tex_id);
    glBindTexture(GL_TEXTURE_2D,global_tf_tex_id);
    set_linear_tex_env();
    set_pixel_store( tf_res );
    glTexImage2D(
        GL_TEXTURE_2D, 0, GL_RGBA8, 
        tf_res, 1, 0, 
        GL_RGBA, GL_UNSIGNED_BYTE, global_tf_tex );


    //load tf index tex
    glGenTextures(1,&tf_index_tex_id);
    glBindTexture(GL_TEXTURE_2D,tf_index_tex_id);
    set_nn_tex_env();
    set_pixel_store( br_tex_ncols );
    glTexImage2D(GL_TEXTURE_2D, 0, GL_LUMINANCE16_ALPHA16, 
        br_tex_ncols, br_tex_nrows, 0, 
        GL_LUMINANCE_ALPHA, GL_UNSIGNED_SHORT, tf_index_tex );

    free(tf_index_tex);
    tf_index_tex=0;

    //load parent tex
    glGenTextures(1,&parent_tex_id);
    glBindTexture(GL_TEXTURE_2D,parent_tex_id);
    set_nn_tex_env();
    set_pixel_store( br_tex_ncols );
    glTexImage2D(   GL_TEXTURE_2D, 0, GL_LUMINANCE16_ALPHA16, 
                    br_tex_ncols, br_tex_nrows, 0, 
                    GL_LUMINANCE_ALPHA, GL_UNSIGNED_SHORT, parent_tex );

    free(parent_tex);
    parent_tex=0;
    
    //load depth tex
    glGenTextures(1,&depth_tex_id);
    glBindTexture(GL_TEXTURE_2D,depth_tex_id);
    set_nn_tex_env();
    set_pixel_store( br_tex_ncols );
    glTexImage2D(   GL_TEXTURE_2D, 0, GL_ALPHA8,
                    br_tex_ncols, br_tex_nrows, 0, 
                    GL_ALPHA, GL_UNSIGNED_BYTE, depth_tex );

    free(depth_tex);
    depth_tex=0;

    //load saddle tex
    glGenTextures(1,&saddle_val_tex_id);
    glBindTexture(GL_TEXTURE_2D,saddle_val_tex_id);
    set_nn_tex_env();
    set_pixel_store( br_tex_ncols );
    glTexImage2D( GL_TEXTURE_2D, 0, GL_ALPHA16,
                    br_tex_ncols, br_tex_nrows, 0, 
                    GL_ALPHA, GL_FLOAT, saddle_val_tex );

    free(saddle_val_tex);
    saddle_val_tex = 0;

    //load data tex
    glGenTextures(1,&scalar_tex_id);
    glBindTexture(GL_TEXTURE_3D,scalar_tex_id);

    glTexParameteri(GL_TEXTURE_3D,GL_TEXTURE_BASE_LEVEL,0);
    glTexParameteri(GL_TEXTURE_3D,GL_TEXTURE_MAX_LEVEL,0);
    glTexParameteri(GL_TEXTURE_3D,GL_TEXTURE_WRAP_S,GL_CLAMP);
    glTexParameteri(GL_TEXTURE_3D,GL_TEXTURE_WRAP_T,GL_CLAMP);
    glTexParameteri(GL_TEXTURE_3D,GL_TEXTURE_WRAP_R,GL_CLAMP);
    glTexParameteri(GL_TEXTURE_3D,GL_TEXTURE_MAG_FILTER,GL_LINEAR);
    glTexParameteri(GL_TEXTURE_3D,GL_TEXTURE_MIN_FILTER,GL_LINEAR);
    glTexEnvi(GL_TEXTURE_ENV,GL_TEXTURE_ENV_MODE,GL_REPLACE);


    glTexImage3D( GL_TEXTURE_3D, 0, GL_ALPHA8, 
                  scalar_tex_size[0],
                  scalar_tex_size[1],
                  scalar_tex_size[2],
                  0, GL_ALPHA, GL_UNSIGNED_BYTE, 0 );

    glPixelStorei(GL_UNPACK_ALIGNMENT,1);
    glPixelStorei(GL_UNPACK_SKIP_PIXELS,0);
    glPixelStorei(GL_UNPACK_ROW_LENGTH, vol_size[0] );
    glPixelStorei(GL_UNPACK_SKIP_ROWS,0);
    glPixelStorei(GL_UNPACK_IMAGE_HEIGHT, vol_size[1] );
    glPixelStorei(GL_UNPACK_SKIP_IMAGES,0);
   
    glTexSubImage3D( GL_TEXTURE_3D, 0, 0,0,0,
                    vol_size[0], vol_size[1], vol_size[2],
                    GL_ALPHA, GL_UNSIGNED_BYTE, voxels );


    //load branch map texture 
    glGenTextures(1,&branch_map_tex_id);
    glBindTexture(GL_TEXTURE_3D,branch_map_tex_id);

    glTexParameteri(GL_TEXTURE_3D,GL_TEXTURE_BASE_LEVEL,0);
    glTexParameteri(GL_TEXTURE_3D,GL_TEXTURE_MAX_LEVEL,0);
    glTexParameteri(GL_TEXTURE_3D,GL_TEXTURE_WRAP_S,GL_CLAMP);
    glTexParameteri(GL_TEXTURE_3D,GL_TEXTURE_WRAP_T,GL_CLAMP);
    glTexParameteri(GL_TEXTURE_3D,GL_TEXTURE_WRAP_R,GL_CLAMP);
    glTexParameteri(GL_TEXTURE_3D,GL_TEXTURE_MAG_FILTER,GL_NEAREST);
    glTexParameteri(GL_TEXTURE_3D,GL_TEXTURE_MIN_FILTER,GL_NEAREST);
    glTexEnvi(GL_TEXTURE_ENV,GL_TEXTURE_ENV_MODE,GL_REPLACE);

    glTexImage3D( GL_TEXTURE_3D, 0, GL_LUMINANCE16_ALPHA16, 
                  scalar_tex_size[0],
                  scalar_tex_size[1],
                  scalar_tex_size[2],
                  0, GL_LUMINANCE_ALPHA, GL_UNSIGNED_SHORT, 0 );

    glPixelStorei(GL_UNPACK_ALIGNMENT,1);
    glPixelStorei(GL_UNPACK_SKIP_PIXELS,0);
    glPixelStorei(GL_UNPACK_ROW_LENGTH, vol_size[0] );
    glPixelStorei(GL_UNPACK_SKIP_ROWS,0);
    glPixelStorei(GL_UNPACK_IMAGE_HEIGHT, vol_size[1] );
    glPixelStorei(GL_UNPACK_SKIP_IMAGES,0);
   
    glTexSubImage3D( GL_TEXTURE_3D, 0, 0,0,0,
                     vol_size[0], vol_size[1], vol_size[2],
                     GL_LUMINANCE_ALPHA, GL_UNSIGNED_SHORT, 
                     branch_map_tex );
    
    free(branch_map_tex);
    branch_map_tex=0;
}


void
ContourTreeVolumeRenderer::update_tf_tex()
{
  cout << "update global tf " << endl;
  glBindTexture(GL_TEXTURE_2D,tf_tex_id);
  set_nn_tex_env();
  set_pixel_store( tf_tex_ncols );

  glTexImage2D(
      GL_TEXTURE_2D, 0, GL_RGBA8, 
      tf_tex_ncols, tf_tex_nrows, 0, 
      GL_RGBA, GL_UNSIGNED_BYTE, tf_tex );
}


void
ContourTreeVolumeRenderer::default_tf()
{
  //RGBA8 *colors = new RGBA8[tf_res];
  //rainbow_colors(colors,tf_res);
  
  for ( uint32_t b=0; b<ct.branches.size(); ++b ) {
    pair<uint32_t,uint32_t> bnds = branch_tf_bounds(b);
    srand(b);
    float h = 360.0f*rand()/float(RAND_MAX);
    uint32_t off = branch_tf_offset[b]+1;
    for ( uint32_t i=bnds.first; i!=bnds.second; ++i ) {
      float r,g,b;
      float t = i/float(tf_res-1);
      hls_to_rgb(h,0.2+0.8*t,1,r,g,b);
      RGBA8* c = tf_tex + off-bnds.first+i;
      c->r = uint8_t(r*255.0f);
      c->g = uint8_t(g*255.0f);
      c->b = uint8_t(b*255.0f);
      //c->a = t > 0.5 ?  uint8_t(t*255.0f) : 0;
      c->a = uint8_t(t*255.0f);
    }
  }

  //for ( uint32_t i=0; i<tf_size; ++i ) tf_tex[i] = (RGBA8){255,255,255,255};
}


void
ContourTreeVolumeRenderer::init_gl ()
{
  load_textures();
  compile_shader();
}



void
ContourTreeVolumeRenderer::enable_gl()
{
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);
    glEnable(GL_DEPTH_TEST);
    glDepthFunc(GL_LESS);
    glUseProgram(program_id);

#define TV_SET_TEX(I,N,X) {\
    glActiveTexture(GL_TEXTURE0+I); \
    glBindTexture(GL_TEXTURE_##N##D,X##_id); \
    glUniform1i(X##_ul,I); }
   
    TV_SET_TEX( 0, 3, scalar_tex );
    TV_SET_TEX( 1, 3, branch_map_tex );
    TV_SET_TEX( 2, 2, tf_tex );
    TV_SET_TEX( 3, 2, tf_index_tex );
    TV_SET_TEX( 4, 2, parent_tex );
    TV_SET_TEX( 5, 2, depth_tex );
    TV_SET_TEX( 6, 2, saddle_val_tex );
    TV_SET_TEX( 7, 2, global_tf_tex );
    
#undef TV_SET_TEX

    glUniform4f( tf_tex_size_ul,
            tf_tex_ncols,     tf_tex_nrows,
        1.0/tf_tex_ncols, 1.0/tf_tex_nrows );

    const GLuint *tz = scalar_tex_size;
    const uint32_t *sz = vol_size;
    glUniform2f( x_size_ul, tz[0], 1.0/(tz[0]) );
    glUniform2f( y_size_ul, tz[1], 1.0/(tz[1]) );
    glUniform2f( z_size_ul, tz[2], 1.0/(tz[2]) );

    glUniform3f( x_inch_ul, 1.0 /tz[0], 0, 0 );
    glUniform3f( y_inch_ul, 0, 1.0 /tz[1], 0 );
    glUniform3f( z_inch_ul, 0, 0, 1.0 /tz[2] );

    { //compute view vector to use as the light direction
        double f[3],b[3];
        GLint vp[4];
        glGetIntegerv(GL_VIEWPORT,vp);
        double mv[16],proj[16];
        glGetDoublev(GL_MODELVIEW_MATRIX,mv);
        glGetDoublev(GL_PROJECTION_MATRIX,proj);
        gluUnProject( vp[2]/2, vp[3]/2, 0, mv,proj,vp, f,f+1,f+2 );
        gluUnProject( vp[2]/2, vp[3]/2, 1, mv,proj,vp, b,b+1,b+2 );
        light_vec = vec3(b[0]-f[0], b[1]-f[1], b[2]-f[2] );
        light_vec.normalize();
        glUniform3f( light_vec_ul, light_vec[0],light_vec[1],light_vec[2] );
    }
    
    { //compute view vector for ray stepping
        double f[3],b[3];
        GLint vp[4];
        glGetIntegerv(GL_VIEWPORT,vp);
        double mv[16],proj[16];
        glGetDoublev(GL_MODELVIEW_MATRIX,mv);
        glGetDoublev(GL_PROJECTION_MATRIX,proj);
        gluUnProject( vp[2]/2, vp[3]/2, 0, mv,proj,vp, f,f+1,f+2 );
        gluUnProject( vp[2]/2, vp[3]/2, 1, mv,proj,vp, b,b+1,b+2 );
        double vx=b[0]-f[0], vy=b[1]-f[1], vz=b[2]-f[2];
        double norm = sqrt(vx*vx+vy*vy+vz*vz);
        double ss = 1.0/sz[0]; //step size
        float v[3] =  { vx*ss/norm, vy*ss/norm, vz*ss/norm };
        glUniform3f( view_vec_ul, v[0],v[1],v[2] );
    }
    
    //  Transform tex coords to go from this...
    //  |-*-|-*-|-*-|
    //  0...........tz
    //                     /tz
    //  0...........1
    //                    0.5
    //    0...........1    
    //                    *(tz-1)/tz
    //    0.......1
    //  |-*-|-*-|-*-|
    //  ... to this

    glActiveTexture(GL_TEXTURE0);
    glMatrixMode(GL_TEXTURE);
    glLoadIdentity();
    
    glScaled( 1.0/tz[0], 1.0/tz[1], 1.0/tz[2] );
    glTranslated( 0.5, 0.5, 0.5 );
    glScaled( (sz[0]-1.0)/sz[0], (sz[1]-1.0)/sz[1], (sz[2]-1.0)/sz[2] );
}



void 
ContourTreeVolumeRenderer::set_global_tf( RGBA8 *colors ) 
{
  copy(colors,colors+tf_res,global_tf_tex);
}


void 
ContourTreeVolumeRenderer::update_global_tf_tex() 
{
  glBindTexture(GL_TEXTURE_2D,global_tf_tex_id);
  set_linear_tex_env();
  set_pixel_store( tf_res );
  glTexImage2D(
      GL_TEXTURE_2D, 0, GL_RGBA8, 
      tf_res, 1, 0, 
      GL_RGBA, GL_UNSIGNED_BYTE, global_tf_tex );
}


void
ContourTreeVolumeRenderer::compute_branch_properties()
{
  uint32_t nbranches = ct.branches.size();
  cout << "computing properties for " << nbranches << " branches" << endl;
  branch_saddle_value.resize(nbranches,0);
  branch_parent.resize(nbranches,uint32_t(-1));
  branch_depth.resize(nbranches,0);
  branch_lo_val.resize(nbranches,HUGE_VAL);
  branch_hi_val.resize(nbranches,-HUGE_VAL);

  typedef pair<uint32_t,uint32_t> Pair; //branch id, depth
  vector<Pair> stack( 1, make_pair(0,0) );

  vector<uint32_t> children;
  while(!stack.empty()) {
    Pair & p = stack.back();
    uint32_t b = p.first; //branch
    uint32_t d = p.second; //depth
    stack.pop_back();
    branch_depth[b] = d;

    pair<ContourTree::Node*,ContourTree::Node*> range = ct.branch_range(b);

    branch_saddle_value[b] = tl.value(range.first->vertex);

    float lo_val = tl.value(range.first->vertex) / 255.0f, 
          hi_val = tl.value(range.second->vertex) / 255.0f;
    if ( lo_val > hi_val ) swap(lo_val,hi_val);
    
    //propagate lo/hi values to parents
    for (uint32_t p=b; p!=uint32_t(-1); p=branch_parent[p]) {
      branch_lo_val[p] = min(lo_val,branch_lo_val[p]); 
      branch_hi_val[p] = max(hi_val,branch_hi_val[p]); 
    }

    children.clear();
    ct.get_branch_children(b,back_inserter(children));
    for ( vector<uint32_t>::iterator itr=children.begin(); itr!=children.end(); ++itr ) {
      uint32_t c = *itr; //child
      assert( c < ct.branches.size() );
      branch_parent[c] = b;
      stack.push_back( make_pair(c,d+1) );
    }
  }
}


void
ContourTreeVolumeRenderer::compute_branch_map()
{
  cout << "computing branch map " << endl;
  branch_map.resize( nvoxels, uint32_t(-1) );
  Trilinear<uint8_t> & tl = ct.tl;

  #pragma omp parallel for
  for ( int i=0; i<int(nvoxels); ++i ) {
    uint32_t up_vert = tl.reachable_max[i];
    uint32_t down_vert = tl.reachable_min[i];
    ContourTree::Node *up_node = ct.node_map[up_vert];
    ContourTree::Node *down_node = ct.node_map[down_vert];
    uint32_t up_branch = up_node->down->branch;
    uint32_t down_branch = down_node->up->branch;
    uint32_t b = select_branch(up_branch,down_branch, voxels[i]);
    while( branch_lo_val[b] == branch_hi_val[b] ) b=branch_parent[b];
    assert(b != uint32_t(-1));
    branch_map[i] = b;
  }

  //free up some memory
  tl.reachable_max = tl.reachable_min = vector<uint32_t>();
}



uint32_t
ContourTreeVolumeRenderer::branch_distance( uint32_t a_, uint32_t b_ )
{
  uint32_t a=a_,b=b_;
  int d=0;
  while (a!=b) {
    assert( a != 0 || b != 0 );
    if (branch_depth[a] > branch_depth[b]) {
      a = branch_parent[a];
    } else {
      b = branch_parent[b];
    }
    ++d;
  }
  return d;
}





void
ContourTreeVolumeRenderer::compute_max_shader_itrs()
{
    max_shader_itrs = 0;
    const uint32_t *size = vol_size;
    uint32_t im_skip = vol_size[0]*vol_size[1];
    
    #pragma omp parallel for
    for (int z=0; z<(int)vol_size[2]-1; ++z) {
      uint32_t max_itrs = 0;
      for (uint32_t y=0; y<vol_size[1]-1; ++y) {
        uint32_t o = (z*vol_size[1]+y)*vol_size[0];

        //vertices of a cell
        uint32_t v[8] = 
          { o, o+1, o+size[0], o+1+size[0],
            o+im_skip, o+im_skip+1, o+im_skip+size[0],
            o+im_skip+size[0]+1 };

        uint32_t b[8]; //branches at vertices
        for (uint32_t x=0; x<size[0]-1; ++x) {
            for (int i=0; i<8; ++i) b[i] = branch_map[v[i]]; 
            for ( int i=0; i<7; ++i) 
            for ( int j=i+1; j<8; ++j ) {
                uint32_t dist = branch_distance(b[i],b[j]);
                max_itrs = dist > max_itrs ? dist : max_itrs;
            }
            
            for (int i=0; i<8; ++i) ++v[i];
        }
      }
      #pragma omp critical 
      max_shader_itrs = max(max_shader_itrs,max_itrs); 
    }
    
    if ( max_shader_itrs > 16 )  {
        cout << "NOTE: cutting off maximum shader iterations at 16. This"
               " may cause some artifacts\n" << endl;
        max_shader_itrs = 16;
    }
}


RGBA8* ContourTreeVolumeRenderer::arc_tf_offset( ContourTree::Arc* a )
{
  uint32_t b = a->branch;
  assert(b<branch_lo_val.size());
  float lo = branch_lo_val[b];
  
  return reinterpret_cast<RGBA8*>(tf_tex) + 
           ( branch_tf_offset[b] - uint32_t(floor(lo*(tf_res-1))) );
}

pair<uint32_t,uint32_t> 
ContourTreeVolumeRenderer::arc_tf_bounds(ContourTree::Arc* arc)
{
  uint32_t b = arc->branch;
  assert(b<branch_lo_val.size());
  float lo = branch_lo_val[b];

  float a_lo = floor( (tl.value(arc->lo->vertex)/255.0f)*(tf_res-1) );
  float a_hi = ceil ( (tl.value(arc->hi->vertex)/255.0f)*(tf_res-1) );

  uint32_t o = branch_tf_offset[b] - uint32_t(floor(lo*(tf_res-1)));;
  return make_pair(o+uint32_t(a_lo), o+uint32_t(a_hi));
}

void ContourTreeVolumeRenderer::test_tf()
{
  for ( uint32_t i=0; i<ct.nodes.size(); ++i ) {
    ContourTree::Node* n = ct.nodes[i]; 
    for ( ContourTree::Arc *a=n->up; a; a=a->next_up ) {
      pair<uint32_t,uint32_t> bounds = arc_tf_bounds(a);

      for ( RGBA8* c = tf_tex+bounds.first; c!=tf_tex+bounds.second; ++c ) {
        assert(c-tf_tex < ptrdiff_t(tf_size));
        c->r = c->g = c->b = c->a = 0; 
      }

      if ( bounds.second - bounds.first > 10 ) {
        RGBA8 *c = tf_tex + (bounds.first+bounds.second)/2;
        c->r = c->g = c->b = c->a = 128; 
      }
    }
  }
}

uint32_t ContourTreeVolumeRenderer::node_cluster_find( uint32_t a )
{
  uint32_t c=node_cluster[a], t;
  while (node_cluster[c]!=c) c=node_cluster[c];
  while ( a != c ) t=node_cluster[a],node_cluster[a]=c,a=t;
  return c;
}

void ContourTreeVolumeRenderer::node_cluster_union( uint32_t a, uint32_t b )
{
  a = node_cluster_find(a);
  b = node_cluster_find(b);
  node_cluster[a] = b; 
}

void ContourTreeVolumeRenderer::read_events_file( const char* filename )
{

  ifstream file(filename,ios::in);
  if (!file) {
    cerr << "could not open events file " << filename << endl;
    abort();
  }
  events.clear();
  int from,to,when;
  for (;;) {
    file >> from >> to >> when;
    if (file.eof()) break;
    ContourTree::Node 
      *anode = ct.node_map[from], 
      *bnode = ct.node_map[to];
    if (!anode) {
      cerr << "events file references extremum at id " << from << " but I don't see one there" << endl; 
      abort();
    }
    if (!bnode) {
      cerr << "events file references extremum at id " << to << " but I don't see one there" << endl; 
      abort();
    }
    events.push_back((Event){from,to,when});
  }
}

void ContourTreeVolumeRenderer::do_merge_events( int until ) 
{
  node_cluster.resize( ct.nodes.size() );
  for ( size_t i=0; i<ct.nodes.size(); ++i ) 
    node_cluster[i] = i; 
  
  int nmerges=0;
  foreach( Event & e, events ) {
    if (e.when >= until ) break;
    ContourTree::Node 
      *anode = ct.node_map[e.from], 
      *bnode = ct.node_map[e.to];
    uint32_t a=anode->id, b=bnode->id;
    if ( node_cluster_find(a) != node_cluster_find(b) ) {
      ++nmerges;
      node_cluster_union(a,b); 
    }
  }
  cout << "merged " << nmerges << " nodes (up to time " << until << ")" << endl;

  propagate_cluster_info();
}


void ContourTreeVolumeRenderer::propagate_cluster_info()
{
  deque<Node*> nodeq;
  foreach( Node *n, ct.nodes ) {
    if (n->is_max()) nodeq.push_back(n->down->lo);
    else if (n->is_min()) nodeq.push_back(n->up->hi);
    else node_cluster[n->id] = uint32_t(-1);
  }
  
  while (!nodeq.empty()) {
    Node *n = nodeq.front();
    nodeq.pop_front();
    
    //look at all the neighbors and see how many -1's there are
    //if there is only one, and the others agree, then mark this node
    //with that label and push the node with the -1 back on the queue
    
    Node *odd_man_out=0;
    uint32_t label=-1;
    bool stop=false;
    for ( Arc *a=n->up; a; a=a->next_up ) {
      if ( node_cluster[a->hi->id] == uint32_t(-1) ) {
        if ( !odd_man_out ) {
          odd_man_out = a->hi;
        } else {
          //stop here because there are two neighbors with -1 labels
          stop=true;
          break;
        }
      } else {
        if (label == uint32_t(-1)) label = node_cluster[a->hi->id];
        else if (label != node_cluster[a->hi->id] ) {
          //stop here because the labels don't agree
          stop=true;
          break;
        }
      }
    }
    for ( Arc *a=n->down; a; a=a->next_down ) {
      if ( node_cluster[a->lo->id] == uint32_t(-1) ) {
        if ( !odd_man_out ) {
          odd_man_out = a->lo;
        } else {
          //stop here because there are two neighbors with -1 labels
          stop=true;
          break;
        }
      } else {
        if (label == uint32_t(-1)) label = node_cluster[a->lo->id];
        else if (label != node_cluster[a->lo->id] ) {
          //stop here because the labels don't agree
          stop=true;
          break;
        }
      }
    }
    if (!stop) {
      node_cluster[n->id] = label;
      nodeq.push_back(odd_man_out);
    }
  }
}

/*
        ....
         \ /
          d
   a     /
    \   b
     \ / \
      c   e 
      |  / \ 
      | f   g
   --------- <- isovalue cut here
      | 
      |
*/


#if 0
void ContourTreeVolumeRenderer::propagate_cluster_info()
{
  vector<int16_t> valence(ct.nodes.size(),0);
  typedef pair<Arc*,Direction> Pair;
  deque<Pair> nodeq;
  set<uint32_t> cluster_ids;
  foreach( Node *n, ct.nodes ) {
    if (n->is_max()) {
      nodeq.push_back(make_pair(n->down,Down));
      cluster_ids.insert(node_cluster[n->id]);
    } else if (n->is_min()) {
      nodeq.push_back(make_pair(n->up,Up));
      cluster_ids.insert(node_cluster[n->id]);
    } else {
      node_cluster[n->id] = uint32_t(-1);
      valence[n->id] = n->up_degree() + n->down_degree();
    }
  }

  cout << nodeq.size() << " extrema have " << cluster_ids.size() << " unique ids" << endl;

  typedef pool_allocator< pair<uint32_t,int>, default_user_allocator_new_delete,
                         details::pool::null_mutex> alloc;
  typedef map<uint32_t,int,less<uint32_t>,alloc > Map; 
          //cluster id,num votes

  while(!nodeq.empty()) {
    Arc *a = nodeq.front().first;
    Direction dir = nodeq.front().second;
    nodeq.pop_front();
    Node *n = dir==Up ? a->hi : a->lo;
    --valence[n->id];
    if ( valence[n->id]==1 && node_cluster[n->id]==uint32_t(-1) ) {
      Map votes;
      for ( Arc *a=n->up; a; a=a->next_up ) 
        ++votes[node_cluster[a->hi->id]];
      for ( Arc *a=n->down; a; a=a->next_down ) 
        ++votes[node_cluster[a->lo->id]];

      uint32_t best=-1;
      int most_votes=-1;
      foreach( Map::value_type v, votes ) {
        if ( v.second > most_votes ) {
          best = v.first, most_votes=v.second;
        }
      }
      node_cluster[n->id] = best;

      for ( Arc *a=n->up; a; a=a->next_up ) 
        if ( node_cluster[a->hi->id]==uint32_t(-1) ) 
          nodeq.push_back(make_pair(a,Up));
      for ( Arc *a=n->down; a; a=a->next_down ) 
        if ( node_cluster[a->lo->id]==uint32_t(-1) ) 
          nodeq.push_back(make_pair(a,Down));
    }
  }

  foreach( uint32_t id, node_cluster ) assert(id!=uint32_t(-1));
}
#endif


void ContourTreeVolumeRenderer::cluster_tf()
{
  cout << "setting automatic transfer function based on cluster" << endl;
  foreach( Arc* a, ct.arcs ) {
    pair<uint32_t,uint32_t> bounds = arc_tf_bounds(a);
    
    if ( node_cluster[a->hi->id] == node_cluster[a->lo->id] ) {
      srand(node_cluster[a->hi->id]);
      double hue = 360* rand()/double(RAND_MAX);
      float r,g,b;
      hls_to_rgb(hue,0.5,1,r,g,b);
      for ( RGBA8* c = tf_tex+bounds.first; c!=tf_tex+bounds.second; ++c ) {
        assert(c-tf_tex < ptrdiff_t(tf_size));
        c->r = uint8_t(255.0*r);
        c->g = uint8_t(255.0*g);
        c->b = uint8_t(255.0*b);
        c->a = 128;
      }
    } else {
      //i forget what this is supposed to be doing
      for ( RGBA8* c = tf_tex+bounds.first; c!=tf_tex+bounds.second; ++c ) {
        assert(c-tf_tex < ptrdiff_t(tf_size));
        c->r = c->g = c->b = c->a = 0; 
      }
    }
  }
}




template <typename T> inline 
double lerp( const T & a, const T & b, const double & t )
{ return (1-t)*double(a) + t*double(b) ; }

static void
maybe_normalize( Vec3d & v )
{
  double n = sqr_norm(v);
  if (n!=0) v/=n;
}


void ContourTreeVolumeRenderer::sw_render
( RGBA8 *image,
  uint32_t img_width,
  uint32_t num_rows,
  int win_left,
  int win_bottom
)
{
  if (!gradients) compute_gradients();

  GLdouble mv[16],proj[16];
  GLint vp[4];
  glGetDoublev(GL_MODELVIEW_MATRIX,mv);
  glGetDoublev(GL_PROJECTION_MATRIX,proj);
  glGetIntegerv(GL_VIEWPORT,vp);

  //#pragma omp parallel for
  #if 0
  for (int j = 0; j < int(num_rows); j++)
  for (uint32_t i = 0; i < img_width; i++) {

    Vec3d fp, bp;
    gluUnProject(i+win_left,j+win_bottom,0,mv,proj,vp,&fp[0],&fp[1],&fp[2]);
    gluUnProject(i+win_left,j+win_bottom,1,mv,proj,vp,&bp[0],&bp[1],&bp[2]);

    RGBAf result = {0,0,0,0};
    
    Vec3d org = fp;
    Vec3d dir = normalize(bp-fp);

    int nsteps=0;
    for( HexMeshRayTracer<uint8_t> tracer(voxels,vol_size,&org[0],&dir[0]); 
          tracer; ++tracer ) 
    {
      ++nsteps;

      //get points at front and back of ray segmen
      Vec3d pFront = org + tracer.tFront()*dir,
            pBack = org + tracer.tBack()*dir;
      
      //put points in cell-local coordinates
      for (int d=0; d<3; ++d) {
        pFront[d] -= floor(pFront[d]);
        pBack[d] -= floor(pBack[d]);
      }

      //do trilinear interpolation
      double tlFront[7], tlBack[7]; 
        //these store the intermediate values of trilinear interpolation
        //which we use for finding a monotone path from the sample points
        //to cell vertices

      tlFront[3] = lerp( tracer.F[0], tracer.F[1], pos[0]);
      tlFront[4] = lerp( tracer.F[2], tracer.F[3], pos[0]);
      tlFront[5] = lerp( tracer.F[4], tracer.F[5], pos[0]);
      tlFront[6] = lerp( tracer.F[6], tracer.F[7], pos[0]);
      tlFront[1] = lerp( tlFront[3] , tlFront[4] , pos[1]);
      tlFront[2] = lerp( tlFront[5] , tlFront[6] , pos[1]);
      tlFront[0] = lerp( tlFront[1] , tlFront[2] , pos[2]);
      double fFront = tlFront[0]; //value at the front point

      tlBack[3] = lerp( tracer.F[0], tracer.F[1], pos[0]);
      tlBack[4] = lerp( tracer.F[2], tracer.F[3], pos[0]);
      tlBack[5] = lerp( tracer.F[4], tracer.F[5], pos[0]);
      tlBack[6] = lerp( tracer.F[6], tracer.F[7], pos[0]);
      tlBack[1] = lerp( tlBack[3] , tlBack[4] , pos[1]);
      tlBack[2] = lerp( tlBack[5] , tlBack[6] , pos[1]);
      tlBack[0] = lerp( tlBack[1] , tlBack[2] , pos[2]);
      double fBack = tlBack[0]; //value at the back point

      //find a monotone path from the sample points to cell vertices,
      //in order to locate the branches containing the samples
      uint32_t maxVert=0, minVert=0;
      uint32_t maxTrace=0, minTrace=0;
      double *tlMin, *tlMax;

      //start the search for the lower vertex from the lower sample point, and
      //likewise for the upper vertex
      if (fFront < fBack) {	
        tlMin = tlFront; 
        tlMax = tlBack;
      } else {
        tlMin = tlBack;
        tlMax = tlFront;
      }
              
      //find vertex below
      if ( tlMin[2] > tlMin[1] ) {
        minTrace = 1;
      } else {
        minVert += 4;
        minTrace = 2;
      }
      if ( tlMin[ (minTrace<<1)+1 ] > tlMin[(minTrace<<1)+2] )  {
        minVert += 2;
      }
      if ( voxels[tracer.verts[minVert+1]] < voxels[tracer.verts[minVert]] ) {
        minVert += 1;
      }

      //find vertex above
      if ( tlMax[2] > tlMax[1] ) {
        maxVert += 4;
        maxTrace = 2;
      } else {
        maxTrace = 1;
      }
      if ( tlMax[ (maxTrace<<1)+1] < tlMin[(maxTrace<<1)+2] ) {
        maxVert += 2;
      }
      if ( voxels[tracer.verts[maxVert+1]] > voxels[tracer.verts[maxVert]] ) {
        maxVert+= 1;
      }

      //get branches
      uint32_t br_up   = branch_map[tracer.verts[maxVert]];
      uint32_t br_down = branch_map[tracer.verts[minVert]];

      //get gradients at front and back points
      Vec3d fGrad, bGrad;

      fGrad[0] = sample_gradient(tracer.verts,0,&pos[0]);
      fGrad[1] = sample_gradient(tracer.verts,1,&pos[0]);
      fGrad[2] = sample_gradient(tracer.verts,2,&pos[0]);

      bGrad[0] = sample_gradient(tracer.verts,0,&pos[0]);
      bGrad[1] = sample_gradient(tracer.verts,1,&pos[0]);
      bGrad[2] = sample_gradient(tracer.verts,2,&pos[0]);

      //maybe_normalize(fGrad);
      //maybe_normalize(bGrad);

      { // Composite

        if ( fabs(fFront-fBack) < 1 ) {
          //segment has only one transfer function entry
        
        }
        //double f = 0.5 * (fFront+fBack);
        int f0 = int(f);
        int f1 = f0 + 1;
        double frac = f - f0;

        
        double length = tracer.tBack() - tracer.tFront();

        uint32_t br = select_branch(up, down, f );
        uint32_t off = branch_tf_offset[br] - uint32_t(branch_lo_val[br]*(tf_res-1));
        const RGBA8 &c0 = tf_tex[off+int32_t(f0)+1];
        const RGBA8 &c1 = tf_tex[off+int32_t(f1)+1];

        float r = lerp( float(c0.r)/255.0f, float(c1.r)/255.0f, frac ); 
        float g = lerp( float(c0.g)/255.0f, float(c1.g)/255.0f, frac ); 
        float b = lerp( float(c0.b)/255.0f, float(c1.b)/255.0f, frac ); 
        float a = lerp( float(c0.a)/255.0f, float(c1.a)/255.0f, frac ); 


        float global = global_tf_tex[int(f0)].a/255.0f;
        a *= global;

        a = 1.0f - powf( 1.0f - a, length);
        float alpha = (1.0f-result.a)*a;

        //cout << "a = " << a << ", result.a = " << result.a << endl;
        result.r += alpha * r;
        result.g += alpha * g;
        result.b += alpha * b;
        result.a += alpha;
      } // Composite

      //if (result.a > 0.95) break;
    }

    //white background
    result.r += 1-result.a;
    result.g += 1-result.a;
    result.b += 1-result.a;

    RGBA8& out = image[j*img_width+i];
    out.r = uint8_t(255.0*result.r);
    out.g = uint8_t(255.0*result.g);
    out.b = uint8_t(255.0*result.b);
    out.a = uint8_t(255.0*result.a);
  }
  #endif
}


//val is in the range 0..255
uint32_t 
ContourTreeVolumeRenderer::select_branch( uint32_t up, uint32_t down, float val )
{
  for (;;) {
    if (up==down) return up;
    if ( branch_depth[up] > branch_depth[down] ) {
      if ( branch_saddle_value[up] < val ) return up; 
      else up=branch_parent[up];
    } else {
      if ( branch_saddle_value[down] > val ) return down; 
      else down=branch_parent[down];
    }
  }
}




void ContourTreeVolumeRenderer::compute_gradients()
{
  gradients = new int16_t[3*nvoxels];
  uint32_t ystride = vol_size[0], zstride = vol_size[0]*vol_size[1];
  #pragma omp parallel for
  for ( int z=0; z<int(vol_size[2]); ++z ) {
    uint32_t i=z*vol_size[0]*vol_size[1];
    for ( uint32_t y=0; y<vol_size[1]; ++y ) 
    for ( uint32_t x=0; x<vol_size[0]; ++x,++i ) {
      
      gradients[3*i+0] = 
        int16_t( x<vol_size[0]-1 ? voxels[i+1] : voxels[i] ) -
        int16_t( x>0 ? voxels[i-1] : voxels[i] );

      gradients[3*i+1] = 
        int16_t( y<vol_size[1]-1 ? voxels[i+ystride] : voxels[i] ) -
        int16_t( y>0 ? voxels[i-ystride] : voxels[i] );

      gradients[3*i+2] = 
        int16_t( uint32_t(z)<vol_size[2]-1 ? voxels[i+zstride]:voxels[i] ) -
        int16_t( uint32_t(z)>0 ? voxels[i-zstride] : voxels[i] );
    }
  }
}


double ContourTreeVolumeRenderer::sample_gradient( uint32_t v[8], int d, double x[3] )
{
  double t[7]; 
  t[3] = lerp( gradients[3*v[0]+d], gradients[3*v[1]+d], x[0]);
  t[4] = lerp( gradients[3*v[2]+d], gradients[3*v[3]+d], x[0]);
  t[5] = lerp( gradients[3*v[4]+d], gradients[3*v[5]+d], x[0]);
  t[6] = lerp( gradients[3*v[6]+d], gradients[3*v[7]+d], x[0]);
  t[1] = lerp( t[3] , t[4] , x[1]);
  t[2] = lerp( t[5] , t[6] , x[1]);
  t[0] = lerp( t[1] , t[2] , x[2]);
  return t[0];
}




