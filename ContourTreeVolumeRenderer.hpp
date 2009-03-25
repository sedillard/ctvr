#ifndef CONTOUR_TREE_VOLUME_RENDERER_HPP_INCLUDED
#define CONTOUR_TREE_VOLUME_RENDERER_HPP_INCLUDED

#include <string>
#include <vector>
#include <GL/gl.h>
#include "ContourTree.hpp"
#include "Trilinear.hpp"
#include "Color.hpp"

struct ContourTreeVolumeRenderer
{
  typedef ContourTree::Node Node;
  typedef ContourTree::Arc Arc;
  ContourTree & ct;
  Trilinear<uint8_t> & tl;

  uint32_t vol_size[3];
  uint32_t nvoxels;
  uint8_t *voxels; //not owned

  ContourTreeVolumeRenderer (
      ContourTree &, uint8_t *voxels,
      uint32_t ncols, uint32_t nrows, uint32_t nstacks );

  void init_gl();
  void enable_gl();


    //^^^^^^^^^^^^^^^^^^^/
    // branch properties /
  //___________________/
  

  std::vector< uint32_t > branch_parent; //branch id -> branch id
  std::vector< uint32_t > branch_depth; //branch id -> int
  std::vector< float > branch_saddle_value; // branch id -> value
  std::vector< uint32_t > branch_map; // voxel -> branch id

  std::vector<float> branch_lo_val, branch_hi_val;
      //the lo and hi values of a branch, in the range [0..1]
      //branch id -> value

  std::vector<uint32_t> branch_tf_offset;
      //branch id -> offset into tf array (tf_tex_master) where that branch's
      //tf starts

  std::vector<uint32_t> node_cluster; // node id -> node id

    //^^^^^^^^^^^^^^^^^^^^^^^^^/
    // transfer function stuff /
  //_________________________/
  uint32_t tf_res; //transfer function resolution: how many entries

  uint32_t tf_size; // total number of transfer function 
                  // entries for all branches
                  
  GLuint tf_tex_nrows, tf_tex_ncols;
      // size of tf_tex, tf_tex_master

  GLuint br_tex_nrows, br_tex_ncols;
      // size of parent_tex, tf_index_tex, depth_tex, saddle_tex
      
  GLushort* parent_tex; 
      // 2 channels, stores x,y texcoords of parent branch

  GLushort* tf_index_tex;
      // 2 channels, stores x,y texcoords of a branch's tf
      // IMPORTANT: the lowest value of the branch is subtracted from 
      // the x-coord so as to yield the correct tf entry when added to
      // the fragments function value

  GLubyte* depth_tex;
      // 1 channel, stores number of parents (0,255)
  
  float* saddle_val_tex;
      // 1 channel, stores branches' saddle values.

  RGBA8* tf_tex;
      // 4 channels (RGBA), stores the per-branch transfer functions
      
  GLushort* branch_map_tex;
      // 2 channels, stores mapping from trilinear mesh vertices to 
      //ranches. 

  GLuint scalar_tex_size[3]; 
      //size of the 3D image, rounded up to 
      //powers of 2

  RGBA8 *global_tf_tex;
      // The global tf
  
  GLubyte *aux_scalar_tex; 
      //the auxiliary scalar texture (eg, original/unsmoothed data)

  GLubyte *surface_val_tex;
      //for each branch, yields the value of the surface there

    //^^^^^^^^^^^^^^/ 
    // Shader stuff / 
  //______________/ 
  GLuint program_id;  // The combined shader's program name, to be used with
                      //   glUseProgram
  GLuint fshader_id;  // The fragment shader's name
  std::string  fshader_src; // The fragment shader source
  uint32_t fshader_len; // The length of the source 

  //Texture objects ids or 'names' as GL calls them
  GLuint scalar_tex_id;
  GLuint aux_scalar_tex_id;
  GLuint global_tf_tex_id;
  GLuint branch_map_tex_id;
  GLuint tf_tex_id;  
  GLuint tf_index_tex_id;
  GLuint parent_tex_id;
  GLuint depth_tex_id;
  GLuint saddle_val_tex_id;
  GLuint surface_val_tex_id;

  //The locations of the shader sampler uniform variables
  GLuint scalar_tex_ul;
  GLuint aux_scalar_tex_ul;
  GLuint global_tf_tex_ul;
  GLuint branch_map_tex_ul;
  GLuint tf_tex_ul;  
  GLuint tf_index_tex_ul;
  GLuint parent_tex_ul;
  GLuint depth_tex_ul;
  GLuint saddle_val_tex_ul;
  GLuint surface_val_tex_ul;

  //Some other uniform variables
  GLuint br_tex_size_ul;
  GLuint tf_tex_size_ul;
  GLuint x_size_ul, y_size_ul, z_size_ul;
  GLuint x_inch_ul, y_inch_ul, z_inch_ul;
  GLuint light_vec_ul;
  GLuint view_vec_ul; 

  uint32_t max_shader_itrs;

    
  //functions to compute things like branch_parent, branch_saddle_val, etc.
  void compute_branch_properties();
  uint32_t branch_distance( uint32_t a, uint32_t b );
  void compute_branch_map();
  void compute_max_shader_itrs();

  void compile_shader();

  //create the branch textures
  void init_branch_textures(); 

  //Upload the branch textures, and free the memory. If you need to do it
  //again call ini_branch_textures.
  void load_textures(); 

  //transform a 1D branch/tf index into a 2D texture coordinate
  GLushort br_tex_x( uint32_t branch_id );
  GLushort br_tex_y( uint32_t branch_id );
  GLushort tf_tex_x( uint32_t offset );
  GLushort tf_tex_y( uint32_t offset );
  
  //the 1D bounds and size of a branch's transfer function
  std::pair<uint32_t,uint32_t> branch_tf_bounds( uint32_t b );
  uint32_t branch_tf_size( uint32_t b );

  //Returns a pointer p such that p + (v*tf_res) yeilds the tf entry for
  //value v along the arc.
  RGBA8* arc_tf_offset( ContourTree::Arc* );
        
  //Returns two offsets into the transfer function array yielding the range
  //of tf entries associated with the arc. End not included, stl-style.
  std::pair<uint32_t,uint32_t> arc_tf_bounds(ContourTree::Arc* arc);

  //node ids to be unioned
  void node_cluster_union(uint32_t, uint32_t); 

  //retrieve the cluster representative, by node id
  uint32_t node_cluster_find( uint32_t ) ; 
  
  //process an extremum tracker file
  void read_tracker_file(const char* filename);
  void propagate_cluster_info();
  void cluster_tf();
  
  void test_tf();

  //after changing the transfer function, re-upload to GL
  void update_tf_tex();

  void set_global_tf( RGBA8 *colors );
  void update_global_tf_tex();
};

#endif
