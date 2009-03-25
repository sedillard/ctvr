#ifndef CONTOUR_TREE_HPP_INCLUDED
#define CONTOUR_TREE_HPP_INCLUDED

#include <stdint.h>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <tr1/unordered_map>

#include "Tourtre.hpp"
#include "Trilinear.hpp"


struct ContourTree
{
  typedef Tourtre::Node<uint32_t> Node;
  typedef Tourtre::Arc<uint32_t> Arc;
  typedef std::tr1::unordered_map<uint32_t,Node*> NodeMap;

  ContourTree ( uint8_t *image, uint32_t ncols, uint32_t nrows, uint32_t nstacks );

  Trilinear<uint8_t> tl;
  uint32_t img_size[3];
  uint32_t nvoxels;
  uint8_t *voxels; //not owned

  //Nodes of the contour tree are given ids in a contiguous range of integers.
 
  //Branches are also given ids, and there are roughly half as many branches.
  //The root branch is always id 0.
  
  //A branch is accessed by its root node, the place where it meets
  //the rest of the tree. You can follow a branch out from that node
  //by walking along arcs whose branch_root field points to that node.

  NodeMap node_map; //voxel -> node pointer
  std::vector< Node* > nodes; //node id -> node pointer
  std::vector< Arc* > arcs; //arc id -> arc pointer
  std::vector< Arc* > branches; //branch id -> pointer to first arc of branch

  void build(); //actually constructs the contour tree

  template <typename OutputItr> 
  void get_branch_children( uint32_t b, OutputItr out );
  
  bool branch_is_ascending( uint32_t b );

  std::pair<Node*,Node*> branch_range( uint32_t b );
    //returns saddle,extremum

  void prune_flat_arcs();

};


template<typename OutputItr>
void get_ascending_branch_children( ContourTree::Arc* arc, OutputItr out )
{
  while( arc ) {
    ContourTree::Arc *next = 0; 
    for ( ContourTree::Arc *a=arc->hi->up; a; a=a->next_up ) {
      if (a->branch==arc->branch) next=a;
      else *(out++) = a->branch;
    }
    for ( ContourTree::Arc *a=arc->hi->down; a; a=a->next_down ) {
      if ( a!=arc ) *(out++) = a->branch;
    }
    arc = next;
  }
}

template<typename OutputItr>
void get_descending_branch_children( ContourTree::Arc* arc, OutputItr out )
{
  while( arc ) {
    ContourTree::Arc *next = 0; 
    for ( ContourTree::Arc *a=arc->lo->down; a; a=a->next_down ) {
      if (a->branch==arc->branch) next=a;
      else *(out++) = a->branch;
    }
    for ( ContourTree::Arc *a=arc->lo->up; a; a=a->next_up ) 
      if ( a!=arc ) *(out++) = a->branch;
    arc = next;
  }
}

template <typename OutputItr>
void ContourTree::get_branch_children( uint32_t b, OutputItr out )
{
  Arc *first=branches[b];
  for ( Arc *a=first->hi->up; a; a=a->next_up ) {
    if (a->branch == b) {
      get_ascending_branch_children(first,out);
      return;
    }
  }
  for ( Arc *a=first->lo->down; a; a=a->next_down ) {
    if (a->branch == b) {
      get_descending_branch_children(first,out);
      return;
    }
  }
}


#endif
