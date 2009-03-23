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

  uint32_t nbranches;

  NodeMap node_map; //voxel -> node pointer
  std::vector< Node* > nodes; //node id -> node pointer
  std::vector< Node* > branches; //branch id -> node pointer
  std::vector< uint32_t > node_branch_id; //node id -> branch id

  void build(); //actually constructs the contour tree
  
  struct BranchChildItr 
  {
    ContourTree & ct;
    Node *broot, *c;
    bool ascending;

    BranchChildItr( ContourTree & ct_, uint32_t b ) : ct(ct_)
    {
      assert( b < ct.branches.size() );
      broot = c = ct.branches[b]; 
      for ( Arc *a=c->up; a; a=a->next_up ) 
        if ( a->branch_root == broot ) { ascending=true; break; } 
    }

    operator bool() const { return c; }
    uint32_t operator*() const { return ct.node_branch_id[c->id]; }

    BranchChildItr & operator++() 
    {
      if ( ascending ) {
        if (c->up) {
          for ( Arc *a=c->up; a; a=a->next_up ) 
            if ( a->branch_root == broot ) { c = a->hi; break; }
        } else {
          c = 0; 
        }
      } else {
        if (c->up) {
          for ( Arc *a=c->down; a; a=a->next_down ) 
            if ( a->branch_root == broot ) { c = a->lo; break; }
        } else {
          c = 0; 
        }
      }
      return *this;
    }
  };

};

#endif
