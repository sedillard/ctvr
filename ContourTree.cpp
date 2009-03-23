
#include "ContourTree.hpp"

using namespace std;
using namespace Tourtre;

#if 0
static
void assert_all_unique( vector<uint32_t> const& v )
{
  vector<uint32_t> x=v;
  sort(x.begin(),x.end());
  for ( size_t i=0; i<x.size()-1; ++i ) assert(x[i] != x[i+1] );
}
#endif


struct JoinSweepClosure
{
  Trilinear<uint8_t> & t;
  typedef uint32_t Vertex;
  int max_link_size;
  JoinSweepClosure( Trilinear<uint8_t> & t_ ) : t(t_),max_link_size(32) {}

  int lower_link( const Vertex & i, Vertex link[] ) const 
  { 
    int nlink = t.lower_link(i,link);
    return nlink;
  }

  Vertex walk_back( const Vertex & v ) const 
  { 
    assert(v < t.nvoxels);
    assert(t.reachable_min[v] != uint32_t(-1));
    return t.reachable_min[v]; 
  }
};

struct SplitSweepClosure
{
  Trilinear<uint8_t> & t;
  typedef uint32_t Vertex;
  int max_link_size;
  SplitSweepClosure( Trilinear<uint8_t> & t_ ) : t(t_),max_link_size(32) {}

  int lower_link( const Vertex & i, Vertex link[] ) const 
  { 
    int nlink = t.upper_link(i,link);
    return nlink;
  }

  Vertex walk_back( const Vertex & v ) const 
  { 
    assert(t.reachable_max[v] != uint32_t(-1));
    return t.reachable_max[v]; 
  }
};


ContourTree::ContourTree 
( uint8_t *image, uint32_t ncols, uint32_t nrows, uint32_t nstacks )
: tl(image,ncols,nrows,nstacks)
{
  voxels = image;
  img_size[0] = ncols; 
  img_size[1] = nrows; 
  img_size[2] = nstacks; 
  nvoxels = ncols*nrows*nstacks;
}

void ContourTree::build ()
{

  tl.find_critical_voxels();
  tl.find_saddles();
  tl.mark_reachable_maxes();
  tl.mark_reachable_mins();
  cout << "found " << tl.maxes.size() << " maxes, " << tl.mins.size() << " mins, "
       << tl.voxel_saddles.size() << " voxel-saddles and " 
       << tl.nonvoxel_saddles.size() << " nonvoxel-saddles" << endl;

  vector<uint32_t> verts;
  append_vector(verts,tl.mins);
  append_vector(verts,tl.voxel_saddles);
  append_vector(verts,tl.maxes);

  tl.maxes.clear();
  tl.mins.clear();
  tl.voxel_saddles.clear();

  if (!tl.place_holders.empty()) {
    cout << "merging in placeholders" << endl;
    append_vector(verts,tl.place_holders);
    sort(verts.begin(),verts.end());
    vector<uint32_t>::iterator e = unique(verts.begin(),verts.end());
    verts.resize( distance(verts.begin(),e) );
    tl.place_holders.clear();
  }

  
  cout << "sorting" << endl;
  #pragma omp parallel sections
  {
    #pragma omp section
    sort( verts.begin(),verts.end(), Trilinear<uint8_t>::compare_voxels(tl) );
    #pragma omp section
    sort( tl.nonvoxel_saddles.begin(), tl.nonvoxel_saddles.end(), 
          Trilinear<uint8_t>::compare_saddles(tl) );
  }

  cout << "merge voxel and saddle lists" << endl;
  tl.merge_sorted_voxels_and_saddles(verts,tl.nonvoxel_saddles);
  tl.nonvoxel_saddles.clear();


  typedef Trilinear<uint8_t>::ComponentMap ComponentMap;
  ComponentMap join_comps(tl,Join), split_comps(tl,Split);
  
  SweepComponent<uint32_t> *join_root, *split_root;

  cout << "sweep" << endl;
  #pragma omp parallel sections 
  {
    #pragma omp section
    {
    JoinSweepClosure join_clos(tl);
    join_root = minimal_sweep(Join,verts.begin(),verts.end(),
                              join_clos,join_comps);
    }
    #pragma omp section
    {
      SplitSweepClosure split_clos(tl);
      split_root= minimal_sweep(Split,verts.rbegin(),verts.rend(),
                                split_clos,split_comps);
    }
  }

  delete[] tl.reachable_max;
  delete[] tl.reachable_min;
  tl.reachable_max = tl.reachable_min = 0;

  remove_regular_points(verts.begin(),verts.end(),join_comps,split_comps);
  verts.clear();

  cout << "merge" << endl;
  Tourtre::merge(join_root,split_root,join_comps,split_comps,node_map);

  cout << "branch decomp" << endl;
  Node *n = node_map.begin()->second;
  
  //get list of nodes
  get_nodes(n,back_inserter(nodes));

  //mark node ids
  for ( uint32_t i=0; i<nodes.size(); ++i ) nodes[i]->id = i;

  Node *root_branch = greedy_branch_decomposition(nodes);
  assert(root_branch);

}
