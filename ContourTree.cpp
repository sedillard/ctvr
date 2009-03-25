
#include "ContourTree.hpp"

#include <utility>

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

  //tl.reachable_max = tl.reachable_min = vector<uint32_t>();

  remove_regular_points(verts.begin(),verts.end(),join_comps,split_comps);
  verts = vector<uint32_t>();

  cout << "merge" << endl;
  Tourtre::merge(join_root,split_root,join_comps,split_comps,node_map);

  tl.join_comps = tl.split_comps = vector<SweepComponent<uint32_t>*>();

  Node *n = node_map.begin()->second;
  
  //get list of nodes
  get_nodes(n,back_inserter(nodes));

  //mark node ids
  for ( uint32_t i=0; i<nodes.size(); ++i ) nodes[i]->id = i;

  //mark arc ids
  uint32_t narcs=0;
  for ( uint32_t i=0; i<nodes.size(); ++i ) {
    for ( Arc* a=nodes[i]->up; a; a=a->next_up ) {
      arcs.push_back(a);
      a->id = ++narcs;
    }
  }
  

  for ( uint32_t i=0; i<nodes.size(); ++i ) {
    assert ( !(nodes[i]->is_max() || nodes[i]->is_min()) || nodes[i]->vertex < nvoxels );
  }

  cout << "branch decomp" << endl;
  greedy_branch_decomposition(nodes,branches);
}




bool ContourTree::branch_is_ascending( uint32_t b )
{
  Arc *first = branches[b];
  if ( !first->hi->up ) return true; //its a max
  if ( !first->lo->down ) return false; //its a min
  for ( Arc *a=first->hi->up; a; a=a->next_up ) 
    if ( a->branch==b ) return true; 
  for ( Arc *a=first->lo->down; a; a=a->next_down ) 
    if ( a->branch==b ) return true; 
  assert(0&&"branch_is_ascending: couldn't find the 2nd arc along the branch");
}


pair<ContourTree::Node*,ContourTree::Node*> 
ContourTree::branch_range( uint32_t b ) 
{
  Arc *first=branches[b], *arc=first, *next;
  if ( branch_is_ascending(b) ) {
    for(;;) {
      next=0;
      for ( Arc*a=arc->hi->up; a; a=a->next_up ) {
        if ( a->branch == b ) {
          next = a;
          break;
        }
      }
      if (!next) return make_pair(first->lo,arc->hi);
      else arc = next;
    }
  } else {
    for(;;) {
      next=0;
      for ( Arc*a=arc->lo->down; a; a=a->next_down ) {
        if ( a->branch == b ) {
          next = a;
          break;
        }
      }
      if (!next) return make_pair(first->hi,arc->lo);
      else arc = next;
    }
  }
  abort();
}


void ContourTree::prune_flat_arcs()
{
  deque<Node*> leafq;
  for ( uint32_t i=0; i<nodes.size(); ++i ) 
    if ( nodes[i]->is_max() || nodes[i]->is_min() ) 
      leafq.push_back(nodes[i]);

  Node *some_node=0;
  while(!leafq.empty()) {
    Node *n = leafq.front();
    leafq.pop_front();
    
    if (n->is_max()) {
      Arc *a = n->down;
      Node *o = a->lo;
      if ( o->up_degree() > 1 ) {
        if ( tl.value(n->vertex) == tl.value(o->vertex) ) {
          o->remove_up_arc(a);    
          delete a;
          delete n;
          if ( o->up_degree()==1 && o->down_degree()==1 ) {
            Arc *u = o->up, *d = o->down;
            Node *h = u->hi, *l = d->lo;
            a = new Arc;
            h->remove_down_arc(u);
            h->add_down_arc(a);
            l->remove_up_arc(d);
            l->add_up_arc(a);
            delete u;
            delete d;
            delete o;
            if (h->is_max()) leafq.push_back(h);
            if (l->is_min()) leafq.push_back(l);
            some_node = l;
          } else {
            some_node = o; 
          }
        }
      }
    } else {
      //TODO 
    }
  }

  nodes.clear();
  get_nodes(some_node,back_inserter(nodes));
}

