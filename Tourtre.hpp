#ifndef TOURTRE_HPP_INCLUDED
#define TOURTRE_HPP_INCLUDED

#include <cassert>
#include <vector>
#include <map>
#include <deque>
#include <iostream>
#include <cstdlib>

namespace Tourtre {

//for type inference
template < typename T > struct pointer_target {};
template < typename T > struct pointer_target<T*> { typedef T type; };

enum SweepType { Join, Split };

template <typename Vertex_>
struct SweepComponent 
{
  typedef Vertex_ Vertex;
  SweepType type;
  Vertex birth;
  SweepComponent *succ, *pred; //successor and predecessor
  SweepComponent *next_pred;   //null-terminated linked list of predecessors
  SweepComponent *uf;          //union-find link

  SweepComponent(SweepType t) 
    : type(t),succ(0),pred(0),next_pred(0),uf(this) {}

  SweepComponent(SweepType t, const Vertex & b) 
    : type(t),birth(b),succ(0),pred(0),next_pred(0),uf(this) {}

  SweepComponent(SweepType t, const Vertex & b,const Vertex & d) 
    : type(t),birth(b),succ(0),pred(0),next_pred(0),uf(this) {}
  
  void unite( SweepComponent* c ) 
  { uf = c->uf; }

  SweepComponent* find()
  { 
    // find parent 
    SweepComponent* c=uf;
    while(c!=c->uf) { c=c->uf; }
    // path compression 
    SweepComponent *s=uf, *t;
    while(s!=c) { t=s->uf; s->uf=c; s=t; }
    uf = c;
    return c;
  }
  
  void add_pred( SweepComponent* c ) 
  { 
    c->succ = this;
    c->next_pred = pred;
    pred = c;
  }

  void remove_pred( SweepComponent* c )
  { 
    assert(pred);
    SweepComponent *p = pred;
    if (p == c) {
      pred = c->next_pred;
      c->next_pred = 0;
    } else {
      while(p && p->next_pred!=c) p=p->next_pred;
      assert(p && p->next_pred==c);
      p->next_pred = c->next_pred;
      c->next_pred = 0;
    }
  }

  void prune() 
  { 
    assert(pred == 0);
    if (succ) succ->remove_pred(this);
    succ = 0;
  }

  bool is_regular() { return pred && pred->next_pred == 0; }
  bool is_leaf() { return !pred; }

  SweepComponent* merge_with_succ()
  { 
    assert( succ && succ->pred == this );
    assert( !next_pred );
    SweepComponent *s = succ;
    SweepComponent *ss = s->succ;
    if (ss) {
      ss->remove_pred(s);
      ss->add_pred(this);
    }
    succ = s->succ;
    s->succ = s->pred = 0;
    s->next_pred = 0;
    unite(this);
    return s;
  }
};


template <typename Vertex_>
struct Node;

template <typename Vertex_>
struct Arc
{
  typedef Vertex_ Vertex;
  typedef Tourtre::Node<Vertex> Node;
  Node *lo,*hi;
  uint32_t id;
  Arc *next_up,*next_down; //null-terminated linked list
  uint32_t branch; //two arcs are on the same branch if this field is the same
  Arc() : lo(0),hi(0),id(-1),next_up(0),next_down(0),branch(-1) {}
};

template <typename Vertex_>
struct Node
{
  typedef Vertex_ Vertex;
  typedef Tourtre::Arc<Vertex> Arc;
  Vertex vertex;
  Arc *up,*down;
  uint32_t id; 
  Node() : up(0),down(0) {}

  Node(Vertex v) : vertex(v),up(0),down(0) {}

  void add_up_arc( Arc *a ) 
  {
    a->next_up = up;
    up = a;
    a->lo = this;
  }

  void add_down_arc( Arc *a ) 
  { 
    a->next_down = down; 
    down = a; 
    a->hi = this;
  }

  void remove_up_arc( Arc *a )
  {
    if ( up == a ) {
      up = a->next_up; 
    } else {
      Arc *u = up;
      while( u->next_up && u->next_up!=a ) u=u->next_up;
      assert(u->next_up && u->next_up==a);
      u->next_up = u->next_up->next_up;
    }
  }

  void remove_down_arc( Arc *a )
  {
    if ( down == a ) {
      down = a->next_down; 
    } else {
      Arc *d = down;
      while( d->next_down && d->next_down!=a ) d=d->next_down;
      assert(d->next_down && d->next_down==a);
      d->next_down = d->next_down->next_down;
    }
  }
  

  bool is_max() const { return !up ;}
  bool is_min() const { return !down ;}

  int up_degree() const 
  {
    int deg=0;
    for (Arc *u=up; u; u=u->next_up) ++deg;
    return deg;
  }

  int down_degree() const 
  {
    int deg=0;
    for(Arc *d = down; d; d=d->next_down) ++deg;
    return deg;
  }

};


template <typename Vertex, typename OutputItr >
void get_nodes( Node<Vertex> *n,  OutputItr out )
{
  typedef std::pair<Node<Vertex>*,Node<Vertex>*> NodePair;
  std::vector< NodePair > stack(1,std::make_pair(n,static_cast<Node<Vertex>*>(0)));
  while(!stack.empty()) {
    NodePair p = stack.back();
    stack.pop_back();
    *out = p.first;
    ++out;
    for ( Arc<Vertex> *a = p.first->up; a; a=a->next_up ) {
      if ( a->hi != p.second ) stack.push_back( std::make_pair(a->hi,p.first)); 
    }
    for ( Arc<Vertex> *a = p.first->down; a; a=a->next_down ) {
      if ( a->lo != p.second ) stack.push_back( std::make_pair(a->lo,p.first)); 
    }
  }
}





// full_sweep takes a Closure argument that must implement the following:
//
//    int link( const Vertex &, Vertex[] )
//      yields the immediate neighbors of a potential critical point (pcp)
//      that preceed it in the sweep
//
//    int max_link_size
//      the minimum size for the 2nd argument to lower_link
//
//  The output of sweep is its final argument, which can be any dictionary-like
//  container implementing the [] operator, taking a Vertex as an index.
//  However this lookup function will be hammered so it should probably be an
//  array.  The sweep components will be stored in that container, keyed by
//  their birth vertices. The container must return a null pointer on the first
//  access of a vertex. (This is the default behavior for stl map-link
//  containers, and for an array or vector it suffices to initialize it with
//  nulls)
//
//  Returns the final component 


template <typename InputItr, typename Closure, typename ComponentMap>
SweepComponent<typename InputItr::value_type>*
full_sweep( 
  SweepType type, //is this a join or split sweep?
  InputItr begin, //the sorted range of vertices
  InputItr end,          
  Closure & closure, //yields lower links and reachable extrema
  ComponentMap & comps ) //output
{
  typedef typename InputItr::value_type Vertex;
  SweepComponent<Vertex> *icomp = 0;
  Vertex last;
  for ( InputItr itr=begin; itr!=end; ++itr ) {
    Vertex i = *itr;
    last = i;
    Vertex link[closure.max_link_size];
    int nlink = closure.link(i,link);
    int num_comps_here=0;
    icomp = 0;
    for ( int l=0; l<nlink; ++l ) {
      Vertex j = link[l];
      SweepComponent<Vertex> *jcomp = comps[j];
      if (jcomp) {
        jcomp = jcomp->find();
        if ( icomp != jcomp ) {
          if (num_comps_here == 0) {
            ++num_comps_here;
            icomp = jcomp;
          } else if (num_comps_here == 1) {
            SweepComponent<Vertex> *new_comp = new SweepComponent<Vertex>(type,i);
            new_comp->add_pred(icomp);
            new_comp->add_pred(jcomp);
            icomp->unite(new_comp);
            jcomp->unite(new_comp);
            icomp = new_comp;
            num_comps_here++;
          } else {
            jcomp->unite(icomp);
            icomp->add_pred(jcomp);
          }
        }
      }
      comps[i] = icomp ;
    }

    if (num_comps_here == 0) {
      icomp = new SweepComponent<Vertex>(type,i);
      comps[i] = icomp;
    }
  }

  if (icomp) { //if we did anything at all
    SweepComponent<Vertex> *inf = new SweepComponent<Vertex>(type,last);
    icomp = comps[last]->find(); 
    inf->add_pred(icomp);
    comps[inf->birth] = inf;
    return inf;
  } else {
    return 0; 
  }
  
}


// augment : ensure that the join and split trees containt the same set of
// nodes. This may change the split tree root; if so, the new split root is
// returned.

template <typename InputItr, typename Vertex, typename ComponentMap>
SweepComponent<Vertex>* 
augment( 
  InputItr begin, //the sorted range of vertices
  InputItr end,          
  ComponentMap & join_comps, 
  ComponentMap & split_comps,
  SweepComponent<Vertex>* split_root )
{
  ++begin;
  --end;
  for ( InputItr itr=begin; itr!=end; ++itr ) {
    Vertex i = *itr;
    SweepComponent<Vertex> *join=join_comps[i], *split=split_comps[i];
    if (join->birth==i && split->birth!=i) {
      SweepComponent<Vertex> *new_comp = new SweepComponent<Vertex>(Split,i);
      if (split->succ) {
          split->succ->remove_pred( split );
          split->succ->add_pred( new_comp );
      }
      new_comp->add_pred(split);
      if (split == split_root) split_root = new_comp;
      split_comps[new_comp->birth] = new_comp;
    } else if ( split->birth==i && join->birth!=i ) {
      SweepComponent<Vertex> *new_comp = new SweepComponent<Vertex>(Join,join->birth);
      join->birth = i;
      while( join->pred ) {
        SweepComponent<Vertex> *p = join->pred;
        join->remove_pred(p);
        new_comp->add_pred(p);
      }
      join->add_pred(new_comp);
      join_comps[join->birth] = join;
      join_comps[new_comp->birth] = new_comp;
    }
  }
  return split_root;
}






// minimal_sweep takes a Closure argument that must implement the following:
//
//    int lower_link( const Vertex &, Vertex[] )
//      yields the immediate neighbors of a potential critical point (pcp)
//      that preceed it in the sweep
//
//    int max_link_size
//      the minimum size for the 2nd argument to lower_link
//
//    Vertex walk_back( const Vertex & )
//      yields an extrema that preceeds the vertex in the sweep and 
//      that is reachable from the vertex by a monotone path
//
//  The output of sweep is its final argument, which can be any dictionary-like
//  container implementing the [] operator, taking a Vertex as an index.  The
//  sweep components will be stored in that container, keyed by their birth
//  vertices.
//
//  Returns the final component 


template <typename InputItr, typename Closure, typename ComponentMap>
SweepComponent<typename InputItr::value_type>*
minimal_sweep( 
  SweepType type, //is this a join or split sweep?
  InputItr begin, //the sorted range of vertices
  InputItr end,          
  Closure & closure, //yields lower links and reachable extrema
  ComponentMap & comps ) //output
{
  typedef typename InputItr::value_type Vertex;
  SweepComponent<Vertex> *icomp = 0;
  for ( InputItr itr=begin; itr!=end; ++itr ) {
    Vertex i = *itr;
    Vertex link[closure.max_link_size];
    int nlink = closure.lower_link(i,link);

    icomp = new SweepComponent<Vertex>(type,i);
    comps[i] = icomp;
    
    //std::cout << "i = " << i << std::endl;
    for ( int l=0; l<nlink; ++l ) {
      Vertex j = closure.walk_back(link[l]);
      //std::cout << "j = " << j << std::endl;
      SweepComponent<Vertex> *jcomp = comps[j];
      if ( !jcomp ) { 
        std::cout << "created componet at " << j << " out of order" << std::endl;
        abort();
      }
      jcomp = jcomp->find();
      if ( icomp != jcomp ) {
        jcomp->unite(icomp);
        icomp->add_pred(jcomp);
      }
    }
  }

  return icomp;
}


template <typename Vertex, typename LeafQ >
void queue_leaves( SweepComponent<Vertex> *root, LeafQ & leafq )
{
  std::vector<SweepComponent<Vertex>*> stack(1,root);
  while(!stack.empty()) {
    SweepComponent<Vertex> *c = stack.back();
    stack.pop_back(); 
    if ( c->is_leaf() ) leafq.push_back(c);
    else {
      for (SweepComponent<Vertex> *p = c->pred; p; p=p->next_pred)
        stack.push_back(p);
    }
  }
}



template <typename Vertex >
void remove_regular_points( 
  std::map<Vertex,SweepComponent<Vertex>*> & join_map,
  std::map<Vertex,SweepComponent<Vertex>*> & split_map )
{
  typedef std::map<Vertex,SweepComponent<Vertex>*> ComponentMap;

  typename ComponentMap::iterator 
    jitr=join_map.begin(), sitr=split_map.begin();
  
  while( jitr != join_map.end() ) {
    assert( jitr->first == sitr->first );
    if ( jitr->second->is_regular() && sitr->second->is_regular() ) {
      typename ComponentMap::iterator j=jitr, s=sitr;
      ++jitr,++sitr;
      
      SweepComponent<Vertex> *dead;
      dead = j->second->pred->merge_with_succ(); 
      delete dead;
      dead = s->second->pred->merge_with_succ(); 
      delete dead;

      join_map.erase(j);
      split_map.erase(s);
    
    } else {
      ++jitr,++sitr;
    }
  }
}


template <typename InputItr, typename ComponentMap  >
void remove_regular_points( 
  InputItr begin, InputItr end,
  ComponentMap & join_map,
  ComponentMap & split_map )
{
  typedef typename InputItr::value_type Vertex;
  for ( InputItr i=begin; i!=end; ++i ) {
    SweepComponent<Vertex> *j=join_map[*i], *s=split_map[*i];

    assert( j->birth==s->birth );
    if ( j->is_regular() && s->is_regular() ) {
      
      SweepComponent<Vertex> *dead;
      dead = j->pred->merge_with_succ(); 
      delete dead;
      dead = s->pred->merge_with_succ(); 
      delete dead;

      join_map[*i] = split_map[*i] = 0;
    }
  }
}












//Merge the join and split trees. The last argument is a value indicating
//a null vertex somehow, such as -1 if Vertex is an integral type or NULL
//if its a pointer.
//
// The output is put in the node map
//
// NB: there is a function named merge in the std namespace. If you
// get strange template errors, try disambiguating

template <typename Vertex, typename ComponentMap, typename NodeMap>
void
merge( 
  SweepComponent<Vertex>* join_root,
  SweepComponent<Vertex>* split_root,
  ComponentMap & join_map, 
  ComponentMap & split_map,
  NodeMap & node_map )
{
  typedef std::deque<SweepComponent<Vertex>*> LeafQueue;
  LeafQueue leafq;
  queue_leaves(join_root,leafq);
  queue_leaves(split_root,leafq);

  ComponentMap *other_map;
  Arc<Vertex> *arc = 0;
  for(;;) {
    assert(!leafq.empty());
    SweepComponent<Vertex> *leaf = leafq.front();
    leafq.pop_front();

    if (!leaf->succ) { // all done 
      break;
    }
    Node<Vertex> *lo,*hi;
    // which tree is this comp from? 
    if ( leaf->type == Join ) {
      // comp is join component 
      other_map = &split_map;
      typename NodeMap::iterator 
        lo_itr = node_map.find(leaf->birth),
        hi_itr = node_map.find(leaf->succ->birth);
      if (lo_itr == node_map.end()) {
        lo = new Node<Vertex>(leaf->birth);
        node_map.insert( std::make_pair(leaf->birth,lo) );
      } else {
        lo = lo_itr->second; 
      }
      if (hi_itr == node_map.end()) {
        hi = new Node<Vertex>(leaf->succ->birth);
        node_map.insert(std::make_pair(leaf->succ->birth,hi));
      } else {
        hi = hi_itr->second; 
      }
    } else { // split component 
      other_map = &join_map;
      typename NodeMap::iterator 
        hi_itr = node_map.find(leaf->birth),
        lo_itr = node_map.find(leaf->succ->birth);
      if (hi_itr == node_map.end()) {
        hi = new Node<Vertex>(leaf->birth);
        node_map.insert(std::make_pair(leaf->birth,hi));
      } else {
        hi = hi_itr->second; 
      }
      if (lo_itr == node_map.end()) {
        lo = new Node<Vertex>(leaf->succ->birth);
        node_map.insert(std::make_pair(leaf->succ->birth,lo));
      } else {
        lo = lo_itr->second; 
      }
    }
    // create arc 
    arc = new Arc<Vertex>;
    lo->add_up_arc(arc);
    hi->add_down_arc(arc);
    // remove leaf 
    SweepComponent<Vertex> *succ = leaf->succ;
    assert(succ);
    leaf->prune();
    // remove leaf's counterpart in other tree 
    SweepComponent<Vertex> 
      *other = (*other_map)[leaf->birth],
      *other_succ = (*other_map)[succ->birth];
    assert(other&&other_succ);
    assert(other->is_regular()) ;
    
    SweepComponent<Vertex>*
     dead = other->pred->merge_with_succ();

    if ( succ->is_leaf() && other_succ->is_regular() ) {
      leafq.push_back(succ);
    } else if ( succ->is_regular() && other_succ->is_leaf() ) {
      leafq.push_back(other_succ);
    }

    delete dead;
    delete leaf;
  }
} 





template <typename Vertex>
void mark_farthest_maxes( std::vector<Node<Vertex>*> & nodes, std::vector<int> & farthest )
{
  typedef std::pair<Node<Vertex>*,int> Pair;
  std::deque<Pair> queue;

  for (size_t i=0; i<nodes.size(); ++i ) {
    Node<Vertex> *n = nodes[i];
    int deg = n->up_degree();
    assert(n->id < farthest.size());
    farthest[n->id] = -(deg-1);
    if (deg==0) { 
      queue.push_back( std::make_pair(n,0) );
    }
  }
  while(!queue.empty()) {
    Pair p = queue.front();
    queue.pop_front();
    if ( farthest[p.first->id] < 0 ) {
      farthest[p.first->id]++; 
    } else {
      farthest[p.first->id] = p.second;
      for ( Arc<Vertex>* d=p.first->down; d; d=d->next_down ) {
        queue.push_back( std::make_pair(d->lo,p.second+1) );
      }
    } 
  }
}

template <typename Vertex>
void mark_farthest_mins( std::vector<Node<Vertex>*> & nodes, std::vector<int> & farthest )
{
  typedef std::pair<Node<Vertex>*,int> Pair;
  std::deque<Pair> queue;

  for (size_t i=0; i<nodes.size(); ++i ) {
    Node<Vertex> *n = nodes[i];
    int deg = n->down_degree();
    assert(n->id < farthest.size());
    farthest[n->id] = -(deg-1);
    if (deg==0) { 
      queue.push_back( std::make_pair(n,0) );
    }
  }
  while(!queue.empty()) {
    Pair p = queue.front();
    queue.pop_front();
    if ( farthest[p.first->id] < 0 ) {
      farthest[p.first->id]++; 
    } else {
      farthest[p.first->id] = p.second;
      for ( Arc<Vertex>* d=p.first->up; d; d=d->next_up ) {
        queue.push_back( std::make_pair(d->hi,p.second+1) );
      }
    } 
  }
}




// Form a branch decomposition where each branch contains as many saddles as
// possible. This requires that the nodes be marked with contiguous ids, e.g.
// after using get_nodes.
//
// Returns the number of branches and the root node. The 1st argument is input,
// a list of nodes. The 2nd argument is output, the root arc of each branch

enum Direction { Up,Down };

template <typename Vertex>
void
greedy_branch_decomposition
( std::vector<Node<Vertex>*> & nodes,  
  std::vector<Arc<Vertex>*> & branches ) 
{
  std::vector<int> farthest_max(nodes.size()), farthest_min(nodes.size());
  mark_farthest_maxes(nodes,farthest_max);
  mark_farthest_mins(nodes,farthest_min);

  //find the max (min) that is furthest from a min (max)
  int best_dist=-1; 
  Node<Vertex> *best_node=0;
  for (size_t i=0; i<nodes.size(); ++i) 
  {
    Node<Vertex> *n = nodes[i];
    if (n->is_max() && farthest_min[n->id] > best_dist ) {
      best_dist = farthest_min[n->id]; 
      best_node = n;
    } else if (n->is_min() && farthest_max[n->id] > best_dist ) {
      best_dist = farthest_max[n->id]; 
      best_node = n;
    }
  }

  std::vector< std::pair<Arc<Vertex>*,Direction> > stack;
  if ( best_node->is_max() ) 
    stack.push_back( std::make_pair(best_node->down,Down) );
  else
    stack.push_back( std::make_pair(best_node->up,Up) );

  //for each root, walk out to the farthest max (min) and call those
  //arcs the branch
  while(!stack.empty()) {
    Arc<Vertex> *first = stack.back().first;  
    Direction dir = stack.back().second;
    stack.pop_back();
    uint32_t branch_id = branches.size();
    branches.push_back(first);
    Arc<Vertex> *arc = first;
    if ( dir == Up ) {
      for(;;) {
        arc->branch = branch_id; 
        Arc<Vertex> *next = arc->hi->up;
        if (!next) break;
        //the next arc along the branch is the one that leads to the farthest max 
        for ( Arc<Vertex> *a=next->next_up; a; a=a->next_up ) {
          if ( farthest_max[a->hi->id] > farthest_max[next->hi->id] )
            next = a;
        }
        //push all the rest of the arcs here onto the stack 
        for ( Arc<Vertex> *a=arc->hi->up; a; a=a->next_up ) 
          if ( a!=next ) {
            stack.push_back( std::make_pair(a,Up) );
          }
        for ( Arc<Vertex> *a=arc->hi->down; a; a=a->next_down ) {
          if ( a!=arc ) {
            stack.push_back( std::make_pair(a,Down) );
          }
        }
        arc = next;
      }
    } else { // dir == Down
      for(;;) {
        arc->branch = branch_id; 
        Arc<Vertex> *next = arc->lo->down;
        if (!next) break;
        for ( Arc<Vertex> *a=next->next_down; a; a=a->next_down ) {
          if ( farthest_min[a->lo->id] > farthest_min[next->lo->id] )
            next = a;
        }
        for ( Arc<Vertex> *a=arc->lo->up; a; a=a->next_up ) {
          if ( a!= arc ) {
            stack.push_back( std::make_pair(a,Up) );
          }
        }
        for ( Arc<Vertex> *a=arc->lo->down; a; a=a->next_down )
          if ( a!=next ) {
            stack.push_back( std::make_pair(a,Down) );
          }
        arc = next;
      }
    }
  }
}



} //namespace Tourtre

#endif

