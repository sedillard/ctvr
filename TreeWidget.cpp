
#include <algorithm>

#include <QMouseEvent>

#include "TreeWidget.hpp"
#include "Geom/Geom.hpp"
#include "Geom/GeomGL.hpp"

#include <boost/foreach.hpp>

using namespace std;
using namespace Geom;
using namespace boost;

#define foreach BOOST_FOREACH

TreeWidget::TreeWidget 
( ContourTreeVolumeRenderer *ctvr_,
  QWidget *parent 
) : 
  QGLWidget(parent),
  ct(&ctvr_->ct),
  tl(&ctvr_->tl),
  ctvr(ctvr_)
{
  xpos = vector<float>(ct->nodes.size(),0.5);
  randomized_tree_layout();
}

void
TreeWidget::initializeGL()
{
  glClearColor(1,1,1,1);
  glEnable(GL_LINE_SMOOTH);
  glEnable(GL_POINT_SMOOTH);
  glEnable(GL_BLEND);
  glBlendFunc(GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);
  
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  glOrtho(0,1,0,255,-1,1);

  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
  glTranslatef( 0.5, 128, 0 );
  glScalef( 0.9,0.9,0.9);
  glTranslatef( -0.5, -128, 0 );
}


void
TreeWidget::resizeGL(int w, int h)
{
  glViewport(0, 0, w, h);
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  glOrtho(0,1,0,255,-1,1);
}


void
TreeWidget::paintGL()
{
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  draw_tree();
  
}




void TreeWidget::mousePressEvent( QMouseEvent *event )
{
  mouse_x = event->x();
  mouse_y = height()-event->y();
  update();
}

void TreeWidget::mouseReleaseEvent( QMouseEvent * )
{
}

void TreeWidget::mouseMoveEvent( QMouseEvent *event )
{
  int event_y = height()-event->y();
  mouse_x = event->x();
  mouse_y = event_y;
  update();
}

void TreeWidget::wheelEvent( QWheelEvent * )
{
  update();
}


void TreeWidget::draw_tree()
{
  glLineWidth(1);
  glBegin(GL_LINES);
  glColor4f(0,0,0,0.5);
  foreach( Node *n, ct->nodes ) {
    for ( Arc *a=n->up; a; a=a->next_up ) {
      glVertex2d( xpos[a->hi->id], tl->value(a->hi->vertex) );
      glVertex2d( xpos[a->lo->id], tl->value(a->lo->vertex) );
    }
  }
  glEnd();

  glPointSize(10);
  glBegin(GL_POINTS);
  foreach (Node *n, ct->nodes ) {
    if ( !ctvr->node_cluster.empty() &&
         ctvr->node_cluster[n->id] != uint32_t(-1) ) 
    {
      srand( ctvr->node_cluster[n->id] );
      double r = rand()/double(RAND_MAX);
      double g = rand()/double(RAND_MAX);
      double b = rand()/double(RAND_MAX);
      double m = sqrt(r*r+g*g+b*b);
      glColor4f(r/m,g/m,b/m,0.5);
      glVertex2d( xpos[n->id], tl->value(n->vertex) );
    } 
  }
  glEnd();
  
}


void TreeWidget::randomized_tree_layout()
{
  foreach( Node* n, ct->nodes ) xpos[n->id] = rand()/double(RAND_MAX);  

  vector<Node*> rnodes = ct->nodes;
  for ( int itr=0; itr<100; ++itr ) {
    random_shuffle(rnodes.begin(),rnodes.end());
    foreach( Node* n, rnodes ) {
      int nnbrs=0, nup=0, ndown=0;
      float avg=0;
      for ( Arc* a=n->up; a; a=a->next_up ) avg+=xpos[a->hi->id],++nnbrs,++nup;
      for ( Arc* a=n->down; a; a=a->next_down ) avg+=xpos[a->lo->id],++nnbrs,++ndown;
      avg /= nnbrs;
      xpos[n->id] = 0.5*(xpos[n->id]+avg);
      //if ( nup == 2 ) {
      //  xpos[n->up->hi->id] -= 0.02; 
      //  xpos[n->up->next_up->hi->id] += 0.02; 
      //}
      //if (ndown == 2 ) {
      //  xpos[n->down->hi->id] -= 0.02; 
      //  xpos[n->down->next_down->hi->id] += 0.02; 
      //}
    }
    
    float max_x = -HUGE_VAL, min_x = HUGE_VAL; 
    foreach( Node* n, ct->nodes ) {
      max_x = max(max_x,xpos[n->id]);
      min_x = min(min_x,xpos[n->id]);
    }
    
    foreach( Node* n, ct->nodes ) {
      xpos[n->id] -= min_x;
      xpos[n->id] /= max_x-min_x;
    }
  }
}

