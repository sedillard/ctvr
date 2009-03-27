

#ifndef TREE_WIDGET_HPP_INCLUDED
#define TREE_WIDGET_HPP_INCLUDED

#define GL_GLEXT_PROTOTYPES
#include <QGLWidget>

#include "Trilinear.hpp"
#include "ContourTree.hpp"
#include "ContourTreeVolumeRenderer.hpp"

#include <vector>


class TreeWidget : public QGLWidget 
{
  Q_OBJECT

  typedef ContourTree::Node Node;
  typedef ContourTree::Arc Arc;

  ContourTree *ct;
  Trilinear<uint8_t> *tl;
  ContourTreeVolumeRenderer *ctvr;

	
  public:
    TreeWidget(
      ContourTreeVolumeRenderer *ctvr, 
      QWidget *parent=0 );
          
    void mousePressEvent ( QMouseEvent* );
    void mouseMoveEvent ( QMouseEvent* );
    void wheelEvent ( QWheelEvent* );

    void mouseReleaseEvent ( QMouseEvent* );
    void mouseDoubleClickEvent ( QMouseEvent* ) {}
    void keyPressEvent ( QKeyEvent* ) {}
    void keyReleaseEvent( QKeyEvent* ) {}
    
    void paintGL();
    void resizeGL(int w, int h);
    void initializeGL();

  private:
    int mouse_x, mouse_y;
    std::vector<float> xpos; //x-position of nodes
    void draw_tree();
    void randomized_tree_layout();

};

#endif
