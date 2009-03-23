
#ifndef VOLUME_RENDERER_WIDGET_HPP_INCLUDED
#define VOLUME_RENDERER_WIDGET_HPP_INCLUDED

#define GL_GLEXT_PROTOTYPES
#include <QGLWidget>

#include "ContourTreeVolumeRenderer.hpp"


class VolumeRendererWidget : public QGLWidget 
{
  Q_OBJECT

  ContourTreeVolumeRenderer *ctvr;
	
  public:
    VolumeRendererWidget(
      ContourTreeVolumeRenderer *ctvr, 
      QWidget *parent = 0 );
          
    void mousePressEvent ( QMouseEvent* );
    void mouseMoveEvent ( QMouseEvent* );
    void wheelEvent ( QWheelEvent* );

    void mouseReleaseEvent ( QMouseEvent* ) {}
    void mouseDoubleClickEvent ( QMouseEvent* ) {}
    void keyPressEvent ( QKeyEvent* ) {}
    void keyReleaseEvent( QKeyEvent* ) {}
    
    void paintGL();
    void resizeGL(int w, int h);
    void initializeGL();

  private:
    void draw_box_faces();

    int mouse_x, mouse_y;
};

#endif
