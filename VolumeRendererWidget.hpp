
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
    void draw_box_faces();
    int mouse_x, mouse_y;

    int sw_rendering; //0=not using sw renderer, 1=rendering, 2=done
    RGBA8 *sw_img; //sw render image
    int32_t sw_win_left, sw_win_bottom, sw_win_width, sw_win_height;
      //where is the sw rendering window located?
    int32_t sw_row; //what was the last row we rendered?
    bool dragging_sw_win; //is the user currently selecting the sw render area?

    void start_sw_rendering();
    void continue_sw_rendering();

};

#endif
