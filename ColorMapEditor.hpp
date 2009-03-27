#ifndef COLOR_MAP_EDITOR_HPP
#define COLOR_MAP_EDITOR_HPP


#include <stdint.h>

#include <QGLWidget>

#include "Color.hpp"


class ColorMapEditor : public QGLWidget   
{
  Q_OBJECT

  public:

  static const uint32_t resolution = 256;

  enum Channel { Hue, Lum, Sat, Alpha };
  enum Mode { Curve, Line };

  ColorMapEditor( QWidget *parent=0 );
  
  void initializeGL();
  void resizeGL(int,int);
  void paintGL();

  Mode mode;
  Channel channel;

  RGBA8 rgba[resolution]; 
  HLSAf hlsa[resolution];

  void update_colors ();
  void stroke( Channel, float, float, float, float );

  void mousePressEvent( QMouseEvent *event );
  void mouseMoveEvent( QMouseEvent *event );
  void mouseReleaseEvent( QMouseEvent *event );

  private:
  void draw_channel(Channel chan);
  void write_hlsa_to_rgba(int ibegin, int iend );
  void write_rgba_to_hlsa( int ibegin, int iend );
  
  double prev_x,prev_y;
  GLuint bg_tex, color_tex;
  bool colors_need_update;


  Q_SIGNALS:
    void edited (RGBA8*);
};



#endif

