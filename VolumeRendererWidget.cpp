

#include <algorithm>

#include <QMouseEvent>

#include "VolumeRendererWidget.hpp"
#include "Geom/Geom.hpp"
#include "Geom/GeomGL.hpp"

using namespace std;
using namespace Geom;

VolumeRendererWidget::VolumeRendererWidget 
( ContourTreeVolumeRenderer *ctvr_,
  QWidget *parent 
) : 
  QGLWidget(parent),
  ctvr(ctvr_)
{
}

void
VolumeRendererWidget::initializeGL()
{
  glClearColor(1,1,1,1);
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  gluPerspective(40.0, 1.0, 1.0, 10.0);

  const uint32_t *size = ctvr->vol_size;

  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
  glTranslatef(-0.5, -0.5, -5.0);
  double scale = 1.0 / max(size[0], max(size[1],size[2]));
  glScaled(scale,scale,scale);

  ctvr->init_gl();
}


void
VolumeRendererWidget::resizeGL(int w, int h)
{
  float ar = float(h) / float(w);
  glViewport(0, 0, w, h);
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  glOrtho(-1,1,-ar,ar,1/256.0,2048);
}


void
VolumeRendererWidget::paintGL()
{
  uint32_t *size = ctvr->vol_size;
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

  GLint vp[4]; 
  glGetIntegerv(GL_VIEWPORT,vp);

  GLdouble mv[16],proj[16];
  glGetDoublev(GL_MODELVIEW_MATRIX,mv);
  glGetDoublev(GL_PROJECTION_MATRIX,proj);

  double bl[3],br[3],tl[3],tr[3]; 
    //bottom left, bottom right, top left, top right
  gluUnProject(vp[0],vp[1],0.9,mv,proj,vp, bl,bl+1,bl+2 );
  gluUnProject(vp[2],vp[1],0.9,mv,proj,vp, br,br+1,br+2 );
  gluUnProject(vp[0],vp[3],0.9,mv,proj,vp, tl,tl+1,tl+2 );
  gluUnProject(vp[2],vp[3],0.9,mv,proj,vp, tr,tr+1,tr+2 );

  glBegin(GL_QUADS);
    glColor3f(0.5,0.5,0.5);
    glVertex3dv(bl);
    glVertex3dv(br);
    glColor3f(0.8,0.8,0.8);
    glVertex3dv(tr);
    glVertex3dv(tl);
  glEnd();

  glBegin(GL_LINES);
    glColor3f(1.0,0,0);
    glVertex3d(      0,  0,  0);
    glVertex3d(size[0],  0,  0);
    glColor3f(0,1.0,0);
    glVertex3d(  0,        0,  0);
    glVertex3d(  0,  size[1],  0);
    glColor3f(0,0,1.0);
    glVertex3d(  0,  0,      0);
    glVertex3d(  0,  0, size[2]);
  glEnd();
  
  glPushAttrib(GL_ALL_ATTRIB_BITS);
  glPushClientAttrib(GL_CLIENT_ALL_ATTRIB_BITS);
  ctvr->enable_gl();
  draw_box_faces();
  glPopClientAttrib();
  glPopAttrib();
  glUseProgram(0);

}


static inline
void plot( double x, double y, double z )
{
  glTexCoord3d(x,y,z);
  glVertex3d(x,y,z);
}

void VolumeRendererWidget::draw_box_faces()
{
  double fp[3],bp[3],vv[3]; //front point, back point, view vector

  GLint vp[4]; 
  glGetIntegerv(GL_VIEWPORT,vp);
  GLdouble mv[16],proj[16];
  glGetDoublev(GL_MODELVIEW_MATRIX,mv);
  glGetDoublev(GL_PROJECTION_MATRIX,proj);

  gluUnProject(vp[2]/2.0,vp[3]/2.0,0,mv,proj,vp,fp,fp+1,fp+2);
  gluUnProject(vp[2]/2.0,vp[3]/2.0,1,mv,proj,vp,bp,bp+1,bp+2);
  for (int i=0; i<3; ++i ) vv[i] = bp[i]-fp[i];

  const uint32_t* s = ctvr->vol_size;

  glBegin(GL_QUADS);
    glColor3f(1,0,0);
    if (vv[0] > 0) {
      plot(0,  0 ,  0 );
      plot(0,  0 ,s[2]);
      plot(0,s[1],s[2]);
      plot(0,s[1],  0 );
    } else {
      plot(s[0],  0 ,  0 );
      plot(s[0],  0 ,s[2]);
      plot(s[0],s[1],s[2]);
      plot(s[0],s[1],  0 );
    }

    glColor3f(0,1,0);
    if (vv[1] > 0) {
      plot(  0 ,0,  0 );
      plot(s[0],0,  0 );
      plot(s[0],0,s[2]);
      plot(  0 ,0,s[2]);
    } else {
      plot(  0 , s[1],   0 );
      plot(s[0], s[1],   0 );
      plot(s[0], s[1], s[2]);
      plot(  0 , s[1], s[2]);
    }

    glColor3f(0,0,1);
    if (vv[2] > 0) {
      plot(  0 ,  0 ,0);
      plot(s[0],  0 ,0);
      plot(s[0],s[1],0);
      plot(  0 ,s[1],0);
    } else {
      plot(  0 ,  0 ,s[2]);
      plot(s[0],  0 ,s[2]);
      plot(s[0],s[1],s[2]);
      plot(  0 ,s[1],s[2]);
    }
  glEnd();
}




void VolumeRendererWidget::mousePressEvent( QMouseEvent *event )
{
  mouse_x = event->x();
  mouse_y = height()-event->y();
  setCursor(QCursor(Qt::ClosedHandCursor));
}

void VolumeRendererWidget::mouseMoveEvent( QMouseEvent *event )
{
  makeCurrent();
  double w=width(),h=height();
  int event_y = height()-event->y();
  double
    x0 = (mouse_x/w)*2-1,
    y0 = (mouse_y/h)*2-1,
    x1 = (event->x()/w)*2-1 ,
    y1 = (event_y   /h)*2-1 ;
  glMatrixMode(GL_MODELVIEW);

  const uint32_t *size = ctvr->vol_size;

  pair< Vec3d, double> ang_ax = trackball( vec2(x0,y0) , vec2(x1,y1) );
  
  Vec3d v = inverse( Matrix33d::from_upper_left( 
               get_gl_matrix(GL_MODELVIEW_MATRIX))) * ang_ax.first;
  glMatrixMode(GL_MODELVIEW); 
  glTranslated(size[0]/2.0, size[1]/2.0, size[2]/2.0 );
  glRotated( ang_ax.second*180.0/M_PI, v[0],v[1],v[2]);
  glTranslated(-(size[0]/2.0), -(size[1]/2.0), -(size[2]/2.0) );

  update();
  mouse_x = event->x();
  mouse_y = event_y;
}

void VolumeRendererWidget::wheelEvent( QWheelEvent *event )
{
  double zoom = event->delta() > 0 ? 1.1 : 0.9;
  const uint32_t *size = ctvr->vol_size;
  glTranslatef(size[0]/2.0, size[1]/2.0, size[2]/2.0 );
  glScaled(zoom,zoom,zoom);
  glTranslatef(-size[0]/2.0, -size[1]/2.0, -size[2]/2.0 );
  update();
}

