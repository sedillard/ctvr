
#include "ColorMapEditor.hpp"

#include <QMouseEvent>


ColorMapEditor::ColorMapEditor( QWidget *parent )
  : QGLWidget(parent)
{

  colors_need_update = true;
  channel = Alpha;
  mode = Curve;
  

  stroke( Hue,   0, 0.60, 1, 0);
  stroke( Lum,   0, 0.5, 1, 0.5);
  stroke( Sat,   0, 1, 1, 1);
  stroke( Alpha, 0, 0, 1, 1);
  
}


void ColorMapEditor::initializeGL()
{
    glDisable( GL_DEPTH_TEST );
    glEnable( GL_BLEND );
    glBlendFunc( GL_SRC_ALPHA , GL_ONE_MINUS_SRC_ALPHA );
    glEnable(GL_LINE_SMOOTH);
    glLineWidth(2);
    glClearColor(1,1,1,1);
    
    //create texture to display color output
    glGenTextures( 1, &color_tex );
    glBindTexture( GL_TEXTURE_1D, color_tex );
    glTexParameteri( GL_TEXTURE_1D, GL_TEXTURE_BASE_LEVEL, 0 );
    glTexParameteri( GL_TEXTURE_1D, GL_TEXTURE_MAX_LEVEL, 0 );
    glTexParameteri( GL_TEXTURE_1D, GL_TEXTURE_WRAP_S, GL_CLAMP );
    glTexParameteri( GL_TEXTURE_1D, GL_TEXTURE_MAG_FILTER, GL_NEAREST );
    glTexParameteri( GL_TEXTURE_1D, GL_TEXTURE_MIN_FILTER, GL_NEAREST );
    glTexEnvi( GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE );
    glPixelStorei( GL_UNPACK_ALIGNMENT, 1 );
    glPixelStorei( GL_UNPACK_SKIP_PIXELS, 0 );
    glPixelStorei( GL_UNPACK_ROW_LENGTH, 0 );
    glPixelStorei( GL_UNPACK_SKIP_ROWS, 0 );
    glPixelStorei( GL_UNPACK_IMAGE_HEIGHT, 0 );
    glPixelStorei( GL_UNPACK_SKIP_IMAGES, 0 );
    glTexImage1D( GL_TEXTURE_1D, 0, GL_RGBA8, resolution , 0, 
                  GL_RGBA, GL_UNSIGNED_BYTE, rgba );
    
}


void ColorMapEditor::resizeGL(int w , int h )
{
  glViewport(0,0,w,h); 
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  glOrtho(0,1,0,1,-1,1);
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
}


void ColorMapEditor::mousePressEvent( QMouseEvent *event )
{
    makeCurrent();
    double x,y,z;
    GLint vp[4];
    GLdouble mv[16],proj[16];
    glGetIntegerv(GL_VIEWPORT,vp);
    glGetDoublev(GL_MODELVIEW_MATRIX,mv);
    glGetDoublev(GL_PROJECTION_MATRIX,proj);
    gluUnProject(event->x(),height()-event->y(),0,mv,proj,vp,&x,&y,&z);
    prev_x = x;
    prev_y = y;
    if (mode == Curve) {
        stroke(channel, x,y, x,y);
        update();
    }
}
        
void ColorMapEditor:: mouseMoveEvent(QMouseEvent *event) 
{
  makeCurrent();
  double x0 = prev_x;
  double y0 = prev_y;
  double x1,y1,z;
  GLint vp[4];
  GLdouble mv[16],proj[16];
  glGetIntegerv(GL_VIEWPORT,vp);
  glGetDoublev(GL_MODELVIEW_MATRIX,mv);
  glGetDoublev(GL_PROJECTION_MATRIX,proj);
  gluUnProject(event->x(),height()-event->y(),0,mv,proj,vp,&x1,&y1,&z);
  stroke(channel,x0,y0,x1,y1);
  if (mode == Curve) {
    prev_x = x1;
    prev_y = y1;
  }
  update();
}

void ColorMapEditor:: mouseReleaseEvent( QMouseEvent* )
{
  Q_EMIT(edited(rgba));
}

void
ColorMapEditor::paintGL () 
{
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  if (colors_need_update) {
      update_colors();
      colors_need_update = false;
  }

  //draw scale
  int scale_delta = 8;
  //while( scale_delta / 255.0 * height() * zoom < 20 ) scale_delta++;
  
  glColor4f(0,0,0,0.5);
  
  //for(size_t i = 0; i < resolution; i += scale_delta) {
  //    glBegin(GL_LINES);
  //    glVertex2f( i/(float)resolution , 0);
  //    glVertex2f( i/(float)resolution , 1);
  //    glEnd();
  //}
  
  //glEnable(GL_LINE_SMOOTH);

  //draw function
  //double res = resolution;
  //glEnable(GL_TEXTURE_1D);
  //glBindTexture(GL_TEXTURE_1D,color_tex);
  //glBegin(GL_QUADS);
  //    glTexCoord1d( 0.5 / (res-1.0)  );
  //    glVertex2d  ( 0, 0 );
  //    glVertex2d  ( 0, 1 );
  //    
  //    glTexCoord1d( 0.5 / (res-1.0) + (res-1.0)/res );
  //    glVertex2d  ( 1, 1 );
  //    glVertex2d  ( 1, 0 );
  //glEnd();
  //glDisable(GL_TEXTURE_1D);

  
  glColor3f(0,0,0);
  glLineWidth(2);
  draw_channel(channel);
}

void
ColorMapEditor::update_colors ()
{
    glBindTexture(GL_TEXTURE_1D,color_tex);
    glTexSubImage1D( 
        GL_TEXTURE_1D, 0, 0, resolution, 
        GL_RGBA, GL_UNSIGNED_BYTE, rgba
    );
}

void 
ColorMapEditor::draw_channel (Channel chan) 
{
    float *x = ((float*)hlsa) + chan;
    glBegin(GL_QUAD_STRIP);
    for (uint32_t i = 0; i < resolution; ++i, x+=4) {
        glVertex2f( i / (float)resolution, *x );
        glVertex2f( i / (float)resolution, 0 );
    }
    glEnd();
}



void
ColorMapEditor::stroke
(   Channel chan,
    float x0,
    float y0,
    float x1,
    float y1 )
{
    int ibegin, iend;
    
    ibegin = x0 * (resolution-1);
    iend  = x1   * (resolution-1);
    
    if (ibegin <  0 ) ibegin =  0 ;
    if (ibegin > (int(resolution)-1)) ibegin = (resolution-1);
    if (iend  <  0 ) iend  =  0 ;
    if (iend  > (int(resolution)-1)) iend  = (resolution-1);
    
    double vbegin = y0, vend = y1;
    if (vbegin < 0) vbegin = 0;
    if (vbegin > 1) vbegin = 1;
    if (vend   < 0) vend = 0;
    if (vend   > 1) vend = 1;

    if (ibegin == iend) {
        switch(chan) {
            case(Hue)   : hlsa[iend].h = vend; break;
            case(Lum)   : hlsa[iend].l = vend; break;
            case(Sat)   : hlsa[iend].s = vend; break;
            case(Alpha) : hlsa[iend].a = vend; break;
            default:;
        }
        write_hlsa_to_rgba(ibegin,iend);
        colors_need_update = true;
        return;
    }
    
    if (ibegin > iend) {
        //swap 
        int ti = ibegin;
        ibegin = iend;
        iend = ti;
        double tv = vbegin;
        vbegin = vend;
        vend = tv;
    }

    for (int i=ibegin; i<=iend; i++) {
        double dif = (i - ibegin) / (double)(iend - ibegin);
        double v = (1-dif)*vbegin + (dif)*vend;
        switch(chan) {
            case(Hue)   : hlsa[i].h = v; break;
            case(Lum)   : hlsa[i].l = v; break;
            case(Sat)   : hlsa[i].s = v; break;
            case(Alpha) : hlsa[i].a = v; break;
            default:;
        }
    }
    
    write_hlsa_to_rgba(ibegin,iend);
    colors_need_update = true;
}



void
ColorMapEditor::write_hlsa_to_rgba( int ibegin, int iend )
{
    HLSAf *x = hlsa+ibegin;
    RGBA8 *y = rgba+ibegin;
    for ( int i=ibegin; i<=iend; ++i,++x,++y ) {
        float h,l,s,r,g,b;
        h = x->h * 360.0;
        l = x->l ;
        s = x->s ;
        hls_to_rgb(h,l,s,r,g,b);
        y->r = (uint8_t)( r*255.0 );
        y->g = (uint8_t)( g*255.0 );
        y->b = (uint8_t)( b*255.0 );
        y->a = (uint8_t)( x->a*255.0 );
    }
}

void
ColorMapEditor::write_rgba_to_hlsa( int ibegin, int iend )
{
    RGBA8 *x = rgba+ibegin;
    HLSAf *y = hlsa+ibegin;
    float h=0,l,s,r,g,b;
    for ( int i=ibegin; i<=iend; ++i,++x,++y ) {
        r = x->r / 255.0;
        g = x->g / 255.0;
        b = x->b / 255.0;
        rgb_to_hls(r,g,b,h,l,s);
        y->h = h/360.0;
        y->l = l;
        y->s = s;
        y->a = x->a / 255.0; 
    }
}

