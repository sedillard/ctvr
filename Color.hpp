
#ifndef COLORS_H
#define COLORS_H

#include <stdint.h>

struct RGBA8{ uint8_t r,g,b,a; };
struct HLSA8{ uint8_t h,l,s,a; };


#define MAX3(a,b,c) (((a)>(c))?(((a)>(b))?(a):(b)):(((b)>(c))?(b):(c)))
#define MIN3(a,b,c) (((a)<(c))?(((a)<(b))?(a):(b)):(((b)<(c))?(b):(c)))

inline
float hls_value(float n1, float n2, float hue) 
{
  if (hue >= 360.0)
    hue -= 360.0;
  else if (hue < 0.0)
    hue += 360.0;

  if (hue < 60.0)
    return n1+(n2-n1)*hue/60.0;
  else if (hue < 180.0)
    return n2;
  else if (hue < 240.0)
    return n1+(n2-n1)*(240.0-hue)/60.0;
  else
    return n1;
} 

inline
void hls_to_rgb(float h, float l, float s, float &r, float &g, float &b) {
/* Given: h in [0,360] or UNDEFINED, l and s in [0,1]
 * Desired: r, g, b each in [0,1]
 */

  float m1,m2;

  m2 = (l<=0.5f)?(l+l*s):(l+s-l*s);
  m1 = (2.0f*l-m2);
  if (s==0.0f)
    r=g=b=l;
  else {
    r = hls_value(m1, m2, h+120.0f);
    g = hls_value(m1, m2, h);
    b = hls_value(m1, m2, h-120.0f);
  }
} 


inline 
void rgb_to_hls(float r, float g, float b, float &h, float &l, float &s) {
/* Given: r, b, g each in [0.1] 
 * Desired: h in [0,360), l and s in [0,1]
 */
 
  float max = MAX3(r,g,b);
  float min = MIN3(r,g,b);
  l = (max+min)/2.0f;		/* This is lightness */
  /* Next calculate saturation */
  if (max == min) {		/* Achromatic, because r = g = b */
    s = 0.0f;
    h = 0.0f; 			/* Actually: UNDEFINED */
  }
  else  {			/* Chromatic case */
    float delta = max-min;

    /* First calculate saturation */
    s = (l<=0.5f) ? (delta / (max+min)) : (delta / (2.0f - (max+min)));
    /* Next calculate the hue */
    if (r == max)
      h = (g-b)/delta;
    else if (g == max)
      h = 2.0f + (b-r)/delta;
    else if (b == max)
      h = 4.0f + (r-g)/delta;
    h *= 60.0f;
    if (h<0.0f)
      h += 360.0f;
  } /*Chromatic case */
} /* RGB_To_HLS */


#endif
