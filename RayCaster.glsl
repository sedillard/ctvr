
//The host must include a #define for MAX_ITRS and TF_RES in this source before
//compiling it. As such, line numbers on the GLSL compiler errors will be off
//by two

#pragma optionNV(fastmath on)

//Because we are reading from 8bit textures, the drivers might want to do math
//at low precision, but we don't want that.

//Convert small if blocks to cmovs 
#pragma optionNV(ifcvt all) 

#pragma optionNV(inline all)
#pragma optionNV(strict on)
//#pragma optionNV(unroll all)

uniform sampler3D scalar_tex;
uniform sampler3D aux_scalar_tex;
uniform sampler3D branch_map_tex;

uniform sampler2D parent_tex;
uniform sampler2D tf_index_tex;
uniform sampler2D depth_tex;
uniform sampler2D saddle_val_tex;

uniform vec4 branch_tex_size;

uniform sampler2D tf_tex;
uniform vec4 tf_tex_size; // { numCols, numRows, 1/numCols, 1/numRows }

uniform vec3 light_vec;
uniform vec3 view_vec;

//these are used for rounding. The x component is 1 fewer than the
//number of voxels in that dimension, and the y component is the inverse
//of that
uniform vec2 x_size; //  { size-1 , 1/(size-1) }
uniform vec2 y_size;
uniform vec2 z_size;

uniform vec3 x_inch; //small offsets for computing gradient
uniform vec3 y_inch;
uniform vec3 z_inch;

uniform sampler2D global_tf_tex;

uniform float isoval;

// x_down,y_up,etc..
// These round a texture coordinate up or down to the nearest 
// texel boundary in a given dimension. Logically:
//   x_down( vec3(3.76,2.34,5.55) ) = vec3(3.0,2.34,5.55)
// The problem is that OpenGL puts the texel sample in the center
// of a voxel, really we're rounding to k+0.5 for integers k.

// |---*---|---*---|---*---|---*---|---*---|---*---|  real tex space
//     |.......|.......|.......|.......|.......|      logical tex space
//      

vec3 x_down( in vec3 p )
{
    p.x *= x_size.r;
    p.x = floor(p.x-0.5)+0.5;
    p.x *= x_size.g;
    return p;
}


vec3 x_up( vec3 p )
{
    p.x *= x_size.r;
    p.x = ceil(p.x-0.5)+0.5;
    p.x *= x_size.g;
    return p;
}


vec3 y_down( vec3 p )
{
    p.y *= y_size.r;
    p.y = floor(p.y-0.5)+0.5;
    p.y *= y_size.g;
    return p;
}


vec3 y_up( vec3 p )
{
    p.y *= y_size.r;
    p.y = ceil(p.y-0.5)+0.5;
    p.y *= y_size.g;
    return p;
}


vec3 z_down( vec3 p )
{
    p.z *= z_size.r;
    p.z = floor(p.z-0.5)+0.5;
    p.z *= z_size.g;
    return p;
}



vec3 z_up( vec3 p )
{
    p.z *= z_size.r;
    p.z = ceil(p.z-0.5)+0.5;
    p.z *= z_size.g;
    return p;
}


void swap( inout vec3 a, inout vec3 b ) 
{
    vec3 t = a;
    a = b;
    b = t;
}

#if 0
vec3 normal(vec3 p)
{
    vec3 grad;
    grad.x = texture3D( scalar_tex, p + x_inch ).a - 
             texture3D( scalar_tex, p - x_inch ).a ;
    grad.y = texture3D( scalar_tex, p + y_inch ).a - 
             texture3D( scalar_tex, p - y_inch ).a ;
    grad.z = texture3D( scalar_tex, p + z_inch ).a - 
             texture3D( scalar_tex, p - z_inch ).a ;
    return normalize(grad);
}
#endif

vec3 normal(vec3 p, float f)
{
    vec3 grad;
    grad.x = texture3D( scalar_tex, p + x_inch ).a - f;
    grad.y = texture3D( scalar_tex, p + y_inch ).a - f;
    grad.z = texture3D( scalar_tex, p + z_inch ).a - f;
    return normalize(grad);
}

vec3 gradient(vec3 p, float f)
{
    vec3 grad;
    grad.x = texture3D( scalar_tex, p + x_inch ).a - f;
    grad.y = texture3D( scalar_tex, p + y_inch ).a - f;
    grad.z = texture3D( scalar_tex, p + z_inch ).a - f;
    return grad;
}


vec2 get_branch(vec3 p, float f)
{
    vec3 lo_pnt = x_down( p );
    vec3 hi_pnt = x_up  ( p );

    vec3 t_pnt;
    
    if ( texture3D(scalar_tex,lo_pnt).a > texture3D(scalar_tex,hi_pnt).a ) 
        swap(lo_pnt,hi_pnt);

    //look for a corner vertex with value no greater than f
    t_pnt  = y_up  (lo_pnt);
    lo_pnt = y_down(lo_pnt);
    if ( texture3D(scalar_tex,lo_pnt).a > texture3D(scalar_tex,t_pnt).a ) 
        swap(lo_pnt,t_pnt);

    t_pnt  = z_up  (lo_pnt);
    lo_pnt = z_down(lo_pnt);
    if ( texture3D(scalar_tex,lo_pnt).a > texture3D(scalar_tex,t_pnt).a ) 
        swap(lo_pnt,t_pnt);

    //look for a corner vertex with value no less than f
    t_pnt  = y_up  (hi_pnt);
    hi_pnt = y_down(hi_pnt);
    if ( texture3D(scalar_tex,hi_pnt).a < texture3D(scalar_tex,t_pnt).a ) 
        swap(hi_pnt,t_pnt);

    t_pnt  = z_up  (hi_pnt);
    hi_pnt = z_down(hi_pnt);
    if ( texture3D(scalar_tex,hi_pnt).a < texture3D(scalar_tex,t_pnt).a ) 
        swap(hi_pnt,t_pnt);

    vec2 hi_branch = texture3D(branch_map_tex, hi_pnt).ra;
    vec2 lo_branch = texture3D(branch_map_tex, lo_pnt).ra;
    float hi_depth = texture2D( depth_tex, hi_branch ).a;
    float lo_depth = texture2D( depth_tex, lo_branch ).a;

    vec2 the_branch;
   
    #if 1
    for(int i=0; i<MAX_ITRS; ++i) {
        float hi_saddle = texture2D( saddle_val_tex, hi_branch ).a;
        float lo_saddle = texture2D( saddle_val_tex, lo_branch ).a;
        vec2 hi_parent  = texture2D( parent_tex,     hi_branch ).ra;
        vec2 lo_parent  = texture2D( parent_tex,     lo_branch ).ra;
        if (hi_depth > lo_depth) {
            if (hi_saddle < f) {
                the_branch = hi_branch;
                break;
            } else {
                hi_branch  = hi_parent;
                hi_depth -= 1.0/256.0;
            }
        } else {
            if (lo_saddle > f) {
                the_branch = lo_branch;
                break;
            } else {
                lo_branch  = lo_parent;
                lo_depth -= 1.0/256.0;
            }
        }
        if (all(hi_branch == lo_branch)) {
            the_branch = hi_branch;
            break;
        }
    }
    #endif
    return the_branch;
}


vec4 sample_branch_tf( vec2 branch, float f )
{
    vec2 tf_index = texture2D( tf_index_tex, branch ).ra;
    
    tf_index.x += f * float(TF_RES) * tf_tex_size.z ;
    if (tf_index.x > 1.0) {
        tf_index.x -= 1.0;
        tf_index.y += tf_tex_size.w;
    }

    return texture2D( tf_tex, tf_index );
}


void main ()
{
  vec3 p = gl_TexCoord[0].xyz;
  vec4 color = vec4(0,0,0,0);

  vec3 box_bounds = vec3 (x_size.x,y_size.x,z_size.x);

  float f0 = texture3D(scalar_tex,p).a;
  f0 = clamp(f0,0.0,255.0/256.0);
  vec3 g0 = gradient(p,f0);

  for ( int i=0; i<256; ++i ) {
    p += view_vec;
    float f1 = texture3D(scalar_tex,p).a;
    f1 = clamp(f1,0.0,255.0/256.0);
    vec3 g1 = gradient(p,f1);
    
    vec4 glob = texture2D(global_tf_tex,vec2(f1,0.0));

    vec4 g = (f0 - isoval) * (f1 - isoval) > 0 ? vec4(0,0,0,0) : vec4(1,1,1,0.2);
    float t = abs(f0-isoval) / (abs(f1-isoval)+abs(f0-isoval));
    
    vec3 nrml = normalize( (1.0-t)*g0 + t*g1 ) ;
    
    g.rgb *= abs(dot(light_vec,nrml));

    vec2 br = get_branch( p, f1 );
    vec4 c = sample_branch_tf(br,f1);
    c.a *= glob.a; 
    c = c+g;
    c = clamp(c,0.0,1.0);


    
    float a = (1.0-color.a)*c.a;
    color.rgb += a*c.rgb;
    color.a += a;
    f0 = f1;
    g0 = g1;

    if (any (lessThan( p, vec3(0.0,0.0,0.0))) ||
        any (greaterThan(p, vec3(1.0,1.0,1.0)))) break;
  }

  gl_FragColor = color;
}



