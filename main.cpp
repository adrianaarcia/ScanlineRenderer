#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <iostream>
#include <cfloat>
#include <math.h>

#include "SETTINGS.h"

using namespace std;

//////////////////////////////////////////////////////////////////////////////////
// Function declarations
//////////////////////////////////////////////////////////////////////////////////
void draw(const string& filename, vector<VEC4>& vertices, vector<VEC3I>& indices, vector<VEC3>& colors, int xres, int yres, float* values, int mode);
void color_vertices(const vector<VEC4>vertices, vector<VEC3I>& indices, vector<VEC3>& colors, int xres, int yres, float* values);

// Geometry Functions
void construct_Mvp(float nx, float ny, MATRIX4& Mvp);
void construct_Mortho(float l, float r, float t, float b, float n, float f, MATRIX4& Mortho);
void construct_Mcam(VEC3 eye, VEC3 lookat, VEC3 up, MATRIX4& Mcam);
void construct_Mper(float fovy, float aspect, float n, float f, MATRIX4& Mper);
float imp_lin_eq(float x, float y, int i, int j, const float* xcoords, const float* ycoords);
void build_big_square(vector<VEC3>& vertices, vector<VEC3I>& indices, vector<VEC3>& colors);
void build_square(vector<VEC3>& vertices, vector<VEC3I>& indices, vector<VEC3>& colors);
void build_cube(vector<VEC3>& vertices, vector<VEC3I>& indices, vector<VEC3>& colors);
void build_cube_vertex_colors(vector<VEC3>& vertices, vector<VEC3I>& indices, vector<VEC3>& colors);

// Vector Functions
VEC3 truncate(const VEC4& v);
VEC4 extend(const VEC3& v);
void transform(vector<VEC4>& vertices, const MATRIX4& M);
void scale_vecs(vector<VEC4>& vertices, float scalar);
void extend_vecs(const vector<VEC3>& vertices, vector<VEC4>& extended);
void truncate_vecs(const vector<VEC4>& vertices, vector<VEC3>& truncated);

// Utility Functions
void clear_floats(int size, float* values);
void read_ppm(const string& filename, int& xres, int& yres, float*& values);
void write_ppm(const string& filename, int& xres, int& yres, const float* values);

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////

/* Viewport Matrix */
void viewport_matrix(vector<VEC4>& extended, vector<VEC3>& vertices, vector<VEC3I>& indices, vector<VEC3>& colors, MATRIX4& Mvp, float nx, float ny, float* values) 
{
  int xres = (int) nx; int yres = (int) ny; int size = xres*yres*3;
  clear_floats(size, values);

  construct_Mvp(nx, ny, Mvp); //construct viewport matrix

  build_square(vertices, indices, colors);//get square
  extend_vecs(vertices, extended);
  
  transform(extended, Mvp); //apply viewport matrix
  color_vertices(extended, indices, colors, xres, yres, values);

  write_ppm("1_viewport_matrix.ppm", xres, yres, values);
}

/* Triangle Rasterization */
void triangle_rasterization(vector<VEC4>& extended, vector<VEC3I>& indices, vector<VEC3>& colors, int xres, int yres, int size, float* values)
{
  clear_floats(size, values); //clear values

  // rasterization occurs in draw()
  draw("2_triangle_rasterization.ppm", extended, indices, colors, xres, yres, values, 0); //draw and output, mode = 0 (simple rasterization)
}

/* Orthographic Projection Matrix */
void ortho_proj_matrix(vector<VEC4>& extended, vector<VEC3>& vertices, vector<VEC3I>& indices, vector<VEC3>& colors,
const MATRIX4 Mvp, MATRIX4& Mortho, int xres, int yres, int size, float* values) 
{
  //reinitialize vars
  clear_floats(size, values);
  extended.clear(); vertices.clear(); indices.clear(); colors.clear(); 
  
  build_big_square(vertices, indices, colors); //get big square
  extend_vecs(vertices, extended);

  //construct Mortho
  float l, r, t, b, n, f; l = b = f = 0.0; r = t = n = 12.0; //planes
  construct_Mortho(l, r, t, b, n, f, Mortho);

  //apply transforms
  transform(extended, Mortho); //apply Mortho
  transform(extended, Mvp); //apply Mvp

  //draw and output
  draw("3_ortho_proj_matrix.ppm", extended, indices, colors, xres, yres, values, 0);
}

/* Camera Matrix */
void camera_matrix(vector<VEC4>& extended, vector<VEC3I>& indices, vector<VEC3>& colors, MATRIX4& Mcam, int xres, int yres, int size, float* values) 
{
  clear_floats(size, values); //reinitialize some vars

  VEC3 eye(0.2, 0.2, 1.0); VEC3 lookat(0,0,0); VEC3 up(0,1,0); //given vectors
  construct_Mcam(eye, lookat, up, Mcam); //construct Mcam
  
  transform(extended, Mcam); //apply Mcam
  draw("4_camera_matrix.ppm", extended, indices, colors, xres, yres, values, 0); //draw and output
}

/* Perspective Projection Matrix */
void persp_proj_matrix(vector<VEC4>& extended, vector<VEC3>& vertices, vector<VEC3I>& indices, vector<VEC3>& colors, MATRIX4& Mvp, MATRIX4& Mcam, MATRIX4& Mper, 
int xres, int yres, int size, float* values)
{
  /* Reinitialize values */
  clear_floats(size, values);
  extended.clear(); vertices.clear(); indices.clear(); colors.clear();
  
  build_cube(vertices, indices, colors); //get cube
  extend_vecs(vertices, extended);
  scale_vecs(extended, 0.5); //scale by factor of 0.5
  
  /* Apply transforms */

  //Mcam
  VEC3 eye(1.0, 1.0, 1.0); VEC3 lookat(0.0,0.0,0.0); VEC3 up(0.0,1.0,0.0); //given vectors
  construct_Mcam(eye, lookat, up, Mcam);
  transform(extended, Mcam);

  //Mper
  float fovy, aspect, n, f; n = 1.0; f = 100.0; fovy = 65.0; aspect = 4.0/3.0; //given values
  construct_Mper(fovy, aspect, n, f, Mper);
  transform(extended, Mper);

  //Ignore w for now -- set to 1.0
  vector<VEC3> truncated; vector<VEC4> verts;
  truncate_vecs(extended, truncated); 
  extend_vecs(truncated, verts);

  //Mvp
  transform(verts, Mvp);

  //draw and output
  draw("5_perspective_proj_matrix.ppm", verts, indices, colors, xres, yres, values, 0);
}
  
/* Perspective Divide */
void p_divide(vector<VEC4>& extended)
{
  int count = extended.size(); float w = 0.0;
  for(int i = 0; i < count; i++) 
  {
    w = extended[i][3];
    extended[i] /= w; //scale each vector by 1/w
  }
}

void persp_divide(vector<VEC4>& extended, vector<VEC3I>& indices, vector<VEC3>& colors, MATRIX4& Mvp, 
int xres, int yres, int size, float* values)
{
  clear_floats(size, values);//reinitialize values array
  
  p_divide(extended);//perspective divide
  transform(extended, Mvp); //apply Mvp

  draw("6_perspective_divide.ppm", extended, indices, colors, xres, yres, values, 0); //draw and output
}

/* Z-buffering */
void z_buff(vector<VEC4>& extended, vector<VEC3I>& indices, vector<VEC3>& colors, int xres, int yres, int size, float* values)
{
  clear_floats(size, values); //clear values array

  // z-buffering occurs in draw()
  draw("7_z_buffering.ppm", extended, indices, colors, xres, yres, values, 1); //draw and output, mode = 1 (z buffering)
}

/* Color Interpolation */
void color_interpolation(vector<VEC4>& extended, vector<VEC3>& vertices, vector<VEC3I>& indices, vector<VEC3>& colors,
int xres, int yres, int size, float* values)
{
  //reinitialize values
  clear_floats(size, values);
  vertices.clear(); indices.clear(); colors.clear(); //extended.clear()

  build_cube_vertex_colors(vertices, indices, colors); //get new colors

  //color interpolation occurs in draw()
  draw("8_color_interpolation.ppm", extended, indices, colors, xres, yres, values, 2); //draw and output, mode = 2 (z-buffering + interpolation)
}

//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
int main(int argc, char* argv[])
{ 
  /* Declare & initialize vars */
  
  //dimensions of output image
  float nx = 800.0;
  float ny = 600.0;
  int xres = (int) nx;
  int yres = (int) ny;
  int size = xres*yres*3;

  //initalize array of values to zero
  float values[size]; clear_floats(size, values);

  //our matrices
  MATRIX4 Mvp; MATRIX4 Mortho; MATRIX4 Mcam; MATRIX4 Mper;
  Mvp.setIdentity(); Mortho.setIdentity(); Mcam.setIdentity(); Mper.setIdentity();

  //vectors
  vector<VEC3> vertices; 
  vector<VEC3I> indices; 
  vector<VEC3> colors; 
  vector<VEC4> extended;
  
  int parameter_count = 9; //expected number of input params
  float cam[parameter_count]; 
  clear_floats(sizeof(float)*parameter_count, cam); //by default camera parameters are set to zero

  if(argc > 1) //check for user input
  { 
    int inputs = argc - 1; //exclude program name
    int lim = min(inputs, parameter_count);
    
    for(int i = 0; i < lim; i++) //convert input to floats -- store in cam array
    {
      cam[i] = atof(argv[i+1]);
    }

    //store user input
    int vec_num = 3;
    vector<VEC3> parameters;

    for(int i = 0; i < vec_num; i++) //for each vector;
    {
      parameters.push_back(VEC3(0,0,0));
      for(int j = 0; j < 3; j++) //for each value
      {
        parameters[i][j] = cam[j+3*i];
      }
    }

    build_cube(vertices, indices, colors); //get cube
    extend_vecs(vertices, extended);
    scale_vecs(extended, 0.5); //scale the vertices by 0.5

    //Mcam
    construct_Mcam(parameters[0], parameters[1], parameters[2], Mcam);
    transform(extended, Mcam);

    //Mper
    float fovy, aspect, n, f; n = 1.0; f = 100.0; fovy = 65.0; aspect = 4.0/3.0; //given values
    construct_Mper(fovy, aspect, n, f, Mper);
    transform(extended, Mper);  

    p_divide(extended); //perspective divide
    
    //Mvp
    construct_Mvp(nx,ny,Mvp);
    transform(extended, Mvp); 

    //draw and output
    draw("custom.ppm", extended, indices, colors, xres, yres, values, 2);
  }
  else
  {
    //Default -- demonstrate process visually, step by step
    viewport_matrix(extended, vertices, indices, colors, Mvp, nx, ny, values); //create viewport matrix
    triangle_rasterization(extended, indices, colors, xres, yres, size, values); //triangle rasterization
    ortho_proj_matrix(extended, vertices, indices, colors, Mvp, Mortho, xres, yres, size, values); // create and apply the orthographic projection matrix
    camera_matrix(extended, indices, colors, Mcam, xres, yres, size, values); // create and apply the camera matrix 
    persp_proj_matrix(extended, vertices, indices, colors, Mvp, Mcam, Mper, xres, yres, size, values); // create and apply the perspective projection matrix
    persp_divide(extended, indices, colors, Mvp, xres, yres, size, values); // perspective divide
    z_buff(extended, indices, colors, xres, yres, size, values); // z buffering
    color_interpolation(extended, vertices, indices, colors, xres, yres, size, values); // color interpolation
  }
  
  return 0;
}

/* draw()
 * Colors and outputs the image. There are three possible modes.
 * 
 * mode = 1 ==> simple triangle rasterization
 * mode = 2 ==> rasterization + z-buffering
 * mode = 3 ==> rasterization + z-buffering + color interpolation
 * 
 */
void draw(const string& filename, vector<VEC4>& vertices, vector<VEC3I>& indices, vector<VEC3>& colors, int xres, int yres, float* values, int mode)
{
  // initialize z-buffer
  int size = xres*yres;
  float depths[size];
  for(int i = 0; i < size; i++)
  {
    depths[i] = FLT_MAX;
  }

  //initialize vars
  VEC3 vert_color; vert_color.setZero();
  float fa, fb, fg, alpha, beta, gamma, curr_coord, xmin, xmax, ymin, ymax, xi, yi, curr_depth;
  fa = fb = fg = alpha = beta = gamma = curr_coord = xmin = xmax = ymin = ymax = xi = yi = curr_depth = 0;
  float barycentric[3];
  int index= 0;
  int triangle_count = indices.size();
  float xcoords[3]; float ycoords[3]; float zcoords[3];
  
  for(int i = 0; i < triangle_count; i++) //for each triangle
  {
    xmin = ymin = FLT_MAX; xmax = ymax = FLT_MIN;
    
    for(int j = 0; j < 3; j++) //get coords of each vertex, store by x,y,z
    {
      index = indices[i][j];
  
      xcoords[j] = xi = vertices[index].x(); 
      ycoords[j] = yi = vertices[index].y();
      zcoords[j] = vertices[index].z();

      //update bounds
      if(xi < xmin) xmin = xi;
      if(xi > xmax) xmax = xi;
      if(yi < ymin) ymin = yi;
      if(yi > ymax) ymax = yi;
    }

    xmin = floor(xmin); ymin = floor(ymin); xmax = ceil(xmax); ymax = ceil(ymax); //update bounds

    //calculate subtriangle areas
    fa = imp_lin_eq(xcoords[0], ycoords[0], 1, 2, xcoords, ycoords);
    fb = imp_lin_eq(xcoords[1], ycoords[1], 2, 0, xcoords, ycoords);
    fg = imp_lin_eq(xcoords[2], ycoords[2], 0, 1, xcoords, ycoords);

    for(int y = ymin; y <= ymax; y++)
    {
      for(int x = xmin; x <= xmax; x++)
      {
        //calculate barycentric coordinates
        barycentric[0] = alpha = imp_lin_eq((float)x, (float)y, 1, 2, xcoords, ycoords)/fa;
        barycentric[1] = beta = imp_lin_eq((float)x, (float)y, 2, 0, xcoords, ycoords)/fb;
        barycentric[2] = gamma = imp_lin_eq((float)x, (float)y, 0, 1, xcoords, ycoords)/fg;

        if(alpha >= 0 && beta >= 0 && gamma >= 0)
        {
          if((alpha > 0 || (fa*imp_lin_eq(-1, -1, 1, 2, xcoords, ycoords)) >0) && 
             (beta > 0 || (fb*imp_lin_eq(-1, -1, 2, 0, xcoords, ycoords))>0) && 
             (gamma > 0 || (fg*imp_lin_eq(-1, -1, 0, 1, xcoords, ycoords))>0))
          {
            index = (yres-y)*xres + x;
            curr_depth = alpha*zcoords[0] + beta*zcoords[1] + gamma*zcoords[2];
            if(depths[index] > curr_depth || mode <= 0) //check depth & mode -- color appropriately
            {
              depths[index] = curr_depth;
              index *= 3;
          
              if(mode < 2)
              {
                for(int j = 0; j < 3; j++) //set each color
                {
                  values[index+j] = 255*colors[i][j];
                }
              }
              else //color interpolation
              {
                for(int j = 0; j < 3; j++) //set each color
                {
                  values[index+j] = 0;

                  for(int k = 0; k < 3; k++) //with info from each vertex & barycentric coords
                  {
                    values[index+j] += barycentric[k]*colors[indices[i][k]][j];
                  }

                  values[index+j] *= 255;
                }
              }
            }
          }
        }
      }
    }
  }

  write_ppm(filename, xres, yres, values);
}

/* color_vertices()
 * A simplified version of draw() that only colors the vertices of the triangles.
 */
void color_vertices(const vector<VEC4> vertices, vector<VEC3I>& indices, vector<VEC3>& colors, int xres, int yres, float* values)
{
  int triangle_count = indices.size();
  VEC4 interm; interm.setZero();
  float x, y, z; int index; x = y = z = 0.0; index = 0;

  for (int i = 0; i < triangle_count; i++) //for each triangle
  {
    for(int j = 0; j < 3; j++) //for each vertex in the triangle
    {
      index = indices[i][j];
      interm = vertices[index];
      
      //get index of pixel
      x = interm[0];
      y = interm[1];
      index = (yres-y)*xres*3 + 3*x;
      
      //color vertex
      values[index] = 255*colors[i][0]; 
      values[index+1] =255*colors[i][1];
      values[index+2] = 255*colors[i][2];
    }
  }
}

////////////////////////////////////////////////////////////////////////////////////////////////
// Geometry Functions
////////////////////////////////////////////////////////////////////////////////////////////////
/* construct_Mvp()
 * Creates viewport matrix.
 */
void construct_Mvp(float nx, float ny, MATRIX4& Mvp)
{
  Mvp.setIdentity();

  Mvp(0,0) = nx/2.0; Mvp(0,3) = nx/2.0;
  Mvp(1,1) = ny/2.0; Mvp(1, 3) = ny/2.0;
}

/* construct_Mortho()
 * Creates orthographic projection matrix.
 */
void construct_Mortho(float l, float r, float t, float b, float n, float f, MATRIX4& Mortho)
{
  Mortho.setIdentity();

  Mortho(0,0) = 2.0/(r-l); Mortho(1,1) = 2.0/(t-b); Mortho(2,2) = -2.0/(n-f);
  Mortho(0,3) = -(r+l)/(r-l); Mortho(1,3) = -(t+b)/(t-b); Mortho(2,3) = -(n+f)/(n-f);
}

/* construct_Mcam()
 * Creates camera matrix.
 */
void construct_Mcam(VEC3 eye, VEC3 lookat, VEC3 up, MATRIX4& Mcam)
{
  //initialize vars
  Mcam.setIdentity(); 
  VEC3 w; w.setZero(); VEC3 u; u.setZero(); VEC3 v; v.setZero();
  
  //calculate basis vectors
  VEC3 g; g = eye - lookat;
  w = g.normalized();
  u = up.cross(w); u.normalize();
  v = w.cross(u);
  
  MATRIX4 Mbasis; Mbasis.setIdentity();
  Mbasis(0,0) = u.x(); Mbasis(0,1) = u.y(); Mbasis(0,2) = u.z();
  Mbasis(1,0) = v.x(); Mbasis(1,1) = v.y(); Mbasis(1,2) = v.z();
  Mbasis(2,0) = w.x(); Mbasis(2,1) = w.y(); Mbasis(2,2) = w.z();
  
  MATRIX4 Meye; Meye.setIdentity();
  Meye(0,3) = -eye.x(); Meye(1,3) = -eye.y(); Meye(2,3) = -eye.z();
  
  Mcam = Mbasis*Meye; //Mcam is the product of the above matrices
}

/* construct_Mper()
 * Creates perspective projection matrix.
 */
void construct_Mper(float fovy, float aspect, float n, float f, MATRIX4& Mper)
{
  Mper.setZero();

  float theta = fovy* M_PI/180; //fovy in radians
  float l,r,t,b; //planes;
  t = n*tanf(theta/2); b = -t;
  r = t*aspect; l = -r;

  float temp1, temp2, temp3, temp4; //store repeated calculations
  temp1 = 2*n; temp2 = r-l; temp3 = n-f; temp4 = t-b;

  Mper(0,0) = temp1/(r-l); Mper(0,2) = (r+l)/(r-l);
  Mper(1,1) = temp1/temp4; Mper(1,2) = (t+b)/temp4; 
  Mper(2, 2) = (n+f)/temp3; Mper(2,3) = (temp1*f)/temp3;
  Mper(3,2) = -1;
}

/* imp_lin_eq()
 * Calculates the result of the implicit line equation given an array of coordinates and the
 * indices of the relevant coordinates.
 */
float imp_lin_eq(float x, float y, int i, int j, const float* xcoords, const float* ycoords)
{
  int temp = 0;
  
  if(xcoords[i] > xcoords[j]) //swap points
  {
    temp = i;
    i = j;
    j = temp;
  }

  return (ycoords[i]-ycoords[j])*x + (xcoords[j]-xcoords[i])*y + xcoords[i]*ycoords[j] - xcoords[j]*ycoords[i];
}

/* build_big_square()
 * build out a single square
 * 
 * provided by Theodore Kim
 */
void build_big_square(vector<VEC3>& vertices, vector<VEC3I>& indices, vector<VEC3>& colors)
{
  vertices.push_back(VEC3(1.0, 1.0, 1.5));
  vertices.push_back(VEC3( 11.0, 1.0, 1.5));
  vertices.push_back(VEC3(1.0,  11.0, 1.5));
  vertices.push_back(VEC3( 11.0,  11.0, 1.5));

  // front face
  indices.push_back(VEC3I(0, 1, 2));
  indices.push_back(VEC3I(2, 1, 3));
  colors.push_back(VEC3(1.0, 0.0, 0.0));
  colors.push_back(VEC3(0.0, 0.0, 1.0));
}

/* build_square()
 * build out a single square
 * 
 * provided by Theodore Kim
 */
void build_square(vector<VEC3>& vertices, vector<VEC3I>& indices, vector<VEC3>& colors)
{
  vertices.push_back(VEC3(-0.5, -0.5,  0.5)); //0 -- (200, 150)
  vertices.push_back(VEC3( 0.5, -0.5,  0.5)); //1 -- (600, 150)
  vertices.push_back(VEC3(-0.5,  0.5,  0.5)); //2 -- (200, 450)
  vertices.push_back(VEC3( 0.5,  0.5,  0.5)); //3 -- (600, 450)

  // front face
  indices.push_back(VEC3I(0, 1, 2));
  indices.push_back(VEC3I(2, 1, 3));
  colors.push_back(VEC3(1.0, 0.0, 0.0));
  colors.push_back(VEC3(0.0, 1.0, 0.0));
}

/* build_cube()
 * build out a cube
 *
 * provided by Theodore Kim
 */
void build_cube(vector<VEC3>& vertices, vector<VEC3I>& indices, vector<VEC3>& colors)
{
  vertices.push_back(VEC3(-0.5, -0.5,  0.5));
  vertices.push_back(VEC3( 0.5, -0.5,  0.5));
  vertices.push_back(VEC3(-0.5,  0.5,  0.5));
  vertices.push_back(VEC3( 0.5,  0.5,  0.5));
  vertices.push_back(VEC3(-0.5, -0.5, -0.5));
  vertices.push_back(VEC3( 0.5, -0.5, -0.5));
  vertices.push_back(VEC3(-0.5,  0.5, -0.5));
  vertices.push_back(VEC3( 0.5,  0.5, -0.5));

  // front face
  indices.push_back(VEC3I(0, 1, 2));
  indices.push_back(VEC3I(2, 1, 3));
  colors.push_back(VEC3(1.0, 0.0, 0.0));
  colors.push_back(VEC3(1.0, 0.0, 0.0));

  // back face
  indices.push_back(VEC3I(5, 4, 7));
  indices.push_back(VEC3I(7, 4, 6));
  colors.push_back(VEC3(0.0, 1.0, 0.0));
  colors.push_back(VEC3(0.0, 1.0, 0.0));

  // left face
  indices.push_back(VEC3I(4, 0, 6));
  indices.push_back(VEC3I(6, 0, 2));
  colors.push_back(VEC3(0.0, 0.0, 1.0));
  colors.push_back(VEC3(0.0, 0.0, 1.0));

  // right face
  indices.push_back(VEC3I(1, 5, 3));
  indices.push_back(VEC3I(3, 5, 7));
  colors.push_back(VEC3(0.0, 1.0, 1.0));
  colors.push_back(VEC3(0.0, 1.0, 1.0));

  // top face
  indices.push_back(VEC3I(2, 3, 6));
  indices.push_back(VEC3I(6, 3, 7));
  colors.push_back(VEC3(1.0, 1.0, 0.0));
  colors.push_back(VEC3(1.0, 1.0, 0.0));

  // bottom face
  indices.push_back(VEC3I(4, 5, 0));
  indices.push_back(VEC3I(0, 5, 1));
  colors.push_back(VEC3(1.0, 0.0, 1.0));
  colors.push_back(VEC3(1.0, 0.0, 1.0));
}

/* build_cube_vertex_colors()
 * build out a cube
 * 
 * provided by Theodore Kim
 */
void build_cube_vertex_colors(vector<VEC3>& vertices, vector<VEC3I>& indices, vector<VEC3>& colors)
{
  vertices.push_back(VEC3(-0.5, -0.5,  0.5));
  vertices.push_back(VEC3( 0.5, -0.5,  0.5));
  vertices.push_back(VEC3(-0.5,  0.5,  0.5));
  vertices.push_back(VEC3( 0.5,  0.5,  0.5));
  vertices.push_back(VEC3(-0.5, -0.5, -0.5));
  vertices.push_back(VEC3( 0.5, -0.5, -0.5));
  vertices.push_back(VEC3(-0.5,  0.5, -0.5));
  vertices.push_back(VEC3( 0.5,  0.5, -0.5));
  colors.push_back(VEC3(1.0, 0.0, 0.0));
  colors.push_back(VEC3(1.0, 0.0, 0.0));
  colors.push_back(VEC3(0.0, 1.0, 0.0));
  colors.push_back(VEC3(0.0, 1.0, 0.0));
  colors.push_back(VEC3(0.0, 0.0, 1.0));
  colors.push_back(VEC3(0.0, 0.0, 1.0));
  colors.push_back(VEC3(1.0, 1.0, 0.0));
  colors.push_back(VEC3(1.0, 1.0, 0.0));

  // front face
  indices.push_back(VEC3I(0, 1, 2));
  indices.push_back(VEC3I(2, 1, 3));

  // back face
  indices.push_back(VEC3I(5, 4, 7));
  indices.push_back(VEC3I(7, 4, 6));

  // left face
  indices.push_back(VEC3I(4, 0, 6));
  indices.push_back(VEC3I(6, 0, 2));

  // right face
  indices.push_back(VEC3I(1, 5, 3));
  indices.push_back(VEC3I(3, 5, 7));

  // top face
  indices.push_back(VEC3I(2, 3, 6));
  indices.push_back(VEC3I(6, 3, 7));

  // bottom face
  indices.push_back(VEC3I(4, 5, 0));
  indices.push_back(VEC3I(0, 5, 1));
}

////////////////////////////////////////////////////////////////////////////////////////////////
// Vector Functions
////////////////////////////////////////////////////////////////////////////////////////////////
/* truncate()
 * Drops the final value in a VEC4 to return a VEC3.
 */
VEC3 truncate(const VEC4& v)
{
  return VEC3(v[0], v[1], v[2]);
}

/* extend()
 * Adds a value of '1' to a VEC3 to return a VEC4.
 */
VEC4 extend(const VEC3& v)
{
  return VEC4(v[0], v[1], v[2], 1.0);
}

/* transform()
 * Applies matrix to a vector of VEC4s.
 */
void transform(vector<VEC4>& vertices, const MATRIX4& M)
{
  int verts_count = vertices.size();

  for(int i = 0; i < verts_count; i++)
  {
    vertices[i] = M*vertices[i];
  }
}

/* scale_vecs()
 * Scales all VEC4s in the vector by the scalar
 */
void scale_vecs(vector<VEC4>& vertices, float scalar)
{
  int verts_count = vertices.size();
  
  for(int i = 0; i < verts_count; i++) 
  {
    vertices[i] *= scalar;
  }
}

/* extend_vecs()
 * Extends all VEC3s in vertices vector and adds the new VEC4 to the
 * extended vector.
 */
void extend_vecs(const vector<VEC3>& vertices, vector<VEC4>& extended)
{
  VEC4 temp; temp.setZero();
  int count = vertices.size();
  
  for(int i = 0; i < count; i++)
  {
    temp = extend(vertices[i]);
    extended.push_back(temp);
  }
}

/* truncate_vecs()
 * Truncates all VEC4s in vertices vector and adds the new VEC3 to the
 * truncated vector.
 */
void truncate_vecs(const vector<VEC4>& vertices, vector<VEC3>& truncated)
{
  VEC3 temp; temp.setZero();
  int count = vertices.size();
  
  for(int i = 0; i < count; i++)
  {
    temp = truncate(vertices[i]);
    truncated.push_back(temp);
  }
}

/////////////////////////////////////////////////////////////////////
// Utility Functions
/////////////////////////////////////////////////////////////////////
/* clear_floats()
 * sets all values in a float array to 0
 */
void clear_floats(int size, float* values)
{
  memset(values,0,sizeof(float)*size);
}

//provided by Theodore Kim
void read_ppm(const string& filename, int& xres, int& yres, float*& values)
{
  // try to open the file
  FILE *fp;
  fp = fopen(filename.c_str(), "rb");
  if (fp == NULL)
  {
    cout << " Could not open file \"" << filename.c_str() << "\" for reading." << endl;
    cout << " Make sure you're not trying to read from a weird location or with a " << endl;
    cout << " strange filename. Bailing ... " << endl;
    exit(0);
  }

  // get the dimensions
  fscanf(fp, "P6\n%d %d\n255\n", &xres, &yres);
  int total_cells = xres * yres;

  // grab the pixel values
  unsigned char* pixels = new unsigned char[3 * total_cells];
  fread(pixels, 1, total_cells * 3, fp);

  // copy to a nicer data type
  values = new float[3 * total_cells];
  for (int i = 0; i < 3 * total_cells; i++)
    values[i] = pixels[i];

  // clean up
  delete[] pixels;
  fclose(fp);
  cout << " Read in file " << filename.c_str() << endl;
}

//provided by Theodore Kim
void write_ppm(const string& filename, int& xres, int& yres, const float* values)
{
  int total_cells = xres * yres;
  unsigned char* pixels = new unsigned char[3 * total_cells];
  for (int i = 0; i < 3 * total_cells; i++)
    pixels[i] = values[i];

  FILE *fp;
  fp = fopen(filename.c_str(), "wb");
  if (fp == NULL)
  {
    cout << " Could not open file \"" << filename.c_str() << "\" for writing." << endl;
    cout << " Make sure you're not trying to write from a weird location or with a " << endl;
    cout << " strange filename. Bailing ... " << endl;
    exit(0);
  }

  fprintf(fp, "P6\n%d %d\n255\n", xres, yres);
  fwrite(pixels, 1, total_cells * 3, fp);
  fclose(fp);
  delete[] pixels;
}
