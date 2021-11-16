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
void draw(const string& filename, vector<VEC4>& vertices, vector<VEC3I>& indices, vector<VEC3>& colors, int xRes, int yRes, float* values, int mode);
void transform(vector<VEC4>& vertices, const MATRIX4& M);
void constructMvp(float nx, float ny, MATRIX4& Mvp);
void constructMortho(float l, float r, float t, float b, float n, float f, MATRIX4& Mortho);
void constructMcam(VEC3 eye, VEC3 lookAt, VEC3 up, MATRIX4& Mcam);
void constructMper(float fovy, float aspect, float n, float f, MATRIX4& Mper);
void scale_vecs(vector<VEC4>& vertices, float scalar);
void extend_vecs(const vector<VEC3>& vertices, vector<VEC4>& extended);
void truncate_vecs(const vector<VEC4>& vertices, vector<VEC3>& truncated);
void color_vertices(const vector<VEC4>vertices, vector<VEC3I>& indices, vector<VEC3>& colors, int xRes, int yRes, float* values);
float imp_lin_eq(float x, float y, int i, int j, const float* xcoords, const float* ycoords);
//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
VEC3 truncate(const VEC4& v)
{
  return VEC3(v[0], v[1], v[2]);
}
VEC4 extend(const VEC3& v)
{
  return VEC4(v[0], v[1], v[2], 1.0);
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void readPPM(const string& filename, int& xRes, int& yRes, float*& values)
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
  fscanf(fp, "P6\n%d %d\n255\n", &xRes, &yRes);
  int totalCells = xRes * yRes;

  // grab the pixel values
  unsigned char* pixels = new unsigned char[3 * totalCells];
  fread(pixels, 1, totalCells * 3, fp);

  // copy to a nicer data type
  values = new float[3 * totalCells];
  for (int i = 0; i < 3 * totalCells; i++)
    values[i] = pixels[i];

  // clean up
  delete[] pixels;
  fclose(fp);
  cout << " Read in file " << filename.c_str() << endl;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void writePPM(const string& filename, int& xRes, int& yRes, const float* values)
{
  int totalCells = xRes * yRes;
  unsigned char* pixels = new unsigned char[3 * totalCells];
  for (int i = 0; i < 3 * totalCells; i++)
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

  fprintf(fp, "P6\n%d %d\n255\n", xRes, yRes);
  fwrite(pixels, 1, totalCells * 3, fp);
  fclose(fp);
  delete[] pixels;
}

//////////////////////////////////////////////////////////////////////////////////
// build out a single square
//////////////////////////////////////////////////////////////////////////////////
void buildBigSquare(vector<VEC3>& vertices, vector<VEC3I>& indices, vector<VEC3>& colors)
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

//////////////////////////////////////////////////////////////////////////////////
// build out a single square
//////////////////////////////////////////////////////////////////////////////////
void buildSquare(vector<VEC3>& vertices, vector<VEC3I>& indices, vector<VEC3>& colors)
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

//////////////////////////////////////////////////////////////////////////////////
// build out a cube
//////////////////////////////////////////////////////////////////////////////////
void buildCube(vector<VEC3>& vertices, vector<VEC3I>& indices, vector<VEC3>& colors)
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

//////////////////////////////////////////////////////////////////////////////////
// build out a cube
//////////////////////////////////////////////////////////////////////////////////
void buildCubePerVertexColors(vector<VEC3>& vertices, vector<VEC3I>& indices, vector<VEC3>& colors)
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

/* Problem 1: Viewport Matrix */
void p1(vector<VEC4>& extended, vector<VEC3>& vertices, vector<VEC3I>& indices, vector<VEC3>& colors, MATRIX4& Mvp, float nx, float ny, float* values) 
{
  int xRes = (int) nx; int yRes = (int) ny; int size = xRes*yRes*3;
  memset(values,0,sizeof(float)*size);

  constructMvp(nx, ny, Mvp); //construct viewport matrix

  buildSquare(vertices, indices, colors);//get square
  extend_vecs(vertices, extended);
  
  transform(extended, Mvp); //apply viewport matrix
  color_vertices(extended, indices, colors, xRes, yRes, values);

  writePPM("1.ppm", xRes, yRes, values);
}

/* Problem 2: Triangle Rasterization */
void p2(vector<VEC4>& extended, vector<VEC3I>& indices, vector<VEC3>& colors, int xRes, int yRes, int size, float* values)
{
  memset(values,0,sizeof(float)*size); //clear values
  draw("2.ppm", extended, indices, colors, xRes, yRes, values, 0); //draw and output, mode = 0 (simple rasterization)
}

/* Problem 3: Orthographic Projection Matrix */
void p3(vector<VEC4>& extended, vector<VEC3>& vertices, vector<VEC3I>& indices, vector<VEC3>& colors,
const MATRIX4 Mvp, MATRIX4& Mortho, int xRes, int yRes, int size, float* values) 
{
  //reinitialize vars
  memset(values,0,sizeof(float)*size);
  extended.clear(); vertices.clear(); indices.clear(); colors.clear(); 
  
  buildBigSquare(vertices, indices, colors); //get big square
  extend_vecs(vertices, extended);

  //construct Mortho
  float l, r, t, b, n, f; l = b = f = 0.0; r = t = n = 12.0; //planes
  constructMortho(l, r, t, b, n, f, Mortho);

  //apply transforms
  transform(extended, Mortho); //apply Mortho
  transform(extended, Mvp); //apply Mvp

  //draw and output
  draw("3.ppm", extended, indices, colors, xRes, yRes, values, 0);
}

/* Problem 4: Camera Matrix */
void p4(vector<VEC4>& extended, vector<VEC3I>& indices, vector<VEC3>& colors, MATRIX4& Mcam, int xRes, int yRes, int size, float* values) 
{
  memset(values,0,sizeof(float)*size); //reinitialize some vars

  VEC3 eye(0.2, 0.2, 1.0); VEC3 lookAt(0,0,0); VEC3 up(0,1,0); //given vectors
  constructMcam(eye, lookAt, up, Mcam); //construct Mcam
  
  transform(extended, Mcam); //apply Mcam
  draw("4.ppm", extended, indices, colors, xRes, yRes, values, 0); //draw and output
}

/* Problem 5: Perspective Projection Matrix */
void p5(vector<VEC4>& extended, vector<VEC3>& vertices, vector<VEC3I>& indices, vector<VEC3>& colors, MATRIX4& Mvp, MATRIX4& Mcam, MATRIX4& Mper, 
int xRes, int yRes, int size, float* values)
{
  /* Reinitialize values */
  memset(values,0,sizeof(float)*size);
  extended.clear(); vertices.clear(); indices.clear(); colors.clear();
  
  buildCube(vertices, indices, colors); //get cube
  extend_vecs(vertices, extended);
  scale_vecs(extended, 0.5); //scale by factor of 0.5
  
  /* Apply transforms */

  //Mcam
  VEC3 eye(1.0, 1.0, 1.0); VEC3 lookAt(0.0,0.0,0.0); VEC3 up(0.0,1.0,0.0); //given vectors
  constructMcam(eye, lookAt, up, Mcam);
  transform(extended, Mcam);

  //Mper
  float fovy, aspect, n, f; n = 1.0; f = 100.0; fovy = 65.0; aspect = 4.0/3.0; //given values
  constructMper(fovy, aspect, n, f, Mper);
  transform(extended, Mper);

  //Ignore w for now -- set to 1.0
  vector<VEC3> truncated; vector<VEC4> verts;
  truncate_vecs(extended, truncated); 
  extend_vecs(truncated, verts);

  //Mvp
  transform(verts, Mvp);

  //draw and output
  draw("5.ppm", verts, indices, colors, xRes, yRes, values, 0);
}
  
/* Problem 6: Perspective Divide */
void p_divide(vector<VEC4>& extended)
{
  int count = extended.size(); float w = 0.0;
  for(int i = 0; i < count; i++) 
  {
    w = extended[i][3];
    extended[i] /= w; //scale each vector by 1/w
  }
}

void p6(vector<VEC4>& extended, vector<VEC3I>& indices, vector<VEC3>& colors, MATRIX4& Mvp, 
int xRes, int yRes, int size, float* values)
{
  memset(values,0,sizeof(float)*size); //reinitialize values array
  
  p_divide(extended);//perspective divide
  transform(extended, Mvp); //apply Mvp

  draw("6.ppm", extended, indices, colors, xRes, yRes, values, 0); //draw and output
}

/* Problem 7: Z-buffering */
void p7(vector<VEC4>& extended, vector<VEC3I>& indices, vector<VEC3>& colors, int xRes, int yRes, int size, float* values)
{
  memset(values,0,sizeof(float)*size); //clear values array
  draw("7.ppm", extended, indices, colors, xRes, yRes, values, 1); //draw and output, mode = 1 (z buffering)
}

/* Problem 8: Color Interpolation */
void p8(vector<VEC4>& extended, vector<VEC3>& vertices, vector<VEC3I>& indices, vector<VEC3>& colors,
int xRes, int yRes, int size, float* values)
{
  //reinitialize values
  memset(values,0,sizeof(float)*size);
  vertices.clear(); indices.clear(); colors.clear(); //extended.clear()

  buildCubePerVertexColors(vertices, indices, colors); //get new colors

  draw("8.ppm", extended, indices, colors, xRes, yRes, values, 2); //draw and output, mode = 2 (z-buffering + interpolation)
}

//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
int main(int argc, char* argv[])
{ 
  /* Declare & initialize vars */
  
  //dimensions of output image
  float nx = 800.0;
  float ny = 600.0;
  int xRes = (int) nx;
  int yRes = (int) ny;
  int size = xRes*yRes*3;

  //initalize array of values to zero
  float values[size]; memset(values,0,sizeof(float)*size); 

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
  memset(cam, 0, sizeof(float) * parameter_count); //by default camera parameters are set to zero

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

    buildCube(vertices, indices, colors); //get cube
    extend_vecs(vertices, extended);
    scale_vecs(extended, 0.5); //scale the vertices by 0.5

    //Mcam
    constructMcam(parameters[0], parameters[1], parameters[2], Mcam);
    transform(extended, Mcam);

    //Mper
    float fovy, aspect, n, f; n = 1.0; f = 100.0; fovy = 65.0; aspect = 4.0/3.0; //given values
    constructMper(fovy, aspect, n, f, Mper);
    transform(extended, Mper);  

    p_divide(extended); //perspective divide
    
    //Mvp
    constructMvp(nx,ny,Mvp);
    transform(extended, Mvp); 

    //draw and output
    draw("custom.ppm", extended, indices, colors, xRes, yRes, values, 2);
  }
  else
  {
    p1(extended, vertices, indices, colors, Mvp, nx, ny, values); //problem 1
    p2(extended, indices, colors, xRes, yRes, size, values); //problem 2
    p3(extended, vertices, indices, colors, Mvp, Mortho, xRes, yRes, size, values); //problem 3
    p4(extended, indices, colors, Mcam, xRes, yRes, size, values); //problem 4 
    p5(extended, vertices, indices, colors, Mvp, Mcam, Mper, xRes, yRes, size, values); //problem 5
    p6(extended, indices, colors, Mvp, xRes, yRes, size, values); //problem 6
    p7(extended, indices, colors, xRes, yRes, size, values); //problem 7
    p8(extended, vertices, indices, colors, xRes, yRes, size, values); //problem 8
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
void draw(const string& filename, vector<VEC4>& vertices, vector<VEC3I>& indices, vector<VEC3>& colors, int xRes, int yRes, float* values, int mode)
{
  // initialize z-buffer
  int size = xRes*yRes;
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
            index = (yRes-y)*xRes + x;
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

  writePPM(filename, xRes, yRes, values);
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

/* constructMvp()
 * Creates viewport matrix.
 */
void constructMvp(float nx, float ny, MATRIX4& Mvp)
{
  Mvp.setIdentity();

  Mvp(0,0) = nx/2.0; Mvp(0,3) = nx/2.0;
  Mvp(1,1) = ny/2.0; Mvp(1, 3) = ny/2.0;
}

/* constructMortho()
 * Creates orthographic projection matrix.
 */
void constructMortho(float l, float r, float t, float b, float n, float f, MATRIX4& Mortho)
{
  Mortho.setIdentity();

  Mortho(0,0) = 2.0/(r-l); Mortho(1,1) = 2.0/(t-b); Mortho(2,2) = -2.0/(n-f);
  Mortho(0,3) = -(r+l)/(r-l); Mortho(1,3) = -(t+b)/(t-b); Mortho(2,3) = -(n+f)/(n-f);
}

/* constructMcam()
 * Creates camera matrix.
 */
void constructMcam(VEC3 eye, VEC3 lookAt, VEC3 up, MATRIX4& Mcam)
{
  //initialize vars
  Mcam.setIdentity(); 
  VEC3 w; w.setZero(); VEC3 u; u.setZero(); VEC3 v; v.setZero();
  
  //calculate basis vectors
  VEC3 g; g = eye - lookAt;
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

/* constructMper()
 * Creates perspective projection matrix.
 */
void constructMper(float fovy, float aspect, float n, float f, MATRIX4& Mper)
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

/* color_vertices()
 * A simplified version of draw() that only colors the vertices of the triangles.
 */
void color_vertices(const vector<VEC4> vertices, vector<VEC3I>& indices, vector<VEC3>& colors, int xRes, int yRes, float* values)
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
      index = (yRes-y)*xRes*3 + 3*x;
      
      //color vertex
      values[index] = 255*colors[i][0]; 
      values[index+1] =255*colors[i][1];
      values[index+2] = 255*colors[i][2];
    }
  }
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

