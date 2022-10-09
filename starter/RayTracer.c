/*
  CSC D18 - RayTracer code.

  Written Dec. 9 2010 - Jan 20, 2011 by F. J. Estrada
  Freely distributable for adacemic purposes only.

  Uses Tom F. El-Maraghi's code for computing inverse
  matrices. You will need to compile together with
  svdDynamic.c

  You need to understand the code provided in
  this file, the corresponding header file, and the
  utils.c and utils.h files. Do not worry about
  svdDynamic.c, we need it only to compute
  inverse matrices.

  You only need to modify or add code in sections
  clearly marked "TO DO" - remember to check what
  functionality is actually needed for the corresponding
  assignment!

  Last updated: Aug. 2017   - F.J.E.
*/

/*****************************************************************************
 * COMPLETE THIS TEXT BOX:
 *
 * 1) Student Name: Bassel Ashi
 * 2) Student Name: Le Keping
 *
 * 1) Student number: 1005116731
 * 2) Student number: 1004948890
 *
 * 1) UtorID: ashipass
 * 2) UtorID: lekeping
 *
 * We hereby certify that the work contained here is our own
 *
 * ________Bassel_______            _________Le__________
 * (sign with your name)            (sign with your name)
 ********************************************************************************/

#include "utils.h"  // <-- This includes RayTracer.h

#define min(A, B) ((A) < (B) ? (A) : (B))

// A couple of global structures and data: An object list, a light list, and the
// maximum recursion depth
struct object3D *object_list;
struct pointLS *light_list;
struct textureNode *texture_list;
int MAX_DEPTH;

void buildScene(void) {
// #include "buildscene.c"  // <-- Import the scene definition!
#include "buildCoolScene.c"
}

void rtShade(struct object3D *obj, struct point3D *p, struct point3D *n,
             struct ray3D *ray, int depth, double a, double b,
             struct colourRGB *col) {
  // This function implements the shading model as described in lecture. It
  // takes
  // - A pointer to the first object intersected by the ray (to get the colour
  // properties)
  // - The coordinates of the intersection point (in world coordinates)
  // - The normal at the point
  // - The ray (needed to determine the reflection direction to use for the
  // global component, as well as for
  //   the Phong specular component)
  // - The current racursion depth
  // - The (a,b) texture coordinates (meaningless unless texture is enabled)
  //
  // Returns:
  // - The colour for this ray (using the col pointer)
  //

  struct colourRGB tmp_col;  // Accumulator for colour components
  double R, G, B;            // Colour for the object in R G and B

  // This will hold the colour as we process all the components of
  // the Phong illumination model
  tmp_col.R = 0;
  tmp_col.G = 0;
  tmp_col.B = 0;

  if (obj->texImg == NULL)  // Not textured, use object colour
  {
    R = obj->col.R;
    G = obj->col.G;
    B = obj->col.B;
  } else {
    // Get object colour from the texture given the texture coordinates (a,b),
    // and the texturing function for the object. Note that we will use textures
    // also for Photon Mapping.
    obj->textureMap(obj->texImg, a, b, &R, &G, &B);
  }

  //////////////////////////////////////////////////////////////
  // TO DO: Implement this function. Refer to the notes for
  // details about the shading model. (DONE)
  //////////////////////////////////////////////////////////////

  // Init some variables
  double lambda = -1;
  double tmp_a, tmp_b;
  struct object3D *tmp_obj;
  struct point3D tmp_p;
  struct point3D tmp_n;

  // Perform shading based on each light
  struct pointLS *curr_light = light_list;
  while (curr_light != NULL) {
    // Calculate p_to_ls (ray from p to light source)
    struct point3D p_to_ls_d;
    copyPoint(&curr_light->p0, &p_to_ls_d);
    subVectors(p, &p_to_ls_d);
    // normalize(&p_to_ls_d);
    struct ray3D p_to_ls;
    initRay(&p_to_ls, p, &p_to_ls_d);

    // Local component
    findFirstHit(&p_to_ls, &lambda, obj, &tmp_obj, &tmp_p, &tmp_n, &tmp_a,
                 &tmp_b);

    // Ambient
    tmp_col.R += R * obj->alb.ra;
    tmp_col.G += G * obj->alb.ra;
    tmp_col.B += B * obj->alb.ra;

    if (lambda <= 0 || lambda >= 1) {      
      // Diffuse
      struct point3D s;
      copyPoint(&p_to_ls_d, &s);
      normalize(&s);
      double n_dot_s = dot(n, &s);
      double n_dot_s_front_back = obj->frontAndBack ? abs(n_dot_s) : n_dot_s; 

      tmp_col.R += R * obj->alb.rd * curr_light->col.R * max(0, n_dot_s_front_back);
      tmp_col.G += G * obj->alb.rd * curr_light->col.G * max(0, n_dot_s_front_back);
      tmp_col.B += B * obj->alb.rd * curr_light->col.B * max(0, n_dot_s_front_back);
      
      // Specular
      struct point3D c;
      copyPoint(&ray->d, &c);
      c.px *= -1;
      c.py *= -1;
      c.pz *= -1;
      normalize(&c);
      
      
      struct point3D m;
      copyPoint(n, &m);
      m.px *= 2 * n_dot_s;
      m.py *= 2 * n_dot_s;
      m.pz *= 2 * n_dot_s;
      subVectors(&s, &m);
      
      double c_dot_m = dot(&c, &m);

      tmp_col.R += obj->alb.rs * curr_light->col.R * pow(max(0, c_dot_m), obj->shinyness);
      tmp_col.G += obj->alb.rs * curr_light->col.G * pow(max(0, c_dot_m), obj->shinyness);
      tmp_col.B += obj->alb.rs * curr_light->col.B * pow(max(0, c_dot_m), obj->shinyness);
    }

    // Global component
    if (depth < MAX_DEPTH) {
      // Specular reflection
      if (obj->alb.rs != 0) {
        // 2n.d
        double two_n_dot_d = 2 * dot(n, &ray->d);

        // (2n.d)n
        struct point3D ms_d;
        copyPoint(n, &ms_d);
        ms_d.px *= two_n_dot_d;
        ms_d.py *= two_n_dot_d;
        ms_d.pz *= two_n_dot_d;
        ms_d.pw = 1;

        // d - (2n.d)n
        subVectors(&ray->d, &ms_d);
        ms_d.px *= -1;
        ms_d.py *= -1;
        ms_d.pz *= -1;
        normalize(&ms_d);

        // Build ms ray
        struct ray3D ms;
        initRay(&ms, p, &ms_d);

        // Perform raytracing using ms
        struct colourRGB c;
        rayTrace(&ms, depth + 1, &c, obj);

        // Update pixel colour
        tmp_col.R += c.R * obj->alb.rg;
        tmp_col.G += c.G * obj->alb.rg;
        tmp_col.B += c.B * obj->alb.rg;
      }

      // Refraction
      if (obj->alpha < 1) {
        // TODO: A3
      }
    }

    curr_light = curr_light->next;
  }

  // Update colour
  col->R = min(tmp_col.R, 1);
  col->G = min(tmp_col.G, 1);
  col->B = min(tmp_col.B, 1);
}

void findFirstHit(struct ray3D *ray, double *lambda, struct object3D *Os,
                  struct object3D **obj, struct point3D *p, struct point3D *n,
                  double *a, double *b) {
  // Find the closest intersection between the ray and any objects in the scene.
  // Inputs:
  //   *ray    -  A pointer to the ray being traced
  //   *Os     -  'Object source' is a pointer toward the object from which the
  //   ray originates. It is used for reflected or refracted rays
  //              so that you can check for and ignore self-intersections as
  //              needed. It is NULL for rays originating at the center of
  //              projection
  // Outputs:
  //   *lambda -  A pointer toward a double variable 'lambda' used to return the
  //   lambda at the intersection point
  //   **obj   -  A pointer toward an (object3D *) variable so you can return a
  //   pointer to the object that has the closest intersection with
  //              this ray (this is required so you can do the shading)
  //   *p      -  A pointer to a 3D point structure so you can store the
  //   coordinates of the intersection point *n      -  A pointer to a 3D point
  //   structure so you can return the normal at the intersection point *a, *b
  //   -  Pointers toward double variables so you can return the texture
  //   coordinates a,b at the intersection point

  /////////////////////////////////////////////////////////////
  // TO DO: Implement this function. See the notes for
  // reference of what to do in here (DONE)
  /////////////////////////////////////////////////////////////

  double curr_lambda = -1;
  struct point3D curr_p, curr_n;
  double curr_a, curr_b;

  struct object3D *curr_obj = object_list;

  // Loop on objects to find closest one for intersection
  while (curr_obj != NULL) {
    if (curr_obj == Os) {
      curr_obj = curr_obj->next;
      continue;
    }

    // Get ray intersection with current object
    curr_obj->intersect(curr_obj, ray, &curr_lambda, &curr_p, &curr_n, &curr_a,
                        &curr_b);

    // Update variables if current object is the new closest
    if ((*lambda < 0 || curr_lambda < *lambda) && curr_lambda > 0) {
      *lambda = curr_lambda;
      *p = curr_p;
      *n = curr_n;
      *a = curr_a;
      *b = curr_b;
      *obj = curr_obj;
    }

    curr_obj = curr_obj->next;
  }
}

void rayTrace(struct ray3D *ray, int depth, struct colourRGB *col,
              struct object3D *Os) {
  // Trace one ray through the scene.
  //
  // Parameters:
  //   *ray   -  A pointer to the ray being traced
  //   depth  -  Current recursion depth for recursive raytracing
  //   *col   - Pointer to an RGB colour structure so you can return the object
  //   colour
  //            at the intersection point of this ray with the closest scene
  //            object.
  //   *Os    - 'Object source' is a pointer to the object from which the ray
  //            originates so you can discard self-intersections due to
  //            numerical errors. NULL for rays originating from the center of
  //            projection.

  double lambda = -1;    // Lambda at intersection
  double a, b;           // Texture coordinates
  struct object3D *obj;  // Pointer to object at intersection
  struct point3D p;      // Intersection point
  struct point3D n;      // Normal at intersection
  struct colourRGB I;    // Colour returned by shading function

  if (depth > MAX_DEPTH)  // Max recursion depth reached. Return invalid colour.
  {
    col->R = -1;
    col->G = -1;
    col->B = -1;
    return;
  }

  ///////////////////////////////////////////////////////
  // TO DO: Complete this function. Refer to the notes
  // if you are unsure what to do here. (DONE)
  ///////////////////////////////////////////////////////

  // Find the first intersection of ray in the scene
  findFirstHit(ray, &lambda, Os, &obj, &p, &n, &a, &b);

  // Perform shading only if ray length is above zero
  if (lambda > 0.0) {
    rtShade(obj, &p, &n, ray, depth, a, b, &I);

    col->R = I.R;
    col->G = I.G;
    col->B = I.B;
  } else {
    // Background colour
    col->R = 0;
    col->G = 0;
    col->B = 0;
  }
}

int main(int argc, char *argv[]) {
  // Main function for the raytracer. Parses input parameters,
  // sets up the initial blank image, and calls the functions
  // that set up the scene and do the raytracing.
  struct image *im;  // Will hold the raytraced image
  struct view *cam;  // Camera and view for this scene
  int sx;            // Size of the raytraced image
  int antialiasing;  // Flag to determine whether antialiaing is enabled or
                     // disabled
  char output_name[1024];  // Name of the output file for the raytraced .ppm
                           // image
  struct point3D e;        // Camera view parameters 'e', 'g', and 'up'
  struct point3D g;
  struct point3D up;
  double du, dv;  // Increase along u and v directions for pixel coordinates
  struct point3D pc, d;  // Point structures to keep the coordinates of a pixel
                         // and the direction or a ray
  struct ray3D ray;      // Structure to keep the ray from e to a pixel
  struct colourRGB col;  // Return colour for raytraced pixels
  struct colourRGB background;  // Background colour
  int i, j;                     // Counters for pixel coordinates
  unsigned char *rgbIm;

  if (argc < 5) {
    fprintf(stderr, "RayTracer: Can not parse input parameters\n");
    fprintf(stderr, "USAGE: RayTracer size rec_depth antialias output_name\n");
    fprintf(stderr, "   size = Image size (both along x and y)\n");
    fprintf(stderr, "   rec_depth = Recursion depth\n");
    fprintf(stderr,
            "   antialias = A single digit, 0 disables antialiasing. Anything "
            "else enables antialiasing\n");
    fprintf(stderr,
            "   output_name = Name of the output file, e.g. MyRender.ppm\n");
    exit(0);
  }
  sx = atoi(argv[1]);
  MAX_DEPTH = atoi(argv[2]);
  if (atoi(argv[3]) == 0)
    antialiasing = 0;
  else
    antialiasing = 1;
  strcpy(&output_name[0], argv[4]);

  fprintf(stderr, "Rendering image at %d x %d\n", sx, sx);
  fprintf(stderr, "Recursion depth = %d\n", MAX_DEPTH);
  if (!antialiasing)
    fprintf(stderr, "Antialising is off\n");
  else
    fprintf(stderr, "Antialising is on\n");
  fprintf(stderr, "Output file name: %s\n", output_name);

  object_list = NULL;
  light_list = NULL;
  texture_list = NULL;

  // Allocate memory for the new image
  im = newImage(sx, sx);
  if (!im) {
    fprintf(stderr, "Unable to allocate memory for raytraced image\n");
    exit(0);
  } else
    rgbIm = (unsigned char *)im->rgbdata;

  ///////////////////////////////////////////////////
  // TO DO: You will need to implement several of the
  //        functions below. For Assignment 2, you can use
  //        the simple scene already provided. But
  //        for Assignment 3 you need to create your own
  //        *interesting* scene.
  ///////////////////////////////////////////////////
  buildScene();  // Create a scene. This defines all the
                 // objects in the world of the raytracer

  //////////////////////////////////////////
  // TO DO: For Assignment 2 you can use the setup
  //        already provided here. For Assignment 3
  //        you may want to move the camera
  //        and change the view parameters
  //        to suit your scene.
  //////////////////////////////////////////

  // Mind the homogeneous coordinate w of all vectors below. DO NOT
  // forget to set it to 1, or you'll get junk out of the
  // geometric transformations later on.

  // Camera center is at (0,0,-1)
  e.px = 0;
  e.py = 0;
  e.pz = -1;
  e.pw = 1;

  // To define the gaze vector, we choose a point 'pc' in the scene that
  // the camera is looking at, and do the vector subtraction pc-e.
  // Here we set up the camera to be looking at the origin.
  g.px = 0 - e.px;
  g.py = 0 - e.py;
  g.pz = 0 - e.pz;
  g.pw = 1;
  // In this case, the camera is looking along the world Z axis, so
  // vector w should end up being [0, 0, -1]

  // Define the 'up' vector to be the Y axis
  up.px = 0;
  up.py = 1;
  up.pz = 0;
  up.pw = 1;

  // Set up view with given the above vectors, a 4x4 window,
  // and a focal length of -1 (why? where is the image plane?)
  // Note that the top-left corner of the window is at (-2, 2)
  // in camera coordinates.
  cam = setupView(&e, &g, &up, -1, -2, 2, 4);

  if (cam == NULL) {
    fprintf(
        stderr,
        "Unable to set up the view and camera parameters. Our of memory!\n");
    cleanup(object_list, light_list, texture_list);
    deleteImage(im);
    exit(0);
  }

  // Set up background colour here
  background.R = 0;
  background.G = 0;
  background.B = 0;

  // Do the raytracing
  //////////////////////////////////////////////////////
  // TO DO: You will need code here to do the raytracing
  //        for each pixel in the image. Refer to the
  //        lecture notes, in particular, to the
  //        raytracing pseudocode, for details on what
  //        to do here. Make sure you undersand the
  //        overall procedure of raytracing for a single
  //        pixel.
  //////////////////////////////////////////////////////
  du = cam->wsize /
       (sx - 1);  // du and dv. In the notes in terms of wl and wr, wt and wb,
  dv = -cam->wsize /
       (sx - 1);  // here we use wl, wt, and wsize. du=dv since the image is
                  // and dv is negative since y increases downward in pixel
                  // coordinates and upward in camera coordinates.

  fprintf(stderr, "View parameters:\n");
  fprintf(stderr, "Left=%f, Top=%f, Width=%f, f=%f\n", cam->wl, cam->wt,
          cam->wsize, cam->f);
  fprintf(stderr,
          "Camera to world conversion matrix (make sure it makes sense!):\n");
  printmatrix(cam->C2W);
  fprintf(stderr, "World to camera conversion matrix:\n");
  printmatrix(cam->W2C);
  fprintf(stderr, "\n");

  fprintf(stderr, "Rendering row: ");
  for (j = 0; j < sx; j++)  // For each of the pixels in the image
  {
    fprintf(stderr, "%d/%d, ", j, sx);
    for (i = 0; i < sx; i++) {
      ///////////////////////////////////////////////////////////////////
      // TO DO - complete the code that should be in this loop to do the
      //         raytracing! (DONE)
      ///////////////////////////////////////////////////////////////////

      // Compute point pij
      pc.px = cam->wl + i * du;
      pc.py = cam->wt + j * dv;
      pc.pz = -1;
      pc.pw = 1;
      matVecMult(cam->C2W, &pc);

      // Compute direction d
      d.px = cam->wl + i * du;
      d.py = cam->wt + j * dv;
      d.pz = -1;
      d.pw = 0;
      matVecMult(cam->C2W, &d);
      normalize(&d);

      // Compute ray
      initRay(&ray, &pc, &d);
      col = background;

      // Perform raytracing
      rayTrace(&ray, 0, &col, NULL);

      // Calculate current pixel index
      int curr_pixel = (i + j * sx) * 3;

      // Update pixels in image
      *(rgbIm + curr_pixel) = (unsigned char)(col.R * 255.0);
      *(rgbIm + curr_pixel + 1) = (unsigned char)(col.G * 255.0);
      *(rgbIm + curr_pixel + 2) = (unsigned char)(col.B * 255.0);
    }  // end for i
  }    // end for j

  fprintf(stderr, "\nDone!\n");

  // Output rendered image
  imageOutput(im, output_name);

  // Exit section. Clean up and return.
  cleanup(object_list, light_list,
          texture_list);  // Object, light, and texture lists
  deleteImage(im);        // Rendered image
  free(cam);              // camera view
  exit(0);
}
