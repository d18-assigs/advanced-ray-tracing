/*
   utils.c - F.J. Estrada, Dec. 9, 2010

   Utilities for the ray tracer. You will need to complete
   some of the functions in this file. Look for the sections
   marked "TO DO". Be sure to read the rest of the file and
   understand how the entire code works.

   HOWEVER: Note that there are a lot of incomplete functions
            that will only be used for the advanced ray tracer!
            So, read the handout carefully and implement only
            the code you need for the corresponding assignment.

   Last updated: Aug. 2017  -  F.J.E.
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

#include "utils.h"

double smallest_cube_size;
// A useful 4x4 identity matrix which can be used at any point to
// initialize or reset object transformations
double eye4x4[4][4] = {{1.0, 0.0, 0.0, 0.0},
                       {0.0, 1.0, 0.0, 0.0},
                       {0.0, 0.0, 1.0, 0.0},
                       {0.0, 0.0, 0.0, 1.0}};

/////////////////////////////////////////////
// Primitive data structure section
/////////////////////////////////////////////
struct point3D *newPoint(double px, double py, double pz) {
    // Allocate a new point structure, initialize it to
    // the specified coordinates, and return a pointer
    // to it.

    struct point3D *pt = (struct point3D *)calloc(1, sizeof(struct point3D));
    if (!pt)
        fprintf(stderr, "Out of memory allocating point structure!\n");
    else {
        pt->px = px;
        pt->py = py;
        pt->pz = pz;
        pt->pw = 1.0;
    }
    return (pt);
}

struct pointLS *newPLS(struct point3D *p0, double r, double g, double b) {
    // Allocate a new point light sourse structure. Initialize the light
    // source to the specified RGB colour
    // Note that this is a point light source in that it is a single point
    // in space, if you also want a uniform direction for light over the
    // scene (a so-called directional light) you need to place the
    // light source really far away.

    struct pointLS *ls = (struct pointLS *)calloc(1, sizeof(struct pointLS));
    if (!ls)
        fprintf(stderr, "Out of memory allocating light source!\n");
    else {
        memcpy(&ls->p0, p0, sizeof(struct point3D)); // Copy light source location

        ls->col.R = r; // Store light source colour and
        ls->col.G = g; // intensity
        ls->col.B = b;
    }
    return (ls);
}

void copyPoint(struct point3D *from, struct point3D *to) {
    to->px = from->px;
    to->py = from->py;
    to->pz = from->pz;
    to->pw = from->pw;
}

/////////////////////////////////////////////
// Ray and normal transforms
/////////////////////////////////////////////
inline void rayTransform(struct ray3D *ray_orig, struct ray3D *ray_transformed,
                         struct object3D *obj) {
    // Transforms a ray using the inverse transform for the specified object. This
    // is so that we can use the intersection test for the canonical object. Note
    // that this has to be done carefully!

    ///////////////////////////////////////////
    // TO DO: Complete this function (DONE)
    ///////////////////////////////////////////

    // Start with original ray
    *ray_transformed = *ray_orig;

    // Set homogeneous coordinate to zero for proper calculations
    ray_transformed->d.pw = 0;

    // Update the initial point and direction of the new ray by the inverse
    // transformation contained in obj
    matVecMult(obj->Tinv, &ray_transformed->p0);
    matVecMult(obj->Tinv, &ray_transformed->d);
}

inline void normalTransform(struct point3D *n_orig,
                            struct point3D *n_transformed,
                            struct object3D *obj) {
    // Computes the normal at an affinely transformed point given the original
    // normal and the object's inverse transformation. From the notes:
    // n_transformed=A^-T*n normalized.

    ///////////////////////////////////////////
    // TO DO: Complete this function (DONE)
    ///////////////////////////////////////////

    // Start with original normal
    *n_transformed = *n_orig;

    // Compute A inverse transpose
    double A_inv_trans[4][4];
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            A_inv_trans[i][j] = obj->Tinv[j][i];
        }
    }

    // Set homogeneous coordinate to zero for proper calculations
    n_transformed->pw = 0;

    // Computer A^-T
    matVecMult(A_inv_trans, n_transformed);

    // Normalize
    normalize(n_transformed);

    // Revert homogeneous coordinates update
    n_transformed->pw = 1;
}

/////////////////////////////////////////////
// Object management section
/////////////////////////////////////////////
void insertObject(struct object3D *o, struct object3D **list) {
    if (o == NULL)
        return;
    // Inserts an object into the object list.
    if (*(list) == NULL) {
        *(list) = o;
        (*(list))->next = NULL;
    } else {
        o->next = (*(list))->next;
        (*(list))->next = o;
    }
}
struct object3D *newPlane(double ra, double rd, double rs, double rg, double r,
                          double g, double b, double alpha, double r_index,
                          double shiny) {
    // Intialize a new plane with the specified parameters:
    // ra, rd, rs, rg - Albedos for the components of the Phong model
    // r, g, b, - Colour for this plane
    // alpha - Transparency, must be set to 1 unless you are doing refraction
    // r_index - Refraction index if you are doing refraction.
    // shiny - Exponent for the specular component of the Phong model
    //
    // The plane is defined by the following vertices (CCW)
    // (1,1,0), (-1,1,0), (-1,-1,0), (1,-1,0)
    // With normal vector (0,0,1) (i.e. parallel to the XY plane)

    struct object3D *plane =
        (struct object3D *)calloc(1, sizeof(struct object3D));

    if (!plane)
        fprintf(stderr, "Unable to allocate new plane, out of memory!\n");
    else {
        plane->alb.ra = ra;
        plane->alb.rd = rd;
        plane->alb.rs = rs;
        plane->alb.rg = rg;
        plane->col.R = r;
        plane->col.G = g;
        plane->col.B = b;
        plane->alpha = alpha;
        plane->r_index = r_index;
        plane->shinyness = shiny;
        plane->intersect = &planeIntersect;
        plane->surfaceCoords = &planeCoordinates;
        plane->randomPoint = &planeSample;
        plane->texImg = NULL;
        plane->photonMap = NULL;
        plane->normalMap = NULL;
        memcpy(&plane->T[0][0], &eye4x4[0][0], 16 * sizeof(double));
        memcpy(&plane->Tinv[0][0], &eye4x4[0][0], 16 * sizeof(double));
        plane->textureMap = &texMap;
        plane->alphaMap = &alphaMap;
        plane->frontAndBack = 1;
        plane->photonMapped = 0;
        plane->normalMapped = 0;
        plane->isCSG = 0;
        plane->isLightSource = 0;
        plane->CSGnext = NULL;
        plane->next = NULL;
    }
    return (plane);
}

void getPMinPmax(struct object3D *obj, struct point3D *p_min, struct point3D *p_max) {
    struct point3D tmp_min, tmp_max;
    assignPoint(&tmp_max, -__DBL_MAX__, -__DBL_MAX__, -__DBL_MAX__);
    assignPoint(&tmp_min, __DBL_MAX__, __DBL_MAX__, __DBL_MAX__);
    for (size_t i = 0; i < NUM_SAMPLE_OCT_TREE; i++) {
        double x, y, z;
        obj->randomPoint(obj, &x, &y, &z);
        assignPoint(&tmp_min, min(tmp_min.px, x), min(tmp_min.py, y), min(tmp_min.pz, z));
        assignPoint(&tmp_max, max(tmp_max.px, x), max(tmp_max.py, y), max(tmp_max.pz, z));
    }

    copyPoint(&tmp_min, p_min);
    copyPoint(&tmp_max, p_max);
}

int haveOverLap(struct object3D *obj, struct point3D *p_min_box, struct point3D *p_max_box) {
    struct point3D tmp_min, tmp_max;
    int x_overlap, y_overlap, z_overlap;
    double x, y, z;
    int obj_in_box;
    int box_in_obj;
    for (int i = 0; i < NUM_SAMPLE_OCT_TREE; i++) {
        obj->randomPoint(obj, &x, &y, &z);
        assignPoint(&tmp_min, min(tmp_min.px, x), min(tmp_min.py, y), min(tmp_min.pz, z));
        assignPoint(&tmp_max, max(tmp_max.px, x), max(tmp_max.py, y), max(tmp_max.pz, z));
        x_overlap = x > p_min_box->px && x < p_max_box->px;
        y_overlap = y > p_min_box->py && y < p_max_box->py;
        z_overlap = z > p_min_box->pz && z < p_max_box->pz;
        obj_in_box = x_overlap && y_overlap && z_overlap;
        if (obj->overlap != NULL) {
            struct point3D tmp_min, tmp_max;
            int x_overlap, y_overlap, z_overlap;
            double x, y, z;
            x = get_random_double(p_min_box->px, p_max_box->px);
            y = get_random_double(p_min_box->py, p_max_box->py);
            z = get_random_double(p_min_box->pz, p_max_box->pz);
            box_in_obj = obj->overlap(obj, x, y, z);
            if (obj_in_box || box_in_obj) {
                return 1;
            }
        } else {
            if (obj_in_box) {
                return 1;
            }
        }
    }
    return 0;
}
struct object3D *newSphere(double ra, double rd, double rs, double rg, double r,
                           double g, double b, double alpha, double r_index,
                           double shiny) {
    // Intialize a new sphere with the specified parameters:
    // ra, rd, rs, rg - Albedos for the components of the Phong model
    // r, g, b, - Colour for this plane
    // alpha - Transparency, must be set to 1 unless you are doing refraction
    // r_index - Refraction index if you are doing refraction.
    // shiny -Exponent for the specular component of the Phong model
    //
    // This is assumed to represent a unit sphere centered at the origin.
    //

    struct object3D *sphere =
        (struct object3D *)calloc(1, sizeof(struct object3D));

    if (!sphere)
        fprintf(stderr, "Unable to allocate new sphere, out of memory!\n");
    else {
        sphere->alb.ra = ra;
        sphere->alb.rd = rd;
        sphere->alb.rs = rs;
        sphere->alb.rg = rg;
        sphere->col.R = r;
        sphere->col.G = g;
        sphere->col.B = b;
        sphere->alpha = alpha;
        sphere->r_index = r_index;
        sphere->shinyness = shiny;
        sphere->intersect = &sphereIntersect;
        sphere->surfaceCoords = &sphereCoordinates;
        sphere->randomPoint = &sphereSample;
        sphere->overlap = &sphereOverlap;
        sphere->texImg = NULL;
        sphere->photonMap = NULL;
        sphere->normalMap = NULL;
        memcpy(&sphere->T[0][0], &eye4x4[0][0], 16 * sizeof(double));
        memcpy(&sphere->Tinv[0][0], &eye4x4[0][0], 16 * sizeof(double));
        sphere->textureMap = &texMap;
        sphere->alphaMap = &alphaMap;
        sphere->frontAndBack = 0;
        sphere->photonMapped = 0;
        sphere->normalMapped = 0;
        sphere->isCSG = 0;
        sphere->isLightSource = 0;
        sphere->CSGnext = NULL;
        sphere->next = NULL;
    }
    return (sphere);
}

struct object3D *newCyl(double ra, double rd, double rs, double rg, double r,
                        double g, double b, double alpha, double r_index,
                        double shiny) {
    ///////////////////////////////////////////////////////////////////////////////////////
    // TO DO (DONE):
    //	Complete the code to create and initialize a new cylinder object.
    ///////////////////////////////////////////////////////////////////////////////////////

    struct object3D *cylinder =
        (struct object3D *)calloc(1, sizeof(struct object3D));

    if (!cylinder)
        fprintf(stderr, "Unable to allocate new cylinder, out of memory!\n");
    else {
        cylinder->alb.ra = ra;
        cylinder->alb.rd = rd;
        cylinder->alb.rs = rs;
        cylinder->alb.rg = rg;
        cylinder->col.R = r;
        cylinder->col.G = g;
        cylinder->col.B = b;
        cylinder->alpha = alpha;
        cylinder->r_index = r_index;
        cylinder->shinyness = shiny;
        cylinder->intersect = &cylIntersect;
        cylinder->surfaceCoords = &cylCoordinates;
        cylinder->randomPoint = &cylSample;
        cylinder->overlap = &cylOverlap;
        cylinder->texImg = NULL;
        cylinder->photonMap = NULL;
        cylinder->normalMap = NULL;
        memcpy(&cylinder->T[0][0], &eye4x4[0][0], 16 * sizeof(double));
        memcpy(&cylinder->Tinv[0][0], &eye4x4[0][0], 16 * sizeof(double));
        cylinder->textureMap = &texMap;
        cylinder->alphaMap = &alphaMap;
        cylinder->frontAndBack = 0;
        cylinder->photonMapped = 0;
        cylinder->normalMapped = 0;
        cylinder->isCSG = 0;
        cylinder->isLightSource = 0;
        cylinder->CSGnext = NULL;
        cylinder->next = NULL;
    }
    return (cylinder);
}

struct object3D *newCone(double ra, double rd, double rs, double rg, double r,
                         double g, double b, double alpha, double r_index,
                         double shiny) {
    ///////////////////////////////////////////////////////////////////////////////////////
    // TO DO (DONE):
    //	Complete the code to create and initialize a new cylinder object.
    ///////////////////////////////////////////////////////////////////////////////////////

    struct object3D *cone =
        (struct object3D *)calloc(1, sizeof(struct object3D));

    if (!cone)
        fprintf(stderr, "Unable to allocate new cylinder, out of memory!\n");
    else {
        cone->alb.ra = ra;
        cone->alb.rd = rd;
        cone->alb.rs = rs;
        cone->alb.rg = rg;
        cone->col.R = r;
        cone->col.G = g;
        cone->col.B = b;
        cone->alpha = alpha;
        cone->r_index = r_index;
        cone->shinyness = shiny;
        cone->intersect = &coneIntersect;
        cone->surfaceCoords = &coneCoordinates;
        cone->randomPoint = &coneSample;
        cone->overlap = &coneOverlap;
        cone->texImg = NULL;
        cone->photonMap = NULL;
        cone->normalMap = NULL;
        memcpy(&cone->T[0][0], &eye4x4[0][0], 16 * sizeof(double));
        memcpy(&cone->Tinv[0][0], &eye4x4[0][0], 16 * sizeof(double));
        cone->textureMap = &texMap;
        cone->alphaMap = &alphaMap;
        cone->frontAndBack = 0;
        cone->photonMapped = 0;
        cone->normalMapped = 0;
        cone->isCSG = 0;
        cone->isLightSource = 0;
        cone->CSGnext = NULL;
        cone->next = NULL;
    }
    return (cone);
}

struct object3D *newOutermostOctCube(struct point3D *p_min, struct point3D *p_max) {

    double dis_x, dis_y, dis_z;
    dis_x = p_min->px - FIRST_CUBE_MARGINS;
    dis_y = p_min->py - FIRST_CUBE_MARGINS;
    dis_z = p_min->pz - FIRST_CUBE_MARGINS;

    double scale;
    scale = (p_max->px - dis_x) / 2;
    scale = max((p_max->py - dis_x) / 2, scale);
    scale = max((p_max->pz - dis_x) / 2, scale);

    smallest_cube_size = scale * pow(0.5,OCT_TREE_DEPTH - 1);

    struct object3D *cube = newCube(.05, .95, .35, .35, .5, .5, .5, .5, .6, 1);
    Translate(cube, 1, 1, 1);
    Scale(cube, scale, scale, scale);
    Translate(cube, dis_x, dis_y, dis_z);
    invert(&cube->T[0][0], &cube->Tinv[0][0]);
    return cube;
}
struct object3D *newCube(double ra, double rd, double rs, double rg, double r, double g, double b, double alpha, double r_index, double shiny) {
    struct object3D *cube = (struct object3D *)calloc(1, sizeof(struct object3D));
    if (!cube)
        fprintf(stderr, "Unable to allocate new sphere, out of memory!\n");
    else {
        cube->alb.ra = ra;
        cube->alb.rd = rd;
        cube->alb.rs = rs;
        cube->alb.rg = rg;
        cube->col.R = r;
        cube->col.G = g;
        cube->col.B = b;
        cube->alpha = alpha;
        cube->r_index = r_index;
        cube->shinyness = shiny;
        cube->intersect = &cubeIntersect;
        cube->surfaceCoords = NULL;
        // cube->randomPoint = &cubeSample;
        cube->texImg = NULL;
        cube->photonMap = NULL;
        cube->normalMap = NULL;
        memcpy(&cube->T[0][0], &eye4x4[0][0], 16 * sizeof(double));
        memcpy(&cube->Tinv[0][0], &eye4x4[0][0], 16 * sizeof(double));
        cube->textureMap = &texMap;
        cube->frontAndBack = 0;
        cube->photonMapped = 0;
        cube->normalMapped = 0;
        cube->isCSG = 0;
        cube->isLightSource = 0;
        cube->CSGnext = NULL;
        cube->next = NULL;
    }
    return (cube);
}

struct object3D *duplicateObj(struct object3D *src_obj, int with_T) {
    // return newCube(src->alb.ra,src->alb.rd,src->alb.rs,src->alb.rg,src->col.R,src->col.G,src->col.B,src->alpha,src->r_index,src->shinyness);
    struct object3D *newObject = (struct object3D *)calloc(1, sizeof(struct object3D));
    if (!newObject)
        fprintf(stderr, "Unable to allocate new object, out of memory!\n");
    else {
        newObject->alb.ra = src_obj->alb.ra;
        newObject->alb.rd = src_obj->alb.rd;
        newObject->alb.rs = src_obj->alb.rs;
        newObject->alb.rg = src_obj->alb.rg;
        newObject->col.R = src_obj->col.R;
        newObject->col.G = src_obj->col.G;
        newObject->col.B = src_obj->col.B;
        newObject->alpha = src_obj->alpha;
        newObject->r_index = src_obj->r_index;
        newObject->shinyness = src_obj->shinyness;
        newObject->intersect = src_obj->intersect;
        newObject->surfaceCoords = src_obj->surfaceCoords;
        newObject->overlap = src_obj->overlap;
        newObject->randomPoint = src_obj->randomPoint;
        newObject->texImg = src_obj->texImg;
        newObject->photonMap = src_obj->photonMap;
        newObject->normalMap = src_obj->normalMap;

        if (with_T) {

            memcpy(&newObject->T[0][0], &src_obj->T[0][0], 16 * sizeof(double));
            memcpy(&newObject->Tinv[0][0], &src_obj->Tinv[0][0], 16 * sizeof(double));
        } else {

            memcpy(&newObject->T[0][0], &eye4x4[0][0], 16 * sizeof(double));
            memcpy(&newObject->Tinv[0][0], &eye4x4[0][0], 16 * sizeof(double));
        }

        newObject->textureMap = src_obj->textureMap;
        newObject->frontAndBack = src_obj->frontAndBack;
        newObject->photonMapped = src_obj->photonMapped;
        newObject->normalMapped = src_obj->normalMapped;
        newObject->isCSG = src_obj->isCSG;
        newObject->isLightSource = src_obj->isLightSource;
        newObject->CSGnext = NULL;
        newObject->next = NULL;
    }
    return (newObject);
}

struct areaLS *newALS(struct object3D *obj, int sample) {
    // construct new als using the shape of obj and the color of obj as light rgb

    struct areaLS *aLS =
        (struct areaLS *)calloc(1, sizeof(struct areaLS));
    aLS->light_shape = obj;
    aLS->k_sample = sample;
    obj->isLightSource = 1;
    return aLS;
}

///////////////////////////////////////////////////////////////////////////////////////
// TO DO:
//	Complete the functions that compute intersections for the canonical
// plane
//      and canonical sphere with a given ray. This is the most fundamental
//      component of the raytracer.
///////////////////////////////////////////////////////////////////////////////////////
void planeIntersect(struct object3D *plane, struct ray3D *ray, double *lambda,
                    struct point3D *p, struct point3D *n, double *a,
                    double *b) {
    // Computes and returns the value of 'lambda' at the intersection
    // between the specified ray and the specified canonical plane.

    /////////////////////////////////
    // TO DO: Complete this function. (DONE)
    /////////////////////////////////

    // The following calculations are based on the intersection defined in p.75
    //    Lambda = ((p1 - a) . n) / (d . n)

    struct ray3D ray_transformed;
    struct point3D normal, p1a;
    rayTransform(ray, &ray_transformed, plane);
    normal.px = 0;
    normal.py = 0;
    normal.pz = 1;
    normal.pw = 1;
    p1a.px = 0;
    p1a.py = 0;
    p1a.pz = 0;
    p1a.pw = 1;

    // Ensure ray is not parallel z-axis
    if (dot(&ray_transformed.d, &normal) <= 0)
        return;

    // Calculate current lambda

    subVectors(&ray_transformed.p0, &p1a);
    double curr_lambda = dot(&p1a, &normal) / dot(&ray_transformed.d, &normal);

    // Ensure current lambda implies intersection
    if (curr_lambda == 0)
        return;

    rayPosition(&ray_transformed, curr_lambda, p);

    // Ensure it occurrs within the plane
    if (p->px < 1 && p->px > -1 && p->py < 1 && p->py > -1) {
        *lambda = curr_lambda;
        *a = - (p->px - 1) / 2;
        *b = - (p->py - 1) / 2;

        normalTransform(&normal, n, plane);
        normalize(n);

        ray->rayPos(ray, *lambda, p);
    }
}

inline void getCubeNormal(struct point3D *p, struct point3D *n) {
    double x = p->px;
    double y = p->py;
    double z = p->pz;

    if (fabs(x) > fabs(y) && fabs(x) > fabs(z)) {
        assignPoint(n, x, 0, 0);
    } else if (fabs(y) > fabs(z) && fabs(y) > fabs(x)) {
        assignPoint(n, 0, y, 0);
    } else {
        assignPoint(n, 0, 0, z);
    }
}
void cubeIntersect(struct object3D *cube, struct ray3D *ray, double *lambda, struct point3D *p, struct point3D *n, double *a, double *b) {

    // cubeIntersect reference is at
    // https://gamedev.stackexchange.com/questions/18436/most-efficient-aabb-vs-ray-collision-algorithms

    struct ray3D ray_transformed;
    rayTransform(ray, &ray_transformed, cube);
    struct point3D dirfrac;
    assignPoint(&dirfrac, 1.0 / (ray_transformed.d.px), 1.0 / (ray_transformed.d.py), 1.0 / (ray_transformed.d.pz));
    double t1, t2, t3, t4, t5, t6, tmin, tmax;

    t1 = (-1 - ray_transformed.p0.px) * dirfrac.px;
    t2 = (1 - ray_transformed.p0.px) * dirfrac.px;
    t3 = (-1 - ray_transformed.p0.py) * dirfrac.py;
    t4 = (1 - ray_transformed.p0.py) * dirfrac.py;
    t5 = (-1 - ray_transformed.p0.pz) * dirfrac.pz;
    t6 = (1 - ray_transformed.p0.pz) * dirfrac.pz;

    tmin = max(max(min(t1, t2), min(t3, t4)), min(t5, t6));
    tmax = min(min(max(t1, t2), max(t3, t4)), max(t5, t6));

    if (tmax > 0 && tmin < tmax) {
        if (tmin < 0) {
            *lambda = tmax;
        } else
            *lambda = tmin;
        struct point3D tmp_p, tmp_n;
        rayPosition(&ray_transformed, *lambda, &tmp_p);

        getCubeNormal(&tmp_p, &tmp_n);

        normalTransform(&tmp_n, n, cube);

        ray->rayPos(ray, *lambda, p);
    }
}

void sphereIntersect(struct object3D *sphere, struct ray3D *ray, double *lambda,
                     struct point3D *p, struct point3D *n, double *a,
                     double *b) {
    // Computes and returns the value of 'lambda' at the intersection
    // between the specified ray and the specified canonical sphere.

    /////////////////////////////////
    // TO DO: Complete this function. (DONE)
    /////////////////////////////////

    // The following calculations are based on the intersection defined in p.77
    //    A = d . d
    //    B = (a-c) . d
    //    C = (a-c) . (a-c) - 1
    //    D = B^2 - AC
    //    Lambda = (-B/A) +- (sqrt(D)/A)

    struct ray3D ray_transformed;
    rayTransform(ray, &ray_transformed, sphere);

    double A = dot(&ray_transformed.d, &ray_transformed.d);
    double B = dot(&ray_transformed.p0, &ray_transformed.d);
    double C = dot(&ray_transformed.p0, &ray_transformed.p0) - 1;
    double D = B * B - A * C;

    // No intersections
    if (D < 0)
        return;

    double lambda1 = -B / A + sqrt(D) / A;
    double lambda2 = -B / A - sqrt(D) / A;

    // One intersection
    if (D == 0) {
        *lambda = lambda1;

        // Two intersections
    } else {
        // Both intersections are behind view-plane and are not visible
        if (lambda1 < 0 && lambda2 < 0)
            return;

        // p(lambda)
        if (lambda1 > 0 && lambda2 < 0)
            *lambda = lambda1;
        else if (lambda1 > lambda2 && lambda2 > 0)
            *lambda = lambda2;
    }

    rayPosition(&ray_transformed, *lambda, p);

    normalTransform(p, n, sphere);
    normalize(n);

    *a = asin(n->px) / PI + 0.5;
    *b = asin(n->py) / PI + 0.5;

    ray->rayPos(ray, *lambda, p);
}

void cylIntersect(struct object3D *cylinder, struct ray3D *ray, double *lambda,
                  struct point3D *p, struct point3D *n, double *a, double *b) {
    // Computes and returns the value of 'lambda' at the intersection
    // between the specified ray and the specified canonical cylinder.

    struct ray3D ray_transformed;
    rayTransform(ray, &ray_transformed, cylinder);

    // The following side intersection calculations are based on the intersection defined in p.77
    //    x = (px + Lambda*dx)
    //    y = (py + Lambda*dy)
    //    A = dx^2 + dy^2
    //    B = 2*px*dx + 2*py*dy
    //    C = -1 + px^2 + py^2
    //    x^2 + y^2 = 1 -> a*Lambda^2 + b*Lambda + c  - 1 0
    //    D = b^2 - 4AC
    //    Lambda = (-B +- sqrt(D)/(2a)
    //    pz>=0 and pz <=1

    // The following calculcates top/bottom intersection
    //    z = pz + Lambda*dz
    //    z = 1 or z = 0
    //    x^2 + y^2 <= 1

    bool intersect_success = false;
    double min_lambda = __DBL_MAX__;
    double lambda1, lambda2, lambda3, lambda4;
    struct point3D test_p_1, test_p_2, test_p_3, test_p_4;
    struct point3D orig_n;
    double A, B, C, D;
    int side  = 0;
    A = pow(ray_transformed.d.px, 2) + pow(ray_transformed.d.py, 2);
    B = 2 * (ray_transformed.p0.px * ray_transformed.d.px + ray_transformed.p0.py * ray_transformed.d.py);
    C = -1 + pow(ray_transformed.p0.px, 2) + pow(ray_transformed.p0.py, 2);
    D = pow(B, 2) - 4 * A * C;

    if (D >= 0) {
        lambda1 = (-B + sqrt(D)) / (2 * A);
        lambda2 = (-B - sqrt(D)) / (2 * A);
    } else {
        lambda1 = -1;
        lambda2 = -1;
    }

    // side surface intersection

    // top/bottom surface lambda

    // top
    lambda3 = (1 - ray_transformed.p0.pz) / ray_transformed.d.pz;
    // bottom
    lambda4 = (-1 - ray_transformed.p0.pz) / ray_transformed.d.pz;

    // caclulate test intersection points
    rayPosition(&ray_transformed, lambda1, &test_p_1);
    rayPosition(&ray_transformed, lambda2, &test_p_2);
    rayPosition(&ray_transformed, lambda3, &test_p_3);
    rayPosition(&ray_transformed, lambda4, &test_p_4);

    // test side 0 <= z <= 1
    if (test_p_1.pz >= -1 && test_p_1.pz <= 1 && lambda1 > 0 && lambda1 < min_lambda) {
        min_lambda = lambda1;
        copyPoint(&test_p_1, &orig_n);
        orig_n.pz = 0;
        intersect_success = true;
        side = 1;
        *a = (atan(test_p_1.px / test_p_1.py) + PI/2 ) / PI;  
        *b = (test_p_1.pz + 1.0 ) / 2.0;
    }
    if (test_p_2.pz >= -1 && test_p_2.pz <= 1 && lambda2 > 0 && lambda2 < min_lambda) {
        min_lambda = lambda2;
        copyPoint(&test_p_2, &orig_n);
        orig_n.pz = 0;
        intersect_success = true;
        side = 1;
        *a = (atan(test_p_2.px / test_p_2.py) + PI/2 ) / PI;  
        *b = (test_p_2.pz + 1.0 ) / 2.0;
    }

    // test top x^2 + y^2 <=1
    if ((pow(test_p_3.px, 2) + pow(test_p_3.py, 2)) <= 1 && lambda3 > 0 && lambda3 < min_lambda) {
        min_lambda = lambda3;
        orig_n.px = 0;
        orig_n.py = 0;
        orig_n.pz = 1;
        orig_n.pw = 1;
        intersect_success = true;
        side = 0;
        *a = (test_p_3.px + 1) / 2;
        *b = (test_p_3.py + 1) / 2;
    }

    // test bottom x^2 + y^2 <=1
    if ((pow(test_p_4.px, 2) + pow(test_p_4.py, 2)) <= 1 && lambda4 > 0 && lambda4 < min_lambda) {
        min_lambda = lambda4;
        orig_n.px = 0;
        orig_n.py = 0;
        orig_n.pz = -1;
        orig_n.pw = 1;
        intersect_success = true;
        side = 0;
        *a = (test_p_4.px + 1) / 2;
        *b = (test_p_4.py + 1) / 2;
    }

    if (intersect_success) {
        // Within cylinder height
        *lambda = min_lambda;

        normalTransform(&orig_n, n, cylinder);
        normalize(n);

        ray->rayPos(ray, *lambda, p);
    }
}

void coneIntersect(struct object3D *cone, struct ray3D *ray, double *lambda,
                   struct point3D *p, struct point3D *n, double *a, double *b) {
    // Computes and returns the value of 'lambda' at the intersection
    // between the specified ray and the specified canonical cylinder.

    struct ray3D ray_transformed;
    rayTransform(ray, &ray_transformed, cone);

    // The following side intersection calculations are based on the intersection defined in p.77
    //    x = (px + Lambda*dx)
    //    y = (py + Lambda*dy)
    //    z = (pz + Lambda*dz)
    //    A = dx^2 + dy^2 - dz^2
    //    B = 2*px*dx + 2*py*dy + 2*dz - 2*pz*dz
    //    C = -1 + px^2 + py^2 + 2*pz - pz^2
    //    x^2 + y^2 = (1-z)^2 -> a*Lambda^2 + b*Lambda + c  - 1 = 0
    //    D = b^2 - 4AC
    //    Lambda = (-B +- sqrt(D)/(2a)

    // The following calculcates bottom intersection
    //    z = pz + Lambda*dz
    //    z = 0
    //    x^2 + y^2 <= 1

    bool intersect_success = false;
    double min_lambda = __DBL_MAX__;
    double lambda1, lambda2, lambda3;
    struct point3D test_p_1, test_p_2, test_p_3;
    struct point3D orig_n;
    double A, B, C, D;
    A = pow(ray_transformed.d.px, 2) + pow(ray_transformed.d.py, 2) - pow(ray_transformed.d.pz, 2);
    B = 2 * (ray_transformed.p0.px * ray_transformed.d.px + ray_transformed.p0.py * ray_transformed.d.py + ray_transformed.d.pz - ray_transformed.p0.pz * ray_transformed.d.pz);
    C = -1 + pow(ray_transformed.p0.px, 2) + pow(ray_transformed.p0.py, 2) + 2 * ray_transformed.p0.pz - pow(ray_transformed.p0.pz, 2);
    D = pow(B, 2) - 4 * A * C;

    // side surface intersection
    if (D >= 0) {
        lambda1 = (-B + sqrt(D)) / (2 * A);
        lambda2 = (-B - sqrt(D)) / (2 * A);
    } else {
        lambda1 = -1;
        lambda2 = -1;
    }

    // bottom
    lambda3 = (-ray_transformed.p0.pz) / ray_transformed.d.pz;

    // caclulate test intersection points
    rayPosition(&ray_transformed, lambda1, &test_p_1);
    rayPosition(&ray_transformed, lambda2, &test_p_2);
    rayPosition(&ray_transformed, lambda3, &test_p_3);

    // test side 0 <= z <= 1
    if (lambda1 > 0 && lambda1 < min_lambda && test_p_1.pz <= 1 && test_p_1.pz >= 0) {
        min_lambda = lambda1;
        copyPoint(&test_p_1, &orig_n);
        orig_n.pz = -orig_n.pz;
        intersect_success = true;
    }
    if (lambda2 > 0 && lambda2 < min_lambda && test_p_2.pz <= 1 && test_p_2.pz >= 0) {
        min_lambda = lambda2;
        copyPoint(&test_p_2, &orig_n);
        orig_n.pz = -orig_n.pz;
        intersect_success = true;
    }

    // test top x^2 + y^2 <=1
    if ((pow(test_p_3.px, 2) + pow(test_p_3.py, 2)) <= 1 && lambda3 > 0 && lambda3 < min_lambda) {
        min_lambda = lambda3;
        orig_n.px = 0;
        orig_n.py = 0;
        orig_n.pz = -1;
        orig_n.pw = 1;
        intersect_success = true;
    }

    if (intersect_success) {
        // Within cylinder height
        *lambda = min_lambda;

        normalTransform(&orig_n, n, cone);
        normalize(n);

        ray->rayPos(ray, *lambda, p);

        *a = atan(p->px / p->py);
        *b = p->pz;
    }
}
/////////////////////////////////////////////////////////////////
// Surface coordinates & random sampling on object surfaces
/////////////////////////////////////////////////////////////////
void planeCoordinates(struct object3D *plane, double a, double b, double *x,
                      double *y, double *z) {
    // Return in (x,y,z) the coordinates of a point on the plane given by the 2
    // parameters a,b in [0,1]. 'a' controls displacement from the left side of
    // the plane, 'b' controls displacement from the bottom of the plane.

    /////////////////////////////////
    // TO DO: Complete this function.
    /////////////////////////////////
}

void sphereCoordinates(struct object3D *sphere, double a, double b, double *x,
                       double *y, double *z) {
    // Return in (x,y,z) the coordinates of a point on the plane given by the 2
    // parameters a,b. 'a' in [0, 2*PI] corresponds to the spherical coordinate
    // theta 'b' in [-PI/2, PI/2] corresponds to the spherical coordinate phi

    /////////////////////////////////
    // TO DO: Complete this function.
    /////////////////////////////////
}

void cylCoordinates(struct object3D *cyl, double a, double b, double *x,
                    double *y, double *z) {
    // Return in (x,y,z) the coordinates of a point on the plane given by the 2
    // parameters a,b. 'a' in [0, 2*PI] corresponds to angle theta around the
    // cylinder 'b' in [0, 1] corresponds to height from the bottom

    /////////////////////////////////
    // TO DO: Complete this function.
    /////////////////////////////////
}

void coneCoordinates(struct object3D *cone, double a, double b, double *x,
                     double *y, double *z) {
    // Return in (x,y,z) the coordinates of a point on the plane given by the 2
    // parameters a,b. 'a' in [0, 2*PI] corresponds to angle theta around the
    // cylinder 'b' in [0, 1] corresponds to height from the bottom

    /////////////////////////////////
    // TO DO: Complete this function.
    /////////////////////////////////
}

void planeSample(struct object3D *plane, double *x, double *y, double *z) {
    // Returns the 3D coordinates (x,y,z) of a randomly sampled point on the plane
    // Sampling should be uniform, meaning there should be an equal change of
    // getting any spot on the plane

    double x_temp = get_random_double(-1.0, 1.0);
    double y_temp = get_random_double(-1.0, 1.0);
    double z_temp = 0;

    struct point3D transformed_p;
    transformed_p.px = x_temp;
    transformed_p.py = y_temp;
    transformed_p.pz = z_temp;
    transformed_p.pw = 1;

    matVecMult(plane->T, &transformed_p);

    *x = transformed_p.px;
    *y = transformed_p.py;
    *z = transformed_p.pz;

    /////////////////////////////////
    // TO DO: Complete this function.
    /////////////////////////////////
}

void sphereSample(struct object3D *sphere, double *x, double *y, double *z) {
    // Returns the 3D coordinates (x,y,z) of a randomly sampled point on the
    // sphere Sampling should be uniform - note that this is tricky for a sphere,
    // do some research and document in your report what method is used to do
    // this, along with a reference to your source.

    double theta = get_random_double(0, 2 * PI);
    double rho = get_random_double(0, 2 * PI);

    struct point3D tmp_p;
    tmp_p.px = sin(rho) * cos(theta);
    tmp_p.py = sin(rho) * sin(theta);
    tmp_p.pz = cos(rho);
    tmp_p.pw = 1;

    matVecMult(sphere->T, &tmp_p);

    *x = tmp_p.px;
    *y = tmp_p.py;
    *z = tmp_p.pz;
}

void cylSample(struct object3D *cyl, double *x, double *y, double *z) {
    // Returns the 3D coordinates (x,y,z) of a randomly sampled point on the
    // cylinder Sampling should be uniform over the cylinder.
    double theta = get_random_double(0, 2 * PI);
    double tmp_z = get_random_double(0, 1);

    struct point3D tmp_p;
    tmp_p.px = sin(theta);
    tmp_p.py = cos(theta);
    tmp_p.pz = tmp_z;
    tmp_p.pw = 1;

    matVecMult(cyl->T, &tmp_p);

    *x = tmp_p.px;
    *y = tmp_p.py;
    *z = tmp_p.pz;
}

int sphereOverlap(struct object3D *sphere, double x, double y, double z) {
    struct point3D tmp_p;
    assignPoint(&tmp_p, x, y, z);
    matVecMult(sphere->Tinv, &tmp_p);

    return pow(x, 2) + pow(y, 2) + pow(z, 2) <= 1;
}

int coneOverlap(struct object3D *cone, double x, double y, double z) {
    struct point3D tmp_p;
    assignPoint(&tmp_p, x, y, z);
    matVecMult(cone->Tinv, &tmp_p);

    return pow(x, 2) + pow(y, 2) <= z && z >= 0 && z <= 1;
}

int cylOverlap(struct object3D *cyl, double x, double y, double z) {
    struct point3D tmp_p;
    assignPoint(&tmp_p, x, y, z);
    matVecMult(cyl->Tinv, &tmp_p);

    return pow(x, 2) + pow(y, 2) <= 1 && z >= -1 && z <= 1;
}

void coneSample(struct object3D *cone, double *x, double *y, double *z) {
    // Returns the 3D coordinates (x,y,z) of a randomly sampled point on the
    // cylinder Sampling should be uniform over the cylinder.

    /////////////////////////////////
    // TO DO: Complete this function.
    /////////////////////////////////
}

/////////////////////////////////
// Texture mapping functions
/////////////////////////////////
void loadTexture(struct object3D *o, const char *filename, int type,
                 struct textureNode **t_list) {
    // Load a texture or normal map image from file and assign it to the
    // specified object.
    // type:   1  ->  Texture map  (RGB, .ppm)
    //         2  ->  Normal map   (RGB, .ppm)
    //         3  ->  Alpha map    (grayscale, .pgm)
    // Stores loaded images in a linked list to avoid replication
    struct image *im;
    struct textureNode *p;

    if (o != NULL) {
        // Check current linked list
        p = *(t_list);
        while (p != NULL) {
            if (strcmp(&p->name[0], filename) == 0) {
                // Found image already on the list
                if (type == 1)
                    o->texImg = p->im;
                else if (type == 2)
                    o->normalMap = p->im;
                else
                    o->alphaImg = p->im;
                return;
            }
            p = p->next;
        }

        // Load this texture image
        if (type == 1 || type == 2)
            im = readPPMimage(filename);
        else if (type == 3)
            im = readPGMimage(filename);

        // Insert it into the texture list
        if (im != NULL) {
            p = (struct textureNode *)calloc(1, sizeof(struct textureNode));
            strcpy(&p->name[0], filename);
            p->type = type;
            p->im = im;
            p->next = NULL;
            // Insert into linked list
            if ((*(t_list)) == NULL)
                *(t_list) = p;
            else {
                p->next = (*(t_list))->next;
                (*(t_list))->next = p;
            }
            // Assign to object
            if (type == 1)
                o->texImg = im;
            else if (type == 2) {
                o->normalMap = im;
                o->normalMapped = 1;
            } else {
                o->alphaImg = im;
                o->alphaMapped = 1;
            }
        }

    } // end if (o != NULL)
}

void texMap(struct image *img, double a, double b, double *R, double *G,
            double *B) {
    /*
     Function to determine the colour of a textured object at
     the normalized texture coordinates (a,b).

     a and b are texture coordinates in [0 1].
     img is a pointer to the image structure holding the texture for
      a given object.

     The colour is returned in R, G, B. Uses bi-linear interpolation
     to determine texture colour.
    */

    //////////////////////////////////////////////////
    // TO DO (Assignment 4 only) (DONE):
    //
    //  Complete this function to return the colour
    // of the texture image at the specified texture
    // coordinates. Your code should use bi-linear
    // interpolation to obtain the texture colour.
    //////////////////////////////////////////////////

    // Bilinear Interpolation reference:
    // https://cseweb.ucsd.edu/classes/wi18/cse167-a/lec9.pdf Slide 11 All
    // variables defined below are corresponding to the reference above

    a = max(0, min(a, 1));
    b = max(0, min(b, 1));

    // Texture image
    double *img_rgb = (double *)img->rgbdata;

    // Current pixel (s, t) in texture coordinates (however, it could be in
    // between actual image pixels)
    double s = a * img->sx - 1;
    double t = b * img->sy - 1;

    // Find surrounding pixels that make up a rectangle around (s, t)
    int s0 = (int)floor(s);
    int s1 = (int)ceil(s);
    int t0 = (int)floor(t);
    int t1 = (int)ceil(t);

    // Precalculation to avoid rerunning this code below. The essence is to end up
    // with the correct pixel in the texture image by adding tx_img to sx
    int t0_img = max(t0 * img->sx - 1, 0);
    int t1_img = max(t1 * img->sx - 1, 0);
    

    // Ratios in s and t directions
    double rs = ((double)(s - s0)) / ((double)(s1 - s0));
    double rt = ((double)(t - t0)) / ((double)(t1 - t0));

    // If bounding pixel is the same on either axis to prevent division by zero
    if (s1 - s0 < 1) rs = 1;
    if (t1 - t0 < 1) rt = 1;

    // printf("a=%f b=%f s=%f t=%f s0=%d s1=%d t0=%d t1=%d t0_img=%d t1_img%d rs=%f rt=%f\n", a, b, s, t, s0, s1, t0, t1, t0_img, t1_img, rs, rt);

    // Calculate top and bottom ratio for all colours
    double ctop_r = *(img_rgb + (t1_img + s0) * 3) * (1.0 - rs) +
                    *(img_rgb + (t1_img + s1) * 3) * rs;
    double cbot_r = *(img_rgb + (t0_img + s0) * 3) * (1.0 - rs) +
                    *(img_rgb + (t0_img + s1) * 3) * rs;

    double ctop_g = *(img_rgb + (t1_img + s0) * 3 + 1) * (1.0 - rs) +
                    *(img_rgb + (t1_img + s1) * 3 + 1) * rs;
    double cbot_g = *(img_rgb + (t0_img + s0) * 3 + 1) * (1.0 - rs) +
                    *(img_rgb + (t0_img + s1) * 3 + 1) * rs;

    double ctop_b = *(img_rgb + (t1_img + s0) * 3 + 2) * (1.0 - rs) +
                    *(img_rgb + (t1_img + s1) * 3 + 2) * rs;
    double cbot_b = *(img_rgb + (t0_img + s0) * 3 + 2) * (1.0 - rs) +
                    *(img_rgb + (t0_img + s1) * 3 + 2) * rs;

    // Get pixel colour based on ratios of surrounding pixels
    *R = cbot_r * (1.0 - rt) + ctop_r * rt;
    *G = cbot_g * (1.0 - rt) + ctop_g * rt;
    *B = cbot_b * (1.0 - rt) + ctop_b * rt;
}

void alphaMap(struct image *img, double a, double b, double *alpha) {
    // Just like texture map but returns the alpha value at a,b,
    // notice that alpha maps are single layer grayscale images, hence
    // the separate function.

    //////////////////////////////////////////////////
    // TO DO (Assignment 4 only):
    //
    //  Complete this function to return the alpha
    // value from the image at the specified texture
    // coordinates. Your code should use bi-linear
    // interpolation to obtain the texture colour.
    //////////////////////////////////////////////////

    // Bilinear Interpolation reference:
    // https://cseweb.ucsd.edu/classes/wi18/cse167-a/lec9.pdf Slide 11 All
    // variables defined below are corresponding to the reference above

    // Texture image
    double *img_rgb = (double *)img->rgbdata;

    // Current pixel (s, t) in texture coordinates (however, it could be in
    // between actual image pixels)
    double s = a * img->sx - 1;
    double t = b * img->sy - 1;

    // Find surrounding pixels that make up a rectangle around (s, t)
    int s0 = (int)floor(s);
    int s1 = (int)ceil(s);
    int t0 = (int)floor(t);
    int t1 = (int)ceil(t);

    // Precalculation to avoid rerunning this code below. The essence is to end up
    // with the correct pixel in the texture image by adding tx_img to sx
    int t0_img = t0 * img->sx - 1;
    int t1_img = t1 * img->sx - 1;

    // Ratios in s and t directions
    double rs = ((double)(s - s0)) / ((double)(s1 - s0));
    double rt = ((double)(t - t0)) / ((double)(t1 - t0));

    // Calculate top and bottom ratio for all colours
    double ctop = *(img_rgb + (t1_img + s0)) * (1.0 - rs) +
                  *(img_rgb + (t1_img + s1)) * rs;
    double cbot = *(img_rgb + (t0_img + s0)) * (1.0 - rs) +
                  *(img_rgb + (t0_img + s1)) * rs;

    // Get alpha value based on ratios of surrounding values
    *alpha = cbot * (1.0 - rt) + ctop * rt;
}

/////////////////////////////
// Light sources
/////////////////////////////
void insertPLS(struct pointLS *l, struct pointLS **list) {
    if (l == NULL)
        return;

    // Inserts a light source into the list of light sources
    if (*(list) == NULL) {
        *(list) = l;
        (*(list))->next = NULL;
    } else {
        l->next = (*(list))->next;
        (*(list))->next = l;
    }
}

void insertALS(struct areaLS *als, struct areaLS **list, struct object3D **obj_list) {
    if (als == NULL || als->light_shape == NULL)
        return;

    insertObject(als->light_shape, obj_list);

    // Inserts a light source into the list of light sources
    if (*(list) == NULL) {
        *(list) = als;
        (*(list))->next = NULL;
    } else {
        als->next = (*(list))->next;
        (*(list))->next = als;
    }
}

void addAreaLight(double sx, double sy, double nx, double ny, double nz,
                  double tx, double ty, double tz, int N, double r, double g,
                  double b, struct object3D **o_list, struct pointLS **l_list) {
    /*
      This function sets up and inserts a rectangular area light source
      with size (sx, sy)
      orientation given by the normal vector (nx, ny, nz)
      centered at (tx, ty, tz)
      consisting of (N) point light sources (uniformly sampled)
      and with colour/intensity (r,g,b)

      Note that the light source must be visible as a uniformly colored rectangle
      which casts no shadows. If you require a lightsource to shade another, you
      must make it into a proper solid box with a back and sides of
      non-light-emitting material
    */

    /////////////////////////////////////////////////////
    // TO DO: (Assignment 4!)
    // Implement this function to enable area light sources
    /////////////////////////////////////////////////////

    // NOTE: The best way to implement area light sources is to random sample from
    // the
    //       light source's object surface within rtShade(). This is a bit more
    //       tricky but reduces artifacts significantly. If you do that, then
    //       there is no need to insert a series of point lightsources in this
    //       function.
}

///////////////////////////////////
// Geometric transformation section
///////////////////////////////////

void invert(double *T, double *Tinv) {
    // Computes the inverse of transformation matrix T.
    // the result is returned in Tinv.

    double *U, *s, *V, *rv1;
    int singFlag, i;

    // Invert the affine transform
    U = NULL;
    s = NULL;
    V = NULL;
    rv1 = NULL;
    singFlag = 0;

    SVD(T, 4, 4, &U, &s, &V, &rv1);
    if (U == NULL || s == NULL || V == NULL) {
        fprintf(
            stderr,
            "Error: Matrix not invertible for this object, returning identity\n");
        memcpy(Tinv, eye4x4, 16 * sizeof(double));
        return;
    }

    // Check for singular matrices...
    for (i = 0; i < 4; i++)
        if (*(s + i) < 1e-9)
            singFlag = 1;
    if (singFlag) {
        fprintf(stderr,
                "Error: Transformation matrix is singular, returning identity\n");
        memcpy(Tinv, eye4x4, 16 * sizeof(double));
        return;
    }

    // Compute and store inverse matrix
    InvertMatrix(U, s, V, 4, Tinv);

    free(U);
    free(s);
    free(V);
}

void RotateXMat(double T[4][4], double theta) {
    // Multiply the current object transformation matrix T in object o
    // by a matrix that rotates the object theta *RADIANS* around the
    // X axis.

    double R[4][4];
    memset(&R[0][0], 0, 16 * sizeof(double));

    R[0][0] = 1.0;
    R[1][1] = cos(theta);
    R[1][2] = -sin(theta);
    R[2][1] = sin(theta);
    R[2][2] = cos(theta);
    R[3][3] = 1.0;

    matMult(R, T);
}

void RotateX(struct object3D *o, double theta) {
    // Multiply the current object transformation matrix T in object o
    // by a matrix that rotates the object theta *RADIANS* around the
    // X axis.

    double R[4][4];
    memset(&R[0][0], 0, 16 * sizeof(double));

    R[0][0] = 1.0;
    R[1][1] = cos(theta);
    R[1][2] = -sin(theta);
    R[2][1] = sin(theta);
    R[2][2] = cos(theta);
    R[3][3] = 1.0;

    matMult(R, o->T);
}

void RotateYMat(double T[4][4], double theta) {
    // Multiply the current object transformation matrix T in object o
    // by a matrix that rotates the object theta *RADIANS* around the
    // Y axis.

    double R[4][4];
    memset(&R[0][0], 0, 16 * sizeof(double));

    R[0][0] = cos(theta);
    R[0][2] = sin(theta);
    R[1][1] = 1.0;
    R[2][0] = -sin(theta);
    R[2][2] = cos(theta);
    R[3][3] = 1.0;

    matMult(R, T);
}

void RotateY(struct object3D *o, double theta) {
    // Multiply the current object transformation matrix T in object o
    // by a matrix that rotates the object theta *RADIANS* around the
    // Y axis.

    double R[4][4];
    memset(&R[0][0], 0, 16 * sizeof(double));

    R[0][0] = cos(theta);
    R[0][2] = sin(theta);
    R[1][1] = 1.0;
    R[2][0] = -sin(theta);
    R[2][2] = cos(theta);
    R[3][3] = 1.0;

    matMult(R, o->T);
}

void RotateZMat(double T[4][4], double theta) {
    // Multiply the current object transformation matrix T in object o
    // by a matrix that rotates the object theta *RADIANS* around the
    // Z axis.

    double R[4][4];
    memset(&R[0][0], 0, 16 * sizeof(double));

    R[0][0] = cos(theta);
    R[0][1] = -sin(theta);
    R[1][0] = sin(theta);
    R[1][1] = cos(theta);
    R[2][2] = 1.0;
    R[3][3] = 1.0;

    matMult(R, T);
}

void RotateZ(struct object3D *o, double theta) {
    // Multiply the current object transformation matrix T in object o
    // by a matrix that rotates the object theta *RADIANS* around the
    // Z axis.

    double R[4][4];
    memset(&R[0][0], 0, 16 * sizeof(double));

    R[0][0] = cos(theta);
    R[0][1] = -sin(theta);
    R[1][0] = sin(theta);
    R[1][1] = cos(theta);
    R[2][2] = 1.0;
    R[3][3] = 1.0;

    matMult(R, o->T);
}

void RotateObjToVec(struct object3D *o, struct point3D *s, struct point3D *d) {
    
    // apply rotation matrix to o such that vector s aligns with vector d
    struct point3D temp_s, temp_d;
    copyPoint(s, &temp_s);
    copyPoint(d, &temp_d);
    temp_s.pw = 0;
    temp_d.pw = 0;
    normalize(&temp_s);
    normalize(&temp_d);

    struct point3D *v = cross(&temp_s, &temp_d);
    double c = dot(&temp_s, &temp_d);

    double R[4][4];
    memset(&R[0][0], 0, 16 * sizeof(double));
    double Vx[4][4];
    memset(&Vx[0][0], 0, 16 * sizeof(double));
    double Vx2[4][4];
    memset(&Vx[0][0], 0, 16 * sizeof(double));

    // Calculation credit:
    // R = I + Vx + Vx^2*(1/(1+c))
    // https://math.stackexchange.com/questions/180418/calculate-rotation-matrix-to-align-vector-a-to-vector-b-in-3d

    // 3x3 on top left corner
    Vx[0][1] = -v->pz;
    Vx[0][2] = v->py;
    Vx[1][0] = v->pz;
    Vx[1][2] = -v->px;
    Vx[2][0] = -v->py;
    Vx[2][1] = v->px;

    // identity R
    R[0][0] = 1;
    R[1][1] = 1;
    R[2][2] = 1;
    R[3][3] = 1.0;

    // Vx2 = Vx
    memcpy(Vx2, Vx, 16 * sizeof(double));

    // Vx2 = Vx2^2
    matMult(Vx2, Vx2);

    // Vx2 = 1/(1+c) * Vx2
    matScaMult(Vx2, 1 / (1 + c));

    // R = R + Vx
    matAdd(Vx, R);

    // R = R + Vx2
    matAdd(Vx2, R);
    matMult(R, o->T);

    free(v);
}



void RotateVecToVec(struct point3D *s, struct point3D *d){
    // apply rotation matrix to o such that vector s aligns with vector d
    struct point3D temp_s, temp_d;
    copyPoint(s, &temp_s);
    copyPoint(d, &temp_d);
    temp_s.pw = 0;
    temp_d.pw = 0;
    normalize(&temp_s);
    normalize(&temp_d);

    struct point3D *v = cross(&temp_s, &temp_d);
    double c = dot(&temp_s, &temp_d);

    double R[4][4];
    memset(&R[0][0], 0, 16 * sizeof(double));
    double Vx[4][4];
    memset(&Vx[0][0], 0, 16 * sizeof(double));
    double Vx2[4][4];
    memset(&Vx[0][0], 0, 16 * sizeof(double));

    // Calculation credit:
    // R = I + Vx + Vx^2*(1/(1+c))
    // https://math.stackexchange.com/questions/180418/calculate-rotation-matrix-to-align-vector-a-to-vector-b-in-3d

    // 3x3 on top left corner
    Vx[0][1] = -v->pz;
    Vx[0][2] = v->py;
    Vx[1][0] = v->pz;
    Vx[1][2] = -v->px;
    Vx[2][0] = -v->py;
    Vx[2][1] = v->px;

    // identity R
    R[0][0] = 1;
    R[1][1] = 1;
    R[2][2] = 1;
    R[3][3] = 1.0;

    // Vx2 = Vx
    memcpy(Vx2, Vx, 16 * sizeof(double));

    // Vx2 = Vx2^2
    matMult(Vx2, Vx2);

    // Vx2 = 1/(1+c) * Vx2
    matScaMult(Vx2, 1 / (1 + c));

    // R = R + Vx
    matAdd(Vx, R);

    // R = R + Vx2
    matAdd(Vx2, R);
    free(v);

    matVecMult(R,d);
}




void getRotationalMatrix(struct point3D *s, struct point3D *d, double *T){
    
    // apply rotation matrix to o such that vector s aligns with vector d
    struct point3D temp_s, temp_d;
    copyPoint(s, &temp_s);
    copyPoint(d, &temp_d);
    temp_s.pw = 0;
    temp_d.pw = 0;
    normalize(&temp_s);
    normalize(&temp_d);

    struct point3D *v = cross(&temp_s, &temp_d);
    double c = dot(&temp_s, &temp_d);

    double R[4][4];
    memset(&R[0][0], 0, 16 * sizeof(double));
    double Vx[4][4];
    memset(&Vx[0][0], 0, 16 * sizeof(double));
    double Vx2[4][4];
    memset(&Vx[0][0], 0, 16 * sizeof(double));

    // Calculation credit:
    // R = I + Vx + Vx^2*(1/(1+c))
    // https://math.stackexchange.com/questions/180418/calculate-rotation-matrix-to-align-vector-a-to-vector-b-in-3d

    // 3x3 on top left corner
    Vx[0][1] = -v->pz;
    Vx[0][2] = v->py;
    Vx[1][0] = v->pz;
    Vx[1][2] = -v->px;
    Vx[2][0] = -v->py;
    Vx[2][1] = v->px;

    // identity R
    R[0][0] = 1;
    R[1][1] = 1;
    R[2][2] = 1;
    R[3][3] = 1.0;

    // Vx2 = Vx
    memcpy(Vx2, Vx, 16 * sizeof(double));

    // Vx2 = Vx2^2
    matMult(Vx2, Vx2);

    // Vx2 = 1/(1+c) * Vx2
    matScaMult(Vx2, 1 / (1 + c));

    // R = R + Vx
    matAdd(Vx, R);

    // R = R + Vx2
    matAdd(Vx2, R);
    memcpy(T,R,16 * sizeof(double));

    free(v);
}

void TranslateMat(double T[4][4], double tx, double ty, double tz) {
    // Multiply the current object transformation matrix T in object o
    // by a matrix that translates the object by the specified amounts.

    double tr[4][4];
    memset(&tr[0][0], 0, 16 * sizeof(double));

    tr[0][0] = 1.0;
    tr[1][1] = 1.0;
    tr[2][2] = 1.0;
    tr[0][3] = tx;
    tr[1][3] = ty;
    tr[2][3] = tz;
    tr[3][3] = 1.0;

    matMult(tr, T);
}

void Translate(struct object3D *o, double tx, double ty, double tz) {
    // Multiply the current object transformation matrix T in object o
    // by a matrix that translates the object by the specified amounts.

    double tr[4][4];
    memset(&tr[0][0], 0, 16 * sizeof(double));

    tr[0][0] = 1.0;
    tr[1][1] = 1.0;
    tr[2][2] = 1.0;
    tr[0][3] = tx;
    tr[1][3] = ty;
    tr[2][3] = tz;
    tr[3][3] = 1.0;

    matMult(tr, o->T);
}

void ScaleMat(double T[4][4], double sx, double sy, double sz) {
    // Multiply the current object transformation matrix T in object o
    // by a matrix that scales the object as indicated.

    double S[4][4];
    memset(&S[0][0], 0, 16 * sizeof(double));

    S[0][0] = sx;
    S[1][1] = sy;
    S[2][2] = sz;
    S[3][3] = 1.0;

    matMult(S, T);
}

void Scale(struct object3D *o, double sx, double sy, double sz) {
    // Multiply the current object transformation matrix T in object o
    // by a matrix that scales the object as indicated.

    double S[4][4];
    memset(&S[0][0], 0, 16 * sizeof(double));

    S[0][0] = sx;
    S[1][1] = sy;
    S[2][2] = sz;
    S[3][3] = 1.0;

    matMult(S, o->T);
}

void printmatrix(double mat[4][4]) {
    fprintf(stderr, "Matrix contains:\n");
    fprintf(stderr, "%f %f %f %f\n", mat[0][0], mat[0][1], mat[0][2], mat[0][3]);
    fprintf(stderr, "%f %f %f %f\n", mat[1][0], mat[1][1], mat[1][2], mat[1][3]);
    fprintf(stderr, "%f %f %f %f\n", mat[2][0], mat[2][1], mat[2][2], mat[2][3]);
    fprintf(stderr, "%f %f %f %f\n", mat[3][0], mat[3][1], mat[3][2], mat[3][3]);
}

/////////////////////////////////////////
// Camera and view setup
/////////////////////////////////////////
struct view *setupView(struct point3D *e, struct point3D *g, struct point3D *up,
                       double f, double wl, double wt, double wsize) {
    /*
      This function sets up the camera axes and viewing direction as discussed in
      the lecture notes. e - Camera center g - Gaze direction up - Up vector fov -
      Fild of view in degrees f - focal length
    */
    struct view *c;
    struct point3D *u, *v;

    u = v = NULL;

    // Allocate space for the camera structure
    c = (struct view *)calloc(1, sizeof(struct view));
    if (c == NULL) {
        fprintf(stderr, "Out of memory setting up camera model!\n");
        return (NULL);
    }

    // Set up camera center and axes
    c->e.px = e->px; // Copy camera center location, note we must make sure
    c->e.py = e->py; // the camera center provided to this function has pw=1
    c->e.pz = e->pz;
    c->e.pw = 1;

    // Set up w vector (camera's Z axis). w=-g/||g||
    c->w.px = -g->px;
    c->w.py = -g->py;
    c->w.pz = -g->pz;
    c->w.pw = 1;
    normalize(&c->w);

    // Set up the horizontal direction, which must be perpenticular to w and up
    u = cross(&c->w, up);
    normalize(u);
    c->u.px = u->px;
    c->u.py = u->py;
    c->u.pz = u->pz;
    c->u.pw = 1;

    // Set up the remaining direction, v=(u x w)  - Mind the signs
    v = cross(&c->u, &c->w);
    normalize(v);
    c->v.px = v->px;
    c->v.py = v->py;
    c->v.pz = v->pz;
    c->v.pw = 1;

    // Copy focal length and window size parameters
    c->f = f;
    c->wl = wl;
    c->wt = wt;
    c->wsize = wsize;

    // Set up coordinate conversion matrices
    // Camera2World matrix (M_cw in the notes)
    // Mind the indexing convention [row][col]
    c->C2W[0][0] = c->u.px;
    c->C2W[1][0] = c->u.py;
    c->C2W[2][0] = c->u.pz;
    c->C2W[3][0] = 0;

    c->C2W[0][1] = c->v.px;
    c->C2W[1][1] = c->v.py;
    c->C2W[2][1] = c->v.pz;
    c->C2W[3][1] = 0;

    c->C2W[0][2] = c->w.px;
    c->C2W[1][2] = c->w.py;
    c->C2W[2][2] = c->w.pz;
    c->C2W[3][2] = 0;

    c->C2W[0][3] = c->e.px;
    c->C2W[1][3] = c->e.py;
    c->C2W[2][3] = c->e.pz;
    c->C2W[3][3] = 1;

    // World2Camera matrix (M_wc in the notes)
    // Mind the indexing convention [row][col]
    c->W2C[0][0] = c->u.px;
    c->W2C[1][0] = c->v.px;
    c->W2C[2][0] = c->w.px;
    c->W2C[3][0] = 0;

    c->W2C[0][1] = c->u.py;
    c->W2C[1][1] = c->v.py;
    c->W2C[2][1] = c->w.py;
    c->W2C[3][1] = 0;

    c->W2C[0][2] = c->u.pz;
    c->W2C[1][2] = c->v.pz;
    c->W2C[2][2] = c->w.pz;
    c->W2C[3][2] = 0;

    c->W2C[0][3] = -dot(&c->u, &c->e);
    c->W2C[1][3] = -dot(&c->v, &c->e);
    c->W2C[2][3] = -dot(&c->w, &c->e);
    c->W2C[3][3] = 1;

    free(u);
    free(v);
    return (c);
}

/////////////////////////////////////////
// Image I/O section
/////////////////////////////////////////
struct image *readPPMimage(const char *filename) {
    // Reads an image from a .ppm file. A .ppm file is a very simple image
    // representation format with a text header followed by the binary RGB data at
    // 24bits per pixel. The header has the following form:
    //
    // P6
    // # One or more comment lines preceded by '#'
    // 340 200
    // 255
    //
    // The first line 'P6' is the .ppm format identifier, this is followed by one
    // or more lines with comments, typically used to inidicate which program
    // generated the .ppm file. After the comments, a line with two integer values
    // specifies the image resolution as number of pixels in x and number of
    // pixels in y. The final line of the header stores the maximum value for
    // pixels in the image, usually 255. After this last header line, binary data
    // stores the RGB values for each pixel in row-major order. Each pixel
    // requires 3 bytes ordered R, G, and B.
    //
    // NOTE: Windows file handling is rather crotchetty. You may have to change
    // the
    //       way this file is accessed if the images are being corrupted on read
    //       on Windows.
    //
    // readPPMdata converts the image colour information to floating point. This
    // is so that the texture mapping function doesn't have to do the conversion
    // every time it is asked to return the colour at a specific location.
    //

    FILE *f;
    struct image *im;
    char line[1024];
    int sizx, sizy;
    int i;
    unsigned char *tmp;
    double *fRGB;
    int tmpi;
    char *tmpc;

    im = (struct image *)calloc(1, sizeof(struct image));
    if (im != NULL) {
        im->rgbdata = NULL;
        f = fopen(filename, "rb+");
        if (f == NULL) {
            fprintf(
                stderr,
                "Unable to open file %s for reading, please check name and path\n",
                filename);
            free(im);
            return (NULL);
        }
        tmpc = fgets(&line[0], 1000, f);
        if (strcmp(&line[0], "P6\n") != 0) {
            fprintf(stderr,
                    "Wrong file format, not a .ppm file or header end-of-line "
                    "characters missing\n");
            free(im);
            fclose(f);
            return (NULL);
        }
        fprintf(stderr, "%s\n", line);
        // Skip over comments
        tmpc = fgets(&line[0], 511, f);
        while (line[0] == '#') {
            fprintf(stderr, "%s", line);
            tmpc = fgets(&line[0], 511, f);
        }
        sscanf(&line[0], "%d %d\n", &sizx, &sizy); // Read file size
        fprintf(stderr, "nx=%d, ny=%d\n\n", sizx, sizy);
        im->sx = sizx;
        im->sy = sizy;

        tmpc = fgets(&line[0], 9, f); // Read the remaining header line
        fprintf(stderr, "%s\n", line);
        tmp = (unsigned char *)calloc(sizx * sizy * 3, sizeof(unsigned char));
        fRGB = (double *)calloc(sizx * sizy * 3, sizeof(double));
        if (tmp == NULL || fRGB == NULL) {
            fprintf(stderr, "Out of memory allocating space for image\n");
            free(im);
            fclose(f);
            return (NULL);
        }

        tmpi = fread(tmp, sizx * sizy * 3 * sizeof(unsigned char), 1, f);
        fclose(f);

        // Conversion to floating point
        for (i = 0; i < sizx * sizy * 3; i++)
            *(fRGB + i) = ((double)*(tmp + i)) / 255.0;
        free(tmp);
        im->rgbdata = (void *)fRGB;

        return (im);
    }

    fprintf(stderr, "Unable to allocate memory for image structure\n");
    return (NULL);
}

struct image *readPGMimage(const char *filename) {
    // Just like readPPMimage() except it is used to load grayscale alpha maps. In
    // alpha maps, a value of 255 corresponds to alpha=1 (fully opaque) and 0
    // correspondst to alpha=0 (fully transparent).
    // A .pgm header of the following form is expected:
    //
    // P5
    // # One or more comment lines preceded by '#'
    // 340 200
    // 255
    //
    // readPGMdata converts the image grayscale data to double floating point in
    // [0,1].

    FILE *f;
    struct image *im;
    char line[1024];
    int sizx, sizy;
    int i;
    unsigned char *tmp;
    double *fRGB;
    int tmpi;
    char *tmpc;

    im = (struct image *)calloc(1, sizeof(struct image));
    if (im != NULL) {
        im->rgbdata = NULL;
        f = fopen(filename, "rb+");
        if (f == NULL) {
            fprintf(
                stderr,
                "Unable to open file %s for reading, please check name and path\n",
                filename);
            free(im);
            return (NULL);
        }
        tmpc = fgets(&line[0], 1000, f);
        if (strcmp(&line[0], "P5\n") != 0) {
            fprintf(stderr,
                    "Wrong file format, not a .pgm file or header end-of-line "
                    "characters missing\n");
            free(im);
            fclose(f);
            return (NULL);
        }
        // Skip over comments
        tmpc = fgets(&line[0], 511, f);
        while (line[0] == '#')
            tmpc = fgets(&line[0], 511, f);
        sscanf(&line[0], "%d %d\n", &sizx, &sizy); // Read file size
        im->sx = sizx;
        im->sy = sizy;

        tmpc = fgets(&line[0], 9, f); // Read the remaining header line
        tmp = (unsigned char *)calloc(sizx * sizy, sizeof(unsigned char));
        fRGB = (double *)calloc(sizx * sizy, sizeof(double));
        if (tmp == NULL || fRGB == NULL) {
            fprintf(stderr, "Out of memory allocating space for image\n");
            free(im);
            fclose(f);
            return (NULL);
        }

        tmpi = fread(tmp, sizx * sizy * sizeof(unsigned char), 1, f);
        fclose(f);

        // Conversion to double floating point
        for (i = 0; i < sizx * sizy; i++)
            *(fRGB + i) = ((double)*(tmp + i)) / 255.0;
        free(tmp);
        im->rgbdata = (void *)fRGB;

        return (im);
    }

    fprintf(stderr, "Unable to allocate memory for image structure\n");
    return (NULL);
}

struct image *newImage(int size_x, int size_y) {
    // Allocates and returns a new image with all zeros. Assumes 24 bit per pixel,
    // unsigned char array.
    struct image *im;

    im = (struct image *)calloc(1, sizeof(struct image));
    if (im != NULL) {
        im->rgbdata = NULL;
        im->sx = size_x;
        im->sy = size_y;
        im->rgbdata = (void *)calloc(size_x * size_y * 3, sizeof(unsigned char));
        if (im->rgbdata != NULL)
            return (im);
    }
    fprintf(stderr, "Unable to allocate memory for new image\n");
    return (NULL);
}

void imageOutput(struct image *im, const char *filename) {
    // Writes out a .ppm file from the image data contained in 'im'.
    // Note that Windows typically doesn't know how to open .ppm
    // images. Use Gimp or any other seious image processing
    // software to display .ppm images.
    // Also, note that because of Windows file format management,
    // you may have to modify this file to get image output on
    // Windows machines to work properly.
    //
    // Assumes a 24 bit per pixel image stored as unsigned chars
    //

    FILE *f;

    if (im != NULL)
        if (im->rgbdata != NULL) {
            f = fopen(filename, "wb+");
            if (f == NULL) {
                fprintf(stderr, "Unable to open file %s for output! No image written\n",
                        filename);
                return;
            }
            fprintf(f, "P6\n");
            fprintf(f, "# Output from RayTracer.c\n");
            fprintf(f, "%d %d\n", im->sx, im->sy);
            fprintf(f, "255\n");
            fwrite((unsigned char *)im->rgbdata,
                   im->sx * im->sy * 3 * sizeof(unsigned char), 1, f);
            fclose(f);
            return;
        }
    fprintf(stderr, "imageOutput(): Specified image is empty. Nothing output\n");
}

void deleteImage(struct image *im) {
    // De-allocates memory reserved for the image stored in 'im'
    if (im != NULL) {
        if (im->rgbdata != NULL)
            free(im->rgbdata);
        free(im);
    }
}

void cleanup(struct object3D *o_list, struct pointLS *l_list,
             struct textureNode *t_list, struct areaLS *als_list) {
    // De-allocates memory reserved for the object list and the point light source
    // list. Note that *YOU* must de-allocate any memory reserved for images
    // rendered by the raytracer.
    struct object3D *p, *q;
    struct pointLS *r, *s;
    struct textureNode *t, *u;
    struct areaLS *als_r, *als_s;

    p = o_list; // De-allocate all memory from objects in the list
    while (p != NULL) {
        q = p->next;
        if (p->photonMap !=
            NULL) // If object is photon mapped, free photon map memory
        {
            if (p->photonMap->rgbdata != NULL)
                free(p->photonMap->rgbdata);
            free(p->photonMap);
        }
        free(p);
        p = q;
    }

    r = l_list; // Delete light source list
    while (r != NULL) {
        s = r->next;
        free(r);
        r = s;
    }

    als_r = als_list; // Delete area light source list
    while (r != NULL) {
        als_s = als_r->next;
        free(als_r->light_shape);
        free(als_r);
        als_r = als_s;
    }

    t = t_list; // Delete texture Images
    while (t != NULL) {
        u = t->next;
        if (t->im->rgbdata != NULL)
            free(t->im->rgbdata);
        free(t->im);
        free(t);
        t = u;
    }
}

void fractal_sphere(struct object3D *prev_o, int depth, double r, struct object3D **list, struct textureNode **textureList) {
    // a hierarchical sphere generation script that generate child shperes around the parent spheres

    if (depth <= 0) {
        // reaches defined recursive depth
        return;
    }

    double step = (PI / 3);
    double zstep = PI / 4;
    struct object3D *cur_object;
    double theta = 0;
    double phi = -PI / 4;
    int side_num = 6;
    for (size_t i = 0; i < 4; i++) {
        phi += zstep;
        theta = step;
        side_num = i == 0 ? 1 : 6;
        for (size_t o = 0; o < side_num; o++) {
            theta += step;

            // uses the parametric equation of a shpere
            // S(theta,phi) = [r*sin(phi)*cos(theta) r*sin(phi)*sin(theta) r*cos(phi)]
            // to generate points around the parent sphere (mother_r + child_r = 1 + r)
            struct point3D p;
            p.px = (1 + r) * sin(phi) * cos(theta);
            p.py = (1 + r) * sin(phi) * sin(theta);
            p.pz = (1 + r) * cos(phi);
            p.pw = 0;
            cur_object = duplicateObj(prev_o, 0);
            struct point3D up;
            up.px = 0;
            up.py = 0;
            up.pz = 1;
            up.pw = 0;
            RotateObjToVec(cur_object, &up, &p);
            Scale(cur_object, r, r, r);
            Translate(cur_object, p.px, p.py, p.pz);
            matMult(prev_o->T, cur_object->T);
            invert(&cur_object->T[0][0], &cur_object->Tinv[0][0]);
            // loadTexture(cur_object,"./Texture/n_wood.ppm",2,textureList);
            // loadTexture(cur_object,"./Texture/t_wood.ppm",1,textureList);
            insertObject(cur_object, list);
            // recursive call
            fractal_sphere(cur_object, depth - 1, r, list, textureList);
        }
    }
}