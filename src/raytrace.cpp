// Source file for raytracing code


// Include files

#include "R2/R2.h"
#include "R3/R3.h"
#include "R3Scene.h"
#include "raytrace.h"

#include <limits>
#include <iostream>

using namespace std;

#define EPSILON 0.001

////////////////////////////////////////////////////////////////////////
// Create image from scene
//
// This is the main ray tracing function called from raypro
// 
// "width" and "height" indicate the size of the ray traced image
//   (keep these small during debugging to speed up your code development cycle)
//
// "max_depth" indicates the maximum number of secondary reflections/transmissions to trace for any ray
//   (i.e., stop tracing a ray if it has already been reflected max_depth times -- 
//   0 means direct illumination, 1 means one bounce, etc.)
//
// "num_primary_rays_per_pixel" indicates the number of random rays to generate within 
//   each pixel during antialiasing.  This argument can be ignored if antialiasing is not implemented.
//
// "num_distributed_rays_per_intersection" indicates the number of secondary rays to generate
//   for each surface intersection if distributed ray tracing is implemented.  
//   It can be ignored otherwise.
// 
////////////////////////////////////////////////////////////////////////

R2Image *RenderImage(R3Scene *scene, int width, int height, int max_depth,
  int num_primary_rays_per_pixel, int num_distributed_rays_per_intersection)
{
  // Allocate  image
  R2Image *image = new R2Image(width, height);
  if (!image) {
    fprintf(stderr, "Unable to allocate image\n");
    return NULL;
  }

  // Calculate each pixel value
  for (int i = 0; i < width; i++)
  {
    cout << i << endl;
    for (int j = 0; j < height; j++)
    {
      // Cast primary ray through the pixel
      R3Ray ray = RayThoughPixel(scene->Camera(), i, j, width, height);
      // Compute the radiance traveling into the camera along that ray
      R3Rgb rgb = ComputeRadiance(scene, ray, 0, max_depth, num_distributed_rays_per_intersection);
      // set image pixel value
      image->SetPixel(i, j, rgb);
    }
  }

  // Return image
  return image;
}

// Construct a ray through a given pixel of a camera
R3Ray RayThoughPixel(R3Camera camera, int i, int j, int width, int height)
{
  R3Vector vx = tan(camera.xfov) * camera.right;
  R3Vector vy = tan(camera.yfov) * camera.up;

  R3Point p = camera.eye + camera.towards 
    + ((i + 0.5) / width  - 0.5) * 2 * vx
    + ((j + 0.5) / height - 0.5) * 2 * vy;

  return R3Ray(camera.eye, p);
}

// Compute the radiance traveling back along a given ray
R3Rgb ComputeRadiance(R3Scene *scene, R3Ray &ray, int depth, int max_depth, int distrib)
{
  // Compute the nearest intersection of the ray with a surface in the scene
  R3Intersection intersection = ComputeIntersection(scene, scene->Root(), ray, numeric_limits<double>::infinity());

  // Compute the radiance emitted back along the ray from the intersection point
  R3Rgb il(0, 0, 0, 1);
  if (intersection.hit) 
  {
    il = ComputeRadiance(scene, ray, intersection, depth, max_depth, distrib);
  }
  else
  {
    // if no intersection found, return scene background color
    il = scene->background;
  }

  return il;
}

// compute the ray reflected from the incoming ray r at intersection i
R3Ray GetReflectedRay(R3Intersection i, R3Ray r)
{
  // add offset to point to ensure the same intersection is not found again due to floating point error
  R3Point offset = i.position + i.normal * EPSILON;
  R3Vector v = -r.Vector();
  return R3Ray(offset, 2 * v.Dot(i.normal) * i.normal - v);
}

// compute ray refracted from incoming ray r using Snell's Law
// NOTE: I assume that refracting objects will always be non-intersecting and that all intersections passed into this function
// will be from empty space with an IOR of 1 to the interior of a shape
R3Ray GetRefractedRay(R3Intersection i, R3Ray r, double ior, R3Scene *scene)
{
  // this time push the point inside the object
  R3Point internal_offset = i.position - i.normal * EPSILON;

  // calculate the ray that will be cast through the object
  R3Vector v = -r.Vector();
  double cos_vn = v.Dot(i.normal);
  R3Vector internal_vec = i.normal * (cos_vn / ior - sqrt(1 - (1 - cos_vn * cos_vn) / (ior * ior))) - v / ior;
  R3Ray trans_ray(internal_offset, internal_vec);

  // cast this ray and determine it's exit point from the refracting object
  R3Intersection exit_i = ComputeIntersection(scene, scene->Root(), trans_ray, numeric_limits<double>::infinity());
  assert (exit_i.hit);

  // calculate the refracted exit ray and 
  R3Point external_offset = exit_i.position - exit_i.normal * EPSILON;
  v = -trans_ray.Vector();
  cos_vn = v.Dot(exit_i.normal);
  R3Vector external_vec = exit_i.normal * (cos_vn * ior - sqrt(1 - (1 - cos_vn * cos_vn) * (ior * ior))) - v * ior;
  return R3Ray(external_offset, external_vec);
}

// Compute the radiance emitted back along a given ray from a given intersection
R3Rgb ComputeRadiance(R3Scene *scene, R3Ray &ray, R3Intersection &intersection, int depth, int max_depth, int distrib)
{
  R3Rgb il(0, 0, 0, 1);
  R3Material *mat = intersection.node->material;
  // add a small epsilon normal to the intersection point so light and secondary rays don't intersect with the same surface
  R3Point offset_intersection_pos = intersection.position + intersection.normal * EPSILON;
  // v = a unit vector pointing from the intersection to the ray origin
  R3Vector v = -ray.Vector();

  // if we are not using distributed ray tracing...
  if (distrib == 0)
  {
    // include lighting from all scene lights
    for (int i = 0; i < scene->NLights(); i++)
    {
      R3Light *light = scene->Light(i);

      // l = a unit vector pointing from the intersection to a light
      // d = the distance between the intersection and a light
      R3Vector l;
      double d;
      if (light->type == R3_DIRECTIONAL_LIGHT)
      {
        l = -light->direction;
        d = numeric_limits<double>::infinity();
      }
      else
      {
        l = (light->position - intersection.position);
        l.Normalize();
        d = (light->position - intersection.position).Length();
      }

      // construct a ray from the intersection toward the light
      R3Ray light_ray(offset_intersection_pos, l);

      // compute a shadow intersection and if one exists (that is between the intersection and the light), ignore the light
      R3Intersection shadow_intersection = ComputeIntersection(scene, scene->Root(), light_ray, d);
      if (shadow_intersection.hit && shadow_intersection.t < d)
        continue;

      // compute the luminance of the light at the intersection position
      R3Rgb light_il = ComputeLuminance(light, intersection.position);

      // evaluate the Phong BRDF equations for diffuse and specular reflection
      double cos_nl = l.Dot(intersection.normal);

      // r = l reflected across the intersection normal
      R3Vector r = 2 * cos_nl * intersection.normal - l; 
      double cos_vr = v.Dot(r);
      // if the light is shining on the front of the surface, include diffuse lighting
      // if the light is also less than 90 degrees from the camera vector, include specular lighting
      if (cos_nl > 0)
      {
        il += light_il * mat->kd * cos_nl;
        if (cos_vr > 0)
          il += light_il * mat->ks * pow(cos_vr, mat->shininess);
      }
    }
  }
  else // if we ARE using distributed...
  {
    if (depth < max_depth)
    {
      int nrays = 0;
      R3Rgb dist_il(0, 0, 0, 1);
      while (nrays < distrib)
      {
        // randomly sample outgoing directions using random variables with normal distributions
        double u1 = (double)rand()/RAND_MAX;
        double u2 = (double)rand()/RAND_MAX;
        double u3 = (double)rand()/RAND_MAX;
        double u4 = (double)rand()/RAND_MAX;

        double s1 = sqrt(-2 * log(u1));
        double rx = s1 * cos(2 * 3.14159 * u2);
        double ry = s1 * sin(2 * 3.14159 * u2);
        double s2 = sqrt(-2 * log(u3));
        double rz = s2 * cos(2 * 3.14159 * u4);

        R3Vector l(rx, ry, rz);
        l.Normalize();

        // check if vector is in the right hemisphere
        if (l.Dot(intersection.normal) < 0)
          continue;

        // generate a random ray and recursively compute the radiance it contributes
        R3Ray rray(offset_intersection_pos, l);
        R3Rgb rgb = ComputeRadiance(scene, rray, depth + 1, max_depth, distrib);
        
        // evaluate the Phong BRDF equations for diffuse and specular reflection
        double cos_nl = l.Dot(intersection.normal);

        // r = l reflected across the intersection normal
        R3Vector r = 2 * cos_nl * intersection.normal - l; 
        double cos_vr = v.Dot(r);
        // if the light is shining on the front of the surface, include diffuse lighting
        // if the light is also less than 90 degrees from the camera vector, include specular lighting
        if (cos_nl > 0)
        {
          dist_il += rgb * mat->kd * cos_nl;
          if (cos_vr > 0)
            dist_il += rgb * mat->ks * pow(cos_vr, mat->shininess);
        }

        nrays++;
      }

      il += dist_il / nrays;
    }
  }

  // if the material has a specular component and the maximum depth hasn't been reached, reflect the camera ray from the
  // intersection surface and recusively evaluate its radiance
  if (!mat->ks.IsBlack() && depth < max_depth)
  {
    R3Ray spec_ray = GetReflectedRay(intersection, ray);
    R3Rgb spec = ComputeRadiance(scene, spec_ray, depth + 1, max_depth, distrib);
    il += spec * mat->ks;
  }

  // add refraction component for both distributed and non-distributed techniques
  if (!mat->kt.IsBlack() && depth < max_depth)
  {
    R3Ray trans_ray = GetRefractedRay(intersection, ray, mat->indexofrefraction, scene);
    R3Rgb trans = ComputeRadiance(scene, trans_ray, depth + 1, max_depth, distrib);
    il += trans * mat->kt;
  }

  // add ambient and emission components
  il += mat->emission;
  il += mat->ka * scene->ambient;

  return il;
}

// Compute the luminance for a light at a specific point
R3Rgb ComputeLuminance(R3Light *light, R3Point &point)
{
  R3Rgb il(0, 0, 0, 1);

  switch (light->type)
  {
    // quadratic attenuation factor
    case R3_POINT_LIGHT:
    {
      double d = (light->position - point).Length();
      double att = light->constant_attenuation + light->linear_attenuation * d + light->quadratic_attenuation * d * d;
      il = light->color / att;
      break;
    }
    // constant luminance
    case R3_DIRECTIONAL_LIGHT:
    {
      il = light->color;
      break;
    }
    // quadratic attenuation with directional falloff
    case R3_SPOT_LIGHT:
    {
      R3Vector l = (point - light->position);
      l.Normalize();
      double cost = l.Dot(light->direction);
      if (acos(cost) <= light->angle_cutoff)
      {
        double d = (light->position - point).Length();
        double att = light->constant_attenuation + light->linear_attenuation * d + light->quadratic_attenuation * d * d;
        il = light->color * pow(cost, light->angle_attenuation) / att;
      }
      break;
    }
    // not implemented
    case R3_AREA_LIGHT:
      break;
    case R3_NUM_LIGHT_TYPES:
      break;
  }

  return il;
}

// test if a point is inside a bounding box
bool InsideBBox(R3Box box, R3Point point)
{
  return point.X() > box.XMin() && point.X() < box.XMax() && point.Y() > box.YMin() && point.Y() < box.YMax() && point.Z() > box.ZMin() && point.Z() < box.ZMax();
}

R3Intersection ComputeIntersection(R3Cylinder *cylinder, R3Ray r) {
  R3Intersection ret;
  R3Point p0 = r.Point(0);
  R3Vector v = r.Vector();
  double h = cylinder->Height();
  
  // check the two caps
  // Bottom First
  R3Plane bottom = R3negxz_plane;
  R3Point bottom_center = R3Point(0, -h/2, 0);
  bottom.Translate(bottom_center - R3null_point);
  
  double denom = r.Vector().Dot(bottom.Normal());
  double t = (denom == 0) ? -1 : -(R3Vector(p0[0], p0[1], p0[2]).Dot(bottom.Normal()) + bottom.D()) / denom;
  double dist = R3Distance(r.Point(t), bottom_center);
  if ((t > EPSILON) && (!ret.hit || t < ret.t) && (dist < cylinder->Radius())) {
    ret.t = t;
    ret.normal = bottom.Normal();
    ret.position = r.Point(t);
    ret.hit = true;
  }
  
  // Now Top
  R3Plane top = R3posxz_plane;
  R3Point top_center = R3Point(0, cylinder->Height()/2, 0);
  top.Translate(top_center - R3null_point);
  
  denom = r.Vector().Dot(top.Normal());
  t = (denom == 0) ? -1 : -(R3Vector(p0[0], p0[1], p0[2]).Dot(top.Normal()) + top.D()) / denom;
  dist = R3Distance(r.Point(t), top_center);
  if ((t > EPSILON) && (!ret.hit || t < ret.t) && (dist < cylinder->Radius())) {
    ret.t = t;
    ret.normal = top.Normal();
    ret.position = r.Point(t);
    ret.hit = true;
  }
  
  // Now check the side
  // we know (p0[0] + v[0]*t)^2 + (p0[0] + v[0]*t)^2 = r^2
  // so we solve the quadratic
  double A = v[0]*v[0] + v[2]*v[2];
  double B = 2*p0[0]*v[0] + 2*p0[2]*v[2];
  double C = p0[0]*p0[0] + p0[2]*p0[2] - cylinder->Radius()*cylinder->Radius();
  double det = B*B - 4 * A * C;
  if (det < 0) return ret;
  double t_plus = (-B + sqrt(det)) / (2*A);
  double t_minus = (-B - sqrt(det)) / (2*A);
  
  R3Point loc = r.Point(t_minus);
  if (((t_minus > EPSILON) && (!ret.hit || t_minus < ret.t)) &&
      ((loc[1] > -h/2) && (loc[1] < h/2)))
  {
    ret.t = t_minus;
    ret.position = loc;
    ret.normal = R3Point(loc[0], 0, loc[2]) - R3null_point;
    ret.normal.Normalize();
    ret.hit = true;
  }
  loc = r.Point(t_plus);
  if (((t_plus > EPSILON) && (!ret.hit || t_plus < ret.t)) &&
      ((loc[1] > -h/2) && (loc[1] < h/2))) {
    ret.t = t_plus;
    ret.position = loc;
    ret.normal = R3Point(loc[0], 0, loc[2]) - R3null_point;
    ret.normal.Normalize();
    ret.hit = true;
  }
  return ret;
}

// compute the nearest intersection between a ray and a scene node
R3Intersection ComputeIntersection(R3Scene *scene, R3Node *node, R3Ray ray, double min_t)
{
  // transform the ray and minimum t from parent node coordinates to this node's coordinates
  R3Vector v_orig = ray.Vector();
  v_orig.InverseTransform(node->transformation);
  double min_t_trans = min_t * v_orig.Length();
  ray.InverseTransform(node->transformation);

  // initialize intersection as infinitely far away
  R3Intersection min_intersection;
  min_intersection.hit = false;
  min_intersection.t = min_t_trans;

  // check this node if it contains a shape
  if (node->shape != NULL && node->is_obstacle)
  {
    R3Intersection i;
    i.hit = false;
    i.t = numeric_limits<double>::infinity();

    // intersect with node shape
    switch (node->shape->type)
    {
      case R3_BOX_SHAPE:
      {
        i = ComputeIntersection(node->shape->box, ray);
      }
      break;
      case R3_SPHERE_SHAPE:
      {
        i = ComputeIntersection(node->shape->sphere, ray);
      }
      break;
      case R3_MESH_SHAPE:
      {
        i = ComputeIntersection(node->shape->mesh, ray, min_intersection.t);
      }
      break;
      // not implemented
      case R3_COIN_SHAPE:
      case R3_CYLINDER_SHAPE:
        i = ComputeIntersection(node->shape->cylinder, ray);
      case R3_CONE_SHAPE:
        break;
      case R3_SEGMENT_SHAPE:
        break;
      case R3_CIRCLE_SHAPE:
        break;
      case R3_NUM_SHAPE_TYPES:
        break;
    }
    if (node->is_coin) {
      printf("I'm a coin");
      i = ComputeIntersection(&(node->bbox), ray);
    }
    // update minimum if found
    if (i.hit && i.t < min_intersection.t)
    {
      i.node = node;
      min_intersection = i;
    }
  }

  // recursively intersect to children nodes
  for (unsigned int i = 0; i < node->children.size(); i++)
  { 
    // check if ray intersects with the child node
    R3Intersection child_bbox_intersection = ComputeIntersection(&(node->children[i]->bbox), ray);
    // if true, intersect with this child node
    if (child_bbox_intersection.hit && child_bbox_intersection.t < min_intersection.t)
    {
      R3Intersection child_intersection = ComputeIntersection(scene, node->children[i], ray, min_intersection.t);
      // update intersection if necessary
      if (child_intersection.hit && child_intersection.t < min_intersection.t)
      {
        min_intersection = child_intersection;
      }
    }
  }

  // transform intersection from this node's coordinates to the parent node coordinates
  if (min_intersection.hit)
  {
    min_intersection.normal.InverseTransform(node->transformation.Transpose());
    min_intersection.normal.Normalize();
    min_intersection.position.Transform(node->transformation);
    R3Vector v2 = ray.Vector();
    v2.Transform(node->transformation);
    min_intersection.t = min_intersection.t * v2.Length();
  }
  return min_intersection;
}

// compute intersection between a sphere and a ray
R3Intersection ComputeIntersection(R3Sphere *sphere, R3Ray &ray)
{
  R3Intersection i;
  bool in_front_of_center = false;
  bool internal = false;

  R3Vector l = sphere->Center() - ray.Start();
  double tca = l.Dot(ray.Vector());
  if (tca < 0) // case where ray originates within sphere and points away from its center
  {
    tca = -tca;
    in_front_of_center = true;
    internal = true;
  }
  double d2 = l.Dot(l) - tca * tca;
  double r2 = sphere->Radius() * sphere->Radius();
  if (d2 > r2) // case where ray misses sphere entirely
  {
    i.hit = false;
    return i;
  }
  double thc = sqrt(r2 - d2);
  double t;
  if (!in_front_of_center)
  {
    t = tca - thc;
    if (t < 0)
    {
      t = tca + thc;
      internal = true;
    }
  }
  else
  {
    t = thc - tca;
  }

  if (t < 0)
  {
    i.hit = false;
    return i;
  }

  R3Point p = ray.Point(t);
  R3Vector n = p - sphere->Center();
  n.Normalize();
  if (internal)
    n = -n;

  // populate intersection data
  i.hit = true;
  i.normal = n;
  i.position = p;
  i.t = t;

  return i;
}

// helper functions for box intersect

bool CheckBoxX(R3Point &p, R3Box *box)
{
  return p.Y() <= box->YMax() && p.Y() >= box->YMin()
      && p.Z() <= box->ZMax() && p.Z() >= box->ZMin();
}

bool CheckBoxY(R3Point &p, R3Box *box)
{
  return p.X() <= box->XMax() && p.X() >= box->XMin()
      && p.Z() <= box->ZMax() && p.Z() >= box->ZMin();
}

bool CheckBoxZ(R3Point &p, R3Box *box)
{
  return p.X() <= box->XMax() && p.X() >= box->XMin()
      && p.Y() <= box->YMax() && p.Y() >= box->YMin();
}

// compute intersection between a box and a ray
R3Intersection ComputeIntersection(R3Box *box, R3Ray &ray)
{
  R3Intersection min_intersection;
  min_intersection.t = numeric_limits<double>::infinity();
  min_intersection.hit = false;

  // check cases where:
  // 1) ray originates on the positive side of the max and points negative
  // 2) ray originates on the negative side of the min and points positive
  // 3) ray originates within the box and points either way
  // for each of X, Y and Z

  // I wrote the function out like this to eliminate most unnecessary checks that would
  // come along with a "cleaner" implementation. This tradeoff between optimization
  // and clean code is a good one because this intersection function is used for bounding boxes,
  // a critical part of ray intersection.

  // X DIRECTION
  if (ray.Start().X() <= box->XMin())
  {
    if (ray.Vector().X() > 0)
    {
      double t = (box->XMin() - ray.Start().X()) / ray.Vector().X();
      R3Point p = ray.Point(t);
      if (t >= 0 && t < min_intersection.t && CheckBoxX(p, box))
      {
        min_intersection.hit = true;
        min_intersection.normal = R3Vector(-1, 0, 0);
        min_intersection.position = p;
        min_intersection.t = t;
      }
    }
    else
      return min_intersection;
  }
  else if (ray.Start().X() >= box->XMax())
  {
    if (ray.Vector().X() < 0)
    {
      double t = (box->XMax() - ray.Start().X()) / ray.Vector().X();
      R3Point p = ray.Point(t);
      if (t >= 0 && t < min_intersection.t && CheckBoxX(p, box))
      {
        min_intersection.hit = true;
        min_intersection.normal = R3Vector(1, 0, 0);
        min_intersection.position = p;
        min_intersection.t = t;
      }
    }
    else
      return min_intersection;
  }
  else
  {
    if (ray.Vector().X() > 0)
    {
      double t = (box->XMax() - ray.Start().X()) / ray.Vector().X();
      R3Point p = ray.Point(t);
      if (t >= 0 && t < min_intersection.t && CheckBoxX(p, box))
      {
        min_intersection.hit = true;
        min_intersection.normal = R3Vector(-1, 0, 0);
        min_intersection.position = p;
        min_intersection.t = t;
      }
    }
    else if (ray.Vector().X() < 0)
    {
      double t = (box->XMin() - ray.Start().X()) / ray.Vector().X();
      R3Point p = ray.Point(t);
      if (t >= 0 && t < min_intersection.t && CheckBoxX(p, box))
      {
        min_intersection.hit = true;
        min_intersection.normal = R3Vector(1, 0, 0);
        min_intersection.position = p;
        min_intersection.t = t;
      }
    }
  }

  // Y DIRECTION
  if (ray.Start().Y() <= box->YMin())
  {
    if (ray.Vector().Y() > 0)
    {
      double t = (box->YMin() - ray.Start().Y()) / ray.Vector().Y();
      R3Point p = ray.Point(t);
      if (t >= 0 && t < min_intersection.t && CheckBoxY(p, box))
      {
        min_intersection.hit = true;
        min_intersection.normal = R3Vector(0, -1, 0);
        min_intersection.position = p;
        min_intersection.t = t;
      }
    }
    else
      return min_intersection;
  }
  else if (ray.Start().Y() >= box->YMax())
  {
    if (ray.Vector().Y() < 0)
    {
      double t = (box->YMax() - ray.Start().Y()) / ray.Vector().Y();
      R3Point p = ray.Point(t);
      if (t >= 0 && t < min_intersection.t && CheckBoxY(p, box))
      {
        min_intersection.hit = true;
        min_intersection.normal = R3Vector(0, 1, 0);
        min_intersection.position = p;
        min_intersection.t = t;
      }
    }
    else
      return min_intersection;
  }
  else
  {
    if (ray.Vector().Y() > 0)
    {
      double t = (box->YMax() - ray.Start().Y()) / ray.Vector().Y();
      R3Point p = ray.Point(t);
      if (t >= 0 && t < min_intersection.t && CheckBoxY(p, box))
      {
        min_intersection.hit = true;
        min_intersection.normal = R3Vector(0, -1, 0);
        min_intersection.position = p;
        min_intersection.t = t;
      }
    }
    else if (ray.Vector().Y() < 0)
    {
      double t = (box->YMin() - ray.Start().Y()) / ray.Vector().Y();
      R3Point p = ray.Point(t);
      if (t >= 0 && t < min_intersection.t && CheckBoxY(p, box))
      {
        min_intersection.hit = true;
        min_intersection.normal = R3Vector(0, 1, 0);
        min_intersection.position = p;
        min_intersection.t = t;
      }
    }
  }

  // Z DIRECTION
  if (ray.Start().Z() <= box->ZMin())
  {
    if (ray.Vector().Z() > 0)
    {
      double t = (box->ZMin() - ray.Start().Z()) / ray.Vector().Z();
      R3Point p = ray.Point(t);
      if (t >= 0 && t < min_intersection.t && CheckBoxZ(p, box))
      {
        min_intersection.hit = true;
        min_intersection.normal = R3Vector(0, 0, -1);
        min_intersection.position = p;
        min_intersection.t = t;
      }
    }
    else
      return min_intersection;
  }
  else if (ray.Start().Z() >= box->ZMax())
  {
    if (ray.Vector().Z() < 0)
    {
      double t = (box->ZMax() - ray.Start().Z()) / ray.Vector().Z();
      R3Point p = ray.Point(t);
      if (t >= 0 && t < min_intersection.t && CheckBoxZ(p, box))
      {
        min_intersection.hit = true;
        min_intersection.normal = R3Vector(0, 0, 1);
        min_intersection.position = p;
        min_intersection.t = t;
      }
    }
    else
      return min_intersection;
  }
  else
  {
    if (ray.Vector().Z() > 0)
    {
      double t = (box->ZMax() - ray.Start().Z()) / ray.Vector().Z();
      R3Point p = ray.Point(t);
      if (t >= 0 && t < min_intersection.t && CheckBoxZ(p, box))
      {
        min_intersection.hit = true;
        min_intersection.normal = R3Vector(0, 0, -1);
        min_intersection.position = p;
        min_intersection.t = t;
      }
    }
    else if (ray.Vector().Z() < 0)
    {
      double t = (box->ZMin() - ray.Start().Z()) / ray.Vector().Z();
      R3Point p = ray.Point(t);
      if (t >= 0 && t < min_intersection.t && CheckBoxZ(p, box))
      {
        min_intersection.hit = true;
        min_intersection.normal = R3Vector(0, 0, 1);
        min_intersection.position = p;
        min_intersection.t = t;
      }
    }
  }

  return min_intersection;
}

// compute intersection between a mesh and a ray
R3Intersection ComputeIntersection(R3Mesh *mesh, R3Ray &ray, double min_t)
{
  R3Intersection min_intersection;
  min_intersection.t = min_t;
  min_intersection.hit = false;

  // check each face
  for (int i = 0; i < mesh->NFaces(); i++)
  {
    R3MeshFace *face = mesh->Face(i);
    assert (face->vertices.size() >= 3);
    
    // if the face has more than 3 vertices, triangulate it
    if (face->vertices.size() > 3)
    {
      for (unsigned int j = 1; j < face->vertices.size()-1; j++)
      {
        // create a virtual face
        vector<R3MeshVertex *> vertices;
        vertices.push_back(face->vertices[0]);
        vertices.push_back(face->vertices[j]);
        vertices.push_back(face->vertices[j + 1]);
        R3MeshFace new_face(vertices);
        new_face.UpdatePlane();

        // compute intersection with triangle
        R3Intersection intersection = ComputeIntersection(&new_face, ray, min_intersection.t);
        if (intersection.hit && intersection.t < min_intersection.t)
          min_intersection = intersection;
      }
    } else {
      // compute intersection with triangle
      R3Intersection intersection = ComputeIntersection(face, ray, min_intersection.t);
      if (intersection.hit && intersection.t < min_intersection.t)
        min_intersection = intersection;
    }
  }

  return min_intersection;
}

// compute intersection between a triangle face and a ray
R3Intersection ComputeIntersection(R3MeshFace *tri, R3Ray &ray, double min_t)
{
  assert (tri->vertices.size() == 3);

  R3Intersection i;

  R3Point t1 = tri->vertices[0]->position;
  R3Point t2 = tri->vertices[1]->position;
  R3Point t3 = tri->vertices[2]->position;

  // intersect ray with triangle's plane
  R3Plane plane = tri->plane;
  double t = -(ray.Start().Vector().Dot(plane.Normal()) + plane.D()) 
    / (ray.Vector().Dot(plane.Normal()));

  // early return if not closer than minimum intersection for the mesh
  if (t > min_t || t < 0)
  {
    i.hit = false;
    return i;
  }


  // check if intersection is within triangle using barycentric coordinate method
  R3Point p = ray.Point(t);

  R3Vector v;

  v = t2 - t1;
  v.Cross(t3 - t1);
  double area = v.Length() / 2;

  v = t2 - t1;
  v.Cross(p - t1);
  if (v.Dot(plane.Normal()) < 0)
  {
    i.hit = false;
    return i;
  }
  double a = v.Length() / (2 * area);

  v = p - t1;
  v.Cross(t3 - t1);
  if (v.Dot(plane.Normal()) < 0)
  {
    i.hit = false;
    return i;
  }
  double b = v.Length() / (2 * area);

  if (a <= 1 && a >= 0 && b <= 1 && b >= 0 && a + b <= 1)
  {
    i.hit = true;
    if (ray.Vector().Dot(plane.Normal()) < 0)
    {
      i.normal = plane.Normal();
    }
    else
    {
      i.normal = -plane.Normal();
    }
    i.position = p;
    i.t = t;
    return i;
  } 
  else
  {
    i.hit = false;
    return i;
  }
}

// return point of intersection between ray and plane
R3Point RayPlaneIntersection(R3Plane plane, R3Ray ray)
{
  double t = -(ray.Start().Vector().Dot(plane.Normal()) + plane.D()) 
    / (ray.Vector().Dot(plane.Normal()));
  return ray.Point(t);
}

// compute area between three points
double Area(R3Point p1, R3Point p2, R3Point p3)
{
  R3Vector v = p2 - p1;
  v.Cross(p3 - p1);
  return v.Length() / 2;
}