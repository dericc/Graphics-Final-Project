// Include file for ray tracing code

struct R3Intersection
{
  bool hit;
  R3Node *node;
  R3Point position;
  R3Vector normal;
  double t;
};

R2Image *RenderImage(R3Scene *scene, int width, int height, int max_depth, 
  int num_primary_rays_per_pixel, int num_distributed_rays_per_intersection);

R3Ray RayThoughPixel(R3Camera camera, int i, int j, int width, int height);

R3Rgb ComputeLuminance(R3Light *light, R3Point &point);

R3Rgb ComputeRadiance(R3Scene *scene, R3Ray &ray, int depth, int max_depth, int distrib);
R3Rgb ComputeRadiance(R3Scene *scene, R3Ray &ray, R3Intersection &intersection, int depth, int max_depth, int distrib);

R3Ray GetReflectedRay(R3Intersection i, R3Ray r);
R3Ray GetRefractedRay(R3Intersection i, R3Ray r, double ior, R3Scene *scene);

bool InsideBBox(R3Box box, R3Point point);
R3Intersection ComputeIntersection(R3Scene *scene, R3Node *node, R3Ray ray, double min_t);
R3Intersection ComputeIntersection(R3Sphere *sphere, R3Ray &ray);
R3Intersection ComputeIntersection(R3Box *box, R3Ray &ray);
R3Intersection ComputeIntersection(R3MeshFace *tri, R3Ray &ray, double min_t);
R3Intersection ComputeIntersection(R3Mesh *mesh, R3Ray &ray, double min_t);

double Area(R3Point p1, R3Point p2, R3Point p3);