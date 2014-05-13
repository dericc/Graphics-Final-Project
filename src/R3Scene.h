// Include file for the R3 scene stuff

#define R3Rgb R2Pixel
#define PI 3.1415926535

// Constant definitions

typedef enum {
  R3_BOX_SHAPE,
  R3_SPHERE_SHAPE,
  R3_CYLINDER_SHAPE,
  R3_CONE_SHAPE,
  R3_MESH_SHAPE,
  R3_SEGMENT_SHAPE,
  R3_CIRCLE_SHAPE,
  R3_NUM_SHAPE_TYPES,
  R3_COIN_SHAPE,
} R3ShapeType;

typedef enum {
  R3_DIRECTIONAL_LIGHT,
  R3_POINT_LIGHT,
  R3_SPOT_LIGHT,
  R3_AREA_LIGHT,
  R3_NUM_LIGHT_TYPES
} R3LightType;



// Scene element definitions

struct R3Shape {
  R3ShapeType type;
  R3Box *box;
  R3Sphere *sphere;
  R3Cylinder *cylinder;
  R3Cone *cone;
  R3Mesh *mesh;
  R3Segment *segment;
  R3Circle *circle;
};  

struct R3Material {
  R3Rgb ka;
  R3Rgb kd;
  R3Rgb ks;
  R3Rgb kt;
  R3Rgb emission;
  double shininess;
  double indexofrefraction;
  R2Image *texture;
  int texture_index;
  int id;
  char texture_name[256]; 
};

struct R3Light {
  R3LightType type;
  R3Point position;
  R3Vector direction;
  double radius;
  R3Rgb color;
  double constant_attenuation;
  double linear_attenuation;
  double quadratic_attenuation;
  double angle_attenuation;
  double angle_cutoff;
};

struct R3Camera {
  R3Point eye;
  R3Vector towards;
  R3Vector right;
  R3Vector up;
  double xfov, yfov;
  double neardist, fardist;
  
  void Rotate(R3Line axis, double angle);
};

struct R3Coin;
struct R3Platform;
struct R3Enemy; 

struct R3Node {
  R3Node(void)
  : is_obstacle(false), is_coin(false), is_platform(false), is_enemy(false), is_visible(true), is_goal(false), is_player(false), del(false) {};
  struct R3Node *parent;
  vector<struct R3Node *> children;
  R3Shape *shape;
  R3Matrix transformation;
  R3Material *material;
  R3Box bbox;
  bool is_obstacle;
  bool is_coin;
  R3Coin *coin;
  bool is_platform;
  R3Platform *platform;
  bool is_enemy; 
  R3Enemy *enemy; 
  bool is_visible;
  bool is_goal;
  bool is_player;
  bool del;
};

struct R3Platform {
  R3Platform(R3Node *node, double speed, R3Point start, R3Point end)
  : node(node), max_speed(speed), center((start+end)/2), start(start), end(end) {  };
  R3Node *node;
  const double max_speed;
  const R3Point center; //center of the path
  const R3Point start; 
  const R3Point end;  
  R3Vector velocity;
  
  R3Vector Forward(void); // normalized forward direction
};

// Particle system definitions

struct R3Particle {
  R3Point position;
  R3Vector velocity;
  double mass;
  bool fixed;
  double drag;
  double elasticity;
  double lifetime;
  R3Material *material;
  vector<struct R3ParticleSpring *> springs;
  vector<R3Point> trail;
};

struct R3ParticleSource {
  R3Shape *shape;
  double rate;
  double velocity;
  double angle_cutoff;
  double mass;
  bool fixed;
  double drag;
  double elasticity;
  double lifetime;
  R3Material *material;
};

struct R3ParticleSink {
  R3Shape *shape;
  double intensity;
  double constant_attenuation;
  double linear_attenuation;
  double quadratic_attenuation;
};

struct R3ParticleSpring {
  R3Particle *particles[2];
  double rest_length;
  double ks;
  double kd;
};

struct R3Fire {
  R3Point position;
};

struct R3Player {
  R3Player(R3Node *node, double max_speed, double mass) :
    node(node), max_speed(max_speed),  mass(mass), is_dead(false), has_won(false), n_coins(0), won_time(0.0f), onPlatform(false) {};
  
  R3Node *node; //
  const double max_speed;
  const double mass;
  
  R3Point Center();
  R3Vector Right();
  R3Vector Towards();
  R3Vector Up();
  
  void Jump();
  void MoveLeft();
  void MoveRight();
  void MoveUp();
  
  R3Vector velocity; // current direction of motion
  bool inAir;
  bool is_dead;
  bool has_won;
  int n_coins;
  double won_time;
  
  bool onPlatform;
  R3Platform *platform;
};

struct R3Coin {
  R3Node *node;
  R3Point position;
  double t;
  bool del;
};

struct R3Enemy {
  R3Enemy(R3Node *node, bool moveLeft, double speed, bool is_jumping, bool is_following, double jumpHeight, double mass) :
    node(node), moveLeft(moveLeft), speed(speed), 
    is_dead(false), del(false), is_jumping(is_jumping), is_following(is_following), jumpHeight(jumpHeight), mass(mass), onPlatform(false) {};
  
  R3Node *node; 
  bool moveLeft; // current direction of motion: left or right
  double speed; 
  R3Vector velocity; // current direction of motion
  
  R3Point Center();
  R3Vector Right();
  R3Vector Towards();
  R3Vector Up();

  bool inAir;
  bool is_dead;
  bool del; 

  bool is_jumping; 
  bool is_following; 

  double jumpHeight; 

  const double mass;
  
  bool onPlatform;
  R3Platform *platform;
};

// goal for player to get to
struct R3Goal {
  R3Goal(R3Node *node) :
    node(node), is_active(true) {};

  R3Point Center();

  R3Node *node;
  bool is_active; // is the goal active?
};

struct R3Sidebar;

// Scene graph definition

struct R3Scene {
 public:
  // Constructor functions
  R3Scene(void);

  // Access functions
  R3Node *Root(void) const;
  int NLights(void) const;
  R3Light *Light(int k) const;
  R3Camera& Camera(void);
  R3Box& BBox(void);

  // Particle stuff
  int NParticles(void) const;
  R3Particle *Particle(int k) const;
  int NParticleSources(void) const;
  R3ParticleSource *ParticleSource(int k) const;
  int NParticleSinks(void) const;
  R3ParticleSink *ParticleSink(int k) const;
  int NParticleSprings(void) const;
  R3ParticleSpring *ParticleSpring(int k) const;
  int NCoins(void) const;
  R3Coin *Coin(int k) const;

  // I/O functions
  int Read(const char *filename, R3Node *root = NULL);

  void WritePlayer(FILE *fp); 
  void WriteFires(FILE *fp); 
  void WriteEnemies(FILE *fp); 
  void WriteGoal(FILE *fp); 
  void WriteMaterials(FILE *fp); 
  void WriteLights(FILE *fp); 
  void WritePlatforms(FILE *fp); 
  void WriteCoins(FILE *fp); 
  void WriteNode(FILE *fp, R3Node *node); 
  void WriteSkybox(FILE *fp);
  void WriteSoundtrack(FILE *fp);
  void WriteNextLevel(FILE *fp);
  int Write(const char *filename, R3Node *node); 

 public:
  int death_y;
  R3Node *root;
  vector<R3Particle *> particles;
  vector<R3Particle *> fire_particles;
  vector<R3ParticleSource *> particle_sources;
  vector<R3ParticleSink *> particle_sinks;
  vector<R3ParticleSpring *> particle_springs;
  vector<R3Coin *> coins;
  vector<R3Platform *> platforms;
  vector<R3Light *> lights;
  vector<R3Enemy *> enemies;
  vector<R3Fire *> fires;
  vector<R3Material *> materials;
  R3Vector gravity;
  R3Camera camera;
  R3Box bbox;
  R3Rgb background;
  R3Rgb ambient;
  R3Player *player;
  R3Goal *goal;
  R3Sidebar *sidebar;
  R3Plane movement_plane;
  R3Shape *coin_shape;
  R3Material *coin_material;
  char skyboxTexture[256]; 
  char soundtrack[256];
  char next_level[256];
};

typedef void (*button_fxn)(void);

struct R3Button {
  // this might be ugly but one of these two guys will be null
  R3Button(int *value, R3Material *material) : value(value), f(NULL), material(material) {};
  R3Button(button_fxn f, R3Material *material) : value(NULL), f(f), material(material) {};
  
  int * value;
  button_fxn f;
  R3Material *material;
};

struct R3Sidebar {
  R3Sidebar(double width, double border) :
  width(width), border(border), button_width(width - 2*border), selected_button(-1) {};
  vector<R3Button *> buttons;
  double width;
  double border;
  double button_width;
  int selected_button;
};

// Inline functions 

inline R3Node *R3Scene::
Root(void) const
{
  // Return root node
  return root;
}



inline int R3Scene::
NLights(void) const
{
  // Return number of lights
  return lights.size();
}



inline R3Light *R3Scene::
Light(int k) const
{
  // Return kth light
  return lights[k];
}



inline R3Camera& R3Scene::
Camera(void) 
{
  // Return camera
  return camera;
}



inline R3Box& R3Scene::
BBox(void) 
{
  // Return bounding box 
  return bbox;
}



inline int R3Scene::
NParticles(void) const
{
  // Return number of particles
  return particles.size();
}



inline R3Particle *R3Scene::
Particle(int k) const
{
  // Return kth particle 
  return particles[k];
}



inline int R3Scene::
NParticleSources(void) const
{
  // Return number of particle sources
  return particle_sources.size();
}



inline R3ParticleSource *R3Scene::
ParticleSource(int k) const
{
  // Return kth particle source
  return particle_sources[k];
}



inline int R3Scene::
NParticleSinks(void) const
{
  // Return number of particle sinks
  return particle_sinks.size();
}



inline R3ParticleSink *R3Scene::
ParticleSink(int k) const
{
  // Return kth particle sink
  return particle_sinks[k];
}



inline int R3Scene::
NParticleSprings(void) const
{
  // Return number of particle springs
  return particle_springs.size();
}



inline R3ParticleSpring *R3Scene::
ParticleSpring(int k) const
{
  // Return kth particle spring
  return particle_springs[k];
}



inline int R3Scene::
NCoins(void) const
{
  // Return number of particle springs
  return coins.size();
}



inline R3Coin *R3Scene::
Coin(int k) const
{
  // Return kth particle spring
  return coins[k];
}




