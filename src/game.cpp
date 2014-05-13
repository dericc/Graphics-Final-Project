// Source file for the scene file viewer


////////////////////////////////////////////////////////////
// INCLUDE FILES
////////////////////////////////////////////////////////////

#include "R3/R3.h"
#include "R3Scene.h"
#include "particle.h"
#include "fglut/fglut.h"

#include <iostream>
#include <unistd.h>
#include <errno.h>
#include "../irrKlang/include/irrKlang.h"
#include "raytrace.h"


using namespace irrklang;

#define TIME_SCALE 2
#define TIME_STEP 0.002

////////////////////////////////////////////////////////////
// GLOBAL CONSTANTS
////////////////////////////////////////////////////////////

static const double VIDEO_FRAME_DELAY = 1./25.; // 25 FPS 

////////////////////////////////////////////////////////////
// GLOBAL VARIABLES
////////////////////////////////////////////////////////////

// Program arguments

static char *input_scene_name = NULL;
static char *output_image_name = NULL;
static const char *video_prefix = "../video-frames/";
static const char *images_path = "../images/";
static const char *skybox_path = "../levels/";
static int integration_type = EULER_INTEGRATION;
static double previous_time = 0;

// Display variables

static R3Scene *scene = NULL;
static R3Camera camera;
static R3Camera minimap_cam;
static R3Node *grabbed = NULL;
static int show_faces = 1;
static int show_edges = 0;
static int show_bboxes = 0;
static int show_lights = 0;
static int show_camera = 0;
static int show_particles = 1;
static int show_particle_springs = 1;
static int show_particle_sources_and_sinks = 1;
static int save_image = 0;
static int save_video = 0;
static int num_frames_to_record = -1; 
static int quit = 0;
static int minimap = 1;
static int level_editor = 0;
static int soundtrack_enabled = 0;
static int camera_mode = 0;
static int blocks_mode = 0;
static int grab_mode = 0;
static int move_mode = 0;
static int coins_mode = 0;
static int paused = 0;
// GLUT variables 

static int GLUTwindow = 0;
static int GLUTwindow_height = 512;
static int GLUTwindow_width = 512;
static int GLUTmouse[2] = { 0, 0 };
static int GLUTbutton[3] = { 0, 0, 0 };
static int GLUTmodifiers = 0;

vector<R3Material *> skyboxMaterials; 

// GLUT command list

enum {
  DISPLAY_FACE_TOGGLE_COMMAND,
  DISPLAY_EDGE_TOGGLE_COMMAND,
  DISPLAY_BBOXES_TOGGLE_COMMAND,
  DISPLAY_LIGHTS_TOGGLE_COMMAND,
  DISPLAY_CAMERA_TOGGLE_COMMAND,
  DISPLAY_PARTICLES_TOGGLE_COMMAND,
  DISPLAY_PARTICLE_SPRINGS_TOGGLE_COMMAND,
  DISPLAY_PARTICLE_SOURCES_AND_SINKS_TOGGLE_COMMAND,
  SAVE_IMAGE_COMMAND,
  SAVE_VIDEO_COMMAND,
  QUIT_COMMAND,
};

// ascii key state

static bool key_state[256] = {false};

// sound engine
static ISoundEngine *sound_engine;
static ISound *soundtrack;

// executable path
static char exec_path[FILENAME_MAX + 1];

void LoadLevel(const char *filename);


////////////////////////////////////////////////////////////
// LEVEL EDITOR STUFF
////////////////////////////////////////////////////////////

R3Camera GetMinimapCam(R3Scene *scene);
void DrawSidebar(R3Scene *);

void CreateShape(R3ShapeType type, R3Scene *s, R3Point p)
{
  if (type == R3_BOX_SHAPE)
  {
    // Create box
    R3Box *box = new R3Box(p + R3Vector(-2.5, -.5, -2), p + R3Vector(2.5, .5, 2));
    
    // Create shape
    R3Shape *shape = new R3Shape();
    shape->type = R3_BOX_SHAPE;
    shape->box = box;
    shape->sphere = NULL;
    shape->cylinder = NULL;
    shape->cone = NULL;
    shape->mesh = NULL;
    shape->segment = NULL;
    
    // Create shape node
    R3Node *node = new R3Node();
    node->transformation = R3identity_matrix;
    node->material = NULL;
    node->shape = shape;
    node->bbox = *box;
    node->is_obstacle = true;
    node->is_visible = true;
    node->is_coin = false;

    // Add to scene
    s->root->bbox.Union(node->bbox);
    s->root->children.push_back(node);
    node->parent = s->root;
  }
  if (type == R3_COIN_SHAPE) {
    R3Coin *coin = new R3Coin();
    coin->position = p;
    coin->t = 0;
    coin->del = false;
    
    R3Matrix tform = R3identity_matrix;
    tform.Translate(p.Vector());
    tform.Rotate(R3_X, PI / 2);
    
    // Create shape node
    R3Node *node = new R3Node();
    node->transformation = tform;
    node->material = s->coin_material;
    node->shape = s->coin_shape;
    node->bbox = s->coin_shape->cylinder->BBox();
    node->is_coin = true;
    node->coin = coin;
    
    coin->node = node;
    
    s->coins.push_back(coin);
    
    s->root->bbox.Union(node->bbox);
    s->root->children.push_back(node);
    node->parent = s->root;
  }
  minimap_cam = GetMinimapCam(scene);
}

////////////////////////////////////////////////////////////
// TIMER CODE
////////////////////////////////////////////////////////////

#ifdef _WIN32
#  include <windows.h>
#else
#  include <sys/time.h>
#endif

#ifdef __MACH__
#  include <mach-o/dyld.h>
#endif

static double GetTime(void)
{
#ifdef _WIN32
  // Return number of seconds since start of execution
  static int first = 1;
  static LARGE_INTEGER timefreq;
  static LARGE_INTEGER start_timevalue;

  // Check if this is the first time
  if (first) {
    // Initialize first time
    QueryPerformanceFrequency(&timefreq);
    QueryPerformanceCounter(&start_timevalue);
    first = 0;
    return 0;
  }
  else {
    // Return time since start
    LARGE_INTEGER current_timevalue;
    QueryPerformanceCounter(&current_timevalue);
    return ((double) current_timevalue.QuadPart - 
            (double) start_timevalue.QuadPart) / 
            (double) timefreq.QuadPart * TIME_SCALE;
  }
#else
  // Return number of seconds since start of execution
  static int first = 1;
  static struct timeval start_timevalue;

  // Check if this is the first time
  if (first) {
    // Initialize first time
    gettimeofday(&start_timevalue, NULL);
    first = 0;
    return 0;
  }
  else {
    // Return time since start
    struct timeval current_timevalue;
    gettimeofday(&current_timevalue, NULL);
    int secs = current_timevalue.tv_sec - start_timevalue.tv_sec;
    int usecs = current_timevalue.tv_usec - start_timevalue.tv_usec;
    return (double) (secs + 1.0E-6F * usecs) * TIME_SCALE;
  }
#endif
}

////////////////////////////////////////////////////////////
// SOUND FUNCTIONS
////////////////////////////////////////////////////////////

void GetExecPath()
{
  #ifdef __linux
    if (readlink("/proc/self/exe", exec_path, FILENAME_MAX + 1) == -1)
    {
      exit(1);
    }
  #else
    uint32_t size = sizeof(exec_path);
    if (_NSGetExecutablePath(exec_path, &size) != 0)
    {
      exit(1);
    }
  #endif
}

void FilePath(char *buf, const char *filename)
{
  strncpy(buf, exec_path, strlen(exec_path) - 5);
  strncpy(buf + strlen(exec_path) - 5, filename, strlen(filename));
  buf[strlen(exec_path) + strlen(filename) - 5] = '\0';
}

void PlaySound(const char *filename, bool looped)
{
  char path[FILENAME_MAX + 1];
  FilePath(path, filename);
  sound_engine->play2D(path, looped);
}



////////////////////////////////////////////////////////////
// SCENE DRAWING CODE
////////////////////////////////////////////////////////////

void KillPlayer(void) {
  R3Player *p = scene->player;
  
  if (!p->is_dead) {
    p->is_dead = true;
    p->node->is_visible = false;
    CreateParticles(scene, p->Center(), 10, scene->materials[0], 5);
    PlaySound("/../sounds/death.wav", false);
  }
}

void KillEnemy(R3Enemy *e) {

  if (!e->is_dead) {
    e->is_dead = true;
    e->node->is_visible = false;
    e->del = true; 
    e->node->del = true;
    // PlaySound("/../sounds/death.wav", false);
  }
}

enum
{
  COLLISION_LEFT,
  COLLISION_TOP,
  COLLISION_RIGHT,
  COLLISION_BOTTOM,
  COLLISION_NONE
};

// returns collision direction in relation to box A
int CollideBoxes(R3Box box_a, R3Box box_b)
{
  // calculate to what depth the player is intersecting in each direction
  double xmin_d = abs(box_b.XMax() - box_a.XMin());
  double xmax_d = abs(box_b.XMin() - box_a.XMax());
  double ymin_d = abs(box_b.YMax() - box_a.YMin());
  double ymax_d = abs(box_b.YMin() - box_a.YMax());

  bool zcoll = (box_a.ZMin() <= box_b.ZMax() && box_a.ZMin() >= box_b.ZMin())
            || (box_a.ZMax() <= box_b.ZMax() && box_a.ZMax() >= box_b.ZMin())
            || (box_a.ZMax() >= box_b.ZMax() && box_a.ZMin() <= box_b.ZMin());

  bool ycoll = (box_a.YMin() <= box_b.YMax() && box_a.YMin() >= box_b.YMin())
            || (box_a.YMax() <= box_b.YMax() && box_a.YMax() >= box_b.YMin())
            || (box_a.YMax() >= box_b.YMax() && box_a.YMin() <= box_b.YMin());

  bool xcoll = (box_a.XMin() <= box_b.XMax() && box_a.XMin() >= box_b.XMin())
            || (box_a.XMax() <= box_b.XMax() && box_a.XMax() >= box_b.XMin())
            || (box_a.XMax() >= box_b.XMax() && box_a.XMin() <= box_b.XMin());

  // determine if there is an intersection for each direction (player edge is inside object)
  bool xmin_coll = zcoll && ycoll && box_a.XMin() <= box_b.XMax() && box_a.XMin() >= box_b.XMin();

  bool xmax_coll = zcoll && ycoll && box_a.XMax() <= box_b.XMax() && box_a.XMax() >= box_b.XMin();

  bool ymin_coll = zcoll && xcoll && box_a.YMin() <= box_b.YMax() && box_a.YMin() >= box_b.YMin();

  bool ymax_coll = zcoll && xcoll && box_a.YMax() <= box_b.YMax() && box_a.YMax() >= box_b.YMin();

  // for each direction, check if it has the minimum depth of all other collisions before correcting player position and setting velocity to 0
  if (xmin_coll && (!ymin_coll || xmin_d < ymin_d) && (!ymax_coll || xmin_d < ymax_d))
    return COLLISION_LEFT;
  if (xmax_coll && (!ymin_coll || xmax_d < ymin_d) && (!ymax_coll || xmax_d < ymax_d))
    return COLLISION_RIGHT;
  if (ymin_coll && (!xmin_coll || ymin_d < xmin_d) && (!xmax_coll || ymin_d < xmax_d))
    return COLLISION_BOTTOM;
  if (ymax_coll && (!xmin_coll || ymax_d < xmin_d) && (!xmax_coll || ymax_d < xmax_d))
    return COLLISION_TOP;

  return COLLISION_NONE;
}

// recursively collide player with the scene
void CollidePlayer(R3Node *node)
{
  // get transformed player box
  R3Player *p = scene->player;
  if (node == p->node)
    return;
  R3Box player_box = *(p->node->shape->box);
  player_box.Transform(p->node->transformation);

  // if the scene node is a box
  if (node->is_obstacle || node->is_coin)
  {
    // get transformed scene box
    R3Box scene_box;
    if (node->is_obstacle)
    {
      scene_box = node->bbox;
      scene_box.Transform(node->transformation);
    }
    else if (node->is_coin)
    {
      scene_box = R3Box(R3Point(-0.5, -0.5, -0.5), R3Point(0.5, 0.5, 0.5));
      scene_box.Translate(node->coin->position.Vector());
    }

    int coll_dir = CollideBoxes(player_box, scene_box);

    if (node->is_obstacle)
    {
      // get player transformation
      R3Matrix tform = p->node->transformation;

      // for each direction, check if it has the minimum depth of all other collisions before correcting player position and setting velocity to 0
      if (coll_dir == COLLISION_LEFT)
      {
        if (node->is_enemy) {
          KillPlayer();  
        }

        tform.Translate(R3Vector(scene_box.XMax() - player_box.XMin(), 0, 0));
        R3Vector v = p->velocity;
        v.SetX(0);
        p->velocity = v;
      }
      if (coll_dir == COLLISION_RIGHT)
      {
        if (node->is_enemy) {
          KillPlayer();  
        }

        tform.Translate(R3Vector(scene_box.XMin() - player_box.XMax(), 0, 0));
        R3Vector v = p->velocity;
        v.SetX(0);
        p->velocity = v;
      }

      if (coll_dir == COLLISION_BOTTOM)
      {
        //If collision with enemy from above, kill the enemy
        if (node->is_enemy) {
          KillEnemy(node->enemy);  
        }

        tform.Translate(R3Vector(0, scene_box.YMax() - player_box.YMin(), 0));
        R3Vector v = p->velocity;
        v.SetY(0);
        
        //Jump up if node is enemy; else, just regular collision
        if (node ->is_enemy) {
          p->velocity = v + 10 * p->Up();
        }
        else {
          p->velocity = v;
        }

        p->inAir = false;
        if (node->is_platform && !node->is_enemy) {
          p->onPlatform = true;
          p->platform = node->platform;
        }
        else {
          p->onPlatform = false;
          p->platform = NULL;
        }
      }
      if (coll_dir == COLLISION_TOP)
      {
        if (node->is_enemy) {
          KillPlayer();  
        }

        tform.Translate(R3Vector(0, scene_box.YMin() - player_box.YMax(), 0));
        R3Vector v = p->velocity;
        v.SetY(0);
        p->velocity = v;
      }
      // update player tform
      p->node->transformation = tform;
    }
    else if (node->is_coin)
    {
      if (coll_dir != COLLISION_NONE)
      {
        PlaySound("/../sounds/coin.wav", false);
        p->n_coins++;
        node->del = true;
        node->coin->del = true;
      }
    }
  }

  for (unsigned int i = 0; i < node->children.size(); i++)
  {
    CollidePlayer(node->children[i]);
  }
}


// Collide enemies with different objects. 
void CollideEnemy(R3Enemy *e, R3Node *node)
{
  if (node == e->node)
    return;
  R3Box enemy_box = *(e->node->shape->box);
  enemy_box.Transform(e->node->transformation);

  // if the scene node is a box
  if (node->is_obstacle || node->is_coin)
  {
    // get transformed scene box
    R3Box scene_box;
    if (node->is_obstacle)
    {
      scene_box = node->bbox;
      scene_box.Transform(node->transformation);
    }
    else if (node->is_coin)
    {
      scene_box = R3Box(R3Point(-0.5, -0.5, -0.5), R3Point(0.5, 0.5, 0.5));
      scene_box.Translate(node->coin->position.Vector());
    }

    int coll_dir = CollideBoxes(enemy_box, scene_box);

    if (node->is_obstacle)
    {
      // get enemy transformation
      R3Matrix tform = e->node->transformation;

      // for each direction, check if it has the minimum depth of all other collisions before correcting enemy position and setting velocity to 0
      if (coll_dir == COLLISION_LEFT)
      {
        tform.Translate(R3Vector(scene_box.XMax() - enemy_box.XMin(), 0, 0));
        R3Vector v = e->velocity;
        v.SetX(0);
        e->velocity = v;
      }

      if (coll_dir == COLLISION_RIGHT)
      {
        tform.Translate(R3Vector(scene_box.XMin() - enemy_box.XMax(), 0, 0));
        R3Vector v = e->velocity;
        v.SetX(0);
        e->velocity = v;
      }

      if (coll_dir == COLLISION_BOTTOM)
      {
        tform.Translate(R3Vector(0, scene_box.YMax() - enemy_box.YMin(), 0));
        R3Vector v = e->velocity;
        v.SetY(0);
        e->velocity = v;
        e->inAir = false;
        if (node->is_platform) {
          e->onPlatform = true;
          e->platform = node->platform;
        }
      }
      if (coll_dir == COLLISION_TOP)
      {
        tform.Translate(R3Vector(0, scene_box.YMin() - enemy_box.YMax(), 0));
        R3Vector v = e->velocity;
        v.SetY(0);
        e->velocity = v;
      }
      // update enemy tform
      e->node->transformation = tform;
    }
  }

  for (unsigned int i = 0; i < node->children.size(); i++)
  {
    CollideEnemy(e, node->children[i]);
  }
}

//Update the player's location based on key presses
void UpdatePlayer(R3Scene *scene, double delta_time) {
  R3Player *p = scene->player;

  if (p == NULL) return;

  // program just started up?
  if (delta_time == 0) {
    p->velocity = R3null_vector;
    p->onPlatform = false;
    p->inAir = true;
  }

  bool up_key = key_state['w'] || key_state['W'];
  bool down_key = key_state['s'] || key_state['S'];
  bool left_key = key_state['a'] || key_state['A'];
  bool right_key = key_state['d'] || key_state['D'];

  // check if they've won
  if (scene->goal)
  {
    R3Goal *goal = scene->goal;
    double goal_dist = R3Distance(p->Center(), goal->Center());
    if (!p->has_won && goal_dist < 0.15f && scene->NCoins() <= 0) {
      PlaySound("/../sounds/victory.wav", false);
      p->won_time = GetTime();
      p->has_won = true;
    }
  }

  // Motion Shit
  // get the forces to move the box
  R3Vector f = R3null_vector;
  if (p->inAir) {
    f += -9.8 * p->Up() * p->mass;
  }
  if (up_key && !p->inAir && !p->has_won) {
    p->velocity += 15 * p->Up();
    PlaySound("/../sounds/jump.wav", false);
    if (p->onPlatform) {
      // only in the horizontal direction;
      R3Vector platformForwardVelocity = p->platform->velocity;
      platformForwardVelocity.Project(p->Towards());
      p->velocity += platformForwardVelocity;
    }
    p->inAir = true;
    p->onPlatform = false;
    p->platform = NULL;
  }

  // side to side
  double TAU = .3; // timescale for velocity relaxation
  const R3Vector forward = p->Towards();
  R3Vector forwardVelocity = p->velocity;
  forwardVelocity.Project(forward);
  if (right_key && !left_key) {
    f += (p->max_speed * forward - forwardVelocity) / TAU;
  }
  else if (left_key && !right_key) {
    f += (-p->max_speed * forward - forwardVelocity) / TAU;
  }
  // Drag
  else if (!p->inAir) {
    const double DRAG_COEFFICIENT = 4;
    f += -1*forwardVelocity*DRAG_COEFFICIENT;// * DRAG_COEFFICIENT;
  }
  
  if (p->has_won)
    p->velocity = R3null_vector;
  else
    p->velocity += (f / p->mass) * delta_time;
  
  // transform the player node
  if (!p->is_dead && !p->has_won)
  {
    R3Matrix tform = p->node->transformation;
    tform.Translate(p->velocity * delta_time);
    p->node->transformation = tform;
  }
  
  // set inair to true: it will be set to false if collision with ground detected

  //p->onPlatform = false;
  p->inAir = true;
  CollidePlayer(scene->root);
  
  if (p->onPlatform) {
    p->node->transformation.Translate(p->platform->velocity * delta_time);
  }

  // Camera Shit
  scene->camera.eye = p->Center() - 40 * p->Right() + 10 * p->Up();
  scene->camera.towards = (p->Center()) - scene->camera.eye;
  scene->camera.towards.Normalize();
  scene->camera.right = p->Towards();
  scene->camera.up = scene->camera.right;
  scene->camera.up.Cross(scene->camera.towards);

  R3Box player_box = *p->node->shape->box;
  player_box.Transform(p->node->transformation);

  if (player_box.Min().Y() <= scene->death_y && !p->is_dead)
  {
    KillPlayer(); 
  }

  if (!(key_state['c'] || key_state['C'])) {
    static double angle = 0;
    double MAX_ROTATION = PI/8;
    double new_angle = -1*MAX_ROTATION * (forwardVelocity.Length() / p->max_speed);
    if (forwardVelocity.Dot(p->Towards()) < 0) { new_angle *= -1; }
    angle = angle * 0.995 + new_angle * 0.005;
    R3Line center_line(p->Center(), p->Up());
    scene->camera.Rotate(center_line, angle);
    camera = scene->camera;
  }
}

void UpdateCoin(R3Coin *coin, double delta_time)
{
  coin->t += delta_time;

  R3Matrix tform = R3identity_matrix;
  tform.Translate(coin->position.Vector());
  tform.Rotate(R3_X, PI / 2);
  tform.Rotate(R3_Z, coin->t * 4);
  coin->node->transformation = tform;
}

void UpdateCoins(R3Scene *scene, double delta_time)
{
  // time passed since starting
  for (int i = 0; i < scene->NCoins(); i++)
  {
    R3Coin *coin = scene->Coin(i);
    UpdateCoin(coin, delta_time);
  }
}


void UpdatePlatform(R3Platform *platform, double delta_time) {
  R3Point pos = platform->node->shape->box->Min();
  pos.Transform(platform->node->transformation);
  
  double K = 2; // spring constant
  R3Vector displacement = platform->center - pos;
  platform->velocity += platform->max_speed * displacement * delta_time;
  platform->node->transformation.Translate(platform->velocity * delta_time);
}


void DeleteNodes(R3Node *node) {
  for (vector<R3Node *>::iterator it = node->children.begin(); it != node->children.end();)
  {
    if ((*it)->del)
    {
      node->children.erase(it);
    }
    else
    {
      it++;
    }
  }
  for (vector<R3Node *>::iterator it = node->children.begin(); it != node->children.end(); it++)
  {
    DeleteNodes(*it);
  }
}

void DeleteCoins() {
  for (vector<R3Coin *>::iterator it = scene->coins.begin(); it != scene->coins.end();)
  {
    if ((*it)->del)
      scene->coins.erase(it);
    else
      it++;
  }
}

void DeleteEnemies() {
  for (vector<R3Enemy *>::iterator it = scene->enemies.begin(); it != scene->enemies.end();)
  {
    if ((*it)->del)
      scene->enemies.erase(it);
    else
      it++;
  }
}

void UpdatePlatforms(R3Scene *scene, double delta_time) {
  int numPlatforms = scene->platforms.size();
  for (int i = 0; i < numPlatforms; ++i) {
    R3Platform *cur = scene->platforms[i];
    if (delta_time == 0) {
      cur->velocity = R3null_vector;
    }
    UpdatePlatform(cur, delta_time);
  }
}

void UpdateEnemies(R3Scene *scene, double delta_time) {
  
  int numEnemies = scene->enemies.size();

  for (int i = 0; i < numEnemies; i++) {
    R3Enemy *p = scene->enemies[i];

    if (p == NULL) return; 
    if (delta_time == 0) {
      p->inAir = true;
      p->onPlatform = false;
    }
    // Motion Shit
    // get the forces to move the box
    R3Vector f = R3null_vector;
    f += -9.8 * p->Up() * p->mass;

    if (!p->inAir && p->is_jumping) {
      p->velocity += 10 * p->Up();
      if (scene->player && R3Distance(scene->player->Center(), p->Center()) < 20.0f)
      {
        char path[FILENAME_MAX + 1];
        FilePath(path, "/../sounds/jump.wav");
        ISound *jump_sound = sound_engine->play2D(path, false, false, true);
        double dist = 1.0f;
        if (scene->player) {
          double r = R3Distance(scene->player->Center(), p->Center()) / 10;
          dist = (1.0f / (1.0f + r * r));
        }
        jump_sound->setVolume(dist);
        jump_sound->setIsPaused(false);
      }
      if (p->onPlatform) {
        p->velocity += p->platform->velocity;
      }
    }

    // side to side
    double TAU = .5; // timescale for velocity relaxation
    const R3Vector forward = p->Towards();
    R3Vector forwardVelocity = p->velocity;
    forwardVelocity.Project(forward);

    R3Box enemy_box = *(p->node->shape->box);
    enemy_box.Transform(p->node->transformation);

    R3Box player_box = *(scene->player->node->shape->box); 
    player_box.Transform(scene->player->node->transformation); 

    double direction = enemy_box.XMin() - player_box.XMin();

    //Only change direction if the enemy is following
    if (p->is_following) {
      if (direction < 0) {p->moveLeft = false;}
      else {p->moveLeft = true;}
    }

    if (!p->moveLeft) {
      f += (p->speed * forward - forwardVelocity) / TAU;
    }
    else if (p->moveLeft) {
      f += (-p->speed * forward - forwardVelocity) / TAU;
    }
    // Drag
    else if (!p->inAir) {
      const double DRAG_COEFFICIENT = 2;
      f += -1*forwardVelocity*DRAG_COEFFICIENT; // * DRAG_COEFFICIENT;
    }
    
    p->velocity += (f / p->mass) * delta_time;
    
    // transform the player node
    if (!p->is_dead)
    {
      R3Matrix tform = p->node->transformation;
      tform.Translate(p->velocity * delta_time);
      p->node->transformation = tform;
    }
    
    // set inair to true: it will be set to false if collision with ground detected
    p->inAir = true;
    p->onPlatform = false;
    CollideEnemy(p, scene->root);

    if (p->onPlatform) {
      p->node->transformation.Translate(p->platform->velocity * delta_time);
    }
  }
}


void ClickSidebar(int x, int y) {
  // did we click a button?
  int sidebarLeftX = GLUTwindow_width - scene->sidebar->width;
  if (x < sidebarLeftX ||
      (x > sidebarLeftX && x < sidebarLeftX + scene->sidebar->border) ||
      (x > GLUTwindow_width - scene->sidebar->border)) {
    return;
  }
  // how many buttons down are we?
  unsigned int button = floor((GLUTwindow_height - y) / (scene->sidebar->button_width + scene->sidebar->border));
  if (button >= scene->sidebar->buttons.size()) {
    return;
  }
  // this code is to determine if we're in the border between buttons
  int remainder = (GLUTwindow_height - y) - button * (scene->sidebar->button_width + scene->sidebar->border);
  if (GLUTwindow_height - remainder < scene->sidebar->border) {
    return;
  }
  
  if (button == scene->sidebar->selected_button) {
    scene->sidebar->selected_button = -1;
    *scene->sidebar->buttons[button]->value = 0;
  }
  else {
    if (scene->sidebar->selected_button != -1) {
      *scene->sidebar->buttons[scene->sidebar->selected_button]->value = 0;
    }
    if (move_mode != 0) {
      grabbed = NULL;
      move_mode = 0;
    }
    scene->sidebar->selected_button = button;
    *scene->sidebar->buttons[button]->value = 1;
  }
  
}

void DrawShape(R3Shape *shape)
{
  // Check shape type
  if (shape->type == R3_BOX_SHAPE) shape->box->Draw();
  else if (shape->type == R3_SPHERE_SHAPE) shape->sphere->Draw();
  else if (shape->type == R3_CYLINDER_SHAPE) shape->cylinder->Draw();
  else if (shape->type == R3_CONE_SHAPE) shape->cone->Draw();
  else if (shape->type == R3_MESH_SHAPE) shape->mesh->Draw();
  else if (shape->type == R3_SEGMENT_SHAPE) shape->segment->Draw();
  else if (shape->type == R3_CIRCLE_SHAPE) shape->circle->Draw();
  else fprintf(stderr, "Unrecognized shape type: %d\n", shape->type);
}

void DrawGoal(R3Shape *shape)
{
  glDisable(GL_LIGHTING);
  glColor3d(1 - scene->background[0], 1 - scene->background[1], 1 - scene->background[2]);
  glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
  shape->box->Draw();
  glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
  glEnable(GL_LIGHTING);
}

void DrawPlayer(R3Shape *shape, bool has_won)
{
  static double win_duration = 3.0f;
  if (!has_won)
    shape->box->Draw();
  else
  {
    R3Box *box = shape->box;
    double shrinkage = (GetTime() - scene->player->won_time) / win_duration;
    if (shrinkage >= 1.0f)
      return;
    // Get box corner points 
    R3Point corners[8];
    corners[0] = box->Corner(0, 0, 0);
    corners[1] = box->Corner(0, 0, 1);
    corners[2] = box->Corner(0, 1, 1);
    corners[3] = box->Corner(0, 1, 0);
    corners[4] = box->Corner(1, 0, 0);
    corners[5] = box->Corner(1, 0, 1);
    corners[6] = box->Corner(1, 1, 1);
    corners[7] = box->Corner(1, 1, 0);

    for (int i = 0; i < 8; i++)
    {
      R3Vector v = box->Centroid() - corners[i];
      v.Normalize();
      v *= shrinkage;
      corners[i].Translate(v);
    }

    // Get normals, texture coordinates, etc.
    static GLdouble normals[6][3] = {
      { -1.0, 0.0, 0.0 },
      { 1.0, 0.0, 0.0 },
      { 0.0, -1.0, 0.0 },
      { 0.0, 1.0, 0.0 },
      { 0.0, 0.0, -1.0 },
      { 0.0, 0.0, 1.0 }
    };
    static GLdouble texcoords[4][2] = {
      { 0.0, 0.0 },
      { 1.0, 0.0 },
      { 1.0, 1.0 },
      { 0.0, 1.0 }
    };
    static int surface_paths[6][4] = {
      { 3, 0, 1, 2 },
      { 4, 7, 6, 5 },
      { 0, 4, 5, 1 },
      { 7, 3, 2, 6 },
      { 3, 7, 4, 0 },
      { 1, 5, 6, 2 }
    };

    // Draw box
    glBegin(GL_QUADS);
    for (int i = 0; i < 6; i++) {
      glNormal3d(normals[i][0], normals[i][1], normals[i][2]);
      for (int j = 0; j < 4; j++) {
        const R3Point& p = corners[surface_paths[i][j]];
        glTexCoord2d(texcoords[j][0], texcoords[j][1]);
        glVertex3d(p[0], p[1], p[2]);
      }
    }
    glEnd();    
  }
}



void LoadMatrix(R3Matrix *matrix)
{
  // Multiply matrix by top of stack
  // Take transpose of matrix because OpenGL represents vectors with 
  // column-vectors and R3 represents them with row-vectors
  R3Matrix m = matrix->Transpose();
  glMultMatrixd((double *) &m);
}



void LoadMaterial(R3Material *material) 
{
  GLfloat c[4];

  // Check if same as current
  static R3Material *current_material = NULL;
  if (material == current_material) return;
  current_material = material;

  // Compute "opacity"
  double opacity = 1 - material->kt.Luminance();

  // Load ambient
  c[0] = material->ka[0];
  c[1] = material->ka[1];
  c[2] = material->ka[2];
  c[3] = opacity;
  glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, c);

  // Load diffuse
  c[0] = material->kd[0];
  c[1] = material->kd[1];
  c[2] = material->kd[2];
  c[3] = opacity;
  glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, c);

  // Load specular
  c[0] = material->ks[0];
  c[1] = material->ks[1];
  c[2] = material->ks[2];
  c[3] = opacity;
  glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, c);

  // Load emission
  c[0] = material->emission.Red();
  c[1] = material->emission.Green();
  c[2] = material->emission.Blue();
  c[3] = opacity;
  glMaterialfv(GL_FRONT_AND_BACK, GL_EMISSION, c);

  // Load shininess
  c[0] = material->shininess;
  glMaterialf(GL_FRONT_AND_BACK, GL_SHININESS, c[0]);

  // Load texture
  if (material->texture) {
    if (material->texture_index <= 0) {
      // Create texture in OpenGL
      GLuint texture_index;
      glGenTextures(1, &texture_index);
      material->texture_index = (int) texture_index;
      glBindTexture(GL_TEXTURE_2D, material->texture_index); 
      R2Image *image = material->texture;
      int npixels = image->NPixels();
      R2Pixel *pixels = image->Pixels();
      GLfloat *buffer = new GLfloat [ 4 * npixels ];
      R2Pixel *pixelsp = pixels;
      GLfloat *bufferp = buffer;
      for (int j = 0; j < npixels; j++) { 
        *(bufferp++) = pixelsp->Red();
        *(bufferp++) = pixelsp->Green();
        *(bufferp++) = pixelsp->Blue();
        *(bufferp++) = pixelsp->Alpha();
        pixelsp++;
      }
      glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MIN_FILTER,GL_LINEAR);
      glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MAG_FILTER,GL_LINEAR);
      glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_WRAP_S,GL_REPEAT);
      glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_WRAP_T,GL_REPEAT);
      glTexEnvi(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);
      glTexImage2D(GL_TEXTURE_2D, 0, 4, image->Width(), image->Height(), 0, GL_RGBA, GL_FLOAT, buffer);
      delete [] buffer;
    }

    // Select texture
    glBindTexture(GL_TEXTURE_2D, material->texture_index);
    glEnable(GL_TEXTURE_2D);
  }
  else {
    glDisable(GL_TEXTURE_2D);
  }

  // Enable blending for transparent surfaces
  if (opacity < 1) {
    glDepthMask(false);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glEnable(GL_BLEND);
  }
  else {
    glDisable(GL_BLEND);
    glBlendFunc(GL_ONE, GL_ZERO);
    glDepthMask(true);
  }
}



void LoadCamera(R3Camera *camera)
{
  // Set projection transformation
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  gluPerspective(2*180.0*camera->yfov/M_PI, (GLdouble) GLUTwindow_width /(GLdouble) GLUTwindow_height, 0.01, 10000);

  // Set camera transformation
  R3Vector t = -(camera->towards);
  R3Vector& u = camera->up;
  R3Vector& r = camera->right;
  GLdouble camera_matrix[16] = { r[0], u[0], t[0], 0, r[1], u[1], t[1], 0, r[2], u[2], t[2], 0, 0, 0, 0, 1 };
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
  glMultMatrixd(camera_matrix);
  glTranslated(-(camera->eye[0]), -(camera->eye[1]), -(camera->eye[2]));
}



void LoadLights(R3Scene *scene)
{
  GLfloat buffer[4];

  // Load ambient light
  static GLfloat ambient[4];
  ambient[0] = scene->ambient[0];
  ambient[1] = scene->ambient[1];
  ambient[2] = scene->ambient[2];
  ambient[3] = 1;
  glLightModelfv(GL_LIGHT_MODEL_AMBIENT, ambient);

  // Load scene lights
  for (int i = 0; i < (int) scene->lights.size(); i++) {
    R3Light *light = scene->lights[i];
    int index = GL_LIGHT0 + i;

    // Temporarily disable light
    glDisable(index);

    // Load color
    buffer[0] = light->color[0];
    buffer[1] = light->color[1];
    buffer[2] = light->color[2];
    buffer[3] = 1.0;
    glLightfv(index, GL_DIFFUSE, buffer);
    glLightfv(index, GL_SPECULAR, buffer);

    // Load attenuation with distance
    buffer[0] = light->constant_attenuation;
    buffer[1] = light->linear_attenuation;
    buffer[2] = light->quadratic_attenuation;
    glLightf(index, GL_CONSTANT_ATTENUATION, buffer[0]);
    glLightf(index, GL_LINEAR_ATTENUATION, buffer[1]);
    glLightf(index, GL_QUADRATIC_ATTENUATION, buffer[2]);

    // Load spot light behavior
    buffer[0] = 180.0 * light->angle_cutoff / M_PI;
    buffer[1] = light->angle_attenuation;
    glLightf(index, GL_SPOT_CUTOFF, buffer[0]);
    glLightf(index, GL_SPOT_EXPONENT, buffer[1]);

    // Load positions/directions
    if (light->type == R3_DIRECTIONAL_LIGHT) {
      // Load direction
      buffer[0] = -(light->direction.X());
      buffer[1] = -(light->direction.Y());
      buffer[2] = -(light->direction.Z());
      buffer[3] = 0.0;
      glLightfv(index, GL_POSITION, buffer);
    }
    else if (light->type == R3_POINT_LIGHT) {
      // Load position
      buffer[0] = light->position.X();
      buffer[1] = light->position.Y();
      buffer[2] = light->position.Z();
      buffer[3] = 1.0;
      glLightfv(index, GL_POSITION, buffer);
    }
    else if (light->type == R3_SPOT_LIGHT) {
      // Load position
      buffer[0] = light->position.X();
      buffer[1] = light->position.Y();
      buffer[2] = light->position.Z();
      buffer[3] = 1.0;
      glLightfv(index, GL_POSITION, buffer);

      // Load direction
      buffer[0] = light->direction.X();
      buffer[1] = light->direction.Y();
      buffer[2] = light->direction.Z();
      buffer[3] = 1.0;
      glLightfv(index, GL_SPOT_DIRECTION, buffer);
    }
    else if (light->type == R3_AREA_LIGHT) {
      // Load position
      buffer[0] = light->position.X();
      buffer[1] = light->position.Y();
      buffer[2] = light->position.Z();
      buffer[3] = 1.0;
      glLightfv(index, GL_POSITION, buffer);

      // Load direction
      buffer[0] = light->direction.X();
      buffer[1] = light->direction.Y();
      buffer[2] = light->direction.Z();
      buffer[3] = 1.0;
      glLightfv(index, GL_SPOT_DIRECTION, buffer);
    }
    else {
      fprintf(stderr, "Unrecognized light type: %d\n", light->type);
      return;
    }

    // Enable light
    glEnable(index);
  }
}



void DrawNode(R3Scene *scene, R3Node *node)
{
  // Push transformation onto stack
  glPushMatrix();
  LoadMatrix(&node->transformation);

  // Load material
  if (node->material) LoadMaterial(node->material);

  // Draw shape
  if (node->is_goal && node->shape && node->is_goal) DrawGoal(node->shape);
  else if (node->is_player && node->shape && scene->player && node->is_visible) DrawPlayer(node->shape, scene->player->has_won);
  else if (node->is_visible && node->shape) DrawShape(node->shape);

  // Draw children nodes
  for (int i = 0; i < (int) node->children.size(); i++) 
    DrawNode(scene, node->children[i]);

  // Restore previous transformation
  glPopMatrix();

  // Show bounding box
  if (show_bboxes) {
    GLboolean lighting = glIsEnabled(GL_LIGHTING);
    glDisable(GL_LIGHTING);
    node->bbox.Outline();
    if (lighting) glEnable(GL_LIGHTING);
  }
}

void DrawLights(R3Scene *scene)
{
  // Check if should draw lights
  if (!show_lights) return;

  // Setup
  GLboolean lighting = glIsEnabled(GL_LIGHTING);
  glDisable(GL_LIGHTING);

  // Draw all lights
  double radius = scene->bbox.DiagonalRadius();
  for (int i = 0; i < scene->NLights(); i++) {
    R3Light *light = scene->Light(i);
    glColor3d(light->color[0], light->color[1], light->color[2]);
    if (light->type == R3_DIRECTIONAL_LIGHT) {
      // Draw direction vector
      glLineWidth(5);
      glBegin(GL_LINES);
      R3Point centroid = scene->bbox.Centroid();
      R3Vector vector = radius * light->direction;
      glVertex3d(centroid[0] - vector[0], centroid[1] - vector[1], centroid[2] - vector[2]);
      glVertex3d(centroid[0] - 1.25*vector[0], centroid[1] - 1.25*vector[1], centroid[2] - 1.25*vector[2]);
      glEnd();
      glLineWidth(1);
    }
    else if (light->type == R3_POINT_LIGHT) {
      // Draw sphere at point light position
      R3Sphere(light->position, 0.1 * radius).Draw();
    }
    else if (light->type == R3_SPOT_LIGHT) {
      // Draw sphere at point light position and line indicating direction
      R3Sphere(light->position, 0.1 * radius).Draw();
  
      // Draw direction vector
      glLineWidth(5);
      glBegin(GL_LINES);
      R3Vector vector = radius * light->direction;
      glVertex3d(light->position[0], light->position[1], light->position[2]);
      glVertex3d(light->position[0] + 0.25*vector[0], light->position[1] + 0.25*vector[1], light->position[2] + 0.25*vector[2]);
      glEnd();
      glLineWidth(1);
    }
    else if (light->type == R3_AREA_LIGHT) {
      // Draw circular area
      R3Vector v1, v2;
      double r = light->radius;
      R3Point p = light->position;
      int dim = light->direction.MinDimension();
      if (dim == 0) { v1 = light->direction % R3posx_vector; v1.Normalize(); v2 = light->direction % v1; }
      else if (dim == 1) { v1 = light->direction % R3posy_vector; v1.Normalize(); v2 = light->direction % v1; }
      else { v1 = light->direction % R3posz_vector; v1.Normalize(); v2 = light->direction % v1; }
      glBegin(GL_POLYGON);
      glVertex3d(p[0] +  1.00*r*v1[0] +  0.00*r*v2[0], p[1] +  1.00*r*v1[1] +  0.00*r*v2[1], p[2] +  1.00*r*v1[2] +  0.00*r*v2[2]);
      glVertex3d(p[0] +  0.71*r*v1[0] +  0.71*r*v2[0], p[1] +  0.71*r*v1[1] +  0.71*r*v2[1], p[2] +  0.71*r*v1[2] +  0.71*r*v2[2]);
      glVertex3d(p[0] +  0.00*r*v1[0] +  1.00*r*v2[0], p[1] +  0.00*r*v1[1] +  1.00*r*v2[1], p[2] +  0.00*r*v1[2] +  1.00*r*v2[2]);
      glVertex3d(p[0] + -0.71*r*v1[0] +  0.71*r*v2[0], p[1] + -0.71*r*v1[1] +  0.71*r*v2[1], p[2] + -0.71*r*v1[2] +  0.71*r*v2[2]);
      glVertex3d(p[0] + -1.00*r*v1[0] +  0.00*r*v2[0], p[1] + -1.00*r*v1[1] +  0.00*r*v2[1], p[2] + -1.00*r*v1[2] +  0.00*r*v2[2]);
      glVertex3d(p[0] + -0.71*r*v1[0] + -0.71*r*v2[0], p[1] + -0.71*r*v1[1] + -0.71*r*v2[1], p[2] + -0.71*r*v1[2] + -0.71*r*v2[2]);
      glVertex3d(p[0] +  0.00*r*v1[0] + -1.00*r*v2[0], p[1] +  0.00*r*v1[1] + -1.00*r*v2[1], p[2] +  0.00*r*v1[2] + -1.00*r*v2[2]);
      glVertex3d(p[0] +  0.71*r*v1[0] + -0.71*r*v2[0], p[1] +  0.71*r*v1[1] + -0.71*r*v2[1], p[2] +  0.71*r*v1[2] + -0.71*r*v2[2]);
      glEnd();
    }
    else {
      fprintf(stderr, "Unrecognized light type: %d\n", light->type);
      return;
    }
  }

  // Clean up
  if (lighting) glEnable(GL_LIGHTING);
}



void DrawCamera(R3Scene *scene)
{
  // Check if should draw lights
  if (!show_camera) return;

  // Setup
  GLboolean lighting = glIsEnabled(GL_LIGHTING);
  glDisable(GL_LIGHTING);
  glColor3d(1.0, 1.0, 1.0);
  glLineWidth(5);

  // Draw view frustum
  R3Camera& c = scene->camera;
  double radius = scene->bbox.DiagonalRadius();
  R3Point org = c.eye + c.towards * radius;
  R3Vector dx = c.right * radius * tan(c.xfov);
  R3Vector dy = c.up * radius * tan(c.yfov);
  R3Point ur = org + dx + dy;
  R3Point lr = org + dx - dy;
  R3Point ul = org - dx + dy;
  R3Point ll = org - dx - dy;
  glBegin(GL_LINE_LOOP);
  glVertex3d(ur[0], ur[1], ur[2]);
  glVertex3d(ul[0], ul[1], ul[2]);
  glVertex3d(ll[0], ll[1], ll[2]);
  glVertex3d(lr[0], lr[1], lr[2]);
  glVertex3d(ur[0], ur[1], ur[2]);
  glVertex3d(c.eye[0], c.eye[1], c.eye[2]);
  glVertex3d(lr[0], lr[1], lr[2]);
  glVertex3d(ll[0], ll[1], ll[2]);
  glVertex3d(c.eye[0], c.eye[1], c.eye[2]);
  glVertex3d(ul[0], ul[1], ul[2]);
  glEnd();

  // Clean up
  glLineWidth(1);
  if (lighting) glEnable(GL_LIGHTING);
}



void DrawScene(R3Scene *scene) 
{
  // Draw nodes recursively
  DrawNode(scene, scene->root);
}


// void DrawParticles(R3Scene *scene)
// {
//   // Get current time (in seconds) since start of execution
//   double current_time = GetTime();
//   static double previous_time = 0;


//   static double time_lost_taking_videos = 0; // for switching back and forth
// 					     // between recording and not
// 					     // recording smoothly

//   // program just started up?
//   if (previous_time == 0) previous_time = current_time;

//   // time passed since starting
//   double delta_time = current_time - previous_time;


//   if (save_video) { // in video mode, the time that passes only depends on the frame rate ...
//     delta_time = VIDEO_FRAME_DELAY;    
//     // ... but we need to keep track how much time we gained and lost so that we can arbitrarily switch back and forth ...
//     time_lost_taking_videos += (current_time - previous_time) - VIDEO_FRAME_DELAY;
//   } else { // real time simulation
//     delta_time = current_time - previous_time;
//   }

//   // Update particles
//   UpdateParticles(scene, current_time - time_lost_taking_videos, delta_time, integration_type);

//   // Generate new particles
//   GenerateParticles(scene, current_time - time_lost_taking_videos, delta_time);

//   // Render particles
//   if (show_particles) RenderParticles(scene, current_time - time_lost_taking_videos, delta_time);

//   // Remember previous time
//   previous_time = current_time;
// }


void DrawParticleSources(R3Scene *scene)
{
  // Check if should draw particle sources
  if (!show_particle_sources_and_sinks) return;

  // Setup
  GLboolean lighting = glIsEnabled(GL_LIGHTING);
  glEnable(GL_LIGHTING);

  // Define source material
  static R3Material source_material;
  if (source_material.id != 33) {
    source_material.ka.Reset(0.2,0.2,0.2,1);
    source_material.kd.Reset(0,1,0,1);
    source_material.ks.Reset(0,1,0,1);
    source_material.kt.Reset(0,0,0,1);
    source_material.emission.Reset(0,0,0,1);
    source_material.shininess = 1;
    source_material.indexofrefraction = 1;
    source_material.texture = NULL;
    source_material.texture_index = -1;
    source_material.id = 33;
  }

  // Draw all particle sources
  glEnable(GL_LIGHTING);
  LoadMaterial(&source_material);
  for (int i = 0; i < scene->NParticleSources(); i++) {
    R3ParticleSource *source = scene->ParticleSource(i);
    DrawShape(source->shape);
  }

  // Clean up
  if (!lighting) glDisable(GL_LIGHTING);
}



void DrawParticleSinks(R3Scene *scene)
{
  // Check if should draw particle sinks
  if (!show_particle_sources_and_sinks) return;

  // Setup
  GLboolean lighting = glIsEnabled(GL_LIGHTING);
  glEnable(GL_LIGHTING);

  // Define sink material
  static R3Material sink_material;
  if (sink_material.id != 33) {
    sink_material.ka.Reset(0.2,0.2,0.2,1);
    sink_material.kd.Reset(1,0,0,1);
    sink_material.ks.Reset(1,0,0,1);
    sink_material.kt.Reset(0,0,0,1);
    sink_material.emission.Reset(0,0,0,1);
    sink_material.shininess = 1;
    sink_material.indexofrefraction = 1;
    sink_material.texture = NULL;
    sink_material.texture_index = -1;
    sink_material.id = 33;
  }

  // Draw all particle sinks
  glEnable(GL_LIGHTING);
  LoadMaterial(&sink_material);
  for (int i = 0; i < scene->NParticleSinks(); i++) {
    R3ParticleSink *sink = scene->ParticleSink(i);
    DrawShape(sink->shape);
  }

  // Clean up
  if (!lighting) glDisable(GL_LIGHTING);
}



void DrawParticleSprings(R3Scene *scene)
{
  // Check if should draw particle springs
  if (!show_particle_springs) return;

  // Setup
  GLboolean lighting = glIsEnabled(GL_LIGHTING);
  glDisable(GL_LIGHTING);

  // Draw all particle sources
  glColor3d(0.5, 0.5, 0.5);
  glBegin(GL_LINES);
  for (unsigned int i = 0; i < scene->particle_springs.size(); i++) {
    R3ParticleSpring *spring = scene->particle_springs[i];
    const R3Point& p0 = spring->particles[0]->position;
    const R3Point& p1 = spring->particles[1]->position;
    glVertex3d(p0[0], p0[1], p0[2]);
    glVertex3d(p1[0], p1[1], p1[2]);
  }
  glEnd();

  // Clean up
  if (lighting) glEnable(GL_LIGHTING);
}

// draw heads up display
void DrawHUD()
{
  // disable lighting and enable orthogonal projection
  glDisable(GL_LIGHTING);
  glMatrixMode(GL_PROJECTION);
  glPushMatrix();
  glLoadIdentity();
  glOrtho(0.0, GLUTwindow_width, GLUTwindow_height, 0.0, -1.0, 10.0);
  glMatrixMode(GL_MODELVIEW);
  glPushMatrix();
  glLoadIdentity();
  glDisable(GL_CULL_FACE);
  glDisable(GL_LIGHTING); 

  glClear(GL_DEPTH_BUFFER_BIT);

  // Level editor stuff
  if (level_editor)
  {
    DrawSidebar(scene);
  }

  // Draw coins as squares in top left
  float spacing = 15.0;
  float size = 50.0;

  //Cancels if no player
  if (scene->player == NULL) return; 

  for (int i = 0; i < scene->player->n_coins; i++)
  {
    float xmin = spacing * (i + 1) + size * i;
    float xmax = spacing * (i + 1) + size * (i + 1);
    float ymin = spacing;
    float ymax = spacing + size;

    // glColor4f(0.0f, 0.0f, 0.0, 0);
    // glBindTexture(GL_TEXTURE_2D, scene->coins[0]->node->material->texture_index); 
    glBindTexture(GL_TEXTURE_2D, 0); 
    glBegin(GL_QUADS);
      glColor3f(.8f, .8f, 0.0);
      // glTexCoord2f(0.0, 0.0); 
      glVertex2f(xmin, ymin);
      // glTexCoord2f(1.0, 0.0); 
      glVertex2f(xmax, ymin);
      // glTexCoord2f(1.0, 1.0); 
      glVertex2f(xmax, ymax);
      // glTexCoord2f(0.0, 1.0); 
      glVertex2f(xmin, ymax);
    glEnd();
  }

  // Making sure we can render 3d again
  glMatrixMode(GL_MODELVIEW);
  glPopMatrix();
  glMatrixMode(GL_PROJECTION);
  glPopMatrix();
}

void DrawSkybox(R3Scene *scene) {
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP); 
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP);

   // Store the current matrix
   glPushMatrix();

   // Reset and transform the matrix.
   glLoadIdentity();
   gluLookAt(
       0,0,0,
       scene->camera.towards.X(),scene->camera.towards.Y(),scene->camera.towards.Z() - 5,
       0,1,0);

   // Enable/Disable features
   glPushAttrib(GL_ENABLE_BIT);
   glEnable(GL_TEXTURE_2D);
   glDisable(GL_DEPTH_TEST);
   glDisable(GL_LIGHTING);
   glDisable(GL_BLEND);

   glColor4f(1,1,1,1);
   // Render the front quad
   glBindTexture(GL_TEXTURE_2D, skyboxMaterials[0]->texture_index);
   glBegin(GL_QUADS);
       glTexCoord2f(0, 0); glVertex3f( -0.5f, -0.5f, -0.5f );
       glTexCoord2f(1, 0); glVertex3f( -0.5f, 0.5f, -0.5f );
       glTexCoord2f(1, 1); glVertex3f( 0.5f,  0.5f, -0.5f );
       glTexCoord2f(0, 1); glVertex3f(  0.5f, -0.5f, -0.5f );
   glEnd();
   glBindTexture(GL_TEXTURE_2D, 0);
   // Restore enable bits wand matrix

   glPopAttrib();
   glPopMatrix();

}

////////////////////////////////////////////////////////////
// GLUT USER INTERFACE CODE
////////////////////////////////////////////////////////////

void GLUTMainLoop(void)
{
  // Run main loop -- never returns 
  glutMainLoop();
}

void GLUTDrawText(const R3Point& p, const char *s)
{
  // Draw text string s and position p
  glRasterPos3d(p[0], p[1], p[2]);
  while (*s) glutBitmapCharacter(GLUT_BITMAP_HELVETICA_12, *(s++));
}

void GLUTSaveImage(const char *filename)
{ 
  // Create image
  R2Image image(GLUTwindow_width, GLUTwindow_height);

  // Read screen into buffer
  GLfloat *pixels = new GLfloat [ 3 * GLUTwindow_width * GLUTwindow_height ];
  glReadPixels(0, 0, GLUTwindow_width, GLUTwindow_height, GL_RGB, GL_FLOAT, pixels);

  // Load pixels from frame buffer
  GLfloat *pixelsp = pixels;
  for (int j = 0; j < GLUTwindow_height; j++) {
    for (int i = 0; i < GLUTwindow_width; i++) {
      double r = (double) *(pixelsp++);
      double g = (double) *(pixelsp++);
      double b = (double) *(pixelsp++);
      R2Pixel pixel(r, g, b, 1);
      image.SetPixel(i, j, pixel);
    }
  }

  // Write image to file
  image.Write(filename);

  // Delete buffer
  delete [] pixels;
}



void GLUTStop(void)
{
  // Destroy window 
  glutDestroyWindow(GLUTwindow);

  // Delete scene
  delete scene;

  // Exit
  exit(0);
}



void GLUTIdle(void)
{
  // Set current window
  if ( glutGetWindow() != GLUTwindow ) 
    glutSetWindow(GLUTwindow);  

  // Redraw
  glutPostRedisplay();
}



void GLUTResize(int w, int h)
{
  // Resize window
  glViewport(0, 0, w, h);

  // Resize camera vertical field of view to match aspect ratio of viewport
  camera.yfov = atan(tan(camera.xfov) * (double) h/ (double) w); 

  // Remember window size 
  GLUTwindow_width = w;
  GLUTwindow_height = h;

  // Redraw
  glutPostRedisplay();
}


void GLUTRedraw(void)
{
  // Initialize OpenGL drawing modes
  glEnable(GL_LIGHTING);
  glDisable(GL_BLEND);
  glBlendFunc(GL_ONE, GL_ZERO);
  glDepthMask(true);

  // Clear window 
  R3Rgb background = scene->background;
  glClearColor(background[0], background[1], background[2], background[3]);
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  
  double current_time = GetTime();

  while (previous_time < current_time) {
    double delta_time;
    if (previous_time == 0) {
      delta_time = 0;
      previous_time = current_time;
      if (level_editor) camera = minimap_cam;
    }
    else {
      delta_time = TIME_STEP;
      previous_time += TIME_STEP;
    }
    
    // Update Player
    if (level_editor != 1 && !scene->player->is_dead) {
      UpdatePlayer(scene, delta_time);
    }

    UpdatePlatforms(scene, delta_time);
    UpdateEnemies(scene, delta_time);

    // Update Coins
    UpdateCoins(scene, delta_time);
    
    // delete objects
    DeleteNodes(scene->root);
    DeleteCoins();
    DeleteEnemies();

    GenerateParticles(scene, current_time, delta_time);
    UpdateParticles(scene, current_time, delta_time, MIDPOINT_INTEGRATION);
  }


  // Load camera
  LoadCamera(&camera);

  // Load scene lights
  LoadLights(scene);

  // Draw scene camera
  DrawCamera(scene);

  // Draw scene lights
  DrawLights(scene);

  // Draw particles
  RenderParticles(scene);

  // Draw particle sources 
  DrawParticleSources(scene);

  // Draw particle sinks 
  DrawParticleSinks(scene);

  // Draw particle springs
  DrawParticleSprings(scene);

  DrawSkybox(scene); 
  

  // Draw scene surfaces
  if (show_faces) {
    glEnable(GL_LIGHTING);
    DrawScene(scene);
  }

  // Draw scene edges
  if (show_edges) {
    glDisable(GL_LIGHTING);
    glColor3d(1 - background[0], 1 - background[1], 1 - background[2]);
    glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
    DrawScene(scene);
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
  }
  
  if (minimap) {
    // Resize window
    glViewport(0, 0, GLUTwindow_width/3,  GLUTwindow_height/3);
    glEnable(GL_LIGHTING);
    LoadCamera(&minimap_cam);
    DrawScene(scene);
    LoadCamera(&camera);
    glViewport(0, 0, GLUTwindow_width, GLUTwindow_height);
  }
  
  DrawHUD();

  // Save image
  if (save_image) {
    char image_name[256];
    static int image_number = 1;
    for (;;) {
      sprintf(image_name, "image%d.jpg", image_number++);
      FILE *fp = fopen(image_name, "r");
      if (!fp) break; 
      else fclose(fp);
    }
    GLUTSaveImage(image_name);
    printf("Saved %s\n", image_name);
    save_image = 0;
  }

  // Save video
  if (save_video) {
    char frame_name[512];
    static int next_frame = 0;
    static int num_frames_recorded = 0;
    for (;;) {
      sprintf(frame_name, "%sframe%04d.jpg", video_prefix, next_frame++);
      FILE *fp = fopen(frame_name, "r");
      if (!fp) break; 
      else fclose(fp);
    }
    GLUTSaveImage(frame_name);
    if (next_frame % 100 == 1) {
      printf("Saved %s\n", frame_name);
    }
    if (num_frames_to_record == ++num_frames_recorded) {
      save_video = 0;
      printf("Recorded %d frames, stopping as instructed.\n", num_frames_recorded);
      quit = 1;
    }
  }

  // Quit here so that can save image before exit
  if (quit) {
    if (output_image_name) GLUTSaveImage(output_image_name);
    GLUTStop();
  }

  if (scene->player && scene->player->has_won)
  {
    if (GetTime() - scene->player->won_time > 10.0f && strlen(scene->next_level))
    {
      LoadLevel(scene->next_level);
    }
  }

  // Swap buffers 
  glutSwapBuffers();
}    



void GLUTMotion(int x, int y)
{
  // Invert y coordinate
  y = GLUTwindow_height - y;
  
  // Compute mouse movement
  int dx = x - GLUTmouse[0];
  int dy = y - GLUTmouse[1];
  
  // Process mouse motion event
  if ((dx != 0) || (dy != 0)) {
    if (camera_mode) {
        R3Point scene_center = scene->bbox.Centroid();
        if ((GLUTbutton[0] && (GLUTmodifiers & GLUT_ACTIVE_SHIFT)) || GLUTbutton[1]) {
          // Scale world 
          double factor = (double) dx / (double) GLUTwindow_width;
          factor += (double) dy / (double) GLUTwindow_height;
          factor = exp(2.0 * factor);
          factor = (factor - 1.0) / factor;
          R3Vector translation = (scene_center - camera.eye) * factor;
          camera.eye += translation;
          glutPostRedisplay();
        }
        else if (GLUTbutton[0] && (GLUTmodifiers & GLUT_ACTIVE_CTRL)) {
          // Rotate world
          double vx = (double) dx / (double) GLUTwindow_width;
          double vy = (double) dy / (double) GLUTwindow_height;
          double theta = 4.0 * (fabs(vx) + fabs(vy));
          R3Vector vector = (camera.right * vx) + (camera.up * vy);
          R3Vector rotation_axis = camera.towards % vector;
          rotation_axis.Normalize();
          camera.eye.Rotate(R3Line(scene_center, rotation_axis), theta);
          camera.towards.Rotate(rotation_axis, theta);
          camera.up.Rotate(rotation_axis, theta);
          camera.right = camera.towards % camera.up;
          camera.up = camera.right % camera.towards;
          camera.towards.Normalize();
          camera.up.Normalize();
          camera.right.Normalize();
          glutPostRedisplay();
        }
        else if (GLUTbutton[0]) {
          // Translate world
          double length = R3Distance(scene_center, camera.eye) * tan(camera.yfov);
          double vx = length * (double) dx / (double) GLUTwindow_width;
          double vy = length * (double) dy / (double) GLUTwindow_height;
          R3Vector translation = -((camera.right * vx) + (camera.up * vy));
          camera.eye += translation;
          glutPostRedisplay();
        }
      }
    else if (move_mode && GLUTbutton[0]) {
      R3Point scene_center = scene->bbox.Centroid();
      double length = R3Distance(scene_center, camera.eye) * tan(camera.yfov);
      double vx = length * (double) dx / (double) GLUTwindow_width;
      double vy = length * (double) dy / (double) GLUTwindow_height;
      R3Vector translation = ((camera.right * vx) + (camera.up * vy));
      // snap the shit to a grid
      grabbed->shape->box->Translate(translation);
      grabbed->bbox = *grabbed->shape->box;
    }
  }
  
  // Remember mouse position
  GLUTmouse[0] = x;
  GLUTmouse[1] = y;
}



void GLUTMouse(int button, int state, int x, int y)
{
  // Invert y coordinate
  y = GLUTwindow_height - y;
  int w = glutGet(GLUT_WINDOW_WIDTH);
  int h = glutGet(GLUT_WINDOW_HEIGHT);
  camera.xfov = atan(tan(camera.yfov) * (double) w/ (double) h);
  
  // Process mouse button event
  if (state == GLUT_DOWN) {
    if (button == GLUT_LEFT_BUTTON && level_editor) {
      if ((x < w - scene->sidebar->width) && (blocks_mode == 1)) {
        R3Ray ray = RayThoughPixel(camera, x, y, w, h);
        R3Point click_location = RayPlaneIntersection(scene->movement_plane, ray);
        CreateShape(R3_BOX_SHAPE, scene, click_location);
      }
      else if ((x < w - scene->sidebar->width) && (grab_mode == 1)) {
        R3Ray ray = RayThoughPixel(camera, x, y, w, h);
        R3Intersection intersection = ComputeIntersection(scene, scene->root, ray, 1000000000);
        if (intersection.hit) {
          grab_mode = 0;
          move_mode = 1;
          grabbed = intersection.node;
          scene->sidebar->selected_button = -1;
        }
      }
      else if ((x < w - scene->sidebar->width) && (coins_mode == 1)) {
        R3Ray ray = RayThoughPixel(camera, x, y, w, h);
        R3Point click_location = RayPlaneIntersection(scene->movement_plane, ray);
        CreateShape(R3_COIN_SHAPE, scene, click_location);
      }
      else if (x >  w - scene->sidebar->width) {
        ClickSidebar(x, y);
      }
    }
    else if (button == GLUT_MIDDLE_BUTTON) {
    }
    else if (button == GLUT_RIGHT_BUTTON) {
    }
  }

  // Remember button state 
  int b = (button == GLUT_LEFT_BUTTON) ? 0 : ((button == GLUT_MIDDLE_BUTTON) ? 1 : 2);
  GLUTbutton[b] = (state == GLUT_DOWN) ? 1 : 0;

  // Remember modifiers 
  GLUTmodifiers = glutGetModifiers();

   // Remember mouse position 
  GLUTmouse[0] = x;
  GLUTmouse[1] = y;

  // Redraw
  glutPostRedisplay();
}



// void GLUTSpecial(int key, int x, int y)
// {
//   // Invert y coordinate
//   y = GLUTwindow_height - y;

//   // Process keyboard button event 
//   switch (key) {
//     case GLUT_KEY_F1:
//       save_image = 1;
//       break;
//     case GLUT_KEY_F2:
//       save_video = save_video ^ 1;
//       break;
//   }

//   // Remember mouse position 
//   GLUTmouse[0] = x;
//   GLUTmouse[1] = y;

//   // Remember modifiers 
//   GLUTmodifiers = glutGetModifiers();

//   // Redraw
//   glutPostRedisplay();
// }

// void GLUTSpecialUp(int key, int x, int y)
// {
//   // Invert y coordinate
//   y = GLUTwindow_height - y;

//   // Process keyboard button event
//   switch (key) {
//     case GLUT_KEY_RIGHT:
//       right = false;
//       break;
//     case GLUT_KEY_LEFT:
//       left = false;
//       break;
//     case GLUT_KEY_UP:
//       up = false;
//       break;
//   }

//   // Remember mouse position
//   GLUTmouse[0] = x;
//   GLUTmouse[1] = y;

//   // Remember modifiers
//   GLUTmodifiers = glutGetModifiers();

//   // Redraw
//   glutPostRedisplay();
// }

void GLUTKeyboard(unsigned char key, int x, int y)
{
  // Invert y coordinate
  y = GLUTwindow_height - y;

  key_state[key] = true;

  // Process keyboard button event 
  switch (key) {
  // case 'L':
  // case 'l':
  //   show_lights = !show_lights;
  //   break;
  case 'c':
  case 'C':
      camera = minimap_cam;
  case 'M':
  case 'm':
      minimap = (minimap == 1) ? 0 : 1;
      break;
  case 'P':
  case 'p':
    scene->Write("DoesNothingRightNow", scene->root); 
    break;

  case 'R':
  case 'r':
      if (scene->player && scene->player->is_dead) {
        LoadLevel(input_scene_name);
        paused = 1;
      }
    break;

  // case 'P':
  // case 'p':
  //   show_particles = !show_particles;
  //   break;

  // case 'R':
  // case 'r':
  //   show_particle_springs = !show_particle_springs;
  //   break;

  // case 'S':
  // case 's':
  //   show_particle_sources_and_sinks = !show_particle_sources_and_sinks;
  //   break;

  // case 'Q':
  // case 'q':
  case 27: // ESCAPE
    quit = 1;
    break;
  // case ' ': {
  //   printf("camera %g %g %g  %g %g %g  %g %g %g  %g  %g %g \n",
  //          camera.eye[0], camera.eye[1], camera.eye[2], 
  //          camera.towards[0], camera.towards[1], camera.towards[2], 
  //          camera.up[0], camera.up[1], camera.up[2], 
  //          camera.xfov, camera.neardist, camera.fardist); 
  //   break; }
  }

  // Remember mouse position 
  GLUTmouse[0] = x;
  GLUTmouse[1] = y;

  // Remember modifiers 
  GLUTmodifiers = glutGetModifiers();

  // Redraw
  glutPostRedisplay();
}

void GLUTKeyboardUp(unsigned char key, int x, int y)
{
  key_state[key] = false;

  // Redraw
  glutPostRedisplay();
}



void GLUTCommand(int cmd)
{
  // Execute command
  switch (cmd) {
  case DISPLAY_PARTICLES_TOGGLE_COMMAND: show_particles = !show_particles; break;
  case DISPLAY_PARTICLE_SPRINGS_TOGGLE_COMMAND: show_particle_springs = !show_particle_springs; break;
  case DISPLAY_PARTICLE_SOURCES_AND_SINKS_TOGGLE_COMMAND: show_particle_sources_and_sinks = !show_particle_sources_and_sinks; break;
  case DISPLAY_FACE_TOGGLE_COMMAND: show_faces = !show_faces; break;
  case DISPLAY_EDGE_TOGGLE_COMMAND: show_edges = !show_edges; break;
  case DISPLAY_BBOXES_TOGGLE_COMMAND: show_bboxes = !show_bboxes; break;
  case DISPLAY_LIGHTS_TOGGLE_COMMAND: show_lights = !show_lights; break;
  case DISPLAY_CAMERA_TOGGLE_COMMAND: show_camera = !show_camera; break;
  case SAVE_IMAGE_COMMAND: save_image = 1; break;
  case SAVE_VIDEO_COMMAND: save_video = save_video ^ 1; break;
  case QUIT_COMMAND: quit = 1; break;
  }

  // Mark window for redraw
  glutPostRedisplay();
}



void GLUTCreateMenu(void)
{
  // Display sub-menu
  int display_menu = glutCreateMenu(GLUTCommand);
  glutAddMenuEntry("Particles (P)", DISPLAY_PARTICLES_TOGGLE_COMMAND);
  glutAddMenuEntry("Particle springs (R)", DISPLAY_PARTICLE_SPRINGS_TOGGLE_COMMAND);
  glutAddMenuEntry("Particle sources and sinks (S)", DISPLAY_PARTICLE_SOURCES_AND_SINKS_TOGGLE_COMMAND);
  glutAddMenuEntry("Faces (F)", DISPLAY_FACE_TOGGLE_COMMAND);
  glutAddMenuEntry("Edges (E)", DISPLAY_EDGE_TOGGLE_COMMAND);
  glutAddMenuEntry("Bounding boxes (B)", DISPLAY_BBOXES_TOGGLE_COMMAND);
  glutAddMenuEntry("Lights (L)", DISPLAY_LIGHTS_TOGGLE_COMMAND);
  glutAddMenuEntry("Camera (C)", DISPLAY_CAMERA_TOGGLE_COMMAND);

  // Main menu
  glutCreateMenu(GLUTCommand);
  glutAddSubMenu("Display", display_menu);
  glutAddMenuEntry("Save Image (F1)", SAVE_IMAGE_COMMAND); 
  glutAddMenuEntry("Capture Video (F2)", SAVE_VIDEO_COMMAND);
 glutAddMenuEntry("Quit", QUIT_COMMAND);

  // Attach main menu to right mouse button
  glutAttachMenu(GLUT_RIGHT_BUTTON);
}



void GLUTInit(int *argc, char **argv)
{
  // Open window 
  glutInit(argc, argv);
  glutInitWindowPosition(100, 100);
  glutInitWindowSize(GLUTwindow_width, GLUTwindow_height);
  glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH); // | GLUT_STENCIL
  GLUTwindow = glutCreateWindow("OpenGL Viewer");
  glutFullScreen();

  // Initialize GLUT callback functions 
  glutIdleFunc(GLUTIdle);
  glutReshapeFunc(GLUTResize);
  glutDisplayFunc(GLUTRedraw);
  glutKeyboardFunc(GLUTKeyboard);
  glutKeyboardUpFunc(GLUTKeyboardUp);
  // glutSpecialFunc(GLUTSpecial);
  // glutSpecialUpFunc(GLUTSpecialUp);
  glutMouseFunc(GLUTMouse);
  glutMotionFunc(GLUTMotion);

  // Initialize graphics modes 
  glEnable(GL_NORMALIZE);
  glEnable(GL_LIGHTING);
  glEnable(GL_DEPTH_TEST);
  glLightModeli(GL_LIGHT_MODEL_LOCAL_VIEWER, GL_TRUE);

  // Create menus
  GLUTCreateMenu();
}




////////////////////////////////////////////////////////////
// SCENE READING
////////////////////////////////////////////////////////////


R3Scene *
ReadScene(const char *filename)
{
  // Allocate scene
  R3Scene *scene = new R3Scene();
  if (!scene) {
    fprintf(stderr, "Unable to allocate scene\n");
    return NULL;
  }

  // Read file
  if (!scene->Read(filename)) {
    fprintf(stderr, "Unable to read scene from %s\n", filename);
    return NULL;
  }

  // scene->Write(filename, scene->root); 

  // Remember initial camera
  camera = scene->camera;

  // Return scene
  return scene;
}



////////////////////////////////////////////////////////////
// PROGRAM ARGUMENT PARSING
////////////////////////////////////////////////////////////

int 
ParseArgs(int argc, char **argv)
{
  // Innocent until proven guilty
  int print_usage = 0;

  // Parse arguments
  argc--; argv++;
  while (argc > 0) {
    if ((*argv)[0] == '-') {
      if (!strcmp(*argv, "-help")) { print_usage = 1; }
      else if (!strcmp(*argv, "-soundtrack_enabled")) { soundtrack_enabled = 1; }
      else if (!strcmp(*argv, "-level_editor")) { level_editor = 1; minimap = 0;}
      else if (!strcmp(*argv, "-exit_immediately")) { quit = 1; }
      else if (!strcmp(*argv, "-output_image")) { argc--; argv++; output_image_name = *argv; }
      else if (!strcmp(*argv, "-video_prefix")) { argc--; argv++; video_prefix = *argv; }
      else if (!strcmp(*argv, "-euler")) integration_type = EULER_INTEGRATION; 
      else if (!strcmp(*argv, "-midpoint")) integration_type = MIDPOINT_INTEGRATION; 
      else if (!strcmp(*argv, "-rk4")) integration_type = RK4_INTEGRATION; 
      else if (!strcmp(*argv, "-adaptive_step_size")) integration_type = ADAPTIVE_STEP_SIZE_INTEGRATION; 
      else if (!strcmp(*argv, "-recordandquit")) { 
        argc--; argv++; num_frames_to_record = atoi(*argv); 
        GLUTwindow_width = 256;
        GLUTwindow_height = 256;
        save_video = 1;
      }
      else { fprintf(stderr, "Invalid program argument: %s", *argv); exit(1); }
      argv++; argc--;
    }
    else {
      if (!input_scene_name) input_scene_name = *argv;
      else { fprintf(stderr, "Invalid program argument: %s", *argv); exit(1); }
      argv++; argc--;
    }
  }

  // Check input_scene_name
  if (!input_scene_name || print_usage) {
    printf("Usage: game <input.scn> [-exit_immediately] [-output_image OUTPUT.JPG]  [-video_prefix VIDEO_DIR/PREFIX_] [-euler] [-midpoint] [-rk4] [-adaptive_step_size]  [-recordandquit NUM_FRAMES] [-v]\n");
    return 0;
  }

  // Return OK status
  return 1;
}




void DrawSidebar(R3Scene *scene) {
  R3Sidebar& sidebar(*scene->sidebar);
  // draw the background
  float xmin = GLUTwindow_width - scene->sidebar->width;
  glBegin(GL_QUADS);
  glColor3f(0, 0, 0);
  glVertex2f(xmin, 0);
  glColor3f(.3, .3, .3);
  glVertex2f(GLUTwindow_width, 0);
  glColor3f(0, 0, 0);
  glVertex2f(GLUTwindow_width, GLUTwindow_height);
  glColor3f(.2, .2, .2);
  glVertex2f(xmin, GLUTwindow_height);
  glEnd();
  
  int numButtons = scene->sidebar->buttons.size();
  double button_width = sidebar.width - 2*sidebar.border;
  for (int i = 0; i < numButtons; ++i) {
    R3Button& cur(*scene->sidebar->buttons[i]);
    {
      float xmin = GLUTwindow_width - sidebar.width + sidebar.border;
      float ymin = (i) * (button_width + sidebar.border) + sidebar.border;
      float ymax = ymin + button_width;
      float xmax = xmin + button_width;


      glBindTexture(GL_TEXTURE_2D, 0); 

      glBegin(GL_QUADS);
      if (i == scene->sidebar->selected_button) {
        glColor3f(1, 0, 0);
      }
      else {
        glColor3f(0, 0, 1);
      }
      glVertex3f(xmin, ymin, .01);
      glVertex3f(xmax, ymin, .01);
      glVertex3f(xmax, ymax, .01);
      glVertex3f(xmin, ymax, .01);
      glEnd();
    }
    {
      float button_width = sidebar.width - 4*sidebar.border;
      float xmin = GLUTwindow_width - sidebar.width + sidebar.border*2;
      float ymin = (i) * (button_width + sidebar.border*3) + sidebar.border*2;
      float ymax = ymin + button_width;
      float xmax = xmin + button_width;
      glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP);
      glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP);
      
      // Store the current matrix
      glPushMatrix();
      // Reset and transform the matrix.
      
      // Enable/Disable features
      glPushAttrib(GL_ENABLE_BIT);
      glEnable(GL_TEXTURE_2D);
      glDisable(GL_DEPTH_TEST);
      glDisable(GL_LIGHTING);
      glDisable(GL_BLEND);
      // Just in case we set all vertices to white.
      glColor4f(1,1,1,1);
      // Render the front quad
      glBindTexture(GL_TEXTURE_2D, cur.material->texture_index);
      glBegin(GL_QUADS);
      glTexCoord2f(1, 0); glVertex3f(xmin, ymin, .04);
      glTexCoord2f(1, 1); glVertex3f(xmax, ymin, .04);
      glTexCoord2f(0, 1); glVertex3f(xmax, ymax, .04);
      glTexCoord2f(0, 0); glVertex3f(xmin, ymax, .04);
      glEnd();
      glPopAttrib();
      glPopMatrix();
    }
  }
  
}

const int NButtons = 4;

int *ButtonVariables[] = {
  &camera_mode,
  &grab_mode,
  &blocks_mode,
  &coins_mode,
};

const char *ButtonIconFiles[] = {
  "camera.jpg",
  "grab.jpg",
  "platform.jpg",
  "coin.jpg"
};

void SetupSkybox(R3Scene *scene) {
  for (int i = 0; i < 5; ++i) {
    // Create material
    R3Material *material = new R3Material();
    *material = *scene->materials[0];
    material->texture_index = -1;

    // Get texture filename
    char buffer[2048];
    memset(buffer, 0, 2048);
    strcpy(buffer, skybox_path);
    strcat(buffer, scene->skyboxTexture);
    strcpy(material->texture_name, buffer);
    
    material->ka = R3Rgb(1, 1, 1, 0);
    material->kd = R3Rgb(1, 1, 1, 0);
    material->ks = R3Rgb(0, 0, 0, 0);
    material->kt = R3Rgb(0, 0, 0, 0);
    material->emission = R3Rgb(0, 0, 0, 0);
    material->shininess = 1;
    material->indexofrefraction = 1;
    
    // Read texture image
    material->texture = new R2Image();
    if (!material->texture->Read(buffer)) {
      fprintf(stderr, "not a good icon file: %s\n", buffer);
    }

    LoadMaterial(material);
    // Insert material
    skyboxMaterials.push_back(material);
  }
  

}

void SetupLevelEditor(R3Scene *scene) {
  for (int i = 0; i < NButtons; ++i) {
    // Create material
    R3Material *material = new R3Material();
    *material = *scene->materials[0];
    material->texture_index = -1;

    // Get texture filename
    char buffer[2048];
    memset(buffer, 0, 2048);
    strcpy(buffer, images_path);
    strcat(buffer, ButtonIconFiles[i]);
    strcpy(material->texture_name, buffer);
    
    material->ka = R3Rgb(1, 1, 1, 0);
    material->kd = R3Rgb(1, 1, 1, 0);
    material->ks = R3Rgb(0, 0, 0, 0);
    material->kt = R3Rgb(0, 0, 0, 0);
    material->emission = R3Rgb(0, 0, 0, 0);
    material->shininess = 1;
    material->indexofrefraction = 1;
    
    // Read texture image
    material->texture = new R2Image();
    if (!material->texture->Read(buffer)) {
      fprintf(stderr, "not a good icon file: %s\n", buffer);
    }
    
    LoadMaterial(material);
    // Insert material
    scene->materials.push_back(material);
    R3Button *button = new R3Button(ButtonVariables[i], material);
    scene->sidebar->buttons.push_back(button);
  }
  scene->sidebar->selected_button = 0;
  *scene->sidebar->buttons[0]->value = 1;
}

R3Camera GetMinimapCam(R3Scene *scene) {
  R3Camera ret(camera);
  ret.towards = R3negz_vector;
  ret.up = R3posy_vector;
  ret.right = R3posx_vector;
  ret.eye = scene->root->bbox.Centroid();
  ret.eye.Translate(-1*ret.towards*(scene->root->bbox.XLength()/2)/tan(ret.xfov));
  return ret;
}

////////////////////////////////////////////////////////////
// MAIN
////////////////////////////////////////////////////////////

void LoadLevel(const char *filename)
{
  if (soundtrack)
  {
    soundtrack->stop();
    soundtrack->drop();
  }

  if (scene) {
    delete scene;
  }
  // Read scene
  scene = ReadScene(filename);
  
  if (!scene) exit(-1);

  previous_time = 0.0f;
  
  SetupSkybox(scene); 
  minimap_cam = GetMinimapCam(scene);

  if (level_editor) {
    SetupLevelEditor(scene);
    scene->camera = minimap_cam;
    camera = minimap_cam;
  }

  if (soundtrack_enabled && scene->soundtrack)
  {
    char path[FILENAME_MAX];
    FilePath(path, scene->soundtrack);
    soundtrack = sound_engine->play2D(path, true, false, true);
    soundtrack->setVolume(0.35);
    soundtrack->setIsPaused(false);
  }
}

int 
main(int argc, char **argv)
{
  // Parse program arguments
  if (!ParseArgs(argc, argv)) exit(1);

  GetExecPath();

  // Initialize GLUT
  GLUTInit(&argc, argv);

  // Initialize sound shit
  sound_engine = createIrrKlangDevice();
  if (!sound_engine)
    return 0; // if there was an error creating the sound engine

  // Read scene
  LoadLevel(input_scene_name);

  // glut loop
  GLUTMainLoop();

  // Return success 
  return 0;
}

