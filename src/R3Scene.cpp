// Source file for the R3 scene graph class



// Include files

#include "R3/R3.h"
#include "R3Scene.h"

#include <iostream>

#include <istream>
#include <fstream>
#include <iterator>

vector<R3Material *> materials;

R3Point R3Player::Center(void) {
  R3Box box = *(node->shape->box);
  box.Transform(node->transformation);
  return (box.Min() + box.Max())/2;
}

// These are directions from the box's persepctive
R3Vector R3Player::Towards(void) {
  R3Box *box = node->shape->box;
  R3Vector towards = R3Point(box->XMax(), box->YMin(), box->ZMin()) - box->Min();
  towards.Normalize();
  return towards;
}

R3Vector R3Player::Right(void) {
  R3Box *box = node->shape->box;
  R3Vector right = box->Min() - R3Point(box->XMin(), box->YMin(), box->ZMax());
  right.Normalize();
  return right;
}

R3Vector R3Player::Up(void) {
  R3Box *box = node->shape->box;
  R3Vector up = R3Point(box->XMin(), box->YMax(), box->ZMin()) - box->Min();
  up.Normalize();
  return up;
}

void R3Camera::Rotate(R3Line axis, double angle) {
  eye.Rotate(axis, angle);
  right.Rotate(axis.Vector(), angle);
  towards.Rotate(axis.Vector(), angle);
  up.Rotate(axis.Vector(), angle);
}



R3Point R3Enemy::Center(void) {
  R3Box box = *(node->shape->box);
  box.Transform(node->transformation);
  return (box.Min() + box.Max())/2;
}

// These are directions from the box's persepctive
R3Vector R3Enemy::Towards(void) {
  R3Box *box = node->shape->box;
  R3Vector towards = R3Point(box->XMax(), box->YMin(), box->ZMin()) - box->Min();
  towards.Normalize();
  return towards;
}

R3Vector R3Enemy::Right(void) {
  R3Box *box = node->shape->box;
  R3Vector right = box->Min() - R3Point(box->XMin(), box->YMin(), box->ZMax());
  right.Normalize();
  return right;
}

R3Vector R3Enemy::Up(void) {
  R3Box *box = node->shape->box;
  R3Vector up = R3Point(box->XMin(), box->YMax(), box->ZMin()) - box->Min();
  up.Normalize();
  return up;
}

R3Point R3Goal::Center(void) {
  R3Box box = *(node->shape->box);
  box.Transform(node->transformation);
  return (box.Min() + box.Max())/2;
}

//R3Sidebar::
//R3Sidebar(void)

R3Scene::
R3Scene(void)
: bbox(R3null_box),
background(0,0,0,1),
ambient(0,0,0,1)
{
  // Setup default camera
  camera.eye = R3zero_point;
  camera.towards = R3negz_vector;
  camera.up = R3posy_vector;
  camera.right = R3posx_vector;
  camera.xfov = 0.0;
  camera.yfov = 0.0;
  camera.neardist = 0.01;
  camera.fardist = 100.0;
  
  // Create root node
  root = new R3Node();
  root->parent = NULL;
  root->transformation = R3identity_matrix;
  root->material = NULL;
  root->shape = NULL;
  root->bbox = R3null_box;
  
  player = NULL;
  death_y = -300;

  // Create sidebar
  sidebar = new R3Sidebar(180, 10);
  
  coin_shape = NULL;
  coin_material = NULL;

  next_level[0] = '\0';
  soundtrack[0] = '\0';
}



static R3Shape *
ReadShape(FILE *fp, int command_number, const char *filename)
{
  // Initialize result
  R3Shape *shape = NULL;
  
  // Read shape type
  char shape_type[1024];
  if (fscanf(fp, "%s", shape_type) != 1) {
    fprintf(stderr, "Unable to read shape type at command %d in file %s\n", command_number, filename);
    return NULL;
  }
  
  // Read shape args
  if (!strcmp(shape_type, "box")) {
    // Read sphere args
    double x1, y1, z1, x2, y2, z2;
    if (fscanf(fp, "%lf%lf%lf%lf%lf%lf", &x1, &y1, &z1, &x2, &y2, &z2) != 6) {
      fprintf(stderr, "Unable to read sphere args at command %d in file %s\n", command_number, filename);
      return NULL;
    }
    
    // Create shape
    shape = new R3Shape();
    shape->type = R3_BOX_SHAPE;
    shape->box = new R3Box(x1, y1, z1, x2, y2, z2);
    shape->sphere = NULL;
    shape->cylinder = NULL;
    shape->cone = NULL;
    shape->mesh = NULL;
    shape->segment = NULL;
    shape->circle = NULL;
  }
  else if (!strcmp(shape_type, "sphere")) {
    // Read sphere args
    double center_x, center_y, center_z, radius;
    if (fscanf(fp, "%lf%lf%lf%lf", &center_x, &center_y, &center_z, &radius) != 4) {
      fprintf(stderr, "Unable to read sphere args at command %d in file %s\n", command_number, filename);
      return NULL;
    }
    
    // Create shape
    shape = new R3Shape();
    shape->type = R3_SPHERE_SHAPE;
    shape->box = NULL;
    shape->sphere = new R3Sphere(R3Point(center_x, center_y, center_z), radius);
    shape->cylinder = NULL;
    shape->cone = NULL;
    shape->mesh = NULL;
    shape->segment = NULL;
    shape->circle = NULL;
  }
  else if (!strcmp(shape_type, "cylinder")) {
    // Read cylinder args
    double center_x, center_y, center_z, radius, height;
    if (fscanf(fp, "%lf%lf%lf%lf%lf", &center_x, &center_y, &center_z, &radius, &height) != 5) {
      fprintf(stderr, "Unable to read cylinder args at command %d in file %s\n", command_number, filename);
      return NULL;
    }
    
    // Create shape
    shape = new R3Shape();
    shape->type = R3_CYLINDER_SHAPE;
    shape->box = NULL;
    shape->sphere = NULL;
    shape->cylinder = new R3Cylinder(R3Point(center_x, center_y, center_z), radius, height);
    shape->cone = NULL;
    shape->mesh = NULL;
    shape->segment = NULL;
    shape->circle = NULL;
  }
  else if (!strcmp(shape_type, "cone")) {
    // Read cylinder args
    double center_x, center_y, center_z, radius, height;
    if (fscanf(fp, "%lf%lf%lf%lf%lf", &center_x, &center_y, &center_z, &radius, &height) != 5) {
      fprintf(stderr, "Unable to read cone args at command %d in file %s\n", command_number, filename);
      return NULL;
    }
    
    // Create shape
    shape = new R3Shape();
    shape->type = R3_CONE_SHAPE;
    shape->box = NULL;
    shape->sphere = NULL;
    shape->cylinder = NULL;
    shape->cone = new R3Cone(R3Point(center_x, center_y, center_z), radius, height);
    shape->mesh = NULL;
    shape->segment = NULL;
    shape->circle = NULL;
  }
  else if (!strcmp(shape_type, "mesh")) {
    // Read mesh args
    char meshname[1024];
    if (fscanf(fp, "%s", meshname) != 1) {
      fprintf(stderr, "Unable to read mesh args at command %d in file %s\n", command_number, filename);
      return NULL;
    }
    
    // Get mesh filename
    char buffer[2048];
    strcpy(buffer, filename);
    char *bufferp = strrchr(buffer, '/');
    if (bufferp) *(bufferp+1) = '\0';
    else buffer[0] = '\0';
    strcat(buffer, meshname);
    
    // Create mesh
    R3Mesh *mesh = new R3Mesh();
    if (!mesh) {
      fprintf(stderr, "Unable to allocate mesh\n");
      return 0;
    }
    
    // Read mesh file
    if (!mesh->Read(buffer)) {
      fprintf(stderr, "Unable to read mesh: %s\n", buffer);
      return 0;
    }
    
    // Create shape
    shape = new R3Shape();
    shape->type = R3_MESH_SHAPE;
    shape->box = NULL;
    shape->sphere = NULL;
    shape->cylinder = NULL;
    shape->cone = NULL;
    shape->mesh = mesh;
    shape->segment = NULL;
    shape->circle = NULL;
  }
  else if (!strcmp(shape_type, "line")) {
    // Read sphere args
    double x1, y1, z1, x2, y2, z2;
    if (fscanf(fp, "%lf%lf%lf%lf%lf%lf", &x1, &y1, &z1, &x2, &y2, &z2) != 6) {
      fprintf(stderr, "Unable to read sphere args at command %d in file %s\n", command_number, filename);
      return NULL;
    }
    
    // Create shape
    shape = new R3Shape();
    shape->type = R3_SEGMENT_SHAPE;
    shape->box = NULL;
    shape->sphere = NULL;
    shape->cylinder = NULL;
    shape->cone = NULL;
    shape->mesh = NULL;
    shape->segment = new R3Segment(R3Point(x1, y1, z1), R3Point(x2, y2, z2));
    shape->circle = NULL;
  }
  else if (!strcmp(shape_type, "circle")) {
    // Read circle args
    double center_x, center_y, center_z, normal_x, normal_y, normal_z, radius;
    if (fscanf(fp, "%lf%lf%lf%lf%lf%lf%lf", &center_x, &center_y, &center_z, &normal_x, &normal_y, &normal_z, &radius) != 7) {
      fprintf(stderr, "Unable to read circle args at command %d in file %s\n", command_number, filename);
      return NULL;
    }
    
    // Create shape
    shape = new R3Shape();
    shape->type = R3_CIRCLE_SHAPE;
    shape->box = NULL;
    shape->sphere = NULL;
    shape->cylinder = NULL;
    shape->cone = NULL;
    shape->mesh = NULL;
    shape->segment = NULL;
    shape->circle = new R3Circle(R3Point(center_x, center_y, center_z), radius, R3Vector(normal_x, normal_y, normal_z));
  }
  else {
    // Unrecognized shape type
    fprintf(stderr, "Invalid shape type (%s) at command %d in file %s\n", shape_type, command_number, filename);
    return NULL;
  }
  
  // Return created shape
  return shape;
}


///////////////////////////////////////////////////////////////////////////////////////////////////////
// WRITER: These files allow you to write a new scene from the existing scenes. ///////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////

void R3Scene:: 
WritePlayer(FILE *fp) {

  //Do nothing if no player
  if (player == NULL) return; 

  R3Node *cNode = player->node;
  //Calculates for material IDs
  R3Material *cMat = cNode->material; 
  R3Box cBox = *(cNode->shape->box); 
  cBox.Transform(cNode->transformation);

  int materialID = -1; 

  for (unsigned int j = 0; j < materials.size(); j++) {
    if (cMat == materials[j]) 
      materialID = j; 
  }

  fprintf(fp, "player %d %lf %lf %lf \n %lf %lf %lf \n %lf %lf \n",
    materialID, cBox.XMin(), cBox.YMin(), cBox.ZMin(), cBox.XMax(), cBox.YMax(), cBox.ZMax(), 
    player->max_speed, player->mass); 

  fprintf(fp, "\n"); 
  
}

void R3Scene:: 
WriteGoal(FILE *fp) {

  if (goal == NULL) return; 

  R3Node *cNode = goal->node; 
  R3Box cBox = *(cNode->shape->box); 
  cBox.Transform(cNode->transformation);

  R3Material *cMat = cNode->material; 
  int materialID = -1; 

  for (unsigned int j = 0; j < materials.size(); j++) {
    if (cMat == materials[j]) 
      materialID = j; 
  }
  
  fprintf(fp, "goal %d %lf %lf %lf %lf %lf %lf \n", 
    materialID, cBox.XMin(), cBox.YMin(), cBox.ZMin(), cBox.XMax(), cBox.YMax(), cBox.ZMax()); 

  fprintf(fp, "\n"); 
}

void R3Scene:: 
WriteEnemies(FILE *fp) {

  for (unsigned int i = 0; i < enemies.size(); i++) {
    R3Enemy *cEnemy = enemies[i]; 
    R3Node *cNode = cEnemy->node; 
    
    R3Material *cMat = cNode->material; 
    int materialID = -1; 

    for (unsigned int j = 0; j < materials.size(); j++) {
      if (cMat == materials[j]) 
        materialID = j; 
    }

    R3Box cBox = *(cNode->shape->box); 
    cBox.Transform(cNode->transformation);

    int moveLeftInt; 
    if (cEnemy->moveLeft)
      moveLeftInt = 1; 
    else 
      moveLeftInt = 0; 

    fprintf(fp, "enemy %d %lf %lf %lf \n %lf %lf %lf \n %d %lf %d %d %lf %lf \n", 
      materialID, cBox.XMin(), cBox.YMin(), cBox.ZMin(), cBox.XMax(), cBox.YMax(), cBox.ZMax(), 
      moveLeftInt, cEnemy->speed, cEnemy->is_jumping, cEnemy->is_following, cEnemy->jumpHeight, cEnemy->mass); 

  }

  fprintf(fp, "\n"); 
}

void R3Scene:: 
WriteFires(FILE *fp) {

  for (unsigned int i = 0; i < fires.size(); i++) {
    R3Fire *cFire = fires[i]; 

    R3Point cPos = cFire->position; 

    fprintf(fp, "fire %lf %lf %lf \n", 
      cPos.X(), cPos.Y(), cPos.Z()); 
  }

  fprintf(fp, "\n"); 
}


void R3Scene:: 
WriteMaterials(FILE *fp) {

  for (unsigned int i = 0; i < materials.size(); i++) {
    R3Material *cMat = materials[i]; 

    R3Rgb ka = cMat->ka; 
    R3Rgb kd = cMat->kd; 
    R3Rgb ks = cMat->ks; 
    R3Rgb kt = cMat->kt; 
    R3Rgb e = cMat->emission; 

    fprintf(fp, "material %lf %lf %lf \n %lf %lf %lf \n %lf %lf %lf \n %lf %lf %lf \n %lf %lf %lf \n %lf %lf %s \n", 
      ka.Red(), ka.Green(), ka.Blue(), 
      kd.Red(), kd.Green(), kd.Blue(), 
      ks.Red(), ks.Green(), ks.Blue(), 
      kt.Red(), kt.Green(), kt.Blue(), 
      e.Red(), e.Green(), e.Blue(), cMat->shininess, cMat->indexofrefraction, cMat->texture_name); 
  }

  fprintf(fp, "\n"); 
}

void R3Scene:: 
WritePlatforms(FILE *fp) {

  for (unsigned int i = 0; i < platforms.size(); i++) {

    R3Platform *cPlatform = platforms[i]; 
    R3Node *cNode = cPlatform->node; 

    R3Material *cMaterial = cNode->material; 
    int materialID = -1; 
    for (unsigned int j = 0; j < materials.size(); j++) {
      if (cMaterial == materials[j]) 
        materialID = j; 
    }

    R3Box *cBox = cNode->shape->box; 

    R3Point p1 = cBox->Min(); 
    R3Point p2 = cBox->Max(); 
    R3Point p3 = cPlatform->end; 

    fprintf(fp, "platform %d %lf %lf %lf \n %lf %lf %lf \n %lf %lf %lf \n %lf \n", 
      materialID, p1.X(), p1.Y(), p1.Z(), 
      p2.X(), p2.Y(), p2.Z(), 
      p3.X(), p3.Y(), p3.Z(), 
      cPlatform->max_speed); 
  }
  fprintf(fp, "\n"); 
}

void R3Scene::
WriteCoins(FILE *fp) {

  if (coin_material != NULL) {
    R3Material *cMaterial = coin_material; 
    int materialID = -1; 
    for (unsigned int j = 0; j < materials.size(); j++) {
      if (cMaterial == materials[j]) 
        materialID = j; 
    }
    fprintf(fp, "coin_material %d \n", materialID); 
  }

  for (unsigned int i = 0; i < coins.size(); i++) {

    R3Coin *cCoin = coins[i]; 
    R3Node *cNode = cCoin->node; 

    R3Material *cMaterial = cNode->material; 
    int materialID = -1; 
    for (unsigned int j = 0; j < materials.size(); j++) {
      if (cMaterial == materials[j]) 
        materialID = j; 
    }

    R3Point cPos = cCoin->position; 

    fprintf(fp, "coin %lf %lf %lf \n", cPos.X(), cPos.Y(), cPos.Z()); 
  }

  fprintf(fp, "\n"); 
}


void R3Scene::
WriteSkybox(FILE *fp) {

  if (skyboxTexture != NULL) {
    fprintf(fp, "skybox %s \n", skyboxTexture); 

    fprintf(fp, "\n"); 
  }
}

void R3Scene::
WriteSoundtrack(FILE *fp) {
  if (soundtrack != NULL) {
    fprintf(fp, "soundtrack %s \n", soundtrack);
  }
}

void R3Scene::
WriteNextLevel(FILE *fp) {
  if (next_level != NULL) {
    fprintf(fp, "next_level %s \n", next_level);
  }
}

void R3Scene::
WriteLights(FILE *fp) {

  for (unsigned int i = 0; i < lights.size(); i++) {
    R3Light *cLight = lights[i]; 

    R3Rgb cColor = cLight->color; 
    R3Vector cDirect = cLight->direction; 
    R3Point cPos = cLight->position; 

    if (cLight->type == R3_DIRECTIONAL_LIGHT) {
      fprintf(fp, "dir_light %lf %lf %lf \n %lf %lf %lf \n", 
        cColor.Red(), cColor.Green(), cColor.Blue(), 
        cDirect.X(), cDirect.Y(), cDirect.Z()); 
    }

    if (cLight->type == R3_POINT_LIGHT) {
      fprintf(fp, "point_light %lf %lf %lf \n %lf %lf %lf \n %lf %lf %lf \n", 
        cColor.Red(), cColor.Green(), cColor.Blue(), 
        cPos.X(), cPos.Y(), cPos.Z(), 
        cLight->constant_attenuation, cLight->linear_attenuation, cLight->quadratic_attenuation); 
    }

    if (cLight->type == R3_SPOT_LIGHT) {

    }

    if (cLight->type == R3_AREA_LIGHT) {

    }
  }

  fprintf(fp, "\n"); 
}

void R3Scene::
WriteNode(FILE *fp, R3Node *node) {


  for (unsigned int i = 0; i < node->children.size(); i++)
  { 
    //Recursively writes all the previous nodes first
    WriteNode(fp, node->children[i]); 

  }

  //Skip redrawing the player node
  if (player != NULL) {
    if (node == player->node) return; 
  }

  //Skip redrawing goal node, coin node, enemy node
  if (node->is_goal || node->is_coin || node->is_enemy) 
    return; 

  //Skip redrawing the platform nodes
  for (unsigned int j = 0; j < platforms.size(); j++) {
    
    R3Node *platNode = platforms[j]->node; 
    if (node == platNode) return; 
  }
    //Skip nodes without shapes
  if (node->shape == NULL) return; 



  if (node->shape->type == R3_BOX_SHAPE) {
    R3Box *cBox = node->shape->box; 

      //Calculates for material IDs
    R3Material *cMaterial = node->material; 
    int materialID = -1; 

    for (unsigned int j = 0; j < materials.size(); j++) {
      if (cMaterial == materials[j]) 
        materialID = j; 
    }

    fprintf(fp, "box %d %lf %lf %lf %lf %lf %lf \n", materialID, 
      cBox->XMin(), cBox->YMin(), cBox->ZMin(), cBox->XMax(), cBox->YMax(), cBox->ZMax()); 
  }

}

int R3Scene::
Write(const char *filename, R3Node *node) {

  const char *newfile = "../levels/output.scn"; 

    // Open file
  FILE *fp;
  if (!(fp = fopen(newfile, "w"))) {
    fprintf(stderr, "Unable to open file %s", newfile);
    return 0;
  }

  WriteMaterials(fp); 
  WriteLights(fp); 

  //Main node loop of objects
  WriteNode(fp, node); 
  fprintf(fp, "\n"); 

  WritePlatforms(fp); 
  WriteCoins(fp); 
  WritePlayer(fp); 
  WriteEnemies(fp);
  WriteGoal(fp); 
  WriteSkybox(fp); 
  WriteFires(fp); 
  WriteSoundtrack(fp);
  WriteNextLevel(fp);

  fclose(fp); 

  return 1; 
}  

int R3Scene::
Read(const char *filename, R3Node *node)
{
  // Open file
  FILE *fp;
  if (!(fp = fopen(filename, "r"))) {
    fprintf(stderr, "Unable to open file %s", filename);
    return 0;
  }
  
  // // Create array of materials
  // vector<R3Material *> materials;
  
  // Create default material
  R3Material *default_material = new R3Material();
  default_material->ka = R3Rgb(0.2, 0.2, 0.2, 1);
  default_material->kd = R3Rgb(0.5, 0.5, 0.5, 1);
  default_material->ks = R3Rgb(0.5, 0.5, 0.5, 1);
  default_material->kt = R3Rgb(0.0, 0.0, 0.0, 1);
  default_material->emission = R3Rgb(0, 0, 0, 1);
  default_material->shininess = 10;
  default_material->indexofrefraction = 1;
  default_material->texture = NULL;
  default_material->id = 0;
  
  // Create stack of group information
  const int max_depth = 1024;
  R3Node *group_nodes[max_depth] = { NULL };
  R3Material *group_materials[max_depth] = { NULL };
  group_nodes[0] = (node) ? node : root;
  group_materials[0] = default_material;
  int depth = 0;
  
  // Read body
  char cmd[128];
  int command_number = 1;
  vector<R3Particle *> particles_for_springs;
  while (fscanf(fp, "%s", cmd) == 1) {
    if (cmd[0] == '#') {
      // Comment -- read everything until end of line
      do { cmd[0] = fgetc(fp); } while ((cmd[0] >= 0) && (cmd[0] != '\n'));
    }
    else if (!strcmp(cmd, "particle")) {
      // Read position and velocity
      R3Point position;
      R3Vector velocity;
      double mass, drag, elasticity, lifetime;
      int fixed, m;
      if (fscanf(fp, "%lf%lf%lf%lf%lf%lf%lf%d%lf%lf%lf%d",
                 &position[0], &position[1], &position[2],
                 &velocity[0], &velocity[1], &velocity[2],
                 &mass, &fixed, &drag, &elasticity, &lifetime, &m) != 12) {
        fprintf(stderr, "Unable to read particle at command %d in file %s\n", command_number, filename);
        return 0;
      }
      
      // Get material
      R3Material *material = group_materials[depth];
      if (m >= 0) {
        if (m < (int) materials.size()) {
          material = materials[m];
        }
        else {
          fprintf(stderr, "Invalid material id at particle command %d in file %s\n", command_number, filename);
          return 0;
        }
      }
      
      // Create particle
      R3Particle *particle = new R3Particle();
      particle->position = position;
      particle->velocity = velocity;
      particle->mass = mass;
      particle->fixed = (fixed) ? true : false;
      particle->drag = drag;
      particle->elasticity = elasticity;
      particle->lifetime = lifetime;
      particle->material = material;
      
      // Add particle to scene
      particles.push_back(particle);
      
      // Update scene bounding box
      bbox.Union(position);
      
      // Add to list of particles available for springs
      particles_for_springs.push_back(particle);
    }
    else if (!strcmp(cmd, "particle_source")) {
      // Read particle source parameters
      double mass, drag, elasticity, lifetime;
      double rate, velocity, angle_cutoff;
      int fixed, m;
      if (fscanf(fp, "%lf%d%lf%lf%lf%d%lf%lf%lf", &mass, &fixed, &drag, &elasticity, &lifetime, &m, &rate, &velocity, &angle_cutoff) != 9) {
        fprintf(stderr, "Unable to read particle source at command %d in file %s\n", command_number, filename);
        return 0;
      }
      
      // Read shape
      R3Shape *shape = ReadShape(fp, command_number, filename);
      if (!shape) {
        fprintf(stderr, "Unable to read particle source at command %d in file %s\n", command_number, filename);
        return 0;
      }
      
      // Get material
      R3Material *material = group_materials[depth];
      if (m >= 0) {
        if (m < (int) materials.size()) {
          material = materials[m];
        }
        else {
          fprintf(stderr, "Invalid material id at particle source command %d in file %s\n", command_number, filename);
          return 0;
        }
      }
      
      // Create particle source
      R3ParticleSource *source = new R3ParticleSource();
      source->mass = mass;
      source->fixed = (fixed) ? true : false;
      source->drag = drag;
      source->elasticity = elasticity;
      source->lifetime = lifetime;
      source->material = material;
      source->rate = rate;
      source->velocity = velocity;
      source->angle_cutoff = angle_cutoff;
      source->shape = shape;
      
      // Add particle source to scene
      particle_sources.push_back(source);
      
      // Update scene bounding box
      if (shape->type == R3_SEGMENT_SHAPE) bbox.Union(shape->segment->BBox());
      else if (shape->type == R3_BOX_SHAPE) bbox.Union(*(shape->box));
      else if (shape->type == R3_CIRCLE_SHAPE) bbox.Union(shape->circle->BBox());
      else if (shape->type == R3_SPHERE_SHAPE) bbox.Union(shape->sphere->BBox());
      else if (shape->type == R3_CYLINDER_SHAPE) bbox.Union(shape->cylinder->BBox());
      else if (shape->type == R3_CONE_SHAPE) bbox.Union(shape->cone->BBox());
      else if (shape->type == R3_MESH_SHAPE) bbox.Union(shape->mesh->bbox);
    }
    else if (!strcmp(cmd, "particle_sink")) {
      // Read sink parameters
      double intensity, ca, la, qa;
      if (fscanf(fp, "%lf%lf%lf%lf", &intensity, &ca, &la, &qa) != 4) {
        fprintf(stderr, "Unable to read particle sink at command %d in file %s\n", command_number, filename);
        return 0;
      }
      
      // Read shape
      R3Shape *shape = ReadShape(fp, command_number, filename);
      if (!shape) {
        fprintf(stderr, "Unable to read particle source at command %d in file %s\n", command_number, filename);
        return 0;
      }
      
      // Create particle sink
      R3ParticleSink *sink = new R3ParticleSink();
      sink->intensity = intensity;
      sink->constant_attenuation = ca;
      sink->linear_attenuation = la;
      sink->quadratic_attenuation = qa;
      sink->shape = shape;
      
      // Add particle sink to scene
      particle_sinks.push_back(sink);
      
      // Update scene bounding box
      if (shape->type == R3_SEGMENT_SHAPE) bbox.Union(shape->segment->BBox());
      else if (shape->type == R3_BOX_SHAPE) bbox.Union(*(shape->box));
      else if (shape->type == R3_CIRCLE_SHAPE) bbox.Union(shape->circle->BBox());
      else if (shape->type == R3_SPHERE_SHAPE) bbox.Union(shape->sphere->BBox());
      else if (shape->type == R3_CYLINDER_SHAPE) bbox.Union(shape->cylinder->BBox());
      else if (shape->type == R3_CONE_SHAPE) bbox.Union(shape->cone->BBox());
      else if (shape->type == R3_MESH_SHAPE) bbox.Union(shape->mesh->bbox);
    }
    else if (!strcmp(cmd, "particle_spring")) {
      // Read gravity parameters
      int id1, id2;
      double rest_length, ks, kd;
      if (fscanf(fp, "%d%d%lf%lf%lf", &id1, &id2, &rest_length, &ks, &kd) != 5) {
        fprintf(stderr, "Unable to read particle spring at command %d in file %s\n", command_number, filename);
        return 0;
      }
      
      // Get particles
      R3Particle *particle1 = particles_for_springs[id1];
      R3Particle *particle2 = particles_for_springs[id2];
      
      // Create particle spring
      R3ParticleSpring *spring = new R3ParticleSpring();
      spring->particles[0] = particle1;
      spring->particles[1] = particle2;
      spring->rest_length = rest_length;
      spring->ks = ks;
      spring->kd = kd;
      
      // Insert spring into particles
      particle1->springs.push_back(spring);
      particle2->springs.push_back(spring);
      
      // Insert spring into scene
      particle_springs.push_back(spring);
    }
    else if (!strcmp(cmd, "particle_gravity")) {
      // Read gravity parameters
      if (fscanf(fp, "%lf%lf%lf", &gravity[0], &gravity[1], &gravity[2]) != 3) {
        fprintf(stderr, "Unable to read particle gravity at command %d in file %s\n", command_number, filename);
        return 0;
      }
    }
    else if (!strcmp(cmd, "tri")) {
      // Read data
      int m;
      R3Point p1, p2, p3;
      if (fscanf(fp, "%d%lf%lf%lf%lf%lf%lf%lf%lf%lf", &m,
                 &p1[0], &p1[1], &p1[2], &p2[0], &p2[1], &p2[2], &p3[0], &p3[1], &p3[2]) != 10) {
        fprintf(stderr, "Unable to read triangle at command %d in file %s\n", command_number, filename);
        return 0;
      }
      
      // Get material
      R3Material *material = group_materials[depth];
      if (m >= 0) {
        if (m < (int) materials.size()) {
          material = materials[m];
        }
        else {
          fprintf(stderr, "Invalid material id at tri command %d in file %s\n", command_number, filename);
          return 0;
        }
      }
      
      // Create mesh
      R3Mesh *mesh = new R3Mesh();
      vector<R3MeshVertex *> vertices;
      vertices.push_back(mesh->CreateVertex(p1, R3zero_vector, R2zero_point));
      vertices.push_back(mesh->CreateVertex(p2, R3zero_vector, R2zero_point));
      vertices.push_back(mesh->CreateVertex(p3, R3zero_vector, R2zero_point));
      mesh->CreateFace(vertices);
      
      // Create shape
      R3Shape *shape = new R3Shape();
      shape->type = R3_MESH_SHAPE;
      shape->box = NULL;
      shape->sphere = NULL;
      shape->cylinder = NULL;
      shape->cone = NULL;
      shape->mesh = mesh;
      shape->segment = NULL;
      
      // Create shape node
      R3Node *node = new R3Node();
      node->transformation = R3identity_matrix;
      node->material = material;
      node->shape = shape;
      node->bbox = R3null_box;
      node->bbox.Union(p1);
      node->bbox.Union(p2);
      node->bbox.Union(p3);
      
      // Insert node
      group_nodes[depth]->bbox.Union(node->bbox);
      group_nodes[depth]->children.push_back(node);
      node->parent = group_nodes[depth];
    }
    else if (!strcmp(cmd, "box")) {
      // Read data
      int m;
      R3Point p1, p2;
      if (fscanf(fp, "%d%lf%lf%lf%lf%lf%lf", &m, &p1[0], &p1[1], &p1[2], &p2[0], &p2[1], &p2[2]) != 7) {
        fprintf(stderr, "Unable to read box at command %d in file %s\n", command_number, filename);
        return 0;
      }
      
      // Get material
      R3Material *material = group_materials[depth];
      if (m >= 0) {
        if (m < (int) materials.size()) {
          material = materials[m];
        }
        else {
          fprintf(stderr, "Invalid material id at box command %d in file %s\n", command_number, filename);
          return 0;
        }
      }
      
      // Create box
      R3Box *box = new R3Box(p1, p2);
      
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
      node->material = material;
      node->shape = shape;
      node->bbox = *box;
      node->is_obstacle = true;
      
      // Insert node
      group_nodes[depth]->bbox.Union(node->bbox);
      group_nodes[depth]->children.push_back(node);
      node->parent = group_nodes[depth];
    }
    else if (!strcmp(cmd, "sphere")) {
      // Read data
      int m;
      R3Point c;
      double r;
      if (fscanf(fp, "%d%lf%lf%lf%lf", &m, &c[0], &c[1], &c[2], &r) != 5) {
        fprintf(stderr, "Unable to read sphere at command %d in file %s\n", command_number, filename);
        return 0;
      }
      
      // Get material
      R3Material *material = group_materials[depth];
      if (m >= 0) {
        if (m < (int) materials.size()) {
          material = materials[m];
        }
        else {
          fprintf(stderr, "Invalid material id at sphere command %d in file %s\n", command_number, filename);
          return 0;
        }
      }
      
      // Create sphere
      R3Sphere *sphere = new R3Sphere(c, r);
      
      // Create shape
      R3Shape *shape = new R3Shape();
      shape->type = R3_SPHERE_SHAPE;
      shape->box = NULL;
      shape->sphere = sphere;
      shape->cylinder = NULL;
      shape->cone = NULL;
      shape->mesh = NULL;
      shape->segment = NULL;
      
      // Create shape node
      R3Node *node = new R3Node();
      node->transformation = R3identity_matrix;
      node->material = material;
      node->shape = shape;
      node->bbox = sphere->BBox();
      
      // Insert node
      group_nodes[depth]->bbox.Union(node->bbox);
      group_nodes[depth]->children.push_back(node);
      node->parent = group_nodes[depth];
    }
    else if (!strcmp(cmd, "cylinder")) {
      // Read data
      int m;
      R3Point c;
      double r, h;
      if (fscanf(fp, "%d%lf%lf%lf%lf%lf", &m, &c[0], &c[1], &c[2], &r, &h) != 6) {
        fprintf(stderr, "Unable to read cylinder at command %d in file %s\n", command_number, filename);
        return 0;
      }
      
      // Get material
      R3Material *material = group_materials[depth];
      if (m >= 0) {
        if (m < (int) materials.size()) {
          material = materials[m];
        }
        else {
          fprintf(stderr, "Invalid material id at cyl command %d in file %s\n", command_number, filename);
          return 0;
        }
      }
      
      // Create cylinder
      R3Cylinder *cylinder = new R3Cylinder(c, r, h);
      
      // Create shape
      R3Shape *shape = new R3Shape();
      shape->type = R3_CYLINDER_SHAPE;
      shape->box = NULL;
      shape->sphere = NULL;
      shape->cylinder = cylinder;
      shape->cone = NULL;
      shape->mesh = NULL;
      shape->segment = NULL;
      
      // Create shape node
      R3Node *node = new R3Node();
      node->transformation = R3identity_matrix;
      node->material = material;
      node->shape = shape;
      node->bbox = cylinder->BBox();
      
      // Insert node
      group_nodes[depth]->bbox.Union(node->bbox);
      group_nodes[depth]->children.push_back(node);
      node->parent = group_nodes[depth];
    }
    else if (!strcmp(cmd, "mesh")) {
      // Read data
      int m;
      char meshname[256];
      if (fscanf(fp, "%d%s", &m, meshname) != 2) {
        fprintf(stderr, "Unable to parse mesh command %d in file %s\n", command_number, filename);
        return 0;
      }
      
      // Get material
      R3Material *material = group_materials[depth];
      if (m >= 0) {
        if (m < (int) materials.size()) {
          material = materials[m];
        }
        else {
          fprintf(stderr, "Invalid material id at cone command %d in file %s\n", command_number, filename);
          return 0;
        }
      }
      
      // Get mesh filename
      char buffer[2048];
      strcpy(buffer, filename);
      char *bufferp = strrchr(buffer, '/');
      if (bufferp) *(bufferp+1) = '\0';
      else buffer[0] = '\0';
      strcat(buffer, meshname);
      
      // Create mesh
      R3Mesh *mesh = new R3Mesh();
      if (!mesh) {
        fprintf(stderr, "Unable to allocate mesh\n");
        return 0;
      }
      
      // Read mesh file
      if (!mesh->Read(buffer)) {
        fprintf(stderr, "Unable to read mesh: %s\n", buffer);
        return 0;
      }
      
      // Create shape
      R3Shape *shape = new R3Shape();
      shape->type = R3_MESH_SHAPE;
      shape->box = NULL;
      shape->sphere = NULL;
      shape->cylinder = NULL;
      shape->cone = NULL;
      shape->mesh = mesh;
      shape->segment = NULL;
      
      // Create shape node
      R3Node *node = new R3Node();
      node->transformation = R3identity_matrix;
      node->material = material;
      node->shape = shape;
      node->bbox = mesh->bbox;
      
      // Insert node
      group_nodes[depth]->bbox.Union(node->bbox);
      group_nodes[depth]->children.push_back(node);
      node->parent = group_nodes[depth];
    }
    else if (!strcmp(cmd, "cone")) {
      // Read data
      int m;
      R3Point c;
      double r, h;
      if (fscanf(fp, "%d%lf%lf%lf%lf%lf", &m, &c[0], &c[1], &c[2], &r, &h) != 6) {
        fprintf(stderr, "Unable to read cone at command %d in file %s\n", command_number, filename);
        return 0;
      }
      
      // Get material
      R3Material *material = group_materials[depth];
      if (m >= 0) {
        if (m < (int) materials.size()) {
          material = materials[m];
        }
        else {
          fprintf(stderr, "Invalid material id at cone command %d in file %s\n", command_number, filename);
          return 0;
        }
      }
      
      // Create cone
      R3Cone *cone = new R3Cone(c, r, h);
      
      // Create shape
      R3Shape *shape = new R3Shape();
      shape->type = R3_CONE_SHAPE;
      shape->box = NULL;
      shape->sphere = NULL;
      shape->cylinder = NULL;
      shape->cone = cone;
      shape->mesh = NULL;
      shape->segment = NULL;
      
      // Create shape node
      R3Node *node = new R3Node();
      node->transformation = R3identity_matrix;
      node->material = material;
      node->shape = shape;
      node->bbox = cone->BBox();
      
      // Insert node
      group_nodes[depth]->bbox.Union(node->bbox);
      group_nodes[depth]->children.push_back(node);
      node->parent = group_nodes[depth];
    }
    else if (!strcmp(cmd, "line")) {
      // Read data
      int m;
      R3Point p1, p2;
      if (fscanf(fp, "%d%lf%lf%lf%lf%lf%lf", &m, &p1[0], &p1[1], &p1[2], &p2[0], &p2[1], &p2[2]) != 7) {
        fprintf(stderr, "Unable to read line at command %d in file %s\n", command_number, filename);
        return 0;
      }
      
      // Get material
      R3Material *material = group_materials[depth];
      if (m >= 0) {
        if (m < (int) materials.size()) {
          material = materials[m];
        }
        else {
          fprintf(stderr, "Invalid material id at line command %d in file %s\n", command_number, filename);
          return 0;
        }
      }
      
      // Create segment
      R3Segment *segment = new R3Segment(p1, p2);
      
      // Create shape
      R3Shape *shape = new R3Shape();
      shape->type = R3_SEGMENT_SHAPE;
      shape->box = NULL;
      shape->sphere = NULL;
      shape->cylinder = NULL;
      shape->cone = NULL;
      shape->mesh = NULL;
      shape->segment = segment;
      
      // Create shape node
      R3Node *node = new R3Node();
      node->transformation = R3identity_matrix;
      node->material = material;
      node->shape = shape;
      node->bbox = segment->BBox();
      
      // Insert node
      group_nodes[depth]->bbox.Union(node->bbox);
      group_nodes[depth]->children.push_back(node);
      node->parent = group_nodes[depth];
    }
    else if (!strcmp(cmd, "begin")) {
      // Read data
      int m;
      double matrix[16];
      if (fscanf(fp, "%d%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf", &m,
                 &matrix[0], &matrix[1], &matrix[2], &matrix[3],
                 &matrix[4], &matrix[5], &matrix[6], &matrix[7],
                 &matrix[8], &matrix[9], &matrix[10], &matrix[11],
                 &matrix[12], &matrix[13], &matrix[14], &matrix[15]) != 17) {
        fprintf(stderr, "Unable to read begin at command %d in file %s\n", command_number, filename);
        return 0;
      }
      
      // Get material
      R3Material *material = group_materials[depth];
      if (m >= 0) {
        if (m < (int) materials.size()) {
          material = materials[m];
        }
        else {
          fprintf(stderr, "Invalid material id at cone command %d in file %s\n", command_number, filename);
          return 0;
        }
      }
      
      // Create new group node
      R3Node *node = new R3Node();
      node->transformation = R3Matrix(matrix);
      node->material = NULL;
      node->shape = NULL;
      node->bbox = R3null_box;
      
      // Push node onto stack
      depth++;
      group_nodes[depth] = node;
      group_materials[depth] = material;
    }
    else if (!strcmp(cmd, "end")) {
      // Pop node from stack
      R3Node *node = group_nodes[depth];
      depth--;
      
      // Transform bounding box
      node->bbox.Transform(node->transformation);
      
      // Insert node
      group_nodes[depth]->bbox.Union(node->bbox);
      group_nodes[depth]->children.push_back(node);
      node->parent = group_nodes[depth];
    }
    else if (!strcmp(cmd, "material")) {
      // Read data
      R3Rgb ka, kd, ks, kt, e;
      double n, ir;
      char texture_name[256];
      if (fscanf(fp, "%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%s",
                 &ka[0], &ka[1], &ka[2], &kd[0], &kd[1], &kd[2], &ks[0], &ks[1], &ks[2], &kt[0], &kt[1], &kt[2],
                 &e[0], &e[1], &e[2], &n, &ir, texture_name) != 18) {
        fprintf(stderr, "Unable to read material at command %d in file %s\n", command_number, filename);
        return 0;
      }
      
      // Create material
      R3Material *material = new R3Material();
      material->ka = ka;
      material->kd = kd;
      material->ks = ks;
      material->kt = kt;
      material->emission = e;
      material->shininess = n;
      material->indexofrefraction = ir;
      material->texture = NULL;

      strcpy(material->texture_name, texture_name); 

      // Read texture
      if (strcmp(texture_name, "0")) {
        // Get texture filename
        char buffer[2048];
        strcpy(buffer, filename);
        char *bufferp = strrchr(buffer, '/');
        if (bufferp) *(bufferp+1) = '\0';
        else buffer[0] = '\0';
        strcat(buffer, texture_name);
        
        // Read texture image
        material->texture = new R2Image();
        if (!material->texture->Read(buffer)) {
          fprintf(stderr, "Unable to read texture from %s at command %d in file %s\n", buffer, command_number, filename);
          return 0;
        }
      }

      // Insert material
      materials.push_back(material);
    }

    else if (!strcmp(cmd, "skybox")) {
      // Read data
      char texture_name[256]; 
      if (fscanf(fp, "%s", texture_name) != 1) {
        fprintf(stderr, "Unable to read material at command %d in file %s\n", command_number, filename);
        return 0;
      }
      
    strcpy(skyboxTexture, texture_name); 

    }

    else if (!strcmp(cmd, "dir_light")) {
      // Read data
      R3Rgb c;
      R3Vector d;
      if (fscanf(fp, "%lf%lf%lf%lf%lf%lf",
                 &c[0], &c[1], &c[2], &d[0], &d[1], &d[2]) != 6) {
        fprintf(stderr, "Unable to read directional light at command %d in file %s\n", command_number, filename);
        return 0;
      }
      
      // Normalize direction
      d.Normalize();
      
      // Create light
      R3Light *light = new R3Light();
      light->type = R3_DIRECTIONAL_LIGHT;
      light->color = c;
      light->position = R3Point(0, 0, 0);
      light->direction = d;
      light->radius = 0;
      light->constant_attenuation = 0;
      light->linear_attenuation = 0;
      light->quadratic_attenuation = 0;
      light->angle_attenuation = 0;
      light->angle_cutoff = M_PI;
      
      // Insert light
      lights.push_back(light);
    }
    else if (!strcmp(cmd, "point_light")) {
      // Read data
      R3Rgb c;
      R3Point p;
      double ca, la, qa;
      if (fscanf(fp, "%lf%lf%lf%lf%lf%lf%lf%lf%lf", &c[0], &c[1], &c[2], &p[0], &p[1], &p[2], &ca, &la, &qa) != 9) {
        fprintf(stderr, "Unable to read point light at command %d in file %s\n", command_number, filename);
        return 0;
      }
      
      // Create light
      R3Light *light = new R3Light();
      light->type = R3_POINT_LIGHT;
      light->color = c;
      light->position = p;
      light->direction = R3Vector(0, 0, 0);
      light->radius = 0;
      light->constant_attenuation = ca;
      light->linear_attenuation = la;
      light->quadratic_attenuation = qa;
      light->angle_attenuation = 0;
      light->angle_cutoff = M_PI;
      
      // Insert light
      lights.push_back(light);
    }
    else if (!strcmp(cmd, "spot_light")) {
      // Read data
      R3Rgb c;
      R3Point p;
      R3Vector d;
      double ca, la, qa, sc, sd;
      if (fscanf(fp, "%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf",
                 &c[0], &c[1], &c[2], &p[0], &p[1], &p[2], &d[0], &d[1], &d[2], &ca, &la, &qa, &sc, &sd) != 14) {
        fprintf(stderr, "Unable to read point light at command %d in file %s\n", command_number, filename);
        return 0;
      }
      
      // Normalize direction
      d.Normalize();
      
      // Create light
      R3Light *light = new R3Light();
      light->type = R3_SPOT_LIGHT;
      light->color = c;
      light->position = p;
      light->direction = d;
      light->radius = 0;
      light->constant_attenuation = ca;
      light->linear_attenuation = la;
      light->quadratic_attenuation = qa;
      light->angle_attenuation = sd;
      light->angle_cutoff = sc;
      
      // Insert light
      lights.push_back(light);
    }
    else if (!strcmp(cmd, "area_light")) {
      // Read data
      R3Rgb c;
      R3Point p;
      R3Vector d;
      double radius, ca, la, qa;
      if (fscanf(fp, "%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf",
                 &c[0], &c[1], &c[2], &p[0], &p[1], &p[2], &d[0], &d[1], &d[2], &radius, &ca, &la, &qa) != 13) {
        fprintf(stderr, "Unable to read area light at command %d in file %s\n", command_number, filename);
        return 0;
      }
      
      // Normalize direction
      d.Normalize();
      
      // Create light
      R3Light *light = new R3Light();
      light->type = R3_AREA_LIGHT;
      light->color = c;
      light->position = p;
      light->direction = d;
      light->radius = radius;
      light->constant_attenuation = ca;
      light->linear_attenuation = la;
      light->quadratic_attenuation = qa;
      light->angle_attenuation = 0;
      light->angle_cutoff = M_PI;
      
      // Insert light
      lights.push_back(light);
    }
    else if (!strcmp(cmd, "camera")) {
      // Read data
      double px, py, pz, dx, dy, dz, ux, uy, uz, xfov, neardist, fardist;
      if (fscanf(fp, "%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf", &px, &py, &pz, &dx, &dy, &dz, &ux, &uy, &uz, &xfov, &neardist, &fardist) != 12) {
        fprintf(stderr, "Unable to read camera at command %d in file %s\n", command_number, filename);
        return 0;
      }
      
      // Assign camera
      camera.eye = R3Point(px, py, pz);
      camera.towards = R3Vector(dx, dy, dz);
      camera.towards.Normalize();
      camera.up = R3Vector(ux, uy, uz);
      camera.up.Normalize();
      camera.right = camera.towards % camera.up;
      camera.right.Normalize();
      camera.up = camera.right % camera.towards;
      camera.up.Normalize();
      camera.xfov = xfov;
      camera.yfov = xfov;
      camera.neardist = neardist;
      camera.fardist = fardist;
    }
    else if (!strcmp(cmd, "include")) {
      // Read data
      char scenename[256];
      if (fscanf(fp, "%s", scenename) != 1) {
        fprintf(stderr, "Unable to read include command %d in file %s\n", command_number, filename);
        return 0;
      }
      
      // Get scene filename
      char buffer[2048];
      strcpy(buffer, filename);
      char *bufferp = strrchr(buffer, '/');
      if (bufferp) *(bufferp+1) = '\0';
      else buffer[0] = '\0';
      strcat(buffer, scenename);
      
      // Read scene from included file
      if (!Read(buffer, group_nodes[depth])) {
        fprintf(stderr, "Unable to read included scene: %s\n", buffer);
        return 0;
      }
    }
    else if (!strcmp(cmd, "background")) {
      // Read data
      double r, g, b;
      if (fscanf(fp, "%lf%lf%lf", &r, &g, &b) != 3) {
        fprintf(stderr, "Unable to read background at command %d in file %s\n", command_number, filename);
        return 0;
      }
      
      // Assign background color
      background = R3Rgb(r, g, b, 1);
    }
    else if (!strcmp(cmd, "ambient")) {
      // Read data
      double r, g, b;
      if (fscanf(fp, "%lf%lf%lf", &r, &g, &b) != 3) {
        fprintf(stderr, "Unable to read ambient at command %d in file %s\n", command_number, filename);
        return 0;
      }
      
      // Assign ambient color
      ambient = R3Rgb(r, g, b, 1);
    }
    else if (!strcmp(cmd, "player")) {
      // Read data
      int m;
      R3Point p1, p2;
      double max_speed;
      double mass;
      if (fscanf(fp, "%d%lf%lf%lf%lf%lf%lf%lf%lf", &m, &p1[0], &p1[1], &p1[2], &p2[0], &p2[1], &p2[2], &max_speed, &mass) != 9) {
        fprintf(stderr, "Unable to read player at command %d in file %s\n", command_number, filename);
        return 0;
      }
      
      // Get material
      R3Material *material = group_materials[depth];
      if (m >= 0) {
        if (m < (int) materials.size()) {
          material = materials[m];
        }
        else {
          fprintf(stderr, "Invalid material id at box command %d in file %s\n", command_number, filename);
          return 0;
        }
      }
      // Create box
      R3Box *box = new R3Box(p1, p2);
      
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
      node->material = material;
      node->shape = shape;
      node->bbox = *box;
      node->is_player = true;
      
      // Insert node
      group_nodes[depth]->bbox.Union(node->bbox);
      group_nodes[depth]->children.push_back(node);
      node->parent = group_nodes[depth];
      
      player = new R3Player(node, max_speed, mass);

      // set movement plane
      movement_plane = R3Plane(player->node->bbox.Centroid(), R3Vector(0, 0, 1));
    }
    else if (!strcmp(cmd, "goal")) {
      // Read data
      int m;
      R3Point p1, p2;
      if (fscanf(fp, "%d%lf%lf%lf%lf%lf%lf", &m, &p1[0], &p1[1], &p1[2], &p2[0], &p2[1], &p2[2]) != 7) {
        fprintf(stderr, "Unable to read goal at command %d in file %s\n", command_number, filename);
        return 0;
      }
      
      // Create box
      R3Box *box = new R3Box(p1, p2);
      
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
      node->is_goal = true;
      
      // Insert node
      group_nodes[depth]->bbox.Union(node->bbox);
      group_nodes[depth]->children.push_back(node);
      node->parent = group_nodes[depth];
      
      goal = new R3Goal(node);
    }
    else if (!strcmp(cmd, "enemy")) {
      // Read data
      int m;
      R3Point p1, p2;
      int moveLeftInt; 
      int isJumpingInt; 
      int isFollowingInt; 
      double jumpHeight; 
      double speed;
      double mass;
      if (fscanf(fp, "%d%lf%lf%lf%lf%lf%lf%d%lf%d%d%lf%lf", &m, &p1[0], &p1[1], &p1[2], &p2[0], &p2[1], &p2[2], &moveLeftInt, &speed, &isJumpingInt, &isFollowingInt, &jumpHeight, &mass) != 13) {
        fprintf(stderr, "Unable to read enemy at command %d in file %s\n", command_number, filename);
        return 0;
      }

      bool moveLeft; 
      if (moveLeftInt == 1) {
        moveLeft = true; 
      }
      else moveLeft = false; 

      bool is_jumping; 
      if (isJumpingInt == 1) {
        is_jumping = true; 
      }
      else is_jumping = false; 

      bool is_following; 
      if (isFollowingInt == 1) {
        is_following = true; 
      }
      else is_following = false; 
      
      // Get material
      R3Material *material = group_materials[depth];
      if (m >= 0) {
        if (m < (int) materials.size()) {
          material = materials[m];
        }
        else {
          fprintf(stderr, "Invalid material id at box command %d in file %s\n", command_number, filename);
          return 0;
        }
      }

      // Create box
      R3Box *box = new R3Box(p1, p2);
      
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
      node->material = material;
      node->shape = shape;
      node->bbox = *box;
      node->is_enemy = true; 
      node->is_obstacle = true;
      
      // Insert node
      group_nodes[depth]->bbox.Union(node->bbox);
      group_nodes[depth]->children.push_back(node);
      node->parent = group_nodes[depth];


      R3Enemy *e = new R3Enemy(node, moveLeft, speed, is_jumping, is_following, jumpHeight, mass);
      
      if (!moveLeft) 
        e->velocity = speed * e->Towards(); 
      else 
        e->velocity = speed * -e->Towards(); 

      e->inAir = true; 
      e->onPlatform = false; 
      
      enemies.push_back(e);
      node->enemy = e;
    }

    else if (!strcmp(cmd, "platform")) {
      // Read data
      int m;
      R3Point p1, p2, p3;
      double speed;
      if (fscanf(fp, "%d%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf", &m, &p1[0], &p1[1], &p1[2], &p2[0], &p2[1], &p2[2], &p3[0], &p3[1], &p3[2], &speed) != 11) {
        fprintf(stderr, "Unable to read platform at command %d in file %s\n", command_number, filename);
        return 0;
      }
      
      // Get material
      R3Material *material = group_materials[depth];
      if (m >= 0) {
        if (m < (int) materials.size()) {
          material = materials[m];
        }
        
        else {
          fprintf(stderr, "Invalid material id at box command %d in file %s\n", command_number, filename);
          return 0;
        }
      }
      
      // Create box
      R3Box *box = new R3Box(p1, p2);
      
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
      node->material = material;
      node->shape = shape;
      node->bbox = *box;
      node->is_platform = true;
      node->is_obstacle = true;
      
      // Insert node
      group_nodes[depth]->bbox.Union(node->bbox);
      group_nodes[depth]->children.push_back(node);
      node->parent = group_nodes[depth];
      
      
      // Create platform
      R3Platform *p = new R3Platform(node, speed, p1, p3);
      platforms.push_back(p);
      node->platform = p;
    }
    else if (!strcmp(cmd, "coin_material")) {
      int m;
      if (fscanf(fp, "%d", &m) != 1) {
        fprintf(stderr, "Unable to read coin_material at command %d in file %s\n", command_number, filename);
        return 0;
      }
      // Get material
      coin_material = group_materials[depth];
      if (m >= 0) {
        if (m < (int) materials.size()) {
          coin_material = materials[m];
        }
        else {
          fprintf(stderr, "Invalid material id at box command %d in file %s\n", command_number, filename);
          return 0;
        }
      }
      // Create box
      R3Cylinder *cyl = new R3Cylinder(R3null_point, 0.5, 0.2);
      
      // Create shape
      coin_shape = new R3Shape();
      coin_shape->type = R3_CYLINDER_SHAPE;
      coin_shape->box = NULL;
      coin_shape->sphere = NULL;
      coin_shape->cylinder = cyl;
      coin_shape->cone = NULL;
      coin_shape->mesh = NULL;
      coin_shape->segment = NULL;
    }
    else if (!strcmp(cmd, "coin")) {
      if (coin_shape == NULL || coin_material == NULL) {
        fprintf(stderr, "Define a coin material first!\n");
        return 0;
      }
      // Read data
      R3Point p;
      if (fscanf(fp, "%lf%lf%lf", &p[0], &p[1], &p[2]) != 3) {
        fprintf(stderr, "Unable to read box at command %d in file %s\n", command_number, filename);
        return 0;
      }
      
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
      node->material = coin_material;
      node->shape = coin_shape;
      node->bbox = coin_shape->cylinder->BBox();
      node->is_coin = true;
      node->coin = coin;
      
      coin->node = node;
      
      coins.push_back(coin);
      
      // Insert node
      group_nodes[depth]->bbox.Union(node->bbox);
      group_nodes[depth]->children.push_back(node);
      node->parent = group_nodes[depth];
    }
    else if (!strcmp(cmd, "death_y")) {
      // Read data
      double y;
      if (fscanf(fp, "%lf", &y) != 1) {
        fprintf(stderr, "Unable to read death_y at command %d in file %s\n", command_number, filename);
        return 0;
      }
      
      death_y = y;
    }
    else if (!strcmp(cmd, "fire")) {
      // Read data
      R3Point p;
      if (fscanf(fp, "%lf%lf%lf", &p[0], &p[1], &p[2]) != 3) {
        fprintf(stderr, "Unable to read fire at command %d in file %s\n", command_number, filename);
        return 0;
      }
      
      R3Fire *fire = new R3Fire();
      fire->position = p;
      fires.push_back(fire);
    }
    else if (!strcmp(cmd, "soundtrack")) {
      // Read data
      char music_file[256];
      if (fscanf(fp, "%s", music_file) != 1) {
        fprintf(stderr, "Unable to read soundtrack at command %d in file %s\n", command_number, filename);
        return 0;
      }
      strcpy(soundtrack, music_file);
    }
    else if (!strcmp(cmd, "next_level")) {
      // Read data
      char level_file[256];
      if (fscanf(fp, "%s", level_file) != 1) {
        fprintf(stderr, "Unable to read next_level at command %d in file %s\n", command_number, filename);
        return 0;
      }
      
      strcpy(next_level, level_file);
    }
    else {
      fprintf(stderr, "Unrecognized command %d in file %s: %s\n", command_number, filename, cmd);
      return 0;
    }
    
    // Increment command number
    command_number++;
  }
  
  // Update bounding box
  bbox.Union(root->bbox);
  
  // Provide default camera
  if (camera.xfov == 0) {
    double scene_radius = bbox.DiagonalRadius();
    R3Point scene_center = bbox.Centroid();
    camera.towards = R3Vector(0, 0, -1);
    camera.up = R3Vector(0, 1, 0);
    camera.right = R3Vector(1, 0, 0);
    camera.eye = scene_center - 3 * scene_radius * camera.towards;
    camera.xfov = 0.25;
    camera.yfov = 0.25;
    camera.neardist = 0.01 * scene_radius;
    camera.fardist = 100 * scene_radius;
  }
  
  // Provide default lights
  if (lights.size() == 0) {
    // Create first directional light
    R3Light *light = new R3Light();
    R3Vector direction(-3,-4,-5);
    direction.Normalize();
    light->type = R3_DIRECTIONAL_LIGHT;
    light->color = R3Rgb(1,1,1,1);
    light->position = R3Point(0, 0, 0);
    light->direction = direction;
    light->radius = 0;
    light->constant_attenuation = 0;
    light->linear_attenuation = 0;
    light->quadratic_attenuation = 0;
    light->angle_attenuation = 0;
    light->angle_cutoff = M_PI;
    lights.push_back(light);
    
    // Create second directional light
    light = new R3Light();
    direction = R3Vector(3,2,3);
    direction.Normalize();
    light->type = R3_DIRECTIONAL_LIGHT;
    light->color = R3Rgb(0.5, 0.5, 0.5, 1);
    light->position = R3Point(0, 0, 0);
    light->direction = direction;
    light->radius = 0;
    light->constant_attenuation = 0;
    light->linear_attenuation = 0;
    light->quadratic_attenuation = 0;
    light->angle_attenuation = 0;
    light->angle_cutoff = M_PI;
    lights.push_back(light);
  }
  
  materials = materials;
  
  // Close file
  fclose(fp);
  
  // Return success
  return 1;
}



