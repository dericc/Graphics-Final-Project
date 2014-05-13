// Include file for the particle system



// Particle system integration types

enum {
  EULER_INTEGRATION,
  MIDPOINT_INTEGRATION,
  ADAPTIVE_STEP_SIZE_INTEGRATION,
  RK4_INTEGRATION
};



// Particle system functions

void UpdateParticles(R3Scene *scene, double current_time, double delta_time, int integration_type);
void EulerStep(R3Scene *scene, double time_step, vector<R3Point> *init_pos=NULL, vector<R3Vector> *init_vel=NULL);
R3Vector ComputeForce(R3Scene *scene, R3Particle *particle);
void GenerateParticles(R3Scene *scene, double current_time, double delta_time);
void RenderParticles(R3Scene *scene);
void CreateParticles(R3Scene *scene, R3Point position, int count, R3Material *mat, double velocity);


