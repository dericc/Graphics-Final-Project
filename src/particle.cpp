// Source file for the particle system



// Include files

#include "R2/R2.h"
#include "R3/R3.h"
#include "R3Scene.h"
#include "particle.h"
#include "raytrace.h"

#include <vector>
#include <iostream>
#include <limits>

using namespace std;

#define EPSILON 0.001
#define ERROR_THRESH 0.001


////////////////////////////////////////////////////////////
// Random Number Generator
////////////////////////////////////////////////////////////

static double 
RandomNumber(void) 
{
#if defined(_WIN32)
  int r1 = rand();
  double r2 = ((double) rand()) / ((double) (RAND_MAX + 1));
  return (r1 + r2) / ((double) (RAND_MAX + 1));
#else
  return drand48();
#endif
}



////////////////////////////////////////////////////////////
// Generating Particles
////////////////////////////////////////////////////////////

void CreateParticles(R3Scene *scene, R3Point position, int count, R3Material *mat, double velocity)
{
  for (int i = 0; i < count; i++)
  {
    double z = (RandomNumber() - 0.5) * 2;
    double phi = RandomNumber() * 2 * 3.14159;
    double d = sqrt(1 - z * z);
    R3Vector n(d * cos(phi), d * sin(phi), z);
    n.Normalize();

    // create new particle
    R3Particle *particle      = new R3Particle();
    particle->drag            = 0;
    particle->elasticity      = 1;
    particle->fixed           = false;
    particle->lifetime        = 10;
    particle->mass            = 1;
    particle->material        = mat;
    particle->position        = position;
    particle->velocity        = n * velocity;
    scene->particles.push_back(particle);
  }
}

void GenerateParticles(R3Scene *scene, double current_time, double delta_time)
{
  // Generate new particles for every source
  for (int i = 0; i < scene->NParticleSources(); i++) 
  {
    R3ParticleSource *source = scene->ParticleSource(i);

    // calculate number of particles to create 
    double start_time = current_time - delta_time;
    int num_particles = floor(current_time * source->rate + 0.5) 
                      - floor(start_time   * source->rate + 0.5);

    for (int j = 0; j < num_particles; j++)
    {
      switch(source->shape->type)
      {
        case R3_SPHERE_SHAPE:
        {
          // pick a random point on the surface
          R3Sphere *sphere = source->shape->sphere;
          double z = (RandomNumber() - 0.5) * 2 * sphere->Radius();
          double phi = RandomNumber() * 2 * 3.14159;
          double d = sqrt(sphere->Radius() * sphere->Radius() - z * z);
          R3Vector n(d * cos(phi), d * sin(phi), z);
          R3Point p = sphere->Center() + n;
          n.Normalize();

          // rotate the velocity vector in a random direction by a random amount less than the cutoff angle
          R3Vector tmp;
          if (n.X() != 0 && n.Y() != 0)
            tmp = R3Vector(0, 0, 1);
          else
            tmp = R3Vector(0, 1, 0);
          R3Vector a = n;
          a.Cross(tmp);
          a.Normalize();
          double t1 = RandomNumber() * 2 * 3.14159;
          double t2 = RandomNumber() * sin(source->angle_cutoff);
          R3Vector v = a;
          v.Rotate(n, t1);
          tmp = v;
          tmp.Cross(n);
          v.Rotate(tmp, acos(t2));

          // create new particle
          R3Particle *particle      = new R3Particle();
          particle->drag            = source->drag;
          particle->elasticity      = source->elasticity;
          particle->fixed           = source->fixed;
          particle->lifetime        = source->lifetime;
          particle->mass            = source->mass;
          particle->material        = source->material;
          particle->position        = p;
          particle->velocity        = v * source->velocity;
          particle->lifetime        = source->lifetime;
          scene->particles.push_back(particle);
        }
        case R3_BOX_SHAPE:
        {
          break;
        }
        case R3_CYLINDER_SHAPE:
        {
          break;
        }
        case R3_CONE_SHAPE:
        {
          break;
        }
        case R3_MESH_SHAPE:
        {
          break;
        }
        case R3_CIRCLE_SHAPE:
        {
          R3Circle *circle = source->shape->circle;

          // find an orthogonal basis for the circle plane
          R3Vector n = circle->Normal();
          R3Vector tmp;
          if (n.X() != 0 || n.Y() != 0)
            tmp = R3Vector(0, 0, 1);
          else
            tmp = R3Vector(0, 1, 0);
          R3Vector b1 = n;
          b1.Cross(tmp);
          b1.Normalize();
          R3Vector b2 = b1;
          b2.Cross(n);

          // pick a random point on the circle using rejection sampling
          double r1;
          double r2;
          while (true)
          {
            r1 = (RandomNumber() - 0.5) * 2 * circle->Radius();
            r2 = (RandomNumber() - 0.5) * 2 * circle->Radius();
            if (r1 * r1 + r2 * r2 <= circle->Radius() * circle->Radius())
              break;
          }

          R3Point p = circle->Center() + r1 * b1 + r2 * b2;

          // rotate velocity in a random direction by a random angle less than the cutoff angle
          double t1 = RandomNumber() * 2 * 3.14159;
          double t2 = RandomNumber() * sin(source->angle_cutoff);
          R3Vector v = b1;
          v.Rotate(b2, acos(t2));
          v.Rotate(n, t1);

          // create new particle
          R3Particle *particle      = new R3Particle();
          particle->drag            = source->drag;
          particle->elasticity      = source->elasticity;
          particle->fixed           = source->fixed;
          particle->lifetime        = source->lifetime;
          particle->mass            = source->mass;
          particle->material        = source->material;
          particle->position        = p;
          particle->velocity        = v * source->velocity;
          particle->lifetime        = source->lifetime;
          scene->particles.push_back(particle);
          break;
        }
        case R3_SEGMENT_SHAPE:
        {
          break;
        }
        case R3_NUM_SHAPE_TYPES:
        {
          break;
        }
      }
    }
  }
}



////////////////////////////////////////////////////////////
// Updating Particles
////////////////////////////////////////////////////////////

void UpdateParticles(R3Scene *scene, double current_time, double delta_time, int integration_type)
{
  switch(integration_type)
  {
    case EULER_INTEGRATION:
    {
      // a single Euler step
      EulerStep(scene, delta_time);
      break;
    }
    case MIDPOINT_INTEGRATION:
    {
      // save original particle properties
      vector<R3Point> ppos(scene->NParticles());
      vector<R3Vector> pvel(scene->NParticles());

      for (int i = 0; i < scene->NParticles(); i++) 
      {
        R3Particle *particle = scene->Particle(i);

        if (!particle->fixed)
        {
          ppos[i] = particle->position;
          pvel[i] = particle->velocity;
        }
      }

      // one normal half step
      EulerStep(scene, delta_time / 2);
      // one full step using force calculated at midpoint values, but overriding initial position/velocity
      EulerStep(scene, delta_time, &ppos, &pvel);
      break;
    }
    case ADAPTIVE_STEP_SIZE_INTEGRATION:
    {
      double error; // average position error
      int nsteps = 1; // number of euler steps
      vector<R3Point> ppos(scene->NParticles());
      vector<R3Point> fpos(scene->NParticles());
      vector<R3Vector> fvel(scene->NParticles());

      // save original particle properties
      for (int i = 0; i < scene->NParticles(); i++) 
      {
        R3Particle *particle = scene->Particle(i);

        if (!particle->fixed)
        {
          fpos[i] = particle->position;
          fvel[i] = particle->velocity;
        }
      }

      // loop w/ increasing number of steps
      do
      {
        for (int i = 0; i < scene->NParticles(); i++) 
        {
          R3Particle *particle = scene->Particle(i);

          if (!particle->fixed)
          {
            particle->position = fpos[i];
            particle->velocity = fvel[i];
          }
        }

        // step n times
        for (int j = 0; j < nsteps; j++)
        {
          EulerStep(scene, delta_time / nsteps);
        }

        // for first step, infinite error
        if (nsteps == 1)
          error = numeric_limits<double>::infinity();
        // otherwise, average position error between this loop and the last one over all non-fixed particles
        else
        {
          error = 0;
          int count = 0;
          for (int i = 0; i < scene->NParticles(); i++)
          {
            R3Particle *particle = scene->Particle(i);
            if (!particle->fixed)
            {
              error += (particle->position - ppos[i]).Length();
              ppos[i] = particle->position;
              count++;
            }
          }
          error /= count;
        }

        // double the number of steps for the next iteration
        nsteps *= 2;
      } while (error > ERROR_THRESH);
      break;
    }
  }

  // delete particles that have reached a sink or have expired
  for (vector<R3Particle *>::iterator i = scene->particles.begin(); i != scene->particles.end(); ) 
  {
    R3Particle *particle = *i;

    // avoids deleting particles that started with a 0 or negative lifetime
    if (particle->lifetime > 0)
    {
      particle->lifetime -= delta_time;
      if (particle->lifetime <= 0)
      {
        scene->particles.erase(i);
        delete particle;
        continue;
      }
    }

    bool del = false;

    for (int j = 0; j < scene->NParticleSinks(); j++)
    {
      if (del)
        break;

      R3ParticleSink *sink = scene->ParticleSink(j);

      switch(sink->shape->type)
      {
        case R3_SPHERE_SHAPE:
        {
          // delete if distance between particle and sphere center is less than the radius
          R3Sphere *sphere = sink->shape->sphere;
          R3Vector diff = sphere->Center() - particle->position;
          if (diff.Length() < sphere->Radius())
          {
            scene->particles.erase(i);
            delete particle;
            del = true;
          }
        }
        case R3_BOX_SHAPE:
        {
          break;
        }
        case R3_CYLINDER_SHAPE:
        {
          break;
        }
        case R3_CONE_SHAPE:
        {
          break;
        }
        case R3_MESH_SHAPE:
        {
          break;
        }
        case R3_CIRCLE_SHAPE:
        {
          break;
        }
        case R3_SEGMENT_SHAPE:
        {
          break;
        }
        case R3_NUM_SHAPE_TYPES:
        {
          break;
        }
      }
    }
    if (del)
      continue;

    i++;
  }
}

void EulerStep(R3Scene *scene, double time_step, vector<R3Point> *init_pos, vector<R3Vector> *init_vel)
{
  // buffers to store updated particle values
  vector<R3Point> npos(scene->NParticles());
  vector<R3Vector> nvel(scene->NParticles());

  for (int i = 0; i < scene->NParticles(); i++) 
  {
    R3Particle *particle = scene->Particle(i);

    if (!particle->fixed)
    {
      // temporary particle values
      R3Point tpos;
      R3Vector tvel;

      // override original position/velocity if arguments provided
      if (init_pos == NULL)
        tpos = particle->position;
      else
        tpos = (*init_pos)[i];

      if (init_vel == NULL)
        tvel = particle->velocity;
      else
        tvel = (*init_vel)[i];

      R3Intersection collision;
      double ttime_step = time_step;
      bool is_coll;
      do
      {
        // cast a ray from the particle along its velocity vector
        R3Ray coll_ray(tpos, tvel);
        collision = ComputeIntersection(scene, scene->Root(), coll_ray, numeric_limits<double>::infinity());
        double coll_t = collision.t / tvel.Length();
        // if the particle is going to collide within the remaining time step...
        is_coll = collision.hit && coll_t < ttime_step;
        if (is_coll)
        {
          // move particle to collision location plus a normal epsilon
          tpos = collision.position + collision.normal * EPSILON;

          // calculate reflected velocity
          R3Vector v = tvel;
          R3Vector parallel = v.Dot(collision.normal) * collision.normal;
          R3Vector tangent = v - parallel;
          tvel = tangent - (parallel * particle->elasticity);

          // subtract simulated time from current time step
          ttime_step -= coll_t;
        }
        // repeat while there was a collision
      } while (is_coll);

      // move the particle along its velocity vector for the remainder of the time step and update velocity using current values
      npos[i] = tpos + tvel * ttime_step;
      nvel[i] = tvel + ComputeForce(scene, particle) * ttime_step / particle->mass;
    }
  } 

  // actually update the particles from the temp array
  for (int i = 0; i < scene->NParticles(); i++)
  {
    R3Particle *particle = scene->Particle(i);
    if (!particle->fixed)
    {
      particle->position = npos[i];
      particle->velocity = nvel[i];
    }
  }
}

R3Vector ComputeForce(R3Scene *scene, R3Particle *particle)
{
  R3Vector f(0, 0, 0);

  // gravitational attractive force between each pair of particles
  // for (int i = 0; i < scene->NParticles(); i++) 
  // {
  //   R3Particle *other = scene->Particle(i);
  //   if (other != particle)
  //   {
  //     R3Vector diff = other->position - particle->position;
  //     double d = diff.Length();
  //     diff.Normalize();
  //     R3Vector fg = diff * 6.67384e-11 * particle->mass * other->mass / (d * d);
  //     f += fg;
  //   }
  // }

  // sink forces
  for (int i = 0; i < scene->NParticleSinks(); i++)
  {
    R3ParticleSink *sink = scene->ParticleSink(i);

    switch(sink->shape->type)
    {
      case R3_SPHERE_SHAPE:
      {
        R3Sphere *sphere = sink->shape->sphere;
        R3Vector diff = sphere->Center() - particle->position;
        double d = diff.Length() - sphere->Radius();
        diff.Normalize();
        if (d < 0)
        {
          diff = -diff;
          d = -d;
        }
        R3Vector fs = diff * sink->intensity / (sink->constant_attenuation + sink->linear_attenuation * d + sink->quadratic_attenuation * d * d);
        f += fs;
      }
      case R3_BOX_SHAPE:
      {
        break;
      }
      case R3_CYLINDER_SHAPE:
      {
        break;
      }
      case R3_CONE_SHAPE:
      {
        break;
      }
      case R3_MESH_SHAPE:
      {
        break;
      }
      case R3_CIRCLE_SHAPE:
      {
        break;
      }
      case R3_SEGMENT_SHAPE:
      {
        break;
      }
      case R3_NUM_SHAPE_TYPES:
      {
        break;
      }
    }
  }

  // spring forces
  for (unsigned int i = 0; i < particle->springs.size(); i++)
  {
    R3ParticleSpring *spring = particle->springs[i];

    // get particle on other side of spring
    R3Particle *other;
    if (particle == spring->particles[0])
      other = spring->particles[1];
    else
      other = spring->particles[0];

    // calculate force
    R3Vector diff = other->position - particle->position;
    double d = diff.Length();
    diff /= d ;
    R3Vector fsp = (spring->ks * (d - spring->rest_length)
                  + spring->kd * (other->velocity - particle->velocity).Dot(diff)) * diff;
    f += fsp;
  }

  // scene gravity
  f += scene->gravity * particle->mass;
  // drag
  f -= particle->velocity * particle->drag;
  return f;
}



////////////////////////////////////////////////////////////
// Rendering Particles
////////////////////////////////////////////////////////////

// void RenderParticles(R3Scene *scene)
// {
//   // Draw every particle
 
//     glPushMatrix();
//   // REPLACE CODE HERE
//   glDisable(GL_LIGHTING);
//   glPointSize(5);


//    glBindTexture(GL_TEXTURE_2D, scene->player->node->material->texture_index);
//   for (int i = 0; i < scene->NParticles(); i++) {
//     R3Particle *particle = scene->Particle(i);
//     // glColor3d(particle->material->kd[0], particle->material->kd[1], particle->material->kd[2]);
//     glBegin(GL_QUADS);
//     glColor3d(1, 0, 0); 
//     const R3Point& position = particle->position;
//     glVertex2f(position[0]-.2, position[1]-.2);
//     glVertex2f(position[0]+.2, position[1]-.2);
//     glVertex2f(position[0]+.2, position[1]+.2);
//     glVertex2f(position[0]-.2, position[1]+.2);
//       glEnd();
//   }   


//    glPopAttrib();
//    glPopMatrix();

// }
void RenderParticles(R3Scene *scene) {
  // glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP); 
  // glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP);

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

   glBindTexture(GL_TEXTURE_2D, scene->player->node->material->texture_index);
  for (int i = 0; i < scene->NParticles(); i++) {
    R3Particle *particle = scene->Particle(i);
    const R3Point& position = particle->position;
   glColor4f(1,1,1,1);
   // Render the front quad

   glBegin(GL_QUADS);
    glVertex2f(position[0]-.2, position[1]-.2);
    glVertex2f(position[0]+.2, position[1]-.2);
    glVertex2f(position[0]+.2, position[1]+.2);
    glVertex2f(position[0]-.2, position[1]+.2);
   glEnd();
   glBindTexture(GL_TEXTURE_2D, 0);
   // Restore enable bits wand matrix

  }   
   glPopAttrib();
   glPopMatrix();
}



