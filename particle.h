#ifndef cmath
#include <cmath>
#endif

#ifndef _particle_h
#define _particle_h

class Particle
{
  public:
    double x;
    double y;
    double z;
    double vx;
    double vy;
    double vz;
    double dist; //Distance from center
    double r; //radius of the particle
    int lastCoIndex;
    int index;
    static int indexer;
    int lastCoStep;
    bool operator<(const Particle &rhs) const {return x < rhs.x;}

  public:
    Particle ()
    {
      x = 0.0l;
      y = 0.0l;
      z = 0.0l;
      vx = 0.0l;
      vy = 0.0l;
      vz = 0.0l;
      dist = 0.0l;
      r = 0.1l;
      lastCoIndex = 0;
      lastCoStep = 0;
      index = ++indexer;
    }

    Particle (double x, double y, double z, double vx, double vy, double vz)
    {
      this->x = x;
      this->y = y;
      this->z = z;
      this->vx = vx;
      this->vy = vy;
      this->vz = vz;
      dist = sqrt(x*x + y*y + z*z);
      r = 0.1l;
      lastCoIndex = 0;
      lastCoStep = 0;
      index = ++indexer;
    }

    Particle (double x, double y, double z)
    {
      this->x = x;
      this->y = y;
      this->z = z;
      vx = 0.0l;
      vy = 0.0l;
      vz = 0.0l;
      dist = sqrt(x*x + y*y + z*z);
      r = 0.1l;
      lastCoIndex = 0;
      lastCoStep = 0;
      index = ++indexer;
    }

    void setParticle (double x, double y, double z, double vx, double vy, double vz)
    {
      this->x = x;
      this->y = y;
      this->z = z;
      this->vx = vx;
      this->vy = vy;
      this->vz = vz;
      dist = sqrt(x*x + y*y + z*z);
      r = 0.1l;
      lastCoIndex = 0;
      lastCoStep = 0;
    }

};

int Particle::indexer = 0;

#endif
