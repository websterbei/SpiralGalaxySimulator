#include <cmath>
#include <iostream>

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
    double mass;
    double r; //radius of the particle
    double TE; //total energy
    bool collided;
    bool trapped;
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
      r = 0.5l;
      lastCoIndex = 0;
      lastCoStep = 0;
      index = ++indexer;
      collided = false;
      trapped = false;
    }

    void setParticle (double x, double y, double z, double vx, double vy, double vz)
    {
      this->x = x;
      this->y = y;
      this->z = z;
      this->vx = vx;
      this->vy = vy;
      this->vz = vz;
      r = 0.5l;
      lastCoIndex = 0;
      lastCoStep = 0;
      collided = false;
      trapped = false;
    }

    void setRadius (double r)
    {
      this->r = r;
    }

    void setMass (double mass)
    {
      this->mass = mass;
    }

    void setTE (double TE)
    {
      this->TE = TE;
    }

    void vCorrection (double (*potential) (double, double, double))
    {
        double PE = (*potential)(x, y, z)*mass;
        double oldKE = TE - PE;
        if(oldKE<0) std::cout<<"oldKE: "<<oldKE<<" "<<x<<" "<<y<<" "<<z<<std::endl;
        double curKE = 0.5*mass*(vx*vx+vy*vy+vz*vz);
        double cFactor = sqrt(oldKE/curKE);
        vx*=cFactor;
        vy*=cFactor;
        vz*=cFactor;
        double rad = sqrt(x*x+y*y+z*z);
        if(rad<0.5) trapped = true;
    }

    void updateTE (double (*potential) (double, double, double))
    {
      double PE = (*potential)(x, y, z)*mass;
      double curKE = 0.5*mass*(vx*vx+vy*vy+vz*vz);
      TE = PE + curKE;
      collided = false;
    }
};

int Particle::indexer = 0;

#endif
