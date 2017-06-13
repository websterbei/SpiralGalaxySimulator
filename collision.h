#ifndef cmath
#include <cmath>
#endif

#ifndef vector
#include <vector>
#endif

#ifndef cstdlib
#include <cstdlib>
#endif

#ifndef iostream
#include <iostream>
#endif

#ifndef _particle_h
#include "particle.h"
#endif

#ifndef algorithm
#include <algorithm>
#endif

#ifndef _collision_h
#define _collision_h
using namespace std;

double angle(double ax, double ay, double az, double bx, double by, double bz)
{
  return acos((ax*bx+ay*by+az*bz)/sqrt((ax*ax+ay*ay+az*az)*(bx*bx+by*by+bz*bz)));
}

double cosAngle3D(Particle a, Particle b)
{
  double dx = b.x-a.x;
  double dy = b.y-a.y;
  double dz = b.z-a.z;
  return (dx*a.vx+dy*a.vy+dz*a.vz)/sqrt((dx*dx+dy*dy+dz*dz)*(a.vx*a.vx+a.vy*a.vy+a.vz*a.vz));
}

double force(Particle a, double cosAngle)
{
  return sqrt(a.vx*a.vx+a.vy*a.vy+a.vz*a.vz)*cosAngle;
}

void forceVector(double TVector[], double NVector[],  double force, Particle a, Particle b)
{
  double dx = b.x-a.x;
  double dy = b.y-a.y;
  double dz = b.z-a.z;
  double length= sqrt(dx*dx+dy*dy+dz*dz);
  TVector[0]=force/length*dx;
  TVector[1]=force/length*dy;
  TVector[2]=force/length*dz;
  NVector[0]=a.vx-TVector[0];
  NVector[1]=a.vy-TVector[1];
  NVector[2]=a.vz-TVector[2];
}

void col(Particle* a, Particle* b, double lambda)
{
  double aTVector[3];
  double aNVector[3];
  double bTVector[3];
  double bNVector[3];
  double ang = cosAngle3D(*a, *b);
  double aforce = force(*a,ang);
  double bforce = force(*b,ang);
  forceVector(aTVector,aNVector,aforce, *a, *b);
  forceVector(bTVector,bNVector,bforce, *b, *a);
  a->vx=bTVector[0]*lambda+aNVector[0];
  a->vy=bTVector[1]*lambda+aNVector[1];
  a->vz=bTVector[2]*lambda+aNVector[2];
  b->vx=aTVector[0]*lambda+bNVector[0];
  b->vy=aTVector[1]*lambda+bNVector[1];
  b->vz=aTVector[2]*lambda+bNVector[2];
}

double distance (Particle a, Particle b)
{
  return (a.x-b.x)*(a.x-b.x)+(a.y-b.y)*(a.y-b.y)+(a.z-b.z)*(a.z-b.z);
}

int collide(vector<Particle> *allParticles, double lambda)
{
  int nCollision = 0;
  vector<Particle> & particles = *allParticles;
  sort(particles.begin(),particles.end());
  for(int i = 0; i<particles.size(); i++)
  {
    int index = i;
    int minIndex=-1;
    double min=100;
    while(true)
    {
      if(index-1<0) break;
      index--;
      if(abs(particles[index].x-particles[i].x) > 4*particles[i].r*particles[i].r) break;
      double dis = distance(particles[i],particles[index]);
      if(dis <= 4*particles[i].r*particles[i].r && dis < min && particles[i].lastCoPar != &particles[index])
      {
        min = dis;
        minIndex = index;
      }
    }

    index=i;
    while(true)
    {
      if(index+1 >= particles.size()) break;
      index++;
      if(abs(particles[index].x-particles[i].x) > 4*particles[i].r*particles[i].r) break;
      double dis = distance(particles[i],particles[index]);
      if(dis <= 4*particles[i].r*particles[i].r && dis < min && particles[i].lastCoPar != &particles[index])
      {
        min = dis;
        minIndex = index;
      }
    }
    if(minIndex==-1) continue;
    nCollision++;
    col(&particles[i], &particles[minIndex], lambda);
    particles[i].lastCoPar = &particles[minIndex];
    particles[minIndex].lastCoPar = &particles[i];
  }
  return nCollision;
}
#endif
