#include <cmath>
#include <vector>
#include <cstdlib>
#include <iostream>
#include <algorithm>
#include "particle.h"

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
  double angA = cosAngle3D(*a, *b);
  double angB = cosAngle3D(*b, *a);
  double aforce = force(*a, angA);
  double bforce = force(*b, angB);
  forceVector(aTVector, aNVector, aforce, *a, *b);
  forceVector(bTVector, bNVector, bforce, *b, *a);

  a->vx=0.5*(aTVector[0]+bTVector[0]+lambda*(bTVector[0]-aTVector[0]))+aNVector[0];
  a->vy=0.5*(aTVector[1]+bTVector[1]+lambda*(bTVector[1]-aTVector[1]))+aNVector[1];
  a->vz=0.5*(aTVector[2]+bTVector[2]+lambda*(bTVector[2]-aTVector[2]))+aNVector[2];
  b->vx=0.5*(aTVector[0]+bTVector[0]+lambda*(aTVector[0]-bTVector[0]))+bNVector[0];
  b->vy=0.5*(aTVector[1]+bTVector[1]+lambda*(aTVector[1]-bTVector[1]))+bNVector[1];
  b->vz=0.5*(aTVector[2]+bTVector[2]+lambda*(aTVector[2]-bTVector[2]))+bNVector[2];

}

double distance2 (Particle a, Particle b)
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
      if(abs(particles[index].x-particles[i].x) > 2*particles[i].r) break;
      double dis = distance2(particles[i], particles[index]);
      if(dis <= 4*particles[i].r*particles[i].r && dis < min &&
        !((particles[i].lastCoIndex == particles[index].index)&&(particles[index].lastCoIndex == particles[i].index)))
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
      if(abs(particles[index].x-particles[i].x) > 2*particles[i].r) break;
      double dis = distance2(particles[i], particles[index]);
      if(dis <= 4*particles[i].r*particles[i].r && dis < min &&
        !((particles[i].lastCoIndex == particles[index].index)&&(particles[index].lastCoIndex == particles[i].index)))
      {
        min = dis;
        minIndex = index;
      }
    }

    if(minIndex==-1) continue;
    nCollision++;
    col(&particles[i], &particles[minIndex], lambda);
    particles[i].lastCoIndex = particles[minIndex].index;
    particles[minIndex].lastCoIndex = particles[i].index;
    particles[i].collided = true;
    particles[minIndex].collided = true;
  }
  return nCollision;
}
#endif
