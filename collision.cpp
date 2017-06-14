#include <cmath>
#include <vector>
#include <cstdlib>
#include <iostream>
#include <algorithm>
#include <fstream>
#include "particle.h"
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

void col(Particle* a, Particle* b)
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
  a->vx=bTVector[0]+aNVector[0];
  a->vy=bTVector[1]+aNVector[1];
  a->vz=bTVector[2]+aNVector[2];
  b->vx=aTVector[0]+bNVector[0];
  b->vy=aTVector[1]+bNVector[1];
  b->vz=aTVector[2]+bNVector[2];
}

double distance (Particle a, Particle b)
{
  return (a.x-b.x)*(a.x-b.x)+(a.y-b.y)*(a.y-b.y)+(a.z-b.z)*(a.z-b.z);
}

void collide(vector<Particle> *allParticles)
{
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
    col(&particles[i], &particles[minIndex]);
    particles[i].lastCoPar = &particles[minIndex];
    particles[minIndex].lastCoPar = &particles[i];
  }

}


int main()
{ 
  ofstream output;
  output.open("output.txt");
  vector<Particle> a(3);
  a[0]=Particle(0,0,0,3,0,0);
  a[1]=Particle(5,5,0,-2,-5,0);
  a[2]=Particle(0,10,0,-3,0,0);
  double step = 0.01;
  double t = 0.0;
  while(t<3) 
  {
    collide(&a);
    output<<a[0].x<<" "<<a[0].y<<" "<<a[0].z<<" "<<a[0].vx<<" "<<a[0].vy<<" "<<a[0].vz<<endl;
    output<<a[1].x<<" "<<a[1].y<<" "<<a[1].z<<" "<<a[1].vx<<" "<<a[1].vy<<" "<<a[1].vz<<endl;
    output<<a[2].x<<" "<<a[2].y<<" "<<a[2].z<<" "<<a[2].vx<<" "<<a[2].vy<<" "<<a[2].vz<<endl;
    a[0].x += a[0].vx*step; a[0].y += a[0].vy*step; a[0].z += a[0].vz*step;
    a[1].x += a[1].vx*step; a[1].y += a[1].vy*step; a[1].z += a[1].vz*step;
    a[2].x += a[2].vx*step; a[2].y += a[2].vy*step; a[2].z += a[2].vz*step;
    cout<<a[0].vx+a[1].vx+a[2].vx<<" "<<a[0].vy+a[1].vy+a[2].vy<<" "<<a[0].vz+a[1].vz+a[2].vz<<endl;
    t+=step;
  }
  output.close();
  return 0;
}
