#include "particle.h"
#include "KentSampler.h"
#include "RKOdeSolver.h"
#include "collision.h"
#include <random>
#include <fstream>
#include <cstdlib>
#include <vector>
#include <cmath>
#include <string>
#include <iostream>
using namespace std;

//Constant declaration
const double r = 0.5; //particle radius default to 1
const double k = 1000.0; //Parameter for kent distribution, assume k = 2
const double density = 0.1;
int n = 100; //Number of particles
int R = 10.0; //Radius of gas cloud
int nPic = 10; //Number of figs generated
double lambda = 1.0; //Coefficient of restitution
double tEnd = 1000;
double mass = density*r*r*r;

//Array of particles
vector<Particle> particles;

//Function declaration
double potential(double x, double y, double z);
normal_distribution<double> getNorm(double sigma);
uniform_real_distribution<double> getUnif();
KentSampler getKent(double x, double y, double z, double k);
void printToFile(string fname);
void singleTest(ofstream *file);
void totAngMom(double z[]);
void optimalTest();
void avgKEPE(double *Kinetic, double *Potential);

//Simulation
int main()
{
  //Initialization of gas cloud
  cout<<"Number of particles: ";
  cin>>n;
  cout<<"Radius of gas cloud: ";
  cin>>R;
  cout<<"Coefficient of restitution: ";
  cin>>lambda;
  cout<<"Ending time: ";
  cin>>tEnd;
  cout<<"Number of snapshots: ";
  cin>>nPic;

  particles = vector<Particle>(n); //Resize vector to contain the particles
  double sigma = ceil(R/sqrt(3));
  mt19937 gen(time(nullptr)); //Random Number generator
  normal_distribution<double> norm = getNorm(sigma); //Normal distribution sampler
  uniform_real_distribution<double> unif = getUnif(); //Uniform real distribution sampler

  for(vector<Particle>::size_type i=0; i<n-1; i++)
  {
    double tmpX = norm(gen)*2;
    double tmpY = norm(gen);
    double tmpZ = norm(gen);
    double pot = potential(tmpX, tmpY, tmpZ); //Calculate potential energy
    double maxV = sqrt(-2*pot/mass); //Maximum velocity, assume TE = 0
    KentSampler kent = getKent(-tmpY, tmpX, 0, k);
    double tmpV[3];
    double scaleV = unif(gen);
    kent.next(tmpV);
    particles[i].setParticle(tmpX, tmpY, tmpZ, maxV*tmpV[0]*scaleV, maxV*tmpV[1]*scaleV, maxV*tmpV[2]*scaleV);
    particles[i].setRadius(r);
    particles[i].setMass(mass);
  }

  particles[n-1].setParticle(1, 0, 0, 0, 14.14, 0);
  particles[n-1].setRadius(r);
  particles[n-1].setMass(mass);

  //Iterative update of particle location
  double t = 0.0;
  double stepSize = 0.01;
  int stepCounter = 0;
  int avgStepSep = (int)(tEnd/stepSize/nPic);
  int nCollision = 0;
  ofstream test;
  test.open("output.txt");

  while(t<tEnd)
  {
    //Single Particle
    //singleTest(&test);
    //Print
    if(stepCounter%avgStepSep==0)
    {
      string fname = "./pics/multiSim" + to_string(stepCounter/avgStepSep);
      printToFile(fname);
      double z[3];
      totAngMom(z);
      cout<<z[0]<<" "<<z[1]<<" "<<z[2]<<endl;
      optimalTest();
      double KE,PE;
      avgKEPE(&KE, &PE);
      cout<<KE*2<<" "<<PE<<endl;
    }
    //Collision
    nCollision += collide(&particles, lambda);
    //RK4
    for(vector<Particle>::size_type i=0; i<n; i++)
    {
      solve(&particles[i], stepSize);
    }

    t+=stepSize;
    stepCounter++;
  }
  test.close();
  cout<<nCollision<<endl;
  return 0;
}

/*double potential(double x, double y, double z) //Potential function
{
  return -exp(-(x*x+y*y+z*z));
}*/

/*double potential(double x, double y, double z) //Potential function
{
  return a*(x*x+y*y+z*z)-1000;
}*/

double potential(double x, double y, double z) //Potential function
{
  return -1.0/sqrt(x*x+y*y+z*z);
}

normal_distribution<double> getNorm(double sigma) //Require the std dev sigma
{
  normal_distribution<double> distribution(0, sigma);
  return distribution;
}

uniform_real_distribution<double> getUnif() //Obtain a uniform real distribution
{
  uniform_real_distribution<double> distribution(0, 1.0);
  return distribution;
}

KentSampler getKent(double x, double y, double z, double k) //k controls the overall angular momentum
{
  KentSampler distribution(x, y, z, k); //2 is a suitable number for k
  return distribution;
}

void printToFile(string fname)
{
  ofstream outputfile;
  outputfile.open(fname);

  for(vector<Particle>::size_type i=0; i<n; i++)
  {
    outputfile<<particles[i].x<<" "<<particles[i].y<<" "<<particles[i].z<<" "\
    <<particles[i].vx<<" "<<particles[i].vy<<" "<<particles[i].vz<<" "<<particles[i].index<<endl;
  }
  outputfile.close();
}

void totAngMom(double z[])
{
  z[0] = 0; z[1] = 0; z[2] = 0;
  for(int i = 0; i < particles.size(); i ++)
  {
    z[0] += (particles[i].y*particles[i].vz)-(particles[i].vy*particles[i].z);
    z[1] += (particles[i].vx*particles[i].z)-(particles[i].vz*particles[i].x);
    z[2] += (particles[i].x*particles[i].vy)-(particles[i].vx*particles[i].y);
  }
}

void singleTest(ofstream *file)
{
  ofstream& test = *file;
  for(vector<Particle>::size_type i=0; i<n; i++)
  {
    if(particles[i].index == n-2)
    {
      Particle *ptc = &particles[i];
      test<<ptc->x<<" "<<ptc->y<<" "<<ptc->z<<" "<<ptc->vx<<" "<<ptc->vy<<" "<<ptc->vz<<endl;
      break;
    }
  }
}

void optimalTest() //For -1/r well
{
  double optimFactor = 0.0;
  int counter = 0;
    for(int i = 0; i < particles.size(); i ++)
    {
      double rad = sqrt(particles[i].x*particles[i].x + particles[i].y*particles[i].y + particles[i].z*particles[i].z);
      double v2 = particles[i].vx*particles[i].vx + particles[i].vy*particles[i].vy + particles[i].vz*particles[i].vz;
      double totE = potential(particles[i].x, particles[i].y, particles[i].z) + 0.5*mass*v2;
      double maxR = -0.5/totE;
      if(maxR>=0) counter++;
      else continue;
      //cout<<maxR<<endl;
      double lZ = particles[i].x*particles[i].vy-particles[i].vx*particles[i].y;
      //cout<<totE-potential(maxR, 0, 0)<<endl;
      optimFactor+=lZ/maxR/sqrt(2*(totE-potential(maxR, 0, 0)/mass));
    }
    cout<<optimFactor/counter<<endl;
}

void avgKEPE(double *Kinetic, double *Potential) //Average kinetic energy
{
  double &KE = *Kinetic;
  double &PE = *Potential;
  for(int i=0;i<particles.size();i++)
  {
    double tmpKE = 0.5*mass*(particles[i].vx*particles[i].vx+particles[i].vy*particles[i].vy+particles[i].vz*particles[i].vz);
    double tmpPE = potential(particles[i].x, particles[i].y, particles[i].z);
    if(tmpKE+tmpPE<0)
    {
      KE+=tmpKE;
      PE+=tmpPE;
    }
  }
}
