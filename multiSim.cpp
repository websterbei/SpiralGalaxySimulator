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

const double r = 0.1; //particle radius default to 1
const double k = 1000.0; //Parameter for kent distribution, assume k = 2
const double p = 1000;
int n = 100; //Number of particles
int R = 10.0; //Radius of gas cloud
int nPic = 10; //Number of figs generated
double lambda = 1.0; //Coefficient of restitution
double tEnd = 1000;

vector<Particle> particles; //Array of particles

double potential(double x, double y, double z) //Potential function
{
  return -exp(-(x*x+y*y+z*z));
}

normal_distribution<double> getNorm(double sigma) //Require the std dev sigma
{
  normal_distribution<double> distribution(0, sigma);
  return distribution;
}

uniform_real_distribution<double> getUnif() //Obtain a uniform real distribution
{
  uniform_real_distribution<double> distribution(0.0, 1.0);
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
    <<particles[i].vx<<" "<<particles[i].vy<<" "<<particles[i].vz<<endl;
  }
  outputfile.close();
}

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

  particles.resize(n); //Resize vector to contain the particles
  double sigma = ceil(R/sqrt(3));
  mt19937 gen(time(nullptr)); //Random Number generator
  normal_distribution<double> norm = getNorm(sigma); //Normal distribution sampler
  uniform_real_distribution<double> unif = getUnif(); //Uniform real distribution sampler

  for(vector<Particle>::size_type i=1; i<n; i++)
  {
    double tmpX = norm(gen);
    double tmpY = norm(gen);
    double tmpZ = norm(gen);
    double pot = potential(tmpX, tmpY, tmpZ); //Calculate potential energy
    double maxV = sqrt(-2*pot/r/r/r/p); //Maximum velocity, assume TE = 0
    KentSampler kent = getKent(-tmpY, tmpX, 0, k);
    double tmpV[3];
    double scaleV = unif(gen);
    //cout<<maxV<<endl;
    kent.next(tmpV);
    Particle newParticle(tmpX, tmpY, tmpZ, maxV*tmpV[0]*scaleV, maxV*tmpV[1]*scaleV, maxV*tmpV[2]*scaleV);
    particles[i] = newParticle;
  }

  particles[0] = Particle(1, 0, 0, 0, 0.8478, 0);
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
    if(stepCounter%avgStepSep==0)
    {
      string fname = "./pics/multiSim" + to_string(stepCounter/avgStepSep);
      printToFile(fname);
      /*double z [3] = {0,0,0};
      for(int i = 0; i < particles.size(); i ++){
        z[0] += (particles[i].y*particles[i].vz)-(particles[i].vy*particles[i].z);
        z[1] += (particles[i].vx*particles[i].z)-(particles[i].vz*particles[i].x);
        z[2] += (particles[i].x*particles[i].vy)-(particles[i].vx*particles[i].y);
      }
      cout<<z[0]<<" "<<z[1]<<" "<<z[2]<<endl;
      */
    }

    int i = 1;
    test<<particles[i].x<<" "<<particles[i].y<<" "<<particles[i].z<<" "<<particles[i].vx<<" "<<particles[i].vy<<" "<<particles[i].vz<<endl;
    nCollision += collide(&particles, lambda);

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
