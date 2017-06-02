#include "particle.h"
#include "KentSampler.h"
#include <random>
#include <fstream>
#include <cstdlib>
#include <vector>
#include <cmath>
using namespace std;

const double mass = 1.0; //Assume mass = 1.0 for now
const double k = 2.0; //Parameter for kent distribution, assume k = 2

vector<Particle> particles;

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

int main(int argc, char** argv)
{
  int n = atoi(argv[1]); //Number of particles to be generated
  double R = atof(argv[2]); //Radius of cloud
  particles.resize(n); //Resize vector to contain the particles

  double sigma = ceil(R/sqrt(3));
  mt19937 gen; //Random Number generator
  normal_distribution<double> norm = getNorm(sigma); //Normal distribution sampler
  uniform_real_distribution<double> unif = getUnif(); //Uniform real distribution sampler

  for(vector<Particle>::size_type i=0; i<n; i++)
  {
    double tmpX = norm(gen);
    double tmpY = norm(gen);
    double tmpZ = norm(gen);
    double pot = potential(tmpX, tmpY, tmpZ); //Calculate potential energy
    double maxV = sqrt(-2*pot/mass); //Maximum velocity, assume TE = 0
    KentSampler kent = getKent(-tmpY, tmpX, 0, k);
    double tmpV[3];
    double scaleV = unif(gen);
    kent.next(tmpV);
    Particle newParticle(tmpX, tmpY, tmpZ, tmpV[0]*scaleV, tmpV[1]*scaleV, tmpV[2]*scaleV);
    particles[i] = newParticle;
  }

  ofstream outputfile;
  outputfile.open("initial.txt");

  for(vector<Particle>::size_type i=0; i<n; i++)
  {
    outputfile<<particles[i].x<<" "<<particles[i].y<<" "<<particles[i].z<<" "\
    <<particles[i].vx<<" "<<particles[i].vy<<" "<<particles[i].vz<<endl;
  }

  outputfile.close();

  return 0;
}
