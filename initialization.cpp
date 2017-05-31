#include "particle.h"
#include <random>
#include <fstream>
#include <cstdlib>
#include <vector>
#include <cmath>
using namespace std;

vector<Particle> particles;

normal_distribution<double> getDist(double sigma) //Require the std dev sigma
{
  normal_distribution<double> distribution(0, sigma);
  return distribution;
}

int main(int argc, char** argv)
{
  int n = atoi(argv[1]); //Number of particles to be generated
  double R = atof(argv[2]); //Radius of cloud
  particles.resize(n); //Resize vector to contain the particles

  double sigma = ceil(R/sqrt(3));
  default_random_engine gen; //Random Number generator
  normal_distribution<double> dist = getDist(sigma); //Normal distribution sampler

  for(vector<Particle>::size_type i=0; i<n; i++)
  {
    double tmpX = dist(gen);
    double tmpY = dist(gen);
    double tmpZ = dist(gen);
    Particle newParticle(tmpX, tmpY, tmpZ);
    particles[i] = newParticle;
  }

  ofstream outputfile;
  outputfile.open("initial.txt");

  for(vector<Particle>::size_type i=0; i<n; i++)
  {
    outputfile<<particles[i].x<<" "<<particles[i].y<<" "<<particles[i].z<<endl;
  }

  outputfile.close();

  return 0;
}
