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

int main()
{
  Particle particle(1, 0, 0, 0, 0.8577, 0);
  double t = 0.0;
  double stepSize = 0.001;
  int stepCounter = 0;
  ofstream test;
  test.open("output.txt");
  while(t<10)
  {
    /*if(stepCounter%avgStepSep==0)
    {
      string fname = "multiSim" + to_string(stepCounter/avgStepSep);
      printToFile(fname);
    }*/
    test<<particle.x<<" "<<particle.y<<" "<<particle.z<<endl;

    solve(&particle, stepSize);
    t+=stepSize;
    stepCounter++;
  }
  test.close();
  return 0;
}
