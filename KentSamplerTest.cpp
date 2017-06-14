#include <fstream>
#include <random>
#include "KentSampler.h"
using namespace std;

int main()
{
  mt19937 gen(time(nullptr));
  uniform_real_distribution<double> unif(0.0, 1.0);
  KentSampler dist(0,1,0,2);
  double x[3];
  ofstream outputFile;
  outputFile.open("KentTest.txt");
  for(int i=0;i<1000;i++)
  {
    //cout<<unif(gen)<<endl;
    dist.next(x);
    outputFile<<x[0]<<" "<<x[1]<<" "<<x[2]<<endl;
  }
  outputFile.close();
  return 0;
}
