#include <iostream>
#include <cmath>
#include <fstream>
#include <cstdlib>
#include "particle.h"
#include <stdio.h>
#include <stdlib.h>
#include <vector>
using namespace std;

/*int cmpfunc (const void * a, const void * b){
  Particle* a1 = (Particle*)a;
  Particle* b1 = (Particle*)b;
return (a1->r -b1->r);
}
void sort(vector<Particle> *particles, int length){
qsort(&(*particles[0]),length,sizeof(Particle),cmpfunc);
}
*/


int main(int argc, char** argv){
  vector<Particle> stuff(3);
  stuff[0]=Particle();
  stuff[1]=Particle();
  stuff[2]=Particle();
  stuff[0].r = 5;
  stuff[1].r = 1;
  stuff[2].r = 3;
  sort(stuff.begin(), stuff.end());
  printf("%f\n",stuff[0].r );
  printf("%f\n",stuff[1].r );
  printf("%f\n",stuff[2].r );

  return 0;
}
