#include <iostream>
#include <cmath>
#include <fstream>
#include <cstdlib>
using namespace std;


int main(int argc, char** argv){
    double lambda = 0.5;
    double step = 0.001;
    double radius = 0.1;
    double Vx [2];
    double Vy [2];
    double Vz [2];
    double x  [2];
    double y  [2];
    double z  [2];
    x[0] = atof(argv[1]);
    y[0] = atof(argv[2]);
    z[0] = atof(argv[3]);
    x[1] = atof(argv[4]);
    y[1] = atof(argv[5]);
    z[1] = atof(argv[6]);
    Vx[0] = atof(argv[7]);
    Vy[0] = atof(argv[8]);
    Vz[0] = atof(argv[9]);
    Vx[1] = atof(argv[10]);
    Vy[1] = atof(argv[11]);
    Vz[1] = atof(argv[12]);
    double t = 0;
    bool a = true;
    int counter = 0;
    ofstream outputFile1;
    outputFile1.open("output1.txt");
    ofstream outputFile2;
    outputFile2.open("output2.txt");
    while(t<20){
      counter++;
      if(counter>=1000){counter = 0;
      outputFile1<<x[0]<<" "<<y[0]<<" "<<z[0]<<endl;
      outputFile1<<x[1]<<" "<<y[1]<<" "<<z[1]<<endl;}
      double r = (x[1]-x[0])*(x[1]-x[0])+(y[1]-y[0])*(y[1]-y[0])+(z[1]-z[0])*(z[1]-z[0]);
      if(r<=4*radius*radius&&a){
        double V1x = (Vx[0]+Vx[1]+lambda*(Vx[1]-Vx[0]))/2;
        double V2x = (Vx[0]+Vx[1]+lambda*(Vx[0]-Vx[1]))/2;
        Vx[0]=V1x;
        Vx[1]=V2x;
        double V1y = (Vy[0]+Vy[1]+lambda*(Vy[1]-Vy[0]))/2;
        double V2y = (Vy[0]+Vy[1]+lambda*(Vy[0]-Vy[1]))/2;
        Vy[0]=V1y;
        Vy[1]=V2y;
        double V1z = (Vz[0]+Vz[1]+lambda*(Vz[1]-Vz[0]))/2;
        double V2z = (Vz[0]+Vz[1]+lambda*(Vz[0]-Vz[1]))/2;
        Vz[0]=V1z;
        Vz[1]=V2z;
        a = false;
      }
      x[0]+=Vx[0]*step;
      x[1]+=Vx[1]*step;
      y[0]+=Vy[0]*step;
      y[1]+=Vy[1]*step;
      z[0]+=Vz[0]*step;
      z[1]+=Vz[1]*step;
      t+=step;

    }
    outputFile1.close();
    outputFile2.close();

  return 0;
}
