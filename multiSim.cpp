#include <random>
#include <fstream>
#include <cstdlib>
#include <vector>
#include <cmath>
#include <string>
#include <iostream>
#include <thread>
#include <cstdio>
using namespace std;

//Constant declaration
//density plot potation curves potential well in x-y plane
#define PI 3.141592653589793238l
const double rmax = 1000000;
int nsteps =  1 + (int)(rmax/1);//nsteps = 1 +(int)(rmax/dr);
double dr = rmax/(nsteps-1);
const double A0 = 1, A2 = -0.5;
const double waveLength0 = 2000, waveLength2 = 1990;
const double k0 = 2 * PI/waveLength0, k2 = 2 * PI/waveLength2;
const double patternPeriod = 50000000;
const double W0 = patternPeriod/8/PI * (k2*k2-k0*k0) - 2 * PI/patternPeriod;
const double W2 = patternPeriod/8/PI * (k2*k2-k0*k0) + 2 * PI/patternPeriod;
const double dw = W2 - W0;
//const double dw = 0;
const double Omega0Squared = (2*PI/patternPeriod)*(2*PI/patternPeriod) + (patternPeriod/8/PI)*(patternPeriod/8/PI)*(k2*k2-k0*k0)*(k2*k2-k0*k0)- (k2*k2+k0*k0)/2;
const double Omega0 = sqrt(Omega0Squared);
const double mu0withSign = 8.7e-13;
const double mu0 = abs(mu0withSign);
const double diskscalelength = 15000;
const double r = 5000; //particle radius default to 1
const double k = 1000.0; //Parameter for kent distribution, assume k = 2
//const double density = 0.001;
int n = 5000; //Number of particles
int nCollision = 0;
int R = 500000; //Radius of gas cloud
int nPic = 10; //Number of figs generated
double lambda = 1.0; //Coefficient of restitution
double tEnd = 100000000;
double stepSize = 10000;
double mass = 1;
double rvalues [1000002];
double f0[1000002];
double f2[1000002];
double u0[1000002];
double u2[1000002];
double u4[1000002];
double uu2[1000002];
double w0[1000002];
double w0p[1000002];
double w2[1000002];
double w2p[1000002];
double w4[1000002];
double w4p[1000002];
double ww2[1000002];
double ww2p[1000002];
double t = 0.0;

#include "RKOdeSolver.h"
#include "particle.h"
#include "KentSampler.h"
#include "collision.h"


//Array of particles
vector<Particle> particles;

//Function declaration
double potential(double x, double y, double z, double t);
normal_distribution<double> getNorm(double sigma);
uniform_real_distribution<double> getUnif();
KentSampler getKent(double x, double y, double z, double k);
void printToFile(string fname);
void RK4Thread(int start, int end, double stepSize);
void potentialGenerator();
void printPotential(double t);
void printSingleTest(ofstream *file);

//Simulation
int main(int argc, char** argv)
{
  if(argc<=1)
  {
    //Initialization of gas cloud
    cout<<"Number of particles: ";
    cin>>n;
    cout<<"Radius of gas cloud: ";
    cin>>R;
    cout<<"Coefficient of restitution: ";
    cin>>lambda;
    cout<<"Number of snapshots: ";
    cin>>nPic;

    particles = vector<Particle>(n); //Resize vector to contain the particles
    double sigma = ceil(R/sqrt(3));
    mt19937 gen(time(nullptr)); //Random Number generator
    normal_distribution<double> norm = getNorm(sigma); //Normal distribution sampler
    uniform_real_distribution<double> unif = getUnif(); //Uniform real distribution sampler
    //Generate Potential Array
    potentialGenerator();
    //printPotential(0);

    for(vector<Particle>::size_type i=0; i<n; i++)
    {
      double tmpX = norm(gen);
      double tmpY = norm(gen);
      double tmpZ = norm(gen);
      double pot = potential(tmpX, tmpY, tmpZ, 0); //Calculate potential energy
      double maxV = sqrt(-2*pot); //Maximum velocity, assume TE = 0
      KentSampler kent = getKent(-tmpY, tmpX, 0, k);
      double tmpV[3];
      double scaleV = unif(gen);
      kent.next(tmpV);
      //cout<<tmpX<<" "<<tmpY<<" "<<tmpZ<<" "<<pot<<endl;
      //r = sqrt(tmpX*tmpX+tmpY*tmpY+tmpZ*tmpZ);
      //double acc[3];
      //accel(acc, tmpX, tmpY, tmpZ, 0);
      particles[i].setParticle(tmpX, tmpY, tmpZ, maxV*tmpV[0]*scaleV, maxV*tmpV[1]*scaleV, maxV*tmpV[2]*scaleV);
      particles[i].setRadius(r);
      particles[i].setMass(mass);
      //double tmpKE = 0.5*mass*(particles[i].vx*particles[i].vx+particles[i].vy*particles[i].vy+particles[i].vz*particles[i].vz);
      //double tmpPE = mass*pot;
      //particles[i].setTE(tmpKE+tmpPE);
    }
  }
  else if(strcmp(argv[1], "resume")==0)
  {
    ifstream configFile;
    configFile.open(argv[2]);
    configFile>>n>>lambda>>t>>nCollision>>nPic;
    particles = vector<Particle>(n);
    uniform_real_distribution<double> unif = getUnif();
    mt19937 gen(time(nullptr));
    double tmp;
    for(int i=0;i<n;i++)
    {
      configFile>>particles[i].x>>particles[i].y>>particles[i].z>>\
      particles[i].vx>>particles[i].vy>>particles[i].vz>>tmp>>particles[i].index;
      //if(i%2==0) particles[i].vx+=particles[i].vx/2;
      //else particles[i].vx-=particles[i].vx/2;
      //if(i%2==0) particles[i].vy+=particles[i].vy/2;
      //else particles[i].vy-=particles[i].vy/2;

      particles[i].setRadius(r);
      particles[i].setMass(mass);
    }
    tEnd=tEnd+t;
    potentialGenerator();
  }

  //Iterative update of particle location
  int stepCounter = 0;
  int avgStepSep = (int)(tEnd/stepSize/nPic);
  int nThread = 4;
  int oldNCollision = 0;
  double oldTE =0;
  thread thrd[nThread];

  // ofstream energyDecrease;
  // ofstream angMomCons;
  // ofstream denergyChange;
  // ofstream diskness;
  // energyDecrease.open("energyDecrease.txt");
  // angMomCons.open("angMomCons.txt");
  // denergyChange.open("denergyChange.txt");
  // diskness.open("diskness.txt");

  int progressCounter = 0;

  ofstream outputfile;
  outputfile.open("output.txt");

  while(t<tEnd)
  {
    printSingleTest(&outputfile);

    //Print
    if(stepCounter%avgStepSep==0)
    {
      string fname = "./pics/multiSim" + to_string(stepCounter/avgStepSep);
      printToFile(fname);
    }

    //Collision
    //nCollision += collide(&particles, lambda);

    for(int i=0;i<n;i++)
    {
      solve(&particles[i], t);
    }

    t+=stepSize;
    stepCounter++;
  }
  outputfile.close();

  // energyDecrease.close();
  // angMomCons.close();
  // denergyChange.close();
  // diskness.close();
  cout<<nCollision<<endl;
  printToFile("./pics/config");

  return 0;
}

void printSingleTest(ofstream *file)
{
  for(int i=0;i<n;++i)
  {
    if(particles[i].index == 1)
      *file<<particles[i].x<<" "<<particles[i].y<<" "<<particles[i].z<<" "<<particles[i].vx<<" "<<particles[i].vy<<" "<<particles[i].vz<<endl;
  }
}

double potential(double x, double y, double z, double t)
{
  double sr = sqrt(x*x+y*y+z*z);
  double w0Value, w2Value, w4Value, ww2Value;
  if(sr <= rmax)
  {
    int ival = round(sr/dr) + 1;
    w0Value = w0[ival];
    w2Value = w2[ival];
    w4Value = w4[ival];
    ww2Value = ww2[ival];
  }
  else
  {
    w0Value = w0[nsteps] * rmax/sr;
    w2Value = w2[nsteps] * pow(rmax/sr, 5);
    w4Value = w4[nsteps] * pow(rmax/sr, 9);
    ww2Value = ww2[nsteps] * pow(rmax/sr, 5);
  }

  double temp = w0Value + w2Value*(3*z*z-sr*sr) + w4Value*(35*z*z*z*z-30*sr*sr*z*z+3*sr*sr*sr*sr);
  double temp2 = temp + ww2Value*(cos(dw*t) * (x*x - y*y) + sin(dw*t) * (2*x*y));
  return 4*PI*mu0*temp2;
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

  outputfile<<n<<" "<<lambda<<" "<<t<<" "<<nCollision<<" "<<nPic<<endl;

  for(vector<Particle>::size_type i=0; i<n; i++)
  {
    outputfile<<particles[i].x<<" "<<particles[i].y<<" "<<particles[i].z<<" "\
    <<particles[i].vx<<" "<<particles[i].vy<<" "<<particles[i].vz<<" "<<particles[i].TE<<" "\
    <<particles[i].index<<endl;
  }
  outputfile.close();
}

void printPotential(double t)
{
  ofstream outputfile;
  string fname = "./potential/potentialat"+to_string(int(t));
  outputfile.open(fname);
  for(int i=-rmax;i<rmax;i+=1000)
  {
    for(int j=-rmax;j<rmax;j+=1000)
    {
      outputfile<<potential(i,j,0,t)<<" ";
    }
    outputfile<<endl;
  }
  outputfile.close();
}

double min(double a[])
{
  double temp = a[1];
  for(int i=2;i<nsteps+1;i++)
  {
    if(a[i]<temp) temp = a[i];
  }
  return temp;
}

double max(double a[])
{
  double temp = a[1];
  for(int i=2;i<nsteps+1;i++)
  {
    if(a[i]>temp) temp = a[i];
  }
  return temp;
}

void potentialGenerator()
{
  //generate rvalues
  for(int i = 2; i<nsteps+1;i++) rvalues[i]=rvalues[i-1]+dr;

  //generate f
  f0[1]=1;
  double f0p = 0;
  double f0pp = -k0*k0*f0[1]/3;
  f2[1]=1;
  double f2p = 0;
  double f2pp = -k2*k2*f2[1]/7;

  for(int i = 2; i<nsteps+1;i++)
  {
    double sr = (i-1)*dr;
    f0[i]=f0[i-1]+f0p*dr + f0pp*dr*dr/2;
    f0p = f0p + f0pp*dr;
    f0pp = -k0*k0 * f0[i] - 2 * f0p / sr;

    f2[i]=f2[i-1]+f2p*dr + f2pp*dr*dr/2;
    f2p = f2p + f2pp*dr;
    f2pp = -k2*k2 * f2[i] - 6 * f2p / sr;
  }

  double *f0temp = (double *)malloc(sizeof(double)*1000002);
  double *f2temp = (double *)malloc(sizeof(double)*1000002);

  for(int i=1;i<nsteps+1;i++)
  {
    // cout<<i<<endl;
    f0temp[i]=f0[i]*rvalues[i];
    f2temp[i]=f2[i]*rvalues[i]*rvalues[i]*rvalues[i];
  }
  double min0 = min(f0temp);
  double min2 = min(f2temp);

  free(f0temp); free(f2temp);

  for(int i=1;i<nsteps+1;i++)
  {
    f0[i]/=-min0;
    f2[i]/=-min2;
  }

  double fac = max(f0);
  for(int i=1;i<nsteps+1;i++)
  {
    f0[i]/=fac;
    f2[i]/=fac;
  }

  //compute u
  for(int i = 1; i<nsteps+1;i++)
  {
    double sr = (i-1)*dr;
    u0[i]= A0 * A0 * f0[i] * f0[i]+ 56.0/105 * A2 * A2 * sr * sr * sr * sr * f2[i] * f2[i];
    u2[i]= -40.0/105 * A2 * A2 * sr * sr * f2[i] * f2[i];
    u4[i]= 3.0/105 * A2 * A2 * f2[i] * f2[i];
    uu2[i]= 2 * A0 * A2 * f0[i] * f2[i];
  }

  //comupte w
  double w0pp = u0[1]/3;
  double w2pp = u2[1]/7;
  double w4pp = u4[1]/11;
  double ww2pp = uu2[1]/7;
  for (int i=2;i<nsteps+1;i++)
  {
    double sr = (i-1)*dr;

    w0[i] = w0[i-1] + w0p[i-1] * dr + w0pp * dr * dr / 2;
    w0p[i] = w0p[i-1] + w0pp * dr;
    w0pp = u0[i] - 2 * w0p[i] / sr;

    w2[i] = w2[i-1] + w2p[i-1] * dr + w2pp * dr * dr / 2;
    w2p[i] = w2p[i-1] + w2pp * dr;
    w2pp = u2[i] - 6 * w2p[i] / sr;

    w4[i] = w4[i-1] + w4p[i-1] * dr + w4pp * dr * dr / 2;
    w4p[i] = w4p[i-1] + w4pp * dr;
    w4pp = u4[i] - 10 * w4p[i] / sr;

    ww2[i] = ww2[i-1] + ww2p[i-1] * dr + ww2pp * dr * dr / 2;
    ww2p[i] = ww2p[i-1] + ww2pp * dr;
    ww2pp = uu2[i] - 6 * ww2p[i] / sr;
  }

  double W0shift = - rmax / 1 * w0p[nsteps] - w0[nsteps];
  double W2shift = - rmax / 5 * w2p[nsteps] - w2[nsteps];
  double W4shift = - rmax / 9 * w4p[nsteps] - w4[nsteps];
  double WW2shift = - rmax / 5 * ww2p[nsteps] - ww2[nsteps];

  for(int i=1;i<nsteps+1;i++)
  {
    w0[i] = w0[i] + W0shift;
    w2[i] = w2[i] + W2shift;
    w4[i] = w4[i] + W4shift;
    ww2[i] = ww2[i] + WW2shift;
  }
}

void accel(double d2[], double x, double y, double z, double t)
{
    double r = sqrt(x*x+y*y+z*z);
    if(r == 0)
    {
      d2[0] = 0; d2[1] = 0; d2[2] = 0;
      return;
    }
    double w0Value, w2Value, w4Value, ww2Value, w0pValue, w2pValue, w4pValue, ww2pValue;
    if(r <= rmax)
    {
      int ival = round(r/dr) + 1;
      w0Value = w0[ival];
      w2Value = w2[ival];
      w4Value = w4[ival];
      ww2Value = ww2[ival];

      w0pValue = w0p[ival];
      w2pValue = w2p[ival];
      w4pValue = w4p[ival];
      ww2pValue = ww2p[ival];
    }
    else
    {
      w0Value = w0[nsteps] * rmax/r;
      w2Value = w2[nsteps] * pow(rmax/r, 5);
      w4Value = w4[nsteps] * pow(rmax/r, 9);
      ww2Value = ww2[nsteps] * pow(rmax/r, 5);

      w0pValue = -w0Value / r;
      w2pValue = -w2Value / r * 5;
      w4pValue = -w4Value / r * 9;
      ww2pValue = -ww2Value / r * 5;
    }
    double term1 = w0pValue + w2pValue * (3*z*z - r*r) + w4pValue * (35*z*z*z*z - 30*r*r*z*z + 3*r*r*r*r);
    double term2 = ww2pValue * (cos(dw*t) * (x*x - y*y) + sin(dw*t) * (2*x*y));
    d2[0] = (term1+term2)*x/r-w2Value*2*x+w4Value*12*x*(r*r-5*z*z)+cos(dw*t)*ww2Value*2*x+sin(dw*t)*ww2Value*2*y;
    d2[1] = (term1+term2)*y/r-w2Value*2*y+w4Value*12*y*(r*r-5*z*z)-cos(dw*t)*ww2Value*2*y+sin(dw*t)*ww2Value*2*x;
    d2[2] = (term1+term2)*z/r+w2Value*4*z+w4Value*16*z*(5*z*z-3*r*r);
    d2[0]*=-4*PI*mu0;
    d2[1]*=-4*PI*mu0;
    d2[2]*=-4*PI*mu0;
}

void RK4Thread(int start, int end, double t)
{
  for(int i=start;i<=end;++i) solve(&particles[i], t);
}
