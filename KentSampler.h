#ifndef cmath
#include <cmath>
#endif

#ifndef random
#include <random>
#endif

#ifndef ctime
#include <ctime>
#endif

#ifndef _KentSampler_h
#define _KentSampler_h
#define PI 3.141592653589793238l
class KentSampler
{
    double u[3];
    double k;
    double Cp;
    std::mt19937 gen;
    std::uniform_real_distribution<double> xdist;
    std::uniform_real_distribution<double> ydist;

  public:
    KentSampler (double ux, double uy, double uz, double k)
    {
      u[0] = ux;
      u[1] = uy;
      u[2] = uz;
      normalize(u);
      this->k = k;
      Cp = k/2/PI/(1-exp(-2*k));
      gen = std::mt19937(time(nullptr));
      xdist = std::uniform_real_distribution<double>(-1.0, 1.0);
      ydist = std::uniform_real_distribution<double>(0, 1.0);
    }

    double VMF(double x[])
    {
      return Cp*exp(k*(x[0]*u[0]+x[1]*u[1]+x[2]*u[2]-1));
    }

    void normalize(double x[])
    {
      double length = sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2]);
      x[0]/=length;
      x[1]/=length;
      x[2]/=length;
    }

    void next(double x[])
    {
      bool success = false;
      while(!success)
      {
        x[0] = xdist(gen);
        x[1] = xdist(gen);
        x[2] = xdist(gen);
        double y = ydist(gen);
        normalize(x);
        double fxMax = Cp;
        double fx = VMF(x);
        if(y<fx/fxMax) success = true;
      }
    }
};
#endif
