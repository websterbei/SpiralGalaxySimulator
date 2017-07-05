#include "particle.h"
#include <iostream>
#include <cmath>
#include <cstdlib>

#ifndef _RKOdeSolver_h
#define _RKOdeSolver_h

const double a = 1.0;

void deri(double d2[], double d[], double x, double y, double z, double vx, double vy, double vz, double t)
{
    //cout<<vx<<endl;
    d[0] = vx;
    d[1] = vy;
    d[2] = vz;

    double sr = sqrt(x*x+y*y+z*z);
    if(sr == 0)
    {
      d[0] = 0; d[1] = 0; d[2] = 0;
      d2[0] = 0; d2[1] = 0; d2[2] = 0;
      return;
    }
    double w0Value, w2Value, w4Value, ww2Value, w0pValue, w2pValue, w4pValue, ww2pValue;
    if(sr <= rmax)
    {
      int ival = round(sr/dr) + 1;
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
      w0Value = w0[nsteps] * rmax/sr;
      w2Value = w2[nsteps] * pow(rmax/sr, 5);
      w4Value = w4[nsteps] * pow(rmax/sr, 9);
      ww2Value = ww2[nsteps] * pow(rmax/sr, 5);

      w0pValue = -w0Value / sr;
      w2pValue = -w2Value / sr * 5;
      w4pValue = -w4Value / sr * 9;
      ww2pValue = -ww2Value / sr * 5;
    }
    double term1 = w0pValue + w2pValue * (3*z*z - sr*sr) + w4pValue * (35*z*z*z*z - 30*sr*sr*z*z + 3*sr*sr*sr*sr);
    double term2 = ww2pValue * (cos(dw*t) * (x*x - y*y) + sin(dw*t) * (2*x*y));
    d2[0] = (term1+term2)*x/sr-w2Value*2*x+w4Value*12*x*(sr*sr-5*z*z)+cos(dw*t)*ww2Value*2*x+sin(dw*t)*ww2Value*2*y;
    d2[1] = (term1+term2)*y/sr-w2Value*2*y+w4Value*12*y*(sr*sr-5*z*z)-cos(dw*t)*ww2Value*2*y+sin(dw*t)*ww2Value*2*x;
    d2[2] = (term1+term2)*z/sr+w2Value*4*z+w4Value*16*z*(5*z*z-3*sr*sr);
    d2[0]*=-4*PI*mu0;
    d2[1]*=-4*PI*mu0;
    d2[2]*=-4*PI*mu0;
}

void solve(Particle *p, double t)
{
    //Initial Condition
    double x = p->x;
    double y = p->y;
    double z = p->z;
    double vx = p->vx;
    double vy = p->vy;
    double vz = p->vz;
    double newX, newY, newZ, newVx, newVy, newVz;

    double d[3];
    double d2[3];

    //Update on k1
    deri(d2, d, x, y, z, vx, vy, vz, t);
    double kx11 = d[0]; //k1 for first ODE
    newX = x + stepSize/2*kx11; //Update x(t)
    double kx12 = d2[0]; //k1 for second ODE
    newVx = vx + stepSize/2*kx12; //Update x'(t)
    double ky11 = d[1]; //k1 for first ODE
    newY = y + stepSize/2*ky11; //Update y(t)
    double ky12 = d2[1]; //k1 for second ODE
    newVy = vy + stepSize/2*ky12; //Update y'(t)
    double kz11 = d[2]; //k1 for first ODE
    newZ = z + stepSize/2*kz11; //Update z(t)
    double kz12 = d2[2]; //k1 for second ODE
    newVz = vz + stepSize/2*kz12; //Update z'(t)

    //Update on k2
    deri(d2, d, newX, newY, newZ, newVx, newVy, newVz, t+stepSize/2);
    double kx21 = d[0]; //k2 for first ODE
    newX = x + stepSize/2*kx21; //Update x(t)
    double kx22 = d2[0]; //k2 for second ODE
    newVx = vx + stepSize/2*kx22; //Update x'(t)
    double ky21 = d[1]; //k2 for first ODE
    newY = y + stepSize/2*ky21; //Update y(t)
    double ky22 = d2[1]; //k2 for second ODE
    newVy = vy + stepSize/2*ky22; //Update y'(t)
    double kz21 = d[2]; //k2 for first ODE
    newZ = z + stepSize/2*kz21; //Update z(t)
    double kz22 = d2[2]; //k2 for second ODE
    newVz = vz + stepSize/2*kz22; //Update z'(t)

    //Update on k3
    deri(d2, d, newX, newY, newZ, newVx, newVy, newVz, t+stepSize/2);
    double kx31 = d[0]; //k3 for first ODE
    newX = x + stepSize*kx31; //Update x(t)
    double kx32 = d2[0]; //k3 for second ODE
    newVx = vx + stepSize*kx32; //Update x'(t)
    double ky31 = d[1]; //k3 for first ODE
    newY = y + stepSize*ky31; //Update y(t)
    double ky32 = d2[1]; //k3 for second ODE
    newVy = vy + stepSize*ky32; //Update y'(t)
    double kz31 = d[2]; //k3 for first ODE
    newZ = z + stepSize*kz31; //Update z(t)
    double kz32 = d2[2]; //k3 for second ODE
    newVz = vz + stepSize*kz32; //Update z'(t)

    //Update on k4
    deri(d2, d, newX, newY, newZ, newVx, newVy, newVz, t+stepSize);
    double kx41 = d[0]; //k4 for first ODE
    double kx42 = d2[0]; //k4 for second ODE
    double ky41 = d[1]; //k4 for first ODE
    double ky42 = d2[1]; //k4 for second ODE
    double kz41 = d[2]; //k4 for first ODE
    double kz42 = d2[2]; //k4 for second ODE

    //Update values
    p->x = x + stepSize*(kx11/6 + kx21/3 + kx31/3 + kx41/6);
    p->vx = vx + stepSize*(kx12/6 + kx22/3 + kx32/3 + kx42/6);
    p->y = y + stepSize*(ky11/6 + ky21/3 + ky31/3 + ky41/6);
    p->vy = vy + stepSize*(ky12/6 + ky22/3 + ky32/3 + ky42/6);
    p->z = z + stepSize*(kz11/6 + kz21/3 + kz31/3 + kz41/6);
    p->vz = vz + stepSize*(kz12/6 + kz22/3 + kz32/3 + kz42/6);
    //cout<<kx11<<" "<<kx21<<" "<<kx31<<" "<<kx41<<endl;
}

#endif
