#include <iostream>
#include <cmath>
#include <fstream>
#include <cstdlib>
using namespace std;

//con[0] first constant 2/m, con[1] con[2] other coordinates
void deri(double dydt[], double y[], double con[]) //dydt[0] y', dydt[1] y'', y[0] old y, y[1] old y'
{
    dydt[0] = y[1];
    dydt[1] = -con[0]*y[0]*exp(-(y[0]*y[0]+con[1]*con[1]+con[2]*con[2]));
    return;
}

int main(int argc, char** argv)
{
    //init*[0] value, init*[1] first order derivative
    double initX[2];
    double initY[2];
    double initZ[2];
    double con[3];
    double step = 0.001;

    //Initialize

    //Initial Condition
    double mass = atof(argv[1]);  //Need to specify
    initX[0] = atof(argv[2]);
    initX[1] = atof(argv[3]);
    initY[0] = atof(argv[4]);
    initY[1] = atof(argv[5]);
    initZ[0] = atof(argv[6]);
    initZ[1] = atof(argv[7]);

    //Fixed coefficient
    con[0] = 2.0/mass;

    //Result array
    double X[110000];
    double Y[110000];
    double Z[110000];

    //Interation array
    X[0] = initX[0];
    Y[0] = initY[0];
    Z[0] = initZ[0];

    double dydt[2];
    double newX[2];
    double newY[2];
    double newZ[2];
    int stepCounter = 0;
    double t = 0;

    while(t<100)
    {
	     //Update on k1
      con[1] = initY[0];
      con[2] = initZ[0];
      deri(dydt, initX, con);
      double kx11 = dydt[0]; //k1 for first ODE
      newX[0] = initX[0] + step/2*kx11; //Update x(t)
      double kx12 = dydt[1]; //k1 for second ODE
      newX[1] = initX[1] + step/2*kx12; //Update x'(t)
      con[1] = initX[0];
      con[2] = initZ[0];
      deri(dydt, initY, con);
      double ky11 = dydt[0]; //k1 for first ODE
      newY[0] = initY[0] + step/2*ky11; //Update y(t)
      double ky12 = dydt[1]; //k1 for second ODE
      newY[1] = initY[1] + step/2*ky12; //Update y'(t)
      con[1] = initY[0];
      con[2] = initZ[0];
      deri(dydt, initZ, con);
      double kz11 = dydt[0]; //k1 for first ODE
      newZ[0] = initZ[0] + step/2*kz11; //Update z(t)
      double kz12 = dydt[1]; //k1 for second ODE
      newZ[1] = initZ[1] + step/2*kz12; //Update z'(t)
      //Update on k2
      con[1] = initY[0];
      con[2] = initZ[0];
      deri(dydt, newX, con);
      double kx21 = dydt[0]; //k2 for first ODE
      newX[0] = initX[0] + step/2*kx21; //Update x(t)
      double kx22 = dydt[1]; //k2 for second ODE
      newX[1] = initX[1] + step/2*kx22; //Update x'(t)
      con[1] = initX[0];
      con[2] = initZ[0];
      deri(dydt, newY, con);
      double ky21 = dydt[0]; //k2 for first ODE
      newY[0] = initY[0] + step/2*ky21; //Update y(t)
      double ky22 = dydt[1]; //k2 for second ODE
      newY[1] = initY[1] + step/2*ky22; //Update y'(t)
      con[1] = initY[0];
      con[2] = initZ[0];
      deri(dydt, newZ, con);
      double kz21 = dydt[0]; //k2 for first ODE
      newZ[0] = initZ[0] + step/2*kz21; //Update z(t)
      double kz22 = dydt[1]; //k2 for second ODE
      newZ[1] = initZ[1] + step/2*kz22; //Update z'(t)
      //Update on k3
      con[1] = initY[0];
      con[2] = initZ[0];
      deri(dydt, newX, con);
      double kx31 = dydt[0]; //k3 for first ODE
      newX[0] = initX[0] + step*kx31; //Update x(t)
      double kx32 = dydt[1]; //k3 for second ODE
      newX[1] = initX[1] + step*kx32; //Update x'(t)
      con[1] = initX[0];
      con[2] = initZ[0];
      deri(dydt, newY, con);
      double ky31 = dydt[0]; //k3 for first ODE
      newY[0] = initY[0] + step*ky31; //Update y(t)
      double ky32 = dydt[1]; //k3 for second ODE
      newY[1] = initY[1] + step*ky32; //Update y'(t)
      con[1] = initY[0];
      con[2] = initZ[0];
      deri(dydt, newZ, con);
      double kz31 = dydt[0]; //k3 for first ODE
      newZ[0] = initZ[0] + step*kz31; //Update z(t)
      double kz32 = dydt[1]; //k3 for second ODE
      newZ[1] = initZ[1] + step*kz32; //Update z'(t)
      //Update on k4
      con[1] = initY[0];
      con[2] = initZ[0];
      deri(dydt, newX, con);
      double kx41 = dydt[0]; //k4 for first ODE
      double kx42 = dydt[1]; //k4 for second ODE
      con[1] = initX[0];
      con[2] = initZ[0];
      deri(dydt, newY, con);
      double ky41 = dydt[0]; //k4 for first ODE
      double ky42 = dydt[1]; //k4 for second ODE
      con[1] = initY[0];
      con[2] = initZ[0];
      deri(dydt, newZ, con);
      double kz41 = dydt[0]; //k4 for first ODE
      double kz42 = dydt[1]; //k4 for second ODE
      initX[0] = initX[0] + step*(kx11/6 + kx21/3 + kx31/3 + kx41/6);
      initX[1] = initX[1] + step*(kx12/6 + kx22/3 + kx32/3 + kx42/6);
      initY[0] = initY[0] + step*(ky11/6 + ky21/3 + ky31/3 + ky41/6);
      initY[1] = initY[1] + step*(ky12/6 + ky22/3 + ky32/3 + ky42/6);
      initZ[0] = initZ[0] + step*(kz11/6 + kz21/3 + kz31/3 + kz41/6);
      initZ[1] = initZ[1] + step*(kz12/6 + kz22/3 + kz32/3 + kz42/6);
      //Update Values
      stepCounter++;
      X[stepCounter] = initX[0];
      Y[stepCounter] = initY[0];
      Z[stepCounter] = initZ[0];
      //cout<<stepCounter<<endl;
      t+= step;
    }

    ofstream outputFile;
    outputFile.open("output.txt");
    for(int i=0; i<stepCounter; i++)
    {
      outputFile<<X[i]<<" "<<Y[i]<<" "<<Z[i]<<endl;
    }
    outputFile.close();

    return 0;
}
