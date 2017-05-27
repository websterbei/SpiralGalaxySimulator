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
    double mass = atof(argv[1]);  //Need to specify

    //init*[0] value, init*[1] first order derivative
    double initX[2];
    double initY[2];
    double initZ[2];
    double con[3];
    con[0] = 2.0/mass;
    double step = 0.01;

    //Initialize
    initX[0] = atof(argv[2]);
    initX[1] = atof(argv[3]);
    initY[0] = atof(argv[4]);
    initY[1] = atof(argv[5]);
    initZ[0] = atof(argv[6]);
    initZ[1] = atof(argv[7]);
    double X[110000];
    double Y[110000];
    double Z[110000];

    X[0] = initX[0];
    Y[0] = initY[0];
    Z[0] = initZ[0];

    double dydt[2];
    double newX[2];
    double newY[2];
    double newZ[2];
    int stepCounter = 0;

    X[0] = initX[0];
    Y[0] = initY[0];
    Z[0] = initZ[0];

    double t = 0;

    while(t<1000)
    {
	     //In X direction
        con[1] = initY[0];
        con[2] = initZ[0];
        deri(dydt, initX, con);
        newX[1] = initX[1] + dydt[1] * step;
        newX[0] = initX[0] + dydt[0] * step;
      //In Y direction
        con[1] = initX[0];
        con[2] = initZ[0];
        deri(dydt, initY, con);
        newY[1] = initY[1] + dydt[1] * step;
        newY[0] = initY[0] + dydt[0] * step;
      //In Z direction
        con[1] = initY[0];
        con[2] = initX[0];
        deri(dydt, initZ, con);
        newZ[1] = initZ[1] + dydt[1] * step;
        newZ[0] = initZ[0] + dydt[0] * step;
      //Update
        initX[0] = newX[0];
        initX[1] = newX[1];
        initY[0] = newY[0];
        initY[1] = newY[1];
        initZ[0] = newZ[0];
        initZ[1] = newZ[1];
        stepCounter++;
        X[stepCounter] = initX[0];
        Y[stepCounter] = initY[0];
        Z[stepCounter] = initZ[0];
        cout<<stepCounter<<endl;
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
