#include <iostream>
#include <cmath>
using namespace std;

//con[0] first constant 2/m, con[1] con[2] other coordinates
void deri(double dydt[], double y[], double con[]) //dydt[0] y', dydt[1] y'', y[0] old y, y[1] old y'
{
    dydt[0] = y[1];
    dydt[1] = con[0]*y[0]*exp(-(y[0]*y[0]+con[1]*con[1]+con[2]*con[2]));
    return;
}

int main()
{
    double mass = //Need to specify

    //init*[0] value, init*[1] first order derivative
    double initX[2];
    double initY[2];
    double initZ[2];
    double con[3];
    con[0] = 2.0/mass;
    double step = 0.0001;

    //Initialize
    double X[11000];
    double Y[11000];
    double Z[11000];

    X[0] = initX[0];
    Y[0] = initY[0];
    Z[0] = initZ[0];

    double dydt[2];
    double newX[2];
    double newY[2];
    double newZ[2];

    while(t<20)
    {
	     //In X direction
        con[1] = initY[0];
        con[2] = initZ[0];
        deri(dydt, initX, con);
        double k11 = dydt[0]; //k1 for first ODE
        newX[0] = initX[0] + step/2*k11; //Update x(t)
        double k12 = dydt[1]; //k1 for second ODE
        newX[1] = initX[1] + step/2*k12; //Update x'(t)
        initX[0] = newX[0];
        initX[1] = newX[1];
        deri(dydt, initX, con);
        double k21 = dydt[0]
    }
}
