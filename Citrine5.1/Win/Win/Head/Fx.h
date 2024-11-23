#pragma once

#include <iostream>
#include <cmath>
#include <Eigen/Dense>
#include <vector>
#include "ReadOut.h"
#include "MNQ.h"


using namespace std;
using namespace Eigen;

class Fx
{
    private:
        VectorXd Yold;
        VectorXd Ynew;
        VectorXd temp;

        double A;
        double rho;
        double d0;
        double E;
        double I;
        // double G;
        // double Ip;

        double Vz;//Velocity of water
        double Vx;
        double Vy;

        double M;//M = RhoCable*▽ = RhoCable*A*diameter
        double ma;//ma = [(Cm*RhoWater*▽) + M]
        double w0;
        double Cdt;
        double Cdn;
        double Cdb;
        double Gx;
        double Gy;
        double Gz;

        double pi;
        double g;

        double Vtx;
        double Vty;
        double Vtz;

        double Vbx;
        double Vby;
        double Vbz;


        double deltaT;
        double deltaS;

        double Gbx;
        double Gby;
        double Gbz;
        double Ax;
        double Ay;
        double Az;

        int k;


    public:
        Fx(VectorXd& arr, VectorXd& brr, int index);
        ~Fx();

        VectorXd fx();
};