#pragma once
#include <iostream>
#include <Eigen/Dense>
#include "ReadOut.h"

using namespace std;
using namespace Eigen;

class Add
{
    private:
        VectorXd Yold;
        VectorXd Ynew;

        double A;
        double rho;
        double d0;
        double E;
        double I;

        double Vz;//Velocity of water
        double Vx;
        double Vy;


        double M;//M = RhoCable*▽ = RhoCable*A*diameter(计算长度)
        double ma;//ma = [(Cm*RhoWater*▽) + M]  此时Cm = 1.00
        double w0;//
        double Cdt;
        double Cdn;
        double Cdb;

        
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

        double Gx;
        double Gy;
        double Gz;
        double Ax;
        double Ay;
        double Az;

        int k;


    public:
            Add(VectorXd &arr, VectorXd &brr, int index);
            ~Add();

            VectorXd Addyold();
            VectorXd Addynew();
};