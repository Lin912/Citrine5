#pragma once
#include <iostream>
#include <Eigen/Dense>
#include "Read&Out.h"

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

        double V1;//Velocity of water
        double V2;
        double V3;


        double M;//M = RhoCable*▽ = RhoCable*A*diameter(计算长度)
        double ma;//ma = [(Cm*RhoWater*▽) + M]  此时Cm = 1.00
        double w0;//
        double Cdt;
        double Cdn;
        double Cdb;

        
        double pi;
        double g;

        double Vt1;
        double Vt2;
        double Vt3;

        double Vb1;
        double Vb2;
        double Vb3;

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