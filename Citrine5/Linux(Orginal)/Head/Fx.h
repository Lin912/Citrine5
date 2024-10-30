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

        double A;//截面面积
        double rho;//介质密度
        double d0;//带缆直径 d0
        double E;//弹性模量 E
        double I;
        // double G;
        // double Ip;

        double Vz;//Velocity of water
        double Vx;
        double Vy;

        //质量特性
        double M;//质量 M = RhoCable*▽ = RhoCable*A*diameter(计算长度)
        double ma;//附加质量 ma = [(Cm*RhoWater*▽) + M]  此时Cm = 1.00
        double w0;//单元重力
        double Cdt;
        double Cdn;
        double Cdb;

        //数学特性
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
        Fx(VectorXd& arr, VectorXd& brr, int index);
        ~Fx();

        VectorXd fx();
};