#pragma once
#include <iostream>
#include <fstream>
#include <cmath>
#include <Eigen/Dense>
#include <vector>
#include "Read&Out.h"

using namespace std;
using namespace Eigen;

class MNQ
{
    private:
        VectorXd Yold;//用于存储临时Yold
        VectorXd Ynew;

        double A;//截面面积
        double rho;//介质密度
        double d0;//带缆直径 d0
        double E;//弹性模量 E
        double I;
        // double G;
        // double Ip;

        double V1;//Velocity of water
        double V2;
        double V3;

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



    public:
        MNQ(VectorXd& arr, VectorXd& brr, int index);
        ~MNQ();

        MatrixXd Mnew();
        MatrixXd Mold();
        MatrixXd Nold();
        MatrixXd Nnew();
        VectorXd Qold();
        VectorXd Qnew();

        void savetxt(Eigen::MatrixXd mat, string filename);

};