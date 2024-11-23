#pragma once
#include <iostream>
#include <cmath>
#include <Eigen/Dense>
#include <vector>
#include "MNQ.h"
#include "ReadOut.h"

class Jacobian
{
    private:
        VectorXd Yold;
        VectorXd Ynew;
        MatrixXd temp;

        double A;
        double rho;
        double d0;
        double E;
        double I;
        // double G;
        // double Ip;

        double Vx;//Velocity of water
        double Vy;
        double Vz;

        double M;
        double ma;
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
        Jacobian(VectorXd& arr, VectorXd& brr, int index);
        ~Jacobian();

        int sign(double a);
        //double D(double a);
        
        MatrixXd jacobian();

};