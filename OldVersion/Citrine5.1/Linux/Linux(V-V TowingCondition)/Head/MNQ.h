#pragma once
#include <iostream>
#include <fstream>
#include <cmath>
#include <Eigen/Dense>
#include <vector>
#include "ReadOut.h"

using namespace std;
using namespace Eigen;

class MNQ
{
    private:
        VectorXd Yold;
        VectorXd Ynew;

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

        double Gx;
        double Gy;
        double Gz;

        double M;
        double ma;
        double w0;
        double Cdt;
        double Cdn;
        double Cdb;

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