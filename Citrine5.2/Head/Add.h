#pragma once
#include <Eigen/Dense>
#include "ReadOut.h"
#include "ParaReader.h"

class Add
{
private:
    Eigen::VectorXd Yold, Ynew;

    double A, rho, d0, E, I, M, ma, Cdt, Cdn, Cdb, pi, g, Gx, Gy, Gz;
    double Vx, Vy, Vz, Vtx, Vty, Vtz, Vbx, Vby, Vbz;
    double deltaT, deltaS;
    double Gbx, Gby, Gbz, Ax, Ay, Az;
    Eigen::Vector3d calculatePoint00(const Eigen::VectorXd& Y);
    Eigen::Vector3d calculatePoint02(const Eigen::VectorXd& Y);

public:

    Add(Eigen::VectorXd &arr, Eigen::VectorXd &brr, int index);
    ~Add();

    Eigen::VectorXd Addyold();
    Eigen::VectorXd Addynew();
};

