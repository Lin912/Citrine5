#pragma once
#include <iostream>
#include <Eigen/Dense>

using namespace std;
using namespace Eigen;

class BC
{
    private:
            VectorXd Yold;
            VectorXd Ynew;

    public:
            BC(VectorXd &arr, VectorXd &brr);
            ~BC();

            VectorXd yold();
            VectorXd ynew();
};