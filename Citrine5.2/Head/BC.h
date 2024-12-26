#pragma once
#include <Eigen/Dense>

class BC
{
private:
    Eigen::VectorXd Yold;
    Eigen::VectorXd Ynew;
    Eigen::VectorXd process(Eigen::VectorXd &Y);

public:
    BC(Eigen::VectorXd &arr, Eigen::VectorXd &brr);
    ~BC();

    Eigen::VectorXd yold();
    Eigen::VectorXd ynew();
};
