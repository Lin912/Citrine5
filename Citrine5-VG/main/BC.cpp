#include <iostream>
#include <Eigen/Dense>
#include "BC.h"

BC::BC(VectorXd &arr, VectorXd &brr)
{
    Yold = arr;
    Ynew = brr;
}

BC::~BC()
{
    // cout << "Bundary condition is Load !" << endl;

}

VectorXd BC::yold()
{
    VectorXd temp(500);

    VectorXd point0(10);//零号边界点
    point0(0) = Yold(0);
    point0(1) = Yold(1);
    point0(2) = Yold(2);
    point0(3) = Yold(3);
    point0(4) = Yold(4);
    point0(5) = Yold(5);
    point0(6) = Yold(6);
    point0(7) = Yold(7);
    point0(8) = Yold(8);
    point0(9) = Yold(9);

    VectorXd point10(10);//尾号边界点
    point10(0) = Yold(490);
    point10(1) = Yold(491);
    point10(2) = Yold(492);
    point10(3) = Yold(493);
    point10(4) = Yold(494);
    point10(5) = Yold(495);
    point10(6) = Yold(496);
    point10(7) = Yold(497);
    point10(8) = Yold(498);
    point10(9) = Yold(499);

    temp.head(10) = point0;
    temp.segment(10, 480) = Yold.segment(10, 480);
    temp.tail(10) = point10;

    return temp;
}

VectorXd BC::ynew()
{
    VectorXd temp(500);

    VectorXd point0(10);//零号边界点
    point0(0) = Ynew(0);
    point0(1) = Ynew(1);
    point0(2) = Ynew(2);
    point0(3) = Ynew(3);
    point0(4) = Ynew(4);
    point0(5) = Ynew(5);
    point0(6) = Ynew(6);
    point0(7) = Ynew(7);
    point0(8) = Ynew(8);
    point0(9) = Ynew(9);

    VectorXd point10(10);//尾号边界点
    point10(0) = Ynew(490);
    point10(1) = Ynew(491);
    point10(2) = Ynew(492);
    point10(3) = Ynew(493);
    point10(4) = Ynew(494);
    point10(5) = Ynew(495);
    point10(6) = Ynew(496);
    point10(7) = Ynew(497);
    point10(8) = Ynew(498);
    point10(9) = Ynew(499);


    temp.head(10) = point0;
    temp.segment(10, 480) = Ynew.segment(10, 480);
    temp.tail(10) = point10;

    return temp;

}