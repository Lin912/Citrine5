#include "../Head/Add.h"

Add::Add(const VectorXd &arr, const VectorXd &brr, int index):Yold(arr), Ynew(brr)
{

    FiberRO a;
    PhysicalData physicalData = ParaReader::ReadAllPhysicalData(a, index);
    Vtx = physicalData.Vtx;
    Vty = physicalData.Vty;
    Vtz = physicalData.Vtz;

    Vbx = physicalData.Vbx;
    Vby = physicalData.Vby;
    Vbz = physicalData.Vbz;
}

Add::~Add() {}

Vector3d Add::calculatePoint00(const VectorXd& Y)
{
    Vector3d point;
    point(0) = Vtx * cos(Y(7)) * cos(Y(6)) + Vty * cos(Y(6)) * sin(Y(7)) - Vtz * sin(Y(6));
    point(1) = Vty * cos(Y(7)) - Vtx * sin(Y(7));
    point(2) = Vtx * cos(Y(7)) * sin(Y(6)) + Vty * sin(Y(6)) * sin(Y(7)) + Vtz * cos(Y(6));
    return point;
}

Vector3d Add::calculatePoint02(const VectorXd& Y)
{
    Vector3d point;
    point(0) = Vbx*cos(Y(497))*cos(Y(496)) + Vby*cos(Y(496))*sin(Y(497)) - Vbz*sin(Y(496));
    point(1) = Vby*cos(Y(497)) - Vbx*sin(Y(497));
    point(2) = Vbx*cos(Y(497))*sin(Y(496)) + Vby*sin(Y(496))*sin(Y(497)) + Vbz*cos(Y(496));
    return point;
}

VectorXd Add::Addyold()
{
    VectorXd temp(500);

    Vector3d point00 = calculatePoint00(Yold);
    VectorXd point01(2);
    point01.setZero();
    Vector3d point02 = calculatePoint02(Yold);
    VectorXd point03(2);
    point03.setZero();

    temp.segment(0, 3) = point00;
    temp.segment(3, 5) = Yold.segment(3, 5);
    temp.segment(8, 2) = point01;
    temp.segment(10, 480) = Yold.segment(10, 480);
    temp.segment(490, 3) = point02;
    temp.segment(493, 5) = Yold.segment(493, 5);
    temp.tail(2) = point03;

    return temp;
}

VectorXd Add::Addynew()
{
    VectorXd temp(500);

    Vector3d point00 = calculatePoint00(Ynew);
    VectorXd point01(2);
    point01.setZero();
    Vector3d point02 = calculatePoint02(Ynew);
    VectorXd point03(2);
    point03.setZero();

    temp.segment(0, 3) = point00;
    temp.segment(3, 5) = Ynew.segment(3, 5);
    temp.segment(8, 2) = point01;
    temp.segment(10, 480) = Ynew.segment(10, 480);
    temp.segment(490, 3) = point02;
    temp.segment(493, 5) = Ynew.segment(493, 5);
    temp.tail(2) = point03;

    return temp;
}
