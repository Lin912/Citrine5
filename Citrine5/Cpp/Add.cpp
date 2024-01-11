#include "../Head/Add.h"

Add::Add(VectorXd &arr, VectorXd &brr, int index)
{
    Yold = arr;
    Ynew = brr;

    k = index;

    FiberRO a;
    A = a.ReadPhysical()[0];
    rho = a.ReadPhysical()[1];
    d0 = a.ReadPhysical()[2];
    E = a.ReadPhysical()[3];
    I = a.ReadPhysical()[4];
    M = a.ReadPhysical()[5];
    ma = a.ReadPhysical()[6];
    w0 = a.ReadPhysical()[7];
    Cdt = a.ReadPhysical()[8];
    Cdn = a.ReadPhysical()[9];
    Cdb = a.ReadPhysical()[10];
    
    pi = a.ReadPhysical()[11];
    g = a.ReadPhysical().back();

    V1 = a.ReadWater(index)[0];
    V2 = a.ReadWater(index)[1];
    V3 = a.ReadWater(index)[2];

    Vt1 = a.ReadTopVel(index)[0];
    Vt2 = a.ReadTopVel(index)[1];
    Vt3 = a.ReadTopVel(index)[2];

    Vb1 = a.ReadBottomVel()[0];
    Vb2 = a.ReadBottomVel()[1];
    Vb3 = a.ReadBottomVel()[2];

    deltaT = a.ReadDelta()[0];
    deltaS = a.ReadDelta()[1];

    Gx = a.ReadBottomG()[0];
    Gy = a.ReadBottomG()[1];
    Gz = a.ReadBottomG()[2];
    Ax = a.ReadBottomG()[3];
    Ay = a.ReadBottomG()[4];
    Az = a.ReadBottomG()[5];
}

Add::~Add()
{

}


VectorXd Add::Addyold()
{
    VectorXd temp(500);

    VectorXd point00(3);
    point00(0) = Vt1*cos(Yold(6))*cos(Yold(7)) - Vt2*sin(Yold(7))*cos(Yold(6)) - Vt3*sin(Yold(6));          //速度u边界条件（上）
    point00(1) = Vt2*cos(Yold(7)) + Vt1*sin(Yold(7));          //速度V边界条件
    point00(2) = Vt1*sin(Yold(6))*cos(Yold(7)) - Vt2*sin(Yold(7))*sin(Yold(6)) + Vt3*cos(Yold(6));            //速度w边界条件
    VectorXd point01(2);

    point01(0) = 0;              //O2mega边界条件（上）
    point01(1) = 0;              //O3mega边界条件

    // VectorXd point02(3);
    // point02(0) = Vb1*cos(Yold(496))*cos(Yold(497)) - Vb2*sin(Yold(497))*cos(Yold(496)) - Vb3*sin(Yold(496));          //速度T边界条件（下）
    // point02(1) = Vb2*cos(Yold(497)) + Vb1*sin(Yold(497));          //速度Sn边界条件
    // point02(2) = Vb1*sin(Yold(496))*cos(Yold(497)) - Vb2*sin(Yold(497))*sin(Yold(496)) + Vb3*cos(Yold(496));             //速度Sb边界条件


    VectorXd point02(3);
    point02(0) =  Gx*cos(Yold(496))*cos(Yold(497)) - Gy*sin(Yold(497))*cos(Yold(496)) - Gz*sin(Yold(496));
    point02(1) =  Gy*cos(Yold(497)) + Gx*sin(Yold(497));
    point02(2) =  Gx*sin(Yold(496))*cos(Yold(497)) - Gy*sin(Yold(497))*sin(Yold(496)) + Gz*cos(Yold(496));

    VectorXd point03(2);
    point03(0) = 0;              //O2mega边界条件（下）
    point03(1) = 0;              //O3mega边界条件

    temp.segment(0,3) = point00;
    temp.segment(3,5) = Yold.segment(3, 5);
    temp.segment(8,2) = point01;


    temp.segment(10, 480) = Yold.segment(10, 480);

    temp.segment(490, 3) = Yold.segment(490, 3);
    temp.segment(493, 3) = point02;
    temp.segment(496, 2) = Yold.segment(496, 2);
    temp.tail(2) = point03;

    return temp;
}

VectorXd Add::Addynew()
{
    VectorXd temp(500);

    VectorXd point00(3);
    point00(0) = Vt1*cos(Ynew(6))*cos(Ynew(7)) - Vt2*sin(Ynew(7))*cos(Ynew(6)) - Vt3*sin(Ynew(6));            //速度u边界条件（上）
    point00(1) = Vt2*cos(Ynew(7)) + Vt1*sin(Ynew(7));            //速度V边界条件
    point00(2) = Vt1*sin(Ynew(6))*cos(Ynew(7)) - Vt2*sin(Ynew(7))*sin(Ynew(6)) + Vt3*cos(Ynew(6));              //速度w边界条件

    VectorXd point01(2);
    point01(0) = 0;              //O2mega边界条件（上）
    point01(1) = 0;              //O3mega边界条件

    // VectorXd point02(3);
    // point02(0) = Vb1*cos(Ynew(496))*cos(Ynew(497)) - Vb2*sin(Ynew(497))*cos(Ynew(496)) - Vb3*sin(Ynew(496));           
    // point02(1) = Vb2*cos(Ynew(497)) + Vb1*sin(Ynew(497));
    // point02(2) = Vb1*sin(Ynew(496))*cos(Ynew(497)) - Vb2*sin(Ynew(497))*sin(Ynew(496)) + Vb3*cos(Ynew(496));

    VectorXd point02(3);
    point02(0) =  Gx*cos(Ynew(496))*cos(Ynew(497)) - Gy*sin(Ynew(497))*cos(Ynew(496)) - Gz*sin(Ynew(496));
    point02(1) =  Gy*cos(Ynew(497)) + Gx*sin(Ynew(497));
    point02(2) =  Gx*sin(Ynew(496))*cos(Ynew(497)) - Gy*sin(Ynew(497))*sin(Ynew(496)) + Gz*cos(Ynew(496));


    VectorXd point03(2);
    point03(0) = 0;              //O2mega边界条件（下）
    point03(1) = 0;              //O3mega边界条件

    temp.segment(0,3) = point00;
    temp.segment(3,5) = Ynew.segment(3, 5);
    temp.segment(8,2) = point01;


    temp.segment(10, 480) = Ynew.segment(10, 480);

    temp.segment(490, 3) = Ynew.segment(490, 3);
    temp.segment(493, 3) = point02;
    temp.segment(496, 2) = Ynew.segment(496, 2);
    temp.tail(2) = point03;

    return temp;
}