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
    Cdt = a.ReadPhysical()[7];
    Cdn = a.ReadPhysical()[8];
    Cdb = a.ReadPhysical()[9];
    pi = a.ReadPhysical()[10];
    g = a.ReadPhysical()[11];
    Gx = a.ReadPhysical()[12];
    Gy = a.ReadPhysical()[13];
    Gz = a.ReadPhysical().back();

    Vx = a.ReadWater(index)[0];
    Vy = a.ReadWater(index)[1];
    Vz = a.ReadWater(index)[2];
    
    //Model(V-G)
    Vtx = a.ReadTopVel(index)[0];
    Vty = a.ReadTopVel(index)[1];
    Vtz = a.ReadTopVel(index)[2];

    //Model(V-V)
    // Vtz = a.ReadTopVel()[0];
    // Vtx = a.ReadTopVel()[1];
    // Vty = a.ReadTopVel()[2];

    Vbx = a.ReadBottomVel()[0];
    Vby = a.ReadBottomVel()[1];
    Vbz = a.ReadBottomVel()[2];

    deltaT = a.ReadDelta()[0];
    deltaS = a.ReadDelta()[1];

    Gbx = a.ReadBottomG()[0];
    Gby = a.ReadBottomG()[1];
    Gbz = a.ReadBottomG()[2];
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
    point00(0) = Vtx*cos(Yold(7))*cos(Yold(6)) + Vty*cos(Yold(6))*sin(Yold(7))- Vtz*sin(Yold(6));
    point00(1) = Vty*cos(Yold(7)) - Vtx*sin(Yold(7));
    point00(2) = Vtx*cos(Yold(7))*sin(Yold(6)) + Vty*sin(Yold(6))*sin(Yold(7)) + Vtz*cos(Yold(6));
    VectorXd point01(2);
    point01(0) = 0;              //O2mega边界条件（上）
    point01(1) = 0;              //O3mega边界条件


    VectorXd point02(3);
    point02(0) =  Gbx*cos(Yold(497))*cos(Yold(496)) + Gby*cos(Yold(496))*sin(Yold(497))- Gbz*sin(Yold(496));
    point02(1) =  Gby*cos(Yold(497)) - Gbx*sin(Yold(497));
    point02(2) =  Gbx*cos(Yold(497))*sin(Yold(496)) + Gby*sin(Yold(496))*sin(Yold(497)) + Gbz*cos(Yold(496));

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
    point00(0) = Vtx*cos(Ynew(7))*cos(Ynew(6)) + Vty*cos(Ynew(6))*sin(Ynew(7))- Vtz*sin(Ynew(6));
    point00(1) = Vty*cos(Ynew(7)) - Vtx*sin(Ynew(7));
    point00(2) = Vtx*cos(Ynew(7))*sin(Ynew(6)) + Vty*sin(Ynew(6))*sin(Ynew(7)) + Vtz*cos(Ynew(6));

    VectorXd point01(2);
    point01(0) = 0;              //O2mega边界条件（上）
    point01(1) = 0;              //O3mega边界条件

    VectorXd point02(3);
    point02(0) =  Gbx*cos(Ynew(497))*cos(Ynew(496)) + Gby*cos(Ynew(496))*sin(Ynew(497))- Gbz*sin(Ynew(496));
    point02(1) =  Gby*cos(Ynew(497)) - Gbx*sin(Ynew(497));
    point02(2) =  Gbx*cos(Ynew(497))*sin(Ynew(496)) + Gby*sin(Ynew(496))*sin(Ynew(497)) + Gbz*cos(Ynew(496));


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