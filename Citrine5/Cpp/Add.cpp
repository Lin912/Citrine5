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

    Vz = a.ReadWater(index)[0];
    Vx = a.ReadWater(index)[1];
    Vy = a.ReadWater(index)[2];

    Vtz = a.ReadTopVel(index)[0];
    Vtx = a.ReadTopVel(index)[1];
    Vty = a.ReadTopVel(index)[2];

    Vbz = a.ReadBottomVel()[0];
    Vbx = a.ReadBottomVel()[1];
    Vby = a.ReadBottomVel()[2];

    deltaT = a.ReadDelta()[0];
    deltaS = a.ReadDelta()[1];

    Gz = a.ReadBottomG()[0];
    Gx = a.ReadBottomG()[1];
    Gy = a.ReadBottomG()[2];
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
    point00(0) = Vtz*cos(Yold(6))*cos(Yold(7)) + Vtx*sin(Yold(7))*cos(Yold(6)) - Vty*sin(Yold(6));         //速度u边界条件（上）
    point00(1) = Vtx*cos(Yold(7)) - Vtz*sin(Yold(7));         //速度V边界条件
    point00(2) = Vtz*cos(Yold(7))*sin(Yold(6)) + Vtx*sin(Yold(6))*sin(Yold(7)) + Vty*cos(Yold(6));          //速度w边界条件
    VectorXd point01(2);
    point01(0) = 0;              //O2mega边界条件（上）
    point01(1) = 0;              //O3mega边界条件


    // VectorXd point02(3);
    // point02(0) = Vb1*cos(Yold(496))*cos(Yold(497)) - Vb2*sin(Yold(497))*cos(Yold(496)) - Vb3*sin(Yold(496));          //速度T边界条件（下）
    // point02(1) = Vb2*cos(Yold(497)) + Vb1*sin(Yold(497));          //速度Sn边界条件
    // point02(2) = Vb1*sin(Yold(496))*cos(Yold(497)) - Vb2*sin(Yold(497))*sin(Yold(496)) + Vb3*cos(Yold(496));             //速度Sb边界条件


    VectorXd point02(3);
    double Sdn = 0.003162;
    double Sdb = 0.003162;
    point02(0) =  Gz*cos(Yold(496))*cos(Yold(497)) + Gx*sin(Yold(497))*cos(Yold(496)) - Gy*sin(Yold(496));
    point02(1) =  Gx*cos(Yold(497)) - Gz*sin(Yold(497));
    point02(2) =  Gz*cos(Yold(497))*sin(Yold(496)) + Gx*sin(Yold(496))*sin(Yold(497)) + Gy*cos(Yold(496));

    VectorXd point03(2);
    point03(0) = 0;              //O2mega边界条件（下）
    point03(1) = 0;              //O3mega边界条件



    temp.segment(0,3) = point00;
    temp.segment(3,5) = Yold.segment(3, 5);
    temp.segment(8,2) = point01;

    temp.segment(10, 480) = Yold.segment(10, 480);

    // temp.segment(490, 3) = point02;
    // temp.segment(493, 5) = Yold.segment(493, 5);
    // temp.tail(2) = point03;

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
    point00(0) = Vtz*cos(Ynew(6))*cos(Ynew(7)) + Vtx*sin(Ynew(7))*cos(Ynew(6)) - Vty*sin(Ynew(6));             //速度u边界条件（上）
    point00(1) = Vtx*cos(Ynew(7)) - Vtz*sin(Ynew(7));           //速度V边界条件
    point00(2) = Vtz*cos(Ynew(7))*sin(Ynew(6)) + Vtx*sin(Ynew(6))*sin(Ynew(7)) + Vty*cos(Ynew(6));           //速度w边界条件

    VectorXd point01(2);
    point01(0) = 0;              //O2mega边界条件（上）
    point01(1) = 0;              //O3mega边界条件

    //Switch (VV) or (VG) 

    //(VV)
    // VectorXd point02(3);
    // point02(0) = Vb1*cos(Ynew(496))*cos(Ynew(497)) - Vb2*sin(Ynew(497))*cos(Ynew(496)) - Vb3*sin(Ynew(496));           
    // point02(1) = Vb2*cos(Ynew(497)) + Vb1*sin(Ynew(497));
    // point02(2) = Vb1*sin(Ynew(496))*cos(Ynew(497)) - Vb2*sin(Ynew(497))*sin(Ynew(496)) + Vb3*cos(Ynew(496));

    //(VG)
    VectorXd point02(3);
    double Sdn = 0.003162;
    double Sdb = 0.003162;
    point02(0) =  Gz*cos(Ynew(496))*cos(Ynew(497)) + Gx*sin(Ynew(497))*cos(Ynew(496)) - Gy*sin(Ynew(496));
    point02(1) =  Gx*cos(Ynew(497)) - Gz*sin(Ynew(497));
    point02(2) =  Gz*cos(Ynew(497))*sin(Ynew(496)) + Gx*sin(Ynew(496))*sin(Ynew(497)) + Gy*cos(Ynew(496));


    VectorXd point03(2);
    point03(0) = 0;              //O2mega边界条件（下）
    point03(1) = 0;              //O3mega边界条件

    temp.segment(0,3) = point00;
    temp.segment(3,5) = Ynew.segment(3, 5);
    temp.segment(8,2) = point01;


    temp.segment(10, 480) = Ynew.segment(10, 480);

    //Switch (VV) or (VG) 

    //(VV)
    // temp.segment(490, 3) = point02;
    // temp.segment(493, 5) = Ynew.segment(493, 5);
    // temp.tail(2) = point03;

    //(VG)
    temp.segment(490, 3) = Ynew.segment(490, 3);
    temp.segment(493, 3) = point02;
    temp.segment(496, 2) = Ynew.segment(496, 2);
    temp.tail(2) = point03;

    return temp;
}