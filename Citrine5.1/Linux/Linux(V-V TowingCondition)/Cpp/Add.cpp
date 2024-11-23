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
    
    Vtx = a.ReadTopVel(index)[0];
    Vty = a.ReadTopVel(index)[1];
    Vtz = a.ReadTopVel(index)[2];

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
    //Yold(7) == phi        Yold(6) == theta
    point00(0) = Vtx*cos(Yold(7))*cos(Yold(6)) + Vty*cos(Yold(6))*sin(Yold(7))- Vtz*sin(Yold(6));
    point00(1) = Vty*cos(Yold(7)) - Vtx*sin(Yold(7));
    point00(2) = Vtx*cos(Yold(7))*sin(Yold(6)) + Vty*sin(Yold(6))*sin(Yold(7)) + Vtz*cos(Yold(6));

    VectorXd point01(2);
    point01(0) = 0; //Omega_2
    point01(1) = 0; //Omega_3

    VectorXd point02(3);
    point02(0) =  Vbx*cos(Yold(497))*cos(Yold(496)) + Vby*cos(Yold(496))*sin(Yold(497))- Vbz*sin(Yold(496));
    point02(1) =  Vby*cos(Yold(497)) - Vbx*sin(Yold(497));
    point02(2) =  Vbx*cos(Yold(497))*sin(Yold(496)) + Vby*sin(Yold(496))*sin(Yold(497)) + Vbz*cos(Yold(496));

    VectorXd point03(2);
    point03(0) = 0;          
    point03(1) = 0;           

    temp.segment(0,3) = point00;
    temp.segment(3,5) = Yold.segment(3, 5);
    temp.segment(8,2) = point01;
    temp.segment(10, 480) = Yold.segment(10, 480);
    temp.segment(490, 3) = point02;
    temp.segment(493, 5) = Yold.segment(493, 5);
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
    point01(0) = 0;             
    point01(1) = 0;             

    VectorXd point02(3);
    point02(0) =  Vbx*cos(Ynew(497))*cos(Ynew(496)) + Vby*cos(Ynew(496))*sin(Ynew(497))- Vbz*sin(Ynew(496));
    point02(1) =  Vby*cos(Ynew(497)) - Vbx*sin(Ynew(497));
    point02(2) =  Vbx*cos(Ynew(497))*sin(Ynew(496)) + Vby*sin(Ynew(496))*sin(Ynew(497)) + Vbz*cos(Ynew(496));

    VectorXd point03(2);
    point03(0) = 0;             
    point03(1) = 0;           

    temp.segment(0,3) = point00;
    temp.segment(3,5) = Ynew.segment(3, 5);
    temp.segment(8,2) = point01;
    temp.segment(10, 480) = Ynew.segment(10, 480);
    temp.segment(490, 3) = point02;
    temp.segment(493, 5) = Ynew.segment(493, 5);
    temp.tail(2) = point03;

    return temp;
}