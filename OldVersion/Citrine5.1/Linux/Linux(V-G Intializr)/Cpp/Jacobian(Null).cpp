#include "../Head/Jacobian.h"

Jacobian::Jacobian(VectorXd& arr, VectorXd& brr, int index)
{
    Ynew = arr;
    Yold = brr;

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

Jacobian::~Jacobian()
{
    // cout << "Jacobian is read in !!" << endl;

}

MatrixXd Jacobian::jacobian()
{

    MatrixXd AA(500, 500);
    MatrixXd temp(500, 500);
    temp = AA.Zero(500, 500);
    
    //BC01
    temp(0, 0) = 1;
    temp(0, 1) = 0;
    temp(0, 2) = 0;
    temp(0, 3) = 0;
    temp(0, 4) = 0;
    temp(0, 5) = 0;
    temp(0, 6) = 0;   
    temp(0, 7) = Vty*cos(Ynew(7)) + Vtx*sin(Ynew(7));
    temp(0, 8) = 0;
    temp(0, 9) = 0;

    temp(1, 0) = 0;
    temp(1, 1) = 1;
    temp(1, 2) = 0;
    temp(1, 3) = 0;
    temp(1, 4) = 0;
    temp(1, 5) = 0;
    temp(1, 6) = Vtz*cos(Ynew(6)) + Vty*cos(Ynew(7))*sin(Ynew(6)) + Vtx*sin(Ynew(7))*sin(Ynew(6));
    temp(1, 7) = Vty*cos(Ynew(6))*sin(Ynew(7)) - Vtx*cos(Ynew(7))*cos(Ynew(6));
    temp(1, 8) = 0;
    temp(1, 9) = 0;

    temp(2, 0) = 0;
    temp(2, 1) = 0;
    temp(2, 2) = 1;
    temp(2, 3) = 0;
    temp(2, 4) = 0;
    temp(2, 5) = 0;
    temp(2, 6) = Vtz*sin(Ynew(6)) - Vty*cos(Ynew(7))*cos(Ynew(6)) - Vtx*cos(Ynew(6))*sin(Ynew(7));
    temp(2, 7) = Vty*sin(Ynew(7))*sin(Ynew(6)) - Vtx*cos(Ynew(7))*sin(Ynew(6));
    temp(2, 8) = 0;
    temp(2, 9) = 0;

    temp(3, 0) = 0;
    temp(3, 1) = 0;
    temp(3, 2) = 0;
    temp(3, 3) = 0;
    temp(3, 4) = 0;
    temp(3, 5) = 0;
    temp(3, 6) = 0;
    temp(3, 7) = 0;
    temp(3, 8) = 1;
    temp(3, 9) = 0;

    temp(4, 0) = 0;
    temp(4, 1) = 0;
    temp(4, 2) = 0;
    temp(4, 3) = 0;
    temp(4, 4) = 0;
    temp(4, 5) = 0;
    temp(4, 6) = 0;
    temp(4, 7) = 0;
    temp(4, 8) = 0;
    temp(4, 9) = 1;


    temp(495, 490) = 0;
    temp(495, 491) = 0;
    temp(495, 492) = 0;
    temp(495, 493) = 1;
    temp(495, 494) = 0;
    temp(495, 495) = 0;
    temp(495, 496) = 0; 
    temp(495, 497) = Gy*cos(Ynew(497)) + Gx*sin(Ynew(497));
    temp(495, 498) = 0;
    temp(495, 499) = 0;

    temp(496, 490) = 0;
    temp(496, 491) = 0;
    temp(496, 492) = 0;
    temp(496, 493) = 0;
    temp(496, 494) = 1;
    temp(496, 495) = 0;
    temp(496, 496) = Gz*cos(Ynew(496)) + Gy*cos(Ynew(497))*sin(Ynew(496)) + Gx*sin(Ynew(497))*sin(Ynew(496));
    temp(496, 497) = Gy*cos(Ynew(496))*sin(Ynew(497)) - Gx*cos(Ynew(497))*cos(Ynew(496));
    temp(496, 498) = 0;
    temp(496, 499) = 0;

    temp(497, 490) = 0;
    temp(497, 491) = 0;
    temp(497, 492) = 0;
    temp(497, 493) = 0;
    temp(497, 494) = 0;
    temp(497, 495) = 1;
    temp(497, 496) = Gz*sin(Ynew(496)) - Gy*cos(Ynew(497))*cos(Ynew(496)) - Gx*cos(Ynew(496))*sin(Ynew(497));
    temp(497, 497) = Gy*sin(Ynew(497))*sin(Ynew(496)) - Gx*cos(Ynew(497))*sin(Ynew(496));
    temp(497, 498) = 0;
    temp(497, 499) = 0;

    temp(498, 490) = 0;
    temp(498, 491) = 0;
    temp(498, 492) = 0;
    temp(498, 493) = 0;
    temp(498, 494) = 0;
    temp(498, 495) = 0;
    temp(498, 496) = 0;
    temp(498, 497) = 0;
    temp(498, 498) = 1;
    temp(498, 499) = 0;

    temp(499, 490) = 0;
    temp(499, 491) = 0;
    temp(499, 492) = 0;
    temp(499, 493) = 0;
    temp(499, 494) = 0;
    temp(499, 495) = 0;
    temp(499, 496) = 0;
    temp(499, 497) = 0;
    temp(499, 498) = 0;
    temp(499, 499) = 1;
    

//Block01    
    for(int i = 0; i < 49; i++)
    {
        temp(i*10 + 5, i*10 + 0) = 
        temp(i*10 + 5, i*10 + 1) = 
        temp(i*10 + 5, i*10 + 2) = 
        temp(i*10 + 5, i*10 + 3) = 
        temp(i*10 + 5, i*10 + 4) = 
        temp(i*10 + 5, i*10 + 5) = 
        temp(i*10 + 5, i*10 + 6) = 
        temp(i*10 + 5, i*10 + 7) = 
        temp(i*10 + 5, i*10 + 8) = 
        temp(i*10 + 5, i*10 + 9) = 
        
        temp(i*10 + 6, i*10 + 0) = 
        temp(i*10 + 6, i*10 + 1) = 
        temp(i*10 + 6, i*10 + 2) = 
        temp(i*10 + 6, i*10 + 3) = 
        temp(i*10 + 6, i*10 + 4) = 
        temp(i*10 + 6, i*10 + 5) = 
        temp(i*10 + 6, i*10 + 6) = 
        temp(i*10 + 6, i*10 + 7) = 
        temp(i*10 + 6, i*10 + 8) = 
        temp(i*10 + 6, i*10 + 9) = 
        
        temp(i*10 + 7, i*10 + 0) = 
        temp(i*10 + 7, i*10 + 1) = 
        temp(i*10 + 7, i*10 + 2) = 
        temp(i*10 + 7, i*10 + 3) = 
        temp(i*10 + 7, i*10 + 4) = 
        temp(i*10 + 7, i*10 + 5) = 
        temp(i*10 + 7, i*10 + 6) = 
        temp(i*10 + 7, i*10 + 7) = 
        temp(i*10 + 7, i*10 + 8) = 
        temp(i*10 + 7, i*10 + 9) = 
        
        temp(i*10 + 8, i*10 + 0) = 
        temp(i*10 + 8, i*10 + 1) = 
        temp(i*10 + 8, i*10 + 2) = 
        temp(i*10 + 8, i*10 + 3) = 
        temp(i*10 + 8, i*10 + 4) = 
        temp(i*10 + 8, i*10 + 5) = 
        temp(i*10 + 8, i*10 + 6) = 
        temp(i*10 + 8, i*10 + 7) = 
        temp(i*10 + 8, i*10 + 8) = 
        temp(i*10 + 8, i*10 + 9) = 
        
        temp(i*10 + 9, i*10 + 0) = 
        temp(i*10 + 9, i*10 + 1) = 
        temp(i*10 + 9, i*10 + 2) = 
        temp(i*10 + 9, i*10 + 3) = 
        temp(i*10 + 9, i*10 + 4) = 
        temp(i*10 + 9, i*10 + 5) = 
        temp(i*10 + 9, i*10 + 6) = 
        temp(i*10 + 9, i*10 + 7) = 
        temp(i*10 + 9, i*10 + 8) = 
        temp(i*10 + 9, i*10 + 9) = 
        
        temp(i*10 + 10, i*10 + 0) = 
        temp(i*10 + 10, i*10 + 1) = 
        temp(i*10 + 10, i*10 + 2) = 
        temp(i*10 + 10, i*10 + 3) = 
        temp(i*10 + 10, i*10 + 4) = 
        temp(i*10 + 10, i*10 + 5) = 
        temp(i*10 + 10, i*10 + 6) = 
        temp(i*10 + 10, i*10 + 7) = 
        temp(i*10 + 10, i*10 + 8) = 
        temp(i*10 + 10, i*10 + 9) = 
        
        temp(i*10 + 11, i*10 + 0) = 
        temp(i*10 + 11, i*10 + 1) = 
        temp(i*10 + 11, i*10 + 2) = 
        temp(i*10 + 11, i*10 + 3) = 
        temp(i*10 + 11, i*10 + 4) = 
        temp(i*10 + 11, i*10 + 5) = 
        temp(i*10 + 11, i*10 + 6) = 
        temp(i*10 + 11, i*10 + 7) = 
        temp(i*10 + 11, i*10 + 8) = 
        temp(i*10 + 11, i*10 + 9) = 
        
        temp(i*10 + 12, i*10 + 0) = 
        temp(i*10 + 12, i*10 + 1) = 
        temp(i*10 + 12, i*10 + 2) = 
        temp(i*10 + 12, i*10 + 3) = 
        temp(i*10 + 12, i*10 + 4) = 
        temp(i*10 + 12, i*10 + 5) = 
        temp(i*10 + 12, i*10 + 6) = 
        temp(i*10 + 12, i*10 + 7) = 
        temp(i*10 + 12, i*10 + 8) = 
        temp(i*10 + 12, i*10 + 9) = 
        
        temp(i*10 + 13, i*10 + 0) = 
        temp(i*10 + 13, i*10 + 1) = 
        temp(i*10 + 13, i*10 + 2) = 
        temp(i*10 + 13, i*10 + 3) = 
        temp(i*10 + 13, i*10 + 4) = 
        temp(i*10 + 13, i*10 + 5) = 
        temp(i*10 + 13, i*10 + 6) = 
        temp(i*10 + 13, i*10 + 7) = 
        temp(i*10 + 13, i*10 + 8) = 
        temp(i*10 + 13, i*10 + 9) = 
        
        temp(i*10 + 14, i*10 + 0) = 
        temp(i*10 + 14, i*10 + 1) = 
        temp(i*10 + 14, i*10 + 2) = 
        temp(i*10 + 14, i*10 + 3) = 
        temp(i*10 + 14, i*10 + 4) = 
        temp(i*10 + 14, i*10 + 5) = 
        temp(i*10 + 14, i*10 + 6) = 
        temp(i*10 + 14, i*10 + 7) = 
        temp(i*10 + 14, i*10 + 8) = 
        temp(i*10 + 14, i*10 + 9) = 


    }

    //Block02
    for(int i = 0; i < 49; i++)
    {
        temp(i*10 + 5, i*10 + 10) = 
        temp(i*10 + 5, i*10 + 11) = 
        temp(i*10 + 5, i*10 + 12) = 
        temp(i*10 + 5, i*10 + 13) = 
        temp(i*10 + 5, i*10 + 14) = 
        temp(i*10 + 5, i*10 + 15) = 
        temp(i*10 + 5, i*10 + 16) = 
        temp(i*10 + 5, i*10 + 17) = 
        temp(i*10 + 5, i*10 + 18) = 
        temp(i*10 + 5, i*10 + 19) = 
        
        temp(i*10 + 6, i*10 + 10) = 
        temp(i*10 + 6, i*10 + 11) = 
        temp(i*10 + 6, i*10 + 12) = 
        temp(i*10 + 6, i*10 + 13) = 
        temp(i*10 + 6, i*10 + 14) = 
        temp(i*10 + 6, i*10 + 15) = 
        temp(i*10 + 6, i*10 + 16) = 
        temp(i*10 + 6, i*10 + 17) = 
        temp(i*10 + 6, i*10 + 18) = 
        temp(i*10 + 6, i*10 + 19) = 
        
        temp(i*10 + 7, i*10 + 10) = 
        temp(i*10 + 7, i*10 + 11) = 
        temp(i*10 + 7, i*10 + 12) = 
        temp(i*10 + 7, i*10 + 13) = 
        temp(i*10 + 7, i*10 + 14) = 
        temp(i*10 + 7, i*10 + 15) = 
        temp(i*10 + 7, i*10 + 16) = 
        temp(i*10 + 7, i*10 + 17) = 
        temp(i*10 + 7, i*10 + 18) = 
        temp(i*10 + 7, i*10 + 19) = 
        
        temp(i*10 + 8, i*10 + 10) = 
        temp(i*10 + 8, i*10 + 11) = 
        temp(i*10 + 8, i*10 + 12) = 
        temp(i*10 + 8, i*10 + 13) = 
        temp(i*10 + 8, i*10 + 14) = 
        temp(i*10 + 8, i*10 + 15) = 
        temp(i*10 + 8, i*10 + 16) = 
        temp(i*10 + 8, i*10 + 17) = 
        temp(i*10 + 8, i*10 + 18) = 
        temp(i*10 + 8, i*10 + 19) = 
        
        temp(i*10 + 9, i*10 + 10) = 
        temp(i*10 + 9, i*10 + 11) = 
        temp(i*10 + 9, i*10 + 12) = 
        temp(i*10 + 9, i*10 + 13) = 
        temp(i*10 + 9, i*10 + 14) = 
        temp(i*10 + 9, i*10 + 15) = 
        temp(i*10 + 9, i*10 + 16) = 
        temp(i*10 + 9, i*10 + 17) = 
        temp(i*10 + 9, i*10 + 18) = 
        temp(i*10 + 9, i*10 + 19) = 
        
        temp(i*10 + 10, i*10 + 10) = 
        temp(i*10 + 10, i*10 + 11) = 
        temp(i*10 + 10, i*10 + 12) = 
        temp(i*10 + 10, i*10 + 13) = 
        temp(i*10 + 10, i*10 + 14) = 
        temp(i*10 + 10, i*10 + 15) = 
        temp(i*10 + 10, i*10 + 16) = 
        temp(i*10 + 10, i*10 + 17) = 
        temp(i*10 + 10, i*10 + 18) = 
        temp(i*10 + 10, i*10 + 19) = 
        
        temp(i*10 + 11, i*10 + 10) = 
        temp(i*10 + 11, i*10 + 11) = 
        temp(i*10 + 11, i*10 + 12) = 
        temp(i*10 + 11, i*10 + 13) = 
        temp(i*10 + 11, i*10 + 14) = 
        temp(i*10 + 11, i*10 + 15) = 
        temp(i*10 + 11, i*10 + 16) = 
        temp(i*10 + 11, i*10 + 17) = 
        temp(i*10 + 11, i*10 + 18) = 
        temp(i*10 + 11, i*10 + 19) = 
        
        temp(i*10 + 12, i*10 + 10) = 
        temp(i*10 + 12, i*10 + 11) = 
        temp(i*10 + 12, i*10 + 12) = 
        temp(i*10 + 12, i*10 + 13) = 
        temp(i*10 + 12, i*10 + 14) = 
        temp(i*10 + 12, i*10 + 15) = 
        temp(i*10 + 12, i*10 + 16) = 
        temp(i*10 + 12, i*10 + 17) = 
        temp(i*10 + 12, i*10 + 18) = 
        temp(i*10 + 12, i*10 + 19) = 
        
        temp(i*10 + 13, i*10 + 10) = 
        temp(i*10 + 13, i*10 + 11) = 
        temp(i*10 + 13, i*10 + 12) = 
        temp(i*10 + 13, i*10 + 13) = 
        temp(i*10 + 13, i*10 + 14) = 
        temp(i*10 + 13, i*10 + 15) = 
        temp(i*10 + 13, i*10 + 16) = 
        temp(i*10 + 13, i*10 + 17) = 
        temp(i*10 + 13, i*10 + 18) = 
        temp(i*10 + 13, i*10 + 19) = 
        
        temp(i*10 + 14, i*10 + 10) = 
        temp(i*10 + 14, i*10 + 11) = 
        temp(i*10 + 14, i*10 + 12) = 
        temp(i*10 + 14, i*10 + 13) = 
        temp(i*10 + 14, i*10 + 14) = 
        temp(i*10 + 14, i*10 + 15) = 
        temp(i*10 + 14, i*10 + 16) = 
        temp(i*10 + 14, i*10 + 17) = 
        temp(i*10 + 14, i*10 + 18) = 
        temp(i*10 + 14, i*10 + 19) = 

    }
    //Orgin JacobianMatrix 
    return temp;
}

int Jacobian::sign(double a)//conj
{
    if(a > 0)
    {
        return 1;
    }
    else
    {
        if(a < 0)
        {
            return -1;
        }
        else
        {
            return 0;
        }
    }
}