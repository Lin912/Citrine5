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

    // Vtz = -a.ReadTopVel()[0];
    // Vtx = -a.ReadTopVel()[1];
    // Vty = -a.ReadTopVel()[2];

    Vbz = a.ReadBottomVel()[0];
    Vbx = -a.ReadBottomVel()[1];
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
    temp(0, 6) = Vty*cos(Ynew(6)) + Vtz*cos(Ynew(7))*sin(Ynew(6)) + Vtx*sin(Ynew(7))*sin(Ynew(6));   
    temp(0, 7) = Vtz*cos(Ynew(6))*sin(Ynew(7)) - Vtx*cos(Ynew(7))*cos(Ynew(6));
    temp(0, 8) = 0;
    temp(0, 9) = 0;

    temp(1, 0) = 0;
    temp(1, 1) = 1;
    temp(1, 2) = 0;
    temp(1, 3) = 0;
    temp(1, 4) = 0;
    temp(1, 5) = 0;
    temp(1, 6) = 0;
    temp(1, 7) = Vtz*cos(Ynew(7)) + Vtx*sin(Ynew(7));
    temp(1, 8) = 0;
    temp(1, 9) = 0;

    temp(2, 0) = 0;
    temp(2, 1) = 0;
    temp(2, 2) = 1;
    temp(2, 3) = 0;
    temp(2, 4) = 0;
    temp(2, 5) = 0;
    temp(2, 6) = Vty*sin(Ynew(6)) - Vtz*cos(Ynew(7))*cos(Ynew(6)) - Vtx*cos(Ynew(6))*sin(Ynew(7));
    temp(2, 7) = Vtz*sin(Ynew(7))*sin(Ynew(6)) - Vtx*cos(Ynew(7))*sin(Ynew(6));
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
    
    

    //BC49
    // temp(495, 490) = 0;
    // temp(495, 491) = 0;
    // temp(495, 492) = 0;
    // temp(495, 493) = 1;
    // temp(495, 494) = 0;
    // temp(495, 495) = 0;
    // temp(495, 496) = Gy*cos(Ynew(496)) + Gz*cos(Ynew(497))*sin(Ynew(496)) + Gx*sin(Ynew(497))*sin(Ynew(496)); 
    // temp(495, 497) = Gz*cos(Ynew(496))*sin(Ynew(497)) - Gx*cos(Ynew(497))*cos(Ynew(496));
    // temp(495, 498) = 0;
    // temp(495, 499) = 0;

    // temp(496, 490) = 0;
    // temp(496, 491) = 0;
    // temp(496, 492) = 0;
    // temp(496, 493) = 0;
    // temp(496, 494) = 1;
    // temp(496, 495) = 0;
    // temp(496, 496) = 0;
    // temp(496, 497) = Gz*cos(Ynew(497)) + Gx*sin(Ynew(497));
    // temp(496, 498) = 0;
    // temp(496, 499) = 0;

    // temp(497, 490) = 0;
    // temp(497, 491) = 0;
    // temp(497, 492) = 0;
    // temp(497, 493) = 0;
    // temp(497, 494) = 0;
    // temp(497, 495) = 1;
    // temp(497, 496) = Gy*sin(Ynew(496)) - Gz*cos(Ynew(497))*cos(Ynew(496)) - Gx*cos(Ynew(496))*sin(Ynew(497));
    // temp(497, 497) = Gz*sin(Ynew(497))*sin(Ynew(496)) - Gx*cos(Ynew(497))*sin(Ynew(496));
    // temp(497, 498) = 0;
    // temp(497, 499) = 0;

    // temp(498, 490) = 0;
    // temp(498, 491) = 0;
    // temp(498, 492) = 0;
    // temp(498, 493) = 0;
    // temp(498, 494) = 0;
    // temp(498, 495) = 0;
    // temp(498, 496) = 0;
    // temp(498, 497) = 0;
    // temp(498, 498) = 1;
    // temp(498, 499) = 0;

    // temp(499, 490) = 0;
    // temp(499, 491) = 0;
    // temp(499, 492) = 0;
    // temp(499, 493) = 0;
    // temp(499, 494) = 0;
    // temp(499, 495) = 0;
    // temp(499, 496) = 0;
    // temp(499, 497) = 0;
    // temp(499, 498) = 0;
    // temp(499, 499) = 1;    


    temp(495, 490) = 1;
    temp(495, 491) = 0;
    temp(495, 492) = 0;
    temp(495, 493) = 0;
    temp(495, 494) = 0;
    temp(495, 495) = 0;
    temp(495, 496) = Vby*cos(Ynew(496)) + Vbz*cos(Ynew(497))*sin(Ynew(496)) + Vbx*sin(Ynew(497))*sin(Ynew(496));   
    temp(495, 497) = Vbz*cos(Ynew(496))*sin(Ynew(497)) - Vbx*cos(Ynew(497))*cos(Ynew(496));
    temp(495, 498) = 0;
    temp(495, 499) = 0;

    temp(496, 490) = 0;
    temp(496, 491) = 1;
    temp(496, 492) = 0;
    temp(496, 493) = 0;
    temp(496, 494) = 0;
    temp(496, 495) = 0;
    temp(496, 496) = 0;
    temp(496, 497) = Vbz*cos(Ynew(497)) + Vbx*sin(Ynew(497));
    temp(496, 498) = 0;
    temp(496, 499) = 0;

    temp(497, 490) = 0;
    temp(497, 491) = 0;
    temp(497, 492) = 1;
    temp(497, 493) = 0;
    temp(497, 494) = 0;
    temp(497, 495) = 0;
    temp(497, 496) = Vby*sin(Ynew(496)) - Vbz*cos(Ynew(497))*cos(Ynew(496)) - Vbx*cos(Ynew(496))*sin(Ynew(497));
    temp(497, 497) = Vbz*sin(Ynew(497))*sin(Ynew(496)) - Vbx*cos(Ynew(497))*sin(Ynew(496));
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
        temp(i*10 + 5, i*10 + 0) = 2*M*deltaS - deltaS*deltaT*(0.5*Cdt*d0*pi*rho*abs(Ynew(i*10+ 0) + Vy*sin(Ynew(i*10+ 6)) - Vz*cos(Ynew(i*10+ 7))*cos(Ynew(i*10+ 6)) - Vx*cos(Ynew(i*10+ 6))*sin(Ynew(i*10+ 7)))*sqrt(Ynew(i*10+ 3)/(A*E) + 1) + 0.5*Cdt*d0*pi*rho*sign(Ynew(i*10+ 0) + Vy*sin(Ynew(i*10+ 6)) - Vz*cos(Ynew(i*10+ 7))*cos(Ynew(i*10+ 6)) - Vx*cos(Ynew(i*10+ 6))*sin(Ynew(i*10+ 7)))*sqrt(Ynew(i*10+ 3)/(A*E) + 1)*(Ynew(i*10+ 0) + Vy*sin(Ynew(i*10+ 6)) - Vz*cos(Ynew(i*10+ 7))*cos(Ynew(i*10+ 6)) - Vx*cos(Ynew(i*10+ 6))*sin(Ynew(i*10+ 7))));
        temp(i*10 + 6, i*10 + 0) = -M*deltaS*cos(Ynew(i*10+ 6))*(Yold(i*10+ 7) - Ynew(i*10+ 7));
        temp(i*10 + 7, i*10 + 0) = M*deltaS*(Yold(i*10+ 6) - Ynew(i*10+ 6));
        temp(i*10 + 10, i*10 + 0) = 2*deltaT;
        temp(i*10 + 11, i*10 + 0) = -Ynew(i*10+ 9 )*deltaS*deltaT;
        temp(i*10 + 12, i*10 + 0) = -Ynew(i*10+ 8 )*deltaS*deltaT;

        temp(i*10 + 5, i*10 + 1) = M*deltaS*cos(Ynew(i*10+ 6))*(Yold(i*10+ 7) - Ynew(i*10+ 7));
        temp(i*10 + 6, i*10 + 1) = deltaS*(2*M - 2*ma) - deltaS*deltaT*(0.5*Cdn*rho*sqrt(Ynew(i*10+ 3)/(A*E) + 1)*sqrt(pow((Ynew(i*10+ 1) - Vx*cos(Ynew(i*10+ 7)) + Vz*sin(Ynew(i*10+ 7))),2) + pow((Vy*cos(Ynew(i*10+ 6)) - Ynew(i*10+ 2) + Vz*cos(Ynew(i*10+ 7))*sin(Ynew(i*10+ 6)) + Vx*sin(Ynew(i*10+ 7))*sin(Ynew(i*10+ 6))),2)) + (0.25*Cdn*rho*sqrt(Ynew(i*10+ 3)/(A*E) + 1)*(2*Ynew(i*10+ 1) - 2*Vx*cos(Ynew(i*10+ 7)) + 2*Vz*sin(Ynew(i*10+ 7)))*(Ynew(i*10+ 1) - Vx*cos(Ynew(i*10+ 7)) + Vz*sin(Ynew(i*10+ 7))))/sqrt(pow((Ynew(i*10+ 1) - Vx*cos(Ynew(i*10+ 7)) + Vz*sin(Ynew(i*10+ 7))),2) + pow((Vy*cos(Ynew(i*10+ 6)) - Ynew(i*10+ 2) + Vz*cos(Ynew(i*10+ 7))*sin(Ynew(i*10+ 6)) + Vx*sin(Ynew(i*10+ 7))*sin(Ynew(i*10+ 6))),2)));
        temp(i*10 + 7, i*10 + 1) = M*deltaS*sin(Ynew(i*10+ 6))*(Yold(i*10+ 7) - Ynew(i*10+ 7)) + (0.25*Cdb*deltaS*deltaT*rho*sqrt(Ynew(i*10+ 3)/(A*E) + 1)*(2*Ynew(i*10+ 1) - 2*Vx*cos(Ynew(i*10+ 7)) + 2*Vz*sin(Ynew(i*10+ 7)))*(Vy*cos(Ynew(i*10+ 6)) - Ynew(i*10+ 2) + Vz*cos(Ynew(i*10+ 7))*sin(Ynew(i*10+ 6)) + Vx*sin(Ynew(i*10+ 7))*sin(Ynew(i*10+ 6))))/sqrt(pow((Ynew(i*10+ 1) - Vx*cos(Ynew(i*10+ 7)) + Vz*sin(Ynew(i*10+ 7))),2) + pow((Vy*cos(Ynew(i*10+ 6)) - Ynew(i*10+ 2) + Vz*cos(Ynew(i*10+ 7))*sin(Ynew(i*10+ 6)) + Vx*sin(Ynew(i*10+ 7))*sin(Ynew(i*10+ 6))),2));
        temp(i*10 + 10, i*10 + 1) = Ynew(i*10+ 9 )*deltaS*deltaT;
        temp(i*10 + 11, i*10 + 1) = 2*deltaT;
        temp(i*10 + 12, i*10 + 1) = -Ynew(i*10+ 9 )*deltaS*deltaT*tan(Ynew(i*10+ 6 ));

        temp(i*10 + 5, i*10 + 2) = -M*deltaS*(Yold(i*10+ 6) - Ynew(i*10+ 6));
        temp(i*10 + 6, i*10 + 2) = (0.25*Cdn*deltaS*deltaT*rho*sqrt(Ynew(i*10+ 3)/(A*E) + 1)*(Ynew(i*10+ 1) - Vx*cos(Ynew(i*10+ 7)) + Vz*sin(Ynew(i*10+ 7)))*(2*Vy*cos(Ynew(i*10+ 6)) - 2*Ynew(i*10+ 2) + 2*Vz*cos(Ynew(i*10+ 7))*sin(Ynew(i*10+ 6)) + 2*Vx*sin(Ynew(i*10+ 7))*sin(Ynew(i*10+ 6))))/sqrt(pow((Ynew(i*10+ 1) - Vx*cos(Ynew(i*10+ 7)) + Vz*sin(Ynew(i*10+ 7))),2) + pow((Vy*cos(Ynew(i*10+ 6)) - Ynew(i*10+ 2) + Vz*cos(Ynew(i*10+ 7))*sin(Ynew(i*10+ 6)) + Vx*sin(Ynew(i*10+ 7))*sin(Ynew(i*10+ 6))),2)) - M*deltaS*sin(Ynew(i*10+ 6))*(Yold(i*10+ 7) - Ynew(i*10+ 7));
        temp(i*10 + 7, i*10 + 2) = deltaS*(2*M - 2*ma) - deltaS*deltaT*(0.5*Cdb*rho*sqrt(Ynew(i*10+ 3)/(A*E) + 1)*sqrt(pow((Ynew(i*10+ 1) - Vx*cos(Ynew(i*10+ 7)) + Vz*sin(Ynew(i*10+ 7))),2) + pow((Vy*cos(Ynew(i*10+ 6)) - Ynew(i*10+ 2) + Vz*cos(Ynew(i*10+ 7))*sin(Ynew(i*10+ 6)) + Vx*sin(Ynew(i*10+ 7))*sin(Ynew(i*10+ 6))),2)) + (0.25*Cdb*rho*sqrt(Ynew(i*10+ 3)/(A*E) + 1)*(Vy*cos(Ynew(i*10+ 6)) - Ynew(i*10+ 2) + Vz*cos(Ynew(i*10+ 7))*sin(Ynew(i*10+ 6)) + Vx*sin(Ynew(i*10+ 7))*sin(Ynew(i*10+ 6)))*(2*Vy*cos(Ynew(i*10+ 6)) - 2*Ynew(i*10+ 2) + 2*Vz*cos(Ynew(i*10+ 7))*sin(Ynew(i*10+ 6)) + 2*Vx*sin(Ynew(i*10+ 7))*sin(Ynew(i*10+ 6))))/sqrt(pow((Ynew(i*10+ 1) - Vx*cos(Ynew(i*10+ 7)) + Vz*sin(Ynew(i*10+ 7))),2) + pow((Vy*cos(Ynew(i*10+ 6)) - Ynew(i*10+ 2) + Vz*cos(Ynew(i*10+ 7))*sin(Ynew(i*10+ 6)) + Vx*sin(Ynew(i*10+ 7))*sin(Ynew(i*10+ 6))),2)));
        temp(i*10 + 10, i*10 + 2) = -Ynew(i*10+ 8 )*deltaS*deltaT;
        temp(i*10 + 11, i*10 + 2) = -Ynew(i*10+ 9 )*deltaS*deltaT*tan(Ynew(i*10+ 6 ));
        temp(i*10 + 12, i*10 + 2) = -2*deltaT;

        temp(i*10 + 5, i*10 + 3) = 2*deltaT - (0.25*Cdt*d0*deltaS*deltaT*pi*rho*abs(Ynew(i*10+ 0) + Vy*sin(Ynew(i*10+ 6)) - Vz*cos(Ynew(i*10+ 7))*cos(Ynew(i*10+ 6)) - Vx*cos(Ynew(i*10+ 6))*sin(Ynew(i*10+ 7)))*(Ynew(i*10+ 0) + Vy*sin(Ynew(i*10+ 6)) - Vz*cos(Ynew(i*10+ 7))*cos(Ynew(i*10+ 6)) - Vx*cos(Ynew(i*10+ 6))*sin(Ynew(i*10+ 7))))/(A*E*sqrt(Ynew(i*10+ 3)/(A*E) + 1));
        temp(i*10 + 6, i*10 + 3) = -deltaS*deltaT*(Ynew(i*10+ 9) + (0.25*Cdn*rho*sqrt(pow((Ynew(i*10+ 1) - Vx*cos(Ynew(i*10+ 7)) + Vz*sin(Ynew(i*10+ 7))),2) + pow((Vy*cos(Ynew(i*10+ 6)) - Ynew(i*10+ 2) + Vz*cos(Ynew(i*10+ 7))*sin(Ynew(i*10+ 6)) + Vx*sin(Ynew(i*10+ 7))*sin(Ynew(i*10+ 6))),2))*(Ynew(i*10+ 1) - Vx*cos(Ynew(i*10+ 7)) + Vz*sin(Ynew(i*10+ 7))))/(A*E*sqrt(Ynew(i*10+ 3)/(A*E) + 1)));
        temp(i*10 + 7, i*10 + 3) = deltaS*deltaT*(Ynew(i*10+ 8) + (0.25*Cdb*rho*sqrt(pow((Ynew(i*10+ 1) - Vx*cos(Ynew(i*10+ 7)) + Vz*sin(Ynew(i*10+ 7))),2) + pow((Vy*cos(Ynew(i*10+ 6)) - Ynew(i*10+ 2) + Vz*cos(Ynew(i*10+ 7))*sin(Ynew(i*10+ 6)) + Vx*sin(Ynew(i*10+ 7))*sin(Ynew(i*10+ 6))),2))*(Vy*cos(Ynew(i*10+ 6)) - Ynew(i*10+ 2) + Vz*cos(Ynew(i*10+ 7))*sin(Ynew(i*10+ 6)) + Vx*sin(Ynew(i*10+ 7))*sin(Ynew(i*10+ 6))))/(A*E*sqrt(Ynew(i*10+ 3)/(A*E) + 1)));
        temp(i*10 + 8, i*10 + 3) = -(3*Ynew(i*10+ 5 )*deltaS*deltaT*pow((Ynew(i*10+ 3 )/(A*E) + 1),2))/(A*E);
        temp(i*10 + 9, i*10 + 3) = (3*Ynew(i*10+ 4 )*deltaS*deltaT*pow((Ynew(i*10+ 3 )/(A*E) + 1),2))/(A*E);
        temp(i*10 + 10, i*10 + 3) = (2*deltaS)/(A*E);
        temp(i*10 + 11, i*10 + 3) = -(deltaS*cos(Ynew(i*10+ 6 ))*(Yold(i*10+ 7 ) - Ynew(i*10+ 7 )))/(A*E);
        temp(i*10 + 12, i*10 + 3) = -(deltaS*(Yold(i*10+ 6 ) - Ynew(i*10+ 6 )))/(A*E);

        temp(i*10 + 5, i*10 + 4) = Ynew(i*10+ 9)*deltaS*deltaT;
        temp(i*10 + 6, i*10 + 4) = 2*deltaT;
        temp(i*10 + 7, i*10 + 4) = -Ynew(i*10+ 9)*deltaS*deltaT*tan(Ynew(i*10+ 6));
        temp(i*10 + 9, i*10 + 4) = deltaS*deltaT*pow((Ynew(i*10+ 3 )/(A*E) + 1),3);

        temp(i*10 + 5, i*10 + 5) = -Ynew(i*10+ 8)*deltaS*deltaT;
        temp(i*10 + 6, i*10 + 5) = -Ynew(i*10+ 9)*deltaS*deltaT*tan(Ynew(i*10+ 6));
        temp(i*10 + 7, i*10 + 5) = 2*deltaT;
        temp(i*10 + 8, i*10 + 5) = -deltaS*deltaT*pow((Ynew(i*10+ 3 )/(A*E) + 1),3);

        temp(i*10 + 5, i*10 + 6) = deltaS*(M*Ynew(i*10+ 2) + M*Yold(i*10+ 2) - M*Ynew(i*10+ 1)*sin(Ynew(i*10+ 6))*(Yold(i*10+ 7) - Ynew(i*10+ 7))) - deltaS*deltaT*(0.5*Cdt*d0*pi*rho*abs(Ynew(i*10+ 0) + Vy*sin(Ynew(i*10+ 6)) - Vz*cos(Ynew(i*10+ 7))*cos(Ynew(i*10+ 6)) - Vx*cos(Ynew(i*10+ 6))*sin(Ynew(i*10+ 7)))*sqrt(Ynew(i*10+ 3)/(A*E) + 1)*(Vy*cos(Ynew(i*10+ 6)) + Vz*cos(Ynew(i*10+ 7))*sin(Ynew(i*10+ 6)) + Vx*sin(Ynew(i*10+ 7))*sin(Ynew(i*10+ 6))) - w0*sin(Ynew(i*10+ 7))*sin(Ynew(i*10+ 6)) + 0.5*Cdt*d0*pi*rho*sign(Ynew(i*10+ 0) + Vy*sin(Ynew(i*10+ 6)) - Vz*cos(Ynew(i*10+ 7))*cos(Ynew(i*10+ 6)) - Vx*cos(Ynew(i*10+ 6))*sin(Ynew(i*10+ 7)))*sqrt(Ynew(i*10+ 3)/(A*E) + 1)*(Vy*cos(Ynew(i*10+ 6)) + Vz*cos(Ynew(i*10+ 7))*sin(Ynew(i*10+ 6)) + Vx*sin(Ynew(i*10+ 7))*sin(Ynew(i*10+ 6)))*(Ynew(i*10+ 0) + Vy*sin(Ynew(i*10+ 6)) - Vz*cos(Ynew(i*10+ 7))*cos(Ynew(i*10+ 6)) - Vx*cos(Ynew(i*10+ 6))*sin(Ynew(i*10+ 7))));
        temp(i*10 + 6, i*10 + 6) = - deltaS*(M*Ynew(i*10+ 2)*cos(Ynew(i*10+ 6)) - M*Ynew(i*10+ 0)*sin(Ynew(i*10+ 6)))*(Yold(i*10+ 7) - Ynew(i*10+ 7)) - deltaS*deltaT*(Ynew(i*10+ 9)*Ynew(i*10+ 5)*(pow(tan(Ynew(i*10+ 6)),2) + 1) + (0.5*Cdn*rho*sqrt(Ynew(i*10+ 3)/(A*E) + 1)*(Vz*cos(Ynew(i*10+ 7))*cos(Ynew(i*10+ 6)) - Vy*sin(Ynew(i*10+ 6)) + Vx*cos(Ynew(i*10+ 6))*sin(Ynew(i*10+ 7)))*(Ynew(i*10+ 1) - Vx*cos(Ynew(i*10+ 7)) + Vz*sin(Ynew(i*10+ 7)))*(Vy*cos(Ynew(i*10+ 6)) - Ynew(i*10+ 2) + Vz*cos(Ynew(i*10+ 7))*sin(Ynew(i*10+ 6)) + Vx*sin(Ynew(i*10+ 7))*sin(Ynew(i*10+ 6))))/sqrt(pow((Ynew(i*10+ 1) - Vx*cos(Ynew(i*10+ 7)) + Vz*sin(Ynew(i*10+ 7))),2) + pow((Vy*cos(Ynew(i*10+ 6)) - Ynew(i*10+ 2) + Vz*cos(Ynew(i*10+ 7))*sin(Ynew(i*10+ 6)) + Vx*sin(Ynew(i*10+ 7))*sin(Ynew(i*10+ 6))),2)));
        temp(i*10 + 7, i*10 + 6) = deltaS*((Yold(i*10+ 7) - Ynew(i*10+ 7))*(M*Ynew(i*10+ 1)*cos(Ynew(i*10+ 6)) - Vx*ma*cos(Ynew(i*10+ 7))*cos(Ynew(i*10+ 6)) + Vz*ma*cos(Ynew(i*10+ 6))*sin(Ynew(i*10+ 7))) - M*Yold(i*10+ 0) - M*Ynew(i*10+ 0) + (Yold(i*10+ 6) - Ynew(i*10+ 6))*(Vy*ma*cos(Ynew(i*10+ 6)) + Vz*ma*cos(Ynew(i*10+ 7))*sin(Ynew(i*10+ 6)) + Vx*ma*sin(Ynew(i*10+ 7))*sin(Ynew(i*10+ 6))) - Vy*ma*sin(Yold(i*10+ 6)) - Vy*ma*sin(Ynew(i*10+ 6)) + Vz*ma*cos(Yold(i*10+ 7))*cos(Yold(i*10+ 6)) + Vz*ma*cos(Ynew(i*10+ 7))*cos(Ynew(i*10+ 6)) + Vx*ma*cos(Yold(i*10+ 6))*sin(Yold(i*10+ 7)) + Vx*ma*cos(Ynew(i*10+ 6))*sin(Ynew(i*10+ 7))) - deltaS*deltaT*(w0*cos(Ynew(i*10+ 6))*sin(Ynew(i*10+ 7)) + Ynew(i*10+ 9)*Ynew(i*10+ 4)*(pow(tan(Ynew(i*10+ 6)),2) + 1) - 0.5*Cdb*rho*sqrt(Ynew(i*10+ 3)/(A*E) + 1)*sqrt(pow((Ynew(i*10+ 1) - Vx*cos(Ynew(i*10+ 7)) + Vz*sin(Ynew(i*10+ 7))),2) + pow((Vy*cos(Ynew(i*10+ 6)) - Ynew(i*10+ 2) + Vz*cos(Ynew(i*10+ 7))*sin(Ynew(i*10+ 6)) + Vx*sin(Ynew(i*10+ 7))*sin(Ynew(i*10+ 6))),2))*(Vz*cos(Ynew(i*10+ 7))*cos(Ynew(i*10+ 6)) - Vy*sin(Ynew(i*10+ 6)) + Vx*cos(Ynew(i*10+ 6))*sin(Ynew(i*10+ 7))) - (0.5*Cdb*rho*sqrt(Ynew(i*10+ 3)/(A*E) + 1)*(Vz*cos(Ynew(i*10+ 7))*cos(Ynew(i*10+ 6)) - Vy*sin(Ynew(i*10+ 6)) + Vx*cos(Ynew(i*10+ 6))*sin(Ynew(i*10+ 7)))*pow((Vy*cos(Ynew(i*10+ 6)) - Ynew(i*10+ 2) + Vz*cos(Ynew(i*10+ 7))*sin(Ynew(i*10+ 6)) + Vx*sin(Ynew(i*10+ 7))*sin(Ynew(i*10+ 6))),2))/sqrt(pow((Ynew(i*10+ 1) - Vx*cos(Ynew(i*10+ 7)) + Vz*sin(Ynew(i*10+ 7))),2) + pow((Vy*cos(Ynew(i*10+ 6)) - Ynew(i*10+ 2) + Vz*cos(Ynew(i*10+ 7))*sin(Ynew(i*10+ 6)) + Vx*sin(Ynew(i*10+ 7))*sin(Ynew(i*10+ 6))),2)));
        temp(i*10 + 8, i*10 + 6) = E*I*pow(Ynew(i*10+ 9 ),2)*deltaS*deltaT*(pow(tan(Ynew(i*10+ 6 )),2) + 1);
        temp(i*10 + 9, i*10 + 6) = -E*I*Ynew(i*10+ 8 )*Ynew(i*10+ 9 )*deltaS*deltaT*(pow(tan(Ynew(i*10+ 6 )),2) + 1);
        temp(i*10 + 11, i*10 + 6) = deltaS*sin(Ynew(i*10+ 6 ))*(Ynew(i*10+ 3 )/(A*E) + 1)*(Yold(i*10+ 7 ) - Ynew(i*10+ 7 )) - Ynew(i*10+ 9 )*deltaS*deltaT*Ynew(i*10+ 2 )*(pow(tan(Ynew(i*10+ 6 )),2) + 1);
        temp(i*10 + 12, i*10 + 6) = deltaS*(Ynew(i*10+ 3 )/(A*E) + Yold(i*10+ 3 )/(A*E) + 2) - Ynew(i*10+ 9 )*deltaS*deltaT*Ynew(i*10+ 1 )*(pow(tan(Ynew(i*10+ 6 )),2) + 1);
        temp(i*10 + 13, i*10 + 6) = 2*deltaT;
        temp(i*10 + 14, i*10 + 6) = -deltaT*sin(Ynew(i*10+ 6 ))*(Ynew(i*10+ 7 ) - Ynew(i*10+ 17 ));

        temp(i*10 + 5, i*10 + 7) = deltaS*deltaT*(0.5*Cdt*d0*pi*rho*abs(Ynew(i*10+ 0) + Vy*sin(Ynew(i*10+ 6)) - Vz*cos(Ynew(i*10+ 7))*cos(Ynew(i*10+ 6)) - Vx*cos(Ynew(i*10+ 6))*sin(Ynew(i*10+ 7)))*sqrt(Ynew(i*10+ 3)/(A*E) + 1)*(Vx*cos(Ynew(i*10+ 7))*cos(Ynew(i*10+ 6)) - Vz*cos(Ynew(i*10+ 6))*sin(Ynew(i*10+ 7))) - w0*cos(Ynew(i*10+ 7))*cos(Ynew(i*10+ 6)) + 0.5*Cdt*d0*pi*rho*sign(Ynew(i*10+ 0) + Vy*sin(Ynew(i*10+ 6)) - Vz*cos(Ynew(i*10+ 7))*cos(Ynew(i*10+ 6)) - Vx*cos(Ynew(i*10+ 6))*sin(Ynew(i*10+ 7)))*sqrt(Ynew(i*10+ 3)/(A*E) + 1)*(Vx*cos(Ynew(i*10+ 7))*cos(Ynew(i*10+ 6)) - Vz*cos(Ynew(i*10+ 6))*sin(Ynew(i*10+ 7)))*(Ynew(i*10+ 0) + Vy*sin(Ynew(i*10+ 6)) - Vz*cos(Ynew(i*10+ 7))*cos(Ynew(i*10+ 6)) - Vx*cos(Ynew(i*10+ 6))*sin(Ynew(i*10+ 7)))) - deltaS*(M*Ynew(i*10+ 1)*cos(Ynew(i*10+ 6)) + M*Yold(i*10+ 1)*cos(Yold(i*10+ 6)));
        temp(i*10 + 6, i*10 + 7) = deltaS*((Vx*ma*cos(Ynew(i*10+ 7)) - Vz*ma*sin(Ynew(i*10+ 7)))*(Yold(i*10+ 7) - Ynew(i*10+ 7)) - Vz*ma*cos(Yold(i*10+ 7)) - Vz*ma*cos(Ynew(i*10+ 7)) + M*Ynew(i*10+ 0)*cos(Ynew(i*10+ 6)) + M*Yold(i*10+ 0)*cos(Yold(i*10+ 6)) - Vx*ma*sin(Yold(i*10+ 7)) - Vx*ma*sin(Ynew(i*10+ 7)) + M*Ynew(i*10+ 2)*sin(Ynew(i*10+ 6)) + M*Yold(i*10+ 2)*sin(Yold(i*10+ 6))) - deltaS*deltaT*(0.5*Cdn*rho*(Vz*cos(Ynew(i*10+ 7)) + Vx*sin(Ynew(i*10+ 7)))*sqrt(Ynew(i*10+ 3)/(A*E) + 1)*sqrt(pow((Ynew(i*10+ 1) - Vx*cos(Ynew(i*10+ 7)) + Vz*sin(Ynew(i*10+ 7))),2) + pow((Vy*cos(Ynew(i*10+ 6)) - Ynew(i*10+ 2) + Vz*cos(Ynew(i*10+ 7))*sin(Ynew(i*10+ 6)) + Vx*sin(Ynew(i*10+ 7))*sin(Ynew(i*10+ 6))),2)) - w0*sin(Ynew(i*10+ 7)) + (0.25*Cdn*rho*(2*(Vx*cos(Ynew(i*10+ 7))*sin(Ynew(i*10+ 6)) - Vz*sin(Ynew(i*10+ 7))*sin(Ynew(i*10+ 6)))*(Vy*cos(Ynew(i*10+ 6)) - Ynew(i*10+ 2) + Vz*cos(Ynew(i*10+ 7))*sin(Ynew(i*10+ 6)) + Vx*sin(Ynew(i*10+ 7))*sin(Ynew(i*10+ 6))) + 2*(Vz*cos(Ynew(i*10+ 7)) + Vx*sin(Ynew(i*10+ 7)))*(Ynew(i*10+ 1) - Vx*cos(Ynew(i*10+ 7)) + Vz*sin(Ynew(i*10+ 7))))*sqrt(Ynew(i*10+ 3)/(A*E) + 1)*(Ynew(i*10+ 1) - Vx*cos(Ynew(i*10+ 7)) + Vz*sin(Ynew(i*10+ 7))))/sqrt(pow((Ynew(i*10+ 1) - Vx*cos(Ynew(i*10+ 7)) + Vz*sin(Ynew(i*10+ 7))),2) + pow((Vy*cos(Ynew(i*10+ 6)) - Ynew(i*10+ 2) + Vz*cos(Ynew(i*10+ 7))*sin(Ynew(i*10+ 6)) + Vx*sin(Ynew(i*10+ 7))*sin(Ynew(i*10+ 6))),2)));
        temp(i*10 + 7, i*10 + 7) = deltaS*deltaT*(0.5*Cdb*rho*sqrt(Ynew(i*10+ 3)/(A*E) + 1)*(Vx*cos(Ynew(i*10+ 7))*sin(Ynew(i*10+ 6)) - Vz*sin(Ynew(i*10+ 7))*sin(Ynew(i*10+ 6)))*sqrt(pow((Ynew(i*10+ 1) - Vx*cos(Ynew(i*10+ 7)) + Vz*sin(Ynew(i*10+ 7))),2) + pow((Vy*cos(Ynew(i*10+ 6)) - Ynew(i*10+ 2) + Vz*cos(Ynew(i*10+ 7))*sin(Ynew(i*10+ 6)) + Vx*sin(Ynew(i*10+ 7))*sin(Ynew(i*10+ 6))),2)) - w0*cos(Ynew(i*10+ 7))*sin(Ynew(i*10+ 6)) + (0.25*Cdb*rho*(2*(Vx*cos(Ynew(i*10+ 7))*sin(Ynew(i*10+ 6)) - Vz*sin(Ynew(i*10+ 7))*sin(Ynew(i*10+ 6)))*(Vy*cos(Ynew(i*10+ 6)) - Ynew(i*10+ 2) + Vz*cos(Ynew(i*10+ 7))*sin(Ynew(i*10+ 6)) + Vx*sin(Ynew(i*10+ 7))*sin(Ynew(i*10+ 6))) + 2*(Vz*cos(Ynew(i*10+ 7)) + Vx*sin(Ynew(i*10+ 7)))*(Ynew(i*10+ 1) - Vx*cos(Ynew(i*10+ 7)) + Vz*sin(Ynew(i*10+ 7))))*sqrt(Ynew(i*10+ 3)/(A*E) + 1)*(Vy*cos(Ynew(i*10+ 6)) - Ynew(i*10+ 2) + Vz*cos(Ynew(i*10+ 7))*sin(Ynew(i*10+ 6)) + Vx*sin(Ynew(i*10+ 7))*sin(Ynew(i*10+ 6))))/sqrt(pow((Ynew(i*10+ 1) - Vx*cos(Ynew(i*10+ 7)) + Vz*sin(Ynew(i*10+ 7))),2) + pow((Vy*cos(Ynew(i*10+ 6)) - Ynew(i*10+ 2) + Vz*cos(Ynew(i*10+ 7))*sin(Ynew(i*10+ 6)) + Vx*sin(Ynew(i*10+ 7))*sin(Ynew(i*10+ 6))),2))) - deltaS*((Yold(i*10+ 6) - Ynew(i*10+ 6))*(Vx*ma*cos(Ynew(i*10+ 7))*cos(Ynew(i*10+ 6)) - Vz*ma*cos(Ynew(i*10+ 6))*sin(Ynew(i*10+ 7))) - (Yold(i*10+ 7) - Ynew(i*10+ 7))*(Vz*ma*cos(Ynew(i*10+ 7))*sin(Ynew(i*10+ 6)) + Vx*ma*sin(Ynew(i*10+ 7))*sin(Ynew(i*10+ 6))) + M*Ynew(i*10+ 1)*sin(Ynew(i*10+ 6)) + M*Yold(i*10+ 1)*sin(Yold(i*10+ 6)) - Vx*ma*cos(Yold(i*10+ 7))*sin(Yold(i*10+ 6)) - Vx*ma*cos(Ynew(i*10+ 7))*sin(Ynew(i*10+ 6)) + Vz*ma*sin(Yold(i*10+ 7))*sin(Yold(i*10+ 6)) + Vz*ma*sin(Ynew(i*10+ 7))*sin(Ynew(i*10+ 6)));
        temp(i*10 + 11, i*10 + 7) = deltaS*(cos(Ynew(i*10+ 6 ))*(Ynew(i*10+ 3 )/(A*E) + 1) + cos(Yold(i*10+ 6 ))*(Yold(i*10+ 3 )/(A*E) + 1));
        temp(i*10 + 14, i*10 + 7) = deltaT*(cos(Ynew(i*10+ 6 )) + cos(Ynew(i*10+ 16 )));
        
        temp(i*10 + 5, i*10 + 8) = -Ynew(i*10+ 5)*deltaS*deltaT;
        temp(i*10 + 7, i*10 + 8) = Ynew(i*10+ 3)*deltaS*deltaT;
        temp(i*10 + 8, i*10 + 8) = -2*E*I*deltaT;
        temp(i*10 + 9, i*10 + 8) = -E*I*Ynew(i*10+ 9 )*deltaS*deltaT*tan(Ynew(i*10+ 6 ));
        temp(i*10 + 10, i*10 + 8) = -deltaS*deltaT*Ynew(i*10+ 2 );
        temp(i*10 + 12, i*10 + 8) = -deltaS*deltaT*Ynew(i*10+ 0 );
        temp(i*10 + 13, i*10 + 8) = deltaS*deltaT;
        
        temp(i*10 + 5, i*10 + 9) = Ynew(i*10+ 4)*deltaS*deltaT;
        temp(i*10 + 6, i*10 + 9) = -deltaS*deltaT*(Ynew(i*10+ 3) + Ynew(i*10+ 5)*tan(Ynew(i*10+ 6)));
        temp(i*10 + 7, i*10 + 9) = -Ynew(i*10+ 4)*deltaS*deltaT*tan(Ynew(i*10+ 6));
        temp(i*10 + 8, i*10 + 9) = 2*E*I*Ynew(i*10+ 9 )*deltaS*deltaT*tan(Ynew(i*10+ 6 ));
        temp(i*10 + 9, i*10 + 9) = - 2*E*I*deltaT - E*I*Ynew(i*10+ 8 )*deltaS*deltaT*tan(Ynew(i*10+ 6 ));
        temp(i*10 + 10, i*10 + 9) = deltaS*deltaT*Ynew(i*10+ 1 );
        temp(i*10 + 11, i*10 + 9) = -deltaS*deltaT*(Ynew(i*10+ 0 ) + Ynew(i*10+ 2 )*tan(Ynew(i*10+ 6 )));
        temp(i*10 + 12, i*10 + 9) = -deltaS*deltaT*Ynew(i*10+ 1 )*tan(Ynew(i*10+ 6 ));
        temp(i*10 + 14, i*10 + 9) = deltaS*deltaT;

    }

    //Block02
    for(int i = 0; i < 49; i++)
    {
        temp(i*10 + 5, i*10 + 10) = 2*M*deltaS - deltaS*deltaT*(0.5*Cdt*d0*pi*rho*abs(Ynew(i*10+ 10) + Vy*sin(Ynew(i*10+ 16)) - Vz*cos(Ynew(i*10+ 17))*cos(Ynew(i*10+ 16)) - Vx*cos(Ynew(i*10+ 16))*sin(Ynew(i*10+ 17)))*sqrt(Ynew(i*10+ 13)/(A*E) + 1) + 0.5*Cdt*d0*pi*rho*sign(Ynew(i*10+ 10) + Vy*sin(Ynew(i*10+ 16)) - Vz*cos(Ynew(i*10+ 17))*cos(Ynew(i*10+ 16)) - Vx*cos(Ynew(i*10+ 16))*sin(Ynew(i*10+ 17)))*sqrt(Ynew(i*10+ 13)/(A*E) + 1)*(Ynew(i*10+ 10) + Vy*sin(Ynew(i*10+ 16)) - Vz*cos(Ynew(i*10+ 17))*cos(Ynew(i*10+ 16)) - Vx*cos(Ynew(i*10+ 16))*sin(Ynew(i*10+ 17))));
        temp(i*10 + 6, i*10 + 10) = -M*deltaS*cos(Ynew(i*10+ 16))*(Yold(i*10+ 17) - Ynew(i*10+ 17));
        temp(i*10 + 7, i*10 + 10) = M*deltaS*(Yold(i*10+ 16) - Ynew(i*10+ 16));
        temp(i*10 + 10, i*10 + 10) = -2*deltaT;
        temp(i*10 + 11, i*10 + 10) = -Ynew(i*10+ 19 )*deltaS*deltaT;
        temp(i*10 + 12, i*10 + 10) = -Ynew(i*10+ 18 )*deltaS*deltaT;
       
        temp(i*10 + 5, i*10 + 11) = M*deltaS*cos(Ynew(i*10+ 16))*(Yold(i*10+ 17) - Ynew(i*10+ 17));
        temp(i*10 + 6, i*10 + 11) = deltaS*(2*M - 2*ma) - deltaS*deltaT*(0.5*Cdn*rho*sqrt(Ynew(i*10+ 13)/(A*E) + 1)*sqrt(pow((Ynew(i*10+ 11) - Vx*cos(Ynew(i*10+ 17)) + Vz*sin(Ynew(i*10+ 17))),2) + pow((Vy*cos(Ynew(i*10+ 16)) - Ynew(i*10+ 12) + Vz*cos(Ynew(i*10+ 17))*sin(Ynew(i*10+ 16)) + Vx*sin(Ynew(i*10+ 17))*sin(Ynew(i*10+ 16))),2)) + (0.25*Cdn*rho*sqrt(Ynew(i*10+ 13)/(A*E) + 1)*(2*Ynew(i*10+ 11) - 2*Vx*cos(Ynew(i*10+ 17)) + 2*Vz*sin(Ynew(i*10+ 17)))*(Ynew(i*10+ 11) - Vx*cos(Ynew(i*10+ 17)) + Vz*sin(Ynew(i*10+ 17))))/sqrt(pow((Ynew(i*10+ 11) - Vx*cos(Ynew(i*10+ 17)) + Vz*sin(Ynew(i*10+ 17))),2) + pow((Vy*cos(Ynew(i*10+ 16)) - Ynew(i*10+ 12) + Vz*cos(Ynew(i*10+ 17))*sin(Ynew(i*10+ 16)) + Vx*sin(Ynew(i*10+ 17))*sin(Ynew(i*10+ 16))),2)));
        temp(i*10 + 7, i*10 + 11) = M*deltaS*sin(Ynew(i*10+ 16))*(Yold(i*10+ 17) - Ynew(i*10+ 17)) + (0.25*Cdb*deltaS*deltaT*rho*sqrt(Ynew(i*10+ 13)/(A*E) + 1)*(2*Ynew(i*10+ 11) - 2*Vx*cos(Ynew(i*10+ 17)) + 2*Vz*sin(Ynew(i*10+ 17)))*(Vy*cos(Ynew(i*10+ 16)) - Ynew(i*10+ 12) + Vz*cos(Ynew(i*10+ 17))*sin(Ynew(i*10+ 16)) + Vx*sin(Ynew(i*10+ 17))*sin(Ynew(i*10+ 16))))/sqrt(pow((Ynew(i*10+ 11) - Vx*cos(Ynew(i*10+ 17)) + Vz*sin(Ynew(i*10+ 17))),2) + pow((Vy*cos(Ynew(i*10+ 16)) - Ynew(i*10+ 12) + Vz*cos(Ynew(i*10+ 17))*sin(Ynew(i*10+ 16)) + Vx*sin(Ynew(i*10+ 17))*sin(Ynew(i*10+ 16))),2));
        temp(i*10 + 10, i*10 + 11) = Ynew(i*10+ 19)*deltaS*deltaT;
        temp(i*10 + 11, i*10 + 11) = -2*deltaT;
        temp(i*10 + 12, i*10 + 11) = -Ynew(i*10+ 19)*deltaS*deltaT*tan(Ynew(i*10+ 16));

        temp(i*10 + 5, i*10 + 12) = -M*deltaS*(Yold(i*10+ 16) - Ynew(i*10+ 16));
        temp(i*10 + 6, i*10 + 12) = (0.25*Cdn*deltaS*deltaT*rho*sqrt(Ynew(i*10+ 13)/(A*E) + 1)*(Ynew(i*10+ 11) - Vx*cos(Ynew(i*10+ 17)) + Vz*sin(Ynew(i*10+ 17)))*(2*Vy*cos(Ynew(i*10+ 16)) - 2*Ynew(i*10+ 12) + 2*Vz*cos(Ynew(i*10+ 17))*sin(Ynew(i*10+ 16)) + 2*Vx*sin(Ynew(i*10+ 17))*sin(Ynew(i*10+ 16))))/sqrt(pow((Ynew(i*10+ 11) - Vx*cos(Ynew(i*10+ 17)) + Vz*sin(Ynew(i*10+ 17))),2) + pow((Vy*cos(Ynew(i*10+ 16)) - Ynew(i*10+ 12) + Vz*cos(Ynew(i*10+ 17))*sin(Ynew(i*10+ 16)) + Vx*sin(Ynew(i*10+ 17))*sin(Ynew(i*10+ 16))),2)) - M*deltaS*sin(Ynew(i*10+ 16))*(Yold(i*10+ 17) - Ynew(i*10+ 17));
        temp(i*10 + 7, i*10 + 12) = deltaS*(2*M - 2*ma) - deltaS*deltaT*(0.5*Cdb*rho*sqrt(Ynew(i*10+ 13)/(A*E) + 1)*sqrt(pow((Ynew(i*10+ 11) - Vx*cos(Ynew(i*10+ 17)) + Vz*sin(Ynew(i*10+ 17))),2) + pow((Vy*cos(Ynew(i*10+ 16)) - Ynew(i*10+ 12) + Vz*cos(Ynew(i*10+ 17))*sin(Ynew(i*10+ 16)) + Vx*sin(Ynew(i*10+ 17))*sin(Ynew(i*10+ 16))),2)) + (0.25*Cdb*rho*sqrt(Ynew(i*10+ 13)/(A*E) + 1)*(Vy*cos(Ynew(i*10+ 16)) - Ynew(i*10+ 12) + Vz*cos(Ynew(i*10+ 17))*sin(Ynew(i*10+ 16)) + Vx*sin(Ynew(i*10+ 17))*sin(Ynew(i*10+ 16)))*(2*Vy*cos(Ynew(i*10+ 16)) - 2*Ynew(i*10+ 12) + 2*Vz*cos(Ynew(i*10+ 17))*sin(Ynew(i*10+ 16)) + 2*Vx*sin(Ynew(i*10+ 17))*sin(Ynew(i*10+ 16))))/sqrt(pow((Ynew(i*10+ 11) - Vx*cos(Ynew(i*10+ 17)) + Vz*sin(Ynew(i*10+ 17))),2) + pow((Vy*cos(Ynew(i*10+ 16)) - Ynew(i*10+ 12) + Vz*cos(Ynew(i*10+ 17))*sin(Ynew(i*10+ 16)) + Vx*sin(Ynew(i*10+ 17))*sin(Ynew(i*10+ 16))),2)));
        temp(i*10 + 10, i*10 + 12) = -Ynew(i*10+ 18 )*deltaS*deltaT;
        temp(i*10 + 11, i*10 + 12) = -Ynew(i*10+ 19 )*deltaS*deltaT*tan(Ynew(i*10+ 16 ));
        temp(i*10 + 12, i*10 + 12) = 2*deltaT;

        temp(i*10 + 5, i*10 + 13) = - 2*deltaT - (0.25*Cdt*d0*deltaS*deltaT*pi*rho*abs(Ynew(i*10+ 10) + Vy*sin(Ynew(i*10+ 16)) - Vz*cos(Ynew(i*10+ 17))*cos(Ynew(i*10+ 16)) - Vx*cos(Ynew(i*10+ 16))*sin(Ynew(i*10+ 17)))*(Ynew(i*10+ 10) + Vy*sin(Ynew(i*10+ 16)) - Vz*cos(Ynew(i*10+ 17))*cos(Ynew(i*10+ 16)) - Vx*cos(Ynew(i*10+ 16))*sin(Ynew(i*10+ 17))))/(A*E*sqrt(Ynew(i*10+ 13)/(A*E) + 1));
        temp(i*10 + 6, i*10 + 13) = -deltaS*deltaT*(Ynew(i*10+ 19) + (0.25*Cdn*rho*sqrt(pow((Ynew(i*10+ 11) - Vx*cos(Ynew(i*10+ 17)) + Vz*sin(Ynew(i*10+ 17))),2) + pow((Vy*cos(Ynew(i*10+ 16)) - Ynew(i*10+ 12) + Vz*cos(Ynew(i*10+ 17))*sin(Ynew(i*10+ 16)) + Vx*sin(Ynew(i*10+ 17))*sin(Ynew(i*10+ 16))),2))*(Ynew(i*10+ 11) - Vx*cos(Ynew(i*10+ 17)) + Vz*sin(Ynew(i*10+ 17))))/(A*E*sqrt(Ynew(i*10+ 13)/(A*E) + 1)));
        temp(i*10 + 7, i*10 + 13) = deltaS*deltaT*(Ynew(i*10+ 18) + (0.25*Cdb*rho*sqrt(pow((Ynew(i*10+ 11) - Vx*cos(Ynew(i*10+ 17)) + Vz*sin(Ynew(i*10+ 17))),2) + pow((Vy*cos(Ynew(i*10+ 16)) - Ynew(i*10+ 12) + Vz*cos(Ynew(i*10+ 17))*sin(Ynew(i*10+ 16)) + Vx*sin(Ynew(i*10+ 17))*sin(Ynew(i*10+ 16))),2))*(Vy*cos(Ynew(i*10+ 16)) - Ynew(i*10+ 12) + Vz*cos(Ynew(i*10+ 17))*sin(Ynew(i*10+ 16)) + Vx*sin(Ynew(i*10+ 17))*sin(Ynew(i*10+ 16))))/(A*E*sqrt(Ynew(i*10+ 13)/(A*E) + 1)));
        temp(i*10 + 8, i*10 + 13) = -(3*Ynew(i*10+ 15 )*deltaS*deltaT*pow((Ynew(i*10+ 13 )/(A*E) + 1),2))/(A*E);
        temp(i*10 + 9, i*10 + 13) = (3*Ynew(i*10+ 14 )*deltaS*deltaT*pow((Ynew(i*10+ 13 )/(A*E) + 1),2))/(A*E);
        temp(i*10 + 10, i*10 + 13) = (2*deltaS)/(A*E);
        temp(i*10 + 11, i*10 + 13) = -(deltaS*cos(Ynew(i*10+ 16 ))*(Yold(i*10+ 17 ) - Ynew(i*10+ 17 )))/(A*E);
        temp(i*10 + 12, i*10 + 13) = -(deltaS*(Yold(i*10+ 16 ) - Ynew(i*10+ 16 )))/(A*E);
        
        temp(i*10 + 5, i*10 + 14) = Ynew(i*10+ 19)*deltaS*deltaT;
        temp(i*10 + 6, i*10 + 14) = -2*deltaT;
        temp(i*10 + 7, i*10 + 14) = -Ynew(i*10+ 19)*deltaS*deltaT*tan(Ynew(i*10+ 16));
        temp(i*10 + 9, i*10 + 14) = deltaS*deltaT*pow((Ynew(i*10+ 13 )/(A*E) + 1),3);
        
        temp(i*10 + 5, i*10 + 15) = -Ynew(i*10+ 18)*deltaS*deltaT;
        temp(i*10 + 6, i*10 + 15) = -Ynew(i*10+ 19)*deltaS*deltaT*tan(Ynew(i*10+ 16));
        temp(i*10 + 7, i*10 + 15) = -2*deltaT;
        temp(i*10 + 8, i*10 + 15) = -deltaS*deltaT*pow((Ynew(i*10+ 13 )/(A*E) + 1),3);
        
        temp(i*10 + 5, i*10 + 16) = deltaS*(M*Ynew(i*10+ 12) + M*Yold(i*10+ 12) - M*Ynew(i*10+ 11)*sin(Ynew(i*10+ 16))*(Yold(i*10+ 17) - Ynew(i*10+ 17))) - deltaS*deltaT*(0.5*Cdt*d0*pi*rho*abs(Ynew(i*10+ 10) + Vy*sin(Ynew(i*10+ 16)) - Vz*cos(Ynew(i*10+ 17))*cos(Ynew(i*10+ 16)) - Vx*cos(Ynew(i*10+ 16))*sin(Ynew(i*10+ 17)))*sqrt(Ynew(i*10+ 13)/(A*E) + 1)*(Vy*cos(Ynew(i*10+ 16)) + Vz*cos(Ynew(i*10+ 17))*sin(Ynew(i*10+ 16)) + Vx*sin(Ynew(i*10+ 17))*sin(Ynew(i*10+ 16))) - w0*sin(Ynew(i*10+ 17))*sin(Ynew(i*10+ 16)) + 0.5*Cdt*d0*pi*rho*sign(Ynew(i*10+ 10) + Vy*sin(Ynew(i*10+ 16)) - Vz*cos(Ynew(i*10+ 17))*cos(Ynew(i*10+ 16)) - Vx*cos(Ynew(i*10+ 16))*sin(Ynew(i*10+ 17)))*sqrt(Ynew(i*10+ 13)/(A*E) + 1)*(Vy*cos(Ynew(i*10+ 16)) + Vz*cos(Ynew(i*10+ 17))*sin(Ynew(i*10+ 16)) + Vx*sin(Ynew(i*10+ 17))*sin(Ynew(i*10+ 16)))*(Ynew(i*10+ 10) + Vy*sin(Ynew(i*10+ 16)) - Vz*cos(Ynew(i*10+ 17))*cos(Ynew(i*10+ 16)) - Vx*cos(Ynew(i*10+ 16))*sin(Ynew(i*10+ 17))));
        temp(i*10 + 6, i*10 + 16) = - deltaS*(M*Ynew(i*10+ 12)*cos(Ynew(i*10+ 16)) - M*Ynew(i*10+ 10)*sin(Ynew(i*10+ 16)))*(Yold(i*10+ 17) - Ynew(i*10+ 17)) - deltaS*deltaT*(Ynew(i*10+ 19)*Ynew(i*10+ 15)*(pow(tan(Ynew(i*10+ 16)),2) + 1) + (0.5*Cdn*rho*sqrt(Ynew(i*10+ 13)/(A*E) + 1)*(Vz*cos(Ynew(i*10+ 17))*cos(Ynew(i*10+ 16)) - Vy*sin(Ynew(i*10+ 16)) + Vx*cos(Ynew(i*10+ 16))*sin(Ynew(i*10+ 17)))*(Ynew(i*10+ 11) - Vx*cos(Ynew(i*10+ 17)) + Vz*sin(Ynew(i*10+ 17)))*(Vy*cos(Ynew(i*10+ 16)) - Ynew(i*10+ 12) + Vz*cos(Ynew(i*10+ 17))*sin(Ynew(i*10+ 16)) + Vx*sin(Ynew(i*10+ 17))*sin(Ynew(i*10+ 16))))/sqrt(pow((Ynew(i*10+ 11) - Vx*cos(Ynew(i*10+ 17)) + Vz*sin(Ynew(i*10+ 17))),2) + pow((Vy*cos(Ynew(i*10+ 16)) - Ynew(i*10+ 12) + Vz*cos(Ynew(i*10+ 17))*sin(Ynew(i*10+ 16)) + Vx*sin(Ynew(i*10+ 17))*sin(Ynew(i*10+ 16))),2)));
        temp(i*10 + 7, i*10 + 16) = deltaS*((Yold(i*10+ 17) - Ynew(i*10+ 17))*(M*Ynew(i*10+ 11)*cos(Ynew(i*10+ 16)) - Vx*ma*cos(Ynew(i*10+ 17))*cos(Ynew(i*10+ 16)) + Vz*ma*cos(Ynew(i*10+ 16))*sin(Ynew(i*10+ 17))) - M*Yold(i*10+ 10) - M*Ynew(i*10+ 10) + (Yold(i*10+ 16) - Ynew(i*10+ 16))*(Vy*ma*cos(Ynew(i*10+ 16)) + Vz*ma*cos(Ynew(i*10+ 17))*sin(Ynew(i*10+ 16)) + Vx*ma*sin(Ynew(i*10+ 17))*sin(Ynew(i*10+ 16))) - Vy*ma*sin(Yold(i*10+ 16)) - Vy*ma*sin(Ynew(i*10+ 16)) + Vz*ma*cos(Yold(i*10+ 17))*cos(Yold(i*10+ 16)) + Vz*ma*cos(Ynew(i*10+ 17))*cos(Ynew(i*10+ 16)) + Vx*ma*cos(Yold(i*10+ 16))*sin(Yold(i*10+ 17)) + Vx*ma*cos(Ynew(i*10+ 16))*sin(Ynew(i*10+ 17))) - deltaS*deltaT*(w0*cos(Ynew(i*10+ 16))*sin(Ynew(i*10+ 17)) + Ynew(i*10+ 19)*Ynew(i*10+ 14)*(pow(tan(Ynew(i*10+ 16)),2) + 1) - 0.5*Cdb*rho*sqrt(Ynew(i*10+ 13)/(A*E) + 1)*sqrt(pow((Ynew(i*10+ 11) - Vx*cos(Ynew(i*10+ 17)) + Vz*sin(Ynew(i*10+ 17))),2) + pow((Vy*cos(Ynew(i*10+ 16)) - Ynew(i*10+ 12) + Vz*cos(Ynew(i*10+ 17))*sin(Ynew(i*10+ 16)) + Vx*sin(Ynew(i*10+ 17))*sin(Ynew(i*10+ 16))),2))*(Vz*cos(Ynew(i*10+ 17))*cos(Ynew(i*10+ 16)) - Vy*sin(Ynew(i*10+ 16)) + Vx*cos(Ynew(i*10+ 16))*sin(Ynew(i*10+ 17))) - (0.5*Cdb*rho*sqrt(Ynew(i*10+ 13)/(A*E) + 1)*(Vz*cos(Ynew(i*10+ 17))*cos(Ynew(i*10+ 16)) - Vy*sin(Ynew(i*10+ 16)) + Vx*cos(Ynew(i*10+ 16))*sin(Ynew(i*10+ 17)))*pow((Vy*cos(Ynew(i*10+ 16)) - Ynew(i*10+ 12) + Vz*cos(Ynew(i*10+ 17))*sin(Ynew(i*10+ 16)) + Vx*sin(Ynew(i*10+ 17))*sin(Ynew(i*10+ 16))),2))/sqrt(pow((Ynew(i*10+ 11) - Vx*cos(Ynew(i*10+ 17)) + Vz*sin(Ynew(i*10+ 17))),2) + pow((Vy*cos(Ynew(i*10+ 16)) - Ynew(i*10+ 12) + Vz*cos(Ynew(i*10+ 17))*sin(Ynew(i*10+ 16)) + Vx*sin(Ynew(i*10+ 17))*sin(Ynew(i*10+ 16))),2)));
        temp(i*10 + 8, i*10 + 16) = E*I*pow(Ynew(i*10+ 19 ),2)*deltaS*deltaT*(pow(tan(Ynew(i*10+ 16 )),2) + 1);
        temp(i*10 + 9, i*10 + 16) = -E*I*Ynew(i*10+ 18 )*Ynew(i*10+ 19 )*deltaS*deltaT*(pow(tan(Ynew(i*10+ 16 )),2) + 1);
        temp(i*10 + 11, i*10 + 16) = deltaS*sin(Ynew(i*10+ 16 ))*(Ynew(i*10+ 13 )/(A*E) + 1)*(Yold(i*10+ 17 ) - Ynew(i*10+ 17 )) - Ynew(i*10+ 19 )*deltaS*deltaT*Ynew(i*10+ 12 )*(pow(tan(Ynew(i*10+ 16 )),2) + 1);
        temp(i*10 + 12, i*10 + 16) = deltaS*(Ynew(i*10+ 13 )/(A*E) + Yold(i*10+ 13 )/(A*E) + 2) - Ynew(i*10+ 19 )*deltaS*deltaT*Ynew(i*10+ 11 )*(pow(tan(Ynew(i*10+ 16 )),2) + 1);
        temp(i*10 + 13, i*10 + 16) = -2*deltaT;
        temp(i*10 + 14, i*10 + 16) = -deltaT*sin(Ynew(i*10+ 16 ))*(Ynew(i*10+ 7 ) - Ynew(i*10+ 17 ));
        
        temp(i*10 + 5, i*10 + 17) = deltaS*deltaT*(0.5*Cdt*d0*pi*rho*abs(Ynew(i*10+ 10) + Vy*sin(Ynew(i*10+ 16)) - Vz*cos(Ynew(i*10+ 17))*cos(Ynew(i*10+ 16)) - Vx*cos(Ynew(i*10+ 16))*sin(Ynew(i*10+ 17)))*sqrt(Ynew(i*10+ 13)/(A*E) + 1)*(Vx*cos(Ynew(i*10+ 17))*cos(Ynew(i*10+ 16)) - Vz*cos(Ynew(i*10+ 16))*sin(Ynew(i*10+ 17))) - w0*cos(Ynew(i*10+ 17))*cos(Ynew(i*10+ 16)) + 0.5*Cdt*d0*pi*rho*sign(Ynew(i*10+ 10) + Vy*sin(Ynew(i*10+ 16)) - Vz*cos(Ynew(i*10+ 17))*cos(Ynew(i*10+ 16)) - Vx*cos(Ynew(i*10+ 16))*sin(Ynew(i*10+ 17)))*sqrt(Ynew(i*10+ 13)/(A*E) + 1)*(Vx*cos(Ynew(i*10+ 17))*cos(Ynew(i*10+ 16)) - Vz*cos(Ynew(i*10+ 16))*sin(Ynew(i*10+ 17)))*(Ynew(i*10+ 10) + Vy*sin(Ynew(i*10+ 16)) - Vz*cos(Ynew(i*10+ 17))*cos(Ynew(i*10+ 16)) - Vx*cos(Ynew(i*10+ 16))*sin(Ynew(i*10+ 17)))) - deltaS*(M*Ynew(i*10+ 11)*cos(Ynew(i*10+ 16)) + M*Yold(i*10+ 11)*cos(Yold(i*10+ 16)));
        temp(i*10 + 6, i*10 + 17) = deltaS*((Vx*ma*cos(Ynew(i*10+ 17)) - Vz*ma*sin(Ynew(i*10+ 17)))*(Yold(i*10+ 17) - Ynew(i*10+ 17)) - Vz*ma*cos(Yold(i*10+ 17)) - Vz*ma*cos(Ynew(i*10+ 17)) + M*Ynew(i*10+ 10)*cos(Ynew(i*10+ 16)) + M*Yold(i*10+ 10)*cos(Yold(i*10+ 16)) - Vx*ma*sin(Yold(i*10+ 17)) - Vx*ma*sin(Ynew(i*10+ 17)) + M*Ynew(i*10+ 12)*sin(Ynew(i*10+ 16)) + M*Yold(i*10+ 12)*sin(Yold(i*10+ 16))) - deltaS*deltaT*(0.5*Cdn*rho*(Vz*cos(Ynew(i*10+ 17)) + Vx*sin(Ynew(i*10+ 17)))*sqrt(Ynew(i*10+ 13)/(A*E) + 1)*sqrt(pow((Ynew(i*10+ 11) - Vx*cos(Ynew(i*10+ 17)) + Vz*sin(Ynew(i*10+ 17))),2) + pow((Vy*cos(Ynew(i*10+ 16)) - Ynew(i*10+ 12) + Vz*cos(Ynew(i*10+ 17))*sin(Ynew(i*10+ 16)) + Vx*sin(Ynew(i*10+ 17))*sin(Ynew(i*10+ 16))),2)) - w0*sin(Ynew(i*10+ 17)) + (0.25*Cdn*rho*(2*(Vx*cos(Ynew(i*10+ 17))*sin(Ynew(i*10+ 16)) - Vz*sin(Ynew(i*10+ 17))*sin(Ynew(i*10+ 16)))*(Vy*cos(Ynew(i*10+ 16)) - Ynew(i*10+ 12) + Vz*cos(Ynew(i*10+ 17))*sin(Ynew(i*10+ 16)) + Vx*sin(Ynew(i*10+ 17))*sin(Ynew(i*10+ 16))) + 2*(Vz*cos(Ynew(i*10+ 17)) + Vx*sin(Ynew(i*10+ 17)))*(Ynew(i*10+ 11) - Vx*cos(Ynew(i*10+ 17)) + Vz*sin(Ynew(i*10+ 17))))*sqrt(Ynew(i*10+ 13)/(A*E) + 1)*(Ynew(i*10+ 11) - Vx*cos(Ynew(i*10+ 17)) + Vz*sin(Ynew(i*10+ 17))))/sqrt(pow((Ynew(i*10+ 11) - Vx*cos(Ynew(i*10+ 17)) + Vz*sin(Ynew(i*10+ 17))),2) + pow((Vy*cos(Ynew(i*10+ 16)) - Ynew(i*10+ 12) + Vz*cos(Ynew(i*10+ 17))*sin(Ynew(i*10+ 16)) + Vx*sin(Ynew(i*10+ 17))*sin(Ynew(i*10+ 16))),2)));
        temp(i*10 + 7, i*10 + 17) = deltaS*deltaT*(0.5*Cdb*rho*sqrt(Ynew(i*10+ 13)/(A*E) + 1)*(Vx*cos(Ynew(i*10+ 17))*sin(Ynew(i*10+ 16)) - Vz*sin(Ynew(i*10+ 17))*sin(Ynew(i*10+ 16)))*sqrt(pow((Ynew(i*10+ 11) - Vx*cos(Ynew(i*10+ 17)) + Vz*sin(Ynew(i*10+ 17))),2) + pow((Vy*cos(Ynew(i*10+ 16)) - Ynew(i*10+ 12) + Vz*cos(Ynew(i*10+ 17))*sin(Ynew(i*10+ 16)) + Vx*sin(Ynew(i*10+ 17))*sin(Ynew(i*10+ 16))),2)) - w0*cos(Ynew(i*10+ 17))*sin(Ynew(i*10+ 16)) + (0.25*Cdb*rho*(2*(Vx*cos(Ynew(i*10+ 17))*sin(Ynew(i*10+ 16)) - Vz*sin(Ynew(i*10+ 17))*sin(Ynew(i*10+ 16)))*(Vy*cos(Ynew(i*10+ 16)) - Ynew(i*10+ 12) + Vz*cos(Ynew(i*10+ 17))*sin(Ynew(i*10+ 16)) + Vx*sin(Ynew(i*10+ 17))*sin(Ynew(i*10+ 16))) + 2*(Vz*cos(Ynew(i*10+ 17)) + Vx*sin(Ynew(i*10+ 17)))*(Ynew(i*10+ 11) - Vx*cos(Ynew(i*10+ 17)) + Vz*sin(Ynew(i*10+ 17))))*sqrt(Ynew(i*10+ 13)/(A*E) + 1)*(Vy*cos(Ynew(i*10+ 16)) - Ynew(i*10+ 12) + Vz*cos(Ynew(i*10+ 17))*sin(Ynew(i*10+ 16)) + Vx*sin(Ynew(i*10+ 17))*sin(Ynew(i*10+ 16))))/sqrt(pow((Ynew(i*10+ 11) - Vx*cos(Ynew(i*10+ 17)) + Vz*sin(Ynew(i*10+ 17))),2) + pow((Vy*cos(Ynew(i*10+ 16)) - Ynew(i*10+ 12) + Vz*cos(Ynew(i*10+ 17))*sin(Ynew(i*10+ 16)) + Vx*sin(Ynew(i*10+ 17))*sin(Ynew(i*10+ 16))),2))) - deltaS*((Yold(i*10+ 16) - Ynew(i*10+ 16))*(Vx*ma*cos(Ynew(i*10+ 17))*cos(Ynew(i*10+ 16)) - Vz*ma*cos(Ynew(i*10+ 16))*sin(Ynew(i*10+ 17))) - (Yold(i*10+ 17) - Ynew(i*10+ 17))*(Vz*ma*cos(Ynew(i*10+ 17))*sin(Ynew(i*10+ 16)) + Vx*ma*sin(Ynew(i*10+ 17))*sin(Ynew(i*10+ 16))) + M*Ynew(i*10+ 11)*sin(Ynew(i*10+ 16)) + M*Yold(i*10+ 11)*sin(Yold(i*10+ 16)) - Vx*ma*cos(Yold(i*10+ 17))*sin(Yold(i*10+ 16)) - Vx*ma*cos(Ynew(i*10+ 17))*sin(Ynew(i*10+ 16)) + Vz*ma*sin(Yold(i*10+ 17))*sin(Yold(i*10+ 16)) + Vz*ma*sin(Ynew(i*10+ 17))*sin(Ynew(i*10+ 16)));
        temp(i*10 + 11, i*10 + 17) = deltaS*(cos(Ynew(i*10+ 16 ))*(Ynew(i*10+ 13 )/(A*E) + 1) + cos(Yold(i*10+ 16 ))*(Yold(i*10+ 13 )/(A*E) + 1));
        temp(i*10 + 14, i*10 + 17) = -deltaT*(cos(Ynew(i*10+ 6 )) + cos(Ynew(i*10+ 16 )));

        temp(i*10 + 5, i*10 + 18) = -Ynew(i*10+ 15)*deltaS*deltaT;
        temp(i*10 + 7, i*10 + 18) = Ynew(i*10+ 13)*deltaS*deltaT;
        temp(i*10 + 8, i*10 + 18) = 2*E*I*deltaT;
        temp(i*10 + 9, i*10 + 18) = -E*I*Ynew(i*10+ 19 )*deltaS*deltaT*tan(Ynew(i*10+ 16 ));
        temp(i*10 + 10, i*10 + 18) = -deltaS*deltaT*Ynew(i*10+ 12 );
        temp(i*10 + 12, i*10 + 18) = -deltaS*deltaT*Ynew(i*10+ 10 );
        temp(i*10 + 13, i*10 + 18) = deltaS*deltaT;
        
        temp(i*10 + 5, i*10 + 19) = Ynew(i*10+ 14)*deltaS*deltaT;
        temp(i*10 + 6, i*10 + 19) = -deltaS*deltaT*(Ynew(i*10+ 13) + Ynew(i*10+ 15)*tan(Ynew(i*10+ 16)));
        temp(i*10 + 7, i*10 + 19) = -Ynew(i*10+ 14)*deltaS*deltaT*tan(Ynew(i*10+ 16));
        temp(i*10 + 8, i*10 + 19) = 2*E*I*Ynew(i*10+ 19 )*deltaS*deltaT*tan(Ynew(i*10+ 16 ));
        temp(i*10 + 9, i*10 + 19) = 2*E*I*deltaT - E*I*Ynew(i*10+ 18 )*deltaS*deltaT*tan(Ynew(i*10+ 16 ));
        temp(i*10 + 10, i*10 + 19) = deltaS*deltaT*Ynew(i*10+ 11 );
        temp(i*10 + 11, i*10 + 19) = -deltaS*deltaT*(Ynew(i*10+ 10 ) + Ynew(i*10+ 12 )*tan(Ynew(i*10+ 16 )));
        temp(i*10 + 12, i*10 + 19) = -deltaS*deltaT*Ynew(i*10+ 11 )*tan(Ynew(i*10+ 16 ));
        temp(i*10 + 14, i*10 + 19) = deltaS*deltaT;
        
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

// double Jacobian::D(double a)
// {
//     temp(i*10 + 9, i*10 + 3) = deltaS*deltaT*((D(Ynew(i*10+ 4))*(Ynew(i*10+ 3)/(A*E) + 1)*(Ynew(i*10+ 3)/(A*E) + 1)^2)/(A*E) + (2*Ynew(i*10+ 4)(Ynew(i*10+ 3)/(A*E) + 1)*(Ynew(i*10+ 3)/(A*E) + 1))/(A*E));
// }