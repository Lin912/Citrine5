#include "../Head/MNQ.h"


MNQ::MNQ(VectorXd& arr, VectorXd& brr, int index)
{
    Yold = arr;
    Ynew = brr;

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

}

MNQ::~MNQ()
{
    //  cout << "MNQ done !!" << endl;
}

MatrixXd MNQ::Mold()
{
    MatrixXd temp = MatrixXd::Zero(500, 500);

    for(int i = 0; i < 50; i++)
    {
        temp(10*i + 0, 10*i + 0) = M;
        temp(10*i + 0, 10*i + 6) = M*Yold(i*10 + 2);
        temp(10*i + 0, 10*i + 7) = -M*Yold(i*10 + 1)*cos(Yold(i*10 + 6));
        temp(10*i + 1, 10*i + 1) = M + ma;
        temp(10*i + 1, 10*i + 7) = M*Yold(i*10 + 0)*cos(Yold(i*10 + 6)) + M*Yold(i*10 + 2)*sin(Yold(i*10 + 6)) + ma*Vy*sin(Yold(i*10 + 7))+ ma*Vx*cos(Yold(i*10 + 7));
        temp(10*i + 2, 10*i + 2) = M + ma;
        temp(10*i + 2, 10*i + 6) = -M*Yold(i*10 + 0) -ma*Vx*cos(Yold(i*10 + 7))*cos(Yold(i*10 + 6)) - ma*Vy*sin(Yold(i*10 + 7))*cos(Yold(i*10 + 6)) + ma*Vz*sin(Yold(i*10 + 6));
        temp(10*i + 2, 10*i + 7) = -M*Yold(i*10 + 1)*sin(Yold(i*10 + 6)) + ma*Vx*sin(Yold(i*10 + 7))*sin(Yold(i*10 + 6)) - ma*Vy*cos(Yold(i*10 + 7))*sin(Yold(i*10 + 6));
        temp(10*i + 3, 10*i + 3) = 1/(A*E);
        temp(10*i + 4, 10*i + 7) = cos(Yold(i*10 + 6))*(Yold(i*10 + 3)/(A*E) + 1);
        temp(10*i + 5, 10*i + 6) = Yold(i*10 + 3)/(A*E) + 1;
    }

    return temp;
}

MatrixXd MNQ::Mnew()
{

    MatrixXd temp = MatrixXd::Zero(500, 500);

    for(int i = 0; i < 50; i++)
    {
        temp(10*i + 0, 10*i + 0) = M;
        temp(10*i + 0, 10*i + 6) = M*Ynew(i*10 + 2);
        temp(10*i + 0, 10*i + 7) = -M*Ynew(i*10 + 1)*cos(Ynew(i*10 + 6));
        temp(10*i + 1, 10*i + 1) = M + ma;
        temp(10*i + 1, 10*i + 7) = M*Ynew(i*10 + 0)*cos(Ynew(i*10 + 6)) + M*Ynew(i*10 + 2)*sin(Ynew(i*10 + 6)) + ma*Vy*sin(Ynew(i*10 + 7))+ ma*Vx*cos(Ynew(i*10 + 7));
        temp(10*i + 2, 10*i + 2) = M + ma;
        temp(10*i + 2, 10*i + 6) = -M*Ynew(i*10 + 0) -ma*Vx*cos(Ynew(i*10 + 7))*cos(Ynew(i*10 + 6)) - ma*Vy*sin(Ynew(i*10 + 7))*cos(Ynew(i*10 + 6)) + ma*Vz*sin(Ynew(i*10 + 6));
        temp(10*i + 2, 10*i + 7) = -M*Ynew(i*10 + 1)*sin(Ynew(i*10 + 6)) + ma*Vx*sin(Ynew(i*10 + 7))*sin(Ynew(i*10 + 6)) - ma*Vy*cos(Ynew(i*10 + 7))*sin(Ynew(i*10 + 6));
        temp(10*i + 3, 10*i + 3) = 1/(A*E);
        temp(10*i + 4, 10*i + 7) = cos(Ynew(i*10 + 6))*(Ynew(i*10 + 3)/(A*E) + 1);
        temp(10*i + 5, 10*i + 6) = Ynew(i*10 + 3)/(A*E) + 1;
    }

    return temp;

}

MatrixXd MNQ::Nold()
{
    MatrixXd temp = MatrixXd::Zero(500, 500);

    for(int i = 0; i < 50; i++)
    {
        temp(i*10 + 0, i*10 + 3) = -1.0;
        temp(i*10 + 1, i*10 + 4) = -1.0;
        temp(i*10 + 2, i*10 + 5) = -1.0;
        temp(i*10 + 3, i*10 + 0) = -1.0;
        temp(i*10 + 4, i*10 + 1) = -1.0;
        temp(i*10 + 5, i*10 + 2) = 1.0;
        temp(i*10 + 6, i*10 + 8) = E*I;
        temp(i*10 + 7, i*10 + 9) = E*I;
        temp(i*10 + 8, i*10 + 6) = 1.0;
        temp(i*10 + 9, i*10 + 7) = cos(Yold(i*10 + 6));
    }

    return temp;

}

MatrixXd MNQ::Nnew()
{

    MatrixXd temp = MatrixXd::Zero(500, 500);

    for(int i = 0; i < 50; i++)
    {   
        temp(i*10 + 0, i*10 + 3) = -1.0;
        temp(i*10 + 1, i*10 + 4) = -1.0;
        temp(i*10 + 2, i*10 + 5) = -1.0;
        temp(i*10 + 3, i*10 + 0) = -1.0;
        temp(i*10 + 4, i*10 + 1) = -1.0;
        temp(i*10 + 5, i*10 + 2) = 1.0;
        temp(i*10 + 6, i*10 + 8) = E*I;
        temp(i*10 + 7, i*10 + 9) = E*I;
        temp(i*10 + 8, i*10 + 6) = 1.0;
        temp(i*10 + 9, i*10 + 7) = cos(Ynew(i*10 + 6));
    }

    return temp;
}

VectorXd MNQ::Qold()
{  
    VectorXd temp(500);

    for(int i = 0; i < 50; i++)
    {
        temp(i*10 + 0) = Yold(i*10 + 9)*Yold(i*10 + 4) - Yold(i*10 + 8)*Yold(i*10 + 5) + Gz*sin(Yold(i*10 + 6)) - Gx*cos(Yold(i*10 + 7))*cos(Yold(i*10 + 6)) - Gy*cos(Yold(i*10 + 6))*sin(Yold(i*10 + 7)) + 0.5*Cdt*d0*pi*rho*abs(Yold(i*10 + 0) + Vz*sin(Yold(i*10 + 6)) - Vx*cos(Yold(i*10 + 7))*cos(Yold(i*10 + 6)) - Vy*cos(Yold(i*10 + 6))*sin(Yold(i*10 + 7)))*sqrt(Yold(i*10 + 3)/(A*E) + 1)*(Yold(i*10 + 0) + Vz*sin(Yold(i*10 + 6)) - Vx*cos(Yold(i*10 + 7))*cos(Yold(i*10 + 6)) - Vy*cos(Yold(i*10 + 6))*sin(Yold(i*10 + 7)));
        temp(i*10 + 1) = Gx*sin(Yold(i*10 + 7)) - Gy*cos(Yold(i*10 + 7)) - Yold(i*10 + 9)*Yold(i*10 + 3) - Yold(i*10 + 9)*Yold(i*10 + 5)*tan(Yold(i*10 + 6)) + 0.5*Cdn*d0*rho*sqrt(Yold(i*10 + 3)/(A*E) + 1)*sqrt(pow((Yold(i*10 + 1) - Vy*cos(Yold(i*10 + 7)) + Vx*sin(Yold(i*10 + 7))),2) + pow((Vz*cos(Yold(i*10 + 6)) - Yold(i*10 + 2) + Vx*cos(Yold(i*10 + 7))*sin(Yold(i*10 + 6)) + Vy*sin(Yold(i*10 + 7))*sin(Yold(i*10 + 6))),2))*(Yold(i*10 + 1) - Vy*cos(Yold(i*10 + 7)) + Vx*sin(Yold(i*10 + 7)));
        temp(i*10 + 2) = Yold(i*10 + 8)*Yold(i*10 + 3) - Gz*cos(Yold(i*10 + 6)) - Gx*cos(Yold(i*10 + 7))*sin(Yold(i*10 + 6)) - Gy*sin(Yold(i*10 + 7))*sin(Yold(i*10 + 6)) + Yold(i*10 + 9)*Yold(i*10 + 4)*tan(Yold(i*10 + 6)) - 0.5*Cdb*d0*rho*sqrt(Yold(i*10 + 3)/(A*E) + 1)*sqrt(pow((Yold(i*10 + 1) - Vy*cos(Yold(i*10 + 7)) + Vx*sin(Yold(i*10 + 7))),2) + pow((Vz*cos(Yold(i*10 + 6)) - Yold(i*10 + 2) + Vx*cos(Yold(i*10 + 7))*sin(Yold(i*10 + 6)) + Vy*sin(Yold(i*10 + 7))*sin(Yold(i*10 + 6))),2))*(Vz*cos(Yold(i*10 + 6)) - Yold(i*10 + 2) + Vx*cos(Yold(i*10 + 7))*sin(Yold(i*10 + 6)) + Vy*sin(Yold(i*10 + 7))*sin(Yold(i*10 + 6)));
        temp(i*10 + 3) = Yold(i*10 + 9)*Yold(i*10 + 1) - Yold(i*10 + 8)*Yold(i*10 + 2);
        temp(i*10 + 4) = - Yold(i*10 + 9)*Yold(i*10 + 0) - Yold(i*10 + 9)*Yold(i*10 + 2)*tan(Yold(i*10 + 6));
        temp(i*10 + 5) = - Yold(i*10 + 8)*Yold(i*10 + 0) - Yold(i*10 + 9)*Yold(i*10 + 1)*tan(Yold(i*10 + 6));
        temp(i*10 + 6) = E*I*pow(Yold(i*10 + 9),2)*tan(Yold(i*10 + 6)) - Yold(i*10 + 5)*pow((Yold(i*10 + 3)/(A*E) + 1),3);
        temp(i*10 + 7) = Yold(i*10 + 4)*pow((Yold(i*10 + 3)/(A*E) + 1),3) - E*I*Yold(i*10 + 8)*Yold(i*10 + 9)*tan(Yold(i*10 + 6));
        temp(i*10 + 8) = -Yold(i*10 + 8);
        temp(i*10 + 9) = -Yold(i*10 + 9);
    }

    return temp;   
}

VectorXd MNQ::Qnew()
{
    VectorXd temp(500);

   for(int i = 0; i < 50; i++)
    {
        temp(i*10 + 0) = Ynew(i*10 + 9)*Ynew(i*10 + 4) - Ynew(i*10 + 8)*Ynew(i*10 + 5) + Gz*sin(Ynew(i*10 + 6)) - Gx*cos(Ynew(i*10 + 7))*cos(Ynew(i*10 + 6)) - Gy*cos(Ynew(i*10 + 6))*sin(Ynew(i*10 + 7)) + 0.5*Cdt*d0*pi*rho*abs(Ynew(i*10 + 0) + Vz*sin(Ynew(i*10 + 6)) - Vx*cos(Ynew(i*10 + 7))*cos(Ynew(i*10 + 6)) - Vy*cos(Ynew(i*10 + 6))*sin(Ynew(i*10 + 7)))*sqrt(Ynew(i*10 + 3)/(A*E) + 1)*(Ynew(i*10 + 0) + Vz*sin(Ynew(i*10 + 6)) - Vx*cos(Ynew(i*10 + 7))*cos(Ynew(i*10 + 6)) - Vy*cos(Ynew(i*10 + 6))*sin(Ynew(i*10 + 7)));
        temp(i*10 + 1) = Gx*sin(Ynew(i*10 + 7)) - Gy*cos(Ynew(i*10 + 7)) - Ynew(i*10 + 9)*Ynew(i*10 + 3) - Ynew(i*10 + 9)*Ynew(i*10 + 5)*tan(Ynew(i*10 + 6)) + 0.5*Cdn*d0*rho*sqrt(Ynew(i*10 + 3)/(A*E) + 1)*sqrt(pow((Ynew(i*10 + 1) - Vy*cos(Ynew(i*10 + 7)) + Vx*sin(Ynew(i*10 + 7))),2) + pow((Vz*cos(Ynew(i*10 + 6)) - Ynew(i*10 + 2) + Vx*cos(Ynew(i*10 + 7))*sin(Ynew(i*10 + 6)) + Vy*sin(Ynew(i*10 + 7))*sin(Ynew(i*10 + 6))),2))*(Ynew(i*10 + 1) - Vy*cos(Ynew(i*10 + 7)) + Vx*sin(Ynew(i*10 + 7)));
        temp(i*10 + 2) = Ynew(i*10 + 8)*Ynew(i*10 + 3) - Gz*cos(Ynew(i*10 + 6)) - Gx*cos(Ynew(i*10 + 7))*sin(Ynew(i*10 + 6)) - Gy*sin(Ynew(i*10 + 7))*sin(Ynew(i*10 + 6)) + Ynew(i*10 + 9)*Ynew(i*10 + 4)*tan(Ynew(i*10 + 6)) - 0.5*Cdb*d0*rho*sqrt(Ynew(i*10 + 3)/(A*E) + 1)*sqrt(pow((Ynew(i*10 + 1) - Vy*cos(Ynew(i*10 + 7)) + Vx*sin(Ynew(i*10 + 7))),2) + pow((Vz*cos(Ynew(i*10 + 6)) - Ynew(i*10 + 2) + Vx*cos(Ynew(i*10 + 7))*sin(Ynew(i*10 + 6)) + Vy*sin(Ynew(i*10 + 7))*sin(Ynew(i*10 + 6))),2))*(Vz*cos(Ynew(i*10 + 6)) - Ynew(i*10 + 2) + Vx*cos(Ynew(i*10 + 7))*sin(Ynew(i*10 + 6)) + Vy*sin(Ynew(i*10 + 7))*sin(Ynew(i*10 + 6)));
        temp(i*10 + 3) = Ynew(i*10 + 9)*Ynew(i*10 + 1) - Ynew(i*10 + 8)*Ynew(i*10 + 2);
        temp(i*10 + 4) = - Ynew(i*10 + 9)*Ynew(i*10 + 0) - Ynew(i*10 + 9)*Ynew(i*10 + 2)*tan(Ynew(i*10 + 6));
        temp(i*10 + 5) = - Ynew(i*10 + 8)*Ynew(i*10 + 0) - Ynew(i*10 + 9)*Ynew(i*10 + 1)*tan(Ynew(i*10 + 6));
        temp(i*10 + 6) = E*I*pow(Ynew(i*10 + 9),2)*tan(Ynew(i*10 + 6)) - Ynew(i*10 + 5)*pow((Ynew(i*10 + 3)/(A*E) + 1),3);
        temp(i*10 + 7) = Ynew(i*10 + 4)*pow((Ynew(i*10 + 3)/(A*E) + 1),3) - E*I*Ynew(i*10 + 8)*Ynew(i*10 + 9)*tan(Ynew(i*10 + 6));
        temp(i*10 + 8) = -Ynew(i*10 + 8);
        temp(i*10 + 9) = -Ynew(i*10 + 9);
    }

    return temp;
}

void MNQ::savetxt(Eigen::MatrixXd mat, string filename)
{
    ofstream outfile(filename, ios::trunc);
    outfile << mat;
    outfile.close();
}