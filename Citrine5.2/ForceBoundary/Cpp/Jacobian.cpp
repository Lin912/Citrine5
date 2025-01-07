#include "../Head/Jacobian.h"

using namespace Eigen;

Jacobian::Jacobian(VectorXd& arr, VectorXd& brr, int index)
    : Ynew(arr), Yold(brr), k(index)
{
    FiberRO a;
    PhysicalData physicalData = ParaReader::ReadAllPhysicalData(a, index);
    A = physicalData.A;
    rho = physicalData.rho;
    d0 = physicalData.d0;
    E = physicalData.E;
    I = physicalData.I;
    M = physicalData.M;
    ma = physicalData.ma;
    Cdt = physicalData.Cdt;
    Cdn = physicalData.Cdn;
    Cdb = physicalData.Cdb;
    pi = physicalData.pi;
    g = physicalData.g;
    Gx = physicalData.Gx;
    Gy = physicalData.Gy;
    Gz = physicalData.Gz;

    Vx = physicalData.Vx;
    Vy = physicalData.Vy;
    Vz = physicalData.Vz;

    Vtx = physicalData.Vtx;
    Vty = physicalData.Vty;
    Vtz = physicalData.Vtz;

    Vbx = physicalData.Vbx;
    Vby = physicalData.Vby;
    Vbz = physicalData.Vbz;

    deltaT = physicalData.deltaT;
    deltaS = physicalData.deltaS;

    Gbx = physicalData.Gbx;
    Gby = physicalData.Gby;
    Gbz = physicalData.Gbz;
    Ax = physicalData.Ax;
    Ay = physicalData.Ay;
    Az = physicalData.Az;
}

Jacobian::~Jacobian() {}

SparseMatrix<double> Jacobian::jacobian() {

    SparseMatrix<double> temp(500, 500);
    double cosY6 = cos(Ynew(6)), cosY7 = cos(Ynew(7)), cosY496 = cos(Ynew(496)), cosY497 = cos(Ynew(497));
    double sinY6 = sin(Ynew(6)), sinY7 = sin(Ynew(7)), sinY496 = sin(Ynew(496)), sinY497 = sin(Ynew(497));

    std::vector<Triplet<double>> triplets;

    triplets.push_back(Triplet<double>(0, 0, 1));
    triplets.push_back(Triplet<double>(0, 6, Vtz * cosY6 + Vtx * cosY7 * sinY6 + Vty * sinY7 * sinY6));
    triplets.push_back(Triplet<double>(0, 7, Vtx * cosY6 * sinY7 - Vty * cosY7 * cosY6));
    triplets.push_back(Triplet<double>(1, 1, 1));
    triplets.push_back(Triplet<double>(1, 7, Vtx * cosY7 + Vty * sinY7));
    triplets.push_back(Triplet<double>(2, 2, 1));
    triplets.push_back(Triplet<double>(2, 6, Vtz * sinY6 - Vtx * cosY7 * cosY6 - Vty * cosY6 * sinY7));
    triplets.push_back(Triplet<double>(2, 7, Vtx * sinY7 * sinY6 - Vty * cosY7 * sinY6));
    triplets.push_back(Triplet<double>(3, 8, 1));
    triplets.push_back(Triplet<double>(4, 9, 1));

    triplets.push_back(Triplet<double>(495, 493, 1));
    triplets.push_back(Triplet<double>(495, 496, Gbx * cosY496 + Gbx * cosY497 * sinY496 + Gby * sinY497 * sinY496));
    triplets.push_back(Triplet<double>(495, 497, Gbx * cosY496 * sinY497 - Gby * cosY497 * cosY496));
    triplets.push_back(Triplet<double>(496, 494, 1));
    triplets.push_back(Triplet<double>(496, 497, Gbx * cosY497 + Gby * sinY497));
    triplets.push_back(Triplet<double>(497, 495, 1));
    triplets.push_back(Triplet<double>(497, 496, Gbz * sinY496 - Gbx * cosY497 * cosY496 - Gby * cosY496 * sinY497));
    triplets.push_back(Triplet<double>(497, 497, Gbx * sinY497 * sinY496 - Gby * cosY497 * sinY496));
    triplets.push_back(Triplet<double>(498, 498, 1));
    triplets.push_back(Triplet<double>(499, 499, 1));

    temp.setFromTriplets(triplets.begin(), triplets.end());

    //Block01
    for (int i = 0; i < 49; i++) {    
        temp.coeffRef(i*10 + 5, i*10 + 0) = 2*M*deltaS + deltaS*deltaT*(0.5*Cdt*d0*pi*rho*abs(Ynew(i*10+ 0) + Vz*sin(Ynew(i*10+ 6)) - Vx*cos(Ynew(i*10+ 7))*cos(Ynew(i*10+ 6)) - Vy*cos(Ynew(i*10+ 6))*sin(Ynew(i*10+ 7)))*sqrt((Ynew(i*10+ 3)/(A*E) + 1)) + 0.5*Cdt*d0*pi*rho*sign(Ynew(i*10+ 0) + Vz*sin(Ynew(i*10+ 6)) - Vx*cos(Ynew(i*10+ 7))*cos(Ynew(i*10+ 6)) - Vy*cos(Ynew(i*10+ 6))*sin(Ynew(i*10+ 7)))*sqrt((Ynew(i*10+ 3)/(A*E) + 1))*(Ynew(i*10+ 0) + Vz*sin(Ynew(i*10+ 6)) - Vx*cos(Ynew(i*10+ 7))*cos(Ynew(i*10+ 6)) - Vy*cos(Ynew(i*10+ 6))*sin(Ynew(i*10+ 7))));
        temp.coeffRef(i*10 + 5, i*10 + 1) = M*deltaS*cos(Ynew(i*10+ 6))*(Yold(i*10+ 7) - Ynew(i*10+ 7));
        temp.coeffRef(i*10 + 5, i*10 + 2) = -M*deltaS*(Yold(i*10+ 6) - Ynew(i*10+ 6));
        temp.coeffRef(i*10 + 5, i*10 + 3) = 2*deltaT + (0.25*Cdt*d0*deltaS*deltaT*pi*rho*abs(Ynew(i*10+ 0) + Vz*sin(Ynew(i*10+ 6)) - Vx*cos(Ynew(i*10+ 7))*cos(Ynew(i*10+ 6)) - Vy*cos(Ynew(i*10+ 6))*sin(Ynew(i*10+ 7)))*(Ynew(i*10+ 0) + Vz*sin(Ynew(i*10+ 6)) - Vx*cos(Ynew(i*10+ 7))*cos(Ynew(i*10+ 6)) - Vy*cos(Ynew(i*10+ 6))*sin(Ynew(i*10+ 7))))/(A*E*sqrt((Ynew(i*10+ 3)/(A*E) + 1)));
        temp.coeffRef(i*10 + 5, i*10 + 4) = Ynew(i*10+ 9)*deltaS*deltaT;
        temp.coeffRef(i*10 + 5, i*10 + 5) = -Ynew(i*10+ 8)*deltaS*deltaT;
        temp.coeffRef(i*10 + 5, i*10 + 6) = deltaS*(M*Ynew(i*10+ 2) + M*Yold(i*10+ 2) - M*Ynew(i*10+ 1)*sin(Ynew(i*10+ 6))*(Yold(i*10+ 7) - Ynew(i*10+ 7))) + deltaS*deltaT*(Gz*cos(Ynew(i*10+ 6)) + Gx*cos(Ynew(i*10+ 7))*sin(Ynew(i*10+ 6)) + Gy*sin(Ynew(i*10+ 7))*sin(Ynew(i*10+ 6)) + 0.5*Cdt*d0*pi*rho*abs(Ynew(i*10+ 0) + Vz*sin(Ynew(i*10+ 6)) - Vx*cos(Ynew(i*10+ 7))*cos(Ynew(i*10+ 6)) - Vy*cos(Ynew(i*10+ 6))*sin(Ynew(i*10+ 7)))*sqrt((Ynew(i*10+ 3)/(A*E) + 1))*(Vz*cos(Ynew(i*10+ 6)) + Vx*cos(Ynew(i*10+ 7))*sin(Ynew(i*10+ 6)) + Vy*sin(Ynew(i*10+ 7))*sin(Ynew(i*10+ 6))) + 0.5*Cdt*d0*pi*rho*sign(Ynew(i*10+ 0) + Vz*sin(Ynew(i*10+ 6)) - Vx*cos(Ynew(i*10+ 7))*cos(Ynew(i*10+ 6)) - Vy*cos(Ynew(i*10+ 6))*sin(Ynew(i*10+ 7)))*sqrt((Ynew(i*10+ 3)/(A*E) + 1))*(Vz*cos(Ynew(i*10+ 6)) + Vx*cos(Ynew(i*10+ 7))*sin(Ynew(i*10+ 6)) + Vy*sin(Ynew(i*10+ 7))*sin(Ynew(i*10+ 6)))*(Ynew(i*10+ 0) + Vz*sin(Ynew(i*10+ 6)) - Vx*cos(Ynew(i*10+ 7))*cos(Ynew(i*10+ 6)) - Vy*cos(Ynew(i*10+ 6))*sin(Ynew(i*10+ 7))));
        temp.coeffRef(i*10 + 5, i*10 + 7) = - deltaS*(M*Ynew(i*10+ 1)*cos(Ynew(i*10+ 6)) + M*Yold(i*10+ 1)*cos(Yold(i*10+ 6))) - deltaS*deltaT*(Gy*cos(Ynew(i*10+ 7))*cos(Ynew(i*10+ 6)) - Gx*cos(Ynew(i*10+ 6))*sin(Ynew(i*10+ 7)) + 0.5*Cdt*d0*pi*rho*abs(Ynew(i*10+ 0) + Vz*sin(Ynew(i*10+ 6)) - Vx*cos(Ynew(i*10+ 7))*cos(Ynew(i*10+ 6)) - Vy*cos(Ynew(i*10+ 6))*sin(Ynew(i*10+ 7)))*sqrt((Ynew(i*10+ 3)/(A*E) + 1))*(Vy*cos(Ynew(i*10+ 7))*cos(Ynew(i*10+ 6)) - Vx*cos(Ynew(i*10+ 6))*sin(Ynew(i*10+ 7))) + 0.5*Cdt*d0*pi*rho*sign(Ynew(i*10+ 0) + Vz*sin(Ynew(i*10+ 6)) - Vx*cos(Ynew(i*10+ 7))*cos(Ynew(i*10+ 6)) - Vy*cos(Ynew(i*10+ 6))*sin(Ynew(i*10+ 7)))*sqrt((Ynew(i*10+ 3)/(A*E) + 1))*(Vy*cos(Ynew(i*10+ 7))*cos(Ynew(i*10+ 6)) - Vx*cos(Ynew(i*10+ 6))*sin(Ynew(i*10+ 7)))*(Ynew(i*10+ 0) + Vz*sin(Ynew(i*10+ 6)) - Vx*cos(Ynew(i*10+ 7))*cos(Ynew(i*10+ 6)) - Vy*cos(Ynew(i*10+ 6))*sin(Ynew(i*10+ 7))));
        temp.coeffRef(i*10 + 5, i*10 + 8) = -Ynew(i*10+ 5)*deltaS*deltaT;
        temp.coeffRef(i*10 + 5, i*10 + 9) = Ynew(i*10+ 4)*deltaS*deltaT;
        
        temp.coeffRef(i*10 + 6, i*10 + 0) = -M*deltaS*cos(Ynew(i*10+ 6))*(Yold(i*10+ 7) - Ynew(i*10+ 7));
        temp.coeffRef(i*10 + 6, i*10 + 1) = deltaS*(2*M + 2*ma) + deltaS*deltaT*(0.5*Cdn*d0*rho*sqrt((Ynew(i*10+ 3)/(A*E) + 1))*sqrt((pow((Ynew(i*10+ 1) - Vy*cos(Ynew(i*10+ 7)) + Vx*sin(Ynew(i*10+ 7))), 2) + pow((Vz*cos(Ynew(i*10+ 6)) - Ynew(i*10+ 2) + Vx*cos(Ynew(i*10+ 7))*sin(Ynew(i*10+ 6)) + Vy*sin(Ynew(i*10+ 7))*sin(Ynew(i*10+ 6))), 2))) + (0.25*Cdn*d0*rho*sqrt((Ynew(i*10+ 3)/(A*E) + 1))*(2*Ynew(i*10+ 1) - 2*Vy*cos(Ynew(i*10+ 7)) + 2*Vx*sin(Ynew(i*10+ 7)))*(Ynew(i*10+ 1) - Vy*cos(Ynew(i*10+ 7)) + Vx*sin(Ynew(i*10+ 7))))/sqrt((pow((Ynew(i*10+ 1) - Vy*cos(Ynew(i*10+ 7)) + Vx*sin(Ynew(i*10+ 7))), 2) + pow((Vz*cos(Ynew(i*10+ 6)) - Ynew(i*10+ 2) + Vx*cos(Ynew(i*10+ 7))*sin(Ynew(i*10+ 6)) + Vy*sin(Ynew(i*10+ 7))*sin(Ynew(i*10+ 6))), 2))));
        temp.coeffRef(i*10 + 6, i*10 + 2) = - M*deltaS*sin(Ynew(i*10+ 6))*(Yold(i*10+ 7) - Ynew(i*10+ 7)) - (0.25*Cdn*d0*deltaS*deltaT*rho*sqrt((Ynew(i*10+ 3)/(A*E) + 1))*(Ynew(i*10+ 1) - Vy*cos(Ynew(i*10+ 7)) + Vx*sin(Ynew(i*10+ 7)))*(2*Vz*cos(Ynew(i*10+ 6)) - 2*Ynew(i*10+ 2) + 2*Vx*cos(Ynew(i*10+ 7))*sin(Ynew(i*10+ 6)) + 2*Vy*sin(Ynew(i*10+ 7))*sin(Ynew(i*10+ 6))))/sqrt((pow((Ynew(i*10+ 1) - Vy*cos(Ynew(i*10+ 7)) + Vx*sin(Ynew(i*10+ 7))), 2) + pow((Vz*cos(Ynew(i*10+ 6)) - Ynew(i*10+ 2) + Vx*cos(Ynew(i*10+ 7))*sin(Ynew(i*10+ 6)) + Vy*sin(Ynew(i*10+ 7))*sin(Ynew(i*10+ 6))), 2)));
        temp.coeffRef(i*10 + 6, i*10 + 3) = -deltaS*deltaT*(Ynew(i*10+ 9) - (0.25*Cdn*d0*rho*sqrt((pow((Ynew(i*10+ 1) - Vy*cos(Ynew(i*10+ 7)) + Vx*sin(Ynew(i*10+ 7))), 2) + pow((Vz*cos(Ynew(i*10+ 6)) - Ynew(i*10+ 2) + Vx*cos(Ynew(i*10+ 7))*sin(Ynew(i*10+ 6)) + Vy*sin(Ynew(i*10+ 7))*sin(Ynew(i*10+ 6))), 2)))*(Ynew(i*10+ 1) - Vy*cos(Ynew(i*10+ 7)) + Vx*sin(Ynew(i*10+ 7))))/(A*E*sqrt((Ynew(i*10+ 3)/(A*E) + 1))));
        temp.coeffRef(i*10 + 6, i*10 + 4) = 2*deltaT;
        temp.coeffRef(i*10 + 6, i*10 + 5) = -Ynew(i*10+ 9)*deltaS*deltaT*tan(Ynew(i*10+ 6));
        temp.coeffRef(i*10 + 6, i*10 + 6) = - deltaS*(M*Ynew(i*10+ 2)*cos(Ynew(i*10+ 6)) - M*Ynew(i*10+ 0)*sin(Ynew(i*10+ 6)))*(Yold(i*10+ 7) - Ynew(i*10+ 7)) - deltaS*deltaT*(Ynew(i*10+ 9)*Ynew(i*10+ 5)*(pow(tan(Ynew(i*10+ 6)), 2) + 1) - (0.5*Cdn*d0*rho*sqrt((Ynew(i*10+ 3)/(A*E) + 1))*(Vx*cos(Ynew(i*10+ 7))*cos(Ynew(i*10+ 6)) - Vz*sin(Ynew(i*10+ 6)) + Vy*cos(Ynew(i*10+ 6))*sin(Ynew(i*10+ 7)))*(Ynew(i*10+ 1) - Vy*cos(Ynew(i*10+ 7)) + Vx*sin(Ynew(i*10+ 7)))*(Vz*cos(Ynew(i*10+ 6)) - Ynew(i*10+ 2) + Vx*cos(Ynew(i*10+ 7))*sin(Ynew(i*10+ 6)) + Vy*sin(Ynew(i*10+ 7))*sin(Ynew(i*10+ 6))))/sqrt((pow((Ynew(i*10+ 1) - Vy*cos(Ynew(i*10+ 7)) + Vx*sin(Ynew(i*10+ 7))), 2) + pow((Vz*cos(Ynew(i*10+ 6)) - Ynew(i*10+ 2) + Vx*cos(Ynew(i*10+ 7))*sin(Ynew(i*10+ 6)) + Vy*sin(Ynew(i*10+ 7))*sin(Ynew(i*10+ 6))), 2))));
        temp.coeffRef(i*10 + 6, i*10 + 7) = deltaS*(Vx*ma*cos(Yold(i*10+ 7)) - (Vy*ma*cos(Ynew(i*10+ 7)) - Vx*ma*sin(Ynew(i*10+ 7)))*(Yold(i*10+ 7) - Ynew(i*10+ 7)) + Vx*ma*cos(Ynew(i*10+ 7)) + M*Ynew(i*10+ 0)*cos(Ynew(i*10+ 6)) + M*Yold(i*10+ 0)*cos(Yold(i*10+ 6)) + Vy*ma*sin(Yold(i*10+ 7)) + Vy*ma*sin(Ynew(i*10+ 7)) + M*Ynew(i*10+ 2)*sin(Ynew(i*10+ 6)) + M*Yold(i*10+ 2)*sin(Yold(i*10+ 6))) + deltaS*deltaT*(Gx*cos(Ynew(i*10+ 7)) + Gy*sin(Ynew(i*10+ 7)) + 0.5*Cdn*d0*rho*(Vx*cos(Ynew(i*10+ 7)) + Vy*sin(Ynew(i*10+ 7)))*sqrt((Ynew(i*10+ 3)/(A*E) + 1))*sqrt((pow((Ynew(i*10+ 1) - Vy*cos(Ynew(i*10+ 7)) + Vx*sin(Ynew(i*10+ 7))), 2) + pow((Vz*cos(Ynew(i*10+ 6)) - Ynew(i*10+ 2) + Vx*cos(Ynew(i*10+ 7))*sin(Ynew(i*10+ 6)) + Vy*sin(Ynew(i*10+ 7))*sin(Ynew(i*10+ 6))), 2))) + (0.25*Cdn*d0*rho*(2*(Vy*cos(Ynew(i*10+ 7))*sin(Ynew(i*10+ 6)) - Vx*sin(Ynew(i*10+ 7))*sin(Ynew(i*10+ 6)))*(Vz*cos(Ynew(i*10+ 6)) - Ynew(i*10+ 2) + Vx*cos(Ynew(i*10+ 7))*sin(Ynew(i*10+ 6)) + Vy*sin(Ynew(i*10+ 7))*sin(Ynew(i*10+ 6))) + 2*(Vx*cos(Ynew(i*10+ 7)) + Vy*sin(Ynew(i*10+ 7)))*(Ynew(i*10+ 1) - Vy*cos(Ynew(i*10+ 7)) + Vx*sin(Ynew(i*10+ 7))))*sqrt((Ynew(i*10+ 3)/(A*E) + 1))*(Ynew(i*10+ 1) - Vy*cos(Ynew(i*10+ 7)) + Vx*sin(Ynew(i*10+ 7))))/sqrt((pow((Ynew(i*10+ 1) - Vy*cos(Ynew(i*10+ 7)) + Vx*sin(Ynew(i*10+ 7))), 2) + pow((Vz*cos(Ynew(i*10+ 6)) - Ynew(i*10+ 2) + Vx*cos(Ynew(i*10+ 7))*sin(Ynew(i*10+ 6)) + Vy*sin(Ynew(i*10+ 7))*sin(Ynew(i*10+ 6))), 2))));
        // temp.coeffRef(i*10 + 6, i*10 + 8) = 0.0;
        temp.coeffRef(i*10 + 6, i*10 + 9) = -deltaS*deltaT*(Ynew(i*10+ 3) + Ynew(i*10+ 5)*tan(Ynew(i*10+ 6)));
        
        temp.coeffRef(i*10 + 7, i*10 + 0) = M*deltaS*(Yold(i*10+ 6) - Ynew(i*10+ 6));
        temp.coeffRef(i*10 + 7, i*10 + 1) = M*deltaS*sin(Ynew(i*10+ 6))*(Yold(i*10+ 7) - Ynew(i*10+ 7)) - (0.25*Cdb*d0*deltaS*deltaT*rho*sqrt((Ynew(i*10+ 3)/(A*E) + 1))*(2*Ynew(i*10+ 1) - 2*Vy*cos(Ynew(i*10+ 7)) + 2*Vx*sin(Ynew(i*10+ 7)))*(Vz*cos(Ynew(i*10+ 6)) - Ynew(i*10+ 2) + Vx*cos(Ynew(i*10+ 7))*sin(Ynew(i*10+ 6)) + Vy*sin(Ynew(i*10+ 7))*sin(Ynew(i*10+ 6))))/sqrt((pow((Ynew(i*10+ 1) - Vy*cos(Ynew(i*10+ 7)) + Vx*sin(Ynew(i*10+ 7))), 2) + pow((Vz*cos(Ynew(i*10+ 6)) - Ynew(i*10+ 2) + Vx*cos(Ynew(i*10+ 7))*sin(Ynew(i*10+ 6)) + Vy*sin(Ynew(i*10+ 7))*sin(Ynew(i*10+ 6))), 2)));
        temp.coeffRef(i*10 + 7, i*10 + 2) = deltaS*(2*M + 2*ma) + deltaS*deltaT*(0.5*Cdb*d0*rho*sqrt((Ynew(i*10+ 3)/(A*E) + 1))*sqrt((pow((Ynew(i*10+ 1) - Vy*cos(Ynew(i*10+ 7)) + Vx*sin(Ynew(i*10+ 7))), 2) + pow((Vz*cos(Ynew(i*10+ 6)) - Ynew(i*10+ 2) + Vx*cos(Ynew(i*10+ 7))*sin(Ynew(i*10+ 6)) + Vy*sin(Ynew(i*10+ 7))*sin(Ynew(i*10+ 6))), 2))) + (0.25*Cdb*d0*rho*sqrt((Ynew(i*10+ 3)/(A*E) + 1))*(Vz*cos(Ynew(i*10+ 6)) - Ynew(i*10+ 2) + Vx*cos(Ynew(i*10+ 7))*sin(Ynew(i*10+ 6)) + Vy*sin(Ynew(i*10+ 7))*sin(Ynew(i*10+ 6)))*(2*Vz*cos(Ynew(i*10+ 6)) - 2*Ynew(i*10+ 2) + 2*Vx*cos(Ynew(i*10+ 7))*sin(Ynew(i*10+ 6)) + 2*Vy*sin(Ynew(i*10+ 7))*sin(Ynew(i*10+ 6))))/sqrt((pow((Ynew(i*10+ 1) - Vy*cos(Ynew(i*10+ 7)) + Vx*sin(Ynew(i*10+ 7))), 2) + pow((Vz*cos(Ynew(i*10+ 6)) - Ynew(i*10+ 2) + Vx*cos(Ynew(i*10+ 7))*sin(Ynew(i*10+ 6)) + Vy*sin(Ynew(i*10+ 7))*sin(Ynew(i*10+ 6))), 2))));
        temp.coeffRef(i*10 + 7, i*10 + 3) = deltaS*deltaT*(Ynew(i*10+ 8) - (0.25*Cdb*d0*rho*sqrt((pow((Ynew(i*10+ 1) - Vy*cos(Ynew(i*10+ 7)) + Vx*sin(Ynew(i*10+ 7))), 2) + pow((Vz*cos(Ynew(i*10+ 6)) - Ynew(i*10+ 2) + Vx*cos(Ynew(i*10+ 7))*sin(Ynew(i*10+ 6)) + Vy*sin(Ynew(i*10+ 7))*sin(Ynew(i*10+ 6))), 2)))*(Vz*cos(Ynew(i*10+ 6)) - Ynew(i*10+ 2) + Vx*cos(Ynew(i*10+ 7))*sin(Ynew(i*10+ 6)) + Vy*sin(Ynew(i*10+ 7))*sin(Ynew(i*10+ 6))))/(A*E*sqrt((Ynew(i*10+ 3)/(A*E) + 1))));
        temp.coeffRef(i*10 + 7, i*10 + 4) = Ynew(i*10+ 9)*deltaS*deltaT*tan(Ynew(i*10+ 6));
        temp.coeffRef(i*10 + 7, i*10 + 5) = 2*deltaT;
        temp.coeffRef(i*10 + 7, i*10 + 6) = - deltaS*(M*Ynew(i*10+ 0) + M*Yold(i*10+ 0) - (Yold(i*10+ 7) - Ynew(i*10+ 7))*(M*Ynew(i*10+ 1)*cos(Ynew(i*10+ 6)) + Vy*ma*cos(Ynew(i*10+ 7))*cos(Ynew(i*10+ 6)) - Vx*ma*cos(Ynew(i*10+ 6))*sin(Ynew(i*10+ 7))) + (Yold(i*10+ 6) - Ynew(i*10+ 6))*(Vz*ma*cos(Ynew(i*10+ 6)) + Vx*ma*cos(Ynew(i*10+ 7))*sin(Ynew(i*10+ 6)) + Vy*ma*sin(Ynew(i*10+ 7))*sin(Ynew(i*10+ 6))) - Vz*ma*sin(Yold(i*10+ 6)) - Vz*ma*sin(Ynew(i*10+ 6)) + Vx*ma*cos(Yold(i*10+ 7))*cos(Yold(i*10+ 6)) + Vx*ma*cos(Ynew(i*10+ 7))*cos(Ynew(i*10+ 6)) + Vy*ma*cos(Yold(i*10+ 6))*sin(Yold(i*10+ 7)) + Vy*ma*cos(Ynew(i*10+ 6))*sin(Ynew(i*10+ 7))) - deltaS*deltaT*(Gx*cos(Ynew(i*10+ 7))*cos(Ynew(i*10+ 6)) - Gz*sin(Ynew(i*10+ 6)) + Gy*cos(Ynew(i*10+ 6))*sin(Ynew(i*10+ 7)) - Ynew(i*10+ 9)*Ynew(i*10+ 4)*(pow(tan(Ynew(i*10+ 6)), 2) + 1) + 0.5*Cdb*d0*rho*sqrt((Ynew(i*10+ 3)/(A*E) + 1))*sqrt((pow((Ynew(i*10+ 1) - Vy*cos(Ynew(i*10+ 7)) + Vx*sin(Ynew(i*10+ 7))), 2) + pow((Vz*cos(Ynew(i*10+ 6)) - Ynew(i*10+ 2) + Vx*cos(Ynew(i*10+ 7))*sin(Ynew(i*10+ 6)) + Vy*sin(Ynew(i*10+ 7))*sin(Ynew(i*10+ 6))), 2)))*(Vx*cos(Ynew(i*10+ 7))*cos(Ynew(i*10+ 6)) - Vz*sin(Ynew(i*10+ 6)) + Vy*cos(Ynew(i*10+ 6))*sin(Ynew(i*10+ 7))) + (0.5*Cdb*d0*rho*sqrt((Ynew(i*10+ 3)/(A*E) + 1))*(Vx*cos(Ynew(i*10+ 7))*cos(Ynew(i*10+ 6)) - Vz*sin(Ynew(i*10+ 6)) + Vy*cos(Ynew(i*10+ 6))*sin(Ynew(i*10+ 7)))*pow((Vz*cos(Ynew(i*10+ 6)) - Ynew(i*10+ 2) + Vx*cos(Ynew(i*10+ 7))*sin(Ynew(i*10+ 6)) + Vy*sin(Ynew(i*10+ 7))*sin(Ynew(i*10+ 6))), 2))/sqrt((pow((Ynew(i*10+ 1) - Vy*cos(Ynew(i*10+ 7)) + Vx*sin(Ynew(i*10+ 7))), 2) + pow((Vz*cos(Ynew(i*10+ 6)) - Ynew(i*10+ 2) + Vx*cos(Ynew(i*10+ 7))*sin(Ynew(i*10+ 6)) + Vy*sin(Ynew(i*10+ 7))*sin(Ynew(i*10+ 6))), 2))));
        temp.coeffRef(i*10 + 7, i*10 + 7) = - deltaS*((Yold(i*10+ 7) - Ynew(i*10+ 7))*(Vx*ma*cos(Ynew(i*10+ 7))*sin(Ynew(i*10+ 6)) + Vy*ma*sin(Ynew(i*10+ 7))*sin(Ynew(i*10+ 6))) - (Yold(i*10+ 6) - Ynew(i*10+ 6))*(Vy*ma*cos(Ynew(i*10+ 7))*cos(Ynew(i*10+ 6)) - Vx*ma*cos(Ynew(i*10+ 6))*sin(Ynew(i*10+ 7))) + M*Ynew(i*10+ 1)*sin(Ynew(i*10+ 6)) + M*Yold(i*10+ 1)*sin(Yold(i*10+ 6)) + Vy*ma*cos(Yold(i*10+ 7))*sin(Yold(i*10+ 6)) + Vy*ma*cos(Ynew(i*10+ 7))*sin(Ynew(i*10+ 6)) - Vx*ma*sin(Yold(i*10+ 7))*sin(Yold(i*10+ 6)) - Vx*ma*sin(Ynew(i*10+ 7))*sin(Ynew(i*10+ 6))) - deltaS*deltaT*(Gy*cos(Ynew(i*10+ 7))*sin(Ynew(i*10+ 6)) - Gx*sin(Ynew(i*10+ 7))*sin(Ynew(i*10+ 6)) + 0.5*Cdb*d0*rho*sqrt((Ynew(i*10+ 3)/(A*E) + 1))*(Vy*cos(Ynew(i*10+ 7))*sin(Ynew(i*10+ 6)) - Vx*sin(Ynew(i*10+ 7))*sin(Ynew(i*10+ 6)))*sqrt((pow((Ynew(i*10+ 1) - Vy*cos(Ynew(i*10+ 7)) + Vx*sin(Ynew(i*10+ 7))), 2) + pow((Vz*cos(Ynew(i*10+ 6)) - Ynew(i*10+ 2) + Vx*cos(Ynew(i*10+ 7))*sin(Ynew(i*10+ 6)) + Vy*sin(Ynew(i*10+ 7))*sin(Ynew(i*10+ 6))), 2))) + (0.25*Cdb*d0*rho*(2*(Vy*cos(Ynew(i*10+ 7))*sin(Ynew(i*10+ 6)) - Vx*sin(Ynew(i*10+ 7))*sin(Ynew(i*10+ 6)))*(Vz*cos(Ynew(i*10+ 6)) - Ynew(i*10+ 2) + Vx*cos(Ynew(i*10+ 7))*sin(Ynew(i*10+ 6)) + Vy*sin(Ynew(i*10+ 7))*sin(Ynew(i*10+ 6))) + 2*(Vx*cos(Ynew(i*10+ 7)) + Vy*sin(Ynew(i*10+ 7)))*(Ynew(i*10+ 1) - Vy*cos(Ynew(i*10+ 7)) + Vx*sin(Ynew(i*10+ 7))))*sqrt((Ynew(i*10+ 3)/(A*E) + 1))*(Vz*cos(Ynew(i*10+ 6)) - Ynew(i*10+ 2) + Vx*cos(Ynew(i*10+ 7))*sin(Ynew(i*10+ 6)) + Vy*sin(Ynew(i*10+ 7))*sin(Ynew(i*10+ 6))))/sqrt((pow((Ynew(i*10+ 1) - Vy*cos(Ynew(i*10+ 7)) + Vx*sin(Ynew(i*10+ 7))), 2) + pow((Vz*cos(Ynew(i*10+ 6)) - Ynew(i*10+ 2) + Vx*cos(Ynew(i*10+ 7))*sin(Ynew(i*10+ 6)) + Vy*sin(Ynew(i*10+ 7))*sin(Ynew(i*10+ 6))), 2))));
        temp.coeffRef(i*10 + 7, i*10 + 8) = Ynew(i*10+ 3)*deltaS*deltaT;
        temp.coeffRef(i*10 + 7, i*10 + 9) = Ynew(i*10+ 4)*deltaS*deltaT*tan(Ynew(i*10+ 6));
        
        temp.coeffRef(i*10 + 8, i*10 + 0) = 2*deltaT;
        temp.coeffRef(i*10 + 8, i*10 + 1) = Ynew(i*10+ 9)*deltaS*deltaT;
        temp.coeffRef(i*10 + 8, i*10 + 2) = -Ynew(i*10+ 8)*deltaS*deltaT;
        temp.coeffRef(i*10 + 8, i*10 + 3) = (2*deltaS)/(A*E);
        // temp.coeffRef(i*10 + 8, i*10 + 4) = 0.0;
        // temp.coeffRef(i*10 + 8, i*10 + 5) = 0.0;
        // temp.coeffRef(i*10 + 8, i*10 + 6) = 0.0;
        // temp.coeffRef(i*10 + 8, i*10 + 7) = 0.0;
        temp.coeffRef(i*10 + 8, i*10 + 8) = -deltaS*deltaT*Ynew(i*10+ 2);
        temp.coeffRef(i*10 + 8, i*10 + 9) = deltaS*deltaT*Ynew(i*10+ 1);
        
        temp.coeffRef(i*10 + 9, i*10 + 0) = -Ynew(i*10+ 9)*deltaS*deltaT;
        temp.coeffRef(i*10 + 9, i*10 + 1) = 2*deltaT;
        temp.coeffRef(i*10 + 9, i*10 + 2) = -Ynew(i*10+ 9)*deltaS*deltaT*tan(Ynew(i*10+ 6));
        temp.coeffRef(i*10 + 9, i*10 + 3) = -(deltaS*cos(Ynew(i*10+ 6))*(Yold(i*10+ 7) - Ynew(i*10+ 7)))/(A*E);
        // temp.coeffRef(i*10 + 9, i*10 + 4) = 0.0;
        // temp.coeffRef(i*10 + 9, i*10 + 5) = 0.0;
        temp.coeffRef(i*10 + 9, i*10 + 6) = deltaS*sin(Ynew(i*10+ 6))*(Ynew(i*10+ 3)/(A*E) + 1)*(Yold(i*10+ 7) - Ynew(i*10+ 7)) - Ynew(i*10+ 9)*deltaS*deltaT*Ynew(i*10+ 2)*(pow(tan(Ynew(i*10+ 6)), 2) + 1);
        temp.coeffRef(i*10 + 9, i*10 + 7) = deltaS*(cos(Ynew(i*10+ 6))*(Ynew(i*10+ 3)/(A*E) + 1) + cos(Yold(i*10+ 6))*(Yold(i*10+ 3)/(A*E) + 1));
        // temp.coeffRef(i*10 + 9, i*10 + 8) = 0.0;
        temp.coeffRef(i*10 + 9, i*10 + 9) = -deltaS*deltaT*(Ynew(i*10+ 0) + Ynew(i*10+ 2)*tan(Ynew(i*10+ 6)));
        
        temp.coeffRef(i*10 + 10, i*10 + 0) = -Ynew(i*10+ 8)*deltaS*deltaT;
        temp.coeffRef(i*10 + 10, i*10 + 1) = -Ynew(i*10+ 9)*deltaS*deltaT*tan(Ynew(i*10+ 6));
        temp.coeffRef(i*10 + 10, i*10 + 2) = -2*deltaT;
        temp.coeffRef(i*10 + 10, i*10 + 3) = -(deltaS*(Yold(i*10+ 6) - Ynew(i*10+ 6)))/(A*E);
        // temp.coeffRef(i*10 + 10, i*10 + 4) = 0.0;
        // temp.coeffRef(i*10 + 10, i*10 + 5) = 0.0;
        temp.coeffRef(i*10 + 10, i*10 + 6) = deltaS*(Ynew(i*10+ 3)/(A*E) + Yold(i*10+ 3)/(A*E) + 2) - Ynew(i*10+ 9)*deltaS*deltaT*Ynew(i*10+ 1)*(pow(tan(Ynew(i*10+ 6)), 2) + 1);
        // temp.coeffRef(i*10 + 10, i*10 + 7) = 0.0;
        temp.coeffRef(i*10 + 10, i*10 + 8) = -deltaS*deltaT*Ynew(i*10+ 0);
        temp.coeffRef(i*10 + 10, i*10 + 9) = -deltaS*deltaT*Ynew(i*10+ 1)*tan(Ynew(i*10+ 6));
        
        // temp.coeffRef(i*10 + 11, i*10 + 0) = 0.0;
        // temp.coeffRef(i*10 + 11, i*10 + 1) = 0.0;
        // temp.coeffRef(i*10 + 11, i*10 + 2) = 0.0;
        temp.coeffRef(i*10 + 11, i*10 + 3) = -(3*Ynew(i*10+ 5)*deltaS*deltaT*pow((Ynew(i*10+ 3)/(A*E) + 1), 2))/(A*E);
        // temp.coeffRef(i*10 + 11, i*10 + 4) = 0.0;
        temp.coeffRef(i*10 + 11, i*10 + 5) = -deltaS*deltaT*pow((Ynew(i*10+ 3)/(A*E) + 1), 3);
        temp.coeffRef(i*10 + 11, i*10 + 6) = E*I*pow(Ynew(i*10+ 9), 2)*deltaS*deltaT*(pow(tan(Ynew(i*10+ 6)), 2) + 1);
        // temp.coeffRef(i*10 + 11, i*10 + 7) = 0.0;
        temp.coeffRef(i*10 + 11, i*10 + 8) = -2*E*I*deltaT;
        temp.coeffRef(i*10 + 11, i*10 + 9) = 2*E*I*Ynew(i*10+ 9)*deltaS*deltaT*tan(Ynew(i*10+ 6));
        
        // temp.coeffRef(i*10 + 12, i*10 + 0) = 0.0;
        // temp.coeffRef(i*10 + 12, i*10 + 1) = 0.0;
        // temp.coeffRef(i*10 + 12, i*10 + 2) = 0.0;
        temp.coeffRef(i*10 + 12, i*10 + 3) = (3*Ynew(i*10+ 4)*deltaS*deltaT*pow((Ynew(i*10+ 3)/(A*E) + 1), 2))/(A*E);
        temp.coeffRef(i*10 + 12, i*10 + 4) = deltaS*deltaT*pow((Ynew(i*10+ 3)/(A*E) + 1), 3);
        // temp.coeffRef(i*10 + 12, i*10 + 5) = 0.0;
        temp.coeffRef(i*10 + 12, i*10 + 6) = -E*I*Ynew(i*10+ 8)*Ynew(i*10+ 9)*deltaS*deltaT*(pow(tan(Ynew(i*10+ 6)), 2) + 1);
        // temp.coeffRef(i*10 + 12, i*10 + 7) = 0.0;
        temp.coeffRef(i*10 + 12, i*10 + 8) = -E*I*Ynew(i*10+ 9)*deltaS*deltaT*tan(Ynew(i*10+ 6));
        temp.coeffRef(i*10 + 12, i*10 + 9) = - 2*E*I*deltaT - E*I*Ynew(i*10+ 8)*deltaS*deltaT*tan(Ynew(i*10+ 6));
        
        // temp.coeffRef(i*10 + 13, i*10 + 0) = 0.0;
        // temp.coeffRef(i*10 + 13, i*10 + 1) = 0.0;
        // temp.coeffRef(i*10 + 13, i*10 + 2) = 0.0;
        // temp.coeffRef(i*10 + 13, i*10 + 3) = 0.0;
        // temp.coeffRef(i*10 + 13, i*10 + 4) = 0.0;
        // temp.coeffRef(i*10 + 13, i*10 + 5) = 0.0;
        temp.coeffRef(i*10 + 13, i*10 + 6) = -2*deltaT;
        // temp.coeffRef(i*10 + 13, i*10 + 7) = 0.0;
        temp.coeffRef(i*10 + 13, i*10 + 8) = -deltaS*deltaT;
        // temp.coeffRef(i*10 + 13, i*10 + 9) = 0.0;
        
        // temp.coeffRef(i*10 + 14, i*10 + 0) = 0.0;
        // temp.coeffRef(i*10 + 14, i*10 + 1) = 0.0;
        // temp.coeffRef(i*10 + 14, i*10 + 2) = 0.0;
        // temp.coeffRef(i*10 + 14, i*10 + 3) = 0.0;
        // temp.coeffRef(i*10 + 14, i*10 + 4) = 0.0;
        // temp.coeffRef(i*10 + 14, i*10 + 5) = 0.0;
        temp.coeffRef(i*10 + 14, i*10 + 6) = deltaT*sin(Ynew(i*10+ 6))*(Ynew(i*10+ 7) - Ynew(i*10+ 17));
        temp.coeffRef(i*10 + 14, i*10 + 7) = -deltaT*(cos(Ynew(i*10+ 6)) + cos(Ynew(i*10+ 16)));
        // temp.coeffRef(i*10 + 14, i*10 + 8) = 0.0;
        temp.coeffRef(i*10 + 14, i*10 + 9) = -deltaS*deltaT;
    }

    //Block02
    for(int i = 0; i < 49; i++)
    {
        temp.coeffRef(i*10 + 5, i*10 + 10) = 2*M*deltaS + deltaS*deltaT*(0.5*Cdt*d0*pi*rho*abs(Ynew(i*10+ 10) + Vz*sin(Ynew(i*10+ 16)) - Vx*cos(Ynew(i*10+ 17))*cos(Ynew(i*10+ 16)) - Vy*cos(Ynew(i*10+ 16))*sin(Ynew(i*10+ 17)))*sqrt((Ynew(i*10+ 13)/(A*E) + 1)) + 0.5*Cdt*d0*pi*rho*sign(Ynew(i*10+ 10) + Vz*sin(Ynew(i*10+ 16)) - Vx*cos(Ynew(i*10+ 17))*cos(Ynew(i*10+ 16)) - Vy*cos(Ynew(i*10+ 16))*sin(Ynew(i*10+ 17)))*sqrt((Ynew(i*10+ 13)/(A*E) + 1))*(Ynew(i*10+ 10) + Vz*sin(Ynew(i*10+ 16)) - Vx*cos(Ynew(i*10+ 17))*cos(Ynew(i*10+ 16)) - Vy*cos(Ynew(i*10+ 16))*sin(Ynew(i*10+ 17))));
        temp.coeffRef(i*10 + 5, i*10 + 11) = M*deltaS*cos(Ynew(i*10+ 16))*(Yold(i*10+ 17) - Ynew(i*10+ 17));
        temp.coeffRef(i*10 + 5, i*10 + 12) = -M*deltaS*(Yold(i*10+ 16) - Ynew(i*10+ 16));
        temp.coeffRef(i*10 + 5, i*10 + 13) = (0.25*Cdt*d0*deltaS*deltaT*pi*rho*abs(Ynew(i*10+ 10) + Vz*sin(Ynew(i*10+ 16)) - Vx*cos(Ynew(i*10+ 17))*cos(Ynew(i*10+ 16)) - Vy*cos(Ynew(i*10+ 16))*sin(Ynew(i*10+ 17)))*(Ynew(i*10+ 10) + Vz*sin(Ynew(i*10+ 16)) - Vx*cos(Ynew(i*10+ 17))*cos(Ynew(i*10+ 16)) - Vy*cos(Ynew(i*10+ 16))*sin(Ynew(i*10+ 17))))/(A*E*sqrt((Ynew(i*10+ 13)/(A*E) + 1))) - 2*deltaT;
        temp.coeffRef(i*10 + 5, i*10 + 14) = Ynew(i*10+ 19)*deltaS*deltaT;
        temp.coeffRef(i*10 + 5, i*10 + 15) = -Ynew(i*10+ 18)*deltaS*deltaT;
        temp.coeffRef(i*10 + 5, i*10 + 16) = deltaS*(M*Ynew(i*10+ 12) + M*Yold(i*10+ 12) - M*Ynew(i*10+ 11)*sin(Ynew(i*10+ 16))*(Yold(i*10+ 17) - Ynew(i*10+ 17))) + deltaS*deltaT*(Gz*cos(Ynew(i*10+ 16)) + Gx*cos(Ynew(i*10+ 17))*sin(Ynew(i*10+ 16)) + Gy*sin(Ynew(i*10+ 17))*sin(Ynew(i*10+ 16)) + 0.5*Cdt*d0*pi*rho*abs(Ynew(i*10+ 10) + Vz*sin(Ynew(i*10+ 16)) - Vx*cos(Ynew(i*10+ 17))*cos(Ynew(i*10+ 16)) - Vy*cos(Ynew(i*10+ 16))*sin(Ynew(i*10+ 17)))*sqrt((Ynew(i*10+ 13)/(A*E) + 1))*(Vz*cos(Ynew(i*10+ 16)) + Vx*cos(Ynew(i*10+ 17))*sin(Ynew(i*10+ 16)) + Vy*sin(Ynew(i*10+ 17))*sin(Ynew(i*10+ 16))) + 0.5*Cdt*d0*pi*rho*sign(Ynew(i*10+ 10) + Vz*sin(Ynew(i*10+ 16)) - Vx*cos(Ynew(i*10+ 17))*cos(Ynew(i*10+ 16)) - Vy*cos(Ynew(i*10+ 16))*sin(Ynew(i*10+ 17)))*sqrt((Ynew(i*10+ 13)/(A*E) + 1))*(Vz*cos(Ynew(i*10+ 16)) + Vx*cos(Ynew(i*10+ 17))*sin(Ynew(i*10+ 16)) + Vy*sin(Ynew(i*10+ 17))*sin(Ynew(i*10+ 16)))*(Ynew(i*10+ 10) + Vz*sin(Ynew(i*10+ 16)) - Vx*cos(Ynew(i*10+ 17))*cos(Ynew(i*10+ 16)) - Vy*cos(Ynew(i*10+ 16))*sin(Ynew(i*10+ 17))));
        temp.coeffRef(i*10 + 5, i*10 + 17) = - deltaS*(M*Ynew(i*10+ 11)*cos(Ynew(i*10+ 16)) + M*Yold(i*10+ 11)*cos(Yold(i*10+ 16))) - deltaS*deltaT*(Gy*cos(Ynew(i*10+ 17))*cos(Ynew(i*10+ 16)) - Gx*cos(Ynew(i*10+ 16))*sin(Ynew(i*10+ 17)) + 0.5*Cdt*d0*pi*rho*abs(Ynew(i*10+ 10) + Vz*sin(Ynew(i*10+ 16)) - Vx*cos(Ynew(i*10+ 17))*cos(Ynew(i*10+ 16)) - Vy*cos(Ynew(i*10+ 16))*sin(Ynew(i*10+ 17)))*sqrt((Ynew(i*10+ 13)/(A*E) + 1))*(Vy*cos(Ynew(i*10+ 17))*cos(Ynew(i*10+ 16)) - Vx*cos(Ynew(i*10+ 16))*sin(Ynew(i*10+ 17))) + 0.5*Cdt*d0*pi*rho*sign(Ynew(i*10+ 10) + Vz*sin(Ynew(i*10+ 16)) - Vx*cos(Ynew(i*10+ 17))*cos(Ynew(i*10+ 16)) - Vy*cos(Ynew(i*10+ 16))*sin(Ynew(i*10+ 17)))*sqrt((Ynew(i*10+ 13)/(A*E) + 1))*(Vy*cos(Ynew(i*10+ 17))*cos(Ynew(i*10+ 16)) - Vx*cos(Ynew(i*10+ 16))*sin(Ynew(i*10+ 17)))*(Ynew(i*10+ 10) + Vz*sin(Ynew(i*10+ 16)) - Vx*cos(Ynew(i*10+ 17))*cos(Ynew(i*10+ 16)) - Vy*cos(Ynew(i*10+ 16))*sin(Ynew(i*10+ 17))));
        temp.coeffRef(i*10 + 5, i*10 + 18) = -Ynew(i*10+ 15)*deltaS*deltaT;
        temp.coeffRef(i*10 + 5, i*10 + 19) = Ynew(i*10+ 14)*deltaS*deltaT;
        
        temp.coeffRef(i*10 + 6, i*10 + 10) = -M*deltaS*cos(Ynew(i*10+ 16))*(Yold(i*10+ 17) - Ynew(i*10+ 17));
        temp.coeffRef(i*10 + 6, i*10 + 11) = deltaS*(2*M + 2*ma) + deltaS*deltaT*(0.5*Cdn*d0*rho*sqrt((Ynew(i*10+ 13)/(A*E) + 1))*sqrt((pow((Ynew(i*10+ 11) - Vy*cos(Ynew(i*10+ 17)) + Vx*sin(Ynew(i*10+ 17))), 2) + pow((Vz*cos(Ynew(i*10+ 16)) - Ynew(i*10+ 12) + Vx*cos(Ynew(i*10+ 17))*sin(Ynew(i*10+ 16)) + Vy*sin(Ynew(i*10+ 17))*sin(Ynew(i*10+ 16))), 2))) + (0.25*Cdn*d0*rho*sqrt((Ynew(i*10+ 13)/(A*E) + 1))*(2*Ynew(i*10+ 11) - 2*Vy*cos(Ynew(i*10+ 17)) + 2*Vx*sin(Ynew(i*10+ 17)))*(Ynew(i*10+ 11) - Vy*cos(Ynew(i*10+ 17)) + Vx*sin(Ynew(i*10+ 17))))/sqrt((pow((Ynew(i*10+ 11) - Vy*cos(Ynew(i*10+ 17)) + Vx*sin(Ynew(i*10+ 17))), 2) + pow((Vz*cos(Ynew(i*10+ 16)) - Ynew(i*10+ 12) + Vx*cos(Ynew(i*10+ 17))*sin(Ynew(i*10+ 16)) + Vy*sin(Ynew(i*10+ 17))*sin(Ynew(i*10+ 16))), 2))));
        temp.coeffRef(i*10 + 6, i*10 + 12) = - M*deltaS*sin(Ynew(i*10+ 16))*(Yold(i*10+ 17) - Ynew(i*10+ 17)) - (0.25*Cdn*d0*deltaS*deltaT*rho*sqrt((Ynew(i*10+ 13)/(A*E) + 1))*(Ynew(i*10+ 11) - Vy*cos(Ynew(i*10+ 17)) + Vx*sin(Ynew(i*10+ 17)))*(2*Vz*cos(Ynew(i*10+ 16)) - 2*Ynew(i*10+ 12) + 2*Vx*cos(Ynew(i*10+ 17))*sin(Ynew(i*10+ 16)) + 2*Vy*sin(Ynew(i*10+ 17))*sin(Ynew(i*10+ 16))))/sqrt((pow((Ynew(i*10+ 11) - Vy*cos(Ynew(i*10+ 17)) + Vx*sin(Ynew(i*10+ 17))), 2) + pow((Vz*cos(Ynew(i*10+ 16)) - Ynew(i*10+ 12) + Vx*cos(Ynew(i*10+ 17))*sin(Ynew(i*10+ 16)) + Vy*sin(Ynew(i*10+ 17))*sin(Ynew(i*10+ 16))), 2)));
        temp.coeffRef(i*10 + 6, i*10 + 13) = -deltaS*deltaT*(Ynew(i*10+ 19) - (0.25*Cdn*d0*rho*sqrt((pow((Ynew(i*10+ 11) - Vy*cos(Ynew(i*10+ 17)) + Vx*sin(Ynew(i*10+ 17))), 2) + pow((Vz*cos(Ynew(i*10+ 16)) - Ynew(i*10+ 12) + Vx*cos(Ynew(i*10+ 17))*sin(Ynew(i*10+ 16)) + Vy*sin(Ynew(i*10+ 17))*sin(Ynew(i*10+ 16))), 2)))*(Ynew(i*10+ 11) - Vy*cos(Ynew(i*10+ 17)) + Vx*sin(Ynew(i*10+ 17))))/(A*E*sqrt((Ynew(i*10+ 13)/(A*E) + 1))));
        temp.coeffRef(i*10 + 6, i*10 + 14) = -2*deltaT;
        temp.coeffRef(i*10 + 6, i*10 + 15) = -Ynew(i*10+ 19)*deltaS*deltaT*tan(Ynew(i*10+ 16));
        temp.coeffRef(i*10 + 6, i*10 + 16) = - deltaS*(M*Ynew(i*10+ 12)*cos(Ynew(i*10+ 16)) - M*Ynew(i*10+ 10)*sin(Ynew(i*10+ 16)))*(Yold(i*10+ 17) - Ynew(i*10+ 17)) - deltaS*deltaT*(Ynew(i*10+ 19)*Ynew(i*10+ 15)*(pow(tan(Ynew(i*10+ 16)), 2) + 1) - (0.5*Cdn*d0*rho*sqrt((Ynew(i*10+ 13)/(A*E) + 1))*(Vx*cos(Ynew(i*10+ 17))*cos(Ynew(i*10+ 16)) - Vz*sin(Ynew(i*10+ 16)) + Vy*cos(Ynew(i*10+ 16))*sin(Ynew(i*10+ 17)))*(Ynew(i*10+ 11) - Vy*cos(Ynew(i*10+ 17)) + Vx*sin(Ynew(i*10+ 17)))*(Vz*cos(Ynew(i*10+ 16)) - Ynew(i*10+ 12) + Vx*cos(Ynew(i*10+ 17))*sin(Ynew(i*10+ 16)) + Vy*sin(Ynew(i*10+ 17))*sin(Ynew(i*10+ 16))))/sqrt((pow((Ynew(i*10+ 11) - Vy*cos(Ynew(i*10+ 17)) + Vx*sin(Ynew(i*10+ 17))), 2) + pow((Vz*cos(Ynew(i*10+ 16)) - Ynew(i*10+ 12) + Vx*cos(Ynew(i*10+ 17))*sin(Ynew(i*10+ 16)) + Vy*sin(Ynew(i*10+ 17))*sin(Ynew(i*10+ 16))), 2))));
        temp.coeffRef(i*10 + 6, i*10 + 17) = deltaS*(Vx*ma*cos(Yold(i*10+ 17)) - (Vy*ma*cos(Ynew(i*10+ 17)) - Vx*ma*sin(Ynew(i*10+ 17)))*(Yold(i*10+ 17) - Ynew(i*10+ 17)) + Vx*ma*cos(Ynew(i*10+ 17)) + M*Ynew(i*10+ 10)*cos(Ynew(i*10+ 16)) + M*Yold(i*10+ 10)*cos(Yold(i*10+ 16)) + Vy*ma*sin(Yold(i*10+ 17)) + Vy*ma*sin(Ynew(i*10+ 17)) + M*Ynew(i*10+ 12)*sin(Ynew(i*10+ 16)) + M*Yold(i*10+ 12)*sin(Yold(i*10+ 16))) + deltaS*deltaT*(Gx*cos(Ynew(i*10+ 17)) + Gy*sin(Ynew(i*10+ 17)) + 0.5*Cdn*d0*rho*(Vx*cos(Ynew(i*10+ 17)) + Vy*sin(Ynew(i*10+ 17)))*sqrt((Ynew(i*10+ 13)/(A*E) + 1))*sqrt((pow((Ynew(i*10+ 11) - Vy*cos(Ynew(i*10+ 17)) + Vx*sin(Ynew(i*10+ 17))), 2) + pow((Vz*cos(Ynew(i*10+ 16)) - Ynew(i*10+ 12) + Vx*cos(Ynew(i*10+ 17))*sin(Ynew(i*10+ 16)) + Vy*sin(Ynew(i*10+ 17))*sin(Ynew(i*10+ 16))), 2))) + (0.25*Cdn*d0*rho*(2*(Vy*cos(Ynew(i*10+ 17))*sin(Ynew(i*10+ 16)) - Vx*sin(Ynew(i*10+ 17))*sin(Ynew(i*10+ 16)))*(Vz*cos(Ynew(i*10+ 16)) - Ynew(i*10+ 12) + Vx*cos(Ynew(i*10+ 17))*sin(Ynew(i*10+ 16)) + Vy*sin(Ynew(i*10+ 17))*sin(Ynew(i*10+ 16))) + 2*(Vx*cos(Ynew(i*10+ 17)) + Vy*sin(Ynew(i*10+ 17)))*(Ynew(i*10+ 11) - Vy*cos(Ynew(i*10+ 17)) + Vx*sin(Ynew(i*10+ 17))))*sqrt((Ynew(i*10+ 13)/(A*E) + 1))*(Ynew(i*10+ 11) - Vy*cos(Ynew(i*10+ 17)) + Vx*sin(Ynew(i*10+ 17))))/sqrt((pow((Ynew(i*10+ 11) - Vy*cos(Ynew(i*10+ 17)) + Vx*sin(Ynew(i*10+ 17))), 2) + pow((Vz*cos(Ynew(i*10+ 16)) - Ynew(i*10+ 12) + Vx*cos(Ynew(i*10+ 17))*sin(Ynew(i*10+ 16)) + Vy*sin(Ynew(i*10+ 17))*sin(Ynew(i*10+ 16))), 2))));
        // temp.coeffRef(i*10 + 6, i*10 + 18) = 0.0;
        temp.coeffRef(i*10 + 6, i*10 + 19) = -deltaS*deltaT*(Ynew(i*10+ 13) + Ynew(i*10+ 15)*tan(Ynew(i*10+ 16)));
        
        temp.coeffRef(i*10 + 7, i*10 + 10) = M*deltaS*(Yold(i*10+ 16) - Ynew(i*10+ 16));
        temp.coeffRef(i*10 + 7, i*10 + 11) = M*deltaS*sin(Ynew(i*10+ 16))*(Yold(i*10+ 17) - Ynew(i*10+ 17)) - (0.25*Cdb*d0*deltaS*deltaT*rho*sqrt((Ynew(i*10+ 13)/(A*E) + 1))*(2*Ynew(i*10+ 11) - 2*Vy*cos(Ynew(i*10+ 17)) + 2*Vx*sin(Ynew(i*10+ 17)))*(Vz*cos(Ynew(i*10+ 16)) - Ynew(i*10+ 12) + Vx*cos(Ynew(i*10+ 17))*sin(Ynew(i*10+ 16)) + Vy*sin(Ynew(i*10+ 17))*sin(Ynew(i*10+ 16))))/sqrt((pow((Ynew(i*10+ 11) - Vy*cos(Ynew(i*10+ 17)) + Vx*sin(Ynew(i*10+ 17))), 2) + pow((Vz*cos(Ynew(i*10+ 16)) - Ynew(i*10+ 12) + Vx*cos(Ynew(i*10+ 17))*sin(Ynew(i*10+ 16)) + Vy*sin(Ynew(i*10+ 17))*sin(Ynew(i*10+ 16))), 2)));
        temp.coeffRef(i*10 + 7, i*10 + 12) = deltaS*(2*M + 2*ma) + deltaS*deltaT*(0.5*Cdb*d0*rho*sqrt((Ynew(i*10+ 13)/(A*E) + 1))*sqrt((pow((Ynew(i*10+ 11) - Vy*cos(Ynew(i*10+ 17)) + Vx*sin(Ynew(i*10+ 17))), 2) + pow((Vz*cos(Ynew(i*10+ 16)) - Ynew(i*10+ 12) + Vx*cos(Ynew(i*10+ 17))*sin(Ynew(i*10+ 16)) + Vy*sin(Ynew(i*10+ 17))*sin(Ynew(i*10+ 16))), 2))) + (0.25*Cdb*d0*rho*sqrt((Ynew(i*10+ 13)/(A*E) + 1))*(Vz*cos(Ynew(i*10+ 16)) - Ynew(i*10+ 12) + Vx*cos(Ynew(i*10+ 17))*sin(Ynew(i*10+ 16)) + Vy*sin(Ynew(i*10+ 17))*sin(Ynew(i*10+ 16)))*(2*Vz*cos(Ynew(i*10+ 16)) - 2*Ynew(i*10+ 12) + 2*Vx*cos(Ynew(i*10+ 17))*sin(Ynew(i*10+ 16)) + 2*Vy*sin(Ynew(i*10+ 17))*sin(Ynew(i*10+ 16))))/sqrt((pow((Ynew(i*10+ 11) - Vy*cos(Ynew(i*10+ 17)) + Vx*sin(Ynew(i*10+ 17))), 2) + pow((Vz*cos(Ynew(i*10+ 16)) - Ynew(i*10+ 12) + Vx*cos(Ynew(i*10+ 17))*sin(Ynew(i*10+ 16)) + Vy*sin(Ynew(i*10+ 17))*sin(Ynew(i*10+ 16))), 2))));
        temp.coeffRef(i*10 + 7, i*10 + 13) = deltaS*deltaT*(Ynew(i*10+ 18) - (0.25*Cdb*d0*rho*sqrt((pow((Ynew(i*10+ 11) - Vy*cos(Ynew(i*10+ 17)) + Vx*sin(Ynew(i*10+ 17))), 2) + pow((Vz*cos(Ynew(i*10+ 16)) - Ynew(i*10+ 12) + Vx*cos(Ynew(i*10+ 17))*sin(Ynew(i*10+ 16)) + Vy*sin(Ynew(i*10+ 17))*sin(Ynew(i*10+ 16))), 2)))*(Vz*cos(Ynew(i*10+ 16)) - Ynew(i*10+ 12) + Vx*cos(Ynew(i*10+ 17))*sin(Ynew(i*10+ 16)) + Vy*sin(Ynew(i*10+ 17))*sin(Ynew(i*10+ 16))))/(A*E*sqrt((Ynew(i*10+ 13)/(A*E) + 1))));
        temp.coeffRef(i*10 + 7, i*10 + 14) = Ynew(i*10+ 19)*deltaS*deltaT*tan(Ynew(i*10+ 16));
        temp.coeffRef(i*10 + 7, i*10 + 15) = -2*deltaT;
        temp.coeffRef(i*10 + 7, i*10 + 16) = - deltaS*(M*Ynew(i*10+ 10) + M*Yold(i*10+ 10) - (Yold(i*10+ 17) - Ynew(i*10+ 17))*(M*Ynew(i*10+ 11)*cos(Ynew(i*10+ 16)) + Vy*ma*cos(Ynew(i*10+ 17))*cos(Ynew(i*10+ 16)) - Vx*ma*cos(Ynew(i*10+ 16))*sin(Ynew(i*10+ 17))) + (Yold(i*10+ 16) - Ynew(i*10+ 16))*(Vz*ma*cos(Ynew(i*10+ 16)) + Vx*ma*cos(Ynew(i*10+ 17))*sin(Ynew(i*10+ 16)) + Vy*ma*sin(Ynew(i*10+ 17))*sin(Ynew(i*10+ 16))) - Vz*ma*sin(Yold(i*10+ 16)) - Vz*ma*sin(Ynew(i*10+ 16)) + Vx*ma*cos(Yold(i*10+ 17))*cos(Yold(i*10+ 16)) + Vx*ma*cos(Ynew(i*10+ 17))*cos(Ynew(i*10+ 16)) + Vy*ma*cos(Yold(i*10+ 16))*sin(Yold(i*10+ 17)) + Vy*ma*cos(Ynew(i*10+ 16))*sin(Ynew(i*10+ 17))) - deltaS*deltaT*(Gx*cos(Ynew(i*10+ 17))*cos(Ynew(i*10+ 16)) - Gz*sin(Ynew(i*10+ 16)) + Gy*cos(Ynew(i*10+ 16))*sin(Ynew(i*10+ 17)) - Ynew(i*10+ 19)*Ynew(i*10+ 14)*(pow(tan(Ynew(i*10+ 16)), 2) + 1) + 0.5*Cdb*d0*rho*sqrt((Ynew(i*10+ 13)/(A*E) + 1))*sqrt((pow((Ynew(i*10+ 11) - Vy*cos(Ynew(i*10+ 17)) + Vx*sin(Ynew(i*10+ 17))), 2) + pow((Vz*cos(Ynew(i*10+ 16)) - Ynew(i*10+ 12) + Vx*cos(Ynew(i*10+ 17))*sin(Ynew(i*10+ 16)) + Vy*sin(Ynew(i*10+ 17))*sin(Ynew(i*10+ 16))), 2)))*(Vx*cos(Ynew(i*10+ 17))*cos(Ynew(i*10+ 16)) - Vz*sin(Ynew(i*10+ 16)) + Vy*cos(Ynew(i*10+ 16))*sin(Ynew(i*10+ 17))) + (0.5*Cdb*d0*rho*sqrt((Ynew(i*10+ 13)/(A*E) + 1))*(Vx*cos(Ynew(i*10+ 17))*cos(Ynew(i*10+ 16)) - Vz*sin(Ynew(i*10+ 16)) + Vy*cos(Ynew(i*10+ 16))*sin(Ynew(i*10+ 17)))*pow((Vz*cos(Ynew(i*10+ 16)) - Ynew(i*10+ 12) + Vx*cos(Ynew(i*10+ 17))*sin(Ynew(i*10+ 16)) + Vy*sin(Ynew(i*10+ 17))*sin(Ynew(i*10+ 16))), 2))/sqrt((pow((Ynew(i*10+ 11) - Vy*cos(Ynew(i*10+ 17)) + Vx*sin(Ynew(i*10+ 17))), 2) + pow((Vz*cos(Ynew(i*10+ 16)) - Ynew(i*10+ 12) + Vx*cos(Ynew(i*10+ 17))*sin(Ynew(i*10+ 16)) + Vy*sin(Ynew(i*10+ 17))*sin(Ynew(i*10+ 16))), 2))));
        temp.coeffRef(i*10 + 7, i*10 + 17) = - deltaS*((Yold(i*10+ 17) - Ynew(i*10+ 17))*(Vx*ma*cos(Ynew(i*10+ 17))*sin(Ynew(i*10+ 16)) + Vy*ma*sin(Ynew(i*10+ 17))*sin(Ynew(i*10+ 16))) - (Yold(i*10+ 16) - Ynew(i*10+ 16))*(Vy*ma*cos(Ynew(i*10+ 17))*cos(Ynew(i*10+ 16)) - Vx*ma*cos(Ynew(i*10+ 16))*sin(Ynew(i*10+ 17))) + M*Ynew(i*10+ 11)*sin(Ynew(i*10+ 16)) + M*Yold(i*10+ 11)*sin(Yold(i*10+ 16)) + Vy*ma*cos(Yold(i*10+ 17))*sin(Yold(i*10+ 16)) + Vy*ma*cos(Ynew(i*10+ 17))*sin(Ynew(i*10+ 16)) - Vx*ma*sin(Yold(i*10+ 17))*sin(Yold(i*10+ 16)) - Vx*ma*sin(Ynew(i*10+ 17))*sin(Ynew(i*10+ 16))) - deltaS*deltaT*(Gy*cos(Ynew(i*10+ 17))*sin(Ynew(i*10+ 16)) - Gx*sin(Ynew(i*10+ 17))*sin(Ynew(i*10+ 16)) + 0.5*Cdb*d0*rho*sqrt((Ynew(i*10+ 13)/(A*E) + 1))*(Vy*cos(Ynew(i*10+ 17))*sin(Ynew(i*10+ 16)) - Vx*sin(Ynew(i*10+ 17))*sin(Ynew(i*10+ 16)))*sqrt((pow((Ynew(i*10+ 11) - Vy*cos(Ynew(i*10+ 17)) + Vx*sin(Ynew(i*10+ 17))), 2) + pow((Vz*cos(Ynew(i*10+ 16)) - Ynew(i*10+ 12) + Vx*cos(Ynew(i*10+ 17))*sin(Ynew(i*10+ 16)) + Vy*sin(Ynew(i*10+ 17))*sin(Ynew(i*10+ 16))), 2))) + (0.25*Cdb*d0*rho*(2*(Vy*cos(Ynew(i*10+ 17))*sin(Ynew(i*10+ 16)) - Vx*sin(Ynew(i*10+ 17))*sin(Ynew(i*10+ 16)))*(Vz*cos(Ynew(i*10+ 16)) - Ynew(i*10+ 12) + Vx*cos(Ynew(i*10+ 17))*sin(Ynew(i*10+ 16)) + Vy*sin(Ynew(i*10+ 17))*sin(Ynew(i*10+ 16))) + 2*(Vx*cos(Ynew(i*10+ 17)) + Vy*sin(Ynew(i*10+ 17)))*(Ynew(i*10+ 11) - Vy*cos(Ynew(i*10+ 17)) + Vx*sin(Ynew(i*10+ 17))))*sqrt((Ynew(i*10+ 13)/(A*E) + 1))*(Vz*cos(Ynew(i*10+ 16)) - Ynew(i*10+ 12) + Vx*cos(Ynew(i*10+ 17))*sin(Ynew(i*10+ 16)) + Vy*sin(Ynew(i*10+ 17))*sin(Ynew(i*10+ 16))))/sqrt((pow((Ynew(i*10+ 11) - Vy*cos(Ynew(i*10+ 17)) + Vx*sin(Ynew(i*10+ 17))), 2) + pow((Vz*cos(Ynew(i*10+ 16)) - Ynew(i*10+ 12) + Vx*cos(Ynew(i*10+ 17))*sin(Ynew(i*10+ 16)) + Vy*sin(Ynew(i*10+ 17))*sin(Ynew(i*10+ 16))), 2))));
        temp.coeffRef(i*10 + 7, i*10 + 18) = Ynew(i*10+ 13)*deltaS*deltaT;
        temp.coeffRef(i*10 + 7, i*10 + 19) = Ynew(i*10+ 14)*deltaS*deltaT*tan(Ynew(i*10+ 16));
        
        temp.coeffRef(i*10 + 8, i*10 + 10) = -2*deltaT;
        temp.coeffRef(i*10 + 8, i*10 + 11) = Ynew(i*10+ 19)*deltaS*deltaT;
        temp.coeffRef(i*10 + 8, i*10 + 12) = -Ynew(i*10+ 18)*deltaS*deltaT;
        temp.coeffRef(i*10 + 8, i*10 + 13) = (2*deltaS)/(A*E);
        // temp.coeffRef(i*10 + 8, i*10 + 14) = 0.0;
        // temp.coeffRef(i*10 + 8, i*10 + 15) = 0.0;
        // temp.coeffRef(i*10 + 8, i*10 + 16) = 0.0;
        // temp.coeffRef(i*10 + 8, i*10 + 17) = 0.0;
        temp.coeffRef(i*10 + 8, i*10 + 18) = -deltaS*deltaT*Ynew(i*10+ 12);
        temp.coeffRef(i*10 + 8, i*10 + 19) = deltaS*deltaT*Ynew(i*10+ 11);
        
        temp.coeffRef(i*10 + 9, i*10 + 10) = -Ynew(i*10+ 19)*deltaS*deltaT;
        temp.coeffRef(i*10 + 9, i*10 + 11) = -2*deltaT;
        temp.coeffRef(i*10 + 9, i*10 + 12) = -Ynew(i*10+ 19)*deltaS*deltaT*tan(Ynew(i*10+ 16));
        temp.coeffRef(i*10 + 9, i*10 + 13) = -(deltaS*cos(Ynew(i*10+ 16))*(Yold(i*10+ 17) - Ynew(i*10+ 17)))/(A*E);
        // temp.coeffRef(i*10 + 9, i*10 + 14) = 0.0;
        // temp.coeffRef(i*10 + 9, i*10 + 15) = 0.0;
        temp.coeffRef(i*10 + 9, i*10 + 16) = deltaS*sin(Ynew(i*10+ 16))*(Ynew(i*10+ 13)/(A*E) + 1)*(Yold(i*10+ 17) - Ynew(i*10+ 17)) - Ynew(i*10+ 19)*deltaS*deltaT*Ynew(i*10+ 12)*(pow(tan(Ynew(i*10+ 16)), 2) + 1);
        temp.coeffRef(i*10 + 9, i*10 + 17) = deltaS*(cos(Ynew(i*10+ 16))*(Ynew(i*10+ 13)/(A*E) + 1) + cos(Yold(i*10+ 16))*(Yold(i*10+ 13)/(A*E) + 1));
        // temp.coeffRef(i*10 + 9, i*10 + 18) = 0.0;
        temp.coeffRef(i*10 + 9, i*10 + 19) = -deltaS*deltaT*(Ynew(i*10+ 10) + Ynew(i*10+ 12)*tan(Ynew(i*10+ 16)));
        
        temp.coeffRef(i*10 + 10, i*10 + 10) = -Ynew(i*10+ 18)*deltaS*deltaT;
        temp.coeffRef(i*10 + 10, i*10 + 11) = -Ynew(i*10+ 19)*deltaS*deltaT*tan(Ynew(i*10+ 16));
        temp.coeffRef(i*10 + 10, i*10 + 12) = 2*deltaT;
        temp.coeffRef(i*10 + 10, i*10 + 13) = -(deltaS*(Yold(i*10+ 16) - Ynew(i*10+ 16)))/(A*E);
        // temp.coeffRef(i*10 + 10, i*10 + 14) = 0.0;
        // temp.coeffRef(i*10 + 10, i*10 + 15) = 0.0;
        temp.coeffRef(i*10 + 10, i*10 + 16) = deltaS*(Ynew(i*10+ 13)/(A*E) + Yold(i*10+ 13)/(A*E) + 2) - Ynew(i*10+ 19)*deltaS*deltaT*Ynew(i*10+ 11)*(pow(tan(Ynew(i*10+ 16)), 2) + 1);
        // temp.coeffRef(i*10 + 10, i*10 + 17) = 0.0;
        temp.coeffRef(i*10 + 10, i*10 + 18) = -deltaS*deltaT*Ynew(i*10+ 10);
        temp.coeffRef(i*10 + 10, i*10 + 19) = -deltaS*deltaT*Ynew(i*10+ 11)*tan(Ynew(i*10+ 16));
        
        // temp.coeffRef(i*10 + 11, i*10 + 10) = 0.0;
        // temp.coeffRef(i*10 + 11, i*10 + 11) = 0.0;
        // temp.coeffRef(i*10 + 11, i*10 + 12) = 0.0;
        temp.coeffRef(i*10 + 11, i*10 + 13) = -(3*Ynew(i*10+ 15)*deltaS*deltaT*pow((Ynew(i*10+ 13)/(A*E) + 1), 2))/(A*E);
        // temp.coeffRef(i*10 + 11, i*10 + 14) = 0.0;
        temp.coeffRef(i*10 + 11, i*10 + 15) = -deltaS*deltaT*pow((Ynew(i*10+ 13)/(A*E) + 1), 3);
        temp.coeffRef(i*10 + 11, i*10 + 16) = E*I*pow(Ynew(i*10+ 19), 2)*deltaS*deltaT*(pow(tan(Ynew(i*10+ 16)), 2) + 1);
        // temp.coeffRef(i*10 + 11, i*10 + 17) = 0.0;
        temp.coeffRef(i*10 + 11, i*10 + 18) = 2*E*I*deltaT;
        temp.coeffRef(i*10 + 11, i*10 + 19) = 2*E*I*Ynew(i*10+ 19)*deltaS*deltaT*tan(Ynew(i*10+ 16));
        
        // temp.coeffRef(i*10 + 12, i*10 + 10) = 0.0;
        // temp.coeffRef(i*10 + 12, i*10 + 11) = 0.0;
        // temp.coeffRef(i*10 + 12, i*10 + 12) = 0.0;
        temp.coeffRef(i*10 + 12, i*10 + 13) = (3*Ynew(i*10+ 14)*deltaS*deltaT*pow((Ynew(i*10+ 13)/(A*E) + 1), 2))/(A*E);
        temp.coeffRef(i*10 + 12, i*10 + 14) = deltaS*deltaT*pow((Ynew(i*10+ 13)/(A*E) + 1), 3);
        // temp.coeffRef(i*10 + 12, i*10 + 15) = 0.0;
        temp.coeffRef(i*10 + 12, i*10 + 16) = -E*I*Ynew(i*10+ 18)*Ynew(i*10+ 19)*deltaS*deltaT*(pow(tan(Ynew(i*10+ 16)), 2) + 1);
        // temp.coeffRef(i*10 + 12, i*10 + 17) = 0.0;
        temp.coeffRef(i*10 + 12, i*10 + 18) = -E*I*Ynew(i*10+ 19)*deltaS*deltaT*tan(Ynew(i*10+ 16));
        temp.coeffRef(i*10 + 12, i*10 + 19) = 2*E*I*deltaT - E*I*Ynew(i*10+ 18)*deltaS*deltaT*tan(Ynew(i*10+ 16));
        
        // temp.coeffRef(i*10 + 13, i*10 + 10) = 0.0;
        // temp.coeffRef(i*10 + 13, i*10 + 11) = 0.0;
        // temp.coeffRef(i*10 + 13, i*10 + 12) = 0.0;
        // temp.coeffRef(i*10 + 13, i*10 + 13) = 0.0;
        // temp.coeffRef(i*10 + 13, i*10 + 14) = 0.0;
        // temp.coeffRef(i*10 + 13, i*10 + 15) = 0.0;
        temp.coeffRef(i*10 + 13, i*10 + 16) = 2*deltaT;
        // temp.coeffRef(i*10 + 13, i*10 + 17) = 0.0;
        temp.coeffRef(i*10 + 13, i*10 + 18) = -deltaS*deltaT;
        // temp.coeffRef(i*10 + 13, i*10 + 19) = 0.0;
        
        // temp.coeffRef(i*10 + 14, i*10 + 10) = 0.0;
        // temp.coeffRef(i*10 + 14, i*10 + 11) = 0.0;
        // temp.coeffRef(i*10 + 14, i*10 + 12) = 0.0;
        // temp.coeffRef(i*10 + 14, i*10 + 13) = 0.0;
        // temp.coeffRef(i*10 + 14, i*10 + 14) = 0.0;
        // temp.coeffRef(i*10 + 14, i*10 + 15) = 0.0;
        temp.coeffRef(i*10 + 14, i*10 + 16) = deltaT*sin(Ynew(i*10+ 16))*(Ynew(i*10+ 7) - Ynew(i*10+ 17));
        temp.coeffRef(i*10 + 14, i*10 + 17) = deltaT*(cos(Ynew(i*10+ 6)) + cos(Ynew(i*10+ 16)));
        // temp.coeffRef(i*10 + 14, i*10 + 18) = 0.0;
        temp.coeffRef(i*10 + 14, i*10 + 19) = -deltaS*deltaT;
    }
    //Orgin JacobianMatrix 
    return temp;  
}


int Jacobian::sign(double a)
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