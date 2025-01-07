#pragma once
#include <iostream>
#include <vector>
#include "../Head/ReadOut.h"

struct PhysicalData {
    double A, rho, d0, E, I, M, ma, Cdt, Cdn, Cdb, pi, g;
    double Gx, Gy, Gz;
    double Vx, Vy, Vz;
    double Vtx, Vty, Vtz;
    double Vbx, Vby, Vbz;
    double deltaT, deltaS;
    double Gbx, Gby, Gbz, Ax, Ay, Az;
};

class ParaReader {
public:
    static PhysicalData ReadAllPhysicalData(FiberRO& a, int index);
    ParaReader();
    ~ParaReader();
};
