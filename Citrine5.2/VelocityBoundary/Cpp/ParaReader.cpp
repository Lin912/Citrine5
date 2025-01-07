#include "../Head/ParaReader.h"

ParaReader::ParaReader(){}
ParaReader::~ParaReader(){}

PhysicalData ParaReader::ReadAllPhysicalData(FiberRO& a, int index)
{
        PhysicalData data;
        
        auto physicalData = a.ReadPhysical();
        data.A = physicalData[0];
        data.rho = physicalData[1];
        data.d0 = physicalData[2];
        data.E = physicalData[3];
        data.I = physicalData[4];
        data.M = physicalData[5];
        data.ma = physicalData[6];
        data.Cdt = physicalData[7];
        data.Cdn = physicalData[8];
        data.Cdb = physicalData[9];
        data.pi = physicalData[10];
        data.g = physicalData[11];
        data.Gx = physicalData[12];
        data.Gy = physicalData[13];
        data.Gz = physicalData.back();

        auto waterData = a.ReadWater(index);
        data.Vx = waterData[0];
        data.Vy = waterData[1];
        data.Vz = waterData[2];

        auto topVelData = a.ReadTopVel(index);
        data.Vtx = topVelData[0];
        data.Vty = topVelData[1];
        data.Vtz = topVelData[2];

        auto bottomVelData = a.ReadBottomVel();
        data.Vbx = bottomVelData[0];
        data.Vby = bottomVelData[1];
        data.Vbz = bottomVelData[2];

        auto deltaData = a.ReadDelta();
        data.deltaT = deltaData[0];
        data.deltaS = deltaData[1];

        auto bottomG = a.ReadBottomG();
        data.Gbx = bottomG[0];
        data.Gby = bottomG[1];
        data.Gbz = bottomG[2];
        data.Ax = bottomG[3];
        data.Ay = bottomG[4];
        data.Az = bottomG[5];
        return data;
}
