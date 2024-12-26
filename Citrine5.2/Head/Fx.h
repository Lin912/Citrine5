#pragma once
#include <iostream>
#include <cmath>
#include <Eigen/Dense>
#include <vector>
#include "ReadOut.h"
#include "MNQ.h"
#include "ParaReader.h"

class Fx {
private:
    VectorXd Yold;    
    VectorXd Ynew;    
    VectorXd temp;    
    
    double A, rho, d0, E, I, M, ma, Cdt, Cdn, Cdb, Gx, Gy, Gz, pi, g;
    double Vx, Vy, Vz;  
    double Vtx, Vty, Vtz; 
    double Vbx, Vby, Vbz; 
    double deltaT, deltaS; 
    double Gbx, Gby, Gbz; 
    double Ax, Ay, Az;
    
    int k;

public:
    Fx(VectorXd& arr, VectorXd& brr, int index);
    ~Fx();
    VectorXd fx();
};
