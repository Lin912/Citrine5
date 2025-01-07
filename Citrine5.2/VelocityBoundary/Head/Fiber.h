#pragma once

#include <iostream>
#include <Eigen/Dense>
#include <fstream>
#include <string>
#include <vector>
#include "Iterator.h"
#include "ReadOut.h"

class FiberMain
{
private:
    int times;                // Iteration times
    double Error;             // Convergence error
    int Nodes;                // Number of nodes
    int TotNoV;               // Total number of variables
    int TimeStep;             // Number of time steps
    double DelTime;           // Time step interval

    VectorXd initializeTransVal(int index);
    MatrixXd initializeZeroMatrix();             

public:
    FiberMain();
    ~FiberMain();
    void Calculation(int index);
};
