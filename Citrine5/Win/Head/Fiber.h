#pragma once
#include <iostream>
#include <Eigen/Dense>
#include <fstream>
#include <stdlib.h>
#include <string>
#include <vector>
#include "Iterator.h"
#include "ReadOut.h"

using namespace std;
using namespace Eigen;

class FiberMain
{
    private:
        int times;
        double Error;
        int Nodes;
        int variable;
        int TotNoV;//Total number of variables
        int TimeStep;
        double DelTime;

        Matrix<double, 500, 1> ans;

        int k;

    public:
        FiberMain();
        ~FiberMain();

        void Calculation(int index);


};