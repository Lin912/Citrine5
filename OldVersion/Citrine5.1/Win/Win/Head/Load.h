#pragma once
#include <iostream>
#include <cmath>
#include <Eigen/Dense>
#include "Fx.h"
#include "Jacobian.h"

using namespace std;
using namespace Eigen;

class Load
{
private:
    VectorXd Yold;
    VectorXd Ynew;


public:
    Load(VectorXd& arr, VectorXd& brr);
    ~Load();

    VectorXd LF(int k);
    MatrixXd LJ(int k);
  
    
};