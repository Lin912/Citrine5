#include <iostream>
#include <cmath>
#include <Eigen/Dense>
#include "../Head/Load.h"

Load::Load(VectorXd& arr, VectorXd& brr)
{
    Yold = arr;
    Ynew = brr;
}

Load::~Load()
{
    // cout << "Load done! Time to get Iteration !!" << endl;
}

VectorXd Load::LF(int k)
{
    Fx A(Yold, Ynew, k);
    return A.fx();
}

MatrixXd Load::LJ(int k)
{
    Jacobian A(Yold, Ynew, k);
    return A.jacobian();   
}