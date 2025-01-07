#pragma once
#include <iostream>
#include <fstream>
#include <cmath>
#include <Eigen/Dense>
#include "Load.h"
#include "BC.h"
#include "Add.h"


using namespace std;
using namespace Eigen;

class Iterator
{
        private:
            int times;//Iteartion Times
		    double Error;//Convergence Error

            VectorXd fx;
            MatrixXd jac;

            VectorXd Yold;
            VectorXd Ynew;

        public:
            Iterator (VectorXd& arr, VectorXd& brr, int a, double b);
            ~Iterator ();

            void begin (int k);
            VectorXd out();

            void savetxt(Eigen::MatrixXd mat, string filename);
};