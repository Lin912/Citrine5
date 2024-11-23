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

            VectorXd fx;//类存储fx
            MatrixXd jac;//类存储jac

            VectorXd Yold;//临时存储Yold
            VectorXd Ynew;//临时存储Ynew

        public:
            Iterator (VectorXd& arr, VectorXd& brr, int a, double b);
            ~Iterator ();

            void begin (int k);
            VectorXd out();

            void savetxt(Eigen::MatrixXd mat, string filename);
};