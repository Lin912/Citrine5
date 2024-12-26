#pragma once
#include <iostream>
#include <fstream>
#include <cmath>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <vector>
#include "ReadOut.h"
#include "ParaReader.h"



using Eigen::SparseMatrix;
using Eigen::VectorXd;


class MNQ
{
    private:
        VectorXd Yold;
        VectorXd Ynew;

        double A, rho, d0, E, I, M, ma, Cdt, Cdn, Cdb, pi, g, Gx, Gy, Gz;
        double Vx, Vy, Vz, Vtx, Vty, Vtz, Vbx, Vby, Vbz;
        double deltaT, deltaS;
        double Gbx, Gby, Gbz, Ax, Ay, Az;

        // 稀疏矩阵的计算方法
        Eigen::SparseMatrix<double> CreateSparseMatrix(const Eigen::VectorXd& Y);
        Eigen::SparseMatrix<double> CreateSparseNMatrix(const Eigen::VectorXd& Y);
        Eigen::VectorXd CreateQVector(const Eigen::VectorXd& Y);



    public:
        MNQ(Eigen::VectorXd& arr, Eigen::VectorXd& brr, int index);
        ~MNQ();

        void savetxt(Eigen::MatrixXd mat, string filename);

        // 使用稀疏矩阵存储
        Eigen::SparseMatrix<double> Mold();
        Eigen::SparseMatrix<double> Mnew();
        Eigen::SparseMatrix<double> Nold();
        Eigen::SparseMatrix<double> Nnew();
        Eigen::VectorXd Qold();
        Eigen::VectorXd Qnew();
        
        void savetxt(Eigen::SparseMatrix<double> mat, string filename);
};