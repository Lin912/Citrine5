#pragma once
#include <iostream>
#include <vector>
#include <cmath>
#include <string>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <Eigen/Dense>
#include <stdlib.h>
#include <algorithm>
#include <stdexcept>
#include <filesystem>
#include <array>



#ifdef _WIN32
#include <windows.h>
#endif

using namespace std;
using namespace Eigen;
using Matrix3x3 = std::array<std::array<double, 3>, 3>;


class FiberRO
{
    private:
        string fileTopVel;
        string fileWater;
        //string fileBottomVelocity;
        string fileBottomVelocityRelative;
        string fileBottomomegaRelative;
        string fileBottomEulerAngle;


        string filePhysical;
        string fileObject;
        string fileDelta;

        string outfile;
        string Outfile;

        int TotNoV;

        string Topforceout;
        string Bottomforceout;


        vector<double> readLastLineData(const string& filePath);
        Matrix3x3 computeRotationMatrix(const vector<double>& Eulerangle);


    public:
        FiberRO();
        ~FiberRO();

        //Model 1
        vector<double> ReadTopVel(int k);
        vector<double> ReadBottomG();

        //Model 2
        vector<double> ReadBottomVel();

        vector<double> ReadWater(int k);  
        vector<double> ReadPhysical();
        

        void Output(const MatrixXd &arr, int rr, int ll);
        double flutov(char* filename);

        vector<double> ReadDelta();

        VectorXd ReadTheLastRow(int index);

        MatrixXd readCSV(int row);

        void OutTopforce(VectorXd v);
        void OutBottomforce(VectorXd v);


};