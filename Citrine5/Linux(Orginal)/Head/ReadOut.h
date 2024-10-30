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



#ifdef _WIN32
#include <windows.h>
#endif

using namespace std;
using namespace Eigen;

class FiberRO
{
    private:
        string fileTopVel;
        string fileWater;
        string fileBottomVelocity;

        string filePhysical;
        string fileObject;
        string fileDelta;

        string outfile;
        string Outfile;

        int TotNoV;

        string Topforceout;
        string Bottomforceout;


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
        

        void Output(MatrixXd &arr, int rr, int ll);
        double flutov(char* filename);

        vector<double> ReadDelta();

        VectorXd ReadTheLastRow(int index);//读取第i行数据output

        MatrixXd readCSV(int row);

        void OutTopforce(VectorXd v);
        void OutBottomforce(VectorXd v);


};