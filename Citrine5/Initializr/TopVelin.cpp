#include <iostream>
#include <Eigen/Dense>
#include <fstream>
#include <stdlib.h>
#include <string>
#include <vector>
#include <cmath>

using namespace std;
using namespace Eigen;

int main()
{
    int Nodes = 50;
    int variable = 10;
    int TV = 500;                                                 //总变量数TV= Nodes * variable
    int TimeStep = 20000;                                          //总时间步数
    double DelTime = 0.001;                                       //时间步长(真实时间步长)
    double pi = 3.1415926;

    MatrixXd a(TimeStep, TV);
    a.Zero(TimeStep, TV);
    //*********************************************************************//
    MatrixXd b(TimeStep, 3);
    b.Zero(TimeStep, 3);
    for(int i = 0; i < TimeStep; i++)
    {
        if(i <= 6000){
            b(i, 0) = 0.00001;
            b(i, 1) = 0;
            b(i, 2) = 0;
        }
        else{
            b(i, 0) = 0.2 * 2*pi/2.0 *cos(2*pi/2.0 * (i - 6000) * DelTime);
            b(i, 1) = 0;
            b(i, 2) = 0;
        }
    }

    ofstream dataFile;
    dataFile.open(".././csv./TopVel.csv", ios::out | ios::trunc);
    for(int i = 0; i < TimeStep; i++)
    {
        for(int j = 0; j < 3; j++)
        {
            dataFile << b(i, j) << ",";
        }
        dataFile << endl;
    }
    dataFile.close();
    return 0;
}