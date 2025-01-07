#include <iostream>
#include <cmath>
#include <Eigen/Dense>
#include "../Head/Iterator.h"

using namespace std;
using namespace Eigen;

Iterator::Iterator (VectorXd& arr, VectorXd& brr, int a, double b)
{
    Yold = arr;
    Ynew = brr;

    times = a;
    Error = b;

    cout << endl;
    // cout << "Newton Iteration strated !!" << endl;
    cout << "Now the Maximum Iteration TIMES is " << times << endl;
    cout << "Now the Maximum Allowable Iteration ERROR is " << Error << endl << endl;

}

Iterator::~Iterator ()
{
    cout << "Iterator is jumped out !!" << endl;
}

void Iterator::begin (int k)
{

    Add bc(Yold, Ynew, k);
    Yold = bc.Addyold();
    Ynew = bc.Addynew();

    Load load(Yold, Ynew);    
    fx = load.LF(k);
    jac = load.LJ(k);

    VectorXd temp;
    VectorXd deltaY;

    for(int i = 0; i < times; i++)
    {
        cout << "Iterator " << i + 1 << " times; " << endl;
       
        double Lambda = 0.5;

        deltaY = jac.inverse() * fx * (-1) * Lambda;

        temp = Ynew;
        Ynew += deltaY;
        Yold = temp;
        
        savetxt(jac, "../Data/loadjac.txt");
        savetxt(fx, "../Data/fx.txt");
        savetxt(deltaY, "../Data/deltaY.txt");   

        // Error 1
        double a;
        double amax = 0;
        for(int i = 10; i < 489; i++)
        {
            if(Yold(i) != 0)
            {
                a = abs(deltaY(i)) / abs(Yold(i));
                if(a > amax){
                    amax = a;
                    cout << "now i: " << i << "     and amax is " << amax <<endl;
                    cout << "new deltaY is " << deltaY(i) << "      and Yold is " << Yold(i) <<endl;
                }
                else{}
            }
        }
        cout << "now Max.incremental Percentage is :   " << amax << endl;

        //Error 2
        double b;
        double bmax = 0;
        for(int  i = 10; i < 489; i++)
        {
            b = abs(fx(i));

            if(bmax > b){
            }
            else{
                bmax = b;
            }

        }
        cout << "now Fx(abs) is :   " << bmax << endl << endl;

        double aError = 0.0000000001;
        double bError = 0.0000000001;

        if(amax < aError || bmax < bError)
        {
            cout << "Iteration convergence !!" << endl;
            break;
        }
        else
        {
            Add tempAdd(Yold, Ynew, k);
            Yold = tempAdd.Addyold();
            Ynew = tempAdd.Addynew();

            savetxt(Ynew, "../Data/Ynew.txt");

            Load templd(Yold, Ynew);
            fx = templd.LF(k);
            jac = templd.LJ(k);   

            BC tempbc(Yold, Ynew);
            Yold = tempbc.yold();
            Ynew = tempbc.ynew();
        }
    }
}

VectorXd Iterator::out()
{
    return Ynew;
}


void Iterator::savetxt(Eigen::MatrixXd mat, string filename)
{
    ofstream outfile(filename, ios::trunc);
    outfile << mat;
    outfile.close();
}