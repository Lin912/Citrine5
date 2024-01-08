#include <iostream>
#include <cmath>
#include <Eigen/Dense>
#include "Iterator.h"

using namespace std;
using namespace Eigen;

Iterator::Iterator (VectorXd& arr, VectorXd& brr, int a, double b)
{
    Yold = arr;
    Ynew = brr;

    times = a;//最大迭代次数
    Error = b;//迭代误差

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

     //********************************************//修剪
    Add bc(Yold, Ynew, k);
    Yold = bc.Addyold();
    Ynew = bc.Addynew();
    //********************************************//

    //*******************************************//装载
    Load load(Yold, Ynew);    
    fx = load.LF(k);//装载 F(x)
    jac = load.LJ(k);//装载 Jacobi(x)

    savetxt(jac, "./Data./loadjac.txt");                //Jacobian输出项，节点的v和w不能同时为0
    //*******************************************//



    VectorXd temp;//用于输出[delta Y]
    VectorXd deltaY;

    for(int i = 0; i < times; i++)
    {
        
        
        //*************************************************************//记步单元
            cout << "Iterator " << i + 1 << " times; " << endl;
        //***************************************************************//


        //*************************************************//求解单元
        double Lambda = 0.5;               //牛顿下山因子(对半取值)

        deltaY = jac.inverse() * fx * (-1) * Lambda;//主要求解单元

        temp = Ynew;
        Ynew += deltaY;//更新单元
        Yold = temp;
        
        savetxt(fx, ".././Data./fx.txt");
        savetxt(deltaY, ".././Data./deltaY.txt");

        //**************************************************//        




        //**************************************************//误差计算单元
        // Error 01: （增量比值判断）
        double a;//单步最大误差
        double amax = 0;
        for(int i = 10; i < 489; i++)//单步最大误差
        {
            if(Yold(i) != 0)
            {
                a = abs(deltaY(i)) / abs(Yold(i));
                if(a > amax)
                {
                    amax = a;
                    cout << "now i: " << i << "     and amax is " << amax <<endl;
                    cout << "new deltaY is " << deltaY(i) << "      and Yold is " << Yold(i) <<endl;
                }
                else
                {
        
                }
            }
        }
        cout << "now Max.incremental Percentage is :   " << amax << endl;//残差输出



        //Error 02: (fx 趋零判断)
        double b;
        double bmax = 0;
        for(int  i = 10; i < 489; i++)
        {
            b = abs(fx(i));

            if(bmax > b)
            {

            }
            else
            {
                bmax = b;
                // cout << "now i: " << i << "     and bmax is " << bmax <<endl;
            }

        }
        cout << "now Fx(abs) is :   " << bmax << endl << endl;//残差输出
        //******************************************************//


        //******************************************************//误差输出及重新装载单元（待修改）
        double aError = 0.001;
        double bError = 0.0000001;

        if(amax < aError || bmax < bError)
        {
            cout << "Iteration convergence !!" << endl;
            break;
        }
        else
        {

            //**************************************************//边界条件装载单元
            Add tempAdd(Yold, Ynew, k);//边界条件装载增加
            Yold = tempAdd.Addyold();
            Ynew = tempAdd.Addynew();

            savetxt(Ynew, ".././Data./Ynew.txt");
            //******************************************************//

            Load templd(Yold, Ynew);
            fx = templd.LF(k);
            jac = templd.LJ(k);   

            //**************************************************//边界条件装载单元
            BC tempbc(Yold, Ynew);//边界条件装载修剪
            Yold = tempbc.yold();
            Ynew = tempbc.ynew();
            //******************************************************//

        }
        //**********************************************************//

 

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