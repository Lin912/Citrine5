#include "Fiber.h"

FiberMain::FiberMain()
{
    times = 20;        //内部迭代次数
    Error = 5e-08;      //内部容许误差
    Nodes = 50;
    variable = 10;
    TotNoV = 500;
    TimeStep = 1000;
    DelTime = 0.01;

    // Mcsv(TimeStep, TotNoV);                                    //初始CSV空表格
    // Mcsv.Zero(TimeStep, TotNoV);  
    


}

FiberMain::~FiberMain()
{}

void FiberMain::Calculation(int index)
{   
    VectorXd TransVal(TotNoV);                                          //计算容器

    if(index == 0)
    {
        MatrixXd aaa(TimeStep, TotNoV);
        MatrixXd az = aaa.Zero(TimeStep, TotNoV);
        FiberRO bbb;
        bbb.Output(az, TimeStep, TotNoV);
    }

    cout << endl;
    cout << "Time Step : " << index << "     Now the real time is " << index*DelTime << "s";

    if(index > 0)
    {
        FiberRO aa;
        TransVal = aa.ReadTheLastRow(index - 1);

    }
    else
    {   
        VectorXd a(TotNoV);
        a.Zero(TotNoV);
        for(int i = 0; i < Nodes; i++)
        {
            a(i*10 + 0) = 0;    //u
            a(i*10 + 1) = 0.002;    //v
            a(i*10 + 2) = 0;    //w

            a(i*10 + 3) = 0;   //T
            a(i*10 + 4) = 0;    //Sn
            a(i*10 + 5) = 0;    //Sb

            a(i*10 + 6) = 0.00001;    //Phi
            a(i*10 + 7) = 0.00001;    //Theta

            a(i*10 + 8) = 0.00001;    //O2mega
            a(i*10 + 9) = 0.00001;    //O3mega
        }
        a(493) = 1.0;
        TransVal = a;
        // cout << TransVal;
    }

    //计算部
    Iterator b(TransVal, TransVal, times, Error);                //迭代设置
    b.begin(index);                                                  //迭代开始(i控制读入流场速度值)
    TransVal = b.out();

    //输运部
    MatrixXd Mcsv(TimeStep, TotNoV);                                    //创建读取CSV表格                                
    MatrixXd zero(TimeStep, TotNoV);
    

    if(index > 0)
    {
        FiberRO bb;
        zero = bb.readCSV(TimeStep);                               //读前index行后，其余行都为0，共TimeStep行
    }else{
        zero = Mcsv.Zero(TimeStep, TotNoV);                        //置空表格"output.csv"
    }
    
    
    for(int j = 0; j < TotNoV; j++)
    {
        zero(index, j) = TransVal(j);
    }


    //储存部
    FiberRO a;
    a.Output(zero, TimeStep, TotNoV);

    a.Outforce(TransVal);                                           //力输出项

}