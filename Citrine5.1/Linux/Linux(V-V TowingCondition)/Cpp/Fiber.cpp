#include "../Head/Fiber.h"

FiberMain::FiberMain()
{
    times = 200;
    Error = 1e-10;
    Nodes = 50;
    variable = 10;
    TotNoV = 500;
    TimeStep = 10000;
    DelTime = 0.001;
}

FiberMain::~FiberMain()
{}

void FiberMain::Calculation(int index)
{   
    VectorXd TransVal(TotNoV);

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
        //Initial Value
        VectorXd a(TotNoV);
        a.Zero(TotNoV);
        for(int i = 0; i < Nodes; i++)
        {
            a(i*10 + 0) = 0.00000000001;    //u
            a(i*10 + 1) = 0.00000000001;    //v
            a(i*10 + 2) = 0.00000000001;    //w

            a(i*10 + 3) = 0;   //T
            a(i*10 + 4) = 0;    //Sn
            a(i*10 + 5) = 0;    //Sb

            a(i*10 + 6) = 0.00000000001;    //Theta
            a(i*10 + 7) = 0.00000000001;    //Phi

            a(i*10 + 8) = 0.00000000001;    //O2mega
            a(i*10 + 9) = 0.00000000001;    //O3mega
        }
        a(493) = 0;
        TransVal = a;
    }

    //Computational Component
    Iterator b(TransVal, TransVal, times, Error);
    b.begin(index);                              
    TransVal = b.out();

    //Data transportation Component
    MatrixXd Mcsv(TimeStep, TotNoV);                                                
    MatrixXd zero(TimeStep, TotNoV);
    

    if(index > 0)
    {
        FiberRO bb;
        zero = bb.readCSV(TimeStep);                           
    }else{
        zero = Mcsv.Zero(TimeStep, TotNoV);               
    }
    
    
    for(int j = 0; j < TotNoV; j++)
    {
        zero(index, j) = TransVal(j);
    }

    //Data storage Component
    FiberRO a;
    a.Output(zero, TimeStep, TotNoV);
    a.OutTopforce(TransVal);
    a.OutBottomforce(TransVal);
}