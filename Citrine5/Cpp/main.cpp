#include <iostream>
#include <time.h>
#include <conio.h>
#include "../Head/Fiber.h"


using namespace std;

void function0();
int main()
{
    int TimeStep = 2000;
    int timeforCitrine5 = 0;

    while(1)
    {    
        for(int i = 0; i < TimeStep; i++)
        {
            int n = -1;
            FILE* file = fopen("../control./judge.txt", "r");
            fscanf(file, "%d", &n);
            fclose(file);

            if(n == 0)
            {
            printf("Data in Citrine\n");  

            FiberMain aa;
            aa.Calculation(timeforCitrine5);
            cout << endl;              
            cout <<"aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa";
            cout << endl;
            function0();

            timeforCitrine5++ ;
            }
            else
            {
                break;
            }
        }

        // for(int i = 0; i < TimeStep; i++)
        // {
        //     int n = -1;
        //     FILE* file = fopen("../control./judge.txt", "r");
        //     fscanf(file, "%d", &n);
        //     fclose(file);

        //     if(n == 0)
        //     {   
        //         printf("Data in Citrine\n");  

        //         FiberMain aa;
        //         aa.Calculation(timeforfiber);//0开始
        //         // aa.Calculation((i - 1) / 2);//1开始

        //         cout << endl;              
        //         cout <<"aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa";
        //         cout << endl;
                
        //         function0();

        //         timeforfiber++ ;
        //     }
        //     else if(n == 1)
        //     {
        //         function1();
        //     }
        //     else if(n ==2)
        //     {
        //         function2();
        //     }
        // }

    }
    return 0;
}

void function0()
{
    //code0

    clock_t start = clock();
    while(1)
    {
        if((clock() - start) >= CLOCKS_PER_SEC) break;
    }
    
    printf("Data in Citrine\n");

    FILE* file = fopen("../control./judge.txt", "w");
    fprintf(file, "0");
    fclose(file);
}

void function1()
{
    //code1
    clock_t start = clock();
    while (1)
    {
        if((clock() - start) >= CLOCKS_PER_SEC) break;
    }
    for(int i = 0; i < 10; i++)
    {
        printf("Data in fluent1\n");
    }

    FILE* file = fopen("../control./judge.txt", "w");
    fprintf(file, "2");
    fclose(file);
}

void function2()
{
    //code1
    clock_t start = clock();
    while (1)
    {
        if((clock() - start) >= CLOCKS_PER_SEC) break;
    }
    for(int i = 0; i < 10; i++)
    {
        printf("Data in fluent2\n");
    }  

    FILE* file = fopen("../control./judge.txt", "w");
    fprintf(file, "0");
    fclose(file);
}