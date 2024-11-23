#include <iostream>
#include <time.h>
#include <conio.h>
#include "../Head/Fiber.h"


using namespace std;

void function0();

int main()
{
    int TimeStep = 10000;
    int timeforCitrine5 = 0;

    while(1)
    {    
        for(int i = 0; i < TimeStep; i++)
        {
            int n = 0;
            if(n == 0)
            {
            printf("Data in Citrine\n");  
            FiberMain aa;
            aa.Calculation(timeforCitrine5);
            cout << endl;              
            cout <<"aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa";
            cout << endl;
            timeforCitrine5++ ;
            }
            else
            {
                break;
            }
        }

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

    FILE* file = fopen("judge.txt", "w");
    fprintf(file, "1");
    fclose(file);
}
