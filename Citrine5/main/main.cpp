#include <iostream>
#include <time.h>
#include <conio.h>
#include "Fiber.h"


using namespace std;

void function0();
int main()
{
    int TimeStep = 100;
    int timeforCitrine5 = 0;

    while(1)
    {    
        for(int i = 0; i < TimeStep; i++)
        {
            int n = -1;
            FILE* file = fopen("../control./judge.txt", "r");
            fscanf(file, "%d", &n);
            fclose(file);

            printf("Data in Citrine\n");  

            FiberMain aa;
            aa.Calculation(timeforCitrine5);
            cout << endl;              
            cout <<"aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa";
            cout << endl;
            function0();

            timeforCitrine5++ ;
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

    FILE* file = fopen("../control./judge.txt", "w");
    fprintf(file, "0");
    fclose(file);
}