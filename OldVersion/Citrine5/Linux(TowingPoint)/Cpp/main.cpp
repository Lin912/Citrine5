#include <iostream>
#include <time.h>
#include "../Head/Fiber.h"


using namespace std;

void function0();
void function1();
void function2();

int main()
{
    int TimeStep = 10000;
    int timeforCitrine5 = 1;

    while(1)
    {    
        for(int i = 0; i < TimeStep; i++)
        {
            int n = -1;
            FILE* file = fopen("../../Star/judge.txt", "r");
            if (file == nullptr) {
                std::cerr << "Failed to open file judge.txt" << std::endl;
                return 1;
            }
            if (fscanf(file, "%d", &n) != 1) {
                std::cerr << "Failed to read integer from file" << std::endl;
                fclose(file);
                return 1;
            }
            fclose(file);

            if(n == 1)
            {   
                printf("Data in Citrine\n");  
                FiberMain aa;
                aa.Calculation(timeforCitrine5);//0开始

                cout << endl;
                cout<< "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA" << endl;

                function0();

                timeforCitrine5++ ;
            }
            else if(n == 0)
            {
                function1();
            }
            else if(n ==2)
            {
                function2();
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

    FILE* file = fopen("../../Star/judge.txt", "w");
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
        printf("Data in CCM\n");
    }

    // FILE* file = fopen("../../flu./judge.txt", "w");
    // fprintf(file, "0");
    // fclose(file);
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

    FILE* file = fopen("../../flu/judge.txt", "w");
    fprintf(file, "0");
    fclose(file);
}
