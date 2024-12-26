#include <iostream>
#include <fstream>
#include <chrono>
#include <thread>
#include "../Head/Fiber.h"

using namespace std;

void logMessage(const string& message);
bool checkExitCondition();

int main()
{
    const int TimeStep = 10000;
    int timeForCitrine = 1;

    while (true) 
    {
        for (int i = 0; i < TimeStep; ++i)
        {
            logMessage("Data in Citrine");

            FiberMain fiberInstance;
            fiberInstance.Calculation(timeForCitrine);

            logMessage("Processing complete");
            timeForCitrine++;

            if (checkExitCondition()) 
            {
                logMessage("Exit condition met. Terminating program.");
                return 0;
            }
        }
    }

    return 0;
}

void logMessage(const string& message)
{
    cout << "[LOG]: " << message << endl;
}

bool checkExitCondition()
{
    ifstream file("judge.txt");
    if (file.is_open())
    {
        string flag;
        file >> flag;
        file.close();
        if (flag == "1")
        {
            return true;
        }
    }
    return false;
}
