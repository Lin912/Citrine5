#include "../Head/Fiber.h"
#include <chrono>
#include <fstream>
#include <iostream>
#include <spdlog/spdlog.h>
#include <thread>

using namespace std;
void function0();
void function1();

void logMessage(const string &message);
bool checkExitCondition();

int main() {
  spdlog::set_level(spdlog::level::info);
  const int TimeStep = 10000;
  int timeForCitrine = 1;

  while (1) {
    for (int i = 0; i < TimeStep; i++) {
      int n = -1;
      FILE *file = fopen("../../Starccm/judge.txt", "r");
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

      if (n == 1) {
        SPDLOG_INFO("Data in Citrine");

        FiberMain fiberInstance;
        fiberInstance.Calculation(timeForCitrine);

        logMessage("Processing complete");
        // function0();

        timeForCitrine++;
      } else if (n == 0) {
        function1();
      }

      // if (checkExitCondition())
      // {
      //     logMessage("Exit condition met. Terminating program.");
      //     return 0;
      // }
    }
  }

  return 0;
}

void logMessage(const string &message) { cout << "[LOG]: " << message << endl; }

bool checkExitCondition() {
  ifstream file("judge.txt");
  if (file.is_open()) {
    string flag;
    file >> flag;
    file.close();
    if (flag == "1") {
      return true;
    }
  }
  return false;
}

void function0() {
  clock_t start = clock();
  while (1) {
    if ((clock() - start) >= CLOCKS_PER_SEC)
      break;
  }

  printf("Data in Citrine\n");

  FILE *file = fopen("../../Starccm/judge.txt", "w");
  fprintf(file, "0");
  fclose(file);
}

void function1() {
  clock_t start = clock();
  while (1) {
    if ((clock() - start) >= CLOCKS_PER_SEC)
      break;
  }
  for (int i = 0; i < 10; i++) {
    printf("Data in CCM\n");
  }
}
