#include "../Head/Fiber.h"
// #include <spdlog/spdlog.h>

FiberMain::FiberMain()
    : times(500), Error(1e-10), Nodes(50), TotNoV(500), TimeStep(10000),
      DelTime(0.0001) {}

FiberMain::~FiberMain() {}

VectorXd FiberMain::initializeTransVal(int index) {
  // SPDLOG_DEBUG("Initializing the transient value for index {}", index);
  VectorXd TransVal(TotNoV);
  if (index > 0) {
    FiberRO aa;
    TransVal = aa.ReadTheLastRow(index - 1);
    // cout << TransVal;
  } else {
    VectorXd a = VectorXd::Zero(TotNoV);
    for (int i = 0; i < Nodes; i++) {
      a(i * 10 + 0) = 1e-11; // u
      a(i * 10 + 1) = 1e-11; // v
      a(i * 10 + 2) = 1e-11; // w
      a(i * 10 + 3) = 0;     // T
      a(i * 10 + 4) = 0;     // Sn
      a(i * 10 + 5) = 0;     // Sb
      a(i * 10 + 6) = 1e-11; // Theta
      a(i * 10 + 7) = 1e-11; // Phi
      a(i * 10 + 8) = 1e-11; // Omega
      a(i * 10 + 9) = 1e-11; // Omega
    }
    a(493) = 0;
    TransVal = a;
  }
  return TransVal;
}

MatrixXd FiberMain::initializeZeroMatrix() {
  return MatrixXd::Zero(TimeStep, TotNoV);
}

void FiberMain::Calculation(int index) {
  VectorXd TransVal = initializeTransVal(index);

  cout << endl;
  cout << "Time Step : " << index << "     Now the real time is "
       << index * DelTime << "s";

  Iterator b(TransVal, TransVal, times, Error);
  b.begin(index);
  TransVal = b.out();

  MatrixXd zero = initializeZeroMatrix();

  if (index > 0) {
    FiberRO bb;
    zero = bb.readCSV(TimeStep);
  }

  for (int j = 0; j < TotNoV; j++) {
    zero(index, j) = TransVal(j);
  }

  FiberRO a;
  a.Output(zero, TimeStep, TotNoV);
}
