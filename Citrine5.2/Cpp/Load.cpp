#include "../Head/Load.h"
#include <Eigen/Dense>
#include <cmath>
#include <iostream>
#include <spdlog/spdlog.h>


Load::Load(const VectorXd &arr, const VectorXd &brr) : Yold(arr), Ynew(brr) {}

Load::~Load() {}

VectorXd Load::LF(int k) {
  SPDLOG_DEBUG("Calculating LF for index {}", k);
  Fx A(Yold, Ynew, k);
  return A.fx();
}

MatrixXd Load::LJ(int k) {
  SPDLOG_DEBUG("Calculating LJ for index {}", k);
  Jacobian A(Yold, Ynew, k);
  return A.jacobian();
}