#include "../Head/ReadOut.h"
#include <algorithm>
#include <spdlog/spdlog.h>

FiberRO::FiberRO()
    : fileTopVel("../csv/TopVel.csv"), fileObject("../csv/TowedObject.csv"),
      fileBottomVelocityRelative("../../Starccm/VelocityRelative.csv"),
      fileBottomomegaRelative("../../Starccm/omegaRelative.csv"),
      fileBottomEulerAngle("../../Starccm/EulerAngle.csv"),
      fileWater("../csv/Water.csv"), filePhysical("../csv/Parameters.csv"),
      fileDelta("../csv/Delta.csv"), Outfile("../csv/output.csv") {
  std::filesystem::create_directory("../csv");
  TotNoV = 500;
  // The .txt used to read into computation through udf
  Topforceout = "../../Starccm/topforce.txt";
  Bottomforceout = "../../Starccm/bottomforce.txt";
}

FiberRO::~FiberRO() {}

vector<double> FiberRO::ReadTopVel(int k) {
  ifstream infile(fileTopVel, ios::in);
  if (!infile) {
    cout << "Can not Open" + fileTopVel + ".csv" << endl;
    exit(1);
  }
  vector<double> arrdata;
  string line;
  string field;
  for (int i = 0; i <= k; i++) {
    getline(infile, line);
  }
  stringstream sin(line);
  // int colindex = 0;
  while (getline(sin, field, ',')) {
    double dvalue = atof(field.c_str());
    arrdata.push_back(dvalue);
  }
  infile.close();
  return arrdata;
}

vector<double> FiberRO::ReadWater(int k) {
  ifstream infile(fileWater, ios::in);
  if (!infile) {
    cout << "Can not Open" + fileWater + ".csv" << endl;
    exit(1);
  }
  vector<double> arrdata;
  string line;
  string field;
  for (int i = 0; i <= k; i++) {
    getline(infile, line);
  }
  stringstream sin(line);
  // int colindex = 0;
  while (getline(sin, field, ',')) {
    double dvalue = atof(field.c_str());
    arrdata.push_back(dvalue);
  }
  infile.close();
  return arrdata;
}

vector<double> FiberRO::readLastLineData(const string &filePath) {
  ifstream file(filePath);
  string line, lastLine;
  while (getline(file, line)) {
    if (!line.empty()) {
      lastLine = line;
    }
  }
  file.close();

  vector<double> data;
  stringstream ss(lastLine);
  string value;
  while (getline(ss, value, ',')) {
    try {
      data.push_back(stod(value));
    } catch (const invalid_argument &) {
    }
  }

  if (data.size() > 3) {
    data.erase(data.begin(), data.end() - 3);
  }
  return data;
}

Matrix3x3 FiberRO::computeRotationMatrix(const vector<double> &Eulerangle) {
  Matrix3x3 Rz = {{{cos(Eulerangle[2]), -sin(Eulerangle[2]), 0},
                   {sin(Eulerangle[2]), cos(Eulerangle[2]), 0},
                   {0, 0, 1}}};

  Matrix3x3 Ry = {{{cos(Eulerangle[1]), 0, sin(Eulerangle[1])},
                   {0, 1, 0},
                   {-sin(Eulerangle[1]), 0, cos(Eulerangle[1])}}};

  Matrix3x3 Rx = {{{1, 0, 0},
                   {0, cos(Eulerangle[0]), -sin(Eulerangle[0])},
                   {0, sin(Eulerangle[0]), cos(Eulerangle[0])}}};

  Matrix3x3 RyRz = {};
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      for (int k = 0; k < 3; ++k) {
        RyRz[i][j] += Ry[i][k] * Rz[k][j];
      }
    }
  }

  Matrix3x3 E = {};
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      for (int k = 0; k < 3; ++k) {
        E[i][j] += Rx[i][k] * RyRz[k][j];
      }
    }
  }
  return E;
}

vector<double> FiberRO::ReadBottomVel() {
  vector<double> velocityRelative =
      readLastLineData(fileBottomVelocityRelative);
  vector<double> omegaRelative = readLastLineData(fileBottomomegaRelative);
  vector<double> Eulerangle = readLastLineData(fileBottomEulerAngle);
  ////////////////////////////////////////////////////Towingpoint Position
  vector<double> R = {-0.0465, 0.0, 0.0};
  ////////////////////////////////////////////////////Towingpoint Position

  Matrix3x3 E = computeRotationMatrix(Eulerangle);

  vector<double> crossProduct = {
      omegaRelative[1] * R[2] - omegaRelative[2] * R[1],
      omegaRelative[2] * R[0] - omegaRelative[0] * R[2],
      omegaRelative[0] * R[1] - omegaRelative[1] * R[0]};

  vector<double> brrdata(3, 0.0);
  for (size_t i = 0; i < 3; ++i) {
    if (i == 1) {
      brrdata[i] = velocityRelative[i] + crossProduct[i];
    } else {
      brrdata[i] = velocityRelative[i] + crossProduct[i];
    }
  }

  vector<double> arrdata(3, 0.0);
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      arrdata[i] += E[i][j] * brrdata[j];
    }
  }
  // swap(arrdata[1], arrdata[2]);
  // swap(arrdata[0], arrdata[1]);

  return arrdata;
}

VectorXd FiberRO::ReadTheLastRow(int index) {
  SPDLOG_DEBUG("Reading the last row of the file with index {}", index);
  Matrix<double, 500, 1> arrcol;
  arrcol.setZero();
  VectorXd TransVal(TotNoV);

  ifstream infile(Outfile);
  if (!infile) {
    cerr << "Open file " + Outfile + " failed" << endl;
    exit(1);
  }

  string line;
  for (int i = 0; i <= index; i++) {
    if (!getline(infile, line)) {
      cerr << "Failed to read line " << i << " from file." << endl;
      infile.close();
      exit(1);
    }
  }

  stringstream sin(line);
  string field;
  int colindex = 0;

  while (getline(sin, field, ',')) {
    if (colindex >= 500) {
      cerr << "Too many columns in the input line." << endl;
      break;
    }
    double dvalue = stod(field);
    arrcol(colindex, 0) = dvalue;
    colindex++;
  }
  infile.close();

  for (int i = 0; i < TotNoV; i++) {
    TransVal(i) = arrcol(i, 0);
  }

  return TransVal;
}

MatrixXd FiberRO::readCSV(int row) {
  SPDLOG_DEBUG("Reading the first {} rows of the file", row);
  ifstream infile(Outfile);
  string line;
  vector<vector<double>> data;

  int currentRow = 0;
  while (getline(infile, line) && currentRow < row) {
    istringstream iss(line);
    vector<double> rowData;

    string valueStr;
    while (getline(iss, valueStr, ',')) {
      double value = stod(valueStr);
      rowData.push_back(value);
    }

    data.push_back(rowData);
    ++currentRow;
  }

  MatrixXd matrix(data.size(), data[0].size());
  for (int i = 0; i < data.size(); ++i) {
    for (int j = 0; j < data[i].size(); ++j) {
      matrix(i, j) = data[i][j];
    }
  }

  return matrix;
}

void FiberRO::OutTopforce(VectorXd v) {
  SPDLOG_DEBUG("Outputting the top force to file");
  ofstream outfileTopforce(Topforceout, ios::trunc);
  MatrixXd aa(1, 3);
  SPDLOG_DEBUG("exmaple log 0");
  aa(0, 0) = v(3); // Forcet
  aa(0, 1) = v(4); // Forcen
  aa(0, 2) = v(5); // Forceb
  SPDLOG_DEBUG("exmaple log 1");
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  MatrixXd tpmat(3, 3);
  tpmat(0, 0) = cos(v(7)) * cos(v(6));
  tpmat(0, 1) = sin(v(7)) * cos(v(6));
  tpmat(0, 2) = -sin(v(6));
  tpmat(1, 0) = -sin(v(7));
  tpmat(1, 1) = cos(v(7));
  SPDLOG_DEBUG("exmaple log 2");
  tpmat(1, 2) = 0;
  tpmat(2, 0) = sin(v(6)) * cos(v(7));
  tpmat(2, 1) = sin(v(6)) * sin(v(7));
  tpmat(2, 2) = cos(v(6));
  SPDLOG_DEBUG("exmaple log 3");
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  MatrixXd bb(1, 3);
  bb = aa * tpmat.inverse();
  outfileTopforce << bb;
  outfileTopforce.close();
}

void FiberRO::OutBottomforce(VectorXd v) {
  SPDLOG_DEBUG("Outputting the bottom force to file");
  ofstream outfileBottomforce(Bottomforceout, ios::trunc);
  MatrixXd aa(1, 3);
  aa(0, 0) = v(493); // Ft
  aa(0, 1) = v(494); // Fn
  aa(0, 2) = v(495); // Fb

  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  MatrixXd tpmat(3, 3);
  tpmat(0, 0) = cos(v(7)) * cos(v(6));
  tpmat(0, 1) = sin(v(7)) * cos(v(6));
  tpmat(0, 2) = -sin(v(6));
  tpmat(1, 0) = -sin(v(7));
  tpmat(1, 1) = cos(v(7));
  tpmat(1, 2) = 0;
  tpmat(2, 0) = sin(v(6)) * cos(v(7));
  tpmat(2, 1) = sin(v(6)) * sin(v(7));
  tpmat(2, 2) = cos(v(6));
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  MatrixXd bb(1, 3);
  bb = aa * tpmat.inverse();
  outfileBottomforce << bb;
  outfileBottomforce.close();
}

vector<double> FiberRO::ReadBottomG() {
  vector<double> arrdata;
  int index = 1;

  ifstream infile(fileObject, ios::in);
  if (!infile) {
    cout << "Could not Open" + fileObject + ".csv" << endl;
    exit(1);
  }

  string line;
  string field;
  for (int i = 0; i <= index; i++) {
    getline(infile, line);
  }
  stringstream sin(line);
  // int colindex = 0;
  while (getline(sin, field, ',')) {
    double dvalue = atof(field.c_str());
    arrdata.push_back(dvalue);
  }
  infile.close();

  return arrdata;
}

vector<double> FiberRO::ReadPhysical() {
  vector<double> arrdata;
  int index = 1;

  ifstream infile(filePhysical, ios::in);
  if (!infile) {
    cout << "Could not Open" + filePhysical + ".csv" << endl;
    exit(1);
  }

  string line;
  string field;
  for (int i = 0; i <= index; i++) {
    getline(infile, line);
  }
  stringstream sin(line);
  // int colindex = 0;
  while (getline(sin, field, ',')) {
    double dvalue = atof(field.c_str());
    arrdata.push_back(dvalue);
  }
  infile.close();

  return arrdata;
}

vector<double> FiberRO::ReadDelta() {
  vector<double> arrdata;
  int index = 1;

  ifstream infile(fileDelta, ios::in);
  if (!infile) {
    cout << "Could not Open" + fileDelta + ".csv" << endl;
    exit(1);
  }

  string line;
  string field;
  for (int i = 0; i <= index; i++) {
    getline(infile, line);
  }
  stringstream sin(line);
  int colindex = 0;
  while (getline(sin, field, ',')) {
    double dvalue = atof(field.c_str());
    arrdata.push_back(dvalue);
  }
  infile.close();

  return arrdata;
}

void FiberRO::Output(const MatrixXd &arr, int rr, int ll) {
  SPDLOG_DEBUG(
      "Outputting the matrix({},{}) to file with {} rows and {} columns",
      arr.rows(), arr.cols(), rr, ll);
  ofstream datafile;
  datafile.open(Outfile.c_str(), ios::out | ios::trunc);
  auto rows = min(static_cast<int>(arr.rows()), rr);
  auto cols = min(static_cast<int>(arr.cols()), ll);
  for (int i = 0; i < rows; i++) {
      SPDLOG_DEBUG("Writing row {}", i);
      for (int j = 0; j < cols; j++) {
        datafile << arr(i, j) << ",";
      }
      for (int j = cols; j < ll; j++) {
        datafile << 0 << ",";
      }
      datafile << endl;
  }
  for (int i = rows; i < rr; i++) {
    for (int j = 0; j < ll; j++) {
      datafile << 0 << ",";
    }
    datafile << endl;
  }
  SPDLOG_DEBUG("Closing the file");
  datafile.close();
}

double FiberRO::flutov(char *filename) {

  double num;
  string s;
  string lastLine;
  ifstream file(filename);
  if (file.is_open()) {
    while (getline(file, s)) {
      lastLine = s;
    }
    file.close();
  } else {
  }
  try {
    int pos = lastLine.find(' ');
    num = stod(lastLine.substr(pos + 1));
  } catch (const std::exception &e) {
  }
  return num;
}