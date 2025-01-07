#include "../Head/ReadOut.h"

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

FiberRO::~FiberRO()
{}

vector<double> FiberRO::ReadTopVel(int k)
{
    ifstream infile(fileTopVel, ios::in);
    if(!infile)
    {
        cout << "Can not Open" + fileTopVel + ".csv" << endl;
        exit(1); 
    }
    vector<double> arrdata;
    string line;
    string field;
    for (int i = 0; i <= k; i++)
    {
        getline(infile, line);
    }
    stringstream sin(line);
    // int colindex = 0;
    while (getline(sin, field, ','))
    {
        double dvalue = atof(field.c_str());
        arrdata.push_back(dvalue);
    }
    infile.close();
    return arrdata;
}

vector<double> FiberRO::ReadWater(int k)
{
    ifstream infile(fileWater, ios::in);
    if(!infile)
    {
        cout << "Can not Open" + fileWater + ".csv" << endl;
        exit(1); 
    }
    vector<double> arrdata;
    string line;
    string field;
    for (int i = 0; i <= k; i++)
    {
        getline(infile, line);
    }
    stringstream sin(line);
    int colindex = 0;
    while (getline(sin, field, ','))
    {
        double dvalue = atof(field.c_str());
        arrdata.push_back(dvalue);
    }
    infile.close();
    return arrdata;
}

vector<double> FiberRO::ReadPhysical()
{
    vector<double> arrdata;
    int index = 1;//The 2nd row

    ifstream infile(filePhysical, ios::in);
    if (!infile)
    {
        cout << "Could not Open" + filePhysical + ".csv" << endl;
        exit(1);
    }

    string line;
    string field;
    for (int i = 0; i <= index; i++)
    {
        getline(infile, line);
    }
    stringstream sin(line);
    int colindex = 0;
    while (getline(sin, field, ','))
    {
        double dvalue = atof(field.c_str());
        arrdata.push_back(dvalue);
    }
    infile.close();

    return arrdata;
}

void FiberRO::Output(const MatrixXd &arr, int rr, int ll) {
//   SPDLOG_DEBUG(
//       "Outputting the matrix({},{}) to file with {} rows and {} columns",
//       arr.rows(), arr.cols(), rr, ll);
  ofstream datafile;
  datafile.open(Outfile.c_str(), ios::out | ios::trunc);
  auto rows = min(static_cast<int>(arr.rows()), rr);
  auto cols = min(static_cast<int>(arr.cols()), ll);
  for (int i = 0; i < rows; i++) {
    //   SPDLOG_DEBUG("Writing row {}", i);
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
//   SPDLOG_DEBUG("Closing the file");
  datafile.close();
}

vector<double> FiberRO::ReadDelta()
{
    vector<double> arrdata;
    int index = 1;//The 2nd row

    ifstream infile(fileDelta, ios::in);
    if (!infile)
    {
        cout << "Could not Open" + fileDelta + ".csv" << endl;
        exit(1);
    }

    string line;
    string field;
    for (int i = 0; i <= index; i++)
    {
        getline(infile, line);
    }
    stringstream sin(line);
    int colindex = 0;
    while (getline(sin, field, ','))
    {
        double dvalue = atof(field.c_str());
        arrdata.push_back(dvalue);
    }
    infile.close();

    return arrdata;

}

VectorXd FiberRO::ReadTheLastRow(int index)
{

    Matrix<double, 500, 1> arrcol;
    arrcol.setZero();
    VectorXd TransVal(TotNoV); 

    ifstream infile(Outfile);  
    if (!infile)
    {
        cerr << "Open file " + Outfile + " failed" << endl;
        exit(1); 
    }

    string line;
    for (int i = 0; i <= index; i++)
    {
        if (!getline(infile, line))
        {
            cerr << "Failed to read line " << i << " from file." << endl;
            infile.close();
            exit(1); 
        }
    }

    stringstream sin(line);
    string field;
    int colindex = 0;
    
    while (getline(sin, field, ','))
    {
        if (colindex >= 500) 
        {
            cerr << "Too many columns in the input line." << endl;
            break;
        }
        double dvalue = stod(field);  
        arrcol(colindex, 0) = dvalue;  
        colindex++;
    }
    infile.close();  

    for (int i = 0; i < TotNoV; i++)
    {
        TransVal(i) = arrcol(i, 0);
    }

    return TransVal;
}

MatrixXd FiberRO::readCSV(int row) {
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

vector<double> FiberRO::ReadBottomG()
{
    vector<double> arrdata;
    int index = 1;//The 2nd row

    ifstream infile(fileObject, ios::in);
    if (!infile)
    {
        cout << "Could not Open" + fileObject + ".csv" << endl;
        exit(1);
    }

    string line;
    string field;
    for (int i = 0; i <= index; i++)
    {
        getline(infile, line);
    }
    stringstream sin(line);
    int colindex = 0;
    while (getline(sin, field, ','))
    {
        double dvalue = atof(field.c_str());
        arrdata.push_back(dvalue);
    }
    infile.close();

    return arrdata;
}