#include "../Head/ReadOut.h"

FiberRO::FiberRO()
{
    //Model 1
    fileTopVel = "../csv./Topvel.csv";
    fileObject = "../csv./TowedObject.csv";

    //Model 2
    fileBottomVelx = (char*)"..\\..\\flu\\Velocity\\sensorvx.out";
    fileBottomVely = (char*)"..\\..\\flu\\Velocity\\sensorvy.out";
    fileBottomVelz = (char*)"..\\..\\flu\\Velocity\\sensorvz.out";//tension

    fileTopVelx = (char*)"..\\..\\flu\\Velocity\\buoyvx.out";//0    z
    fileTopVely = (char*)"..\\..\\flu\\Velocity\\buoyvz.out";//1    x
    fileTopVelz = (char*)"..\\..\\flu\\Velocity\\buoyvy.out";//2    y

    //Model Genrel
    fileWater = "../csv./Water.csv";
    filePhysical = "../csv./Parameters.csv";
    fileDelta = "../csv./Delta.csv";
    

    TotNoV = 500;

    LPCSTR directoryPath = "..\\csv";
    CreateDirectoryA(directoryPath, nullptr);
    Outfile = "..\\csv\\output.csv";

    //The .txt used to read into computation through udf
    Topforceout = "../../flu/topforce.txt";
    Bottomforceout = "../../flu/bottomforce.txt";

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
    int colindex = 0;
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

vector<double> FiberRO::ReadBottomVel()
{

    vector<double> arrdata;

    arrdata.push_back(flutov(fileBottomVelx));
    arrdata.push_back(flutov(fileBottomVelz));
    arrdata.push_back(flutov(fileBottomVely));

    return arrdata;
}

vector<double> FiberRO::ReadTopVel()
{

    vector<double> arrdata;

    arrdata.push_back(flutov(fileTopVelx));
    arrdata.push_back(flutov(fileTopVely));
    arrdata.push_back(flutov(fileTopVelz));

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

void FiberRO::Output(MatrixXd& arr, int rr, int ll)//需要重写
{
    ofstream datafile;
    datafile.open(Outfile.c_str(), ios::out | ios::trunc);

    for(int i = 0; i < rr; i++)
    {
        for(int j = 0; j < ll; j++)
        {
            datafile << arr(i, j) << ",";
        }
        datafile << endl;
    }
    datafile.close();
}

double FiberRO::flutov(char* filename)
{

    double num;
    string s;
    string lastLine;
    ifstream file(filename);
    if (file.is_open())
    {
        while (getline(file, s))
        {
            lastLine = s;
        }
        file.close();
    }
    else
    {}
    try
    {
        int pos = lastLine.find(' ');
        num = stod(lastLine.substr(pos + 1));
    }
    catch (const std::exception& e)
    {}
    return num;
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

    Matrix<double,500,1> arrcol;
    VectorXd TransVal(TotNoV);

    ifstream infile(Outfile, ios::in);
    if (!infile)
    {
      cout << "Open file " + Outfile + " failed" << endl;
      exit(1);
    }

    string line;
    string field;
    for (int i = 0; i <=index; i++)
    {
      getline(infile, line);
    }
    stringstream sin(line);
    int colindex = 0;
    while (getline(sin, field, ','))
    {
      double dvalue = atof(field.c_str());
      arrcol(colindex, 0) = dvalue;
      colindex++;
    }
    infile.close();

    for(int i = 0; i < TotNoV; i++)
    {
      TransVal(i) = arrcol[i];
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


void FiberRO::OutTopforce(VectorXd v)
{
    ofstream outfile(Topforceout, ios::trunc);
    MatrixXd aa(1, 3);
    aa(0, 0) = v(3);//ForceX
    aa(0, 1) = v(4);//ForceY
    aa(0, 2) = v(5);//ForceZ 

    MatrixXd tpmat(3, 3);
    tpmat(0, 0) = cos(v(7))*cos(v(6)); 
    tpmat(0, 1) = -sin(v(7));
    tpmat(0, 2) = sin(v(6))*cos(v(7));
    tpmat(1, 0) = cos(v(6))*sin(v(7));
    tpmat(1, 1) = cos(v(7));
    tpmat(1, 2) = sin(v(6))*sin(v(7));
    tpmat(2, 0) = -sin(v(6));
    tpmat(2, 1) = 0.0;
    tpmat(2, 2) = cos(v(6));

    MatrixXd bb(1, 3);

    bb = aa * tpmat.inverse();
 
    outfile << bb;
    outfile.close();
}


void FiberRO::OutBottomforce(VectorXd v)
{
    ofstream outfile(Bottomforceout, ios::trunc);
    MatrixXd aa(1, 3);
    aa(0, 0) = v(493);//ForceX
    aa(0, 1) = v(494);//ForceY
    aa(0, 2) = v(495);//ForceZ 

    MatrixXd tpmat(3, 3);
    tpmat(0, 0) = cos(v(497))*cos(v(496)); 
    tpmat(0, 1) = -sin(v(497));
    tpmat(0, 2) = sin(v(496))*cos(v(497));
    tpmat(1, 0) = cos(v(496))*sin(v(497));
    tpmat(1, 1) = cos(v(497));
    tpmat(1, 2) = sin(v(496))*sin(v(497));
    tpmat(2, 0) = -sin(v(496));
    tpmat(2, 1) = 0.0;
    tpmat(2, 2) = cos(v(496));

    MatrixXd bb(1, 3);

    bb = aa * tpmat.inverse();
 
    outfile << bb;
    outfile.close();
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