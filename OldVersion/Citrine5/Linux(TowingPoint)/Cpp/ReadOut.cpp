#include "../Head/ReadOut.h"

FiberRO::FiberRO()
{
    //Model 1
    fileTopVel = "../csv/TopVel.csv";
    fileObject = "../csv/TowedObject.csv";

    //Model 2
    //fileBottomVelocity = "../../Star/Velocity.csv";

    fileBottomVelocityRelative = "../../Star/VelocityRelative.csv";
    fileBottomomegaRelative = "../../Star/omegaRelative.csv";
    fileBottomEulerAngle = "../../Star/EulerAngle.csv";


    //Model Genrel
    fileWater = "../csv/Water.csv";
    filePhysical = "../csv/Parameters.csv";
    fileDelta = "../csv/Delta.csv";
    Outfile = "../csv/output.csv";

    std::filesystem::create_directory("../csv");

    TotNoV = 500;

    //The .txt used to read into computation through udf
    Topforceout = "../../Star/topforce.txt";
    Bottomforceout = "../../Star/bottomforce.txt";


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

vector<double> FiberRO::readLastLineData(const string& filePath) {
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
        } catch (const invalid_argument&) {}
    }

    if (data.size() > 3) {
        data.erase(data.begin(), data.end() - 3);
    }
    return data;
}

Matrix3x3 FiberRO::computeRotationMatrix(const vector<double>& Eulerangle) {
    Matrix3x3 Rz = {{
        {cos(Eulerangle[2]), -sin(Eulerangle[2]), 0},
        {sin(Eulerangle[2]), cos(Eulerangle[2]), 0},
        {0, 0, 1}
    }};

    Matrix3x3 Ry = {{
        {cos(Eulerangle[1]), 0, sin(Eulerangle[1])},
        {0, 1, 0},
        {-sin(Eulerangle[1]), 0, cos(Eulerangle[1])}
    }};

    Matrix3x3 Rx = {{
        {1, 0, 0},
        {0, cos(Eulerangle[0]), -sin(Eulerangle[0])},
        {0, sin(Eulerangle[0]), cos(Eulerangle[0])}
    }};

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
    vector<double> velocityRelative = readLastLineData(fileBottomVelocityRelative);
    vector<double> omegaRelative = readLastLineData(fileBottomomegaRelative);
    vector<double> Eulerangle = readLastLineData(fileBottomEulerAngle);
    vector<double> R = {0.0, 0.0, 0.0};

    Matrix3x3 E = computeRotationMatrix(Eulerangle);

    vector<double> crossProduct = {
        omegaRelative[1] * R[2] - omegaRelative[2] * R[1],
        omegaRelative[2] * R[0] - omegaRelative[0] * R[2],
        omegaRelative[0] * R[1] - omegaRelative[1] * R[0]
    };

    vector<double> brrdata(3, 0.0);
    for (size_t i = 0; i < 3; ++i) {
        brrdata[i] = velocityRelative[i] + crossProduct[i];
    }


    vector<double> arrdata(3, 0.0);
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            arrdata[i] += E[i][j] * brrdata[j];
        }
    }
    swap(arrdata[1], arrdata[2]);

    return arrdata;
}




VectorXd FiberRO::ReadTheLastRow(int index)
{
    Matrix<double, 500, 1> arrcol;  // 定义一个 500x1 的矩阵
    VectorXd TransVal(TotNoV);  // 定义一个动态大小的向量，大小为 TotNoV

    ifstream infile(Outfile);  // 打开文件
    if (!infile)
    {
        cerr << "Open file " + Outfile + " failed" << endl;
        exit(1);  // 错误处理
    }

    string line;
    // 跳过前面 index 行，读取第 index 行的数据
    for (int i = 0; i <= index; i++)
    {
        if (!getline(infile, line))
        {
            cerr << "Failed to read line " << i << " from file." << endl;
            infile.close();
            exit(1);  // 如果文件没有足够的行数，则报错退出
        }
    }

    stringstream sin(line);
    string field;
    int colindex = 0;
    
    // 从字符串中解析每个以逗号分隔的数值
    while (getline(sin, field, ','))
    {
        if (colindex >= 500)  // 确保不超过 500 列的限制
        {
            cerr << "Too many columns in the input line." << endl;
            break;
        }
        double dvalue = stod(field);  // 使用 std::stod 进行转换
        arrcol(colindex, 0) = dvalue;  // 将值存入 arrcol 矩阵
        colindex++;
    }
    infile.close();  // 关闭文件

    // 将 arrcol 中的前 TotNoV 个值复制到 TransVal 向量中
    for (int i = 0; i < TotNoV; i++)
    {
        TransVal(i) = arrcol(i, 0);
    }

    return TransVal;
}



//VectorXd FiberRO::ReadTheLastRow(int index)
//{

    // Matrix<double,500,1> arrcol;
    // VectorXd TransVal(TotNoV);

    // ifstream infile(Outfile, ios::in);
    // if (!infile)
    // {
    //   cout << "Open file " + Outfile + " failed" << endl;
    //   exit(1);
    // }

    // string line;
    // string field;
    // for (int i = 0; i <=index; i++)
    // {
    //   getline(infile, line);
    // }
    // stringstream sin(line);
    // int colindex = 0;
    // while (getline(sin, field, ','))
    // {
    //   double dvalue = atof(field.c_str());
    //   arrcol(colindex, 0) = dvalue;
    //   colindex++;
    // }
    // infile.close();

    // for(int i = 0; i < TotNoV; i++)
    // {
    //   TransVal(i) = arrcol[i];
    // }
    // return TransVal;
//}

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