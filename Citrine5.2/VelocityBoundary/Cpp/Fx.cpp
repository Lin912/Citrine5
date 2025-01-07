#include "../Head/Fx.h"
#include <vector>
#include <thread>

Fx::Fx(VectorXd& arr, VectorXd& brr, int index)
    : Yold(arr), Ynew(brr), k(index)
{
    FiberRO a;
    PhysicalData physicalData = ParaReader::ReadAllPhysicalData(a, index);
    A = physicalData.A;
    rho = physicalData.rho;
    d0 = physicalData.d0;
    E = physicalData.E;
    I = physicalData.I;
    M = physicalData.M;
    ma = physicalData.ma;
    Cdt = physicalData.Cdt;
    Cdn = physicalData.Cdn;
    Cdb = physicalData.Cdb;
    pi = physicalData.pi;
    g = physicalData.g;
    Gx = physicalData.Gx;
    Gy = physicalData.Gy;
    Gz = physicalData.Gz;

    Vx = physicalData.Vx;
    Vy = physicalData.Vy;
    Vz = physicalData.Vz;

    Vtx = physicalData.Vtx;
    Vty = physicalData.Vty;
    Vtz = physicalData.Vtz;

    Vbx = physicalData.Vbx;
    Vby = physicalData.Vby;
    Vbz = physicalData.Vbz;

    deltaT = physicalData.deltaT;
    deltaS = physicalData.deltaS;

    Gbx = physicalData.Gbx;
    Gby = physicalData.Gby;
    Gbz = physicalData.Gbz;
    Ax = physicalData.Ax;
    Ay = physicalData.Ay;
    Az = physicalData.Az;
}

Fx::~Fx()
{
}

VectorXd Fx::fx() {
    int numSegments = 50;
    int segmentSize = 10;

    std::vector<VectorXd> Yold_segments(numSegments, VectorXd::Zero(segmentSize));
    std::vector<VectorXd> Ynew_segments(numSegments, VectorXd::Zero(segmentSize));

    //Cutting the Yold and Ynew into segments
    for (int i = 0; i < numSegments; i++) {
        int startIdx = i * segmentSize;
        int length = segmentSize;

        Yold_segments[i] = Yold.segment(startIdx, length);
        Ynew_segments[i] = Ynew.segment(startIdx, length);
    }

    MNQ AA(Yold, Ynew, k);

    std::vector<MatrixXd> Mold(numSegments, MatrixXd(segmentSize, segmentSize));
    std::vector<MatrixXd> Mnew(numSegments, MatrixXd(segmentSize, segmentSize));

    for (int i = 0; i < numSegments; i++) {
        int startIdx = i * segmentSize;
        int length = segmentSize;
        
        Mold[i] = Eigen::SparseMatrix<double>(length, length);
        Mnew[i] = Eigen::SparseMatrix<double>(length, length);    
    
        // 将子矩阵的非零元素插入到新矩阵中
        for (int kold = startIdx; kold < startIdx + length; ++kold) {
            for (Eigen::SparseMatrix<double>::InnerIterator it(AA.Mold(), kold); it; ++it) {
                if (it.row() >= startIdx && it.row() < startIdx + length && it.col() >= startIdx && it.col() < startIdx + length) {
                    Mold[i].coeffRef(it.row() - startIdx, it.col() - startIdx) = it.value();
                }
            }   
        }
        // 对 Mnew 做相同操作
        for (int knew = startIdx; knew < startIdx + length; ++knew) {
            for (Eigen::SparseMatrix<double>::InnerIterator it(AA.Mnew(), knew); it; ++it) {
                if (it.row() >= startIdx && it.row() < startIdx + length && it.col() >= startIdx && it.col() < startIdx + length) {
                    Mnew[i].coeffRef(it.row() - startIdx, it.col() - startIdx) = it.value();
                }
            }
        }
    }


    std::vector<MatrixXd> Nold(numSegments, MatrixXd(segmentSize, segmentSize));
    std::vector<MatrixXd> Nnew(numSegments, MatrixXd(segmentSize, segmentSize));

    for (int i = 0; i < numSegments; i++) {
        int startIdx = i * segmentSize;
        int length = segmentSize;

        // 创建新的稀疏矩阵，存储块
        Nold[i] = Eigen::SparseMatrix<double>(length, length);
        Nnew[i] = Eigen::SparseMatrix<double>(length, length);

        // 通过遍历 Nold 提取块
        for (int knold = startIdx; knold < startIdx + length; ++knold) {
            for (Eigen::SparseMatrix<double>::InnerIterator it(AA.Nold(), knold); it; ++it) {
                if (it.row() >= startIdx && it.row() < startIdx + length && it.col() >= startIdx && it.col() < startIdx + length) {
                    Nold[i].coeffRef(it.row() - startIdx, it.col() - startIdx) = it.value();
                }
            }
        }

        // 通过遍历 Nnew 提取块
        for (int knnew = startIdx; knnew < startIdx + length; ++knnew) {
            for (Eigen::SparseMatrix<double>::InnerIterator it(AA.Nnew(), knnew); it; ++it) {
                if (it.row() >= startIdx && it.row() < startIdx + length && it.col() >= startIdx && it.col() < startIdx + length) {
                    Nnew[i].coeffRef(it.row() - startIdx, it.col() - startIdx) = it.value();
                }
            }
        }
    }

    // 存储 Qold 和 Qnew 向量
    std::vector<VectorXd> Qold(numSegments, VectorXd(segmentSize));
    std::vector<VectorXd> Qnew(numSegments, VectorXd(segmentSize));

    // 填充 Qold 和 Qnew 向量
    for (int i = 0; i < numSegments; ++i) {
        int startIdx = i * segmentSize;
        if (i < numSegments - 1) {
            Qold[i] = AA.Qold().segment(startIdx, segmentSize);
            Qnew[i] = AA.Qnew().segment(startIdx, segmentSize);
        } else {
            Qold[i] = AA.Qold().tail(segmentSize);
            Qnew[i] = AA.Qnew().tail(segmentSize);
        }
    }

    // 创建一个 temp 数组存储每个计算结果
    std::vector<VectorXd> temp(numSegments-1, VectorXd(segmentSize));

    // 计算每个 temp 向量的值
    for (int i = 0; i < numSegments-1; i++) {
            temp[i] = ((Nnew[i] + Nnew[i+1]) * (Ynew_segments[i+1] - Ynew_segments[i]) * deltaT)
                    + ((Nold[i] + Nold[i+1]) * (Yold_segments[i+1] - Yold_segments[i]) * deltaT)
                    + ((Mnew[i+1] + Mold[i+1]) * (Ynew_segments[i+1] - Yold_segments[i+1]) * deltaS)
                    + ((Mnew[i] + Mold[i]) * (Ynew_segments[i] - Yold_segments[i]) * deltaS)
                    + ((Qold[i] + Qold[i+1] + Qnew[i] + Qnew[i+1]) * (deltaT * deltaS));
    }

    // 定义一个计算 BCtemp 的 lambda 函数
        auto calculate_BCtemp = [](const VectorXd& Ynew, double Vx, double Vy, double Vz) {
        VectorXd BCtemp(5);
        double cos_theta = cos(Ynew(6));
        double sin_theta = sin(Ynew(6));
        double cos_phi = cos(Ynew(7));
        double sin_phi = sin(Ynew(7));

        BCtemp(0) = Ynew(0) - (Vx * cos_phi * cos_theta + Vy * cos_theta * sin_phi - Vz * sin_theta);
        BCtemp(1) = Ynew(1) - (Vy * cos_phi - Vx * sin_phi);
        BCtemp(2) = Ynew(2) - (Vx * cos_phi * sin_theta + Vy * sin_theta * sin_phi + Vz * cos_theta);
        BCtemp(3) = Ynew(8);
        BCtemp(4) = Ynew(9);
        return BCtemp;
    };

    // 计算 BCtemp0 和 BCtemp1
    VectorXd BCtemp0 = calculate_BCtemp(Ynew_segments[0], Vtx, Vty, Vtz);
    VectorXd BCtemp1 = calculate_BCtemp(Ynew_segments[numSegments-1], Vbx, Vby, Vbz);

    // 创建最终的返回向量 ret
    VectorXd ret(500);
    
    // 将前五个元素赋值为 BCtemp0，尾部五个元素赋值为 BCtemp1
    ret.head(5) = BCtemp0;
    ret.tail(5) = BCtemp1;

    // 使用循环简化 segment 的赋值操作
    for (int i = 0; i < numSegments-1; i++) {
        int startIdx = 5 + i * segmentSize;
        int length = segmentSize;

        ret.segment(startIdx, length) = temp[i];  // 确保赋值的长度为 length
    }

    return ret;
}



