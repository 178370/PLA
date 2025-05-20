#include <iostream>
#include "ga/GeneticAlgorithmManager.h"
#include <QCoreApplication>
#include <QDebug>
#include <string>
#include <fstream>
#include <sstream>
#include <QDir>


#include <filesystem>
#include <iostream>
#include <fstream>
#include <cerrno>
#include <cstring>






std::vector<std::vector<double>> readCSV_double(const QString &appDir, const QString& filename)
{
    
    QString qpath = appDir.trimmed()
        .remove('\r')
        .remove('\n')
        + QDir::separator()
        + filename;
    QByteArray arr_dev = qpath.toUtf8();
    std::string path(arr_dev.constData());

    std::ifstream in(path);
    if (!in.is_open()) 
    {
        std::cout << "无法打开文件: " << path << endl;
        exit(EXIT_FAILURE);
    }
    std::vector< std::vector<double>> matrix;
    std::string line;
    while (getline(in, line)) 
    {
        if (line.empty()) 
            continue;
        std::stringstream ss(line);

        std::vector<double> row;
        std::string token;
        while (getline(ss, token, ',')) 
        {
            row.push_back(stod(token));
        }
        matrix.push_back(row);
    }
    in.close();
    return matrix;
}

std::vector<double> readCSV_vector(const QString& appDir, const QString& filename)
{
    QString qpath = appDir.trimmed()
        .remove('\r')
        .remove('\n')
        + QDir::separator()
        + filename;
    QByteArray arr_hi0 = qpath.toUtf8();
    std::string path(arr_hi0.constData());
    std::ifstream in(path);
    if (!in.is_open()) 
    {
        std::cout << "无法打开文件: " << path << endl;
        exit(EXIT_FAILURE);
    }
    std::vector<double> vec;
    std::string line;
    while (getline(in, line)) 
    {
        if (line.empty()) 
            continue;
        std::stringstream ss(line);
        std::string token;
        if (line.find(',') == std::string::npos)         // 如果这一行只有一个数字，直接读取；否则只取该行第一个数字
        {
            vec.push_back(stod(line));
        }
        else 
        {
            getline(ss, token, ',');
            vec.push_back(stod(token));
        }
    }
    in.close();
    return vec;
}











int main(int argc, char* argv[])
{
    QCoreApplication a(argc, argv);
    QString appDir = QCoreApplication::applicationDirPath();        // 获取可执行文件所在目录
    
    // 设备尺寸：device_size.csv，每行：[length, width]
    std::vector< std::vector<double>> devSize = readCSV_double(appDir, "device_size.csv");
    int n = static_cast<int>(devSize.size());
    std::vector<double> L(n), DY(n);             // 每个设备的长，宽
    for (int i = 0; i < n; ++i) 
    {
        L[i] = devSize[i][0];
        DY[i] = devSize[i][1];
    }

    // 读取 Pij（n×n），Qij（n×n），hij（n×n）
    std::vector< std::vector<double>> Pij = readCSV_double(appDir, "Pij.csv");
    std::vector< std::vector<double>> Qij = readCSV_double(appDir, "Qij.csv");
    std::vector< std::vector<double>> hij = readCSV_double(appDir, "hij.csv");
    // 读取 hi0（长度 n）、de（长度 n）
    std::vector<double> hi0 = readCSV_vector(appDir, "hi0.csv");
    std::vector<double> de = readCSV_vector(appDir, "de.csv");

    if ((int)hi0.size() != n || (int)de.size() != n
        || (int)Pij.size() != n || (int)Pij[0].size() != n
        || (int)Qij.size() != n || (int)Qij[0].size() != n
        || (int)hij.size() != n || (int)hij[0].size() != n) {
        std::cout << "输入文件尺寸与设备数量不匹配，请检查 CSV 文件格式。" << endl;
        return EXIT_FAILURE;
    }

    // --- 2. GA 参数定义 ---
    int popSize = 80;      // 种群规模
    int maxGen = 200;     // 最大代数
    double mutRate = 0.2;  // 变异概率
    int tourSize = 10;     // 锦标赛规模

    double s = 5000.0;    // 行距
    double s0 = 4000.0;    // 第一行与底部边界距离
    double W = 36000.0;   // 车间长度 (X 轴范围)
    double H = 48000.0;   // 车间宽度 (Y 轴范围)

    // --- 3. 构造并运行遗传算法 ---
    GeneticAlgorithm ga(popSize, maxGen, tourSize, mutRate, s, s0, W, H, L, DY,hij, hi0,  Pij, Qij, de);
    ga.run();

    return a.exec();
}



