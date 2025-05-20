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
        std::cout << "�޷����ļ�: " << path << endl;
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
        std::cout << "�޷����ļ�: " << path << endl;
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
        if (line.find(',') == std::string::npos)         // �����һ��ֻ��һ�����֣�ֱ�Ӷ�ȡ������ֻȡ���е�һ������
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
    QString appDir = QCoreApplication::applicationDirPath();        // ��ȡ��ִ���ļ�����Ŀ¼
    
    // �豸�ߴ磺device_size.csv��ÿ�У�[length, width]
    std::vector< std::vector<double>> devSize = readCSV_double(appDir, "device_size.csv");
    int n = static_cast<int>(devSize.size());
    std::vector<double> L(n), DY(n);             // ÿ���豸�ĳ�����
    for (int i = 0; i < n; ++i) 
    {
        L[i] = devSize[i][0];
        DY[i] = devSize[i][1];
    }

    // ��ȡ Pij��n��n����Qij��n��n����hij��n��n��
    std::vector< std::vector<double>> Pij = readCSV_double(appDir, "Pij.csv");
    std::vector< std::vector<double>> Qij = readCSV_double(appDir, "Qij.csv");
    std::vector< std::vector<double>> hij = readCSV_double(appDir, "hij.csv");
    // ��ȡ hi0������ n����de������ n��
    std::vector<double> hi0 = readCSV_vector(appDir, "hi0.csv");
    std::vector<double> de = readCSV_vector(appDir, "de.csv");

    if ((int)hi0.size() != n || (int)de.size() != n
        || (int)Pij.size() != n || (int)Pij[0].size() != n
        || (int)Qij.size() != n || (int)Qij[0].size() != n
        || (int)hij.size() != n || (int)hij[0].size() != n) {
        std::cout << "�����ļ��ߴ����豸������ƥ�䣬���� CSV �ļ���ʽ��" << endl;
        return EXIT_FAILURE;
    }

    // --- 2. GA �������� ---
    int popSize = 80;      // ��Ⱥ��ģ
    int maxGen = 200;     // ������
    double mutRate = 0.2;  // �������
    int tourSize = 10;     // ��������ģ

    double s = 5000.0;    // �о�
    double s0 = 4000.0;    // ��һ����ײ��߽����
    double W = 36000.0;   // ���䳤�� (X �᷶Χ)
    double H = 48000.0;   // ������ (Y �᷶Χ)

    // --- 3. ���첢�����Ŵ��㷨 ---
    GeneticAlgorithm ga(popSize, maxGen, tourSize, mutRate, s, s0, W, H, L, DY,hij, hi0,  Pij, Qij, de);
    ga.run();

    return a.exec();
}



