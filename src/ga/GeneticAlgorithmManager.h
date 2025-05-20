#include <iostream>
#include <chrono>
#include <vector>
#include <algorithm>
#include <random>


class GeneticAlgorithm
{
public:
    GeneticAlgorithm(int popSize_, int maxGen_, int tourSize_, double mutRate_, double s_, double s0_, double W_, double H_, const std::vector<double>& L_, const std::vector<double>& DY_, const std::vector<std::vector<double>>& hij_, const std::vector<double>& hi0_, const std::vector<std::vector<double>>& Pij_, const std::vector<std::vector<double>>& Qij_, const std::vector<double>& de_);
	~GeneticAlgorithm();

    void initializePopulation();
    int tournamentSelection(const std::vector<double>& costs);
    void crossoverOX(const std::vector<int>& p1, const std::vector<int>& p2,
        std::vector<int>& c1, std::vector<int>& c2);
    void mutate(std::vector<int>& chromo);
    void layoutPositions(const std::vector<int>& order,
        std::vector<std::pair<double, double>>& pos,
        int& rowCount, double& lastRowMaxWidth);
    double fitness(const std::vector<int>& order);
    void run();

private:
    int popSize, maxGen, tourSize;
    double mutRate;
    int n;                                           // �豸����
    double s, s0, W, H;
    std::vector<double> L, DY, hi0, de;
    std::vector<std::vector<double>> hij, Pij, Qij;
    std::vector<std::vector<int>> population;                   // ��ǰ��Ⱥ��popSize �� n
    std::vector<int> bestOrder;                           // �������
    std::vector<std::pair<double, double>> bestPos;             // ��Ѳ����������豸������
    int bestRowCount;
    double bestLastRowMaxWidth;
    std::mt19937 rng;
};

